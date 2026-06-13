#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import dataclasses
import logging
import math

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

from rmgpy.solver.polymer import HybridPolymerSystem, MassTransferConfig, PolymerPoolConfig


def _spc(smiles: str, label: str) -> Species:
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def _two_pool_species():
    """Species + mask for a two-pool system: A (mu 1-3), B (mu 5-7), gas G at 8."""
    sp = {
        "A": _spc("CCCC", "A"),
        "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"), "A_mu2": _spc("C#N", "A_mu2"),
        "B": _spc("CCCCC", "B"),
        "B_mu0": _spc("CCO", "B_mu0"), "B_mu1": _spc("CC=O", "B_mu1"), "B_mu2": _spc("CC#N", "B_mu2"),
        "G": _spc("[CH3]", "G"),
    }
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
            sp["B"], sp["B_mu0"], sp["B_mu1"], sp["B_mu2"], sp["G"]]
    mask = np.array([False] * 8 + [True], dtype=bool)
    return sp, core, mask


def _two_pool_rs(rxn, core, mask, mom_a, mom_b):
    pool_a = PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(1, 2, 3), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    pool_b = PolymerPoolConfig(label="B", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(5, 6, 7), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0}, V_poly=1.0,
        polymer_pools=[pool_a, pool_b], mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments={"A": mom_a, "B": mom_b}, termination=[],
    )
    rs.initialize_model(core, [rxn], [], [])
    return rs

_KIN = dict(kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)

# Reversible kinetics for the thermo reference-state tripwire fixtures (spec
# 2026-06-11): the tripwire keys on rxn.reversible; _KIN is reversible=False.
_REV_KIN = dict(kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
                reversible=True)

# Provenance comment strings, pinned against rmgpy/data/thermo.py (:232/:1818/
# :2237/:2845) and a real EPDM chem_annotated.inp.
_GAV_COMMENT = "Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH)"
_LIB_COMMENT = "Thermo library: primaryThermoLibrary"


def _trivial_nasa(comment=""):
    """Minimal valid NASA thermo (Ne-like constant 2.5R heat capacity) so
    reversible fixtures survive generate_rate_coefficients' Keq path on ANY
    head (identical thermo on both sides -> finite Keq, no explosion). The
    ``comment`` carries the provenance string the tripwire's sensor reads.
    The tripwire's U itself never reads thermo -- it is MW-only by design."""
    poly = NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35],
                          Tmin=(200, "K"), Tmax=(6000, "K"))
    return NASA(polynomials=[poly], Tmin=(200, "K"), Tmax=(6000, "K"), comment=comment)


def _sackur_tetrode_decades(mw_kg_mol, T):
    """Test-local exact Sackur-Tetrode S_trans at P0 = 1e5 Pa, in decades
    (S/(R ln10)). INDEPENDENT recompute for liveness pins -- deliberately not
    imported from production code, so a broken production formula cannot
    vouch for itself."""
    m = mw_kg_mol / constants.Na
    s = constants.R * (math.log((2.0 * math.pi * m * constants.kB * T / constants.h ** 2) ** 1.5
                                * constants.kB * T / 1.0e5) + 2.5)
    return s / (constants.R * math.log(10.0))


def _refstate_pool_species(thermo_comment=_GAV_COMMENT):
    """One-pool fixture for the reference-state tripwire: melt proxy A
    (butane) at 0, mu dummies 1-3, gas products G1 (CH4) and G2 (propene) at
    4-5. A <=> G1 + G2 is mass-balanced (58.12 = 16.04 + 42.08). All species
    carry trivial NASA thermo (see _trivial_nasa)."""
    sp = {
        "A": _spc("CCCC", "A"),
        "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"), "A_mu2": _spc("C#N", "A_mu2"),
        "G1": _spc("C", "G1"), "G2": _spc("C=CC", "G2"),
    }
    for s in sp.values():
        s.thermo = _trivial_nasa(thermo_comment)
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["G1"], sp["G2"]]
    mask = np.array([False] * 4 + [True, True], dtype=bool)
    return sp, core, mask


def _refstate_rs(core, rxns, mask, pools, moments, rs_kwargs=None):
    """Build + initialize a HybridPolymerSystem for the tripwire fixtures.
    core[-1] must be a gas species (it seeds initial_mole_fractions)."""
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[-1]: 0.0}, V_poly=1.0,
        polymer_pools=pools, mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments=moments, termination=[],
        **(rs_kwargs or {}),
    )
    rs.initialize_model(core, rxns, [], [])
    return rs


def _gate_pool_config(monomer_mw_g_mol=28.0):
    """PolymerPoolConfig for the spawn-gate fixtures.

    ``monomer_mw_g_mol`` is added by the mass-flux-spawn-gate change (spec
    2026-06-10 §4.1). The field-presence guard below is deliberate RED-FIRST
    scaffolding, NOT a compatibility shim: it lets the Task-1 integrated
    tripwire run on pre-change HEAD and die at the GROSS-ARRAY assertion
    (the born-dead mechanism, spec §3.1) instead of at fixture construction.
    Once the field lands, the guard always takes the kwargs branch.
    """
    kwargs = dict(label="A", xs=2, explicit_dp_to_species_index={},
                  mu_indices=(1, 2, 3), monomer_poly_index=None)
    if any(f.name == "monomer_mw_g_mol"
           for f in dataclasses.fields(PolymerPoolConfig)):
        kwargs["monomer_mw_g_mol"] = monomer_mw_g_mol
    return PolymerPoolConfig(**kwargs)


def _one_pool_gate_species(rep_smiles="CCO"):
    """Species + mask for the one-pool spawn-gate fixture: canonical proxy A
    at 0, A_mu0/1/2 at 1-3, ordinary POLYMER-PHASE species R at 4 (stands in
    for an absorbed representative — representative status is a python/ledger
    concept; to the solver it is ANY ordinary species produced by
    pool-touching chemistry), gas G at 5."""
    sp = {
        "A": _spc("CCCC", "A"),
        "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"), "A_mu2": _spc("C#N", "A_mu2"),
        "R": _spc(rep_smiles, "R"),
        "G": _spc("[CH3]", "G"),
    }
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["R"], sp["G"]]
    mask = np.array([False] * 5 + [True], dtype=bool)
    return sp, core, mask


def _one_pool_gate_rs(rxn, core, mask, moments, monomer_mw_g_mol=28.0):
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[5]: 0.0}, V_poly=1.0,
        polymer_pools=[_gate_pool_config(monomer_mw_g_mol)], mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments={"A": moments}, termination=[],
    )
    rs.initialize_model(core, [rxn], [], [])
    return rs


class TestHybridPolymerReactor:
    def test_phase_pure_gas_reaction_molar_balance(self):
        """
        For a first-order gas reaction A -> B (irreversible),
        the solver should give dnA/dt = -k*nA and dnB/dt = +k*nA (mol/s).
        """
        A = _spc("C", "A")       # gas
        B = _spc("[CH3]", "B")   # gas

        rxn = Reaction(
            reactants=[A],
            products=[B],
            kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False,
        )

        core_species = [A, B]
        core_reactions = [rxn]

        gas_species_mask = np.array([True, True], dtype=bool)

        T = 1000.0
        P = 1.0e5

        rxn_system = HybridPolymerSystem(
            T=T,
            P=P,
            initial_mole_fractions={A: 1.5, B: 0.0},  # interpreted as moles for this reactor
            V_poly=1.0,  # unused here (no polymer species)
            polymer_pools=[],
            mass_transfer=[],
            gas_species_mask=gas_species_mask,
            constant_gas_volume=False,
            termination=[],
        )

        rxn_system.initialize_model(core_species, core_reactions, [], [])

        # residual returns (dn_dt - dydt); use dydt=0 to get dn_dt directly
        dn_dt = rxn_system.residual(0.0, rxn_system.y, np.zeros_like(rxn_system.y))[0]

        k = rxn.get_rate_coefficient(T, P)
        assert abs(dn_dt[0] - (-k * 1.5)) <= 1e-12
        assert abs(dn_dt[1] - (+k * 1.5)) <= 1e-12

    def test_cross_phase_core_product_disqualifies_reaction(self):
        """
        "Produce-then-Transfer" policy: cross-phase core products disqualify a reaction (rate=0).
        """
        A = _spc("C", "A_gas")     # gas reactant
        C = _spc("O=C=O", "C_pol") # polymer-phase product (core)

        rxn = Reaction(
            reactants=[A],
            products=[C],
            kinetics=Arrhenius(A=(10.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False,
        )

        core_species = [A, C]
        core_reactions = [rxn]

        gas_species_mask = np.array([True, False], dtype=bool)

        rxn_system = HybridPolymerSystem(
            T=1000.0,
            P=1.0e5,
            initial_mole_fractions={A: 1.0},
            V_poly=1.0,
            polymer_pools=[],
            mass_transfer=[],
            gas_species_mask=gas_species_mask,
            constant_gas_volume=False,
            termination=[],
        )
        rxn_system.initialize_model(core_species, core_reactions, [], [])

        dn_dt = rxn_system.residual(0.0, rxn_system.y, np.zeros_like(rxn_system.y))[0]

        # Reaction should be disqualified, so no source/sink from kinetics
        assert abs(dn_dt[0]) <= 1e-20
        assert abs(dn_dt[1]) <= 1e-20

    def test_mass_transfer_sign_convention_poly_to_gas(self):
        """
        Mass transfer: J = kLa * (Cp - K*Cg), dn = J*V_poly.
        Sign convention (per code comment): J > 0 => net poly -> gas (gas gains moles).
        """
        Gg = _spc("N#N", "G_gas")     # gas
        Gp = _spc("O=C=O", "G_poly")  # polymer-dissolved counterpart

        core_species = [Gg, Gp]
        gas_species_mask = np.array([True, False], dtype=bool)

        # Pick K small so Cp - K*Cg > 0 even if Cg is large (ideal gas)
        mt = MassTransferConfig(gas_index=0, poly_index=1, K=0.01, kLa=1.0)

        rxn_system = HybridPolymerSystem(
            T=1000.0,
            P=1.0e5,
            initial_mole_fractions={Gg: 1.0},
            V_poly=1.0,
            polymer_pools=[],
            mass_transfer=[mt],
            gas_species_mask=gas_species_mask,
            constant_gas_volume=False,
            termination=[],
        )
        rxn_system.initialize_model(core_species, [], [], [])

        # Set polymer moles explicitly in the state vector
        y = rxn_system.y.copy()
        y[1] = 2.0  # polymer moles

        dn_dt = rxn_system.residual(0.0, y, np.zeros_like(y))[0]

        # Recompute expected dn based on the same volume conventions
        V_gas = constants.R * rxn_system.T.value_si * (y[0]) / rxn_system.P.value_si
        Cg = y[0] / V_gas
        Cp = y[1] / rxn_system.V_poly
        J = mt.kLa * (Cp - mt.K * Cg)
        expected_dn = J * rxn_system.V_poly

        assert abs(dn_dt[0] - expected_dn) <= 1e-10
        assert abs(dn_dt[1] + expected_dn) <= 1e-10  # polymer loses what gas gains

    def test_tail_handshake_generates_explicit_boundary_species(self):
        """
        If tail is present (mu0 > TAIL_CONC_MIN and mean DP > xs),
        the hybrid handshake should create positive flux into explicit DP=xs species.
        """
        # Gas dummy (keeps gas volume well-defined, but doesn't react)
        Inert = _spc("N#N", "N2")

        # Polymer explicit boundary species (DP = xs)
        P2 = _spc("CC", "P2")  # explicit DP=2

        # Moment "species" placeholders (just indices in solver state)
        Mu0 = _spc("CO", "Mu0")
        Mu1 = _spc("CO", "Mu1")
        Mu2 = _spc("CO", "Mu2")

        core_species = [Inert, P2, Mu0, Mu1, Mu2]

        # gas mask: only Inert is gas
        gas_species_mask = np.array([True, False, False, False, False], dtype=bool)

        pool = PolymerPoolConfig(
            label="poly",
            xs=2,
            explicit_dp_to_species_index={2: 1},  # DP=2 -> P2 index
            mu_indices=(2, 3, 4),                 # Mu0, Mu1, Mu2 indices
            monomer_poly_index=None,
            k_scission=0.0,
            k_unzip=0.1,  # enables handshake
            tail_kinetics=None,
        )

        # Provide tail moments as *moles of moments* (Mu_k = mu_k * V_poly).
        # Choose V_poly=1.0 so mu_k == Mu_k.
        initial_polymer_moments = {"poly": (1.0, 5.0, 30.0)}  # mean=5 > xs=2

        rxn_system = HybridPolymerSystem(
            T=800.0,
            P=1.0e5,
            initial_mole_fractions={Inert: 1.0},
            V_poly=1.0,
            polymer_pools=[pool],
            mass_transfer=[],
            gas_species_mask=gas_species_mask,
            constant_gas_volume=False,
            initial_polymer_moments=initial_polymer_moments,
            initial_explicit_species={"poly": {2: 0.0}},
            termination=[],
        )
        rxn_system.initialize_model(core_species, [], [], [])

        dn_dt = rxn_system.residual(0.0, rxn_system.y, np.zeros_like(rxn_system.y))[0]

        # Handshake should:
        #  - increase explicit DP=xs moles
        #  - decrease Mu0 (only handshake changes Mu0 in the default kinetics)
        assert dn_dt[1] > 0.0
        assert dn_dt[2] < 0.0

    def test_end_group_reaction_scales_by_mu0_not_mu1(self):
        """
        A proxy reaction flagged ``is_end_group_reaction`` must scale by chain-end
        density (mu0); an unflagged proxy reaction scales by monomer-unit density
        (mu1). End-group reactions occur at chain ends, so their rate is set by the
        number of ends (mu0), not the number of repeat units (mu1).

        Verified on a reactant-proxy reaction, where get_reaction_rates yields
        rate = kf * moment / V_poly (the [THE HIJACK] block, polymer.pyx ~1262).
        With mu0 != mu1 the flagged rate is smaller by exactly mu0/mu1. This pins
        the wiring of the previously-dead is_end_group_reaction flag.
        """
        Proxy = _spc("CCCC", "poly")      # pool proxy (label == pool.label)
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        Prod = _spc("[CH3]", "P")         # gas product sink

        core_species = [Proxy, Mu0, Mu1, Mu2, Prod]
        gas_species_mask = np.array([False, False, False, False, True], dtype=bool)

        mu0, mu1, mu2 = 1.0, 5.0, 30.0    # mu0 != mu1 so the two scalings differ
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        initial_polymer_moments = {"poly": (mu0, mu1, mu2)}

        def _proxy_rate(flag):
            rxn = Reaction(
                reactants=[Proxy], products=[Prod],
                kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
                reversible=False,
            )
            rxn.is_end_group_reaction = flag
            rs = HybridPolymerSystem(
                T=800.0, P=1.0e5, initial_mole_fractions={Prod: 0.0}, V_poly=1.0,
                polymer_pools=[pool], mass_transfer=[],
                gas_species_mask=gas_species_mask.copy(), constant_gas_volume=False,
                initial_polymer_moments=initial_polymer_moments, termination=[],
            )
            rs.initialize_model(core_species, [rxn], [], [])
            return rs.get_reaction_rates(rs.y)[0], rxn.get_rate_coefficient(800.0, 1.0e5)

        rate_mu1, kf = _proxy_rate(False)   # default: mu1 scaling
        rate_mu0, _ = _proxy_rate(True)     # flagged: mu0 scaling

        assert np.isclose(rate_mu1, kf * mu1 / 1.0)
        assert np.isclose(rate_mu0, kf * mu0 / 1.0)
        assert np.isclose(rate_mu0 / rate_mu1, mu0 / mu1)

    def test_random_scission_moment_derivatives(self):
        """
        Random backbone scission must satisfy the analytic discrete-bond
        (Ziff-McGrady) moment closure (see docs/multi_pool_design.md §5):

            dmu0/dt = k_s * (mu1 - mu0)     (one new chain per bond broken)
            dmu1/dt = 0                     (monomer units conserved -> mass balance)
            dmu2/dt = (k_s/3) * (mu1 - mu3)

        with the log-Lagrange closure mu3 = mu0 * (mu2/mu1)**3.

        Regression guards:
        - The very first implementation used ``+k_s*mu1*(mu3/mu1 - mu2)`` for
          dmu2/dt (dimensionally inconsistent, wrong sign).
        - The mu0 source MUST be ``k_s*(mu1 - mu0)``, not ``k_s*mu1``. The bare
          +k_s*mu1 form lets mu0 grow past mu1, the moment state leaves the
          realizable cone, and the mu3 closure blows up to a DASSL singularity
          (see ~/Projects/EPDM/polymer_branch_scission_handoff.md, BUG 1).
        """
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "Mu0")
        Mu1 = _spc("C=O", "Mu1")
        Mu2 = _spc("C#N", "Mu2")

        core_species = [Inert, Mu0, Mu1, Mu2]
        gas_species_mask = np.array([True, False, False, False], dtype=bool)

        k_s = 0.3
        pool = PolymerPoolConfig(
            label="poly",
            xs=2,
            explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3),
            monomer_poly_index=None,
            k_scission=k_s,
            k_unzip=0.0,  # isolate scission; no handshake/unzip flux
            tail_kinetics=None,
        )

        # V_poly=1 so moles-of-moment == concentration of moment.
        V_poly = 1.0
        mu0, mu1, mu2 = 1.0, 5.0, 30.0  # mean DP=5, polydisperse
        initial_polymer_moments = {"poly": (mu0, mu1, mu2)}

        rxn_system = HybridPolymerSystem(
            T=800.0,
            P=1.0e5,
            initial_mole_fractions={Inert: 1.0},
            V_poly=V_poly,
            polymer_pools=[pool],
            mass_transfer=[],
            gas_species_mask=gas_species_mask,
            constant_gas_volume=False,
            initial_polymer_moments=initial_polymer_moments,
            termination=[],
        )
        rxn_system.initialize_model(core_species, [], [], [])

        dn_dt = rxn_system.residual(0.0, rxn_system.y, np.zeros_like(rxn_system.y))[0]

        mu3 = mu0 * (mu2 / mu1) ** 3  # log-Lagrange closure used by the solver

        expected_dmu0 = k_s * (mu1 - mu0) * V_poly
        expected_dmu1 = 0.0
        expected_dmu2 = k_s * (mu1 - mu3) / 3.0 * V_poly

        assert np.isclose(dn_dt[1], expected_dmu0, rtol=1e-9)
        assert np.isclose(dn_dt[2], expected_dmu1, atol=1e-12)
        assert np.isclose(dn_dt[3], expected_dmu2, rtol=1e-9)
        # mu0 source is self-limiting: positive while mu0 < mu1, so it can never
        # drive mu0 past mu1 (realizability preserved).
        assert dn_dt[1] > 0.0
        # mu3 dominates mu1 for a polydisperse pool, so mu2 narrows (decreases).
        assert dn_dt[3] < 0.0

    def test_residual_stays_finite_for_extreme_moment_states(self):
        """
        Resurrection robustness (handoff item): the polymer moment RHS must never
        emit NaN/Inf, even for out-of-cone or overflow-broad moment states. A
        non-finite dn_dt surfaces in the solver as 'nans in moles' -> DASxError ->
        model resurrection, which fails UNRECOVERABLY ('invalid_objects could not
        be filled') when the cause is the core moment dynamics rather than a
        promotable edge species (solver/base.pyx ~748-782).

        Two guards keep the RHS bounded: _safe_mu3_from_mu012 returns 0.0 outside
        the realizable cone (mu1<mu0) and for degenerate moments, and Inf on
        log-overflow; the scission dmu2 term is then gated by np.isfinite(mu3)
        (polymer.pyx ~1052). This test pins BOTH so a future refactor can't
        silently reintroduce the unrecoverable failure.
        """
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        core_species = [Inert, Mu0, Mu1, Mu2]
        gas_species_mask = np.array([True, False, False, False], dtype=bool)

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.3, k_unzip=0.1, tail_kinetics=None,  # exercise both moment terms
        )
        rxn_system = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[], gas_species_mask=gas_species_mask,
            constant_gas_volume=False, initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            termination=[],
        )
        rxn_system.initialize_model(core_species, [], [], [])

        # (mu0, mu1, mu2) moles. V_poly=1 so moles == concentration.
        pathological_states = [
            (1.0, 5.0, 30.0),         # normal/realizable
            (21.17, 0.0139, 380.5),   # cone violation mu1 < mu0 (handoff BUG 1 state)
            (0.0, 5.0, 30.0),         # mu0 = 0
            (1.0e-3, 1.0, 1.0e120),   # broad distribution -> mu3 closure overflows to Inf
            (1.0e100, 1.0e100, 1.0e300),  # all-huge near double-max
            (-1.0, -2.0, -3.0),       # negatives (clamped by max(0,.))
            (0.0, 0.0, 0.0),          # fully degenerate
        ]
        for mu0, mu1, mu2 in pathological_states:
            y = rxn_system.y.copy()
            y[1], y[2], y[3] = mu0, mu1, mu2
            dn_dt = rxn_system.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt)), (
                f"Non-finite dn_dt for moment state (mu0,mu1,mu2)=({mu0},{mu1},{mu2}): {dn_dt}. "
                "This would trigger an unrecoverable model resurrection."
            )

    def test_mu3_closure_realizability_guard(self):
        """
        The log-Lagrange mu3 closure must not amplify an out-of-cone moment
        state. For a valid k>=1 distribution mu1 >= mu0 always; if mu1 < mu0 the
        closure returns 0.0 instead of (mu2/mu1)**3 * mu0, which would otherwise
        explode to a DASSL singularity (handoff BUG 1).
        """
        from rmgpy.solver.polymer import _safe_mu3_from_mu012
        # Realizable: returns the finite log-Lagrange value.
        assert _safe_mu3_from_mu012(1.0, 5.0, 30.0) > 0.0
        # Unrealizable (mu1 < mu0): guarded to 0.0 rather than ~4e14.
        assert _safe_mu3_from_mu012(21.17, 0.0139, 380.5) == 0.0

    def test_scission_moment_trajectory_stays_realizable(self):
        """
        Integration regression for the EPDM scission blow-up (handoff BUG 1).

        A scission-dominated pool must keep its moment trajectory inside the
        realizable cone (mu1 >= mu0 >= 0, all finite) for the whole run, and mass
        (mu1) must be conserved. Forward-Euler integration of the moment ODEs
        (no reactions) is enough to expose the original failure: with the buggy
        dmu0 = k_s*mu1 source, mu0 grows without bound and crosses mu1 within ~1
        time constant, after which the mu3 closure explodes (the DASSL IDID=-6
        singularity the EPDM deck hit). The fixed dmu0 = k_s*(mu1 - mu0) form is
        self-limiting and keeps mu0 <= mu1.
        """
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "Mu0")
        Mu1 = _spc("C=O", "Mu1")
        Mu2 = _spc("C#N", "Mu2")
        core_species = [Inert, Mu0, Mu1, Mu2]
        gas_species_mask = np.array([True, False, False, False], dtype=bool)

        k_s = 1.0
        V_poly = 1.0
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=k_s, k_unzip=0.0, tail_kinetics=None,
        )
        mu0_0, mu1_0, mu2_0 = 1.0, 5.0, 30.0
        rxn_system = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=V_poly,
            polymer_pools=[pool], mass_transfer=[], gas_species_mask=gas_species_mask,
            constant_gas_volume=False,
            initial_polymer_moments={"poly": (mu0_0, mu1_0, mu2_0)}, termination=[],
        )
        rxn_system.initialize_model(core_species, [], [], [])

        y = rxn_system.y.copy()
        dt = 2.0e-3
        n_steps = 5000  # t = 10, i.e. 10 scission time constants (1/k_s)
        i0, i1, i2 = 1, 2, 3
        for step in range(n_steps):
            dn_dt = rxn_system.residual(step * dt, y, np.zeros_like(y))[0]
            y = y + dt * dn_dt
            mu0, mu1, mu2 = y[i0] / V_poly, y[i1] / V_poly, y[i2] / V_poly
            assert np.all(np.isfinite(y)), f"non-finite state at step {step}"
            assert mu0 >= -1e-9, f"mu0 negative at step {step}: {mu0}"
            assert mu2 >= -1e-9, f"mu2 negative at step {step}: {mu2}"
            # The realizability invariant the bug violated:
            assert mu1 + 1e-9 >= mu0, f"mu1 < mu0 (unrealizable) at step {step}: mu0={mu0} mu1={mu1}"

        mu0_f, mu1_f, mu2_f = y[i0] / V_poly, y[i1] / V_poly, y[i2] / V_poly
        # Scission conserves mass (mu1).
        assert np.isclose(mu1_f, mu1_0, rtol=1e-6)
        # Ran toward completion: nearly all chains reduced to monomers (mu0 -> mu1).
        assert mu0_f > 0.9 * mu1_0
        # Distribution narrowed: mu2 fell from 30 toward mu1 (=5).
        assert mu2_f < mu2_0 and mu2_f < 10.0

    def test_debug_realizability_check_logs_not_raises(self, caplog):
        """
        The opt-in realizability diagnostic must LOG (never raise) when a pool's
        moment state leaves the realizable cone (here mu1 < mu0), and must stay
        silent when the flag is off. (handoff: a debug realizability assertion to
        localize moment-source bugs in one run.)
        """
        import logging

        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "Mu0")
        Mu1 = _spc("C=O", "Mu1")
        Mu2 = _spc("C#N", "Mu2")
        core_species = [Inert, Mu0, Mu1, Mu2]
        gas_species_mask = np.array([True, False, False, False], dtype=bool)

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )

        def build(flag):
            rs = HybridPolymerSystem(
                T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0,
                polymer_pools=[pool], mass_transfer=[], gas_species_mask=gas_species_mask,
                constant_gas_volume=False,
                # mu1 (0.1) < mu0 (5.0): unrealizable.
                initial_polymer_moments={"poly": (5.0, 0.1, 0.2)}, termination=[],
                debug_check_realizability=flag,
            )
            rs.initialize_model(core_species, [], [], [])
            return rs

        # Flag off: silent through init + residual.
        with caplog.at_level(logging.WARNING):
            rs_off = build(False)
            rs_off.residual(0.0, rs_off.y, np.zeros_like(rs_off.y))
        assert not any("realizable cone" in r.getMessage() for r in caplog.records)

        # Flag on: warns (exactly once, even across repeated residual calls) and
        # never raises. The first warning fires inside initialize_model's own
        # residual call, so capture across the whole build.
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            rs_on = build(True)
            out = rs_on.residual(0.0, rs_on.y, np.zeros_like(rs_on.y))
            rs_on.residual(0.0, rs_on.y, np.zeros_like(rs_on.y))  # still must not re-log
        assert out is not None
        warnings = [r for r in caplog.records if "realizable cone" in r.getMessage()]
        assert len(warnings) == 1, f"expected exactly one realizability warning, got {len(warnings)}"

    def test_flux_archetype_constants_match_enum(self):
        """The solver's mirror constants must equal PolymerFluxArchetype."""
        from rmgpy.polymer import PolymerFluxArchetype
        import rmgpy.solver.polymer as solver_mod
        assert solver_mod.FLUX_NONE == int(PolymerFluxArchetype.NONE) == 0
        assert solver_mod.FLUX_SAME_POOL == int(PolymerFluxArchetype.SAME_POOL) == 1
        assert solver_mod.FLUX_MIGRATION == int(PolymerFluxArchetype.MIGRATION) == 2
        assert solver_mod.FLUX_SCISSION_FRAGMENT == int(PolymerFluxArchetype.SCISSION_FRAGMENT) == 3
        assert solver_mod.FLUX_UNRESOLVED == int(PolymerFluxArchetype.UNRESOLVED) == 4
        assert solver_mod.FLUX_DISCRETE_CHIP == int(PolymerFluxArchetype.DISCRETE_CHIP) == 5

    def test_unstamped_proxy_reaction_remaps_to_unresolved(self):
        """
        A proxy-touching reaction arriving with the default archetype 0 (NONE)
        — e.g. restored from a pickle, which does not serialize the stamp —
        must be remapped to UNRESOLVED at initialize_model so the solver
        applies legacy mu1 flux instead of silently dropping pool moment flux.
        Pure-gas reactions must stay NONE.
        """
        Proxy = _spc("CCCC", "poly")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        A = _spc("C", "A")
        B = _spc("[CH3]", "B")

        core_species = [Proxy, Mu0, Mu1, Mu2, A, B]
        gas_species_mask = np.array([False, False, False, False, True, True], dtype=bool)

        kin = Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K"))
        proxy_rxn = Reaction(reactants=[Proxy], products=[B], kinetics=kin, reversible=False)
        gas_rxn = Reaction(reactants=[A], products=[B], kinetics=kin, reversible=False)
        assert proxy_rxn.polymer_flux_archetype == 0  # unstamped

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={A: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)}, termination=[],
        )
        rs.initialize_model(core_species, [proxy_rxn, gas_rxn], [], [])

        import rmgpy.solver.polymer as solver_mod
        assert rs.reaction_flux_archetype[0] == solver_mod.FLUX_UNRESOLVED
        assert rs.reaction_src_pool[0] == 0
        assert rs.reaction_dst_pool[0] == -1        # gas-only products
        assert rs.reaction_flux_archetype[1] == solver_mod.FLUX_NONE  # pure-gas stays NONE
        assert rs.reaction_src_pool[1] == -1
        assert rs.reaction_dst_pool[1] == -1

    def test_stamped_scission_without_daughter_pool_demotes_to_legacy(self):
        """
        A reaction stamped SCISSION_FRAGMENT whose polymer product is NOT a
        solver pool (e.g. a scission daughter registered as a species before
        its pool spawns) must demote to UNRESOLVED at initialize_model and
        apply the legacy mu1 drain -- otherwise the parent is never drained
        while the explicit daughter species still gains moles (duplication).
        """
        Proxy = _spc("CCCC", "poly")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        Daughter = _spc("CCC", "poly_scission_tail")   # core species, NO pool

        core_species = [Proxy, Mu0, Mu1, Mu2, Daughter]
        gas_species_mask = np.array([False] * 5, dtype=bool)

        rxn = Reaction(reactants=[Proxy], products=[Daughter], **_KIN)
        rxn.polymer_flux_archetype = 3   # stamped SCISSION_FRAGMENT

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        mu0, mu1, mu2 = 1.0, 5.0, 30.0
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (mu0, mu1, mu2)}, termination=[],
        )
        rs.initialize_model(core_species, [rxn], [], [])

        import rmgpy.solver.polymer as solver_mod
        assert rs.reaction_flux_archetype[0] == solver_mod.FLUX_UNRESOLVED  # demoted

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1                       # site-scaled by mu1, V_poly=1
        assert np.isclose(dn_dt[2], -r)    # parent mu1 drained (legacy)
        assert np.isclose(dn_dt[4], +r)    # explicit daughter gains
        # mass moved, not duplicated
        assert np.isclose(dn_dt[2] + dn_dt[4], 0.0, atol=1e-12)

    def test_stamped_chip_without_src_pool_demotes_to_unresolved(self):
        """
        Spec test 15: a reaction stamped DISCRETE_CHIP whose reactant does not
        resolve to a solver pool (src == -1) demotes to UNRESOLVED at
        initialize_model -- the same aggregate unresolvable-pool demotion path
        as MIGRATION/SCISSION_FRAGMENT (chip needs only src; there is no dst:
        the complement folds back and the chip is a gas species). In practice
        src == -1 arises when the chipped pool is a daughter Polymer species
        that is in the core but not yet in the solver's pool configs.
        """
        Proxy = _spc("CCCC", "poly")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        A = _spc("C", "A_gas")
        B = _spc("[CH3]", "B_gas")
        core = [Proxy, Mu0, Mu1, Mu2, A, B]
        mask = np.array([False] * 4 + [True, True], dtype=bool)

        rxn = Reaction(reactants=[A], products=[B], **_KIN)   # no pool reactant
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = 2

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={A: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)}, termination=[],
        )
        rs.initialize_model(core, [rxn], [], [])

        import rmgpy.solver.polymer as sp
        assert rs.reaction_src_pool[0] == -1
        assert rs.reaction_flux_archetype[0] == sp.FLUX_UNRESOLVED   # demoted
        assert rs.reaction_chip_units[0] == 2    # array filled regardless

    def test_validate_configuration_rejects_moment_in_stoichiometry(self):
        """
        validate_configuration should fail if a moment index appears in reaction stoichiometry.
        """
        A = _spc("C", "A_gas")
        Mu0 = _spc("CO", "Mu0")  # moment placeholder
        P2 = _spc("CC", "P2")    # polymer explicit placeholder

        core_species = [A, Mu0, P2]
        gas_species_mask = np.array([True, False, False], dtype=bool)

        pool = PolymerPoolConfig(
            label="poly",
            xs=2,
            explicit_dp_to_species_index={2: 2},
            mu_indices=(1, 1, 1),  # intentionally nonsense but won't be reached if stoich check triggers
            k_unzip=0.1,
        )

        # Put Mu0 in a reaction -> should be caught by moment isolation check
        rxn = Reaction(
            reactants=[Mu0],
            products=[P2],
            kinetics=Arrhenius(A=(1.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False,
        )

        rxn_system = HybridPolymerSystem(
            T=1000.0,
            P=1.0e5,
            initial_mole_fractions={A: 1.0},
            V_poly=1.0,
            polymer_pools=[pool],
            mass_transfer=[],
            gas_species_mask=gas_species_mask,
            constant_gas_volume=False,
            termination=[],
        )

        with pytest.raises(ValueError):
            rxn_system.initialize_model(core_species, [rxn], [], [])

    def test_initialize_model_accepts_two_pools(self):
        """Synthetic multi-pool: HybridPolymerSystem must accept and resolve
        two structurally-distinct PolymerPoolConfig objects in one solver.

        The single-pool tests above prove the per-pool plumbing works; this
        test covers the multi-pool path that the dynamic spawning feature
        (see ~/Code/RMG-Py/docs/multi_pool_design.md) targets: by the end of
        an RMG run with spawn-detection enabled, the solver may hold N>1
        pools. Validate that initialize_model resolves both pools cleanly
        and that pool_mu0_indices is populated for each.
        """
        Inert = _spc("N#N", "N2")

        # Pool A — three moment dummies (any non-isomorphic SMILES works).
        A_mu0 = _spc("CO", "PoolA_mu0")
        A_mu1 = _spc("C=O", "PoolA_mu1")
        A_mu2 = _spc("C#N", "PoolA_mu2")
        for spc in (A_mu0, A_mu1, A_mu2):
            spc.reactive = False

        # Pool B — three moment dummies with distinct structures.
        B_mu0 = _spc("CC", "PoolB_mu0")
        B_mu1 = _spc("C=C", "PoolB_mu1")
        B_mu2 = _spc("C#C", "PoolB_mu2")
        for spc in (B_mu0, B_mu1, B_mu2):
            spc.reactive = False

        core_species = [Inert, A_mu0, A_mu1, A_mu2, B_mu0, B_mu1, B_mu2]
        # Indices: Inert=0, PoolA mu0/1/2 = 1/2/3, PoolB mu0/1/2 = 4/5/6
        gas_mask = np.array(
            [True, False, False, False, False, False, False], dtype=bool
        )

        pool_a = PolymerPoolConfig(
            label="PoolA", xs=3,
            explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3),
        )
        pool_b = PolymerPoolConfig(
            label="PoolB", xs=3,
            explicit_dp_to_species_index={},
            mu_indices=(4, 5, 6),
        )

        rxn_system = HybridPolymerSystem(
            T=900.0, P=1.0e5,
            initial_mole_fractions={Inert: 1.0},
            V_poly=1.0e-3,
            polymer_pools=[pool_a, pool_b],
            gas_species_mask=gas_mask,
            initial_polymer_moments={
                "PoolA": (1.0e-3, 1.0e-2, 1.0e-1),
                "PoolB": (5.0e-4, 5.0e-3, 5.0e-2),
            },
        )

        rxn_system.initialize_model(core_species, [], [], [])

        assert len(rxn_system.polymer_pools) == 2
        assert rxn_system.polymer_pools[0].label == "PoolA"
        assert rxn_system.polymer_pools[1].label == "PoolB"
        # mu0 indices populated per-pool from the configured tuples.
        assert int(rxn_system.pool_mu0_indices[0]) == 1
        assert int(rxn_system.pool_mu0_indices[1]) == 4
        # Initial moments land at the correct state-vector slots.
        assert abs(rxn_system.y[1] - 1.0e-3) < 1e-12
        assert abs(rxn_system.y[4] - 5.0e-4) < 1e-12

    def test_spawn_rebuild_round_trip_preserves_moments(self):
        """
        Integration round-trip: a daughter pool produced by drain_spawn_intents
        (the Python iteration-boundary hook) must survive the solver rebuild. The
        rebuilt HybridPolymerSystem resolves BOTH the parent and the spawned
        daughter by label, carries the parent's moments across unchanged (moment
        continuity), and places the daughter's honest-empty moments ([0,0,0],
        item #14a) at its own state-vector slots.

        Closes the gap flagged in TestDrainSpawnIntents ("the Cython solver-reinit
        half follows separately and consumes the Polymer objects produced here")
        and the spawn->rebuild round-trip promised by docs/multi_pool_design.md.
        There is no in-place CVodeReInit: spawning works by registering the
        daughter's _muN dummy species and rebuilding the solver, which resolves
        moment indices by label.
        """
        from rmgpy.polymer import Polymer, SpawnIntent, drain_spawn_intents

        # --- Parent pool, pre-spawn solver (v1) -----------------------------
        parent = Polymer(label="PE", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"],
                         cutoff=3, Mn=1000.0, Mw=2500.0, initial_mass=1.0)
        parent.mu_indices = (1, 2, 3)

        Inert = _spc("N#N", "N2")
        P_mu0 = _spc("CO", "PE_mu0"); P_mu1 = _spc("C=O", "PE_mu1"); P_mu2 = _spc("C#N", "PE_mu2")
        for s in (P_mu0, P_mu1, P_mu2):
            s.reactive = False
        core_v1 = [Inert, P_mu0, P_mu1, P_mu2]
        gas_v1 = np.array([True, False, False, False], dtype=bool)

        parent_moments = (1.0e-3, 8.0e-3, 7.0e-2)
        rs1 = HybridPolymerSystem(
            T=900.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0e-3,
            polymer_pools=[PolymerPoolConfig(label="PE", xs=3,
                                             explicit_dp_to_species_index={}, mu_indices=(1, 2, 3))],
            gas_species_mask=gas_v1.copy(),
            initial_polymer_moments={"PE": parent_moments},
        )
        rs1.initialize_model(core_v1, [], [], [])
        # Parent state carried out of the pre-spawn solver.
        carried = (float(rs1.y[1]), float(rs1.y[2]), float(rs1.y[3]))
        assert abs(carried[0] - parent_moments[0]) < 1e-12

        # --- Spawn a daughter at the iteration boundary ---------------------
        # Honest-empty seeding (item #14a uniform-t=0): daughters start
        # mu = [0, 0, 0]; pools[].moments in the artifact are t=0 initial
        # conditions. triggering_dp survives as spawn METADATA only.
        DP = 5
        intent = SpawnIntent(parent_pool=parent, monomer=parent.backbone_group,
                             end_groups=["[H]", "[H]"], triggering_dp=DP)
        daughter = drain_spawn_intents([intent], iteration=1, existing_pools=[parent])[0]
        assert daughter.label == "PE_d1"
        assert np.allclose(daughter.moments, [0.0, 0.0, 0.0])
        assert daughter.spawn_metadata["triggering_dp"] == DP
        assert "triggering_moles" not in daughter.spawn_metadata

        # --- Rebuild solver (v2) with parent + daughter ---------------------
        # The daughter registers _muN dummy species labelled "<label>_muK";
        # the rebuilt solver resolves them by label.
        D_mu0 = _spc("CC", "PE_d1_mu0"); D_mu1 = _spc("C=C", "PE_d1_mu1"); D_mu2 = _spc("C#C", "PE_d1_mu2")
        for s in (D_mu0, D_mu1, D_mu2):
            s.reactive = False
        core_v2 = [Inert, P_mu0, P_mu1, P_mu2, D_mu0, D_mu1, D_mu2]
        gas_v2 = np.array([True, False, False, False, False, False, False], dtype=bool)

        rs2 = HybridPolymerSystem(
            T=900.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0e-3,
            polymer_pools=[
                PolymerPoolConfig(label="PE", xs=3, explicit_dp_to_species_index={}, mu_indices=(1, 2, 3)),
                PolymerPoolConfig(label="PE_d1", xs=3, explicit_dp_to_species_index={}, mu_indices=(4, 5, 6)),
            ],
            gas_species_mask=gas_v2.copy(),
            initial_polymer_moments={"PE": carried, "PE_d1": tuple(float(m) for m in daughter.moments)},
        )
        rs2.initialize_model(core_v2, [], [], [])

        # Both pools resolved by label.
        assert int(rs2.pool_mu0_indices[0]) == 1
        assert int(rs2.pool_mu0_indices[1]) == 4
        # Parent moments unchanged across the rebuild (continuity).
        assert abs(rs2.y[1] - carried[0]) < 1e-12
        assert abs(rs2.y[2] - carried[1]) < 1e-12
        assert abs(rs2.y[3] - carried[2]) < 1e-12
        # Daughter's honest-empty moments land at its own slots.
        assert rs2.y[4] == 0.0
        assert rs2.y[5] == 0.0
        assert rs2.y[6] == 0.0
        # The rebuilt RHS evaluates cleanly (finite, no NaN).
        dn = rs2.residual(0.0, rs2.y, np.zeros_like(rs2.y))[0]
        assert np.all(np.isfinite(dn))

    def test_initialization_subtracts_explicit_from_total_moments(self):
        """
        If initial_polymer_moments are provided (Total) along with explicit species,
        the solver should subtract the explicit contribution from the tail moments stored in y0.
        """
        # Define species
        Inert = _spc("N#N", "N2")
        P1 = _spc("C", "P1")  # DP=1
        Mu0 = _spc("CO", "Mu0")
        Mu1 = _spc("C=O", "Mu1")
        Mu2 = _spc("C#N", "Mu2")

        core_species = [Inert, P1, Mu0, Mu1, Mu2]
        # Indices: Inert=0, P1=1, Mu0=2, Mu1=3, Mu2=4
        gas_mask = np.array([True, False, False, False, False], dtype=bool)

        # Explicitly set P1 = 10.0 moles
        initial_explicit = {"poly": {1: 10.0}}

        # Provide TOTAL moments that *include* P1
        # Total Mu0 = 20.0 (10 chains of P1 + 10 chains of Tail)
        # Total Mu1 = 10*1 + 10*10 = 110.0 (assuming tail has DP=10)
        # Total Mu2 = 10*1^2 + 10*10^2 = 1010.0
        initial_total_moments = {"poly": (20.0, 110.0, 1010.0)}

        pool = PolymerPoolConfig(
            label="poly", xs=1,
            explicit_dp_to_species_index={1: 1},
            mu_indices=(2, 3, 4)
        )

        rxn_system = HybridPolymerSystem(
            T=300, P=1e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0,
            polymer_pools=[pool], gas_species_mask=gas_mask,
            initial_polymer_moments=initial_total_moments,
            initial_explicit_species=initial_explicit
        )

        # Initialize (triggers set_initial_conditions -> subtraction logic)
        rxn_system.initialize_model(core_species, [], [], [])

        # Check y0.
        # Expected Tail Mu0 = Total(20) - Explicit(10) = 10.0
        # Expected Tail Mu1 = Total(110) - Explicit(10*1) = 100.0
        y0_mu0 = rxn_system.y[2]
        y0_mu1 = rxn_system.y[3]

        assert abs(y0_mu0 - 10.0) < 1e-9
        assert abs(y0_mu1 - 100.0) < 1e-9

    def test_edge_reaction_fluxes_are_diagnostic_only(self):
        """
        Edge reactions (all-core reactants, edge product) must not perturb the
        integrated core state: dn_dt and consumption/production stay zero, while
        edge_reaction_rates / edge_species_rates carry the diagnostic flux.
        Matches simple.pyx semantics (core ODE = core reactions only).
        """
        A = _spc("C", "A")        # core gas reactant
        E = _spc("[CH3]", "E")    # edge product

        kin = Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K"))
        edge_rxn = Reaction(reactants=[A], products=[E], kinetics=kin, reversible=False)

        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={A: 1.0}, V_poly=1.0,
            polymer_pools=[], mass_transfer=[],
            gas_species_mask=np.array([True], dtype=bool), constant_gas_volume=False,
            termination=[],
            # Item 17 A5-2: direct (no-blueprint-phase) build with a non-empty
            # edge -- a legitimate last-resort fallback; flag it so R1-EDGE does
            # not raise on the default-filled edge suffix.
            allow_default_prospective_edge=True,
        )
        rs.initialize_model([A], [], [E], [edge_rxn])

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        assert dn_dt[0] == 0.0                          # A untouched by edge rxn
        assert rs.core_species_consumption_rates[0] == 0.0
        assert rs.core_species_production_rates[0] == 0.0
        assert rs.edge_reaction_rates[0] > 0.0          # diagnostics still flow
        assert rs.edge_species_rates[0] > 0.0

    def test_migration_moves_whole_chain_bundle(self):
        """
        MIGRATION (archetype 2): one event moves a whole length-biased chain
        from pool A to pool B: bundle (1, mu2/mu1, mu3/mu1) with the
        log-Lagrange closure mu3 = mu0*(mu2/mu1)**3. A loses exactly what B
        gains (conservation by construction).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1a                       # site-scaled by A's mu1, V_poly=1
        mu3a = mu0a * (mu2a / mu1a) ** 3    # 216.0
        b1 = mu2a / mu1a                    # 6.0
        b2 = mu3a / mu1a                    # 43.2

        assert np.isclose(dn_dt[1], -r * 1.0)        # A mu0
        assert np.isclose(dn_dt[2], -r * b1)         # A mu1
        assert np.isclose(dn_dt[3], -r * b2)         # A mu2
        assert np.isclose(dn_dt[5], +r * 1.0)        # B mu0
        assert np.isclose(dn_dt[6], +r * b1)         # B mu1
        assert np.isclose(dn_dt[7], +r * b2)         # B mu2
        assert np.isclose(dn_dt[1] + dn_dt[5], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[2] + dn_dt[6], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[3] + dn_dt[7], 0.0, atol=1e-14)

    def test_migration_reverse_leg_uses_target_pool_statistics(self):
        """
        Per-direction MIGRATION bundles: the reverse (rr) leg must move
        B-statistics chains B->A, not A-statistics. Pool stats are chosen
        distinguishable (b_A=(1,6,43.2) vs b_B=(1,2.5,7.8125)) so a
        wrong-source bundle fails loudly. Also pins the rf/rr volume-factor
        identity: at b0=1 on both legs, the net mu0 flux must equal the
        legacy net molar rate (rf-rr)*V_rxn.
        """
        sp, core, mask = _two_pool_species()
        # reversible=False so generate_rate_coefficients needs no thermo
        # (kb=0); the reverse leg is then driven by overriding rs.kb directly.
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        mu_a = (1.0, 5.0, 30.0)
        mu_b = (2.0, 4.0, 10.0)
        rs = _two_pool_rs(rxn, core, mask, mu_a, mu_b)
        rs.kb[0] = 0.6

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        rf = kf * mu_a[1]               # 2.0 * 5 = 10 (site-scaled by A mu1)
        rr = 0.6 * mu_a[1]              # kb * C(proxyB)=1, then *= site -> 3.0
        mu3_a = mu_a[0] * (mu_a[2] / mu_a[1]) ** 3   # 216.0
        mu3_b = mu_b[0] * (mu_b[2] / mu_b[1]) ** 3   # 31.25
        bA1, bA2 = mu_a[2] / mu_a[1], mu3_a / mu_a[1]    # 6.0, 43.2
        bB1, bB2 = mu_b[2] / mu_b[1], mu3_b / mu_b[1]    # 2.5, 7.8125

        assert np.isclose(dn_dt[1], -(rf - rr))              # == legacy net (b0=1)
        assert np.isclose(dn_dt[2], -rf * bA1 + rr * bB1)    # -52.5
        assert np.isclose(dn_dt[3], -rf * bA2 + rr * bB2)    # -408.5625
        assert np.isclose(dn_dt[5], +(rf - rr))
        assert np.isclose(dn_dt[6], +rf * bA1 - rr * bB1)
        assert np.isclose(dn_dt[7], +rf * bA2 - rr * bB2)

    def test_end_group_migration_uses_uniform_chain_bundle(self):
        """
        A mu0-scaled (is_end_group_reaction) MIGRATION picks chains uniformly:
        bundle (1, mu1/mu0, mu2/mu0), and the rate itself scales by mu0.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        rxn.is_end_group_reaction = True
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu0a                    # site-scaled by mu0, not mu1
        assert np.isclose(dn_dt[1], -r * 1.0)
        assert np.isclose(dn_dt[2], -r * mu1a / mu0a)        # -10
        assert np.isclose(dn_dt[3], -r * mu2a / mu0a)        # -60
        assert np.isclose(dn_dt[5], +r * 1.0)
        assert np.isclose(dn_dt[6], +r * mu1a / mu0a)
        assert np.isclose(dn_dt[7], +r * mu2a / mu0a)

    def test_scission_fragment_complement_stays_in_parent(self):
        """
        SCISSION_FRAGMENT (archetype 3), complement-stays accounting: parent
        (0, -r*mu2/(2 mu1), -r*(2/3)*mu3/mu1); daughter (+r, +r*mu2/(2 mu1),
        +r*mu3/(3 mu1)). mu1 conserves exactly; total mu0 +r per event; total
        mu2 drops by r*mu3/(3 mu1).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1a
        mu3a = mu0a * (mu2a / mu1a) ** 3    # 216.0
        e_n = mu2a / mu1a                   # 6.0
        e_n2 = mu3a / mu1a                  # 43.2

        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)             # parent mu0
        assert np.isclose(dn_dt[2], -r * e_n / 2.0)              # parent mu1
        assert np.isclose(dn_dt[3], -r * (2.0 / 3.0) * e_n2)     # parent mu2
        assert np.isclose(dn_dt[5], +r)                          # daughter mu0
        assert np.isclose(dn_dt[6], +r * e_n / 2.0)              # daughter mu1
        assert np.isclose(dn_dt[7], +r * e_n2 / 3.0)             # daughter mu2
        assert np.isclose(dn_dt[2] + dn_dt[6], 0.0, atol=1e-14)  # mu1 conserved
        assert np.isclose(dn_dt[3] + dn_dt[7], -r * e_n2 / 3.0)  # mu2 destroyed
        assert np.isclose(dn_dt[8], +r)                          # gas co-product G

    def test_same_pool_reaction_leaves_moments_untouched(self):
        """
        SAME_POOL (archetype 1): fold-back means net-zero pool moment flux —
        the dispatch skips pool writes entirely. Gas co-species still flow.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        assert np.allclose(dn_dt[1:4], 0.0, atol=1e-14)   # A moments untouched
        assert np.allclose(dn_dt[5:8], 0.0, atol=1e-14)   # B untouched
        assert np.isclose(dn_dt[8], kf * 5.0)             # gas G produced

    def test_unresolved_applies_legacy_mu1_flux(self):
        """UNRESOLVED (archetype 4) replicates the pre-apportionment behavior:
        whole event flux on mu1 only (reactant pool -r, product pool +r)."""
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 4
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * 5.0
        assert np.isclose(dn_dt[2], -r)                   # A mu1 only
        assert np.isclose(dn_dt[6], +r)                   # B mu1 only
        assert np.allclose([dn_dt[1], dn_dt[3], dn_dt[5], dn_dt[7]], 0.0, atol=1e-14)

    def test_bundle_guard_empty_source_pool_no_nan(self):
        """
        An empty source pool (all moments zero) must produce zero flux and no
        NaN: the site scaling already zeroes the rate, and the bundle guard
        prevents the 0/0 in mu2/mu1.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        rs = _two_pool_rs(rxn, core, mask, (0.0, 0.0, 0.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(np.isfinite(dn_dt))
        assert np.allclose(dn_dt[1:8], 0.0, atol=1e-14)

    def test_bundle_guard_mu3_overflow_skips_mu2_component_only(self):
        """
        When the mu3 closure overflows to inf (huge mu2), MIGRATION still
        applies the mu0/mu1 bundle components and skips ONLY mu2 (mirrors the
        scission-ODE finiteness gate).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        mu0a, mu1a, mu2a = 1.0, 5.0, 1.0e120
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1a
        assert np.all(np.isfinite(dn_dt))
        assert np.isclose(dn_dt[1], -r)                       # mu0 applied
        assert np.isclose(dn_dt[2], -r * mu2a / mu1a)         # mu1 applied
        assert dn_dt[3] == 0.0                                # mu2 skipped
        assert dn_dt[7] == 0.0

    def test_discrete_chip_monodisperse_closed_form_both_picks(self):
        """
        Spec test 10: monodisperse pool (mu_j = N*L^j) -> E[n] = L under BOTH
        picks (uniform mu1/mu0 == length-biased mu2/mu1 == L on monodisperse),
        so the chip drain has the same closed form for either flag value:
        dmu0 = 0, dmu1 = -r*a, dmu2 = -r*(2aL - a^2). Only the rate's site
        scaling differs (mu0 vs mu1). The chip species flows through the
        standard gas path (+r).
        """
        N, L, a = 2.0, 10.0, 3
        for end_group in (False, True):
            sp, core, mask = _two_pool_species()
            rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
            rxn.polymer_flux_archetype = 5
            rxn.polymer_chip_units = a
            rxn.is_end_group_reaction = end_group
            rs = _two_pool_rs(rxn, core, mask, (N, N * L, N * L * L), (0.1, 0.2, 0.5))

            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

            kf = rxn.get_rate_coefficient(800.0, 1.0e5)
            r = kf * (N if end_group else N * L)          # site scaling per flag
            assert np.isclose(dn_dt[1], 0.0, atol=1e-14), end_group   # mu0
            assert np.isclose(dn_dt[2], -r * a), end_group            # mu1
            assert np.isclose(dn_dt[3], -r * (2.0 * a * L - a * a)), end_group
            assert np.allclose(dn_dt[5:8], 0.0, atol=1e-14)           # pool B idle
            assert np.isclose(dn_dt[8], +r)                           # chip gas

    def test_discrete_chip_uses_scaling_consistent_e_n(self):
        """
        Spec test 11 (decision D5): polydisperse pool (1, 5, 30) separates the
        picks -- uniform E[n] = mu1/mu0 = 5, length-biased E[n] = mu2/mu1 = 6.
        The mu2 drain must use the E[n] matching the reaction's rate-scaling
        flag; pairing a mu0 rate with length-biased E[n] (or vice versa)
        fails this test loudly.
        """
        a = 1
        for end_group, e_n, site in ((False, 6.0, 5.0), (True, 5.0, 1.0)):
            sp, core, mask = _two_pool_species()
            rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
            rxn.polymer_flux_archetype = 5
            rxn.polymer_chip_units = a
            rxn.is_end_group_reaction = end_group
            rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.1, 0.2, 0.5))

            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

            kf = rxn.get_rate_coefficient(800.0, 1.0e5)
            r = kf * site
            assert np.isclose(dn_dt[2], -r * a), end_group
            assert np.isclose(dn_dt[3], -r * (2.0 * a * e_n - a * a)), end_group

    def test_discrete_chip_clamp_regime_skips_mu2_write(self):
        """
        Spec test 12: a >= 2*E[n] makes 2aE[n] - a^2 <= 0 -- impossible
        per-chain (n >= a) but reachable in expectation for a pool whose mean
        length decayed toward chip size. The forward mu2 decrement is clamped
        (write skipped), mu1 still drains, RHS finite.
        """
        a = 13   # length-biased E[n] = 6 -> 2*13*6 - 169 = -13 < 0
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = a
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * 5.0                                  # mu1-scaled
        assert np.all(np.isfinite(dn_dt))
        assert np.isclose(dn_dt[2], -r * a)           # mu1 still drains
        assert dn_dt[3] == 0.0                        # mu2 write skipped
        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)  # mu0 untouched

    def test_discrete_chip_reverse_leg_exact_extension_form(self):
        """
        Spec test 14: the reverse leg is the EXACT extension form
        Delta(n^2) = (n+a)^2 - n^2 = +(2aE[n] + a^2): dmu1 = +rr*a,
        dmu2 = +rr*(2aE[n] + a^2) -- PLUS a^2, not the forward sign-flip
        (which would subtract it) -- and never clamps, even at a >= 2*E[n]
        where the forward leg would. Driven via the rs.kb override with
        kf = 0 and injected chip moles. The site magnitude is the spec 5
        exhaustion-throttle value min(mu0, mu1/a) (amendment 2026-06-10),
        derived from the fixture moments -- so this test also pins that the
        throttle multiplies the reverse leg.
        """
        a = 13                                # forward would clamp here
        mu_a = (1.0, 5.0, 30.0)
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = a
        rxn.is_end_group_reaction = True      # uniform pick: E[n] = mu1/mu0 = 5
        rs = _two_pool_rs(rxn, core, mask, mu_a, (0.1, 0.2, 0.5))
        rs.kf[0] = 0.0                        # silence the forward leg
        rs.kb[0] = 0.4                        # drive the reverse leg directly

        y = rs.y.copy()
        y[8] = 1.0                            # inject chip (G) moles

        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]

        # G is the only gas species with moles -> C_G = P/(R*T) exactly.
        c_g = 1.0e5 / (constants.R * 800.0)
        # Exhaustion throttle (spec 5 amendment 2026-06-10): the mu0-scaled
        # chip site is min(mu0, mu1/a), and it multiplies BOTH directions --
        # here mu1/a = 5/13 < mu0 = 1, so the reverse leg is throttled too.
        site = min(mu_a[0], mu_a[1] / a)
        rr = 0.4 * 1.0 * c_g * site           # kb * C(fold-back proxy)=1 * C(G) * site
        e_n = 5.0
        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)              # mu0 untouched
        assert np.isclose(dn_dt[2], +rr * a)                      # exact form
        assert np.isclose(dn_dt[3], +rr * (2.0 * a * e_n + a * a))  # +a^2, no clamp
        assert np.isclose(dn_dt[8], -rr)      # chip consumed via standard path

    def test_discrete_chip_exhaustion_trajectory_throttles(self):
        """
        Spec test 12b (exhaustion-throttle amendment 2026-06-10): a mu0-scaled
        chip rate is constant in mu1 unthrottled (rate ~ mu0, which chip
        events never drain), so past unit exhaustion mu1 would run linearly
        negative while chip moles keep being created. With the throttle,
        mu1 decays at worst exponentially and never crosses zero, chip
        production slows in lockstep, the FULL cone (mu1 >= 0 AND
        mu0*mu2 >= mu1^2) holds throughout, and Sum(mu1) + a*n_chip is
        conserved at every step.
        """
        a = 3
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = a
        rxn.is_end_group_reaction = True
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.0, 0.0, 0.0))

        y = rs.y.copy()
        invariant0 = y[2] + a * y[8]          # Sum(mu1) + a*n_chip
        # kf = 2/s, so UNTHROTTLED dynamics drain mu1 at kf*mu0*a = 6/s
        # CONSTANT, crossing zero at t ~ 0.83 s and reaching -19 by t = 4 s.
        # Forward Euler well past exhaustion: dt*kf = 0.01 (stable).
        dt = 0.005
        chip_increments = []
        for _ in range(800):                  # t = 4 s
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
            chip_increments.append(dt * dn_dt[8])
            assert y[2] >= 0.0                              # mu1 never crosses zero
            if y[2] > 1e-12:
                assert y[1] * y[3] >= y[2] ** 2 * (1.0 - 1e-9)   # full cone
            assert dn_dt[8] >= 0.0                          # chip moles nondecreasing
        # Production decays in lockstep with mu1 (exponential, not constant).
        assert chip_increments[-1] < 1e-2 * chip_increments[0]
        assert np.isclose(y[2] + a * y[8], invariant0, rtol=1e-9, atol=1e-12)

    def test_discrete_chip_throttle_direct_rhs_regimes(self):
        """
        Spec 5 amendment, direct-RHS pins. Throttled regime (mu1 < a*mu0):
        site = mu1/a, so dmu1/dt = -kf*(mu1/a)*a = -kf*mu1 EXACTLY -- not
        -kf*mu0*a, which would be 3x larger here. Healthy regime
        (mu1 >> a*mu0): site = mu0, byte-identical to the pre-throttle path.
        """
        a = 3
        kf = None
        for mom_a, site in (((1.0, 1.0, 2.0), 1.0 / 3.0),    # throttled: min(1, 1/3)
                            ((1.0, 50.0, 3000.0), 1.0)):     # healthy:  min(1, 50/3)
            sp, core, mask = _two_pool_species()
            rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
            rxn.polymer_flux_archetype = 5
            rxn.polymer_chip_units = a
            rxn.is_end_group_reaction = True
            rs = _two_pool_rs(rxn, core, mask, mom_a, (0.1, 0.2, 0.5))

            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            kf = rxn.get_rate_coefficient(800.0, 1.0e5)
            assert np.isclose(dn_dt[2], -kf * site * a), mom_a   # mu1 drain
            assert np.isclose(dn_dt[8], +kf * site), mom_a       # chip production
            if site < mom_a[0]:
                # Throttled identity: dmu1/dt = -kf*mu1 exactly.
                assert np.isclose(dn_dt[2], -kf * mom_a[1])
                assert not np.isclose(dn_dt[2], -kf * mom_a[0] * a)

    def test_discrete_chip_zero_unit_chip_exempt_from_throttle(self):
        """
        Spec 5 amendment: a = 0 chips drain nothing, so they are exempt from
        the throttle (no mu1/a division, no rate reduction): the chip species
        rate stays kf*mu0 even when mu1 is small relative to mu0.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = 0
        rxn.is_end_group_reaction = True
        rs = _two_pool_rs(rxn, core, mask, (2.0, 1.0, 1.0), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        assert np.all(np.isfinite(dn_dt))
        assert np.isclose(dn_dt[8], +kf * 2.0)        # chip rate kf*mu0, unthrottled
        assert np.isclose(dn_dt[2], 0.0, atol=1e-14)  # a = 0 drains no units
        assert np.isclose(dn_dt[3], 0.0, atol=1e-14)  # and no mu2

    def test_discrete_chip_throttle_diagnostic_rate_parity(self):
        """
        Spec 5 amendment, diagnostic-path parity: get_reaction_rates (the
        [THE HIJACK] block feeding get_edge_reaction_rates) must apply the
        SAME exhaustion throttle as the residual's site scaling. mu0 = 1,
        mu1 = 1, a = 3 is throttled (mu1/a = 1/3 < mu0), so the diagnostic
        rate is kf*min(mu0, mu1/a) = kf/3 -- NOT the unthrottled kf*mu0 --
        and equals what the residual's dmu1/(-a) implies.
        """
        a = 3
        mom_a = (1.0, 1.0, 2.0)               # throttled: min(1, 1/3) = 1/3
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = a
        rxn.is_end_group_reaction = True
        rs = _two_pool_rs(rxn, core, mask, mom_a, (0.1, 0.2, 0.5))

        rate = rs.get_reaction_rates(rs.y)[0]
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        site = min(mom_a[0], mom_a[1] / a)
        assert np.isclose(rate, kf * site)                # throttled, not kf*mu0
        assert not np.isclose(rate, kf * mom_a[0])
        assert np.isclose(rate, dn_dt[2] / (-a))          # parity with residual

    def test_discrete_chip_trajectory_conserves_units_and_cone(self):
        """
        Spec test 13 -- the exact conservation invariant (closed system,
        chip reactions only): d/dt [ Sigma_pools mu1 + Sigma_chips a_i*n_i ]
        = 0, chip moles weighted by the STAMPED unit count (not raw moles).
        Multi-chip generalization of the single-reaction invariant already
        pinned by test_discrete_chip_exhaustion_trajectory_throttles: TWO
        chip reactions with DIFFERENT unit counts on DIFFERENT pools (a = 3,
        mu0-scaled on A; a = 5, mu1-scaled on B) ejecting DISTINCT chip
        species -- the invariant only closes if each chip is weighted by its
        own a_i. Over a chip-only forward-Euler trajectory: the weighted
        invariant is constant to roundoff, mu0 never changes for either pool
        (chain counts preserved), and both pools stay in the realizable cone
        (mu1 >= mu0 >= 0).
        """
        sp, core, mask = _two_pool_species()
        g2 = _spc("[CH2]C", "G2")
        sp["G2"] = g2
        core = core + [g2]
        mask = np.append(mask, True)
        a1, a2 = 3, 5
        rxn1 = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn1.polymer_flux_archetype = 5
        rxn1.polymer_chip_units = a1
        rxn1.is_end_group_reaction = True     # mu0-scaled, uniform pick
        rxn2 = Reaction(reactants=[sp["B"]], products=[sp["B"], sp["G2"]], **_KIN)
        rxn2.polymer_flux_archetype = 5
        rxn2.polymer_chip_units = a2
        rxn2.is_end_group_reaction = False    # mu1-scaled, length-biased pick

        pool_a = PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={},
                                   mu_indices=(1, 2, 3), monomer_poly_index=None,
                                   k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
        pool_b = PolymerPoolConfig(label="B", xs=2, explicit_dp_to_species_index={},
                                   mu_indices=(5, 6, 7), monomer_poly_index=None,
                                   k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={sp["G"]: 0.0, sp["G2"]: 0.0},
            V_poly=1.0, polymer_pools=[pool_a, pool_b], mass_transfer=[],
            gas_species_mask=mask, constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 50.0, 3000.0),
                                     "B": (1.0, 50.0, 3000.0)},
            termination=[],
        )
        rs.initialize_model(core, [rxn1, rxn2], [], [])

        y = rs.y.copy()
        invariant0 = y[2] + y[6] + a1 * y[8] + a2 * y[9]   # Sum(mu1) + Sum(a_i*n_i)
        # rxn1: r = kf*mu0 = 2/s -> dmu1_A = -6/s; rxn2: dmu1_B/dt =
        # -kf*a2*mu1_B = -10*mu1_B/s (self-limiting). 200 * 1e-4 s drains
        # ~0.12 of A and ~18% of B: far from depletion, forward Euler safe.
        dt = 1e-4
        for _ in range(200):
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
            for mu0_i, mu1_i in ((1, 2), (5, 6)):
                assert y[mu0_i] >= -1e-12              # mu0 >= 0
                assert y[mu1_i] - y[mu0_i] >= -1e-9    # cone: mu1 >= mu0
        assert np.isclose(y[2] + y[6] + a1 * y[8] + a2 * y[9], invariant0,
                          rtol=1e-12)
        # DISCRETE_CHIP writes no mu0 leg, so chain counts are bit-exact.
        assert np.isclose(y[1], 1.0, rtol=0, atol=1e-14)
        assert np.isclose(y[5], 1.0, rtol=0, atol=1e-14)
        assert y[8] > 0.0 and y[9] > 0.0      # both chips accumulated

    def test_scission_monodisperse_limit_closed_form(self):
        """
        Monodisperse pool (PDI=1, mu_j = N*k^j): the mu3 closure is exact
        (mu3 = N*k^3), so SCISSION_FRAGMENT has closed-form derivatives:
        parent (0, -r*k/2, -r*(2/3)*k^2), daughter (+r, +r*k/2, +r*k^2/3).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        N, k = 2.0, 5.0
        rs = _two_pool_rs(rxn, core, mask, (N, N * k, N * k * k), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * (N * k)                                  # site-scaled by mu1

        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[2], -r * k / 2.0)
        assert np.isclose(dn_dt[3], -r * (2.0 / 3.0) * k * k)
        assert np.isclose(dn_dt[5], +r)
        assert np.isclose(dn_dt[6], +r * k / 2.0)
        assert np.isclose(dn_dt[7], +r * k * k / 3.0)

    def test_apportionment_trajectory_conserves_mu1_and_stays_realizable(self):
        """
        Forward-Euler trajectory with a SCISSION_FRAGMENT reaction A -> B + G:
        total monomer units (mu1_A + mu1_B) stay constant (the gas co-product
        G tracks events, not units), total chain count mu0_A + mu0_B grows,
        and both pools stay in the realizable cone (mu1 >= mu0 >= 0).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        rs = _two_pool_rs(rxn, core, mask, (1.0, 50.0, 3000.0), (0.0, 0.0, 0.0))

        y = rs.y.copy()
        mu1_total0 = y[2] + y[6]
        mu0_total0 = y[1] + y[5]
        # Step sizing: parent mu1 drains at ~kf*mu2/2 = 3000/s initially, so
        # keep t_total = 2e-3 s (~12% of parent mu1) to stay far from the
        # depletion overshoot regime that forward Euler handles poorly.
        dt = 1e-5
        for _ in range(200):
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
            for i0, i1 in ((1, 2), (5, 6)):
                assert y[i1] >= -1e-9                  # mu1 >= 0
                assert y[i1] - y[i0] >= -1e-6          # mu1 >= mu0 (cone)
        assert np.isclose(y[2] + y[6], mu1_total0, rtol=1e-9)   # units conserved
        assert y[1] + y[5] > mu0_total0                          # chains created

    def test_scission_mu2_overdrain_stays_finite_and_warns(self, caplog):
        """
        A high-PDI parent makes the closure-estimated mu3 huge; the resulting
        parent mu2 drain can overshoot in an explicit step. The residual must
        stay finite when evaluated at the resulting negative-mu2 state (reads
        are max(0,.)-clamped), and debug_check_realizability must log the cone
        violation rather than raise.
        """
        import logging as _logging
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 1000.0), (0.1, 0.2, 0.5))
        rs.debug_check_realizability = True

        y = rs.y.copy()
        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
        assert np.all(np.isfinite(dn_dt))
        # Overshoot mu2 negative with one huge explicit step, then re-evaluate.
        y2 = y + 1e3 * dn_dt
        assert y2[3] < 0.0, "test setup: expected a mu2 overshoot"
        with caplog.at_level(_logging.WARNING):
            dn_dt2 = rs.residual(0.0, y2, np.zeros_like(y2))[0]
        assert np.all(np.isfinite(dn_dt2))
        assert any("realizable cone" in rec.message for rec in caplog.records)

    def test_scission_net_reverse_guard_on_empty_daughter(self):
        """
        SCISSION_FRAGMENT with net-reverse flux (r < 0) depletes the DAUGHTER;
        with an empty daughter the dispatch must skip entirely (no negative
        moments manufactured) and stay finite.

        Net-reverse is achieved by zeroing kf[0] and setting kb[0] while
        injecting G moles (y[8]=1.0) so that rr = kb*_C(B)*_C(G)*site > 0.
        With no G present the gas _C(G)=0 and rr=0 — the G injection is the
        minimal change needed to reach the r<0 branch. Note _KIN builds the
        reaction with reversible=False (so initialization leaves kb=0); the
        kb array is patched directly on the solver to bypass that without
        needing thermo for Keq.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.0, 0.0, 0.0))
        # Force net-reverse: zero forward (kf) and inject a reverse rate (kb).
        rs.kf[0] = 0.0
        rs.kb[0] = 0.6   # arbitrary nonzero -- only the sign of the net rate matters
        # Inject G moles so _C(G) = n_G/V_gas > 0 and rr > 0.
        y = rs.y.copy()
        y[8] = 1.0
        # Verify net-reverse before the main assertion.
        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
        assert np.all(np.isfinite(dn_dt))
        assert np.allclose(dn_dt[1:8], 0.0, atol=1e-14)   # all pool moments untouched

    def test_scission_mu3_overflow_skips_mu2_components(self):
        """
        The SCISSION_FRAGMENT branch independently guards its mu2 components
        on mu3 finiteness (the migration branch has its own gate inside
        _chain_bundle): with a mu2 huge enough to overflow the closure, the
        mu0/mu1 scission terms still apply and BOTH mu2 components are
        skipped.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        mu0a, mu1a, mu2a = 1.0, 5.0, 1.0e120
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1a
        e_n = mu2a / mu1a
        assert np.all(np.isfinite(dn_dt))
        assert np.isclose(dn_dt[2], -r * e_n / 2.0)   # parent mu1 applied
        assert np.isclose(dn_dt[5], +r)               # daughter mu0 applied
        assert np.isclose(dn_dt[6], +r * e_n / 2.0)   # daughter mu1 applied
        assert dn_dt[3] == 0.0                        # parent mu2 skipped
        assert dn_dt[7] == 0.0                        # daughter mu2 skipped


class TestIntegratedSpawnGateTripwire:
    """Spec §7.1 — INTEGRATED live-path tripwire, two halves, RED-FIRST.

    A real HybridPolymerSystem integrated solve (the apportionment-trajectory
    forward-Euler idiom), NOT a fabricated snapshot. The SAME_POOL reaction
    A -> A + R is apportioned (non-UNRESOLVED) pool-touching chemistry whose
    product R is an ordinary (non-canonical-proxy) core species. The
    numerator half MUST fail on pre-change-(a) HEAD: gross arrays are
    maintained only in is_pool_proxy branches today (spec §3.1) — that red
    run is the proof the fix is real and the gate to executing the rest of
    the plan. The polymer.pyx line citations in spec §3 are informational;
    these tests assert the MECHANISM, not the address.
    """

    @staticmethod
    def _solve():
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1  # SAME_POOL: apportioned, non-UNRESOLVED
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0), monomer_mw_g_mol=28.0)
        # Short integrated solve (the trajectory-test idiom): evolve the
        # engine state, then evaluate the residual AT the evolved state so
        # the gross arrays hold the last evaluation — exactly what
        # spawn_gate_flux_snapshot() reads after a production simulate().
        y = rs.y.copy()
        dt = 1e-4
        for _ in range(50):
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
        rs.y[:] = y
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        return sp, rs, dn_dt

    def test_numerator_half_ordinary_species_gross_is_real(self):
        """Numerator half (the regression that would have caught the
        born-dead class): the ordinary product R has a NONZERO gross entry
        in the snapshot, equal to max(0, core_species_production_rates[R])
        recomputed independently from the engine arrays, and
        g_R = gross * E[n] * monomer_MW under its parent pool's pool_stats.
        """
        sp, rs, dn_dt = self._solve()
        i_r = 4
        # LIVENESS PIN — must come BEFORE the red assertion. R must be
        # chemically ALIVE (nonzero net production) so the red below can
        # only mean "alive but no gross record" — never "fixture dead".
        # A dead fixture (archetype eats the reaction, phase gate, wrong
        # index) would zero prod_r too and launder the wrong failure as
        # the born-dead proof. xfail(strict) cannot tell those apart;
        # this assertion can.
        assert float(dn_dt[i_r]) > 0.0, (
            "FIXTURE BROKEN, not a valid red: representative R has zero "
            "net production - fix the fixture before trusting the red"
        )
        prod_r = float(rs.core_species_production_rates[i_r])
        # THE red assertion: fails on pre-change-(a) HEAD because ordinary
        # core species get dn_dt writes only — no gross record exists.
        assert prod_r > 0.0, (
            "ordinary core species R has no gross production record: the "
            "gross arrays are proxy-only (spec §3.1 born-dead hole; "
            "change (a) pending)"
        )
        # Independent recompute: irreversible reaction, V_poly = V_rxn = 1
        # -> R's gross production must equal its net dn_dt exactly.
        assert prod_r == pytest.approx(float(dn_dt[i_r]), rel=1e-12)

        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()
        assert gross["R"] == pytest.approx(max(0.0, prod_r), rel=1e-12)
        e_n, mw = pool_stats["A"]
        assert e_n == pytest.approx(float(rs.y[2]) / float(rs.y[1]), rel=1e-12)
        assert mw == pytest.approx(28.0)
        g_r = gross["R"] * e_n * mw
        assert g_r == pytest.approx(
            max(0.0, prod_r) * (float(rs.y[2]) / float(rs.y[1])) * 28.0, rel=1e-12)
        assert g_r > 0.0

    def test_denominator_half_proxy_net_rerouted_gross_nonzero(self):
        """Denominator half: the canonical proxy's net dn_dt contribution is
        ~= 0 (the apportionment reroutes proxy flux to pool moments) while
        its gross entry is nonzero — an assertion that is only true of the
        GROSS array and dies if the denominator path is ever rewired to net
        rates."""
        sp, rs, dn_dt = self._solve()
        assert dn_dt[0] == pytest.approx(0.0, abs=1e-14)
        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()
        assert gross["A"] > 0.0
        e_n, mw = pool_stats["A"]
        assert proxy_total == pytest.approx(gross["A"] * e_n * mw, rel=1e-12)
        assert proxy_total > 0.0


class TestSpawnGateFluxSnapshot:
    """spawn_gate_flux_snapshot() unit pins (spec §4.1): 3-tuple shape,
    all-core gross coverage, engine-attributed proxy event-mass total, and
    the E[n]*MW calibration (spec §7.6 / decision 3)."""

    def test_snapshot_three_tuple_covers_all_core_species(self):
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0), monomer_mw_g_mol=28.0)

        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()

        # gross: EVERY core species by label, max(0, production) — moment
        # dummies and the untouched gas species carry explicit zeros.
        assert set(gross.keys()) == {"A", "A_mu0", "A_mu1", "A_mu2", "R", "G"}
        assert gross["A_mu0"] == 0.0 and gross["A_mu1"] == 0.0 and gross["A_mu2"] == 0.0
        assert gross["G"] == 0.0
        assert gross["R"] == pytest.approx(
            max(0.0, float(rs.core_species_production_rates[4])), rel=1e-12)
        assert gross["R"] > 0.0
        assert gross["A"] > 0.0
        # pool_stats: pool label -> (E[n], monomer MW), live E[n] = mu1/mu0.
        assert set(pool_stats.keys()) == {"A"}
        e_n, mw = pool_stats["A"]
        assert e_n == pytest.approx(5.0)
        assert mw == pytest.approx(28.0)
        # proxy_event_mass_total: engine-attributed CANONICAL PROXIES only
        # (species_to_pool_indices + is_pool_proxy); attributing the
        # ordinary R is the python ledger's job (spec §4.1).
        assert proxy_total == pytest.approx(gross["A"] * 5.0 * 28.0, rel=1e-12)

    def test_snapshot_e_n_calibration_dominates_fragment_mw(self):
        """Spec §7.6 (decision 3): one mole of representative production is
        one mole of EVENTS; the mass entering the motif class per event is a
        chain's worth (E[n]*monomer_MW), not the representative fragment's
        own MW. With parent-pool E[n]=60 and a ~3-monomer-sized
        representative (hexane vs C2H4 repeat unit) the calibrated
        event-mass must read ~20x a fragment-MW accounting."""
        sp, core, mask = _one_pool_gate_species(rep_smiles="CCCCCC")
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        # E[n] = 60; realizable: mu0*mu2 = 3700 >= mu1^2 = 3600.
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 60.0, 3700.0), monomer_mw_g_mol=28.0)

        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        gross, pool_stats, _ = rs.spawn_gate_flux_snapshot()

        prod = float(rs.core_species_production_rates[4])
        assert prod > 0.0
        e_n, mw = pool_stats["A"]
        g_r = gross["R"] * e_n * mw  # the gate's representative g_i (spec §3)
        assert g_r == pytest.approx(prod * 60.0 * 28.0, rel=1e-12)
        rep_mw_g_mol = sp["R"].molecule[0].get_molecular_weight() * 1000.0  # ~86.18
        ratio = g_r / (prod * rep_mw_g_mol)
        assert ratio == pytest.approx(60.0 * 28.0 / rep_mw_g_mol, rel=1e-12)
        assert 15.0 < ratio < 25.0, "E[n]-calibrated mass must read ~20x the fragment-MW accounting"

    def test_snapshot_mu0_exhaustion_defers_not_inflates(self):
        """Spec §7.7: mu0 <= SMALL_EPS with tiny mu1 -> pool_stats E[n]
        clamps to 0 -> g_i = 0 -> the gate DEFERS. Asserts the deferral
        DIRECTION, not just finiteness: the naive mu1/mu0 would explode
        toward +inf and wave the motif through. Note the amended split:
        gross itself stays nonzero (it is the raw production record); the
        zeroing lives in pool_stats.
        """
        from rmgpy.polymer import (MassFluxAccumulator, MotifLedgerEntry,
                                   Polymer, discover_repeat_motif,
                                   process_polymer_candidates_multipool)

        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        # mu0 exhausted (0 <= SMALL_EPS), mu1 tiny but nonzero.
        rs = _one_pool_gate_rs(rxn, core, mask, (0.0, 1e-25, 1e-20))

        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        snapshot = rs.spawn_gate_flux_snapshot()
        gross, pool_stats, proxy_total = snapshot

        assert gross["A"] > 0.0, (
            "gross production is nonzero — only the E[n] clamp zeroes g_i"
        )
        assert pool_stats["A"][0] == 0.0, "E[n] must clamp to 0 under SMALL_EPS"
        assert pool_stats["A"][1] == pytest.approx(28.0)
        assert proxy_total == 0.0

        # Feed the engine snapshot to the gate: the exhausted pool must DEFER.
        parent = Polymer(label="PE", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"],
                         cutoff=3, Mn=1000.0, Mw=2500.0, initial_mass=1.0)
        cand = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        cand.label = "phenolic_arrival"
        motif = discover_repeat_motif(cand.molecule[0])
        assert motif is not None

        class _Model:
            pass

        model = _Model()
        model.polymer_motif_ledger = [MotifLedgerEntry(
            motif=motif, accumulator_key="motif-0",
            representatives=[("R", "A")],  # parent pool "A" recorded at absorption
        )]
        model.polymer_flux_snapshot = snapshot

        _, intents = process_polymer_candidates_multipool(
            candidates=[cand],
            reaction_model=model,
            pool_registry=[parent],
            iteration=2,
            flux_accumulator=MassFluxAccumulator(window=3),
        )
        assert intents == [], (
            "a mu0-exhausted pool must defer the spawn (epsilon errs toward deferral)"
        )

    def test_duplicate_core_label_raises_not_masks(self):
        """Label uniqueness is NOT enforced on the standard non-RMS path
        (model.py's dedup loop is gated on edge.phase_system) — a duplicate
        core label must RAISE, never silently overwrite: an overwritten
        gross entry would misattribute spawn-gate flux without any signal.
        The main.py stash turns the raise into warn + None snapshot, so the
        gate defers honestly (spec section 4.5)."""
        sp, core, mask = _one_pool_gate_species()
        sp["G"].label = "R"  # collide the gas species with the ordinary R
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0))
        with pytest.raises(ValueError, match="duplicate core species label"):
            rs.spawn_gate_flux_snapshot()


class TestAttributionTrustFloor:
    """Attribution trust floor (item #14a, spec 2026-06-11 §3): the spawn-gate
    snapshot distrusts E[n] = mu1/mu0 for pools whose mu0 sits in the
    integrator-noise band (SMALL_EPS, max(SMALL_EPS, 100 * atol_mu0)].
    SMALL_EPS keeps its realizability job everywhere else — two constants,
    two jobs (T3 exhibits the divergence directly).
    """

    # _one_pool_gate_rs initializes with the base default atol=1e-16
    # (vector atol_array = np.ones(neq) * 1e-16, base.pyx:390), so the
    # trust floor is max(1e-30, 100 * 1e-16) = 1e-14 on every slot.
    TRUST_FLOOR = 100.0 * 1e-16

    def test_band_mu0_zeroes_snapshot_e_n(self):
        """T1 (band explosion, spec §6): a pool with mu0 inside
        (SMALL_EPS, K*atol] MUST report E[n] = 0.0 in pool_stats."""
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        # mu0 = 1e-20 (in the band); mu1 = 1e-12 -> raw E[n] = 1e8 (absurd).
        # Realizable: mu0*mu2 = 1e-23 >= mu1^2 = 1e-24.
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0e-20, 1.0e-12, 1.0e-3))

        rs.residual(0.0, rs.y, np.zeros_like(rs.y))

        # LIVENESS PINS — BEFORE the red assertion (tripwire discipline):
        # (1) mu0 genuinely ABOVE SMALL_EPS: it passes the OLD guard, so the
        #     red below can only mean "trust floor absent", never
        #     "fixture dead / mu0 exhausted".
        mu0 = float(rs.y[1])
        assert mu0 == pytest.approx(1.0e-20) and mu0 > 1.0e-30, (
            "FIXTURE BROKEN, not a valid red: mu0 must sit ABOVE SMALL_EPS"
        )
        # (2) the raw ratio is genuinely huge — far beyond any physical
        #     chain length — so zeroing it is a TRUST verdict, not noise.
        raw_ratio = float(rs.y[2]) / mu0
        assert raw_ratio == pytest.approx(1.0e8) and raw_ratio > 1.0e6, (
            "FIXTURE BROKEN, not a valid red: raw mu1/mu0 must be absurdly "
            "large for the band-explosion red to mean anything"
        )
        # (3) the band sits below the tolerance-anchored floor.
        assert mu0 < self.TRUST_FLOOR

        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()
        # THE red assertion: pre-change HEAD returns the laundered 1e8.
        assert pool_stats["A"][0] == 0.0, (
            "noise-band mu0 passed the attribution guard: E[n] = mu1/mu0 "
            "laundered into pool_stats (trust floor absent)"
        )
        # The zeroed pool contributes nothing to the proxy denominator
        # either, and gross stays a raw production record.
        assert proxy_total == 0.0
        assert gross["A"] > 0.0

    def test_two_constants_two_jobs_divergence(self):
        """T3 (spec §6): on the SAME band state the snapshot reports
        E[n] = 0 while the solver's own E[n] machinery (and the residual)
        still compute with SMALL_EPS realizability semantics — the gate
        distrusts what the solver still legally computes."""
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0e-20, 1.0e-12, 1.0e-3))
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(np.isfinite(dn_dt))  # the solver legally evaluates

        # Solver-side E[n]: _chain_bundle's uniform pick keeps the SMALL_EPS
        # realizability guard — mu0 = 1e-20 is NOT "too empty to move a
        # chain" (b0 = 1) and b1 is the raw 1e8 ratio, legally computed.
        b0, b1, b2, mu2_ok = rs._chain_bundle(0, rs.y, 1.0, True)
        assert b0 == 1.0
        assert b1 == pytest.approx(1.0e8)

        # Gate-side: the snapshot distrusts the SAME state.
        _, pool_stats, _ = rs.spawn_gate_flux_snapshot()
        assert pool_stats["A"][0] == 0.0

    def test_all_pools_in_band_fraction_zero_no_nan(self):
        """T4 (spec §6): every pool in the band -> motif numerator AND the
        whole denominator collapse to zero -> _spawn_gate_fraction returns
        0.0 via the polymer.py:2387 `denominator <= 0.0` guard — fraction
        0.0, defer, no NaN. The deep-conversion TGA-tail regime, pinned
        end-to-end."""
        from rmgpy.polymer import (MotifLedgerEntry, _spawn_gate_fraction,
                                   discover_repeat_motif)
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0e-20, 1.0e-12, 1.0e-3))
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        snapshot = rs.spawn_gate_flux_snapshot()
        gross, pool_stats, proxy_total = snapshot
        assert gross["R"] > 0.0           # flux exists...
        assert pool_stats["A"][0] == 0.0  # ...but the only pool is distrusted
        assert proxy_total == 0.0

        cand = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        motif = discover_repeat_motif(cand.molecule[0])
        assert motif is not None
        entry = MotifLedgerEntry(motif=motif, accumulator_key="motif-0",
                                 representatives=[("R", "A")])
        fraction = _spawn_gate_fraction(entry, [entry], snapshot)
        assert not math.isnan(fraction)
        assert fraction == 0.0

    def test_dying_pool_trajectory_phases_to_zero_before_small_eps(self):
        """T5 (spec §6): a pool draining monotonically through the band
        records fractions that collapse to 0 while mu0 is still decades
        ABOVE SMALL_EPS, and the sum/3 window statistic shows phased decay
        (zero records mix in as the pool crosses the floor) — decay through
        the window, not a cliff at SMALL_EPS."""
        from rmgpy.polymer import (MassFluxAccumulator, MotifLedgerEntry,
                                   _spawn_gate_fraction, discover_repeat_motif)
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        # Cycle 1 starts trusted: mu0 = 1e-6 >> trust floor 1e-14, E[n] = 50.
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0e-6, 5.0e-5, 3.0e-3))

        cand = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        motif = discover_repeat_motif(cand.molecule[0])
        entry = MotifLedgerEntry(motif=motif, accumulator_key="motif-0",
                                 representatives=[("R", "A")])
        acc = MassFluxAccumulator(window=3)

        # mu0 ladder: trusted, then three band states — ALL above SMALL_EPS.
        mu0_ladder = [1.0e-6, 1.0e-16, 1.0e-18, 1.0e-20]
        fractions, stats = [], []
        for it, mu0 in enumerate(mu0_ladder, start=1):
            rs.y[1] = mu0
            # band mu1 = 1e-12 -> raw ratio up to 1e8 (would have exploded
            # pre-floor); trusted cycle keeps the physical E[n] = 50.
            rs.y[2] = mu0 * 50.0 if it == 1 else 1.0e-12
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
            assert float(rs.y[1]) > 1.0e-30  # NEVER reaches SMALL_EPS
            snapshot = rs.spawn_gate_flux_snapshot()
            f = _spawn_gate_fraction(entry, [entry], snapshot)
            acc.record(entry.accumulator_key, f, it)
            fractions.append(f)
            stats.append(acc.gate_statistic(entry.accumulator_key))

        # Cycle 1 (trusted): a genuine nonzero fraction in (0, 1].
        assert 0.0 < fractions[0] <= 1.0
        # Cycles 2-4 (band): fractions collapse to EXACTLY 0 long BEFORE
        # mu0 <= SMALL_EPS.
        assert fractions[1:] == [0.0, 0.0, 0.0]
        # Window statistic: phased decay — monotone nonincreasing, holding
        # f1/3 while the trusted record remains in the window, reaching 0
        # only when it ages out. No cliff at the floor: the trusted record
        # must SURVIVE in the window through cycles 2-3 (a band entry that
        # wiped the history would show [f1/3, 0, 0, 0] and is ruled out).
        assert stats[0] == pytest.approx(fractions[0] / 3.0)
        assert stats[1] == pytest.approx(stats[0])
        assert stats[2] == pytest.approx(stats[0])
        assert all(s2 <= s1 + 1e-15 for s1, s2 in zip(stats, stats[1:]))
        assert stats[-1] == 0.0

    def test_census_line_once_per_snapshot_and_silent_when_trusted(self, caplog):
        """T6 (spec §6): the band pool emits EXACTLY ONE
        SPAWN-GATE ATTRIBUTION CENSUS warning per snapshot, carrying pool
        label / mu0 / trust floor / motif count; a trusted-mu0 build emits
        none; an EXHAUSTED pool (mu0 <= SMALL_EPS) emits none either —
        distinguishable from genuinely-zero flux."""
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1

        rs = _one_pool_gate_rs(rxn, core, mask, (1.0e-20, 1.0e-12, 1.0e-3))
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        with caplog.at_level(logging.WARNING):
            caplog.clear()
            rs.spawn_gate_flux_snapshot(motif_counts_by_pool={"A": 2})
        census = [r for r in caplog.records
                  if "SPAWN-GATE ATTRIBUTION CENSUS" in r.getMessage()]
        assert len(census) == 1  # once per snapshot per pool
        msg = census[0].getMessage()
        assert "pool A" in msg
        assert "1.000000e-20" in msg  # mu0
        assert "1.000000e-14" in msg  # trust floor = 100 * atol 1e-16
        assert "motifs attributed to this pool: 2" in msg

        # Trusted mu0: silent.
        rs2 = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0))
        rs2.residual(0.0, rs2.y, np.zeros_like(rs2.y))
        with caplog.at_level(logging.WARNING):
            caplog.clear()
            rs2.spawn_gate_flux_snapshot()
        assert not any("SPAWN-GATE ATTRIBUTION CENSUS" in r.getMessage()
                       for r in caplog.records)

        # Exhausted mu0 (<= SMALL_EPS): silent too — deferral was already
        # honest there; the census marks ONLY the distrust band.
        rs3 = _one_pool_gate_rs(rxn, core, mask, (0.0, 1e-25, 1e-20))
        rs3.residual(0.0, rs3.y, np.zeros_like(rs3.y))
        with caplog.at_level(logging.WARNING):
            caplog.clear()
            rs3.spawn_gate_flux_snapshot()
        assert not any("SPAWN-GATE ATTRIBUTION CENSUS" in r.getMessage()
                       for r in caplog.records)


class TestThermoReferenceStateTripwire:
    """Build-time thermo reference-state tripwire (spec 2026-06-11 §§5-8).

    RED-FIRST discipline: the refusal and provenance tests below were written
    and confirmed FAILING on pre-change HEAD for their PINNED reasons (refusal:
    initialize_model completes without ValueError; provenance: warning absent)
    before the check landed. The liveness pins come BEFORE each red assertion
    so the red can only mean "check absent", never "fixture dead".
    """

    def test_refusal_on_genuine_unpaired_reference_state(self):
        """Spec §8.1: a genuine unpaired reaction (melt chain <=> all-gas
        products, reversible) must REFUSE at initialize_model with the
        cliff-sign ValueError."""
        sp, core, mask = _refstate_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]], **_REV_KIN)

        # LIVENESS PIN -- BEFORE the red assertion. The fixture must actually
        # cross the phase boundary with chain-scale U: reversible, melt
        # reactant, all-gas products, and an INDEPENDENTLY recomputed
        # U = S_trans(A)/(R ln10) + log10(P0/(R*T*C0)) above the 3.0-decade
        # refuse bound (computed: 9.404 + 1.177 = 10.58). A failure HERE
        # means the fixture is dead, not a valid red.
        assert rxn.reversible
        assert not mask[0] and mask[4] and mask[5]
        mw_a = sp["A"].molecule[0].get_molecular_weight()
        u_expected = (_sackur_tetrode_decades(mw_a, 800.0)
                      + math.log10(1.0e5 / (constants.R * 800.0 * 1.0)))
        assert u_expected > 3.0, (
            "FIXTURE BROKEN, not a valid red: independently recomputed U "
            f"({u_expected:.2f}) is not above the refuse bound"
        )

        # THE red assertion: on pre-change HEAD initialize_model COMPLETES
        # (verified live: the Keq path survives on the trivial NASA thermo,
        # kb finite) -- no tripwire exists to refuse. The match string pins
        # the pinned reason; an unrelated ValueError cannot launder a pass.
        with pytest.raises(ValueError, match="unpaired reference-state"):
            _refstate_rs(core, [rxn], mask,
                         [_gate_pool_config()], {"A": (1.0, 5.0, 30.0)})

    def test_mixed_provenance_chain_counterparty_warns(self, caplog):
        """Spec §8.4: one melt-class species takes library thermo while its
        chain-scale counterparty takes GAV -> the mixed-provenance warning
        must fire. The counterparty is the §2 decoupling shape: a
        gas-classified but proxy-TAGGED same-mass radical (butyl 57.11 vs
        butane 58.12 g/mol, inside the one-monomer window 28+10), so the pair
        is mass-paired (U = 1.5*log10(58.12/57.11) = 0.011 -- no census, no
        refusal) and ONLY the sensor can speak. The small gas H also takes
        library thermo and must NOT matter (outside the window)."""
        import logging
        sp = {
            "A": _spc("CCCC", "A"),
            "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"), "A_mu2": _spc("C#N", "A_mu2"),
            "B": _spc("[CH2]CCC", "B"),
            "H": _spc("[H]", "H"),
        }
        sp["A"].thermo = _trivial_nasa(_LIB_COMMENT)   # melt proxy: library
        sp["B"].thermo = _trivial_nasa(_GAV_COMMENT)   # chain counterparty: GAV
        sp["H"].thermo = _trivial_nasa(_LIB_COMMENT)   # small gas: library
        for k in ("A_mu0", "A_mu1", "A_mu2"):
            sp[k].thermo = _trivial_nasa(_LIB_COMMENT)
        # gas-CLASSIFIED, proxy-TAGGED, chain-scale (57.11 >= window
        # 28 + 10 = 38 from the pool config below): physically-melt via the
        # tag branch of the C3-AMENDED spec-§5.1 class (tag AND MW >= window)
        sp["B"].is_polymer_proxy = True
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["B"], sp["H"]]
        mask = np.array([False] * 4 + [True, True], dtype=bool)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["H"]], **_REV_KIN)

        # LIVENESS PINS: the pair is mass-paired (inside the counterparty
        # window) and the build must SUCCEED on any head -- only the warning
        # distinguishes pre/post-change.
        mw_a = sp["A"].molecule[0].get_molecular_weight() * 1000.0
        mw_b = sp["B"].molecule[0].get_molecular_weight() * 1000.0
        assert abs(mw_a - mw_b) <= 28.0 + 10.0
        # amended-class pin (spec §5.1 C3): B must clear the chain-scale
        # window, or it would not be physically-melt at all and the sensor
        # would have no melt-vs-counterparty pair to warn on
        assert mw_b >= 28.0 + 10.0
        with caplog.at_level(logging.WARNING):
            _refstate_rs(core, [rxn], mask,
                         [_gate_pool_config()], {"A": (1.0, 5.0, 30.0)})

        # THE red assertion: pre-change HEAD emits no provenance warning.
        assert any("THERMO REFERENCE-STATE PROVENANCE" in r.getMessage()
                   for r in caplog.records), (
            "mixed library-vs-GAV provenance among chain-scale counterparties "
            "went unwarned (the spec-§5.3 decoupling-fingerprint sensor is "
            "absent)"
        )

    def test_melt_sum_leak_guard_raises_classification_error(self):
        """Spec §5.1 C3 amendment: the cannot-happen leak guard inside the
        melt-sum accumulation. Under the amended class a tagged below-window
        species fails the MW conjunct and is EXCLUDED by the gate -- expected
        and silent (the family.py:1657 over-tagging fingerprint, H2 on every
        proxy-touching reaction); the raise is only for such a species
        REACHING the melt sum. Because the gate and the guard share ONE
        window definition, the violation is structurally unreachable through
        public paths -- so the raise is pinned by calling the helper directly
        with a hand-built violating member (documented bypass, per the C3
        amendment)."""
        from rmgpy.solver.polymer import _assert_chain_scale_melt_member
        # Valid members never raise: a condensed-branch member (gate-exempt,
        # pool-configured by input -- any MW, even below the window) and a
        # chain-scale tag-branch member (MW above the window).
        _assert_chain_scale_melt_member("M1", 0.016043, False, 0.038)
        _assert_chain_scale_melt_member("Erad", 0.211407, True, 0.080)
        # The hand-built violation: a gas-classified (tag-branch) member
        # below the window inside the melt sum. The message must steer the
        # operator to CLASSIFICATION, never to reference states.
        with pytest.raises(ValueError, match="classification leak, NOT a thermo problem"):
            _assert_chain_scale_melt_member("H2", 0.002016, True, 0.038)
        with pytest.raises(ValueError, match="non-chain species in the melt sum"):
            _assert_chain_scale_melt_member("H2", 0.002016, True, 0.038)

    def test_unpaired_decades_formula_pins_spec_numbers(self):
        """Spec §5.1/§2 pinned numerically against the production helper:
        the C15H32 EPDM proxy gives S_trans/(R ln10) = 10.49 decades at
        1000 K plus the C0 term 1.08 -> a genuine melt-chain => all-gas
        reaction carries U = 11.57 (the §2 measurement, exact); the
        same-mass C15H32/C15H31 melt pair collapses to 0.0031 decades (the
        structural cancellation, <= the spec's 0.03)."""
        from rmgpy.solver.polymer import _unpaired_reference_decades
        mw_chain = Molecule().from_smiles("CCC(C)CCCC(C)CCCC(C)C").get_molecular_weight()
        mw_rad = Molecule().from_smiles("CCC(C)CCCC(C)CCCC(C)[CH2]").get_molecular_weight()
        u_unpaired = _unpaired_reference_decades([mw_chain], [], 1000.0)
        assert u_unpaired == pytest.approx(11.571, abs=0.01)
        u_paired = _unpaired_reference_decades([mw_chain], [mw_rad], 1000.0)
        assert u_paired == pytest.approx(0.0031, abs=0.001)
        assert u_paired < 0.03

    def test_census_fires_above_half_decade_not_below(self, caplog):
        """Spec §8.5/§6: census at U > 0.5 and not below. Two melt<=>melt
        reactions in one build (Dn_melt = 0, so U = 1.5*log10(MW ratio),
        T-independent): CH4 => C2H6 gives U = 0.409 (below; silent) and
        CH4 => C4H10 gives U = 0.839 (above; census, no refusal)."""
        import logging
        m1, m2, m4, g = (_spc("C", "M1"), _spc("CC", "M2"),
                         _spc("CCCC", "M4"), _spc("[H][H]", "G"))
        for s in (m1, m2, m4, g):
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [m1, m2, m4, g]
        mask = np.array([False, False, False, True], dtype=bool)
        rxn_below = Reaction(reactants=[m1], products=[m2], **_REV_KIN)
        rxn_above = Reaction(reactants=[m1], products=[m4], **_REV_KIN)
        with caplog.at_level(logging.WARNING):
            rs = _refstate_rs(core, [rxn_below, rxn_above], mask, [], {})
        census = [r for r in caplog.records
                  if "THERMO REFERENCE-STATE CENSUS" in r.getMessage()]
        assert len(census) == 1
        assert "M4" in census[0].getMessage()
        assert len(rs.reference_state_census) == 1
        assert rs.reference_state_census[0][1] == pytest.approx(0.839, abs=0.002)
        # reference_state_max_decades tracks ALL reversible melt-touching
        # reactions (census bound or not); here the max is the 0.839 case.
        assert rs.reference_state_max_decades == pytest.approx(0.839, abs=0.002)

    def test_refusal_fires_above_three_decades_not_below(self, caplog):
        """Spec §8.5/§6: refusal at U > 3.0 and not below. melt CH4 => melt
        C91H184 (MW 1278.4) gives U = 1.5*log10(1278.42/16.043) = 2.852
        (census only, builds); melt CH4 => melt C120H242 (MW 1685.2) gives
        U = 1.5*log10(1685.21/16.043) = 3.032 (> 3.0, refuses). The census
        is emitted even on the refusing build (spec §7: census regardless
        of pass/fail)."""
        import logging

        def _build(n_carbons):
            m1 = _spc("C", "M1")
            big = _spc("C" * n_carbons, "BIG")
            g = _spc("[H][H]", "G")
            for s in (m1, big, g):
                s.thermo = _trivial_nasa(_GAV_COMMENT)
            core = [m1, big, g]
            mask = np.array([False, False, True], dtype=bool)
            rxn = Reaction(reactants=[m1], products=[big], **_REV_KIN)
            return _refstate_rs(core, [rxn], mask, [], {})

        with caplog.at_level(logging.WARNING):
            rs = _build(91)
        assert rs.reference_state_max_decades == pytest.approx(2.852, abs=0.005)
        assert any("THERMO REFERENCE-STATE CENSUS" in r.getMessage()
                   for r in caplog.records)

        caplog.clear()
        with caplog.at_level(logging.WARNING):
            with pytest.raises(ValueError, match="unpaired reference-state"):
                _build(120)
        assert any("THERMO REFERENCE-STATE CENSUS" in r.getMessage()
                   for r in caplog.records)

    def test_override_knob_builds_and_still_logs_census(self, caplog):
        """Spec §8.1/§7: the SAME genuine-unpaired fixture as the refusal
        test builds with allow_unpaired_reference_state=True; the census and
        the explicit bypass warning are still emitted (the override silences
        ONLY the refusal)."""
        import logging
        sp, core, mask = _refstate_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]], **_REV_KIN)
        with caplog.at_level(logging.WARNING):
            rs = _refstate_rs(core, [rxn], mask, [_gate_pool_config()],
                              {"A": (1.0, 5.0, 30.0)},
                              rs_kwargs={"allow_unpaired_reference_state": True})
        assert rs.reference_state_max_decades > 3.0
        assert any("THERMO REFERENCE-STATE CENSUS" in r.getMessage()
                   for r in caplog.records)
        assert any("allow_unpaired_reference_state=True" in r.getMessage()
                   for r in caplog.records)


class TestAllowUnpairedReferenceStateKnobPlumb:
    """Spec §7/§10 -- the override knob's exact path: input.py
    hybrid_polymer_reactor -> HybridPolymerReactor -> to_solver_object ->
    HybridPolymerSystem kwarg. Mirrors the constant_gas_volume pattern (the
    spec's named 'strict_phase_check pattern' does not exist in the code --
    plan contradiction C1; debug_check_realizability has no input plumb)."""

    @staticmethod
    def _reactor(**kwargs):
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import HybridPolymerReactor, PolymerPhase
        a = _spc("CCCC", "A")
        phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"), initial_moments={},
                             initial_explicit={a: 1.0}, pools=[], mass_transfer=[])
        return a, HybridPolymerReactor(
            temperature=(800.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles={a: 1.0}, polymerPhase=phase,
            terminationTime=(1.0, "s"), **kwargs)

    def test_to_solver_object_passes_knob(self):
        a, reactor = self._reactor(allow_unpaired_reference_state=True)
        solver = reactor.to_solver_object([a], [], [], [])
        assert solver.allow_unpaired_reference_state is True

    def test_knob_default_false_and_input_function_carries_it(self):
        import inspect
        a, reactor = self._reactor()
        solver = reactor.to_solver_object([a], [], [], [])
        assert solver.allow_unpaired_reference_state is False
        from rmgpy.rmg.input import hybrid_polymer_reactor
        sig = inspect.signature(hybrid_polymer_reactor)
        assert "allow_unpaired_reference_state" in sig.parameters
        assert sig.parameters["allow_unpaired_reference_state"].default is False
        # Static pin of the forwarding line: executing hybrid_polymer_reactor
        # needs module-global rmg/species_dict state, so the source pin keeps
        # this test cheap while making a silent kwarg drop impossible.
        src = inspect.getsource(hybrid_polymer_reactor)
        assert "allow_unpaired_reference_state=allow_unpaired_reference_state" in src


class TestThermoReferenceStateEpdmShaped:
    """Spec §8.2/§8.3 -- the EPDM shape stays quantitatively clean.

    Fixture geometry mirrors the real deck (T = 1000 K): proxy C15H32
    (212.41 g/mol, GAV), same-length radical C15H31 (211.41 g/mol, GAV,
    gas-classified + proxy-tagged), H/H2 library. Counterparty window =
    70 (pool monomer) + 10 = 80 g/mol: Erad is inside (|dMW| = 1.0), H/H2
    far outside (|dMW| >= 209)."""

    @staticmethod
    def _build():
        c15h32 = "CCC(C)CCCC(C)CCCC(C)C"
        c15h31 = "CCC(C)CCCC(C)CCCC(C)[CH2]"
        sp = {
            "E": _spc(c15h32, "E"),
            "E_mu0": _spc("CO", "E_mu0"), "E_mu1": _spc("C=O", "E_mu1"), "E_mu2": _spc("C#N", "E_mu2"),
            "Erad": _spc(c15h31, "Erad"),
            "H": _spc("[H]", "H"), "H2": _spc("[H][H]", "H2"),
        }
        for k in ("E", "Erad"):
            sp[k].thermo = _trivial_nasa(_GAV_COMMENT)
        for k in ("E_mu0", "E_mu1", "E_mu2", "H", "H2"):
            sp[k].thermo = _trivial_nasa(_LIB_COMMENT)
        # The §2 shape: the same-length abstraction radical is gas-CLASSIFIED
        # but proxy-TAGGED (the chain-variant judgment the spawn-pass
        # machinery stamps in production; family.py:1657 -> model.py:486 ->
        # multipool re-classification).
        sp["Erad"].is_polymer_proxy = True
        core = [sp["E"], sp["E_mu0"], sp["E_mu1"], sp["E_mu2"],
                sp["Erad"], sp["H"], sp["H2"]]
        mask = np.array([False] * 4 + [True, True, True], dtype=bool)
        rxn = Reaction(
            reactants=[sp["E"], sp["H"]], products=[sp["Erad"], sp["H2"]],
            kinetics=Arrhenius(A=(1.0e3, "m^3/(mol*s)"), n=0.0,
                               Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=True)
        pool = dataclasses.replace(_gate_pool_config(monomer_mw_g_mol=70.0),
                                   label="E")
        rs = HybridPolymerSystem(
            T=1000.0, P=1.0e5,
            initial_mole_fractions={sp["H"]: 0.01, sp["H2"]: 0.0},
            V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"E": (1.0, 5.0, 30.0)}, termination=[],
        )
        rs.initialize_model(core, [rxn], [], [])
        return rs

    def test_max_u_below_benign_ceiling(self, caplog):
        """Spec §8.2: the boundary-crossing H-abstraction pair comes in at
        the paired-cancellation scale -- max U <= 0.33 + margin ASSERTED
        (not just 'no exception'); in fact 1.5*log10(212.41/211.41) =
        0.0031 decades. The census is empty."""
        import logging
        with caplog.at_level(logging.WARNING):
            rs = self._build()
        assert rs.reference_state_max_decades <= 0.33 + 0.1
        assert rs.reference_state_max_decades == pytest.approx(0.0031, abs=0.002)
        assert rs.reference_state_census == []
        assert not any("THERMO REFERENCE-STATE CENSUS" in r.getMessage()
                       for r in caplog.records)

    def test_zero_provenance_warnings(self, caplog):
        """Spec §8.3 green assertion that PINS the narrow counterparty
        scope: the reaction genuinely mixes library (H, H2) and GAV (E,
        Erad) provenance, but H/H2 sit far outside the one-monomer MW
        window, so the sensor must stay SILENT. The broad 'all
        co-participants' definition would warn here -- and on all 26 EPDM
        reactions on day one (the false-positive storm spec §5.3 forbids)."""
        import logging
        from rmgpy.solver.polymer import _thermo_provenance
        with caplog.at_level(logging.WARNING):
            rs = self._build()
        # LIVENESS PINS -- BEFORE the silence assertion, so silence cannot
        # be the silence of a dead pass: (a) the tripwire ran and visited
        # the reaction (the same 0.0031 the sibling test pins); (b) the
        # provenance classifier genuinely sees the library/GAV mix on the
        # fixture's comment strings -- the sensor stays silent on SCOPE
        # (H/H2 outside the MW window), not on blindness.
        assert rs.reference_state_max_decades == pytest.approx(0.0031, abs=0.002)
        probe_lib, probe_gav = _spc("C", "PL"), _spc("CC", "PG")
        probe_lib.thermo = _trivial_nasa(_LIB_COMMENT)
        probe_gav.thermo = _trivial_nasa(_GAV_COMMENT)
        assert _thermo_provenance(probe_lib) == "library"
        assert _thermo_provenance(probe_gav) == "gav"
        assert not any("THERMO REFERENCE-STATE PROVENANCE" in r.getMessage()
                       for r in caplog.records)


# ---------------------------------------------------------------------------
# Item 17 (spec 2026-06-12): phase-gate / enlargement consistency fixtures.
# Ported from /tmp/census_probe.py (the Q0 census probes) — same species,
# same kinetics, same moments; two solver builds per case (E = product in
# edge, C = product promoted to core).
# ---------------------------------------------------------------------------

def _gate17_species():
    """Census-probe species set. Pool A (proxy + mu0/1/2), inert gas seed X,
    gas driver target Y, product-under-test G, would-be pool-G mu dummies,
    ordinary condensed R."""
    return {
        "A": _spc("CCCC", "A"),
        "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"),
        "A_mu2": _spc("C#N", "A_mu2"),
        "X": _spc("N#N", "X"),
        "Y": _spc("C", "Y"),
        "G": _spc("[CH3]", "G"),
        "G_mu0": _spc("CCO", "G_mu0"), "G_mu1": _spc("CC=O", "G_mu1"),
        "G_mu2": _spc("CC#N", "G_mu2"),
        "R17": _spc("CCCO", "R17"),
    }


def _stage1_classifier(species_list):
    """Item 17 A5-2 fixture stage-1 classifier: f(species_list) -> bool array
    (True = gas). Stands in for the bound polymerPhase.get_gas_mask on the
    production live-edge path -- it classifies the product-under-test G as
    CONDENSED by label (exactly as get_gas_mask condenses a registered pool
    member), and everything else GAS. The KEY property the armed-row mutation
    proof rests on: G's condensed verdict here comes from STAGE 1 (this
    classifier), NOT from a stage-2 pool-label override -- so wiring None
    (forced fallback) genuinely flips edge G to GAS and breaks the umbrella
    parity, the way the stale-seed fallback silently would in production."""
    return np.array([s.label != "G" for s in species_list], dtype=bool)


def _gate17_rs(core, mask, rxns_core, edge_spcs=(), rxns_edge=(),
               pools=(("A", (1, 2, 3)),), moments=None):
    """Build + initialize a HybridPolymerSystem for the item-17 fixtures.

    ``pools`` is a tuple of (label, mu_indices) — the §5 config-state axis:
    adding ("G", ...) is the solver-level expression of an item-16 daughter
    config (spec §3(a) stage-2 labels)."""
    moments = moments if moments is not None else {
        lbl: (1.0, 5.0, 30.0) for lbl, _ in pools}
    mask_arr = np.array(mask, dtype=bool)
    seed_idx = int(np.where(mask_arr)[0][0])
    pool_cfgs = [PolymerPoolConfig(label=lbl, xs=2,
                                   explicit_dp_to_species_index={},
                                   mu_indices=mu, monomer_poly_index=None)
                 for lbl, mu in pools]
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[seed_idx]: 1.0},
        V_poly=1.0, polymer_pools=pool_cfgs, mass_transfer=[],
        gas_species_mask=mask_arr, constant_gas_volume=False,
        initial_polymer_moments=moments, termination=[],
        # Item 17 A5-2: a direct-test build with no blueprint phase object is
        # a legitimate last-resort fallback -- flag it so R1-EDGE (the
        # edge-suffix provenance guard) does not raise on its default-filled
        # edge. (Fixtures express prospectively-condensed edge species via
        # configured pool labels in stage 2; the production live-edge stage-1
        # path is exercised by the classifier-wired fixtures in T12/T13.)
        allow_default_prospective_edge=True)
    rs.initialize_model(list(core), list(rxns_core), list(edge_spcs),
                        list(rxns_edge))
    return rs


class TestUmbrellaPhaseGateParity:
    """Spec 2026-06-12 §5 — THE umbrella invariant, one test:

        promotion-time flux must equal post-promotion flux under FULL
        post-promotion semantics (mask + config).

    Item 17 delivers the mask projection; item 16 RE-RUNS this exact
    parameterized test under engine-created configs (the parameterization
    axis IS 16's activation stage). xfail(strict=True) rows are the
    promoted-then-zeroed census shapes, RED at pre-17 HEAD for their pinned
    reasons; the §3(c) gate rewrite (Task 3) removes the marks."""

    # Post-17 expected common rate per case (edge == core == this).
    EXPECTED = {
        "B1_head": 0.0,                       # parity via zero
        "B1_configured": 10.0,                # parity via flux (probe Cp)
        "B2_allgas": 0.0,                     # zero, census-loud (Task 4)
        "A_head": 2.0e5 / (constants.R * 800.0),   # 30.068... ungated
        "A_armed": 0.0,                       # Gate A, armed shape
    }

    def _build_pair(self, case):
        sp = _gate17_species()
        pools = [("A", (1, 2, 3))]
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = [False, False, False, False, True]
        g_configured = case in ("B1_configured", "A_armed")
        if g_configured:
            core = core + [sp["G_mu0"], sp["G_mu1"], sp["G_mu2"]]
            mask = mask + [False, False, False]
            pools.append(("G", (5, 6, 7)))
        if case.startswith("B1"):
            rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        elif case == "B2_allgas":
            # Probe 3b: fold-back-less DISCRETE_CHIP — genuinely all-gas
            # products. Throttled site = min(mu0, mu1/2)/V_poly = 1.0.
            rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
            rxn.is_end_group_reaction = True
            rxn.polymer_flux_archetype = 5  # DISCRETE_CHIP
            rxn.polymer_chip_units = 2
        else:  # A-shaped: gas event
            rxn = Reaction(reactants=[sp["X"]], products=[sp["G"]], **_KIN)
        rs_edge = _gate17_rs(core, mask, [], edge_spcs=[sp["G"]],
                             rxns_edge=[rxn], pools=pools)
        # Post-promotion build: G in core. Condensed iff a G config exists
        # (the prospective verdict item 16 will create; HEAD daughters run
        # as ordinary gas — the probed B1 truth).
        rs_core = _gate17_rs(core + [sp["G"]], mask + [not g_configured],
                             [rxn], pools=pools)
        return sp, rs_edge, rs_core

    @pytest.mark.parametrize("case", [
        "B1_head",
        "B1_configured",
        "B2_allgas",
        "A_head",
        "A_armed",
    ])
    def test_umbrella_parity_edge_rate_equals_core_rate(self, case):
        sp, rs_edge, rs_core = self._build_pair(case)
        # Liveness pins FIRST (the red can only mean "edge evaluation uses
        # different semantics than core", never "dead fixture"): kinetics
        # alive at T, and the pool site factor alive for poly events.
        assert float(rs_edge.kf[0]) > 0.0
        if case.startswith(("B1", "B2")):
            assert float(rs_edge.y[2]) > 0.0  # pool A mu1 site factor
        else:
            assert float(rs_edge.y[4]) > 0.0  # gas reactant X
        rs_edge.residual(0.0, rs_edge.y, np.zeros_like(rs_edge.y))
        rs_core.residual(0.0, rs_core.y, np.zeros_like(rs_core.y))
        edge_rate = float(np.asarray(rs_edge.edge_reaction_rates)[0])
        core_rate = float(np.asarray(rs_core.core_reaction_rates)[0])
        # THE UMBRELLA INVARIANT — this assertion dies first, by design.
        assert edge_rate == pytest.approx(core_rate, abs=1e-12), (
            f"umbrella parity broken for {case}: edge rate {edge_rate} vs "
            f"post-promotion core rate {core_rate}")
        assert edge_rate == pytest.approx(self.EXPECTED[case],
                                          rel=1e-9, abs=1e-12)


class TestProspectiveMask:
    """Spec 2026-06-12 §3(a)/(b)/(d) — the prospective mask is the real mask
    evaluated early, and R1 makes that claim self-verifying every build."""

    def test_prospective_mask_prefix_and_length_contract(self):
        """T3: mixed core+edge build -- prospective_gas_mask has length
        n_core + n_edge, its core prefix equals gas_species_mask elementwise
        (R1's green half), edge species default prospectively-GAS in
        fallback mode (A3), and a configured pool label in the EDGE flips
        the prospective verdict while gas_species_mask stays core-sized and
        untouched (R3: second array, never a resize)."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["G_mu0"], sp["G_mu1"], sp["G_mu2"]]
        mask = [False, False, False, False, True, False, False, False]
        edge = [sp["G"], sp["Y"]]
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        rs = _gate17_rs(core, mask, [], edge_spcs=edge, rxns_edge=[rxn],
                        pools=(("A", (1, 2, 3)), ("G", (5, 6, 7))))
        n_core = rs.num_core_species
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        gm = np.asarray(rs.gas_species_mask, dtype=bool)
        assert pm.shape[0] == n_core + len(edge)
        assert gm.shape[0] == n_core  # R3: never a resize
        assert np.array_equal(pm[:n_core], gm)
        # stage 2 over the combined list: configured pool label G lives in
        # the EDGE -- prospectively condensed; unconfigured edge Y stays gas.
        assert pm[n_core + 0] == False  # G: configured -> condensed
        assert pm[n_core + 1] == True   # Y: fallback default GAS (A3)

    def test_stale_stage1_seed_falls_back_loudly(self, caplog):
        """Engine-reuse pin (plan-level decision): a constructor seed sized
        for a DIFFERENT edge list must not crash initialize_model
        (polymer_input.py:176-180 reuses engines across simulate calls);
        it logs PROSPECTIVE-MASK SEED STALE: and takes the documented
        fallback, and R1 still proves the prefix."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = np.array([False, False, False, False, True], dtype=bool)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        # Seed sized for n_core + 3 edge species; build with ONE edge species.
        stale_seed = np.concatenate([mask, np.ones(3, dtype=bool)])
        # A5-2 re-pin: with a stale seed, no classifier, and a non-empty edge,
        # the build takes the edge-defaults-GAS fallback (branch 3) -- a
        # default-filled edge suffix that R1-EDGE would RAISE on for a
        # production build. This is a legitimate direct-test fallback, so it
        # is flagged allow_default_prospective_edge=True; the back-compat
        # PROSPECTIVE-MASK SEED STALE warning is still emitted (a stale seed
        # was present) and the fallback is taken (no provenance raise).
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_gas_mask=stale_seed,
            allow_default_prospective_edge=True)
        with caplog.at_level(logging.WARNING):
            rs.initialize_model(list(core), [], [sp["G"]], [rxn])
        assert any("PROSPECTIVE-MASK SEED STALE:" in r.getMessage()
                   for r in caplog.records)
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert pm.shape[0] == rs.num_core_species + 1
        assert np.array_equal(pm[:rs.num_core_species],
                              np.asarray(rs.gas_species_mask, dtype=bool))

    def test_divergent_seed_prefix_raises_tripwire(self):
        """T4: a doctored stage-1 seed whose core prefix disagrees with
        gas_species_mask on ONE index (a non-pool species, so stage 2 cannot
        repair it) must make initialize_model RAISE with the
        PROSPECTIVE-MASK TRIPWIRE: sentinel naming that species -- the raise
        is live, not decorative (R1 is raise-never-warn)."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = np.array([False, False, False, False, True], dtype=bool)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        doctored = np.concatenate([mask, np.ones(1, dtype=bool)])
        doctored[4] = False  # X: gas in the real mask, condensed in the seed
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_gas_mask=doctored)
        with pytest.raises(ValueError, match="PROSPECTIVE-MASK TRIPWIRE:"):
            rs.initialize_model(list(core), [], [sp["G"]], [rxn])
        # the diverging species is NAMED
        with pytest.raises(ValueError, match=r"X"):
            rs.initialize_model(list(core), [], [sp["G"]], [rxn])


class TestPhaseGateNonRegression:
    """Spec 2026-06-12 SS3(c) precision pins: the rewrite kills ONLY the
    has_edge_prod phase bypass. The reverse-rate concentration-availability
    hole (:1494-1499 -- vanilla simple.pyx parity, Z6) and the chip
    exhaustion throttle (mask-free, edge-inclusive by construction) must be
    provably untouched."""

    def test_reverse_rate_hole_for_edge_products_untouched(self):
        """An edge product has no concentration in y, so the reverse rate is
        UNCOMPUTABLE (not phase-forbidden): rr stays 0 however large kb is.
        Control: the same doctored kb IS live when the product sits in core
        with state -- proving the pin is not vacuous."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = [False, False, False, False, True]
        rxn = Reaction(reactants=[sp["X"]], products=[sp["G"]], **_KIN)
        rs = _gate17_rs(core, mask, [], edge_spcs=[sp["G"]], rxns_edge=[rxn])
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        rate_before = float(np.asarray(rs.edge_reaction_rates)[0])
        assert rate_before > 0.0  # liveness: ungated gas->gas edge flux
        rs.kb[0] = 1.0e6  # doctor the reverse coefficient
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert float(np.asarray(rs.edge_reaction_rates)[0]) == \
            pytest.approx(rate_before, rel=1e-12)
        # control: product in CORE (gas) with nonzero state -> kb is live
        rs2 = _gate17_rs(core + [sp["G"]], mask + [True], [rxn])
        y2 = rs2.y.copy()
        y2[5] = 0.5  # give G concentration so rr = kb*C(G) > 0
        rs2.residual(0.0, y2, np.zeros_like(y2))
        fwd_only = float(np.asarray(rs2.core_reaction_rates)[0])
        rs2.kb[0] = 1.0e6
        rs2.residual(0.0, y2, np.zeros_like(y2))
        assert float(np.asarray(rs2.core_reaction_rates)[0]) < fwd_only

    def test_chip_throttle_edge_core_symmetry_pinned(self):
        """T8: probe 3's HEAD values pinned -- the canonical DISCRETE_CHIP
        (fold-back proxy product) is throttled identically at edge and core
        (2.0 == 2.0), and mu0 = 0 exhausts both to 0.0. The throttle never
        reads any mask; item 17 must not change it (its config-flip cousin
        Z4 is item 16's)."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = [False, False, False, False, True]

        def chip_rxn():
            r = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]],
                         **_KIN)
            r.is_end_group_reaction = True
            r.polymer_flux_archetype = 5  # DISCRETE_CHIP
            r.polymer_chip_units = 2
            return r

        for moments, expected in (((1.0, 5.0, 30.0), 2.0),
                                  ((0.0, 5.0, 30.0), 0.0)):
            rs_e = _gate17_rs(core, mask, [], edge_spcs=[sp["G"]],
                              rxns_edge=[chip_rxn()],
                              moments={"A": moments})
            rs_c = _gate17_rs(core + [sp["G"]], mask + [True], [chip_rxn()],
                              moments={"A": moments})
            rs_e.residual(0.0, rs_e.y, np.zeros_like(rs_e.y))
            rs_c.residual(0.0, rs_c.y, np.zeros_like(rs_c.y))
            e = float(np.asarray(rs_e.edge_reaction_rates)[0])
            c = float(np.asarray(rs_c.core_reaction_rates)[0])
            assert e == pytest.approx(expected, abs=1e-12)
            assert c == pytest.approx(expected, abs=1e-12)


def _gate17_simulate(rs, core, rxns_core, edge_spcs, rxns_edge,
                     tol_move_to_core=1.0e-3, t_end=1.0e-6):
    """Drive a gate-17 fixture through a REAL simulate() -- the hook's only
    habitat (tol_move_to_core is a simulate local, base.pyx:635; char_rate
    exists only per accepted snapshot)."""
    from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
    from rmgpy.solver.base import TerminationTime
    rs.termination.append(TerminationTime((t_end, "s")))
    ms = ModelSettings(tol_keep_in_edge=0.0,
                       tol_move_to_core=tol_move_to_core,
                       tol_interrupt_simulation=1.0e8)
    rs.simulate(list(core), list(rxns_core), list(edge_spcs),
                list(rxns_edge), [], [],
                model_settings=ms, simulator_settings=SimulatorSettings())


def _census_lines(caplog):
    return [r.getMessage() for r in caplog.records
            if "PHASE-GATE FLUX CENSUS:" in r.getMessage()]


class TestPhaseGateFluxCensus:
    """Spec 2026-06-12 SS3(e) dynamic half -- any flux the gates zero at edge
    must emit a census line when it would have cleared the enlargement bar:
    correct-but-loud. The hook reads the same edge_species_rates-snapshot
    staleness as the enlargement read it audits (amendment A2 -- a feature;
    do not 'fix' onto accepted-state-only plumbing)."""

    def _b1_with_driver(self, archetype=None):
        """Gated B1 channel (A -> G, edge) + a slow pure-gas core driver
        X -> Y so char_rate > 0 (without it the base.pyx:839 singularity
        path owns the step and the hook correctly defers)."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["Y"]]
        mask = [False, False, False, False, True, True]
        gated = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        if archetype is not None:
            gated.polymer_flux_archetype = archetype
        driver = Reaction(
            reactants=[sp["X"]], products=[sp["Y"]],
            kinetics=Arrhenius(A=(1.0e-3, "1/s"), n=0.0,
                               Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)
        rs = _gate17_rs(core, mask, [driver], edge_spcs=[sp["G"]],
                        rxns_edge=[gated])
        return sp, core, driver, gated, rs

    def test_census_fires_with_payload_and_warn_once_key(self, caplog):
        """T5 (base): exactly one census line per build however many
        accepted steps; gate code B; the ungated ratio matches the
        hand-computed species ratio; a second build re-announces (warn-once
        is per engine rebuild = per RMG iteration, deliberately)."""
        import re
        sp, core, driver, gated, rs = self._b1_with_driver()
        with caplog.at_level(logging.WARNING):
            _gate17_simulate(rs, core, [driver], [sp["G"]], [gated])
        lines = _census_lines(caplog)
        assert len(lines) == 1, lines
        msg = lines[0]
        assert "gate=B" in msg
        assert "no prospectively-condensed product" in msg
        assert "tol_move_to_core" in msg
        # hand-computed t~0 ratio: ungated G rate = kf*mu1/V_poly = 10.0;
        # char_rate = sqrt(2) * (1e-3 / V_gas), V_gas = R*800/1e5.
        v_gas = constants.R * 800.0 / 1.0e5
        expected_ratio = 10.0 / (np.sqrt(2.0) * 1.0e-3 / v_gas)
        ratio = float(re.search(r"ungated_ratio=([0-9.eE+-]+)", msg).group(1))
        assert ratio == pytest.approx(expected_ratio, rel=1e-2)
        # the species carrying the ratio is named
        assert "G" in msg
        # warn-once across steps, re-announce per rebuild:
        caplog.clear()
        sp2, core2, driver2, gated2, rs2 = self._b1_with_driver()
        with caplog.at_level(logging.WARNING):
            _gate17_simulate(rs2, core2, [driver2], [sp2["G"]], [gated2])
        assert len(_census_lines(caplog)) == 1

    def test_census_payload_carries_predemotion_stamp(self, caplog):
        """T5 (stamp-divergence variant): a reaction stamped
        SCISSION_FRAGMENT(3) with unresolvable dst is demoted to
        UNRESOLVED(4) by the init pass AND Gate-B zeroed; the census payload
        shows pre-demotion 3 / post-demotion 4 -- the stamp thread proven
        end-to-end on the live path (spec SS7(iv))."""
        sp, core, driver, gated, rs = self._b1_with_driver(archetype=3)
        # liveness: the solver did demote (edge reaction index n_rxn + 0)
        assert int(np.asarray(rs.reaction_flux_archetype)[1]) == 4
        with caplog.at_level(logging.WARNING):
            _gate17_simulate(rs, core, [driver], [sp["G"]], [gated])
        lines = _census_lines(caplog)
        assert len(lines) == 1, lines
        assert "pre-demotion=3" in lines[0]
        assert "post-demotion=4" in lines[0]

    def test_census_gate_a_variant_names_decisive_product(self, caplog):
        """T5 (Gate-A variant, armed shape A): X(gas) -> G(edge) with a
        configured G pool -- prospective Gate A zeroes at edge; the census
        line carries gate=A and the prospectively-condensed product label."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["Y"], sp["G_mu0"], sp["G_mu1"], sp["G_mu2"]]
        mask = [False, False, False, False, True, True, False, False, False]
        gated = Reaction(reactants=[sp["X"]], products=[sp["G"]], **_KIN)
        driver = Reaction(
            reactants=[sp["X"]], products=[sp["Y"]],
            kinetics=Arrhenius(A=(1.0e-3, "1/s"), n=0.0,
                               Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)
        rs = _gate17_rs(core, mask, [driver], edge_spcs=[sp["G"]],
                        rxns_edge=[gated],
                        pools=(("A", (1, 2, 3)), ("G", (6, 7, 8))))
        with caplog.at_level(logging.WARNING):
            _gate17_simulate(rs, core, [driver], [sp["G"]], [gated])
        lines = _census_lines(caplog)
        assert len(lines) == 1, lines
        assert "gate=A" in lines[0]
        assert "G" in lines[0]

    def test_census_gate_a_core_product_decisive_label(self, caplog):
        """T5 (Gate-A variant, decisive product CORE-resident): X(gas) ->
        A(core pool proxy) + G(edge) -- the prospectively-condensed product
        that decides Gate A lives in CORE, so the census must resolve its
        label through core_species (threaded via the base hook) instead of
        printing an opaque core_index=N."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["Y"]]
        mask = [False, False, False, False, True, True]
        gated = Reaction(reactants=[sp["X"]], products=[sp["A"], sp["G"]],
                         **_KIN)
        driver = Reaction(
            reactants=[sp["X"]], products=[sp["Y"]],
            kinetics=Arrhenius(A=(1.0e-3, "1/s"), n=0.0,
                               Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)
        rs = _gate17_rs(core, mask, [driver], edge_spcs=[sp["G"]],
                        rxns_edge=[gated])
        with caplog.at_level(logging.WARNING):
            _gate17_simulate(rs, core, [driver], [sp["G"]], [gated])
        lines = _census_lines(caplog)
        assert len(lines) == 1, lines
        msg = lines[0]
        assert "gate=A" in msg
        # the decisive core product is named by LABEL, never by index
        assert "decisive=A;" in msg
        assert "core_index=" not in msg

    def test_census_silent_on_clean_shapes(self, caplog):
        """T6: canonical chip (probe 3), demoted-but-ungated (probe 4) and
        shape A @ HEAD config (probe 2) emit ZERO census lines through a
        real simulate(), and their edge rates keep the HEAD values
        (2.0 / 10.0 / 30.068) -- the consistent-silent-at-HEAD pin for shape
        A that item 16 will flip to the armed expectation. Silence here is
        owned by gate_code == 0 (clean shapes), NOT by a dead char_rate:
        the proxy-touching rows feed proxy_activity into core_species_rates
        (measured chip-case char_rate = 4.0); only the pure-gas shape-A row
        has char_rate == 0. The fourth case is the converse pin: a GATED B1
        row whose only flux the gate zeroes leaves char_rate == 0, and the
        hook defers silently to the base.pyx singularity path -- the
        char_rate == 0 guard is load-bearing, not decorative."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = [False, False, False, False, True]
        chip = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]],
                        **_KIN)
        chip.is_end_group_reaction = True
        chip.polymer_flux_archetype = 5
        chip.polymer_chip_units = 2
        demoted = Reaction(reactants=[sp["A"]],
                           products=[sp["R17"], sp["G"]], **_KIN)
        demoted.polymer_flux_archetype = 3  # dst unresolvable -> UNRESOLVED
        gas_a = Reaction(reactants=[sp["X"]], products=[sp["G"]], **_KIN)
        gated_b1 = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        cases = [
            (chip, core, mask, 2.0, 0),
            (demoted,
             [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
              sp["R17"]],
             [False, False, False, False, True, False], 10.0, 0),
            (gas_a, core, mask, 2.0e5 / (constants.R * 800.0), 0),
            # gated row, char_rate == 0: gate B zeroes the only flux, the
            # suppressed proxy_activity leaves core_species_rates all-zero,
            # so the hook must defer WITHOUT a census line even though the
            # gate code is live and the counterfactual is huge.
            (gated_b1, core, mask, 0.0, 2),
        ]
        for rxn, c, m, expected, expected_gate in cases:
            caplog.clear()
            rs = _gate17_rs(list(c), list(m), [], edge_spcs=[sp["G"]],
                            rxns_edge=[rxn])
            with caplog.at_level(logging.WARNING):
                _gate17_simulate(rs, list(c), [], [sp["G"]], [rxn])
            assert _census_lines(caplog) == []
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
            assert float(np.asarray(rs.edge_reaction_rates)[0]) == \
                pytest.approx(expected, rel=1e-3)
            assert int(np.asarray(rs.edge_reaction_gate_code)[0]) == \
                expected_gate
            if expected_gate:
                # liveness for the deference pin: gated, counterfactual
                # alive, and char_rate genuinely 0 (no core flux at all).
                assert float(np.asarray(
                    rs.edge_reaction_rates_ungated)[0]) > 0.0
                assert not np.any(np.asarray(rs.core_species_rates))

    def test_counterfactual_purity(self):
        """T9: R2 observes, never leaks into state. On the gated B1 fixture
        the consistency point (edge_reaction_rates == 0, no
        edge_species_rates contribution) holds, the ungated arrays carry
        exactly the counterfactual, and the parent proxy's
        core_species_rates (proxy_activity-derived) excludes the gated
        |flux| -- ghost flux must not feed spawn/similarity diagnostics."""
        sp, core, driver, gated, rs = self._b1_with_driver()
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert float(np.asarray(rs.edge_reaction_rates)[0]) == 0.0
        assert float(np.asarray(rs.edge_species_rates)[0]) == 0.0
        assert int(np.asarray(rs.edge_reaction_gate_code)[0]) == 2
        assert float(np.asarray(rs.edge_reaction_rates_ungated)[0]) == \
            pytest.approx(10.0)
        assert float(np.asarray(rs.edge_species_rates_ungated)[0]) == \
            pytest.approx(10.0)
        # ghost-flux suppression: proxy A (index 0) carries NO activity from
        # the gated edge row (residual section 11 reports
        # proxy_activity/V_poly for proxies).
        assert float(np.asarray(rs.core_species_rates)[0]) == 0.0

    def test_counterfactual_purity_r1_site_bimolecular(self):
        """T9 (r1 site): X(gas) + A(pool proxy) -> Y + G, all products
        prospectively gas, is Gate-B zeroed at edge. The proxy sits at the
        SECOND reactant slot, so this pins the r1-site proxy_activity
        suppression in residual section 3 -- the r0-site test above never
        reaches it. A's core_species_rates (proxy_activity-derived) must
        stay 0 while the counterfactual proves the row alive."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["Y"]]
        mask = [False, False, False, False, True, True]
        gated = Reaction(reactants=[sp["X"], sp["A"]],
                         products=[sp["Y"], sp["G"]], **_KIN)
        rs = _gate17_rs(core, mask, [], edge_spcs=[sp["G"]],
                        rxns_edge=[gated])
        # fixture-shape pins: proxy A really is r1, gate really is B
        assert int(np.asarray(rs.reactant_indices)[0, 1]) == 0
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert int(np.asarray(rs.edge_reaction_gate_code)[0]) == 2
        assert float(np.asarray(rs.edge_reaction_rates)[0]) == 0.0
        assert float(np.asarray(rs.edge_species_rates)[0]) == 0.0
        # liveness: the suppression is judged on a LIVE counterfactual,
        # not a dead row (kf * C_X * mu1-site > 0).
        assert float(np.asarray(rs.edge_reaction_rates_ungated)[0]) > 0.0
        # r1-site ghost-flux suppression: no spawn/similarity activity on A.
        assert float(np.asarray(rs.core_species_rates)[0]) == 0.0

    def test_counterfactual_purity_product_site_core_proxy(self):
        """T9 (product site): X(gas) -> A(CORE pool proxy) + G(edge) is
        Gate-A zeroed at edge (gas event, prospectively-condensed core
        product). The proxy appears ONLY as a product, so this pins the
        p_idx_tmp-site proxy_activity suppression in residual section 4."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = [False, False, False, False, True]
        gated = Reaction(reactants=[sp["X"]], products=[sp["A"], sp["G"]],
                         **_KIN)
        rs = _gate17_rs(core, mask, [], edge_spcs=[sp["G"]],
                        rxns_edge=[gated])
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert int(np.asarray(rs.edge_reaction_gate_code)[0]) == 1
        assert float(np.asarray(rs.edge_reaction_rates)[0]) == 0.0
        assert float(np.asarray(rs.edge_species_rates)[0]) == 0.0
        # liveness: counterfactual = kf * C_X > 0 (pure-gas rate law).
        assert float(np.asarray(rs.edge_reaction_rates_ungated)[0]) > 0.0
        # product-site ghost-flux suppression: no activity lands on A.
        assert float(np.asarray(rs.core_species_rates)[0]) == 0.0


class TestPhaseGateStaticCensus:
    """Spec 2026-06-12 SS3(e) STATIC half (amendment A1): a reaction whose
    species each reached core on OTHER channels' flux sits core-resident and
    gate-zeroed without ever crossing the edge on its own flux (the third
    route -- proven live by the 13 EPDM Gate-A recombinations). Gate
    verdicts for core reactions are STATIC (masks + event type, no rates),
    so the census enumerates them once per initialize_model, zero residual
    cost."""

    def test_static_census_announces_core_gated_reaction_per_rebuild(
            self, caplog):
        """T11: B2-shaped reaction placed directly in core_reactions with
        every participant in core (the third-route end state by
        construction). Liveness pins: kinetics alive (kf > 0) and the gate
        verdict genuinely zero under the masks -- the census line can only
        mean the init-time enumeration ran."""
        sp = _gate17_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["G"]]
        mask = [False, False, False, False, True, True]
        with caplog.at_level(logging.WARNING):
            rs = _gate17_rs(core, mask, [rxn])
        lines = _census_lines(caplog)
        assert len(lines) == 1, lines
        msg = lines[0]
        assert "gate=B" in msg
        assert "static (core, init-time)" in msg
        assert "no prospectively-condensed product" in msg
        # pre-demotion NONE(0) -> post-demotion UNRESOLVED(4): the proxy-
        # touching unstamped reaction is remapped by the init pass; the
        # static payload carries both ends of the stamp thread.
        assert "pre-demotion=0" in msg
        assert "post-demotion=4" in msg
        # liveness + ownership: kinetics alive, the GATE owns the zero
        assert float(rs.kf[0]) > 0.0
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert float(np.asarray(rs.core_reaction_rates)[0]) == 0.0
        # once per rebuild: a SECOND initialize_model re-announces; the same
        # build never repeats (the enumeration runs exactly once per init).
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            rs.initialize_model(list(core), [rxn], [], [])
        assert len(_census_lines(caplog)) == 1


# ---------------------------------------------------------------------------
# Item 17 Task 5A (spec 2026-06-12 A5-2/A5-3): live-edge stage-1 + R1-EDGE
# edge-suffix provenance guard. T12 (production fallback RAISES) is the
# live-path tripwire confirmed RED on HEAD f13e4ce52 (where the guard does
# not exist and the production build silently falls back); T13 is the armed-
# rows fail-under-forced-fallback mutation proof. Skeletons per plan 5A.1 --
# the `...` halves are filled by the implementer against the 5A.2 knob.
# ---------------------------------------------------------------------------


class TestProspectiveEdgeProvenance:
    """Spec 2026-06-12 SS3(d) R1-EDGE (A5-2): R1's core-prefix check CANNOT
    see the staging defect -- under the edge-defaults-GAS fallback the core
    PREFIX still matches gas_species_mask exactly (the fallback copies it
    verbatim), so R1 passes green while the edge SUFFIX is silently
    default-filled. A separate guard fires on the part R1 is structurally
    blind to: a PRODUCTION build (not flagged as a legitimate test-fallback)
    whose prospective edge classification came from defaults must RAISE."""

    def _prod_build_kwargs(self, core, mask, cfg):
        # A "production-shaped" build is one NOT flagged
        # allow_default_prospective_edge=True (i.e. built from a blueprint
        # phase handle, the live path). Tests express that by passing a phase
        # classifier handle (set below) and leaving the flag default-False.
        # initial_mole_fractions is a MANDATORY positional on the HEAD
        # constructor (polymer.pyx:391) -- seed the first gas-masked core
        # species, exactly as the landed _gate17_rs fixture does, so the build
        # reaches initialize_model (the provenance raise site) instead of
        # dying at __init__.
        seed_idx = int(np.where(np.asarray(mask, dtype=bool))[0][0])
        return dict(
            T=800.0, P=1.0e5, initial_mole_fractions={core[seed_idx]: 1.0},
            V_poly=1.0, polymer_pools=[cfg],
            mass_transfer=[], gas_species_mask=mask.copy(),
            constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[])

    def test_production_fallback_raises_provenance(self):
        """T12: a production-shaped build (live-edge classifier handle present,
        allow_default_prospective_edge default-False) that is forced to
        default-fill the edge suffix -- the exact shape that fired
        PROSPECTIVE-MASK SEED STALE on the live Task-6 run -- RAISES with the
        PROSPECTIVE-MASK PROVENANCE: sentinel naming the default-filled edge
        count. Liveness pin before the raise: the edge list is genuinely
        non-empty AND carries a prospectively-condensed edge species (the G
        product under a configured G pool), so the raise can only mean 'the
        production path consumed defaults', never 'empty edge / nothing to
        classify'."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["G_mu0"], sp["G_mu1"], sp["G_mu2"]]
        mask = np.array([False, False, False, False, True,
                         False, False, False], dtype=bool)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        # Build with a CLASSIFIER HANDLE present (production path), but force
        # the fallback (no live-edge stage-1 result available) -- e.g. by
        # constructing with a stale/absent seed AND no live classifier wired,
        # the implementer's chosen "force fallback on a production build" knob
        # (see 5A.2). The edge carries G (prospectively condensed via the G
        # pool config in stage 2), so default-filling it is a real provenance
        # violation.
        # Liveness: edge non-empty, G prospectively condensed once classified.
        rs = HybridPolymerSystem(
            **self._prod_build_kwargs(core, mask, cfg),
            # production marker: NOT allow_default_prospective_edge
        )
        # the implementer wires the production classifier handle + the forced-
        # fallback knob per 5A.2; with it absent on a production build the
        # guard must raise.
        with pytest.raises(ValueError, match="PROSPECTIVE-MASK PROVENANCE:"):
            rs.initialize_model(list(core), [], [sp["G"]], [rxn])

    def test_live_edge_stage1_does_not_raise(self):
        """T12 green counterpart: a production build whose prospective edge
        suffix came from the LIVE-EDGE stage-1 classifier (provenance all
        stage-1-classified) does NOT raise -- the guard fires only on
        default-fill, not on a legitimately-classified edge."""
        # Production live-edge path: a classifier handle is wired (the bound-
        # method stand-in _stage1_classifier, classifying G condensed by label
        # exactly as polymerPhase.get_gas_mask classifies a registered pool
        # member), allow_default_prospective_edge default-False. The solver
        # re-runs stage 1 over the live chain(core, edge) at initialize_model.
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = np.array([False, False, False, False, True], dtype=bool)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        rs = HybridPolymerSystem(
            **self._prod_build_kwargs(core, mask, cfg),
            prospective_classifier=_stage1_classifier)  # production marker
        # No raise: the edge suffix was classified by the live-edge stage-1.
        rs.initialize_model(list(core), [], [sp["G"]], [rxn])
        n_core = rs.num_core_species
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        # provenance all stage-1-classified, no default-fill
        assert np.all(np.asarray(rs._prospective_edge_provenance) == 1)
        # the edge suffix equals a fresh classifier call over chain(core, edge)
        fresh = np.asarray(_stage1_classifier(list(core) + [sp["G"]]),
                           dtype=bool)
        assert np.array_equal(pm[n_core:], fresh[n_core:])
        assert bool(pm[-1]) is False  # G classified CONDENSED by stage-1

    def test_flagged_test_fallback_does_not_raise(self):
        """T12 green counterpart: a build FLAGGED
        allow_default_prospective_edge=True (direct-test/runner construction
        with no blueprint phase object -- e.g. polymer_moments_runner) does
        NOT raise even with default-filled edges. This is the legitimate
        last-resort fallback (spec SS3(a)/SS3(d))."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = np.array([False, False, False, False, True], dtype=bool)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg],
            mass_transfer=[], gas_species_mask=mask.copy(),
            constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            allow_default_prospective_edge=True)
        # No raise: the flag marks this build as legitimately fallback-permitted.
        rs.initialize_model(list(core), [], [sp["G"]], [rxn])
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert pm.shape[0] == rs.num_core_species + 1
        assert bool(pm[-1]) is True  # edge G default GAS, permitted under flag


class TestArmedRowsExerciseStage1:
    """Spec 2026-06-12 SS5/A5-3 acceptance gate: the two ARMED parity rows
    (B1_configured, A_armed) PASS with the real A5-2 live-edge stage-1 mask AND
    go RED under forced fallback (edge-defaults-GAS). Without this proof the
    armed rows could pass vacuously (parity-via-zero by accident) the way the
    live EPDM deck 'passed' by luck.

    Design note (why NOT TestUmbrellaPhaseGateParity._build_pair): there the
    edge G is condensed by a STAGE-2 pool-label override (a configured "G"
    pool), which runs identically whether stage-1 fell back or not -- so
    forcing the fallback would NOT flip G and the mutation proof would be
    vacuous. These fixtures instead source G's condensed verdict from STAGE 1
    (the wired classifier), with NO G pool, so wiring the classifier to None
    genuinely flips edge G to GAS. The post-promotion core build condenses G
    via a direct mask bit (the verdict item 16's promotion would record),
    matching the EXPECTED armed values (B1_configured 10.0, A_armed 0.0)."""

    def _build_armed_pair(self, case, live):
        """Build the (edge, core) pair for an armed case. ``live`` selects the
        EDGE build's stage-1 source: True wires the live-edge classifier
        (production path, G classified condensed by stage 1); False forces the
        edge-defaults-GAS fallback (classifier None) -- flagged permitted so
        the build itself does not raise, isolating the parity break to the
        lost gate verdict, not a provenance raise."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
        mask = [False, False, False, False, True]
        if case == "B1_configured":  # poly event A -> G
            rxn = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        else:  # A_armed: gas event X -> G
            rxn = Reaction(reactants=[sp["X"]], products=[sp["G"]], **_KIN)
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        rs_edge = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=np.array(mask, dtype=bool),
            constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_classifier=(_stage1_classifier if live else None),
            allow_default_prospective_edge=(not live))
        rs_edge.initialize_model(list(core), [], [sp["G"]], [rxn])
        # Post-promotion core build: G promoted to core CONDENSED (the
        # prospective verdict item 16 records). No G pool; the mask bit carries
        # the verdict. allow_default_prospective_edge via _gate17_rs.
        rs_core = _gate17_rs(core + [sp["G"]], mask + [False], [rxn])
        return rs_edge, rs_core

    @pytest.mark.parametrize("case", ["B1_configured", "A_armed"])
    def test_armed_row_red_under_forced_fallback(self, case):
        """Force the edge-defaults-GAS fallback on the armed parity pair: the
        edge G, classified condensed by stage 1 on the live path, instead
        defaults GAS, so Gate B (B1_configured) / Gate A (A_armed) never fires
        at edge and the umbrella parity (edge rate == post-promotion core rate)
        BREAKS.

        B1_configured (poly event A -> G): in the post-promotion core build G
        is condensed, so Gate B PASSES and the core rate is the armed flux
        (10.0). Under the forced fallback the edge G defaults GAS, so Gate B
        ZEROES the edge rate (0.0). 10.0 != 0.0 -- parity broken.

        A_armed (gas event X -> G): in the core build G is condensed, so Gate A
        ZEROES the core rate (0.0). Under the forced fallback the edge G
        defaults GAS, so Gate A does NOT fire and the edge rate is the ungated
        30.068... 30.068 != 0.0 -- parity broken."""
        rs_edge, rs_core = self._build_armed_pair(case, live=False)
        # liveness: the forced-fallback edge G really did default GAS (stage-1
        # would have condensed it), so the break is the lost gate verdict.
        pm = np.asarray(rs_edge.prospective_gas_mask, dtype=bool)
        assert bool(pm[-1]) is True  # edge G default GAS (fallback)
        rs_edge.residual(0.0, rs_edge.y, np.zeros_like(rs_edge.y))
        rs_core.residual(0.0, rs_core.y, np.zeros_like(rs_core.y))
        edge_rate = float(np.asarray(rs_edge.edge_reaction_rates)[0])
        core_rate = float(np.asarray(rs_core.core_reaction_rates)[0])
        # THE PROOF: under forced fallback the armed flux verdict is lost ->
        # parity is BROKEN (the vacuous-pass guard).
        assert edge_rate != pytest.approx(core_rate, abs=1e-9), (
            f"armed row {case} did NOT break under forced fallback "
            f"(edge {edge_rate} vs core {core_rate}) -- the mutation proof is "
            f"vacuous; G's condensed verdict is not stage-1-sourced")

    @pytest.mark.parametrize("case", ["B1_configured", "A_armed"])
    def test_armed_row_green_with_live_stage1(self, case):
        """Counterpart: with the real A5-2 live-edge stage-1 mask (classifier
        wired), the armed rows PASS -- edge rate == post-promotion core rate ==
        the expected armed value (B1_configured 10.0, A_armed 0.0). Closes the
        green-RED-green cycle in ONE class and proves the verdict is genuinely
        stage-1-sourced (provenance all 1, no default-fill)."""
        rs_edge, rs_core = self._build_armed_pair(case, live=True)
        assert np.all(np.asarray(rs_edge._prospective_edge_provenance) == 1)
        pm = np.asarray(rs_edge.prospective_gas_mask, dtype=bool)
        assert bool(pm[-1]) is False  # edge G classified CONDENSED by stage-1
        rs_edge.residual(0.0, rs_edge.y, np.zeros_like(rs_edge.y))
        rs_core.residual(0.0, rs_core.y, np.zeros_like(rs_core.y))
        edge_rate = float(np.asarray(rs_edge.edge_reaction_rates)[0])
        core_rate = float(np.asarray(rs_core.core_reaction_rates)[0])
        assert edge_rate == pytest.approx(core_rate, abs=1e-12)
        assert edge_rate == pytest.approx(
            TestUmbrellaPhaseGateParity.EXPECTED[case], rel=1e-9, abs=1e-12)


class TestLiveEdgeRebuildWidensR1RaiseSurface:
    """Spec 2026-06-12 §3(d) R1 raise-surface change (A5-review refinement):
    the A5-2 live-edge stage-1 rebuild WIDENS R1's core-prefix tripwire surface
    to label-fallback + edge-collision decks, where the OLD stale path (which
    copied gas_species_mask verbatim into the prospective prefix) could never
    raise. Intended loudness; EPDM-clean (0 tripwire lines on the live Task-6
    run -- every poly member is id-registered, no label-fallback core members,
    no edge label collisions).

    The reviewer's divergence vector: a CORE poly member condensed ONLY via
    label fallback (id mismatch -- a copied species whose id() is not in
    get_gas_mask's registry but whose label is) PLUS a duplicate of that label
    in the EDGE. get_gas_mask's label_fallback_safe (polymer_input.py:634-637)
    is computed PER-LIST, so the edge duplicate disables the fallback for the
    combined chain(core, edge) call: the core member that condensed in the
    core-only call (gas_species_mask) flips GAS in the combined prefix, the
    prospective core prefix diverges, and R1 RAISES. The old stale path was
    structurally incapable of this -- it copied gas_species_mask verbatim."""

    def _phase_and_pieces(self):
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import PolymerPhase, PolymerPool
        mono = _spc("CCCC", "MON")        # the registered pool monomer
        mono_copy = _spc("CCCC", "MON")   # SAME label, DIFFERENT id -> label
        #                                   fallback only (id not in registry)
        mu0 = _spc("CO", "A_mu0")
        mu1 = _spc("C=O", "A_mu1")
        mu2 = _spc("C#N", "A_mu2")
        x = _spc("N#N", "X")              # inert gas seed
        edge_dup = _spc("C", "MON")       # duplicate label ONLY in the edge
        g = _spc("[CH3]", "G")
        pool_in = PolymerPool(label="A", xs=2, monomer=mono, explicit_map={},
                              mu_species=[mu0, mu1, mu2])
        phase = PolymerPhase(density=Quantity(900.0, "kg/m^3"),
                             initial_moments={"A": (1.0, 5.0, 30.0)},
                             initial_explicit={}, pools=[pool_in])
        core = [mono_copy, mu0, mu1, mu2, x]
        return phase, core, edge_dup, g

    def test_get_gas_mask_diverges_per_list_under_collision(self):
        """The narrow invariant the reviewer probed: get_gas_mask(core) and
        get_gas_mask(chain(core, edge))[:n_core] DIVERGE under the collision --
        the core member that condensed via label fallback in the core-only call
        flips GAS in the combined call (fallback disabled by the edge dup)."""
        phase, core, edge_dup, g = self._phase_and_pieces()
        core_only = np.asarray(phase.get_gas_mask(core), dtype=bool)
        combined = np.asarray(phase.get_gas_mask(core + [edge_dup, g]),
                              dtype=bool)
        n_core = len(core)
        # core-only: index 0 (MON copy) condensed via label fallback
        assert bool(core_only[0]) is False
        # combined: edge dup disables the fallback -> index 0 flips GAS
        assert bool(combined[0]) is True
        assert not np.array_equal(core_only, combined[:n_core])

    def test_r1_raises_on_label_fallback_plus_edge_collision(self):
        """The live-edge rebuild makes R1 capable of raising on this deck: the
        production prospective_classifier (the real bound phase.get_gas_mask)
        recomputes stage-1 over chain(core, edge), the core prefix flips at MON,
        and R1 raises PROSPECTIVE-MASK TRIPWIRE: naming the diverging member.
        Liveness: gas_species_mask is the genuine core-only get_gas_mask (MON
        condensed), so the raise is the lost-condensation divergence, not a
        doctored seed (the way T4 injects one)."""
        phase, core, edge_dup, g = self._phase_and_pieces()
        gas_mask = np.asarray(phase.get_gas_mask(core), dtype=bool)
        assert bool(gas_mask[0]) is False  # MON condensed in the static mask
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        rxn = Reaction(reactants=[core[0]], products=[g], **_KIN)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=gas_mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            # production marker: real bound classifier, fallback NOT permitted
            prospective_classifier=phase.get_gas_mask)
        with pytest.raises(ValueError, match="PROSPECTIVE-MASK TRIPWIRE:"):
            rs.initialize_model(list(core), [], [edge_dup, g], [rxn])
        # the diverging member (MON) is NAMED in the raise
        with pytest.raises(ValueError, match=r"MON"):
            rs.initialize_model(list(core), [], [edge_dup, g], [rxn])


class TestProspectiveClassifierLengthGuard:
    """Spec 2026-06-12 §3(a): the live-edge stage-1 rebuild defensively checks
    that the plumbed prospective_classifier returns an array sized to
    chain(core, edge). A classifier returning the wrong length RAISES
    PROSPECTIVE-MASK CLASSIFIER: -- and production get_gas_mask never trips it
    (it returns len(input) by construction). Both halves pinned so the guard is
    proven live AND proven free of false-positives on the real classifier."""

    def _core_pieces(self):
        a = _spc("CCCC", "A")
        mu0 = _spc("CO", "A_mu0")
        mu1 = _spc("C=O", "A_mu1")
        mu2 = _spc("C#N", "A_mu2")
        x = _spc("N#N", "X")
        g = _spc("[CH3]", "G")
        core = [a, mu0, mu1, mu2, x]
        mask = np.array([False, False, False, False, True], dtype=bool)
        cfg = PolymerPoolConfig(label="A", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None)
        rxn = Reaction(reactants=[a], products=[g], **_KIN)
        return core, mask, cfg, rxn, g

    def test_wrong_length_classifier_raises(self):
        """A classifier returning len(input) - 1 (a wrong-length array) makes
        the live-edge rebuild RAISE PROSPECTIVE-MASK CLASSIFIER:, not silently
        mis-size the prospective mask."""
        core, mask, cfg, rxn, g = self._core_pieces()

        def bad_classifier(species_list):
            return np.ones(len(species_list) - 1, dtype=bool)

        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_classifier=bad_classifier)
        with pytest.raises(ValueError, match="PROSPECTIVE-MASK CLASSIFIER:"):
            rs.initialize_model(list(core), [], [g], [rxn])

    def test_production_get_gas_mask_returns_input_length_no_false_positive(
            self):
        """The real production classifier (bound phase.get_gas_mask) returns an
        array of len(input) for any chain(core, edge), so the length guard is a
        no-op on the production path -- it never false-positives. Pinned via the
        classifier directly AND through a full build that does NOT raise the
        length guard."""
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import PolymerPhase, PolymerPool
        core, mask, cfg, rxn, g = self._core_pieces()
        pool_in = PolymerPool(label="A", xs=2, monomer=core[0],
                              explicit_map={},
                              mu_species=[core[1], core[2], core[3]])
        phase = PolymerPhase(density=Quantity(900.0, "kg/m^3"),
                             initial_moments={"A": (1.0, 5.0, 30.0)},
                             initial_explicit={}, pools=[pool_in])
        combined_input = list(core) + [g]
        out = np.asarray(phase.get_gas_mask(combined_input))
        assert out.shape[0] == len(combined_input)  # exactly len(input)
        # full build through the live-edge rebuild: no CLASSIFIER length raise
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=np.asarray(phase.get_gas_mask(core), dtype=bool),
            constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_classifier=phase.get_gas_mask)
        rs.initialize_model(list(core), [], [g], [rxn])
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert pm.shape[0] == rs.num_core_species + 1
