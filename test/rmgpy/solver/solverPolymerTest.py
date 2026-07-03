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


def _qssa_triplet(A=1.0e13, n=0.0, Ea=1.0e5):
    """Arrhenius triplet for the radical_qssa_unzip channel (SI convention:
    A [s^-1] for initiation/depropagation AND for transfer -- transfer is
    pseudo-first-order as implemented, ktr multiplies R directly -- and
    [m^3 mol^-1 s^-1] for the bimolecular termination only; Ea [J/mol])."""
    return dict(A=A, n=n, Ea=Ea)


def _qssa_channel(**overrides):
    """A valid minimal radical_qssa_unzip channel config (mandatory blocks
    only; efficiency/monomer_yield/basis/transfer left to their defaults)."""
    ch = dict(
        initiation=_qssa_triplet(A=1.0e15, Ea=3.0e5),
        depropagation=_qssa_triplet(A=1.0e13, Ea=8.0e4),
        termination=_qssa_triplet(A=1.0e8, Ea=1.0e4),
    )
    ch.update(overrides)
    return ch


def _exact_triplet(A):
    """Arrhenius triplet with n=0, Ea=0 so k(T) == A EXACTLY (T**0.0 == 1.0,
    exp(-0.0) == 1.0 bitwise): lets weak-link pins assert bitwise/analytic
    equalities instead of tolerance soup."""
    return dict(A=float(A), n=0.0, Ea=0.0)


def _weaklink_channel(**overrides):
    """A valid weak-link allyl/U-state channel (schema-2.2 vocabulary): the
    four weak-link keys present, legacy summed 'termination' absent."""
    ch = dict(
        initiation=_qssa_triplet(A=1.0e15, Ea=3.0e5),
        depropagation=_qssa_triplet(A=1.0e13, Ea=8.0e4),
        initiation_allyl=_qssa_triplet(A=2.0e14, Ea=2.4e5),
        termination_recombination=_qssa_triplet(A=6.0e7, Ea=8.0e3),
        termination_disproportionation=_qssa_triplet(A=4.0e7, Ea=1.2e4),
        unsaturated_tail_ends_initial=0.0,
    )
    ch.update(overrides)
    return ch


def _weaklink_exact_channel(ki_A=1.0e-3, kdp_A=1.0e2, kia_A=2.0e4,
                            ktrec_A=5.0e7, ktdisp_A=5.0e7, u0=0.0):
    """Weak-link channel built entirely from n=0/Ea=0 triplets so every
    k(T) == A exactly (see _exact_triplet)."""
    return dict(
        initiation=_exact_triplet(ki_A),
        depropagation=_exact_triplet(kdp_A),
        initiation_allyl=_exact_triplet(kia_A),
        termination_recombination=_exact_triplet(ktrec_A),
        termination_disproportionation=_exact_triplet(ktdisp_A),
        unsaturated_tail_ends_initial=float(u0),
    )


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

        # Released-monomer target for the lumped unzip channel. k_unzip > 0
        # requires a monomer_poly_index (solver invariant: the drain is
        # unconditional, the emission is gated on the index) -- a SEPARATE
        # species from P2 so the handshake assertion dn_dt[P2] > 0 keeps
        # proving handshake flux, not monomer emission. GAS since the
        # 2026-07-03 monomer-gas fix (release = devolatilization).
        M = _spc("C", "M")

        core_species = [Inert, P2, Mu0, Mu1, Mu2, M]

        # gas mask: Inert and the released-monomer target M are gas
        gas_species_mask = np.array([True, False, False, False, False, True],
                                    dtype=bool)

        pool = PolymerPoolConfig(
            label="poly",
            xs=2,
            explicit_dp_to_species_index={2: 1},  # DP=2 -> P2 index
            mu_indices=(2, 3, 4),                 # Mu0, Mu1, Mu2 indices
            monomer_poly_index=5,                 # unzip releases into M
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
        # GAS released-monomer target (2026-07-03 fix): k_unzip > 0 requires
        # a monomer_poly_index (solver invariant). The emission k_unzip*mu0
        # is part of the finiteness surface under test.
        M = _spc("C", "M")
        core_species = [Inert, Mu0, Mu1, Mu2, M]
        gas_species_mask = np.array([True, False, False, False, True],
                                    dtype=bool)

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
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
        assert solver_mod.FLUX_VOLATILE_EJECTION == int(PolymerFluxArchetype.VOLATILE_EJECTION) == 6

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
        # NOTE (item 18): this invariant is MOLAR, not MW-weighted, so it cannot
        # see the FEATURE-radical mass leak (pool debited 1 monomer-unit vs gas
        # gaining full chain MW). See
        # test_refused_feature_abstraction_conserves_backbone_mass for the
        # MW-weighted conservation gate.
        assert np.isclose(dn_dt[2] + dn_dt[4], 0.0, atol=1e-12)

    def test_refused_feature_abstraction_conserves_backbone_mass(self):
        import numpy as np
        from rmgpy.reaction import Reaction
        from rmgpy.kinetics import Arrhenius
        from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
        Proxy = _spc("CCC(C)CCCC(C)CCCC(C)C", "epdm")
        Mu0 = _spc("CO", "epdm_mu0"); Mu1 = _spc("C=O", "epdm_mu1"); Mu2 = _spc("C#N", "epdm_mu2")
        Macro = _spc("CCC(C)CCC[C](C)CCCC(C)C", "epdm_macroradical")
        H = _spc("[H]", "H"); H2 = _spc("[H][H]", "H2")
        macro_mw = Macro.molecule[0].get_molecular_weight() * 1000.0
        core = [Proxy, Mu0, Mu1, Mu2, Macro, H, H2]
        mask = np.array([False, False, False, False, False, True, True], dtype=bool)  # macro condensed -> past Gate B
        mono_mw = Proxy.molecule[0].get_molecular_weight() * 1000.0 / 3
        rxn = Reaction(reactants=[Proxy, H], products=[Macro, H2],
                       kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
                       reversible=False)
        rxn.polymer_flux_archetype = 4          # UNRESOLVED
        rxn.polymer_refused = True              # set by Task 3 in production; direct here to isolate the suppression
        pool = PolymerPoolConfig(label="epdm", xs=2, explicit_dp_to_species_index={},
                                 mu_indices=(1, 2, 3), monomer_poly_index=None,
                                 k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
                                 monomer_mw_g_mol=mono_mw)
        rs = HybridPolymerSystem(T=800.0, P=1.0e5, initial_mole_fractions={H: 1.0}, V_poly=1.0,
                                 polymer_pools=[pool], mass_transfer=[], gas_species_mask=mask.copy(),
                                 constant_gas_volume=False, initial_polymer_moments={"epdm": (1.0, 50.0, 3000.0)},
                                 termination=[])
        rs.initialize_model(core, [rxn], [], [])
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        i = {s.label: k for k, s in enumerate(core)}
        net = dn[i["epdm_mu1"]] * mono_mw + dn[i["epdm_macroradical"]] * macro_mw
        assert np.isclose(net, 0.0, atol=1e-6), f"backbone mass not conserved: {net:.4e} g/s fabricated"

    def test_refused_census_distinguishes_eliminating_and_accumulating(self):
        import numpy as np
        from rmgpy.reaction import Reaction
        from rmgpy.kinetics import Arrhenius
        from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
        Proxy = _spc("CCC(C)CCCC(C)CCCC(C)C", "epdm")
        Mu0 = _spc("CO", "epdm_mu0"); Mu1 = _spc("C=O", "epdm_mu1"); Mu2 = _spc("C#N", "epdm_mu2")
        SatRad = _spc("CCC(C)CCC[C](C)CCCC(C)C", "sat_macro")
        AllylRad = _spc("CCC(C)CCCC=C[C](C)CCC(C)CC", "allyl_macro")
        H = _spc("[H]", "H"); H2 = _spc("[H][H]", "H2")
        core = [Proxy, Mu0, Mu1, Mu2, SatRad, AllylRad, H, H2]
        mask = np.array([False, False, False, False, False, False, True, True], dtype=bool)
        kin = Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K"))
        r1 = Reaction(reactants=[Proxy, H], products=[SatRad, H2], kinetics=kin, reversible=False)
        r1.polymer_flux_archetype = 4; r1.polymer_refused = True; r1.polymer_refused_accumulating = False
        r2 = Reaction(reactants=[Proxy, H], products=[AllylRad, H2], kinetics=kin, reversible=False)
        r2.polymer_flux_archetype = 4; r2.polymer_refused = True; r2.polymer_refused_accumulating = True
        pool = PolymerPoolConfig(label="epdm", xs=2, explicit_dp_to_species_index={}, mu_indices=(1, 2, 3),
                                 monomer_poly_index=None, k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
        rs = HybridPolymerSystem(T=800.0, P=1.0e5, initial_mole_fractions={H: 1.0}, V_poly=1.0,
                                 polymer_pools=[pool], mass_transfer=[], gas_species_mask=mask.copy(),
                                 constant_gas_volume=False, initial_polymer_moments={"epdm": (1.0, 50.0, 3000.0)},
                                 termination=[])
        rs.initialize_model(core, [r1, r2], [], [])
        census = rs.refused_reaction_census
        assert len(census) == 2
        assert any(c["radical_class"] == "eliminating" and c["reason"] == "conduit-deferred" for c in census)
        assert any(c["radical_class"] == "accumulating" and c["reason"] == "qssa-invalid" for c in census)

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

    def test_double_count_tripwire_fires_on_copresence_and_dormant_otherwise(self):
        """
        item 18 (T2): a pool carrying BOTH an explicit, surviving beta-scission/
        chip reaction sourced from it AND a nonzero phenomenological k_scission/
        k_unzip double-counts chain degradation. initialize_model must census it
        (correct-but-loud, warn-once; severity calibrated by item 19) and must be
        DORMANT when the phenomenological constant is present but no surviving
        explicit scission/chip is sourced from the pool.

        DISCRETE_CHIP (archetype 5) is used for the explicit reaction because it
        needs only src resolved (not dst) and so survives the demotion loop when
        its reactant maps to the pool.
        """
        import rmgpy.solver.polymer as sp

        def _build(pool_a_k_scission, with_chip):
            spx, core, mask = _two_pool_species()
            # Pool A's proxy species is matched by base_label == pool.label, so
            # rename both proxy and pool to 'epdm' for them to bind.
            spx["A"].label = "epdm"
            # Chip reaction: reactant sp["A"] resolves to pool A (src == 0).
            # Complement folds back to A; chip ejects gas G. DISCRETE_CHIP needs
            # only src, so this survives demotion.
            rxn = Reaction(reactants=[spx["A"]], products=[spx["A"], spx["G"]], **_KIN)
            rxn.polymer_flux_archetype = 5   # DISCRETE_CHIP
            rxn.polymer_chip_units = 1
            pool_a = PolymerPoolConfig(
                label="epdm", xs=2, explicit_dp_to_species_index={},
                mu_indices=(1, 2, 3), monomer_poly_index=None,
                k_scission=pool_a_k_scission, k_unzip=0.0, tail_kinetics=None)
            pool_b = PolymerPoolConfig(
                label="B", xs=2, explicit_dp_to_species_index={},
                mu_indices=(5, 6, 7), monomer_poly_index=None,
                k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
            rxns = [rxn] if with_chip else []
            rs = HybridPolymerSystem(
                T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0}, V_poly=1.0,
                polymer_pools=[pool_a, pool_b], mass_transfer=[],
                gas_species_mask=mask.copy(), constant_gas_volume=False,
                initial_polymer_moments={"epdm": (1.0, 5.0, 30.0), "B": (1.0, 5.0, 30.0)},
                termination=[])
            rs.initialize_model(core, rxns, [], [])
            return rs

        # CO-PRESENCE: pool 'epdm' has k_scission=1.0 AND a surviving DISCRETE_CHIP
        rs_copresent = _build(1.0, with_chip=True)
        # Verify the chip survived demotion (still archetype 5, src == pool epdm).
        assert rs_copresent.reaction_flux_archetype[0] == sp.FLUX_DISCRETE_CHIP
        assert rs_copresent.reaction_src_pool[0] == 0   # pool 'epdm' index
        assert len(rs_copresent.double_count_census) >= 1
        assert any(c["pool"] == "epdm" for c in rs_copresent.double_count_census)

        # DORMANT: same pool k_scission=1.0 but NO explicit scission/chip reaction.
        rs_clean = _build(1.0, with_chip=False)
        assert rs_clean.double_count_census == []

    def test_double_count_tripwire_dormant_when_scission_reaction_demoted(self):
        """
        item 18 (T6 hardening): the EPDM-faithful complement to the co-presence
        case. A pool ('epdm') carries k_scission=1.0 AND an explicit DISCRETE_CHIP
        reaction IS present -- but the chip's reactant does NOT resolve to the pool
        (src == -1, an unmapped gas reactant), so the src/dst-resolution loop
        DEMOTES it to UNRESOLVED (archetype 4) BEFORE the census is built. The
        demoted reaction's pool index therefore never enters scission_src_pools,
        so the tripwire must stay DORMANT even though k_scission > 0.

        This pins the load-bearing post-demotion placement of the census against
        regression: the sibling co-presence test proves "survived demotion ->
        fires"; this proves the complementary "demoted -> does NOT fire". (If the
        census were built BEFORE the demotion loop, this would falsely fire.)
        """
        import rmgpy.solver.polymer as sp

        spx, core, mask = _two_pool_species()
        # Rename pool A's proxy + config to 'epdm' (proxy binds by base_label).
        spx["A"].label = "epdm"
        # Chip reaction whose reactant is the UNMAPPED gas 'G' (mask True, no
        # pool config): species_to_pool_indices -> -1, so src stays -1 and the
        # DISCRETE_CHIP demotes to UNRESOLVED per polymer.pyx:1107-1122. Same
        # demotion mechanism as test_stamped_chip_without_src_pool_demotes_to_unresolved.
        rxn = Reaction(reactants=[spx["G"]], products=[spx["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5   # DISCRETE_CHIP
        rxn.polymer_chip_units = 1
        # Pool 'epdm' still carries k_scission=1.0 so the k-guard would pass --
        # the ONLY reason the census stays empty is that the chip demoted.
        pool_a = PolymerPoolConfig(
            label="epdm", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=1.0, k_unzip=0.0, tail_kinetics=None)
        pool_b = PolymerPoolConfig(
            label="B", xs=2, explicit_dp_to_species_index={},
            mu_indices=(5, 6, 7), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0}, V_poly=1.0,
            polymer_pools=[pool_a, pool_b], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"epdm": (1.0, 5.0, 30.0), "B": (1.0, 5.0, 30.0)},
            termination=[])
        rs.initialize_model(core, [rxn], [], [])

        # Demotion-proof: assert the fixture GENUINELY demoted so the test isn't
        # vacuous (it must not stay empty merely because the chip never resolved
        # to a scission/chip archetype in the first place).
        assert rs.reaction_src_pool[0] == -1                          # unresolved src
        assert rs.reaction_flux_archetype[0] == sp.FLUX_UNRESOLVED    # demoted (archetype 4)

        # Payoff: even though pool 'epdm' has k_scission=1.0, the DEMOTED chip
        # does NOT trip the double-count tripwire -- census stays dormant.
        assert rs.double_count_census == []

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
            # k_unzip stays 0.0: it is incidental here, and k_unzip > 0 with
            # monomer_poly_index=None now trips the unzip invariant FIRST --
            # also a ValueError, which would let this test pass without ever
            # reaching the moment-isolation check it pins.
            k_unzip=0.0,
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

    def test_initialize_model_rejects_unzip_pool_without_monomer_index(self):
        """SOLVER-LEVEL INVARIANT (last line of defense behind the deck helper,
        PolymerPool.to_config and the artifact runner, which all guard the same
        shape upstream): a pool with k_unzip > 0 and monomer_poly_index=None
        makes the residual drain the condensed moments unconditionally
        (polymer.pyx: dmu1_dt -= k_unzip*mu0) while the released-monomer
        emission is gated on monomer_poly_index is not None -- mass would leave
        the condensed phase un-conserved. A directly-constructed
        PolymerPoolConfig bypasses every upstream guard, so
        validate_configuration (invoked by initialize_model) must refuse it,
        naming the pool."""
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        core_species = [Inert, Mu0, Mu1, Mu2]
        gas_species_mask = np.array([True, False, False, False], dtype=bool)

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_unzip=0.5,
        )
        rxn_system = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask, constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            termination=[],
        )
        with pytest.raises(ValueError, match=r"poly.*k_unzip.*un-conserved"):
            rxn_system.initialize_model(core_species, [], [], [])

    def test_initialize_model_rejects_negative_k_unzip(self):
        """A negative k_unzip is not a valid rate constant and must raise at the
        same solver invariant layer -- NOT silently become an inert channel
        (every k_unzip consumer in the residual is gated on k_unzip > 0, so a
        negative value would masquerade as a frozen pool while the config is
        nonsense). A wired monomer_poly_index must not dodge the check."""
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")  # GAS released-monomer target (2026-07-03 fix)
        core_species = [Inert, Mu0, Mu1, Mu2, M]
        gas_species_mask = np.array([True, False, False, False, True],
                                    dtype=bool)

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
            k_unzip=-1.0,
        )
        rxn_system = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask, constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            termination=[],
        )
        with pytest.raises(ValueError,
                           match=r"poly.*not a valid rate constant"):
            rxn_system.initialize_model(core_species, [], [], [])

    @pytest.mark.parametrize("field", ["k_unzip", "k_scission"])
    @pytest.mark.parametrize("bad", [float("nan"), float("inf"),
                                     float("-inf")],
                             ids=["nan", "+inf", "-inf"])
    def test_initialize_model_rejects_non_finite_rates(self, field, bad):
        """NaN passes BOTH the `< 0` and `> 0` gates as False, so a
        non-finite k_unzip/k_scission on a directly-constructed
        PolymerPoolConfig would make the channel SILENTLY INERT (NaN) or
        poison the residual (inf) -- a laundered no-op behind every upstream
        guard. validate_configuration (invoked by initialize_model) is the
        last line of defense and must refuse it, naming the pool."""
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")  # GAS released-monomer target (2026-07-03 fix)
        core_species = [Inert, Mu0, Mu1, Mu2, M]
        gas_species_mask = np.array([True, False, False, False, True],
                                    dtype=bool)

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
            **{field: bad},
        )
        rxn_system = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask, constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            termination=[],
        )
        with pytest.raises(ValueError,
                           match=rf"poly.*{field}.*not a valid rate constant"):
            rxn_system.initialize_model(core_species, [], [], [])

    # ------------------------------------------------------------------
    # radical_qssa_unzip channel (M1: config + validation only, NO RHS)
    # ------------------------------------------------------------------

    @staticmethod
    def _qssa_system(pool):
        """(rxn_system, core_species) for the radical_qssa_unzip solver-invariant
        fixtures: gas N2, one pool with mu at 1-3 and a condensed released-
        monomer slot M at 4."""
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")  # GAS released-monomer target (2026-07-03 fix)
        core_species = [Inert, Mu0, Mu1, Mu2, M]
        gas_species_mask = np.array([True, False, False, False, True],
                                    dtype=bool)
        rxn_system = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask, constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            termination=[],
        )
        return rxn_system, core_species

    def test_pool_config_radical_qssa_unzip_field_defaults_to_none(self):
        """PolymerPoolConfig grows a radical_qssa_unzip field (plain dict,
        default None). Channel-absent pools are completely unaffected."""
        pool = PolymerPoolConfig(label="poly", xs=2,
                                 explicit_dp_to_species_index={},
                                 mu_indices=(1, 2, 3))
        assert pool.radical_qssa_unzip is None

    def test_initialize_model_accepts_valid_radical_qssa_unzip_config(self):
        """GREEN path: a valid channel dict (mandatory Arrhenius triplets,
        defaults for efficiency/monomer_yield/basis/transfer) with the monomer
        routing wired and k_unzip=0 passes validate_configuration and is
        stored on the pool config."""
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=_qssa_channel(),
        )
        rxn_system, core_species = self._qssa_system(pool)

        rxn_system.initialize_model(core_species, [], [], [])

        assert rxn_system.polymer_pools[0].radical_qssa_unzip["initiation"] == \
            dict(A=1.0e15, n=0.0, Ea=3.0e5)

    def test_initialize_model_rejects_qssa_channel_without_monomer_index(self):
        """SOLVER-LEVEL INVARIANT (last line of defense; a directly-constructed
        PolymerPoolConfig bypasses the deck helper and to_config): a channel
        without a released-monomer emission target would, once the M2 rate law
        lands, let depropagated mass leave the condensed phase un-conserved.
        Refuse it now, naming the pool."""
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None, k_unzip=0.0,
            radical_qssa_unzip=_qssa_channel(),
        )
        rxn_system, core_species = self._qssa_system(pool)

        with pytest.raises(ValueError,
                           match=r"poly.*radical_qssa_unzip.*un-conserved"):
            rxn_system.initialize_model(core_species, [], [], [])

    def test_initialize_model_rejects_qssa_channel_with_positive_k_unzip(self):
        """Double-counting guard at the solver invariant: radical_qssa_unzip
        and k_unzip > 0 are two representations of the same chain-end
        depropagation channel -- mutually exclusive on a pool."""
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.1,
            radical_qssa_unzip=_qssa_channel(),
        )
        rxn_system, core_species = self._qssa_system(pool)

        with pytest.raises(ValueError, match=r"poly.*mutually exclusive"):
            rxn_system.initialize_model(core_species, [], [], [])

    @pytest.mark.parametrize("channel, pattern", [
        pytest.param({k: v for k, v in _qssa_channel().items()
                      if k != "depropagation"},
                     r"poly.*radical_qssa_unzip.*missing.*depropagation",
                     id="missing-depropagation-block"),
        pytest.param(_qssa_channel(initiation=_qssa_triplet(A=float("nan"))),
                     r"poly.*initiation.*A.*not finite", id="nan-A"),
        pytest.param(_qssa_channel(termination=_qssa_triplet(Ea=float("inf"))),
                     r"poly.*termination.*Ea.*not finite", id="inf-Ea"),
        pytest.param(_qssa_channel(depropagation=_qssa_triplet(A=0.0)),
                     r"poly.*depropagation.*A.*> 0", id="zero-A"),
        pytest.param(_qssa_channel(efficiency=0.0),
                     r"poly.*efficiency.*\(0, 1\]", id="efficiency-zero"),
        pytest.param(_qssa_channel(monomer_yield=1.5),
                     r"poly.*monomer_yield.*\(0, 1\]",
                     id="monomer-yield-above-one"),
        pytest.param(_qssa_channel(basis="chain_ends_mu0"),
                     r"poly.*basis.*backbone_bonds_mu1_minus_mu0",
                     id="bad-basis"),
    ])
    def test_initialize_model_rejects_malformed_qssa_channel(self, channel,
                                                             pattern):
        """The solver re-validates the channel dict itself (finite A/n/Ea,
        A > 0, Ea >= 0, efficiency/monomer_yield in (0, 1], pinned basis):
        a directly-constructed PolymerPoolConfig bypasses every upstream
        guard."""
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=channel,
        )
        rxn_system, core_species = self._qssa_system(pool)

        with pytest.raises(ValueError, match=pattern):
            rxn_system.initialize_model(core_species, [], [], [])

    def test_initialize_model_weaklink_channel_initializes_cleanly(self):
        """Milestones ii+iii REPLACE the milestone-i anti-silent-no-op guard:
        a weak-link allyl/U-state channel (initiation_allyl, split
        terminations, unsaturated_tail_ends_initial) now has real solver
        support -- initialization must succeed and expose the flattened
        weak-link state (the guard test this replaces pinned the refusal)."""
        channel = dict(
            initiation=_qssa_triplet(A=1.0e15, Ea=3.0e5),
            depropagation=_qssa_triplet(A=1.0e13, Ea=8.0e4),
            initiation_allyl=_qssa_triplet(A=2.0e14, Ea=2.4e5),
            termination_recombination=_qssa_triplet(A=6.0e7, Ea=8.0e3),
            termination_disproportionation=_qssa_triplet(A=4.0e7, Ea=1.2e4),
            unsaturated_tail_ends_initial=0.02,
        )
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=channel,
        )
        rxn_system, core_species = self._qssa_system(pool)

        rxn_system.initialize_model(core_species, [], [], [])

        rs = rxn_system
        assert rs.qssa_enabled[0] == 1
        assert rs.qssa_weaklink[0] == 1
        assert rs.qssa_kia_A[0] == 2.0e14
        assert rs.qssa_kia_n[0] == 0.0
        assert rs.qssa_kia_Ea[0] == 2.4e5
        assert rs.qssa_ktrec_A[0] == 6.0e7
        assert rs.qssa_ktrec_Ea[0] == 8.0e3
        assert rs.qssa_ktdisp_A[0] == 4.0e7
        assert rs.qssa_ktdisp_Ea[0] == 1.2e4
        assert rs.qssa_u0[0] == 0.02
        # legacy summed-kt slots stay zero on a weak-link pool: the RHS must
        # gate on qssa_weaklink, never on kt_A != 0.
        assert rs.qssa_kt_A[0] == 0.0

    def test_initialize_model_rejects_weaklink_with_legacy_termination(self):
        """The solver is the last line of defense for the mutual-exclusion
        rule too: a directly-constructed config carrying BOTH the legacy
        summed 'termination' and the weak-link vocabulary is refused by the
        shared validator before the not-implemented guard is even reached."""
        channel = dict(
            _qssa_channel(),  # legacy, includes summed termination
            initiation_allyl=_qssa_triplet(A=2.0e14, Ea=2.4e5),
            termination_recombination=_qssa_triplet(A=6.0e7, Ea=8.0e3),
            termination_disproportionation=_qssa_triplet(A=4.0e7, Ea=1.2e4),
            unsaturated_tail_ends_initial=0.02,
        )
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=channel,
        )
        rxn_system, core_species = self._qssa_system(pool)

        with pytest.raises(ValueError,
                           match=r"poly.*'termination'.*mutually exclusive"):
            rxn_system.initialize_model(core_species, [], [], [])

    def test_radical_qssa_unzip_footprint_confined_to_signature(self):
        """M2 supersedes the M1 zero-RHS pin on the channel-ON side: the
        channel is now LIVE, and its residual footprint must be confined to
        the documented signature -- mu1 (drain), mu2 (drain) and the released
        monomer slot (emission). Everything else (mu0, the gas inert, every
        shared k_scission contribution) must stay BITWISE identical between
        two systems that differ only by the channel; the channel-OFF path is
        untouched."""
        pool_kwargs = dict(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
            k_scission=0.3, k_unzip=0.0,
        )
        rs_off, core_off = self._qssa_system(
            PolymerPoolConfig(**pool_kwargs))
        rs_on, core_on = self._qssa_system(
            PolymerPoolConfig(radical_qssa_unzip=_qssa_channel(),
                              **pool_kwargs))
        rs_off.initialize_model(core_off, [], [], [])
        rs_on.initialize_model(core_on, [], [], [])

        dn_off = rs_off.residual(0.0, rs_off.y, np.zeros_like(rs_off.y))[0]
        dn_on = rs_on.residual(0.0, rs_on.y, np.zeros_like(rs_on.y))[0]

        assert np.any(dn_off != 0.0)  # the fixture has live dynamics
        diff = dn_on - dn_off
        # M2: the channel is live -- it must move mu1, mu2 and the monomer.
        assert diff[2] < 0.0   # mu1 drained
        assert diff[3] < 0.0   # mu2 drained
        assert diff[4] > 0.0   # monomer emitted
        # ... and NOTHING else: mu0 and the gas inert are bitwise-shared.
        assert dn_on[0] == dn_off[0]
        assert dn_on[1] == dn_off[1]

    def test_direct_construction_qssa_channel_is_normalized_in_storage(self):
        """NORMALIZED CONSUMPTION (review round 21, finding 1): a directly
        constructed PolymerPoolConfig carrying a minimal-but-valid channel
        (mandatory blocks only) must come OUT of initialize_model with the
        validator's normalized form stored back on the config -- defaults
        filled (efficiency=1.0, monomer_yield=1.0, transfer=None, basis
        pinned), numerics coerced to float. Without the store-back, future
        RHS code doing q['efficiency'] would KeyError on exactly this
        (validation-passing) config."""
        minimal = dict(
            initiation=dict(A=1e15, n=0, Ea=3e5),
            depropagation=dict(A=1e13, n=0, Ea=8e4),
            termination=dict(A=1e8, n=0, Ea=1e4),
        )
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=minimal,
        )
        rxn_system, core_species = self._qssa_system(pool)
        rxn_system.initialize_model(core_species, [], [], [])

        q = rxn_system.polymer_pools[0].radical_qssa_unzip
        assert q["efficiency"] == 1.0
        assert q["monomer_yield"] == 1.0
        assert q["transfer"] is None
        assert q["basis"] == "backbone_bonds_mu1_minus_mu0"
        assert q["initiation"] == dict(A=1.0e15, n=0.0, Ea=3.0e5)
        assert isinstance(q["initiation"]["n"], float)

    def test_qssa_stored_channel_does_not_alias_caller_dict(self):
        """NORMALIZED CONSUMPTION, aliasing half: the stored normalized
        structure must be defensively deep-copied -- mutating the caller's
        original dict (or its nested triplets) after initialize_model must
        not reach the stored config."""
        ch = _qssa_channel()
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=ch,
        )
        rxn_system, core_species = self._qssa_system(pool)
        rxn_system.initialize_model(core_species, [], [], [])

        ch["initiation"]["A"] = float("nan")
        ch["efficiency"] = -1.0

        q = rxn_system.polymer_pools[0].radical_qssa_unzip
        assert q["initiation"]["A"] == 1.0e15
        assert q["efficiency"] == 1.0

    def test_qssa_channel_flattened_into_solver_state(self):
        """M2 prep (review round 21, finding 1, flattening half): the
        validated+normalized channel is flattened into solver-owned per-pool
        flat arrays at initialization. The M2 rate law will read ONLY these
        arrays, never the dict."""
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=_qssa_channel(
                transfer=_qssa_triplet(A=2.0e7, n=0.5, Ea=5.0e4),
                efficiency=0.8, monomer_yield=0.9),
        )
        rxn_system, core_species = self._qssa_system(pool)
        rxn_system.initialize_model(core_species, [], [], [])

        rs = rxn_system
        assert rs.qssa_enabled.shape == (1,)
        assert rs.qssa_enabled[0] == 1
        assert rs.qssa_ki_A[0] == 1.0e15
        assert rs.qssa_ki_n[0] == 0.0
        assert rs.qssa_ki_Ea[0] == 3.0e5
        assert rs.qssa_kdp_A[0] == 1.0e13
        assert rs.qssa_kdp_n[0] == 0.0
        assert rs.qssa_kdp_Ea[0] == 8.0e4
        assert rs.qssa_kt_A[0] == 1.0e8
        assert rs.qssa_kt_n[0] == 0.0
        assert rs.qssa_kt_Ea[0] == 1.0e4
        assert rs.qssa_efficiency[0] == 0.8
        assert rs.qssa_monomer_yield[0] == 0.9
        assert rs.qssa_has_transfer[0] == 1
        assert rs.qssa_ktr_A[0] == 2.0e7
        assert rs.qssa_ktr_n[0] == 0.5
        assert rs.qssa_ktr_Ea[0] == 5.0e4

    def test_qssa_flattened_state_channel_absent_pool(self):
        """A channel-absent pool flattens to the disabled/default row:
        enabled=0, all Arrhenius slots 0, efficiency/monomer_yield 1 (inert
        defaults), has_transfer=0."""
        pool = PolymerPoolConfig(label="poly", xs=2,
                                 explicit_dp_to_species_index={},
                                 mu_indices=(1, 2, 3))
        rxn_system, core_species = self._qssa_system(pool)
        rxn_system.initialize_model(core_species, [], [], [])

        rs = rxn_system
        assert rs.qssa_enabled[0] == 0
        assert rs.qssa_ki_A[0] == 0.0
        assert rs.qssa_kdp_A[0] == 0.0
        assert rs.qssa_kt_A[0] == 0.0
        assert rs.qssa_has_transfer[0] == 0
        assert rs.qssa_ktr_A[0] == 0.0
        assert rs.qssa_efficiency[0] == 1.0
        assert rs.qssa_monomer_yield[0] == 1.0

    def test_qssa_dict_mutation_after_initialize_cannot_reach_solver_state(self):
        """MUTATION BYPASS (review round 21, finding 2): the nested dict is
        mutable behind the frozen dataclass. With the flattened solver state,
        post-initialization dict mutation is harmless -- corrupt the stored
        dict aggressively and PROVE the flattened arrays (the only thing the
        M2 rate law will read) are unchanged."""
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=_qssa_channel(),
        )
        rxn_system, core_species = self._qssa_system(pool)
        rxn_system.initialize_model(core_species, [], [], [])
        rs = rxn_system

        snapshot = {name: getattr(rs, name).copy() for name in (
            "qssa_enabled", "qssa_ki_A", "qssa_ki_n", "qssa_ki_Ea",
            "qssa_kdp_A", "qssa_kdp_n", "qssa_kdp_Ea",
            "qssa_kt_A", "qssa_kt_n", "qssa_kt_Ea",
            "qssa_efficiency", "qssa_monomer_yield",
            "qssa_has_transfer", "qssa_ktr_A", "qssa_ktr_n", "qssa_ktr_Ea")}

        q = rs.polymer_pools[0].radical_qssa_unzip
        q["termination"]["A"] = float("nan")
        q["initiation"]["Ea"] = -1.0e9
        q["efficiency"] = float("inf")
        q["transfer"] = dict(A=float("nan"), n=0.0, Ea=0.0)

        for name, before in snapshot.items():
            assert np.array_equal(getattr(rs, name), before), name
        assert rs.qssa_kt_A[0] == 1.0e8
        assert math.isfinite(rs.qssa_efficiency[0])

    # ------------------------------------------------------------------
    # radical_qssa_unzip channel (M2: solver RHS -- rate law + moments +
    # handshake + census)
    # ------------------------------------------------------------------

    @staticmethod
    def _qssa_m2_pool(channel, explicit=None, k_unzip=0.0, k_scission=0.0):
        """One-pool config for the M2 RHS fixtures: mu at 1-3, condensed
        released-monomer slot M at 4, optional explicit DP=2 species at 5."""
        return PolymerPoolConfig(
            label="poly", xs=2,
            explicit_dp_to_species_index=dict(explicit) if explicit else {},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
            k_scission=k_scission, k_unzip=k_unzip, tail_kinetics=None,
            radical_qssa_unzip=channel)

    @staticmethod
    def _qssa_m2_system(pool, moments=(1.0, 5.0, 30.0), V_poly=1.0,
                        explicit_moles=0.0):
        """Initialized HybridPolymerSystem for the M2 RHS fixtures: gas N2 at
        0, mu dummies at 1-3, condensed monomer M at 4, explicit DP=2 chain
        P2 at 5. ``moments`` are MOLES of moments (mu_k = Mu_k / V_poly);
        ``explicit_moles`` seeds the P2 slot (a pool built with
        explicit={2: 5} then carries explicit-species moment contributions
        that the step-6 consistency pass subtracts from these totals)."""
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")
        P2 = _spc("CC", "P2")
        core_species = [Inert, Mu0, Mu1, Mu2, M, P2]
        # M (released-monomer target) is GAS since the 2026-07-03 monomer-gas
        # fix; the mu dummies and explicit P2 stay condensed.
        gas_species_mask = np.array(
            [True, False, False, False, True, False], dtype=bool)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0},
            V_poly=V_poly, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask, constant_gas_volume=False,
            initial_polymer_moments={"poly": tuple(moments)},
            initial_explicit_species={"poly": {2: explicit_moles}},
            termination=[],
        )
        rs.initialize_model(core_species, [], [], [])
        return rs

    @staticmethod
    def _qssa_oracle_rate(channel, mu0, mu1, T=800.0):
        """Independent analytic recompute of the QSSA monomer-release rate
        r_mono [mol m^-3 s^-1 of condensed volume] -- deliberately NOT
        imported from production code. Concentration basis mol/m^3 (SI),
        matching the solver's C_poly = y/V_poly convention; kt is
        m^3 mol^-1 s^-1 on that basis, ktr is s^-1 (pseudo-first-order:
        it multiplies R directly in ktr*R + 2*kt*R^2 = G_R).
        R_gas = 8.314 J/mol/K (M1 SI pin)."""
        R_gas = 8.314

        def k(tri):
            return tri["A"] * T ** tri["n"] * math.exp(-tri["Ea"] / (R_gas * T))

        f = channel.get("efficiency", 1.0)
        y_m = channel.get("monomer_yield", 1.0)
        B = max(mu1 - mu0, 0.0)
        if B <= 0.0:
            return 0.0
        ki, kdp, kt = (k(channel["initiation"]), k(channel["depropagation"]),
                       k(channel["termination"]))
        transfer = channel.get("transfer")
        if transfer is None:
            R_ss = math.sqrt(f * ki * B / kt)
        else:
            ktr = k(transfer)
            G_R = 2.0 * f * ki * B
            R_ss = (math.sqrt(ktr * ktr + 8.0 * kt * G_R) - ktr) / (4.0 * kt)
        return y_m * kdp * R_ss

    def _qssa_channel_diff(self, channel, moments=(1.0, 5.0, 30.0),
                           V_poly=1.0, state=None):
        """Residual difference (channel-ON minus channel-OFF) of two systems
        identical but for the channel: isolates the QSSA contribution."""
        rs_on = self._qssa_m2_system(self._qssa_m2_pool(channel),
                                     moments, V_poly)
        rs_off = self._qssa_m2_system(self._qssa_m2_pool(None),
                                      moments, V_poly)
        y = rs_on.y.copy()
        if state is not None:
            y[1], y[2], y[3] = state
        dn_on = rs_on.residual(0.0, y, np.zeros_like(y))[0]
        dn_off = rs_off.residual(0.0, y.copy(), np.zeros_like(y))[0]
        return dn_on - dn_off

    def test_qssa_rate_law_analytic_oracle_no_transfer(self):
        """Analytic rate oracle, no transfer: the net dmu1 contribution of the
        channel is exactly -y_m*kdp(T)*sqrt(f*ki(T)*B/kt(T)) with
        B = mu1 - mu0 (breakable-bond basis), at the solver temperature."""
        channel = _qssa_channel(efficiency=0.6, monomer_yield=0.9)
        mu0, mu1, mu2 = 1.0, 5.0, 30.0
        r = self._qssa_oracle_rate(channel, mu0, mu1)
        assert r > 0.0  # the oracle itself is live

        diff = self._qssa_channel_diff(channel, moments=(mu0, mu1, mu2))
        assert diff[2] == pytest.approx(-r, rel=1e-10)

    def test_qssa_rate_law_analytic_oracle_with_transfer(self):
        """Analytic rate oracle WITH transfer: R_ss solves the quadratic
        active-end balance ktr*R + 2*kt*R^2 = G_R, i.e.
        R_ss = (sqrt(ktr^2 + 8*kt*G_R) - ktr)/(4*kt). Transfer is a
        first-order active-end LOSS, so it must strictly REDUCE the rate --
        pinning that the accepted transfer block is implemented, not
        silently ignored.

        UNITS PIN: ktr multiplies R DIRECTLY in the implemented balance, so
        the configured transfer triplet yields a PSEUDO-FIRST-ORDER rate
        constant, ktr(T) [s^-1] -- NOT a bimolecular m^3 mol^-1 s^-1 one. A
        literature bimolecular k_tr [L mol^-1 s^-1 / m^3 mol^-1 s^-1] must be
        premultiplied by the relevant substrate concentration [mol/m^3, SI]
        BEFORE entering the config. This oracle recomputes the documented
        pseudo-first-order law exactly; if the solver ever applied a hidden
        concentration factor to ktr the rtol-1e-10 match below would break."""
        base = _qssa_channel()
        with_tr = _qssa_channel(transfer=_qssa_triplet(A=1.0e5, Ea=5.0e4))
        mu0, mu1, mu2 = 1.0, 5.0, 30.0
        r_tr = self._qssa_oracle_rate(with_tr, mu0, mu1)
        r_no = self._qssa_oracle_rate(base, mu0, mu1)
        assert 0.0 < r_tr < r_no  # transfer visibly suppresses the rate

        diff = self._qssa_channel_diff(with_tr, moments=(mu0, mu1, mu2))
        assert diff[2] == pytest.approx(-r_tr, rel=1e-10)

    def test_qssa_mass_conservation_monomer_yield_scales_both_sides(self):
        """monomer_yield y_m scales the moment drain AND the emission
        TOGETHER: the monomer emission (mol/s) must equal the mu1 drain
        (mol/s) exactly -- monomer_mw * drain == monomer_mw * emission, so
        net condensed+gas mass is conserved by the channel. y_m=0.7 pins
        that the yield cannot be applied to one side only."""
        channel = _qssa_channel(monomer_yield=0.7)
        diff = self._qssa_channel_diff(channel)
        assert diff[2] < 0.0
        assert diff[4] == pytest.approx(-diff[2], rel=1e-12)
        # ... and y_m scaled BOTH: the drain is 0.7x the y_m=1 drain.
        diff_full = self._qssa_channel_diff(_qssa_channel())
        assert diff[2] == pytest.approx(0.7 * diff_full[2], rel=1e-10)

    def test_qssa_moment_signature(self):
        """Chain-END monomer release signature: dmu0 == 0 (no chain created
        or destroyed), dmu1 == -r, dmu2 == -r*max(2*E[n] - 1, 0) with
        E[n] = mu1/mu0 (same-pool-VE clamp idiom: the drain must never make
        mu2 increase)."""
        channel = _qssa_channel()
        mu0, mu1, mu2 = 1.0, 5.0, 30.0
        r = self._qssa_oracle_rate(channel, mu0, mu1)
        diff = self._qssa_channel_diff(channel, moments=(mu0, mu1, mu2))
        assert diff[1] == 0.0                                    # dmu0
        assert diff[2] == pytest.approx(-r, rel=1e-10)           # dmu1
        assert diff[3] == pytest.approx(-r * (2.0 * mu1 / mu0 - 1.0),
                                        rel=1e-10)               # dmu2
        # mu2 must NEVER increase from this drain, across pathological
        # states too (mu0 ~ 0 guarded by the eps clamp).
        for state in [(1.0, 5.0, 30.0), (0.0, 5.0, 30.0),
                      (1.0e-3, 1.0, 1.0e3), (2.0, 3.0, 30.0)]:
            d = self._qssa_channel_diff(channel, state=state)
            assert np.all(np.isfinite(d)), state
            assert d[3] <= 0.0, state

    def test_qssa_dp1_and_inverted_states_rate_zero_no_nan(self):
        """DP->1 pools (mu1 == mu0) and cone-inverted states (mu1 < mu0):
        B = max(mu1 - mu0, 0) -> 0, so the rate goes to zero SMOOTHLY -- no
        sqrt-of-negative, no NaN, bitwise-zero channel contribution."""
        channel = _qssa_channel()
        for state in [(5.0, 5.0, 30.0),    # DP = 1 exactly
                      (5.0, 4.0, 30.0),    # inverted: mu1 < mu0
                      (0.0, 0.0, 0.0)]:    # fully degenerate
            diff = self._qssa_channel_diff(channel, state=state)
            assert np.all(np.isfinite(diff)), state
            assert np.all(diff == 0.0), state

    def test_qssa_emission_volume_factor(self):
        """The emission lands on monomer_poly_index with the SAME volume
        factor the surrounding small_src code uses: dn_dt[M] = r_mono*V_poly
        with r_mono computed on mu_k = Mu_k/V_poly (mol/m^3 of condensed
        volume). V_poly=2 with doubled moment MOLES keeps mu identical, so
        dn_dt[M] must be exactly 2x the V_poly=1 emission."""
        channel = _qssa_channel()
        mu0, mu1 = 1.0, 5.0
        r = self._qssa_oracle_rate(channel, mu0, mu1)
        diff_v2 = self._qssa_channel_diff(channel,
                                          moments=(2.0, 10.0, 60.0),
                                          V_poly=2.0)
        assert diff_v2[4] == pytest.approx(r * 2.0, rel=1e-10)
        assert diff_v2[2] == pytest.approx(-r * 2.0, rel=1e-10)

    def test_qssa_handshake_equivalence_with_k_unzip(self):
        """Hybrid handshake: a QSSA pool with an explicit-tail boundary must
        produce F = (r_mono/mu0) * N_boundary flux into the explicit DP=xs
        species -- i.e. at a frozen state the ENTIRE residual must match a
        k_unzip pool with k_unzip = r_mono/mu0 (the instantaneous per-chain
        unzip frequency): same drain, same emission, same handshake. QSSA
        pools must not strand low-DP condensed residue."""
        channel = _qssa_channel()
        mu0, mu1, mu2 = 1.0, 5.0, 30.0   # mean DP 5 > xs 2 -> tail valid
        r = self._qssa_oracle_rate(channel, mu0, mu1)
        k_eq = r / mu0

        rs_qssa = self._qssa_m2_system(
            self._qssa_m2_pool(channel, explicit={2: 5}),
            moments=(mu0, mu1, mu2))
        rs_kz = self._qssa_m2_system(
            self._qssa_m2_pool(None, explicit={2: 5}, k_unzip=k_eq),
            moments=(mu0, mu1, mu2))

        dn_qssa = rs_qssa.residual(0.0, rs_qssa.y,
                                   np.zeros_like(rs_qssa.y))[0]
        dn_kz = rs_kz.residual(0.0, rs_kz.y, np.zeros_like(rs_kz.y))[0]

        assert dn_qssa[5] > 0.0    # explicit DP=2 boundary species fed
        assert dn_qssa[1] < 0.0    # mu0 drained by the handshake
        np.testing.assert_allclose(dn_qssa, dn_kz, rtol=1e-9, atol=0.0)

    # ------------------------------------------------------------------
    # weak-link allyl/U-state channel (milestones ii+iii: per-pool U state
    # slot + allyl initiation RHS, nu=1, disproportionation-sourced U)
    # ------------------------------------------------------------------

    def test_legacy_state_layout_and_rhs_golden_frozen(self):
        """REGRESSION PIN 1 (contractual): a legacy config (no weak-link
        keys) keeps the EXACT pre-milestone state layout (neq == n_core, no
        U slot) and the EXACT pre-milestone RHS derivatives. The golden
        values below were captured from the pre-implementation build
        (commit 7e1fa0671) on this fixed fixture: legacy _qssa_channel() +
        k_scission=0.3, moments (1, 5, 30) mol, V_poly=1, T=800 K. Bitwise
        equality via float.fromhex -- any drift in the legacy path fails
        here."""
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_qssa_channel(), k_scission=0.3))
        assert rs.neq == 6              # n_core exactly: NO U slot appended
        assert len(rs.y0) == 6
        assert rs.num_qssa_u == 0
        assert rs.qssa_u_slot[0] == -1

        dn = rs.residual(0.0, rs.y0.copy(), np.zeros(6))[0]
        golden = [
            0.0,                                     # gas inert
            float.fromhex("0x1.3333333333333p+0"),   # dmu0 (k_scission)
            float.fromhex("-0x1.015b879ec323fp+7"),  # dmu1 (QSSA drain)
            float.fromhex("-0x1.26cd5ef901eedp+10"), # dmu2
            float.fromhex("0x1.015b879ec323fp+7"),   # monomer emission
            0.0,                                     # explicit P2 slot
        ]
        assert dn.shape == (6,)
        for i, g in enumerate(golden):
            assert dn[i] == g, (i, dn[i].hex() if dn[i] else dn[i], g)

    def test_weaklink_state_layout_appends_u_slot(self):
        """A weak-link pool gets EXACTLY ONE extra ODE slot, appended after
        the core-species block, initialized to unsaturated_tail_ends_initial
        [mol -- same amount basis as mu0]. Tolerance arrays follow neq."""
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_channel(
                unsaturated_tail_ends_initial=0.5)))
        assert rs.num_qssa_u == 1
        assert rs.neq == 7
        assert len(rs.y0) == 7
        assert rs.qssa_u_slot[0] == 6          # first slot after n_core=6
        assert rs.y0[6] == 0.5
        assert len(rs.atol_array) == 7
        assert len(rs.rtol_array) == 7

    def test_weaklink_sensitivity_layout_not_supported(self):
        """The DASPK sensitivity layout interleaves per-species blocks; a
        trailing U slot does not fit it. Refuse loudly instead of silently
        corrupting the sensitivity state vector."""
        pool = self._qssa_m2_pool(_weaklink_channel())
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")
        P2 = _spc("CC", "P2")
        core = [Inert, Mu0, Mu1, Mu2, M, P2]
        mask = np.array([True] + [False] * 5, dtype=bool)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0},
            V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask, constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            initial_explicit_species={"poly": {2: 0.0}}, termination=[])
        with pytest.raises(ValueError, match=r"sensitivity"):
            rs.initialize_model(core, [], [], [], sensitivity=True)

    def test_weaklink_bridge_matches_legacy_at_u_zero(self):
        """REGRESSION PIN 7 (bridge) + PIN 2a: a weak-link pool with
        kt_rec(T) + kt_disp(T) == the legacy summed kt(T) EXACTLY (n=0/Ea=0
        triplets, A_rec = A_disp = A_legacy/2, both exactly representable)
        and U = 0 must reproduce the legacy pool's moment derivatives
        BITWISE at t=0: same G_R (the allyl term contributes exactly 0 at
        U=0), same R_ss, same drains, same emission. This proves the kt
        split did not change the base algebra. The U slot itself grows
        (disproportionation sources U) -- pinned separately."""
        legacy = dict(initiation=_exact_triplet(1.0e-3),
                      depropagation=_exact_triplet(1.0e2),
                      termination=_exact_triplet(1.0e8))
        weak = _weaklink_exact_channel(ki_A=1.0e-3, kdp_A=1.0e2,
                                       kia_A=1.0e3, ktrec_A=5.0e7,
                                       ktdisp_A=5.0e7, u0=0.0)
        rs_leg = self._qssa_m2_system(self._qssa_m2_pool(legacy))
        rs_weak = self._qssa_m2_system(self._qssa_m2_pool(weak))

        dn_leg = rs_leg.residual(0.0, rs_leg.y0.copy(), np.zeros(6))[0]
        dn_weak = rs_weak.residual(0.0, rs_weak.y0.copy(), np.zeros(7))[0]

        assert dn_leg[2] < 0.0  # the fixture's channel is live
        for i in range(6):
            assert dn_weak[i] == dn_leg[i], i  # bitwise bridge
        # PIN 2b/4a: kt_disp legally > 0 and R_ss > 0 -> U grows.
        assert dn_weak[6] > 0.0
        # ... by exactly the disproportionation event rate kt_disp*R_ss^2
        # (halved-radical-disappearance convention: kt_disp*R^2 IS the
        # event rate; 1 unsaturated end per event), amount basis *V_poly.
        B = 4.0  # mu1 - mu0 = 5 - 1 [mol/m^3, V_poly = 1]
        R_ss = math.sqrt(1.0e-3 * B / 1.0e8)  # sqrt(f*ki*B/kt_total), f=1
        assert dn_weak[6] == pytest.approx(5.0e7 * R_ss * R_ss, rel=1e-12)

    def test_weaklink_all_zero_state_is_inert(self):
        """PIN 2b (exact construction): B = 0 (mu1 == mu0) AND U = 0 -> the
        whole weak-link channel contributes exactly nothing: every residual
        entry, including dU/dt, is bitwise 0. U stays exactly 0 when there
        is no radical source (validator forbids A=0, so this zero-source
        state is the exact pin; production requires R_ss > 0)."""
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel(u0=0.0)),
            moments=(1.0, 1.0, 1.0))  # all-monomer: B = mu1 - mu0 = 0
        dn = rs.residual(0.0, rs.y0.copy(), np.zeros(7))[0]
        assert np.all(dn == 0.0)

    def test_weaklink_allyl_nu_is_one(self):
        """PIN 3: the allyl channel enters G_R with stoichiometry nu = 1
        (ONE unzipping radical per weak-link fission; the allylic
        co-fragment does not unzip): dG_R/du == f*ki_allyl EXACTLY, not
        2*f*ki_allyl. Constructed on a B > 0 state with u varied BELOW both
        B (the r36 active-site clamp u_active = min(u, B) stays inactive)
        and 2*mu0 (capacity); the allyl contribution is isolated as
        G_R(u) - G_R(0). G_R is recovered from the no-transfer law
        r = y_m*kdp*sqrt(G_R/(2*kt)) => G_R = 2*kt*(r/(y_m*kdp))^2."""
        ki, kia, kdp, kt = 1.0e-3, 2.0e4, 1.0e2, 1.0e8  # exact, f = y_m = 1
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel(
                ki_A=ki, kia_A=kia, kdp_A=kdp,
                ktrec_A=kt / 2.0, ktdisp_A=kt / 2.0)),
            moments=(1.0, 5.0, 30.0))  # B = 4; capacity 2*mu0 = 2
        u_slot = rs.qssa_u_slot[0]

        def g_r_at(u_mol):
            y = rs.y0.copy()
            y[u_slot] = u_mol
            dn = rs.residual(0.0, y, np.zeros(7))[0]
            r = -dn[2]  # dmu1 = -r_mono, V_poly = 1
            return 2.0 * kt * (r / kdp) ** 2

        g0 = g_r_at(0.0)
        assert g0 == pytest.approx(2.0 * ki * 4.0, rel=1e-10)  # chain term
        u1, u2 = 0.02, 0.04  # << B = 4 and << capacity 2: clamp inactive
        g1, g2 = g_r_at(u1), g_r_at(u2)
        # levels: G_R(u) - G_R(0) == 1 * kia * u (u = U/V_poly, V_poly = 1)
        assert g1 - g0 == pytest.approx(kia * u1, rel=1e-10)
        assert g2 - g0 == pytest.approx(kia * u2, rel=1e-10)
        # derivative: dG_R/du == f*ki_allyl exactly (nu = 1, never 2)
        dg_du = (g2 - g1) / (u2 - u1)
        assert dg_du == pytest.approx(kia, rel=1e-9)
        assert abs(dg_du - 2.0 * kia) > 0.5 * kia  # nu = 2 is far away

    def test_weaklink_u_clamped_negative_state(self):
        """u enters rates as max(U, 0)/V_poly: a (solver-transient) negative
        U state contributes exactly like U = 0 -- bitwise-identical residual
        (chain initiation still live, allyl term zero, production throttle
        at its U=0 value), no negative G_R, no NaN."""
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel()),
            moments=(1.0, 5.0, 30.0))  # B = 4 > 0: the channel is live
        u_slot = rs.qssa_u_slot[0]
        y_neg = rs.y0.copy()
        y_neg[u_slot] = -1.0e-3
        dn_neg = rs.residual(0.0, y_neg, np.zeros(7))[0]
        y_zero = rs.y0.copy()
        y_zero[u_slot] = 0.0
        dn_zero = rs.residual(0.0, y_zero, np.zeros(7))[0]
        assert np.all(np.isfinite(dn_neg))
        assert dn_zero[2] < 0.0  # live fixture, not vacuous
        assert np.array_equal(dn_neg, dn_zero)

    def test_weaklink_allyl_inert_when_no_backbone_bonds(self):
        """r36 P1-2: the weak channel fissions a backbone bond ALLYLIC to
        the unsaturated end -- with B = 0 there is nothing to fission and
        nothing to depropagate. B = 0 AND U > 0 must be completely inert:
        no G_R, no monomer drain (which could push mu1 < mu0 and defeat the
        DP->1 self-termination floor), no U consumption, no U production."""
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel()),
            moments=(1.0, 1.0, 1.0))  # all-monomer: B = 0
        u_slot = rs.qssa_u_slot[0]
        y = rs.y0.copy()
        y[u_slot] = 0.5  # U > 0 but no backbone bonds to fission
        dn = rs.residual(0.0, y, np.zeros(7))[0]
        assert np.all(dn == 0.0)

    def test_weaklink_active_site_clamp_u_above_b(self):
        """r36 P1-2: rates use the ACTIVE u, u_active = min(max(U,0)/V_poly,
        B): unsaturated ends beyond the remaining backbone-bond count have
        no bond to fission. With u > B both the G_R term and the U sink
        must saturate at B: the moment derivatives at U = 3.5 are BITWISE
        those at U = B = 3 (fixture: moments (2, 5, 30) -> B = 3, capacity
        2*mu0 = 4, so U in (3, 4) is above B yet below capacity -- the
        clamp is exercised independently of the capacity throttle)."""
        ki, kia, kdp = 1.0e-3, 2.0e4, 1.0e2
        ktrec, ktdisp = 5.0e7, 5.0e7
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel(
                ki_A=ki, kia_A=kia, kdp_A=kdp,
                ktrec_A=ktrec, ktdisp_A=ktdisp)),
            moments=(2.0, 5.0, 30.0))
        u_slot = rs.qssa_u_slot[0]

        def dn_at(u_mol):
            y = rs.y0.copy()
            y[u_slot] = u_mol
            return rs.residual(0.0, y, np.zeros(7))[0]

        dn_at_b = dn_at(3.0)      # u == B exactly
        dn_above = dn_at(3.5)     # u > B, still below capacity 4
        for i in range(6):        # clamped: same G_R -> bitwise moments
            assert dn_above[i] == dn_at_b[i], i
        # ... and the analytic dU at U = 3.5: sink uses u_active = B = 3,
        # production throttled by 1 - U/(2*mu0) = 1 - 3.5/4:
        G_R = 2.0 * ki * 3.0 + kia * 3.0
        R_sq = G_R / (2.0 * (ktrec + ktdisp))
        expected = ktdisp * R_sq * (1.0 - 3.5 / 4.0) - kia * 3.0
        assert dn_above[6] == pytest.approx(expected, rel=1e-12)
        # a u-not-clamped sink law would differ by kia*0.5:
        assert abs(dn_above[6] - (ktdisp * R_sq * (1.0 - 3.5 / 4.0)
                                  - kia * 3.5)) > 0.4 * kia

    def test_weaklink_efficiency_symmetric_on_allyl_channel(self):
        """r36 P1-1: f multiplies BOTH allyl terms. f is the escaped-
        radical-pair efficiency: a caged recombination restores the allylic
        bond, so U is NOT consumed by the caged fraction --
        dU/dt sink = f*ki_allyl*u_active, never the bare ki_allyl*u_active
        (which would destroy U without producing radicals for f < 1).
        Production keeps NO f (R_ss already contains the escape efficiency;
        production is a termination event)."""
        f, kia = 0.5, 2.0e4
        ktrec = ktdisp = 0.5e-2  # kt_total = 1e-2
        channel = _weaklink_exact_channel(
            ki_A=1.0e-9, kdp_A=1.0e-6, kia_A=kia,
            ktrec_A=ktrec, ktdisp_A=ktdisp)
        channel["efficiency"] = f
        rs = self._qssa_m2_system(self._qssa_m2_pool(channel),
                                  moments=(1.0, 5.0, 30.0))  # B=4, cap=2
        u_slot = rs.qssa_u_slot[0]
        u = 1.0  # < B (clamp inactive), = cap/2 (throttle = 0.5)
        y = rs.y0.copy()
        y[u_slot] = u
        dn = rs.residual(0.0, y, np.zeros(7))[0]
        G_R = 2.0 * f * 1.0e-9 * 4.0 + f * kia * u
        R_sq = G_R / (2.0 * (ktrec + ktdisp))
        expected = ktdisp * R_sq * (1.0 - u / 2.0) - f * kia * u
        assert dn[u_slot] == pytest.approx(expected, rel=1e-12)
        # the f-asymmetric law (bare kia*u sink) is far away:
        assert abs(dn[u_slot] - (ktdisp * R_sq * (1.0 - u / 2.0)
                                 - kia * u)) > 0.4 * kia * u

    def test_weaklink_u_capacity_saturation(self):
        """r36 P1-3: U(t) must never exceed the chain-end capacity 2*mu0
        mid-integration. Production is throttled by the linear factor
        max(0, 1 - U/(2*mu0)) -- exactly zero at capacity -- so in a
        strong-kt_disp regime U saturates at 2*mu0 from below:
        U(t) ~ cap*(1 - exp(-p0*t/cap)) with p0 = kt_disp*R_ss^2 (sink
        negligible: kia tiny). The channel never touches mu0, so
        cap = 2*mu0(0) = 2 exactly throughout."""
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel(
                ki_A=1.0e-3, kdp_A=1.0e-3, kia_A=1.0e-9,
                ktrec_A=50.0, ktdisp_A=50.0)),
            moments=(1.0, 5.0, 30.0))
        u_slot = rs.qssa_u_slot[0]
        cap = 2.0  # 2 * mu0 amount
        p0 = 50.0 * (1.0e-3 * 4.0 / 100.0)  # kt_disp * R_ss^2 = 2e-3 mol/s
        u_traj = []
        for t in (500.0, 1000.0, 2000.0, 5000.0):
            rs.advance(t)
            y = np.asarray(rs.y)
            assert np.all(np.isfinite(y)), t
            u_traj.append(float(y[u_slot]))
            assert y[u_slot] <= cap + 1.0e-6, t   # NEVER above capacity
        assert all(b > a for a, b in zip(u_traj, u_traj[1:]))
        assert u_traj[-1] > 1.9        # genuinely saturating, not stalled
        # hand check against the closed-form throttled growth (B drifts
        # ~0.8% over 5000 s from the mu1 drain -> few-% tolerance):
        expected = cap * (1.0 - math.exp(-p0 * 5000.0 / cap))
        assert u_traj[-1] == pytest.approx(expected, rel=3e-2)

    def test_weaklink_du_production_scales_with_ktdisp_not_kttotal(self):
        """PIN 4: ONLY disproportionation sources U. Two channels with the
        SAME kt_total (same R_ss, bitwise-same moment drains) but swapped
        rec/disp splits (3:1 vs 1:3) must produce dU/dt in exactly the
        disp-ratio 1:3 -- production scales with kt_disp, never kt_total
        (recombination and random initiation create no U). Fixture
        constraints keeping the ratio EXACT under the r36 law: U0 = 0, so
        u_active = 0 (clamp trivially inactive) and the capacity throttle
        1 - U/(2*mu0) is exactly 1.0 on both systems."""
        rs_lo = self._qssa_m2_system(self._qssa_m2_pool(
            _weaklink_exact_channel(kia_A=1.0e-9,
                                    ktrec_A=3.0e7, ktdisp_A=1.0e7)))
        rs_hi = self._qssa_m2_system(self._qssa_m2_pool(
            _weaklink_exact_channel(kia_A=1.0e-9,
                                    ktrec_A=1.0e7, ktdisp_A=3.0e7)))
        dn_lo = rs_lo.residual(0.0, rs_lo.y0.copy(), np.zeros(7))[0]
        dn_hi = rs_hi.residual(0.0, rs_hi.y0.copy(), np.zeros(7))[0]

        # same kt_total -> bitwise-identical moment drains and emission
        for i in range(6):
            assert dn_hi[i] == dn_lo[i], i
        assert dn_lo[6] > 0.0
        assert dn_hi[6] == pytest.approx(3.0 * dn_lo[6], rel=1e-12)

    # ------------------------------------------------------------------
    # milestone (vi) equivalence pins: inert-by-values / allyl-off /
    # disp-vs-recomb directionality. The validator pins A > 0 (0.0 is
    # ILLEGAL for every Arrhenius block), so the legal minimum for a
    # "switched-off" constant is the smallest positive double 5e-324;
    # these pins lean on two exact float facts, probed before writing:
    # big + 5e-324 == big bitwise (absorbed below the ulp) and
    # 5e-324 * x == 0.0 exactly for x < 0.25 (underflow past the
    # smallest subnormal).
    # ------------------------------------------------------------------

    def test_weaklink_inert_by_values_u_slot_frozen_at_ktdisp_floor(self):
        """PIN P2 (inert-by-values): kt_disp at its legal minimum (A =
        5e-324; the validator rejects A <= 0) with kt_rec == the FULL
        legacy summed kt and U0 = 0 makes the weak-link pool inert-by-
        values in every slot: kt_total = kt_rec + 5e-324 absorbs the
        addend bitwise, so all six core derivatives are bitwise-identical
        to the legacy pool, AND the production arm kt_disp*R_ss^2
        underflows to exactly 0.0 while U0 = 0 zeroes the sink -- the U
        slot is FROZEN (dU/dt == 0.0 exactly), not merely slow. This is
        the deck-authoring recipe for running the 2.2 vocabulary as a
        pure legacy channel."""
        legacy = dict(initiation=_exact_triplet(1.0e-3),
                      depropagation=_exact_triplet(1.0e2),
                      termination=_exact_triplet(1.0e8))
        weak = _weaklink_exact_channel(ki_A=1.0e-3, kdp_A=1.0e2,
                                       kia_A=1.0e3, ktrec_A=1.0e8,
                                       ktdisp_A=5e-324, u0=0.0)
        rs_leg = self._qssa_m2_system(self._qssa_m2_pool(legacy))
        rs_weak = self._qssa_m2_system(self._qssa_m2_pool(weak))

        dn_leg = rs_leg.residual(0.0, rs_leg.y0.copy(), np.zeros(6))[0]
        dn_weak = rs_weak.residual(0.0, rs_weak.y0.copy(), np.zeros(7))[0]

        assert dn_leg[2] < 0.0  # the fixture's channel is live
        for i in range(6):
            assert dn_weak[i] == dn_leg[i], i  # bitwise bridge
        assert dn_weak[6] == 0.0  # U slot frozen exactly, not just tiny

    def test_weaklink_moments_invariant_to_u_at_kia_floor(self):
        """PIN P3 (allyl-off): with ki_allyl at its legal minimum (A =
        5e-324; A = 0 is illegal) U is inert AS A SOURCE -- the allyl
        G_R term f*kia*u_active is absorbed below the ulp of 2*f*ki*B,
        so the moment/species derivatives are BITWISE independent of the
        U value. U itself is still live (the capacity throttle differs
        between the twins), which is exactly the pin: U may grow via
        disproportionation without ever feeding back into the moments."""
        u_on = self._qssa_m2_system(self._qssa_m2_pool(
            _weaklink_exact_channel(kia_A=5e-324, u0=1.0)))
        u_off = self._qssa_m2_system(self._qssa_m2_pool(
            _weaklink_exact_channel(kia_A=5e-324, u0=0.0)))

        dn_on = u_on.residual(0.0, u_on.y0.copy(), np.zeros(7))[0]
        dn_off = u_off.residual(0.0, u_off.y0.copy(), np.zeros(7))[0]

        assert dn_on[2] < 0.0  # live channel, not a vacuous comparison
        for i in range(6):
            assert dn_on[i] == dn_off[i], i  # bitwise: moments blind to U
        # non-vacuity: the U state itself IS different dynamics (capacity
        # throttle 1 - U/(2*mu0): 0.5 vs 1.0), production still positive.
        assert dn_on[6] > 0.0 and dn_off[6] > 0.0
        assert dn_on[6] != dn_off[6]

    def test_weaklink_recombination_never_sources_u(self):
        """PIN P4 (directional, recombination arm): with kt_disp at the
        legal floor (production arm underflows to exactly 0.0) and
        kt_rec = 1e8 carrying ALL the termination, a live channel
        (B > 0, R_ss > 0, throttle > 0) must produce NO U: dU/dt is the
        pure allyl-fission sink -f*kia*u_active*V_poly EXACTLY.
        Recombination does not create unsaturated tail ends. The
        production arm (kt_disp > 0 grows U) is pinned by
        test_weaklink_du_production_scales_with_ktdisp_not_kttotal and
        the bridge test."""
        rs = self._qssa_m2_system(self._qssa_m2_pool(
            _weaklink_exact_channel(kia_A=2.0e4, ktrec_A=1.0e8,
                                    ktdisp_A=5e-324, u0=1.0)))
        dn = rs.residual(0.0, rs.y0.copy(), np.zeros(7))[0]
        assert dn[2] < 0.0   # moments live: radicals exist, chains unzip
        # pure sink, bitwise: production contributes exactly nothing
        # (f = 1, u_active = min(U/V_poly, B) = 1.0, V_poly = 1)
        assert dn[6] == -(1.0 * 2.0e4 * 1.0) * 1.0
        assert dn[6] < 0.0   # U strictly shrinks: no recombination source

    def test_weaklink_u_consumption_is_allyl_fission(self):
        """dU/dt sink = f*ki_allyl(T) * u_active * V_poly (r36: f-symmetric
        with the G_R term, u_active clamped -- inactive here, u < B): pins
        the sink law and its amount/concentration volume convention on a
        sink-dominated exact balance (f = 1)."""
        ki, kia = 1.0e-9, 2.0e4
        ktrec = ktdisp = 0.5e-2  # kt_total = 1e-2
        rs = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel(
                ki_A=ki, kia_A=kia, kdp_A=1.0e-6,
                ktrec_A=ktrec, ktdisp_A=ktdisp)),
            moments=(1.0, 5.0, 30.0))  # B = 4, capacity 2*mu0 = 2
        u_slot = rs.qssa_u_slot[0]
        u = 1.0  # < B (clamp inactive), = cap/2 (throttle = 0.5)
        y = rs.y0.copy()
        y[u_slot] = u
        dn = rs.residual(0.0, y, np.zeros(7))[0]
        # exact balance (production NOT negligible with these constants):
        G_R = 2.0 * ki * 4.0 + kia * u
        R_ss_sq = G_R / (2.0 * 1.0e-2)
        expected = ktdisp * R_ss_sq * (1.0 - u / 2.0) - kia * u
        assert dn[u_slot] == pytest.approx(expected, rel=1e-12)
        assert dn[u_slot] < 0.0  # sink-dominated regime: U shrinks

    def test_weaklink_u_is_massless(self):
        """PIN 5: U is MASSLESS bookkeeping -- it must never enter condensed
        mass/TGA, pool mass, or Mn/Mw/PDI. Structural half: the U slot lives
        OUTSIDE the core-species block (>= n_core), so every mass consumer
        (which reads species amounts + mu slots, all < n_core) is blind to
        it by construction; y0[:n_core] is bitwise identical between U0=0
        and U0=2 systems. Dynamic half: at U huge, the channel still moves
        mass ONLY mu1 -> monomer, exactly balanced (emission == drain), and
        creates/destroys no chains (dmu0 == 0)."""
        rs0 = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel(u0=0.0)))
        rs2 = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_exact_channel(u0=2.0)))
        n_core = rs0.num_core_species
        assert rs0.qssa_u_slot[0] >= n_core
        # identical species/moment state regardless of U0: any mass or
        # Mn/Mw/PDI computed from it is identical.
        assert np.array_equal(rs0.y0[:n_core], rs2.y0[:n_core])
        # Mn/Mw read the mu slots at their unchanged indices:
        assert rs0.y0[2] / rs0.y0[1] == rs2.y0[2] / rs2.y0[1]  # Mn ~ mu1/mu0
        assert rs0.y0[3] / rs0.y0[2] == rs2.y0[3] / rs2.y0[2]  # Mw ~ mu2/mu1

        y_huge = rs0.y0.copy()
        y_huge[rs0.qssa_u_slot[0]] = 1.0e6  # mid-run state, not an IC
        dn = rs0.residual(0.0, y_huge, np.zeros(7))[0]
        assert np.all(np.isfinite(dn))
        assert dn[2] < 0.0                # huge U drives a live drain...
        assert dn[4] == -dn[2]            # ...mass-balanced exactly
        assert dn[1] == 0.0               # no chains created or destroyed
        assert dn[0] == 0.0               # gas untouched

    def test_weaklink_u0_census_trap(self):
        """PIN 6 (r37: TAIL-ONLY basis): at initialization U0 must fit the
        TAIL-DISTRIBUTION chain-end capacity 2*mu0_tail (each chain carries
        at most 2 tail ends; a mol vs mol/L typo must not become a hidden
        initiation source). U is a tail-distribution state -- B and the RHS
        capacity throttle count the tail moments only -- so explicit-species
        ends must NOT back U0 (census and throttle share ONE basis).
        U0 = 3*mu0_tail rejected naming both numbers; U0 == 2*mu0_tail
        accepted; U0 = 0 with mu0 = 0 accepted."""
        # tail mu0 amount = 1.0 mol in the fixture -> capacity = 2.0 mol
        with pytest.raises(ValueError,
                           match=r"poly.*unsaturated_tail_ends_initial=3"
                                 r".*2\*mu0 = 2\b"):
            self._qssa_m2_system(
                self._qssa_m2_pool(_weaklink_channel(
                    unsaturated_tail_ends_initial=3.0)))

        rs_ok = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_channel(
                unsaturated_tail_ends_initial=2.0)))
        assert rs_ok.y0[rs_ok.qssa_u_slot[0]] == 2.0

        rs_empty = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_channel(
                unsaturated_tail_ends_initial=0.0)),
            moments=(0.0, 0.0, 0.0))
        assert rs_empty.y0[rs_empty.qssa_u_slot[0]] == 0.0

    def test_weaklink_u0_census_is_tail_only_not_explicit_backed(self):
        """r37: the census basis must MATCH the RHS throttle basis
        (tail-only). Discriminating fixture: total mu0 = 2 mol of which
        1 mol is explicit DP=2 chains (1 mol P2 at slot 5, subtracted by
        the step-6 consistency pass) -> tail mu0 = 1, tail capacity = 2;
        the OLD (tail+explicit) law would have allowed U0 up to 4.
        U0 = 2.5 sits between the two bounds: it MUST be rejected (a U0
        justified by explicit-species ends would otherwise evolve against
        tail-only capacity -- mixed semantics). U0 = 2.0 (== tail capacity)
        stays accepted in the same fixture."""
        totals = (2.0, 7.0, 34.0)  # tail (1, 5, 30) + explicit (1, 2, 4)
        with pytest.raises(ValueError,
                           match=r"poly.*unsaturated_tail_ends_initial=2\.5"
                                 r".*tail-distribution.*2\*mu0 = 2\b"):
            self._qssa_m2_system(
                self._qssa_m2_pool(_weaklink_channel(
                    unsaturated_tail_ends_initial=2.5), explicit={2: 5}),
                moments=totals, explicit_moles=1.0)

        rs_ok = self._qssa_m2_system(
            self._qssa_m2_pool(_weaklink_channel(
                unsaturated_tail_ends_initial=2.0), explicit={2: 5}),
            moments=totals, explicit_moles=1.0)
        assert rs_ok.y0[rs_ok.qssa_u_slot[0]] == 2.0
        assert rs_ok.y0[1] == 1.0  # step 6 left tail mu0 = 1 (basis proof)
        assert rs_ok.y0[5] == 1.0  # the explicit chains are really there

    def test_weaklink_two_pool_integration_smoke(self):
        """PIN 8: a two-pool system (pool A legacy channel, pool B
        weak-link) integrates without error; the legacy pool gets no U slot
        (state layout [-1, n_core]); U(t) is strictly increasing in the
        hand-checkable production-dominated regime (kia tiny, kt_disp
        production ~ constant while the mu1 drain is negligible), matching
        U(t) ~ kt_disp * R_ss^2 * t."""
        Inert = _spc("N#N", "N2")
        Amu0, Amu1, Amu2 = _spc("CO", "A_mu0"), _spc("C=O", "A_mu1"), _spc("C#N", "A_mu2")
        MA = _spc("C", "MA")
        Bmu0, Bmu1, Bmu2 = _spc("CCO", "B_mu0"), _spc("CC=O", "B_mu1"), _spc("CC#N", "B_mu2")
        MB = _spc("CC", "MB")
        core = [Inert, Amu0, Amu1, Amu2, MA, Bmu0, Bmu1, Bmu2, MB]
        # MA/MB (released-monomer targets) are GAS (2026-07-03 fix)
        mask = np.array([True, False, False, False, True,
                         False, False, False, True], dtype=bool)
        pool_a = PolymerPoolConfig(
            label="A", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
            radical_qssa_unzip=dict(initiation=_exact_triplet(1.0e-3),
                                    depropagation=_exact_triplet(1.0e2),
                                    termination=_exact_triplet(1.0e8)))
        pool_b = PolymerPoolConfig(
            label="B", xs=2, explicit_dp_to_species_index={},
            mu_indices=(5, 6, 7), monomer_poly_index=8,
            radical_qssa_unzip=_weaklink_exact_channel(
                ki_A=1.0e-3, kdp_A=1.0e2, kia_A=1.0e-9,
                ktrec_A=5.0e7, ktdisp_A=5.0e7, u0=0.0))
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0},
            V_poly=1.0, polymer_pools=[pool_a, pool_b], mass_transfer=[],
            gas_species_mask=mask, constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0),
                                     "B": (1.0, 5.0, 30.0)},
            termination=[])
        rs.initialize_model(core, [], [], [])

        assert rs.neq == 10
        assert list(rs.qssa_u_slot) == [-1, 9]
        assert rs.y0[9] == 0.0

        u_traj, mu1_traj = [], []
        for t in (1.0e-3, 1.0e-2, 1.0e-1, 0.5, 1.0):
            rs.advance(t)
            y = np.asarray(rs.y)
            assert np.all(np.isfinite(y)), t
            u_traj.append(float(y[9]))
            mu1_traj.append(float(y[6]))
        assert all(b > a for a, b in zip(u_traj, u_traj[1:]))       # U up
        assert all(b < a for a, b in zip(mu1_traj, mu1_traj[1:]))   # mu1 down
        # hand check: dU/dt ~ kt_disp*R_ss^2 with R_ss = sqrt(ki*B/kt_total)
        R_ss = math.sqrt(1.0e-3 * 4.0 / 1.0e8)
        assert u_traj[-1] == pytest.approx(5.0e7 * R_ss * R_ss * 1.0,
                                           rel=1e-2)

    @staticmethod
    def _qssa_census_system(channel, with_ve):
        """System for the QSSA/scission-VE double-count census: pool proxy
        bound by label, optional surviving same-pool VOLATILE_EJECTION
        (src == dst survives the demotion loop, polymer.pyx ~1531)."""
        Proxy = _spc("CCCC", "poly")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")   # GAS released-monomer target (2026-07-03 fix)
        G = _spc("[CH3]", "G")
        core = [Proxy, Mu0, Mu1, Mu2, M, G]
        mask = np.array([False, False, False, False, True, True], dtype=bool)
        rxns = []
        if with_ve:
            rxn = Reaction(reactants=[Proxy], products=[Proxy, G], **_KIN)
            rxn.polymer_flux_archetype = 6   # VOLATILE_EJECTION
            rxn.polymer_eject_units = 1.0
            rxns = [rxn]
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
            radical_qssa_unzip=channel)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={G: 0.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            termination=[])
        rs.initialize_model(core, rxns, [], [])
        return rs

    def test_qssa_scission_ve_copresence_census_warns_once(self, caplog):
        """Double-count census: QSSA initiation is backbone homolysis, so a
        pool carrying BOTH the channel AND a surviving generated scission/VE
        reaction sourced from it may double-count the same physics. That
        must be censused (warn-once, NEVER refuse -- generated scission may
        cover different bonds)."""
        import rmgpy.solver.polymer as sp
        # Warn-once is process-global and keyed (pool, overlap kind); an
        # earlier fixture in this file may legitimately warm this key, so
        # reset exactly the key under test to make the assertion
        # deterministic (mirrors the k_scission co-presence twin).
        from rmgpy.polymer import _qssa_double_count_warned
        _qssa_double_count_warned.discard(("poly", "generated_scission_ve"))
        with caplog.at_level(logging.WARNING):
            rs = self._qssa_census_system(_qssa_channel(), with_ve=True)
            # The VE survived demotion (src == dst same-pool shape).
            assert rs.reaction_flux_archetype[0] == sp.FLUX_VOLATILE_EJECTION
            assert rs.reaction_src_pool[0] == 0
            assert len(rs.qssa_double_count_census) == 1
            assert rs.qssa_double_count_census[0]["pool"] == "poly"
            # Warn-once: a SECOND build censuses again but does not re-log.
            rs2 = self._qssa_census_system(_qssa_channel(), with_ve=True)
            assert len(rs2.qssa_double_count_census) == 1
        qssa_warnings = [rec for rec in caplog.records
                         if "QSSA" in rec.getMessage()
                         and "poly" in rec.getMessage()
                         and "double-count" in rec.getMessage().lower()]
        assert len(qssa_warnings) == 1

    def test_qssa_only_pool_census_silent(self):
        """A QSSA pool with NO scission/VE reactions must NOT be censused;
        neither must a scission/VE pool without the channel."""
        rs_qssa_only = self._qssa_census_system(_qssa_channel(),
                                                with_ve=False)
        assert rs_qssa_only.qssa_double_count_census == []
        rs_ve_only = self._qssa_census_system(None, with_ve=True)
        assert rs_ve_only.qssa_double_count_census == []

    def test_qssa_disabled_row_arrhenius_garbage_is_inert(self):
        """The RHS must gate on qssa_enabled, NEVER on A != 0: poisoning a
        disabled pool's flattened Arrhenius slots with live-looking values
        must leave the residual bitwise unchanged (the channel-off path is
        untouched by M2)."""
        rs = self._qssa_m2_system(self._qssa_m2_pool(None, k_unzip=0.1))
        dn_before = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0].copy()
        assert np.any(dn_before != 0.0)  # k_unzip keeps the fixture live

        assert rs.qssa_enabled[0] == 0
        rs.qssa_ki_A[0] = 1.0e15
        rs.qssa_ki_Ea[0] = 3.0e5
        rs.qssa_kdp_A[0] = 1.0e13
        rs.qssa_kdp_Ea[0] = 8.0e4
        rs.qssa_kt_A[0] = 1.0e8
        rs.qssa_kt_Ea[0] = 1.0e4
        rs.qssa_has_transfer[0] = 1
        rs.qssa_ktr_A[0] = 1.0e5

        dn_after = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.array_equal(dn_after, dn_before)

    def test_qssa_runtime_guard_kt_underflow_raises(self):
        """Runtime degenerate-rate guard: a huge-but-FINITE termination Ea
        passes M1 config validation but underflows kt(T) to exactly 0.0 at
        the solver temperature -- R_ss = sqrt(f*ki*B/kt) would divide by
        zero into inf. The RHS must raise ValueError naming the pool, T and
        the termination constant (fail-loud), never emit inf/NaN flux."""
        # exp(-5e6 / (8.314*800)) = exp(-751.7) -> 0.0 (double underflow)
        channel = _qssa_channel(termination=_qssa_triplet(A=1.0e8, Ea=5.0e6))
        with pytest.raises(ValueError, match=r"poly.*800.*termination"):
            # initialize_model evaluates the residual once at setup, so the
            # guard may fire at construction (even better: t=0 fail-loud);
            # if it survives construction the direct RHS call must raise.
            rs = self._qssa_m2_system(self._qssa_m2_pool(channel))
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))

    def test_qssa_runtime_guard_overflow_to_inf_raises(self):
        """Runtime degenerate-rate guard: an initiation block finite in
        config (A=1e300, n=10 -- every field passes the M1 finiteness rules)
        overflows to A*T**n = inf at evaluation. The RHS must raise
        ValueError naming the pool, T and the initiation constant instead
        of propagating inf through R_ss into the residual."""
        channel = _qssa_channel(
            initiation=_qssa_triplet(A=1.0e300, n=10.0, Ea=0.0))
        with pytest.raises(ValueError, match=r"poly.*800.*initiation"):
            rs = self._qssa_m2_system(self._qssa_m2_pool(channel))
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))

    def test_qssa_k_scission_copresence_census_warns_once(self, caplog):
        """Census gap closure: QSSA initiation and the pool's own k_scission
        BOTH represent random backbone homolysis -- the most direct
        initiation double-count. A pool with the channel AND k_scission > 0
        must be censused (warn-once, never refuse), mirroring the
        archetype-reaction co-presence census."""
        # Warn-once is process-global and keyed (pool, overlap kind); an
        # earlier fixture in this file (the M2 footprint test: 'poly' with
        # k_scission=0.3 + channel) legitimately warms this key, so reset
        # exactly the key under test to make the assertion deterministic.
        from rmgpy.polymer import _qssa_double_count_warned
        _qssa_double_count_warned.discard(("poly", "k_scission"))
        with caplog.at_level(logging.WARNING):
            rs = self._qssa_m2_system(
                self._qssa_m2_pool(_qssa_channel(), k_scission=0.02))
            ks_entries = [e for e in rs.qssa_double_count_census
                          if e.get("overlap") == "k_scission"]
            assert len(ks_entries) == 1
            assert ks_entries[0]["pool"] == "poly"
            assert ks_entries[0]["k_scission"] == pytest.approx(0.02)
            # Warn-once: a SECOND build censuses again but does not re-log.
            rs2 = self._qssa_m2_system(
                self._qssa_m2_pool(_qssa_channel(), k_scission=0.02))
            assert len([e for e in rs2.qssa_double_count_census
                        if e.get("overlap") == "k_scission"]) == 1
        ks_warnings = [rec for rec in caplog.records
                       if "QSSA" in rec.getMessage()
                       and "k_scission" in rec.getMessage()
                       and "poly" in rec.getMessage()]
        assert len(ks_warnings) == 1

    def test_qssa_without_k_scission_census_silent(self):
        """A QSSA pool with k_scission == 0 must NOT produce a k_scission
        census entry; a k_scission-only pool without the channel must stay
        entirely out of the QSSA census."""
        rs_qssa_only = self._qssa_m2_system(
            self._qssa_m2_pool(_qssa_channel(), k_scission=0.0))
        assert [e for e in rs_qssa_only.qssa_double_count_census
                if e.get("overlap") == "k_scission"] == []
        rs_ks_only = self._qssa_m2_system(
            self._qssa_m2_pool(None, k_scission=0.02))
        assert rs_ks_only.qssa_double_count_census == []

    # ------------------------------------------------------------------
    # radical_qssa_unzip daughter-pool inheritance (M5): a scission event
    # spawns a daughter that must KEEP depolymerizing, or the TGA S-curve
    # freezes after the parent's first generation.
    # ------------------------------------------------------------------

    @staticmethod
    def _qssa_scission_daughter():
        """Produce a daughter Polymer through a REAL scission event on a QSSA
        parent (create_reacted_copy on a head-wing + labeled-methyl fragment),
        exactly the production path that registers scission daughters."""
        from rmgpy.polymer import Polymer, stitch_molecules_by_labeled_atoms

        channel = _qssa_channel()
        parent = Polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                         end_groups=["[CH3]", "[H]"], cutoff=2,
                         Mn=5000.0, Mw=6000.0, initial_mass=1.0,
                         radical_qssa_unzip=channel)
        styrene = _spc("C=Cc1ccccc1", "styrene")
        parent.monomer_product_species = styrene

        methyl_star2 = Molecule().from_adjacency_list("""multiplicity 2
            1 *2 C u1 p0 c0 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 H u0 p0 c0 {1,S}""")
        frag = stitch_molecules_by_labeled_atoms(
            parent._stitch_wing("head"), methyl_star2)
        daughter = parent.create_reacted_copy(frag)
        assert daughter is not None
        assert daughter.label == "PS_scission_tail"
        return daughter, styrene, channel

    def _qssa_daughter_system(self, daughter, styrene, moments=(1.0, 5.0, 30.0),
                              reactions=None):
        """Solver built the production way for a spawned daughter: pool config
        comes from derive_daughter_pool_configs over the registered core
        species (daughter + _muN dummies), NOT hand-built -- the exact path
        HybridPolymerReactor.to_solver_object runs on every rebuild."""
        from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

        Inert = _spc("N#N", "N2")
        d_mu0 = _spc("CO", "PS_scission_tail_mu0")
        d_mu1 = _spc("C=O", "PS_scission_tail_mu1")
        d_mu2 = _spc("C#N", "PS_scission_tail_mu2")
        core = [Inert, daughter, d_mu0, d_mu1, d_mu2, styrene]
        spc_map = {s: i for i, s in enumerate(core)}

        cfgs = derive_daughter_pool_configs(core, spc_map,
                                            existing_pool_labels={"PS"})
        assert len(cfgs) == 1 and cfgs[0].label == "PS_scission_tail"

        # styrene (the routed monomer_product at index 5) is GAS
        # (2026-07-03 monomer-gas fix); the daughter pool members stay
        # condensed.
        mask = np.array([True, False, False, False, False, True], dtype=bool)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0},
            V_poly=1.0, polymer_pools=[cfgs[0]], mass_transfer=[],
            gas_species_mask=mask, constant_gas_volume=False,
            initial_polymer_moments={"PS_scission_tail": tuple(moments)},
            termination=[])
        rs.initialize_model(core, list(reactions or []), [], [])
        return rs, cfgs[0]

    def test_spawned_scission_daughter_ejects_monomer_through_inherited_qssa(self):
        """M5 integration: a scission EVENT produces a daughter Polymer; the
        daughter's derived pool config carries the inherited channel + monomer
        routing; a fresh solver build FLATTENS it (qssa_enabled -- the only
        signal the RHS trusts, killing the enabled-in-dict/absent-in-arrays
        failure) and the RHS actually DRAINS the daughter's mu1 and EMITS
        monomer at the analytic QSSA rate. Companion: stripping the channel
        from the daughter yields a solver whose residual shows NO channel
        contribution (channel-free daughters stay honest-inert)."""
        daughter, styrene, channel = self._qssa_scission_daughter()
        mu0, mu1, mu2 = 1.0, 5.0, 30.0

        rs_on, cfg_on = self._qssa_daughter_system(daughter, styrene,
                                                   moments=(mu0, mu1, mu2))
        assert cfg_on.radical_qssa_unzip is not None
        assert cfg_on.monomer_poly_index == 5
        assert rs_on.qssa_enabled[0] == 1

        # Channel-free twin: same daughter with the channel stripped.
        daughter_off = daughter.copy()
        daughter_off.radical_qssa_unzip = None
        rs_off, cfg_off = self._qssa_daughter_system(daughter_off, styrene,
                                                     moments=(mu0, mu1, mu2))
        assert cfg_off.radical_qssa_unzip is None
        assert rs_off.qssa_enabled[0] == 0

        y = rs_on.y.copy()
        dn_on = rs_on.residual(0.0, y, np.zeros_like(y))[0]
        dn_off = rs_off.residual(0.0, y.copy(), np.zeros_like(y))[0]
        diff = dn_on - dn_off

        r = self._qssa_oracle_rate(channel, mu0, mu1)
        assert r > 0.0                                   # the oracle is live
        assert diff[3] == pytest.approx(-r, rel=1e-10)   # daughter mu1 drains
        assert diff[5] == pytest.approx(+r, rel=1e-10)   # monomer emitted

    def test_spawned_daughter_inherited_qssa_seen_by_double_count_census(self):
        """The M2 double-count census keys on qssa_enabled, which is rebuilt
        (with the flattening) on EVERY initialize_model -- so a daughter that
        inherits the channel is censused exactly like a deck pool when a
        surviving same-pool scission/VE reaction sources from it."""
        import rmgpy.solver.polymer as sp
        from rmgpy.polymer import _qssa_double_count_warned
        _qssa_double_count_warned.discard(
            ("PS_scission_tail", "generated_scission_ve"))

        daughter, styrene, _ = self._qssa_scission_daughter()
        G = _spc("[CH3]", "G_ve")
        rxn = Reaction(reactants=[daughter], products=[daughter, G], **_KIN)
        rxn.polymer_flux_archetype = 6   # VOLATILE_EJECTION (same-pool)
        rxn.polymer_eject_units = 1.0

        # G must be in core for the reaction indices to resolve; extend the
        # daughter system's core via the reactions hook.
        from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
        Inert = _spc("N#N", "N2")
        d_mu0 = _spc("CO", "PS_scission_tail_mu0")
        d_mu1 = _spc("C=O", "PS_scission_tail_mu1")
        d_mu2 = _spc("C#N", "PS_scission_tail_mu2")
        core = [Inert, daughter, d_mu0, d_mu1, d_mu2, styrene, G]
        spc_map = {s: i for i, s in enumerate(core)}
        cfgs = derive_daughter_pool_configs(core, spc_map,
                                            existing_pool_labels={"PS"})
        # styrene (routed monomer_product at 5) GAS -- 2026-07-03 fix
        mask = np.array([True, False, False, False, False, True, True],
                        dtype=bool)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0},
            V_poly=1.0, polymer_pools=[cfgs[0]], mass_transfer=[],
            gas_species_mask=mask, constant_gas_volume=False,
            initial_polymer_moments={"PS_scission_tail": (1.0, 5.0, 30.0)},
            termination=[])
        rs.initialize_model(core, [rxn], [], [])

        assert rs.reaction_flux_archetype[0] == sp.FLUX_VOLATILE_EJECTION
        assert rs.reaction_src_pool[0] == 0
        entries = [e for e in rs.qssa_double_count_census
                   if e.get("overlap") == "generated_scission_ve"]
        assert len(entries) == 1
        assert entries[0]["pool"] == "PS_scission_tail"

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

    # ------------------------------------------------------------------
    # VOLATILE_EJECTION (archetype 6): a MIGRATION mirror that debits `a`
    # (= polymer_eject_units, Sigma volatile MW / source monomer MW) on the
    # destination (to-pool) leg. The from-leg loses the FULL bundle
    # (identical to MIGRATION); the to-leg gains the SHIFTED bundle, so the
    # net condensed mu1 change is -a*ev forward / +a*ev reverse (backbone
    # mass ejected as / restored from the discrete gas volatile, which is
    # itself credited by the standard section-4 net-rate product path).
    # ------------------------------------------------------------------
    def test_volatile_ejection_forward_one_event_shifts_dest_leg(self):
        """
        Forward VE one-event: src (A) loses its full length-biased bundle;
        dst (B) gains the a-shifted bundle. Net condensed Delta mu1 == -a*ev
        (mass leaves the melt). mu2 legs: src -ev*E[n^2]; dst
        +ev*(E[n^2] - 2a*E[n] + a^2).
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        ev = kf * mu1a                       # site-scaled by A's mu1, V_poly=1
        mu3a = mu0a * (mu2a / mu1a) ** 3     # 216.0
        b1 = mu2a / mu1a                     # E[n]  = 6.0
        b2 = mu3a / mu1a                     # E[n^2]= 43.2

        # src (A) loses the FULL bundle
        assert np.isclose(dn_dt[1], -ev * 1.0)                 # A mu0
        assert np.isclose(dn_dt[2], -ev * b1)                  # A mu1
        assert np.isclose(dn_dt[3], -ev * b2)                  # A mu2 == -ev*E[n^2]
        # dst (B) gains the a-SHIFTED bundle
        assert np.isclose(dn_dt[5], +ev * 1.0)                 # B mu0
        assert np.isclose(dn_dt[6], +ev * (b1 - a))            # B mu1 == ev*(E[n]-a)
        assert np.isclose(dn_dt[7], +ev * (b2 - 2.0 * a * b1 + a * a))  # B mu2
        # net condensed conservation: mu0 conserved, mu1 loses exactly a*ev
        assert np.isclose(dn_dt[1] + dn_dt[5], 0.0, atol=1e-12)
        assert np.isclose(dn_dt[2] + dn_dt[6], -a * ev)

    def test_volatile_ejection_reverse_one_event_restores_parent(self):
        """
        Reverse VE one-event: dst (B) loses its full bundle; src (A, the
        re-formed parent) gains the +a-shifted bundle. Net condensed Delta
        mu1 == +a*ev (parent backbone mass restored). Reverse-only isolation:
        kf zeroed, kb driven directly (mirrors the MIGRATION reverse harness).
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        mu_a = (1.0, 5.0, 30.0)
        mu_b = (2.0, 4.0, 10.0)
        rs = _two_pool_rs(rxn, core, mask, mu_a, mu_b)
        rs.kf[0] = 0.0        # forward off: isolate the reverse leg
        rs.kb[0] = 0.6

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        # rr = kb * C(proxyB=1), then site-scaled by the REACTANT-pool (A) mu1
        # (mirrors MIGRATION: reverse site scaling keys on the reactant pool;
        # the deferred reverse-site question is out of scope here).
        ev = 0.6 * mu_a[1]                              # 3.0
        mu3_b = mu_b[0] * (mu_b[2] / mu_b[1]) ** 3      # 31.25
        bB1 = mu_b[2] / mu_b[1]                         # E[n_dst] = 2.5
        bB2 = mu3_b / mu_b[1]                           # E[n_dst^2] = 7.8125

        # dst (B) loses its FULL bundle
        assert np.isclose(dn_dt[5], -ev * 1.0)                 # B mu0
        assert np.isclose(dn_dt[6], -ev * bB1)                 # B mu1
        assert np.isclose(dn_dt[7], -ev * bB2)                 # B mu2
        # src (A, re-formed parent) gains the +a-SHIFTED bundle
        assert np.isclose(dn_dt[1], +ev * 1.0)                 # A mu0
        assert np.isclose(dn_dt[2], +ev * (bB1 + a))           # A mu1 == ev*(E[n_dst]+a)
        assert np.isclose(dn_dt[3], +ev * (bB2 + 2.0 * a * bB1 + a * a))  # A mu2
        # net condensed conservation: mu1 GAINS exactly a*ev (parent restored)
        assert np.isclose(dn_dt[1] + dn_dt[5], 0.0, atol=1e-12)
        assert np.isclose(dn_dt[2] + dn_dt[6], +a * ev)

    def test_volatile_ejection_fractional_a_not_rounded(self):
        """
        The stamped fractional a (e.g. 1.135 for alpha-methylstyrene off a
        styrene pool) must flow verbatim -- NOT rounded to 1. The net
        condensed mu1 change must equal -1.135*ev, distinguishable from the
        integer-rounded -1.0*ev.
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        ev = kf * 5.0
        net_mu1 = dn_dt[2] + dn_dt[6]
        assert np.isclose(net_mu1, -a * ev)          # exact fractional a
        assert not np.isclose(net_mu1, -1.0 * ev)    # NOT rounded to 1

    def test_migration_unchanged_by_volatile_ejection_dispatch(self):
        """
        Regression guard: an equivalent MIGRATION reaction (archetype 2, no
        eject_units) still conserves mu1 across pools (net Delta mu1 == 0).
        The VE branch must not leak an `a`-shift into MIGRATION.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2   # MIGRATION, eject_units defaults 0.0
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.isclose(dn_dt[1] + dn_dt[5], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[2] + dn_dt[6], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[3] + dn_dt[7], 0.0, atol=1e-14)

    # ------------------------------------------------------------------
    # SAME-POOL volatile ejection (signed a, spec 2026-06-2x round-13):
    # unzip / depropagation write on ONE pool (src == dst). This is a
    # DIRECT chip-style drain, NOT the cross-pool per-direction cancellation
    # (which can push mu2 the wrong way near exhaustion). a > 0 sheds mass
    # (mu1 drops by a*ev); a < 0 GAINS mass (mu1 rises by |a|*ev, chain
    # growth). The mu2 decrement is CLAMPED at > 0 exactly like DISCRETE_CHIP.
    # ------------------------------------------------------------------
    def test_same_pool_ve_forward_drains_one_pool(self):
        """
        Same-pool VE (src == dst) forward, a > 0: a single pool's mu1 drops
        by a*ev; mu2 by ev*(2a*E[n] - a^2) (clamp holds, decrement > 0).
        No second pool is touched (B moments stay zero) -- distinct from the
        cross-pool VE which splits the write across A and B.
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        ev = kf * mu1a                       # site-scaled by A's mu1, V_poly=1
        b1 = mu2a / mu1a                     # E[n] = 6.0 (length-biased pick)
        mu2_dec = 2.0 * a * b1 - a * a       # 12.331775 > 0 -> clamp holds
        assert mu2_dec > 0.0
        assert np.isclose(dn_dt[2], -a * ev)             # A mu1 drained by a*ev
        assert np.isclose(dn_dt[3], -ev * mu2_dec)       # A mu2 clamped decrement
        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)     # mu0 unchanged (no chain lost)
        # the OTHER pool (B) is untouched -- this is a single-pool write
        assert np.allclose(dn_dt[5:8], 0.0, atol=1e-14)
        # net condensed mu1 change == -a*ev (mass ejected to the gas volatile)
        assert np.isclose(dn_dt[2] + dn_dt[6], -a * ev)

    def test_same_pool_ve_reverse_restores_one_pool(self):
        """
        Same-pool VE reverse, a > 0: mu1 rises by a*ev on the single pool;
        mu2 by ev*(2a*E[n] + a^2) -- the EXACT extension form (no clamp,
        unconditionally positive). Reverse-only isolation: kf zeroed, kb
        driven directly (mirrors the cross-pool VE reverse harness).
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))
        rs.kf[0] = 0.0        # forward off: isolate the reverse leg
        rs.kb[0] = 0.6

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        ev = 0.6 * mu1a                      # rr = kb*C(proxyA=1), site-scaled by mu1
        b1 = mu2a / mu1a                     # E[n] = 6.0
        assert np.isclose(dn_dt[2], +a * ev)                       # A mu1 restored
        assert np.isclose(dn_dt[3], +ev * (2.0 * a * b1 + a * a))  # A mu2 exact extension
        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)               # mu0 unchanged
        assert np.allclose(dn_dt[5:8], 0.0, atol=1e-14)            # B untouched
        assert np.isclose(dn_dt[2] + dn_dt[6], +a * ev)            # net mu1 GAINS a*ev

    def test_same_pool_ve_negative_a_adds_mass(self):
        """
        Signed a < 0 (monomer/radical ADDITION back onto the same pool):
        forward mu1 -= rf_mol*a with a < 0 ADDS mass -> mu1 RISES by |a|*ev
        (chain growth). The mu2 decrement 2a*E[n] - a^2 is negative for a < 0,
        so the `> 0` clamp SKIPS the mu2 write (documented approximation: the
        forward same-pool leg does not credit mu2 growth for a < 0; the exact
        +extension lives on the reverse leg). mu1 growth is the load-bearing
        signed behavior and is asserted exactly.
        """
        a = -1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        ev = kf * mu1a
        b1 = mu2a / mu1a
        assert 2.0 * a * b1 - a * a < 0.0                # decrement negative -> clamp
        assert dn_dt[2] > 0.0                            # mu1 RISES (mass gained)
        assert np.isclose(dn_dt[2], -a * ev)             # == +|a|*ev
        assert np.isclose(dn_dt[3], 0.0, atol=1e-14)     # mu2 write clamped out
        assert np.allclose(dn_dt[5:8], 0.0, atol=1e-14)  # B untouched

    def test_same_pool_ve_throttle_direct_rhs_regimes(self):
        """
        Exhaustion throttle (spec s5 amendment) extended to same-pool a > 0
        VE: a mu0-scaled (end-group) same-pool VE drains mu1 but never mu0,
        so unthrottled mu1 would run linearly negative past exhaustion.
        Throttled regime (mu1 < a*mu0): site = mu1/a, dmu1/dt = -kf*mu1
        EXACTLY. Healthy regime (mu1 >> a*mu0): site = mu0 (byte-identical to
        the pre-throttle path). Mirrors the DISCRETE_CHIP direct-RHS test.
        """
        a = 3.0
        for mom_a, site in (((1.0, 1.0, 2.0), 1.0 / 3.0),    # throttled: min(1, 1/3)
                            ((1.0, 50.0, 3000.0), 1.0)):     # healthy:  min(1, 50/3)
            sp, core, mask = _two_pool_species()
            rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
            rxn.polymer_flux_archetype = 6
            rxn.polymer_eject_units = a
            rxn.is_end_group_reaction = True
            rs = _two_pool_rs(rxn, core, mask, mom_a, (0.1, 0.2, 0.5))

            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            kf = rxn.get_rate_coefficient(800.0, 1.0e5)
            assert np.isclose(dn_dt[2], -kf * site * a), mom_a   # mu1 drain throttled
            if site < mom_a[0]:
                # Throttled identity: dmu1/dt = -kf*mu1 exactly.
                assert np.isclose(dn_dt[2], -kf * mom_a[1])
                assert not np.isclose(dn_dt[2], -kf * mom_a[0] * a)

    def test_same_pool_ve_negative_a_exempt_from_throttle(self):
        """
        a < 0 same-pool VE GROWS the chain (mu1 rises), so there is no
        exhaustion to throttle: the guard is a > 0. The rate must stay the
        unthrottled mu0-scaled value even when mu1 is small relative to mu0
        (no mu1/a division, which for a < 0 would be a spurious negative site).
        """
        a = -3.0
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        rxn.is_end_group_reaction = True
        mom_a = (2.0, 1.0, 1.0)               # mu1 small vs mu0 -- would throttle if a>0
        rs = _two_pool_rs(rxn, core, mask, mom_a, (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        assert np.all(np.isfinite(dn_dt))
        # site == mu0 (unthrottled), so mu1 drain leg = -kf*mu0*a = +kf*mu0*|a|
        assert np.isclose(dn_dt[2], -kf * mom_a[0] * a)   # mu1 RISES, unthrottled

    def test_same_pool_ve_exhaustion_trajectory_no_negative(self):
        """
        Same-pool a > 0 VE forward-Euler trajectory past exhaustion: with the
        throttle, mu1 decays at worst exponentially and never crosses zero,
        the FULL cone holds, and Sum(mu1) + a*n_gas is conserved at every step
        (each ejected unit leaves the melt as gas mass). Mirrors the
        DISCRETE_CHIP exhaustion trajectory.
        """
        a = 3.0
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        rxn.is_end_group_reaction = True
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.0, 0.0, 0.0))

        y = rs.y.copy()
        invariant0 = y[2] + a * y[8]          # Sum(mu1) + a*n_gas
        dt = 0.005
        gas_increments = []
        for _ in range(800):                  # t = 4 s, well past unthrottled zero-cross
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
            gas_increments.append(dt * dn_dt[8])
            assert y[2] >= 0.0                              # mu1 never crosses zero
            if y[2] > 1e-12:
                assert y[1] * y[3] >= y[2] ** 2 * (1.0 - 1e-9)   # full cone
        assert gas_increments[-1] < 1e-2 * gas_increments[0]     # production decays
        assert np.isclose(y[2] + a * y[8], invariant0, rtol=1e-9, atol=1e-12)

    def test_same_pool_ve_throttle_diagnostic_rate_parity(self):
        """
        Diagnostic-path parity: get_reaction_rates (the [THE HIJACK] block)
        must apply the SAME same-pool VE exhaustion throttle as the residual.
        mu0=1, mu1=1, a=3 is throttled (mu1/a = 1/3 < mu0), so the diagnostic
        rate is kf*min(mu0, mu1/a) = kf/3, NOT the unthrottled kf*mu0, and
        equals what the residual's dmu1/(-a) implies.
        """
        a = 3.0
        mom_a = (1.0, 1.0, 2.0)               # throttled: min(1, 1/3) = 1/3
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        rxn.is_end_group_reaction = True
        rs = _two_pool_rs(rxn, core, mask, mom_a, (0.1, 0.2, 0.5))

        rate = rs.get_reaction_rates(rs.y)[0]
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        site = min(mom_a[0], mom_a[1] / a)
        assert np.isclose(rate, kf * site)                # throttled, not kf*mu0
        assert not np.isclose(rate, kf * mom_a[0])
        assert np.isclose(rate, dn_dt[2] / (-a))          # parity with residual

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

    def test_gas_volatile_veto_excludes_from_melt_tag_branch(self, caplog):
        """Spec §5.1 (durable gas-volatile veto): a genuine discrete gas
        volatile (alpha-methylstyrene, C9H10, ~118 g/mol) reaches the solver
        gas-masked (gas_mask=True) but proxy-TAGGED and chain-scale, so the
        melt-gate TAG branch (proxy AND MW>=window) would wrongly count it as a
        melt reference-state participant. Because is_polymer_proxy is a
        monotonic multi-writer sticky cache with no gas authority, the durable
        verdict lives in ``Species.props['polymer_reference_state_gas_veto']``
        (set once at the discrete-product creation point, copied by
        Species.copy, never touched by the proxy stamping machinery). The gate
        must honor it: TAG branch = proxy AND MW>=window AND NOT veto.

        Live shape (the PS-run crash): melt reactant [chain] vs melt products
        [volatile, tail] -> dn_melt=+1 -> one uncancelled Sackur-Tetrode term
        -> U~11 decades -> REFUSE. With the volatile vetoed it drops out of the
        product melt set -> dn_melt=0 -> benign, builds.

        RED before the gate honors the veto: the vetoed build still REFUSES
        (the veto is ignored). The no-veto control refuses in BOTH pre/post
        (it is the liveness pin proving the scenario genuinely refuses).
        """
        import logging
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO_KEY

        def _build(set_veto):
            # one-pool core: melt proxy A + mu0/1/2 dummies (indices 0-3), a
            # condensed melt chain TAIL (index 4, gas_mask=False), and the
            # discrete gas volatile VOL (index 5, gas_mask=True, proxy-TAGGED,
            # chain-scale MW 0.118 >= window -> only the TAG branch can make it
            # melt).
            a = _spc("CCCC", "A")
            a_mu0, a_mu1, a_mu2 = (_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"),
                                   _spc("C#N", "A_mu2"))
            tail = _spc("CCCCCCCC", "TAIL")
            vol = _spc("C=C(C)c1ccccc1", "VOL")
            vol.is_polymer_proxy = True
            if set_veto:
                vol.props[VETO_KEY] = True
            for s in (a, a_mu0, a_mu1, a_mu2, tail, vol):
                s.thermo = _trivial_nasa(_GAV_COMMENT)
            core = [a, a_mu0, a_mu1, a_mu2, tail, vol]
            mask = np.array([False, False, False, False, False, True], dtype=bool)
            rxn = Reaction(reactants=[a], products=[vol, tail], **_REV_KIN)
            # pool config -> chain_window = (104 + 10)/1000 = 0.114 kg/mol,
            # just below the volatile's 0.118 so the TAG branch admits it.
            return _refstate_rs(core, [rxn], mask,
                                [_gate_pool_config(104.0)], {"A": (1.0, 5.0, 30.0)})

        # LIVENESS PIN: without the veto the scenario genuinely refuses
        # (proves the volatile-as-melt asymmetry is the cause).
        with pytest.raises(ValueError, match="unpaired reference-state"):
            _build(set_veto=False)

        # THE red assertion: WITH the durable veto the volatile is excluded
        # from the melt sum -> the build SUCCEEDS (no refusal).
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            rs = _build(set_veto=True)
        assert rs.reference_state_max_decades < 3.0, (
            "durable gas-volatile veto must exclude the volatile from the melt "
            "tag branch so the reference-state term stays paired/benign"
        )

    def test_gas_veto_suppression_of_chain_scale_member_is_logged(self, caplog):
        """Backstop (code-review IMPORTANT #1): the veto silences the
        reference-state tripwire for a CHAIN-SCALE (MW >= window) product,
        because create_reacted_copy returns None both for a genuine gas
        volatile AND for a wing-match FAILURE on a real chain-scale fragment.
        The gate cannot distinguish them (alpha-MS is itself above the window),
        so a handshake false-None on a genuine chain would be silently dropped
        from the melt sum instead of loudly refused. To preserve the backstop,
        the gate must emit a visible census/warning naming any chain-scale
        member whose melt classification the veto suppressed -- so a human can
        catch a mis-vetoed real chain even though the build no longer refuses.

        RED before the census logic: no such warning is emitted.
        """
        import logging
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO_KEY

        a = _spc("CCCC", "A")
        a_mu0, a_mu1, a_mu2 = (_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"),
                               _spc("C#N", "A_mu2"))
        tail = _spc("CCCCCCCC", "TAIL")
        vol = _spc("C=C(C)c1ccccc1", "VOL")       # chain-scale MW 0.118 >= window
        vol.is_polymer_proxy = True
        vol.props[VETO_KEY] = True
        for s in (a, a_mu0, a_mu1, a_mu2, tail, vol):
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [a, a_mu0, a_mu1, a_mu2, tail, vol]
        mask = np.array([False, False, False, False, False, True], dtype=bool)
        rxn = Reaction(reactants=[a], products=[vol, tail], **_REV_KIN)

        with caplog.at_level(logging.WARNING):
            rs = _refstate_rs(core, [rxn], mask,
                              [_gate_pool_config(104.0)], {"A": (1.0, 5.0, 30.0)})

        # a visible notice naming the suppressed chain-scale species
        assert any("GAS VETO" in r.getMessage() and "VOL" in r.getMessage()
                   for r in caplog.records), (
            "veto suppression of a chain-scale melt member must be logged so a "
            "handshake false-None on a genuine chain surfaces for a human"
        )
        # and the census is recorded on the engine for programmatic inspection
        assert any("VOL" in str(entry)
                   for entry in getattr(rs, "gas_veto_census", [])), (
            "the suppressed chain-scale member must be recorded in gas_veto_census"
        )

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


class TestDaughterPoolRegistration:
    """Stage 1 (proxy_reaction_reality_rules.md Layer 2): a scission/spawn
    daughter Polymer registered as a core species (with its _mu0/_mu1/_mu2
    dummies) must be auto-registered as a solver pool by to_solver_object, so a
    stamped SCISSION_FRAGMENT/MIGRATION reaction targeting it resolves its pool
    instead of demoting to UNRESOLVED ("could not resolve their solver pool(s)").
    The deck list polymerPhase.pools never grows at runtime, so to_solver_object
    derives the daughter config from the live core species themselves."""

    @staticmethod
    def _moment_dummy(label):
        s = _spc("[Ne]", label)
        s.reactive = False
        s.is_moment_dummy = True
        return s

    def test_to_solver_object_registers_core_daughter_as_pool(self):
        """STAGE 1 / CYCLE 3 (integration): a daughter Polymer in core_species
        but absent from the static deck pool list must appear as a solver pool
        after to_solver_object."""
        from rmgpy.quantity import Quantity
        from rmgpy.polymer import Polymer
        from rmgpy.rmg.polymer_input import HybridPolymerReactor, PolymerPhase

        a = _spc("CCCC", "A")
        daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                           end_groups=["[CH3]", "[H]"], cutoff=3,
                           Mn=5000.0, Mw=6000.0, initial_mass=0.001)
        mu0 = self._moment_dummy("PS_d1_mu0")
        mu1 = self._moment_dummy("PS_d1_mu1")
        mu2 = self._moment_dummy("PS_d1_mu2")
        core = [a, daughter, mu0, mu1, mu2]

        phase = PolymerPhase(density=Quantity(1050.0, "kg/m^3"), initial_moments={},
                             initial_explicit={a: 1.0}, pools=[], mass_transfer=[])
        reactor = HybridPolymerReactor(
            temperature=(1000.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles={a: 1.0}, polymerPhase=phase, terminationTime=(0.1, "s"))

        solver = reactor.to_solver_object(core, [], [], [])

        labels = {p.label for p in solver.polymer_pools}
        assert "PS_d1" in labels

    def test_derived_daughter_config_resolves_scission_instead_of_demoting(self):
        """STAGE 1 / CYCLE 4 (payoff): the positive counterpart to
        test_stamped_scission_without_daughter_pool_demotes_to_legacy. A stamped
        SCISSION_FRAGMENT parent->daughter reaction whose daughter pool is
        DERIVED by derive_daughter_pool_configs must resolve BOTH src and dst and
        stay SCISSION_FRAGMENT -- not demote to UNRESOLVED. Closes the '2
        reactions could not resolve their solver pool(s)' gap."""
        import rmgpy.solver.polymer as solver_mod
        from rmgpy.polymer import Polymer
        from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

        Proxy = _spc("CCCC", "PS")                      # static parent pool proxy
        PMu0 = _spc("CO", "PS_mu0"); PMu1 = _spc("C=O", "PS_mu1"); PMu2 = _spc("C#N", "PS_mu2")
        daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                           end_groups=["[CH3]", "[H]"], cutoff=3,
                           Mn=5000.0, Mw=6000.0, initial_mass=0.001)
        DMu0 = self._moment_dummy("PS_d1_mu0")
        DMu1 = self._moment_dummy("PS_d1_mu1")
        DMu2 = self._moment_dummy("PS_d1_mu2")
        core = [Proxy, PMu0, PMu1, PMu2, daughter, DMu0, DMu1, DMu2]
        spc_map = {s: i for i, s in enumerate(core)}
        gas_mask = np.array([False] * 8, dtype=bool)

        rxn = Reaction(reactants=[Proxy], products=[daughter], **_KIN)
        rxn.polymer_flux_archetype = solver_mod.FLUX_SCISSION_FRAGMENT  # stamped 3

        parent = PolymerPoolConfig(
            label="PS", xs=3, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        # The daughter config is DERIVED (not hand-fed) -- the unit under test.
        derived = derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})
        assert [c.label for c in derived] == ["PS_d1"]

        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={}, V_poly=1.0,
            polymer_pools=[parent] + derived, mass_transfer=[],
            gas_species_mask=gas_mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"PS": (1.0, 5.0, 30.0), "PS_d1": (0.0, 0.0, 0.0)},
            termination=[],
        )
        rs.initialize_model(core, [rxn], [], [])

        # Resolved, NOT demoted: both the parent (src) and derived daughter (dst)
        # pools resolve, so the scission kernel applies instead of legacy mu1-only.
        assert rs.reaction_flux_archetype[0] == solver_mod.FLUX_SCISSION_FRAGMENT


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
               pools=(("A", (1, 2, 3)),), moments=None, k_unzip=0.0,
               monomer_poly_index=None):
    """Build + initialize a HybridPolymerSystem for the item-17 fixtures.

    ``pools`` is a tuple of (label, mu_indices) — the §5 config-state axis:
    adding ("G", ...) is the solver-level expression of an item-16 daughter
    config (spec §3(a) stage-2 labels).

    ``k_unzip`` (default 0.0, no-op for the existing callers) arms the lumped
    chain-end unzip channel on every pool; it drains the moment coordinates
    (dμ1/dt -= k_unzip·μ0, dμ2/dt -= k_unzip·(2μ1−μ0)) — i.e. it loads the
    moment-dummy core positions with large, known fluxes. That is exactly the
    char_rate dimensional bug. Arming it REQUIRES ``monomer_poly_index`` (the
    released monomer's core index, applied to every pool): the solver
    invariant refuses k_unzip > 0 with no emission target, because that shape
    drains condensed mass to nowhere."""
    moments = moments if moments is not None else {
        lbl: (1.0, 5.0, 30.0) for lbl, _ in pools}
    mask_arr = np.array(mask, dtype=bool)
    seed_idx = int(np.where(mask_arr)[0][0])
    pool_cfgs = [PolymerPoolConfig(label=lbl, xs=2,
                                   explicit_dp_to_species_index={},
                                   mu_indices=mu,
                                   monomer_poly_index=monomer_poly_index,
                                   k_unzip=k_unzip)
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


class TestCharRateMomentDummyExclusion:
    """Fix #2a (2026-06-27): the enlargement characteristic flux
    ``char_rate = sqrt(Σ core_species_rates²)`` (base.pyx:822) must be an L2
    norm over REAL species production rates only. Moment-dummy core positions
    (PS_mu0/_mu1/_mu2, ``is_moment_dummy``) carry moment-COORDINATE derivatives
    (e.g. μ2 under k_unzip), not molar species fluxes, and mixing them into the
    norm is a dimensional bug: the lumped channel scales char_rate ~linearly in
    k_unzip and buries family chemistry below tol_move_to_core, leaving the core
    empty. char_rate must exclude moment dummies; the real monomer-release flux
    stays (2a, NOT 2b — we do not delete real species chemistry to force
    promotion)."""

    def _b1_with_unzip(self, k_unzip):
        """The census B1 fixture (A→G gated edge + slow gas driver X→Y so
        char_rate>0) with the chain-end unzip armed on pool A. Unzip drains
        μ1/μ2 → the A_mu1/A_mu2 core positions carry large moment-coordinate
        fluxes — and releases the monomer into R17 (GAS real species since
        the 2026-07-03 monomer-gas fix, rate k_unzip·μ0 = 1.0·μ0), the
        solver-invariant-mandated emission target. R17's flux is REAL
        chemistry and belongs IN char_rate (2a, not 2b); the
        moment-coordinate fluxes do not."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["Y"], sp["R17"]]
        mask = [False, False, False, False, True, True, True]
        gated = Reaction(reactants=[sp["A"]], products=[sp["G"]], **_KIN)
        driver = Reaction(
            reactants=[sp["X"]], products=[sp["Y"]],
            kinetics=Arrhenius(A=(1.0e-3, "1/s"), n=0.0,
                               Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)
        rs = _gate17_rs(core, mask, [driver], edge_spcs=[sp["G"]],
                        rxns_edge=[gated], k_unzip=k_unzip,
                        monomer_poly_index=6)  # unzip releases into R17
        return sp, core, driver, gated, rs

    def test_moment_dummy_flux_does_not_deflate_enlargement_ratio(self, caplog):
        """RED→GREEN tracer: arm k_unzip so the moment-dummy μ2 derivative
        (≈ −9·k_unzip) dwarfs the real gas-driver flux. The census ungated_ratio
        for the gated product G is edge_rate/char_rate; with the bug, the moment
        dummies inflate char_rate and crush the ratio ~400×. After 2a, char_rate
        is the real-species-only L2 norm and the ratio matches the hand-computed
        value:
            ungated G rate = kf·μ1/V_poly = 10.0
            char_rate (real only) = sqrt(emission² + 2·driver²) with
              emission = k_unzip·μ0·V_poly / V_gas (monomer release into the
              GAS R17, 2026-07-03 monomer-gas fix — REAL chemistry, kept in
              the yardstick: 2a, not 2b)
              driver = 1e-3 / V_gas per X/Y, V_gas = R·800/1e5.
        The moment-coordinate fluxes (dμ1 = −1, dμ2 = −9) stay excluded: were
        they mixed back in, char_rate would shift well beyond the assertion's
        2e-2 tolerance."""
        import re
        sp, core, driver, gated, rs = self._b1_with_unzip(k_unzip=1.0)
        with caplog.at_level(logging.WARNING):
            _gate17_simulate(rs, core, [driver], [sp["G"]], [gated])
        lines = _census_lines(caplog)
        assert len(lines) == 1, lines
        v_gas = constants.R * 800.0 / 1.0e5
        emission = 1.0 * 1.0 / v_gas  # k_unzip·μ0·V_poly [mol/s] / V_gas
        char_real = np.sqrt(emission ** 2 + 2.0 * (1.0e-3 / v_gas) ** 2)
        expected_ratio = 10.0 / char_real
        ratio = float(re.search(r"ungated_ratio=([0-9.eE+-]+)",
                                lines[0]).group(1))
        assert ratio == pytest.approx(expected_ratio, rel=2e-2), (
            f"char_rate still mixes moment-coordinate flux: ungated_ratio="
            f"{ratio:.4g}, expected real-only {expected_ratio:.4g}")

    def test_include_mask_excludes_pool_mu_indices(self):
        """Populate guard (the seam that silently breaks as the core grows):
        after initialize_model the mask spans the core index space and is False
        at exactly the pool moment-coordinate positions, True for real species.
        Uses mu_indices (the copy-proof authority) -- NOT is_moment_dummy, which
        Species.copy() drops (species.py:784), so a flag-only mask would go
        live-inert on copied core dummies."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["Y"]]
        mask = [False, False, False, False, True, True]
        rs = _gate17_rs(core, mask, [], pools=(("A", (1, 2, 3)),))
        m = np.asarray(rs._char_rate_include_mask, dtype=bool)
        assert m.shape == (len(core),)
        # moment-coordinate positions excluded; real species (proxy A, gas X/Y)
        # included -- the monomer-release yardstick survives (2a, not 2b).
        assert list(m) == [True, False, False, False, True, True]

    def test_include_mask_honors_is_moment_dummy_flag(self):
        """The mask also drops any position flagged is_moment_dummy even when it
        is NOT a registered pool mu_index -- the union honors the contract's
        named identifier. (mu_indices is the robust primary; the flag is the
        belt-and-suspenders for a flagged dummy outside an active pool.)"""
        sp = _gate17_species()
        # R17 is an ordinary condensed species, not in any pool's mu_indices;
        # flag it a moment dummy and it must drop out of the yardstick.
        sp["R17"].is_moment_dummy = True
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["R17"]]
        mask = [False, False, False, False, True, False]
        rs = _gate17_rs(core, mask, [], pools=(("A", (1, 2, 3)),))
        m = np.asarray(rs._char_rate_include_mask, dtype=bool)
        # index 5 (R17) excluded via the flag; indices 1-3 via mu_indices.
        assert list(m) == [True, False, False, False, True, False]


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


class TestGasMaskInvariantToProxyTagContamination:
    """Regression (polystyrene reactive deck, 2026-06-26): a CORE gas solvent
    (N2) must keep the SAME phase verdict from get_gas_mask before AND after
    reaction generation stamps is_polymer_proxy onto it.

    Mechanism the deck actually hits (NOT label-duplication -- the live deck has
    zero duplicate core/edge labels and label_fallback stays enabled): family.py
    :1657 blanket-tags EVERY structure of a proxy-touching reaction with
    is_polymer_proxy, and model.py:486 propagates the tag to the Species. The
    solvent N2 reacts with the C25 PS proxy during initial-reaction generation,
    so the SAME core N2 object that was is_polymer_proxy=False at build time
    (gas_species_mask built GAS) becomes is_polymer_proxy=True by the time the
    solver rebuilds the prospective mask. If get_gas_mask registers
    initial_explicit species on the mutable is_polymer_proxy tag, N2 flips
    CONDENSED at simulate time, the prospective core prefix diverges from
    gas_species_mask, and RIDER R1 raises. The verdict for a core species must
    be invariant to this tag flip; is_moment_dummy (set once at creation, never
    propagated) is the stable polymer-phase marker.
    """

    def _phase_and_core(self):
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import PolymerPhase, PolymerPool
        proxy = _spc("CCCC", "PS")          # pool proxy / monomer
        mu0 = _spc("CO", "PS_mu0")
        mu1 = _spc("C=O", "PS_mu1")
        mu2 = _spc("C#N", "PS_mu2")
        for m in (mu0, mu1, mu2):
            m.is_moment_dummy = True
        n2 = _spc("N#N", "N2")              # gas solvent in initial_explicit
        pool = PolymerPool(label="PS", xs=3, monomer=proxy, explicit_map={},
                           mu_species=[mu0, mu1, mu2])
        phase = PolymerPhase(density=Quantity(1050.0, "kg/m^3"),
                             initial_moments={"PS": (1.0, 5.0, 30.0)},
                             initial_explicit={n2: 0.99}, pools=[pool])
        core = [proxy, mu0, mu1, mu2, n2]   # N2 at index 4
        return phase, core, n2

    def test_solvent_verdict_invariant_to_proxy_tag(self):
        """N2 is GAS at build (tag False) and MUST stay GAS after the proxy tag
        is stamped on it at simulate time. On the buggy code the tagged N2 is
        registered by get_gas_mask section A and flips CONDENSED."""
        phase, core, n2 = self._phase_and_core()
        clean = np.asarray(phase.get_gas_mask(core), dtype=bool)
        assert bool(clean[4]) is True  # N2 classified GAS at build

        # Simulate-time contamination (family.py:1657 -> model.py:486).
        n2.is_polymer_proxy = True
        contaminated = np.asarray(phase.get_gas_mask(core), dtype=bool)
        assert bool(contaminated[4]) is True  # N2 MUST stay GAS
        assert np.array_equal(clean, contaminated)

    def test_core_prefix_invariant_with_contaminated_edge_appended(self):
        """The R1 invariant directly: get_gas_mask(core)[N2] == get_gas_mask(
        chain(core, edge))[N2] even after the proxy tag is stamped on N2 and
        edge species (incl. a duplicate label) are appended. The combined-list
        core prefix must not diverge from the core-only mask for the solvent."""
        phase, core, n2 = self._phase_and_core()
        n2.is_polymer_proxy = True  # contaminated, as at simulate time
        edge_dup_a = _spc("CC", "FRAG")
        edge_dup_b = _spc("CCC", "FRAG")  # duplicate edge label (red herring)
        core_only = np.asarray(phase.get_gas_mask(core), dtype=bool)
        combined = np.asarray(
            phase.get_gas_mask(core + [edge_dup_a, edge_dup_b]), dtype=bool)
        assert bool(core_only[4]) is True
        assert np.array_equal(core_only, combined[:len(core)])


class TestCondensedEdgeDaughterBases:
    """Spec 2026-06-29 -- the pure-Python qualifying predicate that decides
    which edge Polymer daughters get prospectively condensed. Policy lives in
    PolymerPhase (testable with no Cython rebuild); the solver only applies it
    (Task 2). The moment-triplet + is_moment_dummy check is load-bearing;
    is_polymer_proxy is over-stamped and must never qualify a species alone."""

    @staticmethod
    def _phase(pool_labels=()):
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import PolymerPhase, PolymerPool
        pools = [PolymerPool(label=lbl, xs=2, monomer=_spc("CCCC", lbl),
                             explicit_map={}, mu_species=[])
                 for lbl in pool_labels]
        return PolymerPhase(density=Quantity(1000.0, "kg/m^3"),
                            initial_moments={}, initial_explicit={},
                            pools=pools, mass_transfer=[])

    @staticmethod
    def _daughter(label):
        from rmgpy.polymer import Polymer
        d = Polymer(label=label, monomer="[CH2][CH]c1ccccc1", Mn=100.0, Mw=200.0)
        d.is_polymer_proxy = True
        return d

    @staticmethod
    def _dummy(owner_label, k):
        # production convention: dummy label = {owner.label}_mu{k}
        s = _spc("CO", f"{owner_label}_mu{k}")
        s.is_moment_dummy = True
        return s

    def _triplet(self, owner_label):
        return [self._dummy(owner_label, k) for k in (0, 1, 2)]

    def test_qualifies_paren_less_production_daughter(self):
        # real PS shape: daughter PS_scission_tail_2 + PS_scission_tail_2_mu0/1/2
        phase = self._phase()
        d = self._daughter("PS_scission_tail_2")
        combined = [d] + self._triplet("PS_scission_tail_2")
        assert phase.get_condensed_edge_daughter_bases(combined) == {"PS_scission_tail_2"}

    def test_qualifies_parenthesized_daughter_as_base(self):
        # production dummy convention for a parenthesized daughter: D(2) +
        # D(2)_mu0/1/2 -> both normalize to base 'D' -> qualifies as 'D'
        phase = self._phase()
        d = self._daughter("D(2)")
        combined = [d] + self._triplet("D(2)")
        assert phase.get_condensed_edge_daughter_bases(combined) == {"D"}

    def test_excludes_static_proxy_by_base(self):
        # static_pool_labels = {"PS"}; candidate PS(2) + PS(2)_mu0/1/2 ->
        # base 'PS' is a static deck pool -> excluded
        phase = self._phase(pool_labels=("PS",))
        d = self._daughter("PS(2)")
        combined = [d] + self._triplet("PS(2)")
        assert phase.get_condensed_edge_daughter_bases(combined) == set()

    def test_rejects_incomplete_triplet(self):
        phase = self._phase()
        d = self._daughter("PS_scission_tail_2")
        combined = [d, self._dummy("PS_scission_tail_2", 0),
                    self._dummy("PS_scission_tail_2", 1)]   # missing _mu2
        assert phase.get_condensed_edge_daughter_bases(combined) == set()

    def test_rejects_dummy_missing_moment_flag(self):
        phase = self._phase()
        d = self._daughter("PS_scission_tail_2")
        m0 = _spc("CO", "PS_scission_tail_2_mu0")   # NO is_moment_dummy set
        combined = [d, m0, self._dummy("PS_scission_tail_2", 1),
                    self._dummy("PS_scission_tail_2", 2)]
        assert phase.get_condensed_edge_daughter_bases(combined) == set()

    def test_rejects_non_polymer_even_with_sticky_proxy_stamp(self):
        phase = self._phase()
        # (a) a plain Species sharing the daughter shape -- not a Polymer
        non_poly = _spc("C=Cc1ccccc1", "PS_scission_tail_2")
        combined_a = [non_poly] + self._triplet("PS_scission_tail_2")
        assert phase.get_condensed_edge_daughter_bases(combined_a) == set()
        # (b) the reachable over-stamp hazard: an ORDINARY Species carrying a
        # sticky is_polymer_proxy=True plus a complete moment triplet must NOT
        # qualify, because it is not a Polymer. A non-proxy Polymer is not a
        # meaningful fixture here: Polymer construction stamps the proxy flag and
        # Species.is_polymer_proxy latches true. Sticky ordinary species are the
        # reachable over-stamp hazard. This pins the safety property: raw
        # is_polymer_proxy alone never condenses an edge daughter -- the daughter
        # must be a Polymer AND have the full moment triplet.
        sticky = _spc("CCO", "OTH")
        sticky.is_polymer_proxy = True
        combined_b = [sticky] + self._triplet("OTH")
        assert phase.get_condensed_edge_daughter_bases(combined_b) == set()


def _condensed_edge_stub(bases):
    """Stand-in for the bound PolymerPhase.get_condensed_edge_daughter_bases on
    the solver application path (Task 2): returns a fixed base-label set,
    ignoring the live list. Isolates the .pyx APPLY logic from the predicate
    (predicate correctness is pinned in TestCondensedEdgeDaughterBases)."""
    def _classify(combined_species):
        return set(bases)
    return _classify


class TestEdgeDaughterCondensedMaskApply:
    """Spec 2026-06-29 -- the solver-side APPLY pass. A passed-in classifier
    (callable, never frozen) is re-run over the live combined list every
    initialize_model; qualifying EDGE slots flip CONDENSED in
    prospective_gas_mask, after the stage-2 override and before the R1
    core-prefix tripwire. Edge slots only -- core prefix parity (R1) holds."""

    @staticmethod
    def _build(core, mask, edge_spcs, classifier):
        cfg = PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None,
                                k_unzip=0.0)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[-1]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=np.asarray(mask, dtype=bool).copy(),
            constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_condensed_edge_daughter_classifier=classifier,
            # direct-test build: no blueprint phase -> permit the default edge
            # fill so R1-EDGE does not raise; the condensed-mask pass then flips
            # the qualifying edge slot CONDENSED on top of that GAS default.
            allow_default_prospective_edge=True)
        rs.initialize_model(list(core), [], list(edge_spcs), [])
        return rs

    def _core_and_mask(self):
        # core: A (melt proxy), A_mu0..2 (dummies), G (gas seed). mask: gas=True
        core = [_spc("CCCC", "A"), _spc("CO", "A_mu0"), _spc("C=O", "A_mu1"),
                _spc("C#N", "A_mu2"), _spc("[CH3]", "G")]
        mask = [False, False, False, False, True]
        return core, mask

    def test_qualifying_edge_daughter_flips_condensed_others_gas(self):
        core, mask = self._core_and_mask()
        d = _spc("C=Cc1ccccc1", "D(2)")     # qualifying edge daughter, base "D"
        other = _spc("CCO", "OTH")          # non-qualifying edge gas species
        rs = self._build(core, mask, [d, other], _condensed_edge_stub({"D"}))
        n_core = len(core)
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert bool(pm[n_core]) is False        # D(2) edge slot -> CONDENSED
        assert bool(pm[n_core + 1]) is True     # OTH edge slot -> still GAS
        # core prefix bit-identical to gas_species_mask (R1 holds, no core write)
        assert np.array_equal(pm[:n_core],
                              np.asarray(rs.gas_species_mask, dtype=bool))

    def test_classifier_is_live_not_frozen_across_edge_growth(self):
        """MANDATORY staleness regression. Build/init with edge_A (no daughter
        'D'), then re-run initialize_model with edge_B that introduces D(2). A
        constructor-FROZEN set computed at build time would miss D; the live
        callable, re-run over the larger combined list, flips it CONDENSED.
        This is the bug the authors already hit with the frozen prospective
        seed (polymer_input.py:283-291)."""
        core, mask = self._core_and_mask()
        d = _spc("C=Cc1ccccc1", "D(2)")
        rs = self._build(core, mask, [], _condensed_edge_stub({"D"}))  # edge_A: empty
        # edge_A had no D -> nothing condensed beyond the core
        pm_a = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert pm_a.shape[0] == len(core)               # no edge yet
        # edge_B introduces D(2); same solver object, larger edge
        rs.initialize_model(list(core), [], [d], [])
        pm_b = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert bool(pm_b[len(core)]) is False           # D(2) now CONDENSED

    def test_classifier_failure_is_loud_no_gas_default(self):
        """A raising classifier must propagate, never silently leave the
        daughter GAS (which would re-hide Gate B)."""
        core, mask = self._core_and_mask()
        d = _spc("C=Cc1ccccc1", "D(2)")

        def _boom(combined_species):
            raise RuntimeError("classifier exploded")

        with pytest.raises(RuntimeError, match="classifier exploded"):
            self._build(core, mask, [d], _boom)


class TestCondensedEdgeDaughterClassifierPlumb:
    """Spec 2026-06-29 -- the wiring path: to_solver_object must pass the bound
    PolymerPhase.get_condensed_edge_daughter_bases into the solver as the
    prospective_condensed_edge_daughter_classifier, so production builds run
    the real predicate over the live edge (not a stub)."""

    @staticmethod
    def _reactor():
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import HybridPolymerReactor, PolymerPhase
        a = _spc("CCCC", "A")
        phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"), initial_moments={},
                             initial_explicit={a: 1.0}, pools=[], mass_transfer=[])
        return a, phase, HybridPolymerReactor(
            temperature=(800.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles={a: 1.0}, polymerPhase=phase,
            terminationTime=(1.0, "s"))

    def test_to_solver_object_passes_bound_predicate(self):
        a, phase, reactor = self._reactor()
        solver = reactor.to_solver_object([a], [], [], [])
        clf = solver.prospective_condensed_edge_daughter_classifier
        assert clf is not None
        # bound method of THIS phase's get_condensed_edge_daughter_bases
        assert clf.__self__ is phase
        assert clf.__func__.__name__ == "get_condensed_edge_daughter_bases"


class TestMonomerProductGasEmission:
    """Incident 2026-07-03 (PS regeneration refusal): the deck-declared
    monomer_product (the unzip release target, e.g. styrene) is the
    mechanism's principal GAS volatile, but the solver's pool-override pass
    forced it CONDENSED -- so every reversible gas core reaction producing it
    carried an ~11.17-decade unpaired reference-state term and the (correct)
    tripwire refused the run. Ratified fix (design B-prime): the release
    target is validated GAS and the unzip/QSSA release deposits into the gas
    species amount basis; the condensed pool-moment drain is unchanged."""

    @staticmethod
    def _ps_incident_system():
        """Mirror of the real PS incident: pool PS (proxy + mu dummies,
        condensed) with monomer_poly_index -> styrene, PLUS a reversible
        gas-phase core reaction producing styrene against gas partners
        (1-phenylethyl <=> H + styrene). Styrene is deck-declared GAS."""
        sp = {
            "PS": _spc("CCCC", "PS"),
            "PS_mu0": _spc("CO", "PS_mu0"), "PS_mu1": _spc("C=O", "PS_mu1"),
            "PS_mu2": _spc("C#N", "PS_mu2"),
            "R": _spc("C[CH]c1ccccc1", "R"),     # 1-phenylethyl, gas
            "H": _spc("[H]", "H"),               # gas
            "STY": _spc("C=Cc1ccccc1", "STY"),   # styrene: monomer_product AND gas volatile
        }
        for s in sp.values():
            s.thermo = _trivial_nasa(_LIB_COMMENT)
        core = [sp["PS"], sp["PS_mu0"], sp["PS_mu1"], sp["PS_mu2"],
                sp["R"], sp["H"], sp["STY"]]
        mask = np.array([False] * 4 + [True, True, True], dtype=bool)
        pool = PolymerPoolConfig(
            label="PS", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=6,
            monomer_mw_g_mol=104.15, k_scission=0.0, k_unzip=0.1,
            tail_kinetics=None)
        return sp, core, mask, pool

    def test_incident_pin_gas_monomer_product_with_reversible_gas_chemistry(self):
        """THE incident pin (red-first). On pre-fix HEAD the pool override
        flips styrene CONDENSED and initialize_model REFUSES via the
        reference-state tripwire (U ~ 11 decades, one melt participant facing
        all-gas partners). Post-fix the build must initialize cleanly with
        styrene GAS and NOT pool-mapped (its gas chemistry runs on the gas
        concentration basis)."""
        sp, core, mask, pool = self._ps_incident_system()
        rxn = Reaction(reactants=[sp["R"]], products=[sp["H"], sp["STY"]],
                       **_REV_KIN)

        # LIVENESS PIN: were styrene condensed, its independently recomputed
        # unpaired reference-state magnitude would be far above the 3.0-decade
        # refuse bound -- i.e. this fixture genuinely reproduces the incident
        # shape, and pre-fix HEAD refuses it (verified red).
        mw_sty = sp["STY"].molecule[0].get_molecular_weight()
        u_if_condensed = (_sackur_tetrode_decades(mw_sty, 800.0)
                          + math.log10(1.0e5 / (constants.R * 800.0 * 1.0)))
        assert u_if_condensed > 3.0, (
            "FIXTURE BROKEN, not a valid red: recomputed U for a condensed "
            f"styrene ({u_if_condensed:.2f}) is not above the refuse bound")
        assert rxn.reversible and mask[4] and mask[5] and mask[6]

        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={sp["H"]: 0.0},
            V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"PS": (1.0, 5.0, 30.0)}, termination=[],
        )
        # Pre-fix: raises "THERMO REFERENCE-STATE TRIPWIRE ... unpaired
        # reference-state" here. Post-fix: completes.
        rs.initialize_model(core, [rxn], [], [])

        final_mask = np.asarray(rs.gas_species_mask, dtype=bool)
        assert bool(final_mask[6]) is True, (
            "monomer_product target was condensed by the pool override -- "
            "the incident defect")
        # Not pool-mapped: its concentration in gas reactions must be the
        # real gas concentration, never the pooled 1.0 site convention.
        assert int(np.asarray(rs.species_to_pool_indices)[6]) == -1
        # The unzip release lands ON the gas styrene (amount basis, mol/s).
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert dn_dt[6] > 0.0

    def test_unzip_emission_deposits_to_gas_and_conserves_mass(self):
        """Mass conservation across the unzip release with the GAS deposit:
        the condensed mu1 drain (repeat units) must equal the gas monomer
        deposit mol-for-mol -- dn_dt[mu1] == -dn_dt[monomer], both equal to
        k_unzip * Mu0 in magnitude. Guards against a regression to the
        May-era silent mass leak (drained repeat units going nowhere) AND
        against double-deposit (condensed + gas)."""
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")  # GAS released-monomer target (post-fix contract)
        core = [Inert, Mu0, Mu1, Mu2, M]
        mask = np.array([True, False, False, False, True], dtype=bool)

        k_unzip = 0.5
        V_poly = 2.0
        Mu0_moles = 1.0
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4,
            k_scission=0.0, k_unzip=k_unzip, tail_kinetics=None)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0},
            V_poly=V_poly, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (Mu0_moles, 5.0, 30.0)},
            termination=[],
        )
        rs.initialize_model(core, [], [], [])

        # The release target must STAY gas (red now: the override flips it).
        assert bool(np.asarray(rs.gas_species_mask)[4]) is True

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        # moments dict is moles-of-moments: mu0_conc = Mu0/V_poly, and the
        # emission r_events*V_poly == k_unzip*Mu0 [mol/s] lands on the gas M.
        expected = k_unzip * Mu0_moles
        assert np.isclose(dn_dt[4], expected, rtol=1e-12)
        assert np.isclose(dn_dt[2], -expected, rtol=1e-12)  # mu1 drain
        # repeat-unit moles conserved across the phase boundary:
        assert abs(dn_dt[2] + dn_dt[4]) <= 1e-15 * max(1.0, expected)

    def test_qssa_release_deposits_to_gas_and_conserves_mass(self):
        """The radical_qssa_unzip release writes the SAME small_src path as
        legacy unzip and must therefore also deposit to the GAS species
        amount basis, one gas monomer mole per drained mu1 repeat unit."""
        Inert = _spc("N#N", "N2")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        M = _spc("C", "M")  # GAS released-monomer target (post-fix contract)
        core = [Inert, Mu0, Mu1, Mu2, M]
        mask = np.array([True, False, False, False, True], dtype=bool)
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=4, k_unzip=0.0,
            radical_qssa_unzip=_qssa_channel(),
        )
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={Inert: 1.0},
            V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)},
            termination=[],
        )
        rs.initialize_model(core, [], [], [])

        assert bool(np.asarray(rs.gas_species_mask)[4]) is True

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert dn_dt[4] > 0.0                       # release flows to gas M
        assert np.isclose(dn_dt[2], -dn_dt[4], rtol=1e-12)  # mu1 drain pairs it
