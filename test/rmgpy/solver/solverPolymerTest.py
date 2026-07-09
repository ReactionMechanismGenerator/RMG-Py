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
    # r71 FIX 4 escape (documented direct-test posture): these fixtures
    # deliberately build unstamped proxy rows to exercise the legacy
    # mu1-only demotion; production keeps the constructor default (False,
    # hard-fail) -- pinned by TestRefusalAdjudicationSurvivesRebuild.
    kwargs = dict(allow_unstamped_proxy_rows=True)
    kwargs.update(rs_kwargs or {})
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[-1]: 0.0}, V_poly=1.0,
        polymer_pools=pools, mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments=moments, termination=[],
        **kwargs,
    )
    rs.initialize_model(core, rxns, [], [])
    return rs


def _gate_pool_config(monomer_mw_g_mol=28.0, chain_mass_defect_g_mol=0.0,
                      monomer_heavy_atoms=0):
    """PolymerPoolConfig for the spawn-gate fixtures.

    ``monomer_mw_g_mol`` is added by the mass-flux-spawn-gate change (spec
    2026-06-10 §4.1). The field-presence guard below is deliberate RED-FIRST
    scaffolding, NOT a compatibility shim: it lets the Task-1 integrated
    tripwire run on pre-change HEAD and die at the GROSS-ARRAY assertion
    (the born-dead mechanism, spec §3.1) instead of at fixture construction.
    Once the field lands, the guard always takes the kwargs branch.

    ``monomer_heavy_atoms`` is the r89 dual-axis heavy denominator (same
    field-presence guard posture: pre-r89 HEAD builds the config without it
    and the dual-axis tests die at their own assertions, not at fixture
    construction). Default 0 = heavy axis uncomputable.
    """
    kwargs = dict(label="A", xs=2, explicit_dp_to_species_index={},
                  mu_indices=(1, 2, 3), monomer_poly_index=None)
    if any(f.name == "monomer_mw_g_mol"
           for f in dataclasses.fields(PolymerPoolConfig)):
        kwargs["monomer_mw_g_mol"] = monomer_mw_g_mol
    if monomer_heavy_atoms and any(f.name == "monomer_heavy_atoms"
                                   for f in dataclasses.fields(PolymerPoolConfig)):
        kwargs["monomer_heavy_atoms"] = monomer_heavy_atoms
    if chain_mass_defect_g_mol:
        kwargs["chain_mass_defect_g_mol"] = chain_mass_defect_g_mol
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


def _one_pool_gate_rs(rxn, core, mask, moments, monomer_mw_g_mol=28.0,
                      chain_mass_defect_g_mol=0.0):
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[5]: 0.0}, V_poly=1.0,
        polymer_pools=[_gate_pool_config(monomer_mw_g_mol,
                                         chain_mass_defect_g_mol)],
        mass_transfer=[],
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
            rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
                allow_unstamped_proxy_rows=True,
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
            (0.0, 0.0, 0.0),          # fully degenerate (r86: census only, RHS live from raw y)
        ]
        for mu0, mu1, mu2 in pathological_states:
            y = rxn_system.y.copy()
            y[1], y[2], y[3] = mu0, mu1, mu2
            dn_dt = rxn_system.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt)), (
                f"Non-finite dn_dt for moment state (mu0,mu1,mu2)=({mu0},{mu1},{mu2}): {dn_dt}. "
                "This would trigger an unrecoverable model resurrection."
            )
        # Negatives beyond the atol-derived exhaustion floor are NO LONGER
        # silently max(0,.)-clamped: r81's negative-moment rule makes them a
        # HARD error (integrator corruption must surface, not integrate).
        # The pre-r81 pin ((-1,-2,-3) stays finite) is superseded.
        y = rxn_system.y.copy()
        y[1], y[2], y[3] = -1.0, -2.0, -3.0
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            rxn_system.residual(0.0, y, np.zeros_like(y))

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
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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

    def test_refused_zero_while_volatile_ejection_contributes(self):
        """Cross-solver agreement pin (schema 2.4 refused-row contract), RMG
        half, numbers not strings: a refused reaction contributes EXACTLY
        zero to the residual (elementwise: the two-reaction system equals
        the VE-only system), while a non-refused VOLATILE_EJECTION reaction
        on the same pool still contributes (regression guard). TA pins the
        same two numbers on its side (tests/test_polymer_simulator.py)."""
        import numpy as np
        from rmgpy.reaction import Reaction
        from rmgpy.kinetics import Arrhenius
        from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
        Proxy = _spc("CCC(C)CCCC(C)CCCC(C)C", "epdm")
        Mu0 = _spc("CO", "epdm_mu0"); Mu1 = _spc("C=O", "epdm_mu1"); Mu2 = _spc("C#N", "epdm_mu2")
        Macro = _spc("CCC(C)CCC[C](C)CCCC(C)C", "epdm_macroradical")
        H = _spc("[H]", "H"); H2 = _spc("[H][H]", "H2")
        Vol = _spc("C=C", "C2H4")
        core = [Proxy, Mu0, Mu1, Mu2, Macro, H, H2, Vol]
        mask = np.array([False, False, False, False, False, True, True, True], dtype=bool)
        refused = Reaction(reactants=[Proxy, H], products=[Macro, H2],
                           kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0,
                                              Ea=(0.0, "J/mol"), T0=(1.0, "K")),
                           reversible=False)
        refused.polymer_flux_archetype = 4      # UNRESOLVED
        refused.polymer_refused = True

        def _ve():
            ve = Reaction(reactants=[Proxy], products=[Proxy, Vol],
                          kinetics=Arrhenius(A=(0.5, "1/s"), n=0.0,
                                             Ea=(0.0, "J/mol"), T0=(1.0, "K")),
                          reversible=False)
            ve.polymer_flux_archetype = 6       # VOLATILE_EJECTION
            ve.polymer_eject_units = 0.5
            return ve

        def _system(rxns):
            pool = PolymerPoolConfig(label="epdm", xs=2,
                                     explicit_dp_to_species_index={},
                                     mu_indices=(1, 2, 3),
                                     monomer_poly_index=None,
                                     k_scission=0.0, k_unzip=0.0,
                                     tail_kinetics=None)
            rs = HybridPolymerSystem(
                T=800.0, P=1.0e5, initial_mole_fractions={H: 1.0},
                V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
                gas_species_mask=mask.copy(), constant_gas_volume=False,
                initial_polymer_moments={"epdm": (1.0, 50.0, 3000.0)},
                termination=[])
            rs.initialize_model(core, rxns, [], [])
            return rs

        rs_both = _system([refused, _ve()])
        dn_both = rs_both.residual(0.0, rs_both.y, np.zeros_like(rs_both.y))[0]
        rs_ve = _system([_ve()])
        dn_ve = rs_ve.residual(0.0, rs_ve.y, np.zeros_like(rs_ve.y))[0]
        assert np.any(dn_ve != 0.0)              # the VE row is live
        assert np.array_equal(dn_both, dn_ve)    # refused adds exactly zero

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

        UPDATED for the adjudicated direction-specific availability law
        (run-5 DASPK forensics): the reverse leg's site now comes from the
        DEBITED pool B's mu1 (rr = kb * mu1_B = 2.4), no longer from the
        forward reactant pool A's (kb * mu1_A = 3.0, the run-5 defect this
        test previously pinned).
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
        rr = 0.6 * mu_b[1]              # kb * C(proxyB)=1, then *= B site -> 2.4
        mu3_a = mu_a[0] * (mu_a[2] / mu_a[1]) ** 3   # 216.0
        mu3_b = mu_b[0] * (mu_b[2] / mu_b[1]) ** 3   # 31.25
        bA1, bA2 = mu_a[2] / mu_a[1], mu3_a / mu_a[1]    # 6.0, 43.2
        bB1, bB2 = mu_b[2] / mu_b[1], mu3_b / mu_b[1]    # 2.5, 7.8125

        assert np.isclose(dn_dt[1], -(rf - rr))              # == legacy net (b0=1)
        assert np.isclose(dn_dt[2], -rf * bA1 + rr * bB1)    # -54.0
        assert np.isclose(dn_dt[3], -rf * bA2 + rr * bB2)    # -413.25
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

        # rr = kb * C(proxyB=1), then site-scaled by the DEBITED pool (B) mu1
        # (direction-specific availability, adjudicated Part C -- this
        # resolves the previously deferred reverse-site question: the
        # reverse leg drains B, so its availability comes from B).
        ev = 0.6 * mu_b[1]                              # 2.4
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

    # ------------------------------------------------------------------
    # Direction-specific source availability (run-5 DASPK IDID=-7
    # forensics, adjudicated Part C): each per-direction exchange leg's
    # event flux must scale with the moments of the pool ACTUALLY BEING
    # DEBITED in that direction. The reverse leg debits the dst pool, so
    # its site factor comes from the dst pool's OWN moments -- previously
    # the forward reactant pool's (the "deferred reverse-site question"),
    # which let a near-empty spawned pool be drained at a rate set by the
    # healthy source pool: mu2 negative in ~1e-25 s, corrector divergence,
    # h -> ~4e-21, IDID=-7, resurrection failure.
    # ------------------------------------------------------------------
    def test_cross_pool_reverse_availability_scales_with_debited_pool(self):
        """
        Run-5 shape: near-zero dst pool B (the run-5 forensic
        polypropylene_mod_3 moments rescaled x1e4 to sit just ABOVE the r81
        exhaustion floor of 1e-14 mol, preserving the pool's shape) with a
        finite reverse driving force on a cross-pool VOLATILE_EJECTION row.
        The outbound B fluxes must scale with B's own moments
        (ev = kb * mu1_B ~ 1.6e-13), not with healthy pool A's
        (ev = kb * mu1_A = 3.0, which pre-fix drained B's mu2 at -2831 mol/s
        -- 1e12 times B's entire mu2 content per second). r86
        re-adjudication (Codex-adjudicated, spar r86 -- r81 read-projection
        reverted): BELOW the floor the SAME availability law keeps holding
        from the raw state -- the ORIGINAL in-band forensic values yield a
        tiny LIVE outbound flux scaled by B's own moments, not an
        exactly-zero projection (second phase below).
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        mu_a = (1.0, 5.0, 30.0)
        mu_b = (1.9e-14, 2.6e-13, 6.1e-12)
        rs = _two_pool_rs(rxn, core, mask, mu_a, mu_b)
        rs.kf[0] = 0.0        # forward off: isolate the reverse leg
        rs.kb[0] = 0.6

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        ev = 0.6 * mu_b[1]                          # kb * mu1 of the DEBITED pool
        mu3_b = mu_b[0] * (mu_b[2] / mu_b[1]) ** 3
        bB1 = mu_b[2] / mu_b[1]
        bB2 = mu3_b / mu_b[1]
        # Outbound legs scale with B's own content (exact law, atol=0 so the
        # tiny magnitudes are pinned relatively, not swallowed by a default
        # absolute tolerance).
        assert np.isclose(dn_dt[5], -ev * 1.0, rtol=1e-10, atol=0.0)
        assert np.isclose(dn_dt[6], -ev * bB1, rtol=1e-10, atol=0.0)
        assert np.isclose(dn_dt[7], -ev * bB2, rtol=1e-10, atol=0.0)
        # No absolute-scale drain of a ~1e-14-scale pool (pre-fix: -2831).
        assert abs(dn_dt[7]) < 1e-9
        # Mass conservation holds EXACTLY across the fixed leg with the same
        # event flux: chains conserved; backbone units gain exactly a*ev
        # (the re-absorbed volatile's units), as for any reverse VE.
        assert np.isclose(dn_dt[1] + dn_dt[5], 0.0, atol=1e-30)
        assert np.isclose(dn_dt[2] + dn_dt[6], +a * ev, rtol=1e-12, atol=0.0)
        # r86 re-adjudication (r81 read-projection reverted): at the
        # ORIGINAL run-5 forensic moments (1.9e-18/2.6e-17/6.1e-16 mol --
        # all inside the max(SMALL_EPS, 100*atol) = 1e-14 census band) the
        # exhaustion census announces the pool, but the reverse leg keeps
        # the SAME mu_B-scaled availability law from the raw state -- a
        # tiny LIVE outbound flux, never an exactly-zero projection.
        mu_b2 = (1.9e-18, 2.6e-17, 6.1e-16)
        rs2 = _two_pool_rs(rxn, core, mask, mu_a, mu_b2)
        rs2.kf[0] = 0.0
        rs2.kb[0] = 0.6
        dn_dt2 = rs2.residual(0.0, rs2.y, np.zeros_like(rs2.y))[0]
        ev2 = 0.6 * mu_b2[1]
        mu3_b2 = mu_b2[0] * (mu_b2[2] / mu_b2[1]) ** 3
        bB1_2 = mu_b2[2] / mu_b2[1]
        bB2_2 = mu3_b2 / mu_b2[1]
        assert np.isclose(dn_dt2[5], -ev2 * 1.0, rtol=1e-10, atol=0.0)
        assert np.isclose(dn_dt2[6], -ev2 * bB1_2, rtol=1e-10, atol=0.0)
        assert np.isclose(dn_dt2[7], -ev2 * bB2_2, rtol=1e-10, atol=0.0)

    def test_cross_pool_reverse_fix_migration_mass_conservation_exact(self):
        """
        MIGRATION with BOTH legs live and distinguishable pools: all three
        pool-moment sums conserve to machine precision under the
        direction-specific availability law (source debit + destination
        credit use the same per-direction event flux).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))
        rs.kb[0] = 0.6

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        assert dn_dt[1] + dn_dt[5] == pytest.approx(0.0, abs=1e-13)
        assert dn_dt[2] + dn_dt[6] == pytest.approx(0.0, abs=1e-12)
        assert dn_dt[3] + dn_dt[7] == pytest.approx(0.0, abs=1e-11)
        # Reverse leg actually live (fixture liveness).
        assert dn_dt[5] != 0.0

    def test_cross_pool_reverse_fix_volatile_balance_exact(self):
        """
        Cross-pool VE with a discrete gas volatile, BOTH legs live:
        A(poolA) <=> B(poolB) + G. Backbone-unit balance must close EXACTLY
        with the SAME event fluxes that move the pool moments:
        d(mu1_A + mu1_B)/dt + a * dn_G/dt == 0 and chain count conserves
        (d(mu0_A + mu0_B)/dt == 0).
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))
        rs.kb[0] = 0.6
        y = rs.y.copy()
        y[8] = 1.0   # G moles so C(G) > 0 and the reverse leg is live

        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]

        assert dn_dt[8] != 0.0   # fixture liveness: net gas flux flows
        assert dn_dt[1] + dn_dt[5] == pytest.approx(0.0, abs=1e-12)
        assert (dn_dt[2] + dn_dt[6]) + a * dn_dt[8] == pytest.approx(0.0, abs=1e-11)

    def test_cross_pool_reverse_fix_healthy_pools_regression_pin(self):
        """
        Numeric regression pin captured at db76881b9 (PRE-fix head): with
        equal availability moments (mu1_A == mu1_B) the directional law
        coincides with the old shared-site law, so the full dn_dt vector is
        byte-identical. Guards the fix against collateral factor changes.
        """
        a = 1.135
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        rs = _two_pool_rs(rxn, core, mask, (1.0e-3, 4.0e-3, 2.0e-2),
                          (2.0e-3, 4.0e-3, 1.0e-2))
        rs.kb[0] = 0.6

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        expected = np.array([
            0.0, -0.0056, -0.031276, -0.21453825999999995,
            0.0, +0.0056, +0.02492, +0.15075579999999994, 0.0,
        ])
        assert np.allclose(dn_dt[:9], expected, rtol=1e-12, atol=0.0)

    def test_cross_pool_reverse_flux_vanishes_continuously(self):
        """
        Continuity: as the debited pool's moments -> 0 (fixed distribution
        shape, amplitude s = 1e-6, 1e-9, 1e-12, 1e-18 -- spanning the live
        region above the 1e-14 mol census floor AND the census band inside
        it, r86 re-adjudication: the r81 read-projection is reverted, so
        the linear law continues THROUGH the floor with no step at any
        threshold, no cliff -- the run-9d H-early wake-up shape must see a
        continuous RHS across the floor crossing) the outbound mu2 drain
        decreases monotonically and LINEARLY with s. Pre-fix the drain was
        CONSTANT (-2831.2 mol/s at every s: availability taken from the
        healthy source pool).
        """
        a = 1.135
        mu_a = (1.0, 5.0, 30.0)
        drains = []
        for s in (1.0e-6, 1.0e-9, 1.0e-12, 1.0e-18):
            sp, core, mask = _two_pool_species()
            rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
            rxn.polymer_flux_archetype = 6
            rxn.polymer_eject_units = a
            rs = _two_pool_rs(rxn, core, mask, mu_a,
                              (1.9 * s, 26.0 * s, 610.0 * s))
            rs.kf[0] = 0.0
            rs.kb[0] = 0.6
            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            drains.append(-dn_dt[7])
        # Monotone decrease toward 0, all finite and positive (still live --
        # INCLUDING the in-band amplitude s = 1e-18: census only, no
        # projected-to-zero cliff at the floor).
        assert drains[0] > drains[1] > drains[2] > drains[3] > 0.0
        # Linear in the pool amplitude: successive ratios == amplitude ratio.
        assert drains[1] / drains[0] == pytest.approx(1.0e-3, rel=1e-9)
        assert drains[2] / drains[1] == pytest.approx(1.0e-3, rel=1e-9)
        # Through the floor: the SAME line continues into the census band.
        assert drains[3] / drains[2] == pytest.approx(1.0e-6, rel=1e-9)

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
    # Detailed balance under the direction-specific availability law
    # (round-64 adjudication, P2 pin on 53e92e1b0): the forward event flux
    # scales with the SOURCE pool's availability moment and the reverse
    # event flux with the DESTINATION pool's same-order moment ("reverse
    # availability from the debited pool"). Consequence, pinned here as
    # the intended post-53e92e1b0 equilibrium semantics: for a reversible
    # cross-pool VE row A(poolA) <=> B(poolB) + G with kb = kf/Keq
    # (generate_rate_coefficients, polymer.pyx:2205/:2214), the net event
    # flux vanishes at
    #     kf * (mu1_A/V_poly) = kb * C_G * (mu1_B/V_poly)
    #     =>  C_G* = Keq * (mu1_A / mu1_B)
    # -- the balance point depends on the pool moments BY DESIGN ("the old
    # pairing was the bug"): under the pre-fix shared-site pairing the
    # balance point was C_G = Keq regardless of which pool the reverse leg
    # debited.
    # ------------------------------------------------------------------
    def _detailed_balance_rs(self, mu_a, mu_b, a):
        """Reversible cross-pool VE deck A(poolA) <=> B(poolB) + G with the
        solver's OWN thermo-paired reverse rate (kb = kf/Keq, not a manual
        rs.kb poke). All species carry the Ne-like trivial NASA except G,
        whose enthalpy offset (a6 = 15663) puts Keq(800 K) at ~2 mol/m^3 so
        the balance gas concentration is physically scaled; the tests read
        rs.Keq back rather than assuming the value. constant_gas_volume +
        V_gas0 = 1 m^3 makes C_G == y[G] directly tunable (with the default
        constant-P gas volume a lone gas species is pinned at P/RT)."""
        sp, core, mask = _two_pool_species()
        g_poly = NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, 15663.0, 3.35],
                                Tmin=(200, "K"), Tmax=(6000, "K"))
        for lbl, s in sp.items():
            s.thermo = (NASA(polynomials=[g_poly], Tmin=(200, "K"),
                             Tmax=(6000, "K"))
                        if lbl == "G" else _trivial_nasa())
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       **_REV_KIN)
        rxn.polymer_flux_archetype = 6
        rxn.polymer_eject_units = a
        pool_a = PolymerPoolConfig(label="A", xs=2,
                                   explicit_dp_to_species_index={},
                                   mu_indices=(1, 2, 3),
                                   monomer_poly_index=None, k_scission=0.0,
                                   k_unzip=0.0, tail_kinetics=None)
        pool_b = PolymerPoolConfig(label="B", xs=2,
                                   explicit_dp_to_species_index={},
                                   mu_indices=(5, 6, 7),
                                   monomer_poly_index=None, k_scission=0.0,
                                   k_unzip=0.0, tail_kinetics=None)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0},
            V_poly=1.0, polymer_pools=[pool_a, pool_b], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=True,
            V_gas0=1.0, initial_polymer_moments={"A": mu_a, "B": mu_b},
            termination=[],
        )
        rs.initialize_model(core, [rxn], [], [])
        return rxn, rs

    def test_cross_pool_ve_detailed_balance_at_new_law_equilibrium(self):
        """
        Detailed-balance regression (round-64 P2 pin): BOTH pools populated
        with UNEQUAL moments (every component differs), reversible VE row,
        kb generated by the solver from Keq. At the new law's equilibrium
        gas concentration C_G* = (kf/kb) * (mu1_A/mu1_B) = Keq * mu1_A/mu1_B
        the residual's net exchange flux is zero on EVERY affected row
        (both pools' mu0/mu1/mu2 and the gas row) to rtol 1e-9, and
        perturbing the gas away from C_G* produces net flux of the correct
        sign toward equilibrium with the exact magnitude -delta*kf*mu1_A.

        A full moment-space fixed point needs, beyond event-flux balance,
        consistency of the exchanged chain-length statistics: the forward
        leg moves A's length-biased bundle (1, bA1, bA2) and the reverse
        leg B's (1, bB1, bB2) with the +a shift, so the mu1/mu2 rows cancel
        iff bA1 = bB1 + a and bA2 = bB2 + 2*a*bB1 + a^2 (mu3 closure
        mu3 = mu0*(mu2/mu1)^3). The moments below satisfy exactly that
        while staying componentwise unequal (mu1_A/mu1_B = 2), so C_G* =
        2*Keq -- visibly DIFFERENT from the naive C_G = Keq point, which is
        precisely the adjudicated semantics (reverse availability from the
        debited pool).
        """
        a = 1.135
        mu_b = (4.0e-4, 5.0e-3, 0.07)
        bB1 = mu_b[2] / mu_b[1]                    # 14.0
        bB2 = mu_b[0] * bB1 ** 3 / mu_b[1]         # 219.52
        bA1 = bB1 + a                              # 15.135
        bA2 = bB2 + 2.0 * a * bB1 + a * a          # 252.588225
        mu1_a = 1.0e-2
        mu_a = (mu1_a * bA2 / bA1 ** 3, mu1_a, bA1 * mu1_a)  # (~7.286e-4, 1e-2, 0.15135)
        rxn, rs = self._detailed_balance_rs(mu_a, mu_b, a)

        # Probe pin: the reverse rate constant is thermo-paired kb = kf/Keq
        # with Keq the solver's own get_equilibrium_constant(T) (Kc, SI).
        assert rs.kb[0] > 0.0
        assert rs.kb[0] == pytest.approx(rs.kf[0] / rs.Keq[0], rel=1e-12)
        assert rs.Keq[0] == pytest.approx(
            rxn.get_equilibrium_constant(800.0), rel=1e-12)

        c_g_star = (rs.kf[0] / rs.kb[0]) * (mu_a[1] / mu_b[1])
        assert c_g_star == pytest.approx(2.0 * rs.Keq[0], rel=1e-12)
        y = rs.y.copy()
        y[8] = c_g_star * 1.0  # V_gas0 = 1 m^3: moles == concentration

        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]

        ev = rs.kf[0] * mu_a[1]  # one-way event flux at balance (V_poly=1)
        assert ev > 0.0          # fixture liveness: both legs actually flow
        # Every affected row zero, pinned RELATIVE to its own one-way flux
        # magnitude (rtol 1e-9), not swallowed by an absolute tolerance.
        refs = ev * np.array([1.0, bA1, bA2, 1.0, bA1, bA2, 1.0])
        assert np.all(np.abs(dn_dt[[1, 2, 3, 5, 6, 7, 8]]) <= 1e-9 * refs)

        # Perturbation: correct sign toward equilibrium, exact magnitude.
        # Net event flux = kf*mu1_A - kb*C_G*mu1_B = -delta*kf*mu1_A at
        # C_G = (1+delta)*C_G*.
        for delta in (+0.05, -0.05):
            y2 = rs.y.copy()
            y2[8] = c_g_star * (1.0 + delta)
            d2 = rs.residual(0.0, y2, np.zeros_like(y2))[0]
            assert d2[8] == pytest.approx(-delta * ev, rel=1e-9)
            # Excess G (+delta): reverse dominates, G consumed, chains
            # restored to A; deficit (-delta): mirror image.
            assert np.sign(d2[8]) == -np.sign(delta)
            assert np.sign(d2[1]) == +np.sign(delta)   # A mu0
            assert np.sign(d2[5]) == -np.sign(delta)   # B mu0

    def test_cross_pool_ve_event_balance_statistics_residue(self):
        """
        Companion pin: with NON-bundle-matched pools (the raw round-64
        example moments) at the SAME new-law balance point C_G* =
        Keq*mu1_A/mu1_B, the EVENT fluxes still cancel -- chain count and
        gas rows are exactly zero -- while the mu1/mu2 rows carry the pure
        statistics-exchange residue (forward removes A-statistics chains,
        reverse returns B-statistics chains): a fixed point of the counting
        rows, not of the full moment state. Pinned values at ev = kf*mu1_A:
        mu1_A: +ev*(bB1 + a - bA1); mu2_A: +ev*(bB2 + 2a*bB1 + a^2 - bA2);
        mu1_B/mu2_B mirror with the -a shift.
        """
        a = 1.135
        mu_a = (1.0e-3, 1.0e-2, 0.12)
        mu_b = (4.0e-4, 5.0e-3, 0.07)
        rxn, rs = self._detailed_balance_rs(mu_a, mu_b, a)

        bA1 = mu_a[2] / mu_a[1]                    # 12.0
        bA2 = mu_a[0] * bA1 ** 3 / mu_a[1]         # 172.8
        bB1 = mu_b[2] / mu_b[1]                    # 14.0
        bB2 = mu_b[0] * bB1 ** 3 / mu_b[1]         # 219.52

        c_g_star = (rs.kf[0] / rs.kb[0]) * (mu_a[1] / mu_b[1])
        y = rs.y.copy()
        y[8] = c_g_star
        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]

        ev = rs.kf[0] * mu_a[1]
        # Event balance: chain-count and gas rows zero.
        assert abs(dn_dt[1]) <= 1e-9 * ev
        assert abs(dn_dt[5]) <= 1e-9 * ev
        assert abs(dn_dt[8]) <= 1e-9 * ev
        # Statistics-exchange residue, exact (probe: +0.0627, +1.5957645,
        # -0.0627, -1.4534355 at ev = 0.02).
        assert dn_dt[2] == pytest.approx(+ev * (bB1 + a - bA1), rel=1e-12)
        assert dn_dt[3] == pytest.approx(
            +ev * (bB2 + 2.0 * a * bB1 + a * a - bA2), rel=1e-12)
        assert dn_dt[6] == pytest.approx(-ev * (bB1 + a - bA1), rel=1e-12)
        assert dn_dt[7] == pytest.approx(
            -ev * (bB2 - (bA2 - 2.0 * a * bA1 + a * a)), rel=1e-12)
        # mu1 residue is a pure A<->B transfer (conserved); mu2 is not.
        assert dn_dt[2] + dn_dt[6] == pytest.approx(0.0, abs=1e-12)

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
        # Overshoot mu2 negative with one huge explicit step, then
        # re-evaluate. Pre-r81 this pinned "finite + cone warning"; under
        # the r81 negative-moment rule a moment this far beyond -floor
        # (-1e10 mol vs -1e-14) is integrator corruption and HARD-errors
        # instead of being silently max(0,.)-clamped.
        y2 = y + 1e3 * dn_dt
        assert y2[3] < 0.0, "test setup: expected a mu2 overshoot"
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            rs.residual(0.0, y2, np.zeros_like(y2))
        # The cone diagnostic itself (debug_check_realizability) survives
        # for NONNEGATIVE out-of-cone states: mu0*mu2 < mu1^2 with all
        # moments >= 0 still warns once and stays finite.
        y3 = y.copy()
        y3[1], y3[2], y3[3] = 1.0, 5.0, 1.0   # 1*1 < 25: outside the cone
        with caplog.at_level(_logging.WARNING):
            dn_dt3 = rs.residual(0.0, y3, np.zeros_like(y3))[0]
        assert np.all(np.isfinite(dn_dt3))
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
        e_n, mw, defect = pool_stats["A"]
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
        e_n, mw, defect = pool_stats["A"]
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
        e_n, mw, defect = pool_stats["A"]
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
        e_n, mw, defect = pool_stats["A"]
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
        # mu0 exhausted (0 <= SMALL_EPS), mu1/mu2 small but ABOVE the r81
        # exhaustion floor (1e-14 mol): the pool is NOT fully exhausted, so
        # the row still carries flux and the deferral tested here is the
        # spawn-gate E[n] clamp itself (a fully exhausted pool would defer
        # even earlier -- zero RHS -> zero gross -- by r81 conditioning).
        rs = _one_pool_gate_rs(rxn, core, mask, (0.0, 1e-10, 1e-5))

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

    def test_refused_row_with_unpaired_reference_state_passes_initialize(self):
        """PP v1 gas-association refusal, pin pair (adjudicated round 63,
        design constraint 3): a ``polymer_refused`` row (stamp-but-keep, its
        WHOLE flux suppressed via ``reaction_refused`` in the residual) with a
        genuinely unpaired reference state must PASS ``initialize_model`` --
        a refused row's reference-state pairing is moot. The skip is keyed on
        the refused stamp ONLY: the SAME row unrefused must still trip (second
        half of the pair, below). RED at b917becd7: the tripwire census loop
        reads every reversible core row regardless of refused status, so the
        refused build raises the cliff-sign ValueError."""
        sp, core, mask = _refstate_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]],
                       **_REV_KIN)
        rxn.polymer_refused = True
        rxn.polymer_refused_accumulating = False   # conduit-deferred

        # LIVENESS PIN -- BEFORE the red assertion: the row genuinely carries
        # chain-scale unpaired U above the refuse bound (independently
        # recomputed, same recipe as the refusal test). A failure HERE means
        # the fixture is dead, not a valid red.
        assert rxn.reversible
        mw_a = sp["A"].molecule[0].get_molecular_weight()
        u_expected = (_sackur_tetrode_decades(mw_a, 800.0)
                      + math.log10(1.0e5 / (constants.R * 800.0 * 1.0)))
        assert u_expected > 3.0, (
            "FIXTURE BROKEN, not a valid red: independently recomputed U "
            f"({u_expected:.2f}) is not above the refuse bound"
        )

        # THE red assertion: the refused build must COMPLETE (no tripwire
        # ValueError) -- and the refusal is not a free pass: the row's flux
        # is suppressed and censused as conduit-deferred.
        rs = _refstate_rs(core, [rxn], mask,
                          [_gate_pool_config()], {"A": (1.0, 5.0, 30.0)})
        assert any(c["reason"] == "conduit-deferred"
                   for c in rs.refused_reaction_census)

        # Pin-pair second half: the SAME row UNREFUSED still trips -- the
        # skip is via the refused stamp only, no other relaxation.
        sp2, core2, mask2 = _refstate_pool_species()
        rxn2 = Reaction(reactants=[sp2["A"]], products=[sp2["G1"], sp2["G2"]],
                        **_REV_KIN)
        assert not getattr(rxn2, "polymer_refused", False)
        with pytest.raises(ValueError, match="unpaired reference-state"):
            _refstate_rs(core2, [rxn2], mask2,
                         [_gate_pool_config()], {"A": (1.0, 5.0, 30.0)})

    @staticmethod
    def _split_isomerization_fixture():
        """r87 shape-B fixture (FR1 run-3 ``(117) <=> (111)``): a real
        Polymer pool A in core (the restamp's pool registry is collected
        from the SPECIES lists, so the monomer scale must be discoverable
        there) plus two chain-scale proxy-derived C21H42Br2 adduct isomers:
        D_melt melt-classified (tagged, 454.37 g/mol, unvetoed) and D_gas
        veto-suppressed (tagged AND gas-vetoed, run-3's (117) posture). Both
        are genuine multi-repeat chain-scale defects (454.37 g/mol / 23 heavy)
        clearing the r95 absolute floor on both r85 axes -- and isomers, so
        they pair off exactly in the tripwire U."""
        from rmgpy.polymer import Polymer, set_polymer_gas_veto
        pool = Polymer(label="A", monomer="[CH2][CH]C",
                       Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        sp = {
            "A": pool,
            "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"),
            "A_mu2": _spc("C#N", "A_mu2"),
            "D_gas": _spc("BrCC(C)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", "D_gas"),
            "D_melt": _spc("CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", "D_melt"),
        }
        sp["D_melt"].is_polymer_proxy = True
        sp["D_gas"].is_polymer_proxy = True
        set_polymer_gas_veto(sp["D_gas"])
        for s in sp.values():
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
                sp["D_gas"], sp["D_melt"]]
        mask = np.array([False] * 4 + [True, True], dtype=bool)
        return sp, core, mask

    def test_reference_state_split_isomerization_refused_at_restamp(self, caplog):
        """r87 (FR1 run-3 shape B, RED-first): the ``(117) <=> (111)``
        reference-state-SPLIT isomerization -- one side melt-classified,
        the other veto-suppressed, U ~ 13 decades unpaired (run-3 RMG.out
        line 5140) -- arrives UNSTAMPED at initialize_model and must be
        refused (conduit-deferred) by the r71 rebuild restamp BEFORE the
        thermo reference-state tripwire scans core rows, so the build
        COMPLETES with the row quarantined instead of dying on the
        cliff-sign ValueError. RED at 9b22bba8d: no shape-B classifier
        exists, the live row reaches the tripwire and the build raises.

        Census pin (r87 'gas-veto census unchanged'): the veto-suppressed
        chain-scale adduct still surfaces through the THERMO
        REFERENCE-STATE GAS VETO warning -- refusal must not eat it."""
        import logging as _logging
        sp, core, mask = self._split_isomerization_fixture()
        rxn = Reaction(reactants=[sp["D_gas"]], products=[sp["D_melt"]],
                       **_REV_KIN)

        # LIVENESS PIN -- BEFORE the red assertion: the melt-classified side
        # genuinely carries chain-scale unpaired U above the 3.0-decade
        # refuse bound (independently recomputed, same recipe as the
        # refusal test). A failure HERE means the fixture is dead, not a
        # valid red.
        assert rxn.reversible
        mw_d = sp["D_melt"].molecule[0].get_molecular_weight()
        u_expected = (_sackur_tetrode_decades(mw_d, 800.0)
                      + math.log10(1.0e5 / (constants.R * 800.0 * 1.0)))
        assert u_expected > 3.0, (
            "FIXTURE BROKEN, not a valid red: independently recomputed U "
            f"({u_expected:.2f}) is not above the refuse bound"
        )
        assert not getattr(rxn, "polymer_refused", False)

        # THE red assertion: pre-change the unstamped split row reaches the
        # tripwire and initialize_model raises; post-change the restamp
        # refuses it first and the build completes.
        with caplog.at_level(_logging.WARNING):
            rs = _refstate_rs(core, [rxn], mask,
                              [_gate_pool_config(monomer_mw_g_mol=42.08,
                                                 monomer_heavy_atoms=3)],
                              {"A": (1.0, 5.0, 30.0)})
        assert rxn.polymer_refused is True
        assert rxn.polymer_refused_accumulating is False  # conduit-deferred
        assert any(c["reason"] == "conduit-deferred"
                   for c in rs.refused_reaction_census)
        # Gas-veto census unchanged: the vetoed chain-scale adduct is still
        # announced (run-3 RMG.out line 5113 posture).
        assert any("THERMO REFERENCE-STATE GAS VETO" in r.getMessage()
                   for r in caplog.records)

    def test_no_split_isomerization_stays_live_through_tripwire(self):
        """r87 negative control (pin): the SAME adduct isomerization with
        BOTH sides melt-classified (no veto anywhere) does NOT split the
        reference-state classification -- U pairs off exactly (isomers,
        same MW) -- so the row stays LIVE and the tripwire passes it."""
        from rmgpy.polymer import Polymer
        pool = Polymer(label="A", monomer="[CH2][CH]C",
                       Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        sp = {
            "A": pool,
            "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"),
            "A_mu2": _spc("C#N", "A_mu2"),
            "D1": _spc("BrCC(C)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", "D1"),
            "D2": _spc("CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", "D2"),
        }
        sp["D1"].is_polymer_proxy = True
        sp["D2"].is_polymer_proxy = True
        for s in sp.values():
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
                sp["D1"], sp["D2"]]
        mask = np.array([False] * 4 + [True, True], dtype=bool)
        rxn = Reaction(reactants=[sp["D1"]], products=[sp["D2"]], **_REV_KIN)
        rs = _refstate_rs(core, [rxn], mask,
                          [_gate_pool_config(monomer_mw_g_mol=42.08,
                                                 monomer_heavy_atoms=3)],
                          {"A": (1.0, 5.0, 30.0)})
        assert not getattr(rxn, "polymer_refused", False)
        assert not any(c["reason"] == "conduit-deferred"
                       for c in rs.refused_reaction_census)

    @staticmethod
    def _two_pool_shape_d_fixture():
        """r93 shape-D fixture (FR1 run-5 ``FR1 + FR1_sidegrp <=> (5) + FR1``,
        U = 13.05, RMG.out): a real PP-scale Polymer pool A in core (the
        restamp's pool registry is collected from the SPECIES lists, so the
        general branch sees a real Polymer participant) plus a second
        chain-scale melt participant B (mask-False, the FR1_sidegrp stand-in)
        and a chain-scale proxy-derived GAS-VETOED discrete D5 (the (5)
        stand-in, C21H42Br2 454.37 g/mol / 23 heavy, a genuine multi-repeat
        chain-scale defect clearing the r95 absolute floor). D5 is
        gas-vetoed, so it is EXCLUDED from the melt sum -- the two-melt-vs-one
        imbalance (B unpaired) is what the tripwire scores, exactly run-5."""
        from rmgpy.polymer import Polymer, set_polymer_gas_veto
        pool = Polymer(label="A", monomer="[CH2][CH]C",
                       Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        sp = {
            "A": pool,
            "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"),
            "A_mu2": _spc("C#N", "A_mu2"),
            # B: chain-scale melt participant (mask-False), the second pool
            # of the two-pool abstraction; distinct MW so it is UNPAIRED.
            "B": _spc("CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)C", "B"),
            "D5": _spc("CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", "D5"),
        }
        sp["D5"].is_polymer_proxy = True
        set_polymer_gas_veto(sp["D5"])
        for s in sp.values():
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        # The refused row still has its Keq computed by generate_rate_coefficients
        # (refusal only zeroes flux), so the pool proxy needs usable thermo.
        pool.get_proxy_species().thermo = _trivial_nasa(_GAV_COMMENT)
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
                sp["B"], sp["D5"]]
        # A, A_mu*, B are condensed/melt (mask False); D5 is gas (mask True).
        mask = np.array([False] * 5 + [True], dtype=bool)
        return sp, core, mask

    def test_two_pool_shape_d_row_refused_at_restamp(self, caplog):
        """r93 RED matrix #1 (FR1 run-5 shape D, RED-first, LIVE path): the
        two-pool abstraction ``A + B <=> D5 + A`` (D5 a vetoed chain-scale
        proxy-derived discrete) arrives UNSTAMPED at initialize_model and must
        be refused (conduit-deferred) by the r71 rebuild restamp BEFORE the
        thermo reference-state tripwire scans core rows, so the build COMPLETES
        with the row quarantined instead of dying on the cliff-sign ValueError.

        RED at f271af2ce: r74 same-proxy carves the row out (net gas mass
        change >= 0.5 monomer) and r87 shape A skips it (equal 1v1 pool count),
        so the live row reaches the tripwire and initialize_model raises with
        U = Sackur-Tetrode(B) unpaired.

        Census pin (run-5 'gas-veto census unchanged'): the vetoed chain-scale
        discrete still surfaces through the THERMO REFERENCE-STATE GAS VETO
        warning -- refusal must not eat it."""
        import logging as _logging
        sp, core, mask = self._two_pool_shape_d_fixture()
        rxn = Reaction(reactants=[sp["A"], sp["B"]],
                       products=[sp["D5"], sp["A"]], **_REV_KIN)

        # LIVENESS PIN -- BEFORE the red assertion: the unpaired melt
        # participant B genuinely carries chain-scale unpaired U above the
        # 3.0-decade refuse bound (independently recomputed). A failure HERE
        # means the fixture is dead, not a valid red.
        assert rxn.reversible
        mw_b = sp["B"].molecule[0].get_molecular_weight()
        u_expected = (_sackur_tetrode_decades(mw_b, 800.0)
                      + math.log10(1.0e5 / (constants.R * 800.0 * 1.0)))
        assert u_expected > 3.0, (
            "FIXTURE BROKEN, not a valid red: independently recomputed U "
            f"({u_expected:.2f}) is not above the refuse bound"
        )
        assert not getattr(rxn, "polymer_refused", False)

        # THE red assertion: pre-change the unstamped shape-D row reaches the
        # tripwire and initialize_model raises; post-change the restamp
        # refuses it first and the build completes.
        with caplog.at_level(_logging.WARNING):
            rs = _refstate_rs(core, [rxn], mask,
                              [_gate_pool_config(monomer_mw_g_mol=42.08,
                                                 monomer_heavy_atoms=3)],
                              {"A": (1.0, 5.0, 30.0)})
        assert rxn.polymer_refused is True
        assert rxn.polymer_refused_accumulating is False  # conduit-deferred
        assert any(c["reason"] == "conduit-deferred"
                   for c in rs.refused_reaction_census)
        # Gas-veto census unchanged: the vetoed chain-scale discrete is still
        # announced (run-5 RMG.out gas-veto census posture).
        assert any("THERMO REFERENCE-STATE GAS VETO" in r.getMessage()
                   for r in caplog.records)

    def test_two_pool_shape_d_row_unrefused_still_trips(self):
        """r93 RED matrix #1 pin-pair second half: the SAME shape-D row, if the
        restamp is prevented from seeing a Polymer participant (both melt
        participants plain gas-referenced Species, no pool object in the row),
        still trips the tripwire -- the build's refusal is via the general
        branch's Polymer-participant conjunct ONLY, no other relaxation."""
        sp, core, mask = self._two_pool_shape_d_fixture()
        # Swap the Polymer pool A OUT of the row for a plain melt Species with
        # the same reference state, so conjunct (ii) (a Polymer participant)
        # fails and the general branch cannot fire -- the row stays live and
        # the unpaired melt term still trips the tripwire.
        a_plain = _spc("CC(C)CC(C)CC(C)C", "A_plain")
        a_plain.thermo = _trivial_nasa(_GAV_COMMENT)
        core2 = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
                 sp["B"], sp["D5"], a_plain]
        mask2 = np.array([False] * 5 + [True, False], dtype=bool)
        rxn = Reaction(reactants=[a_plain, sp["B"]],
                       products=[sp["D5"], a_plain], **_REV_KIN)
        assert not getattr(rxn, "polymer_refused", False)
        with pytest.raises(ValueError, match="unpaired reference-state"):
            _refstate_rs(core2, [rxn], mask2,
                         [_gate_pool_config(monomer_mw_g_mol=42.08,
                                            monomer_heavy_atoms=3)],
                         {"A": (1.0, 5.0, 30.0)})

    def test_mixed_provenance_chain_counterparty_warns(self, caplog):
        """Spec §8.4: one melt-class species takes library thermo while its
        chain-scale counterparty takes GAV -> the mixed-provenance warning
        must fire. The counterparty is the §2 decoupling shape: a
        gas-classified but proxy-TAGGED same-mass radical (n-C22H45 309.60 vs
        n-C22H46 310.61 g/mol), so the pair is mass-paired
        (U = 1.5*log10(310.61/309.60) = 0.0021 -- no census, no refusal) and
        ONLY the sensor can speak. The small gas H also takes library thermo
        and must NOT matter (outside the counterparty window). r89/r95 lockstep:
        both A and B are GENUINE multi-repeat chains (C22, 22 heavy) clearing
        the r95 absolute chain-scale floor (ABS_CHAIN_SCALE_MW/HEAVY), so B
        clears the dual-axis polymer-sized melt gate -- a below-floor tagged
        radical would be conservative-gas and A would be UNPAIRED (refusal),
        which is a different test."""
        import logging
        sp = {
            "A": _spc("CCCCCCCCCCCCCCCCCCCCCC", "A"),
            "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"), "A_mu2": _spc("C#N", "A_mu2"),
            "B": _spc("[CH2]CCCCCCCCCCCCCCCCCCCCC", "B"),
            "H": _spc("[H]", "H"),
        }
        sp["A"].thermo = _trivial_nasa(_LIB_COMMENT)   # melt proxy: library
        sp["B"].thermo = _trivial_nasa(_GAV_COMMENT)   # chain counterparty: GAV
        sp["H"].thermo = _trivial_nasa(_LIB_COMMENT)   # small gas: library
        for k in ("A_mu0", "A_mu1", "A_mu2"):
            sp[k].thermo = _trivial_nasa(_LIB_COMMENT)
        # gas-CLASSIFIED, proxy-TAGGED, dual-axis polymer-sized against the
        # CH2 pool (57.11 / 14.03 = 4.07 mass-units, 4 / 1 = 4.0
        # heavy-units, both >= 2.5): physically-melt via the tag branch of
        # the r89-amended spec-§5.1 class (tag AND NOT veto AND dual-axis)
        sp["B"].is_polymer_proxy = True
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["B"], sp["H"]]
        mask = np.array([False] * 4 + [True, True], dtype=bool)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["H"]], **_REV_KIN)

        # LIVENESS PINS: the pair is mass-paired (inside the counterparty
        # window) and the build must SUCCEED on any head -- only the warning
        # distinguishes pre/post-change.
        mw_a = sp["A"].molecule[0].get_molecular_weight() * 1000.0
        mw_b = sp["B"].molecule[0].get_molecular_weight() * 1000.0
        assert abs(mw_a - mw_b) <= 14.03 + 10.0
        # amended-class pin (spec §5.1 C3, r89 dual-axis + r95 absolute floor):
        # B must clear the polymer-sized threshold on BOTH axes, or it would
        # not be physically-melt at all and the sensor would have no
        # melt-vs-counterparty pair to warn on
        from rmgpy.polymer import ABS_CHAIN_SCALE_MW, ABS_CHAIN_SCALE_HEAVY
        assert mw_b >= ABS_CHAIN_SCALE_MW
        heavy_b = (sp["B"].molecule[0].get_num_atoms()
                   - sp["B"].molecule[0].get_num_atoms('H'))
        assert heavy_b >= ABS_CHAIN_SCALE_HEAVY
        with caplog.at_level(logging.WARNING):
            _refstate_rs(core, [rxn], mask,
                         [_gate_pool_config(monomer_mw_g_mol=14.03,
                                            monomer_heavy_atoms=1)],
                         {"A": (1.0, 5.0, 30.0)})

        # THE red assertion: pre-change HEAD emits no provenance warning.
        assert any("THERMO REFERENCE-STATE PROVENANCE" in r.getMessage()
                   for r in caplog.records), (
            "mixed library-vs-GAV provenance among chain-scale counterparties "
            "went unwarned (the spec-§5.3 decoupling-fingerprint sensor is "
            "absent)"
        )

    def test_melt_sum_leak_guard_raises_classification_error(self):
        """Spec §5.1 C3 amendment (r89 dual-axis form): the cannot-happen
        leak guard inside the melt-sum accumulation. Under the amended class
        a tagged below-threshold species fails the dual-axis polymer-sized
        conjunct and is EXCLUDED by the gate -- expected and silent (the
        family.py:1657 over-tagging fingerprint, H2 on every proxy-touching
        reaction; the PP run-9 DP-2 hexadiene); the raise is only for such a
        species REACHING the melt sum. Because the gate and the guard share
        ONE dual-axis predicate, the violation is structurally unreachable
        through public paths -- so the raise is pinned by calling the helper
        directly with a hand-built violating member (documented bypass, per
        the C3 amendment)."""
        from rmgpy.solver.polymer import _assert_chain_scale_melt_member
        pp_axes = [(42.08, 3)]   # propene pool: monomer MW [g/mol], heavy atoms
        units = 2.5
        # Valid members never raise: a condensed-branch member (gate-exempt,
        # pool-configured by input -- any size, even below threshold) and a
        # genuine chain-scale tag-branch member (C22 backbone radical: 309.6
        # g/mol / 22 heavy -- clears the r95 absolute floor on both axes).
        _assert_chain_scale_melt_member("M1", 0.016043, 1, False, pp_axes, units)
        _assert_chain_scale_melt_member("Erad", 0.309602, 22, True, pp_axes, units)
        # The hand-built violations: gas-classified (tag-branch) members
        # below the dual-axis threshold inside the melt sum. The message must
        # steer the operator to CLASSIFICATION, never to reference states.
        with pytest.raises(ValueError, match="classification leak, NOT a thermo problem"):
            _assert_chain_scale_melt_member("H2", 0.002016, 0, True, pp_axes, units)
        with pytest.raises(ValueError, match="non-chain species in the melt sum"):
            _assert_chain_scale_melt_member("H2", 0.002016, 0, True, pp_axes, units)
        # r89: the run-9 DP-2 shape (hexadiene, 82.15 g/mol = 1.95
        # mass-units) violates too -- under the pre-r89 window form
        # (>= monomer + 10 slack) it slipped THROUGH this guard.
        with pytest.raises(ValueError, match="classification leak"):
            _assert_chain_scale_melt_member(
                "C=CCCC=C", 0.08215, 6, True, pp_axes, units)
        # And an axis-undecidable member (heavy unknown: -1) may never reach
        # the sum either -- undecidable is conservative-gas, not melt.
        with pytest.raises(ValueError, match="classification leak"):
            _assert_chain_scale_melt_member(
                "NOAXIS", 0.211407, -1, True, pp_axes, units)

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

    @staticmethod
    def _volatile_vs_pool_build(vol_smiles, pool_mw, pool_heavy, set_veto,
                                set_tag=True):
        """Shared r89 fixture: one-pool core (melt proxy A + mu dummies,
        indices 0-3), a condensed melt chain TAIL (index 4, gas_mask=False)
        and the discrete gas species VOL (index 5, gas_mask=True) in the
        unpaired shape ``A <=> VOL + TAIL`` (dn_melt=+1 if VOL classifies
        melt -> one uncancelled Sackur-Tetrode term -> U ~ 11 decades ->
        REFUSE; VOL classified gas -> dn_melt=0 -> benign, builds)."""
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO_KEY
        a = _spc("CCCC", "A")
        a_mu0, a_mu1, a_mu2 = (_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"),
                               _spc("C#N", "A_mu2"))
        tail = _spc("CCCCCCCC", "TAIL")
        vol = _spc(vol_smiles, "VOL")
        if set_tag:
            vol.is_polymer_proxy = True
        if set_veto:
            vol.props[VETO_KEY] = True
        for s in (a, a_mu0, a_mu1, a_mu2, tail, vol):
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [a, a_mu0, a_mu1, a_mu2, tail, vol]
        mask = np.array([False, False, False, False, False, True], dtype=bool)
        rxn = Reaction(reactants=[a], products=[vol, tail], **_REV_KIN)
        return _refstate_rs(core, [rxn], mask,
                            [_gate_pool_config(monomer_mw_g_mol=pool_mw,
                                               monomer_heavy_atoms=pool_heavy)],
                            {"A": (1.0, 5.0, 30.0)})

    @pytest.mark.parametrize("vol_smiles,pool_mw,pool_heavy,mass_units,heavy_units", [
        # r89 RED 1 sibling shape at the tripwire seam: PP run-9's
        # 1,5-hexadiene (C6H10, 82.15 g/mol) vs the propene pool -- above
        # the pre-r89 window (42.08 + 10 = 52.1) but a genuine DP-2 gas
        # volatile at 1.95 mass-units / 2.0 heavy-units < 2.5.
        ("C=CCCC=C", 42.08, 3, 82.15 / 42.08, 6 / 3.0),
        # r89 RED 2 (PS sibling): alpha-methylstyrene (118.18 g/mol) vs the
        # styrene pool -- above the pre-r89 window (104.15 + 10 = 114.2) but
        # 1.13 mass-units / 1.125 heavy-units < 2.5.
        ("C=C(C)c1ccccc1", 104.15, 8, 118.18 / 104.15, 9 / 8.0),
    ])
    def test_dp2_and_ams_volatiles_classify_gas_without_needing_veto(
            self, vol_smiles, pool_mw, pool_heavy, mass_units, heavy_units,
            caplog):
        """r89 dual-axis melt gate (PP run-9 fix, RED-first at 5bba9e3eb): a
        proxy-TAGGED, gas-masked, UN-vetoed genuine volatile whose MW clears
        the pre-r89 chain window but which sits below
        _IMPOSTOR_DISCRETE_MONOMER_UNITS (2.5) monomer-equivalents on the
        mass and/or heavy-atom axis must classify GAS -- the tag branch is
        proxy AND NOT veto AND dual-axis polymer-sized, so the unpaired
        shape ``A <=> VOL + TAIL`` stays benign (VOL's gas reference state
        is simply CORRECT). RED at 5bba9e3eb: the window-based tag branch
        melt-classifies VOL and initialize_model dies on the cliff-sign
        ValueError (run-9's death on ``allyl + allyl <=> 1,5-hexadiene``,
        RMG.log:2564-2587)."""
        import logging
        # LIVENESS PINS: VOL is genuinely the r89 shape -- above the
        # pre-r89 window, below the dual-axis threshold.
        mw_vol = Molecule().from_smiles(vol_smiles).get_molecular_weight() * 1000.0
        assert mw_vol >= pool_mw + 10.0, "not the window-leak shape"
        assert mass_units < 2.5 or heavy_units < 2.5, "not a DP-2 volatile"
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            rs = self._volatile_vs_pool_build(vol_smiles, pool_mw, pool_heavy,
                                              set_veto=False)
        assert rs.reference_state_max_decades < 3.0, (
            "the dual-axis gate must classify a below-threshold tagged "
            "volatile GAS (no veto needed) so the reference-state term stays "
            "paired/benign"
        )
        # A decided-False verdict is NOT undecidable: no axis census.
        assert not any("AXIS UNDECIDABLE" in r.getMessage()
                       for r in caplog.records)
        # And not veto-censused either (nothing was vetoed).
        assert not getattr(rs, "gas_veto_census", [])

    def test_gas_volatile_veto_excludes_from_melt_tag_branch(self, caplog):
        """Spec §5.1 (durable gas-volatile veto), r89 form: a POLYMER-SIZED
        proxy-derived discrete (the FR1 run-3 adduct shape, C21H42Br2,
        454.37 g/mol / 23 heavy vs the propene pool -- a genuine multi-repeat
        chain-scale defect clearing the r95 absolute floor) reaches the solver gas-masked but
        proxy-TAGGED. Because is_polymer_proxy is a monotonic multi-writer
        sticky cache with no gas authority, the durable verdict lives in
        ``Species.props['polymer_reference_state_gas_veto']`` (set once at
        the discrete-product creation point, copied by Species.copy, never
        touched by the proxy stamping machinery). The gate must honor it
        with UNCHANGED precedence (r89 RED 4): TAG branch = proxy AND NOT
        veto AND dual-axis polymer-sized.

        The no-veto control refuses in BOTH pre/post r89 (r89 RED 3 at the
        tripwire seam: a genuinely chain-scale tagged unvetoed discrete
        REMAINS melt -- the dual-axis amendment must not blanket-gas the
        FR1 adduct cohort)."""
        import logging

        # LIVENESS PIN + r89 RED 3 (negative control): without the veto the
        # polymer-sized adduct is MELT and the scenario genuinely refuses
        # (proves the volatile-as-melt asymmetry is the cause and pins that
        # r87 tripwire behavior on chain-scale adducts is unchanged).
        with pytest.raises(ValueError, match="unpaired reference-state"):
            self._volatile_vs_pool_build(
                "CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", 42.08, 3,
                set_veto=False)

        # r89 RED 4 (veto precedence unchanged): WITH the durable veto the
        # adduct is excluded from the melt sum -> the build SUCCEEDS.
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            rs = self._volatile_vs_pool_build(
                "CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", 42.08, 3,
                set_veto=True)
        assert rs.reference_state_max_decades < 3.0, (
            "durable gas-volatile veto must exclude the discrete from the "
            "melt tag branch so the reference-state term stays paired/benign"
        )

    def test_gas_veto_suppression_of_chain_scale_member_is_logged(self, caplog):
        """Backstop (code-review IMPORTANT #1, r89-aligned): the veto
        silences the reference-state tripwire for a genuinely CHAIN-SCALE
        (dual-axis polymer-sized) product, because create_reacted_copy
        returns None both for a genuine gas volatile AND for a wing-match
        FAILURE on a real chain-scale fragment. A handshake false-None on a
        genuine chain would be silently dropped from the melt sum instead of
        loudly refused, so the gate must emit a visible census/warning
        naming any POLYMER-SIZED member whose melt classification the veto
        suppressed. r89 alignment: the census keys on the SAME dual-axis
        chain-scale notion as the gate -- a vetoed genuine volatile below
        scale (DP-2 hexadiene, 1.95 propene-monomer-equivalents) is NO
        LONGER announced (pre-r89 its MW cleared the window, 82.15 > 52.1,
        and it was censused despite the warning text itself calling such
        exclusions correct)."""
        import logging
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO_KEY

        a = _spc("CCCC", "A")
        a_mu0, a_mu1, a_mu2 = (_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"),
                               _spc("C#N", "A_mu2"))
        tail = _spc("CCCCCCCC", "TAIL")
        # chain-scale vetoed member: the FR1 adduct shape vs propene pool
        # (C21H42Br2, 454.37 g/mol / 23 heavy -- clears the r95 absolute floor)
        vol = _spc("CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", "VOL")
        vol.is_polymer_proxy = True
        vol.props[VETO_KEY] = True
        # below-scale vetoed volatile: DP-2 hexadiene (1.95 mass-units, 2.0
        # heavy-units vs the propene pool -- decided False on both axes)
        hexd = _spc("C=CCCC=C", "HEXD")
        hexd.is_polymer_proxy = True
        hexd.props[VETO_KEY] = True
        for s in (a, a_mu0, a_mu1, a_mu2, tail, vol, hexd):
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [a, a_mu0, a_mu1, a_mu2, tail, vol, hexd]
        mask = np.array([False] * 5 + [True, True], dtype=bool)
        rxn = Reaction(reactants=[a], products=[vol, tail], **_REV_KIN)

        with caplog.at_level(logging.WARNING):
            rs = _refstate_rs(core, [rxn], mask,
                              [_gate_pool_config(monomer_mw_g_mol=42.08,
                                                 monomer_heavy_atoms=3)],
                              {"A": (1.0, 5.0, 30.0)})

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
        # r89 alignment pin: the below-scale vetoed volatile is NOT censused
        # (its exclusion is a decided gas verdict, not a suppression worth a
        # human's attention). Pre-r89 RED: HEXD's MW cleared the window and
        # it polluted the census.
        assert not any("HEXD" in str(entry)
                       for entry in getattr(rs, "gas_veto_census", [])), (
            "a below-polymer-sized vetoed volatile must not pollute the "
            "gas-veto census (dual-axis alignment, r89)"
        )

    def test_run9_allyl_recombination_dp2_product_passes_tripwire(self, caplog):
        """r89 RED 1, the EXACT PP run-9 death shape
        (/home/alon/Projects/polymer/PP/rmg/run9/RMG.log:2564-2587): the
        all-gas reversible recombination ``[CH2]C=C + [CH2]C=C <=>
        C=CCCC=C`` where the 1,5-hexadiene product (82.15 g/mol, a genuine
        DP-2 gas volatile at 2.0 propene-monomer-equivalents) is
        proxy-tagged (family.py:1657 blanket over-tagging from some OTHER
        proxy-touching reaction) and un-vetoed (no polymer reactant in THIS
        row, so _handshake_structures never stamped it). RED at 5bba9e3eb:
        the window conjunct (82.15 >= 42.08 + 10) melt-classifies hexadiene,
        the allyls (41.07 g/mol) stay gas, and the tripwire dies on a false
        11.01-decade unpaired term. Post-r89 the row has NO melt participant
        at all (U contribution 0) and the build completes silently."""
        import logging
        a = _spc("CCCC", "A")
        a_mu0, a_mu1, a_mu2 = (_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"),
                               _spc("C#N", "A_mu2"))
        allyl = _spc("[CH2]C=C", "allyl")
        hexadiene = _spc("C=CCCC=C", "hexadiene")
        hexadiene.is_polymer_proxy = True   # blanket over-tagging fingerprint
        for s in (a, a_mu0, a_mu1, a_mu2, allyl, hexadiene):
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [a, a_mu0, a_mu1, a_mu2, allyl, hexadiene]
        mask = np.array([False] * 4 + [True, True], dtype=bool)
        rxn = Reaction(reactants=[allyl, allyl], products=[hexadiene],
                       **_REV_KIN)

        # LIVENESS PINS: the run-9 numbers. Hexadiene clears the pre-r89 PP
        # window (chain window 42.08 + 10 = 52.1) but is 1.95 mass-units /
        # 2.0 heavy-units < 2.5; un-vetoed; the row is reversible all-gas.
        mw_hex = hexadiene.molecule[0].get_molecular_weight() * 1000.0
        assert mw_hex == pytest.approx(82.15, abs=0.1)
        assert mw_hex >= 42.08 + 10.0
        assert mw_hex < 2.5 * 42.08
        assert not hexadiene.props.get("polymer_reference_state_gas_veto",
                                       False)
        assert rxn.reversible

        with caplog.at_level(logging.WARNING):
            rs = _refstate_rs(core, [rxn], mask,
                              [_gate_pool_config(monomer_mw_g_mol=42.08,
                                                 monomer_heavy_atoms=3)],
                              {"A": (1.0, 5.0, 30.0)})
        # No melt participant in the row: U identically 0, no census, and
        # the row is untouched (not refused).
        assert rs.reference_state_max_decades == 0.0
        assert not rs.reference_state_census
        assert not getattr(rxn, "polymer_refused", False)
        assert not any("THERMO REFERENCE-STATE" in r.getMessage()
                       for r in caplog.records)

    def test_condensed_branch_melt_below_threshold_unchanged(self):
        """r89 RED 5 (negative control): a genuinely-melt species that is
        condensed via gas_species_mask=False stays MELT even though it sits
        BELOW 2.5 monomer-equivalents -- the dual-axis amendment touches
        ONLY the tag branch; the not-gas-masked branch is size-blind
        (pool-configured by input). Fixture: the §8.1 refusal shape (butane
        proxy A, 58.12 g/mol = 1.38 mass-units / 1.33 heavy-units vs the
        propene pool) must STILL refuse on the unpaired melt <=> all-gas
        row."""
        sp, core, mask = _refstate_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]],
                       **_REV_KIN)
        # LIVENESS PIN: A is below the dual-axis threshold on both axes --
        # only the condensed branch can classify it melt.
        mw_a = sp["A"].molecule[0].get_molecular_weight() * 1000.0
        assert mw_a < 2.5 * 42.08
        assert not mask[0]
        with pytest.raises(ValueError, match="unpaired reference-state"):
            _refstate_rs(core, [rxn], mask,
                         [_gate_pool_config(monomer_mw_g_mol=42.08,
                                            monomer_heavy_atoms=3)],
                         {"A": (1.0, 5.0, 30.0)})

    def test_axis_undecidable_candidate_classifies_gas_and_warns(self, caplog):
        """r89 undecidable-axis contract (mirrors rmgpy.polymer
        ._warn_impostor_axis_undecidable): a proxy-tagged gas-masked
        candidate whose polymer-sized verdict cannot be established on BOTH
        axes (here: mass axis clears the threshold but the pool config
        carries no monomer_heavy_atoms -> heavy denominator uncomputable) is
        NEVER melt-classified blind -- it degrades to conservative-gas
        (build completes; a single-axis melt call is forbidden) AND the
        degradation is announced through the census warning + the
        reference_state_axis_undecidable record."""
        import logging
        # The chain-scale adduct (454.37 g/mol, clears the r95 mass floor) vs a
        # pool WITHOUT a heavy denominator: mass axis computable and above the
        # floor (True), heavy axis uncomputable -> undecidable.
        with caplog.at_level(logging.WARNING):
            rs = self._volatile_vs_pool_build(
                "CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", 42.08,
                0, set_veto=False)
        # conservative-gas: the build completes, no refusal
        assert rs.reference_state_max_decades < 3.0
        # ...and the degradation is announced, naming species and axis
        assert any("AXIS UNDECIDABLE" in r.getMessage()
                   and "VOL" in r.getMessage()
                   and "structure" in r.getMessage()
                   for r in caplog.records), (
            "an undecidable dual-axis verdict must be censused, never silent"
        )
        assert any(lbl == "VOL" and "structure" in miss
                   for lbl, miss in
                   getattr(rs, "reference_state_axis_undecidable", [])), (
            "the undecidable candidate must be recorded on the engine"
        )

    def test_all_gas_no_proxy_reactions_stay_ignored(self):
        """r89 RED 7 (negative control): an ordinary reversible all-gas
        reaction with NO proxy tag anywhere has no melt participant and is
        ignored by the tripwire -- no census, no refusal, max U untouched."""
        a = _spc("CCCC", "A")
        a_mu0, a_mu1, a_mu2 = (_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"),
                               _spc("C#N", "A_mu2"))
        m1, m4 = _spc("C", "M1"), _spc("CCCC", "M4")
        for s in (a, a_mu0, a_mu1, a_mu2, m1, m4):
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [a, a_mu0, a_mu1, a_mu2, m1, m4]
        mask = np.array([False] * 4 + [True, True], dtype=bool)
        # Wildly mass-unpaired all-gas row: only melt membership could make
        # the tripwire care.
        rxn = Reaction(reactants=[m1], products=[m4], **_REV_KIN)
        rs = _refstate_rs(core, [rxn], mask,
                          [_gate_pool_config(monomer_mw_g_mol=42.08,
                                             monomer_heavy_atoms=3)],
                          {"A": (1.0, 5.0, 30.0)})
        assert rs.reference_state_max_decades == 0.0
        assert not rs.reference_state_census

    def test_solver_dual_axis_gate_agrees_with_polymer_module_classifier(self):
        """r89 agreement pin: the solver seam's _dual_axis_polymer_sized
        (numbers-only mirror) and rmgpy.polymer._discrete_is_polymer_sized
        (structure-based original) must return the SAME verdict on a shared
        fixture matrix, with the SAME threshold constant -- the refusal
        mirror (shape B) and the tripwire's is_melt tag branch are built on
        these two implementations respectively, and they MUST agree or the
        refused/live sets drift."""
        from rmgpy.polymer import (Polymer, _IMPOSTOR_DISCRETE_MONOMER_UNITS,
                                   _discrete_is_polymer_sized)
        from rmgpy.solver.polymer import _dual_axis_polymer_sized

        pp = Polymer(label="polypropylene", monomer="[CH2][CH]C",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        ps = Polymer(label="polystyrene", monomer="[CH2][CH]c1ccccc1",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        matrix = [
            ("C=CC", pp), ("C=CCCCC", pp), ("C=CCCC=C", pp),   # 1.0/2.0/1.95u
            ("C=C(C)c1ccccc1", ps),                            # AMS 1.13u
            ("CC(CBr)CC(C)CC(C)Br", pp),                       # 286 g/mol: below r95 floor
            ("[H][H]", pp), ("BrBr", pp),                      # heavy-vs-mass split
            ("CCC(C)CCCC(C)CCCC(C)C", pp),                     # C15 212 g/mol: below r95 floor
            # Genuine chain-scale (454 g/mol / 23 heavy): clears the r95
            # absolute floor -> both gates must AGREE on a TRUE verdict too.
            ("CC(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br", pp),
        ]
        for smiles, pool in matrix:
            spc = Species(molecule=[Molecule().from_smiles(smiles)])
            spc.label = smiles
            expected = _discrete_is_polymer_sized(spc, pool)
            mol = spc.molecule[0]
            mono = pool.monomer
            axes = [(mono.get_molecular_weight() * 1000.0,
                     mono.get_num_atoms() - mono.get_num_atoms('H'))]
            got, missing = _dual_axis_polymer_sized(
                mol.get_molecular_weight() * 1000.0,
                mol.get_num_atoms() - mol.get_num_atoms('H'),
                axes, _IMPOSTOR_DISCRETE_MONOMER_UNITS)
            assert got == expected, (
                f"{smiles} vs {pool.label}: solver gate {got} != polymer "
                f"module classifier {expected} -- refusal mirror and "
                f"tripwire is_melt drift"
            )
            assert missing is None  # fully-computable matrix: all decided

    def test_heavy_atom_key_literal_matches_solver_gate(self):
        """Literal-drift guard (same posture as
        test_gas_veto_key_literal_matches_solver_gate): the Cython solver
        melt gate reads the consumer-world heavy-atom props key as a
        HARDCODED string literal. If the Python constant is renamed without
        updating the pyx, structureless consumer species silently lose their
        heavy axis and every tag-branch candidate degrades to
        conservative-gas (unpairing genuine chains). Pin the exact literal
        so that rename fails loudly here."""
        from rmgpy.polymer import POLYMER_HEAVY_ATOM_COUNT_KEY
        assert POLYMER_HEAVY_ATOM_COUNT_KEY == "polymer_heavy_atom_count", (
            "the heavy-atom props key literal is hardcoded in polymer.pyx's "
            "melt gate; update both together"
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

    Fixture geometry mirrors the real deck (T = 1000 K): proxy C22H46
    (310.61 g/mol, GAV), same-length radical C22H45 (309.60 g/mol, GAV,
    gas-classified + proxy-tagged; genuine multi-repeat chain clearing the r95
    absolute floor -- 22 heavy atoms, 309.6 g/mol -- vs the C5H10 repeat unit),
    H/H2 library. Counterparty window = 70 (pool monomer) + 10 = 80 g/mol: Erad
    is inside (|dMW| = 1.0), H/H2 far outside (|dMW| >= 308)."""

    @staticmethod
    def _build():
        # C22 same-length pair (r95): both clear the absolute chain-scale floor
        # (310.61 / 309.60 g/mol, 22 heavy), still mass-paired (|dMW| = 1.0).
        c22h46 = "CCCCCCCCCCCCCCCCCCCCCC"
        c22h45 = "[CH2]CCCCCCCCCCCCCCCCCCCCC"
        sp = {
            "E": _spc(c22h46, "E"),
            "E_mu0": _spc("CO", "E_mu0"), "E_mu1": _spc("C=O", "E_mu1"), "E_mu2": _spc("C#N", "E_mu2"),
            "Erad": _spc(c22h45, "Erad"),
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
        pool = dataclasses.replace(_gate_pool_config(monomer_mw_g_mol=70.0,
                                                     monomer_heavy_atoms=5),
                                   label="E")
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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
        (not just 'no exception'); in fact 1.5*log10(310.61/309.60) =
        0.0021 decades. The census is empty."""
        import logging
        with caplog.at_level(logging.WARNING):
            rs = self._build()
        assert rs.reference_state_max_decades <= 0.33 + 0.1
        assert rs.reference_state_max_decades == pytest.approx(0.0021, abs=0.002)
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
        assert rs.reference_state_max_decades == pytest.approx(0.0021, abs=0.002)
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
        allow_default_prospective_edge=True,
        # r71 FIX 4 escape (direct-test posture): item-17 gate fixtures
        # deliberately run unstamped proxy rows on the legacy mu1-only
        # path; the production hard-fail default is pinned by
        # TestRefusalAdjudicationSurvivesRebuild.
        allow_unstamped_proxy_rows=True)
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
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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
        rs_edge = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
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
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=np.asarray(phase.get_gas_mask(core), dtype=bool),
            constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_classifier=phase.get_gas_mask)
        rs.initialize_model(list(core), [], [g], [rxn])
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert pm.shape[0] == rs.num_core_species + 1


class TestGasMaskSmilesBaseAliasParity:
    """PP run-5 crash pin (PROSPECTIVE-MASK TRIPWIRE, 2026-07-04): _base_label
    truncated labels at the FIRST '(' (label.partition('(')[0]), a convention
    meant for the RMG index suffix ('PS(2)' -> 'PS'). SMILES-derived labels
    legitimately contain '(' as BRANCHING syntax, so truncation aliases
    structurally different species onto one base:

      edge C9H19 daughter 'C[CH]CC(C)CC(C)C'  -> base 'C[CH]CC'
      core C6H13 tail     'C[CH]CC(C)C'       -> base 'C[CH]CC'  (COLLIDES)

    A genuine H-loss radical daughter of the C9H20 proxy in the EDGE therefore
    registers a truncated base that the get_gas_mask section-D application
    loop matches against unrelated CORE scission-tail radicals -- but ONLY in
    the combined chain(core, edge) call (core-only has no qualifying C9H19),
    so prospective says CONDENSED where gas_species_mask says GAS and RIDER
    R1 raises (run 5: indices 24-26, C[CH]CC(C)C / [CH2]C(C)C /
    [CH2]C(C)CCC). The core mask is RIGHT (the C6/C4 tails are not H-loss
    daughters of any proxy); the prospective mask lied via the alias.

    Fix under test: base labels strip ONLY a trailing '(<int>)' index suffix.
    """

    def _phase_and_lists(self):
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import PolymerPhase, PolymerPool

        mono = _spc("C=CC", "MON")
        mu0 = _spc("CO", "polypropylene_mu0")
        mu1 = _spc("C=O", "polypropylene_mu1")
        mu2 = _spc("C#N", "polypropylene_mu2")
        for m in (mu0, mu1, mu2):
            m.is_moment_dummy = True
        proxy = _spc("CC(C)CC(C)CCC", "polypropylene")   # C9H20 trimer proxy
        pool = PolymerPool(label="polypropylene", xs=3, monomer=mono,
                           explicit_map={}, mu_species=[mu0, mu1, mu2],
                           proxy_species=proxy)
        phase = PolymerPhase(density=Quantity(905.0, "kg/m^3"),
                             initial_moments={"polypropylene": (1.0, 5.0, 30.0)},
                             initial_explicit={}, pools=[pool])
        n2 = _spc("N#N", "N2")
        # CORE scission-tail radicals, labelled by their SMILES (the run-5
        # shape: labels contain structural parentheses, no index suffix).
        c6a = _spc("C[CH]CC(C)C", "C[CH]CC(C)C")         # C6H13
        c4 = _spc("[CH2]C(C)C", "[CH2]C(C)C")            # C4H9
        c6b = _spc("[CH2]C(C)CCC", "[CH2]C(C)CCC")       # C6H13
        core = [n2, proxy, mono, mu0, mu1, mu2, c6a, c4, c6b]
        # EDGE: genuine C9H19 H-loss radical daughters of the C9H20 proxy
        # whose SMILES labels share a first-'(' prefix with the core tails.
        d1 = _spc("C[CH]CC(C)CC(C)C", "C[CH]CC(C)CC(C)C")     # base 'C[CH]CC'
        d2 = _spc("[CH2]C(C)CC(C)CCC", "[CH2]C(C)CC(C)CCC")   # base '[CH2]C'
        edge = [d1, d2]
        return phase, core, edge

    def test_edge_daughters_qualify_and_tails_do_not(self):
        """Sanity of the fixture: the C9H19 edge daughters DO qualify under
        the narrow H-loss predicate (they are the real Gate-B beneficiaries);
        the C6/C4 core tails do NOT."""
        from rmgpy.rmg.polymer_input import _base_label
        phase, core, edge = self._phase_and_lists()
        bases = phase.get_h_loss_radical_daughter_bases(core + edge)
        # Every returned base must be the base of a QUALIFYING species (the
        # two C9H19 daughters), never of a core tail.
        daughter_bases = {_base_label(s.label) for s in edge}
        assert bases, "fixture lost its qualifying C9H19 daughters"
        assert bases <= daughter_bases
        # core-only: nothing qualifies
        assert phase.get_h_loss_radical_daughter_bases(core) == set()

    def test_gas_mask_core_prefix_parity_under_smiles_alias(self):
        """THE run-5 pin: get_gas_mask(chain(core, edge))[:n_core] must equal
        get_gas_mask(core). On the buggy _base_label the C9H19 edge daughters'
        truncated bases alias the three core tails CONDENSED in the combined
        call only (exactly the PROSPECTIVE-MASK TRIPWIRE divergence)."""
        phase, core, edge = self._phase_and_lists()
        core_only = np.asarray(phase.get_gas_mask(core), dtype=bool)
        combined = np.asarray(phase.get_gas_mask(core + edge), dtype=bool)
        # the three tails are GAS in the core mask (the RIGHT verdict)
        assert bool(core_only[6]) and bool(core_only[7]) and bool(core_only[8])
        # core-prefix parity (R1's invariant) -- RED before the fix
        assert np.array_equal(combined[:len(core)], core_only), (
            "combined get_gas_mask core prefix diverged: %s vs %s"
            % (combined[:len(core)].tolist(), core_only.tolist()))

    def test_genuine_edge_daughters_stay_condensed(self):
        """Guard against 'fixing' the alias by killing branch 2: the genuine
        C9H19 daughters must classify CONDENSED in the combined call, and
        their bases must survive an RMG index suffix on the label."""
        from rmgpy.rmg.polymer_input import _base_label
        phase, core, edge = self._phase_and_lists()
        combined = np.asarray(phase.get_gas_mask(core + edge), dtype=bool)
        assert not bool(combined[len(core)])      # d1 condensed
        assert not bool(combined[len(core) + 1])  # d2 condensed
        # index-suffixed copies (RMG relabels on promotion: 'label(31)')
        bases = phase.get_condensed_edge_daughter_bases(core + edge)
        for s in edge:
            assert _base_label(s.label + "(31)") in bases

    def test_base_label_strips_only_trailing_index_suffix(self):
        """The corrected base-label convention: strip ONLY a trailing
        '(<int>)' RMG index suffix; structural SMILES parentheses survive."""
        from rmgpy.rmg.polymer_input import _base_label
        assert _base_label("PS(2)") == "PS"
        assert _base_label("polypropylene(2)") == "polypropylene"
        assert _base_label("polypropylene_mod_2") == "polypropylene_mod_2"
        assert _base_label("C[CH]CC(C)C(6)") == "C[CH]CC(C)C"
        assert _base_label("C[CH]CC(C)C") == "C[CH]CC(C)C"
        assert _base_label("[CH2]C(C)CCC") == "[CH2]C(C)CCC"


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


class TestHLossRadicalDaughterQualifier:
    """Ratified 2026-07-03 (PP thin-core fix): Gate B zeroed every proxy
    H-abstraction (H + C9H20-proxy <=> H2 + C9H19*) because the edge-daughter
    classifier only recognized spawned Polymer proxies with mu-triplets --
    so the PP core could never grow (thin-core incident, 14 spc / 2 rxn).
    NARROW qualifier: a plain Species edge daughter is prospectively
    condensed ONLY when ALL hold: (i) plain Species (not a Polymer proxy);
    (ii) radical-bearing and neutral; (iii) same non-H element composition
    as a configured condensed pool proxy; (iv) differs from that proxy by
    H-loss only; (v) not durably gas-vetoed and not a smaller
    scission/volatile fragment. MW-window-alone qualification is REJECTED
    (near-monomer volatiles would misclassify)."""

    PP_PROXY_SMILES = "CCCC(C)CC(C)C"        # C9H20 PP trimer proxy
    C9H19_SMILES = "CCCC(C)CC(C)[CH2]"       # H-loss radical daughter

    @staticmethod
    def _phase_with_pp():
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import PolymerPhase, PolymerPool
        pool = PolymerPool(
            label="PP", xs=2,
            monomer=_spc("C=CC", "PP_repeat"),
            explicit_map={}, mu_species=[])
        return PolymerPhase(density=Quantity(1000.0, "kg/m^3"),
                            initial_moments={}, initial_explicit={},
                            pools=[pool], mass_transfer=[])

    def _combined(self, *extra):
        # the configured pool's core proxy (base label == pool label) plus
        # the candidate daughters under test
        proxy = _spc(self.PP_PROXY_SMILES, "PP(2)")
        return [proxy] + list(extra)

    def test_h_loss_radical_daughter_qualifies(self):
        """THE PP red pin: the C9H19 H-abstraction daughter (plain Species,
        radical, neutral, same C9 skeleton, proxy-H minus one) qualifies."""
        phase = self._phase_with_pp()
        d = _spc(self.C9H19_SMILES, "C9H19a(15)")
        bases = phase.get_condensed_edge_daughter_bases(self._combined(d))
        assert "C9H19a" in bases

    def test_gas_product_propene_does_not_qualify(self):
        """A gas product (propene: closed-shell, C3 != C9) must stay
        unqualified -- Gate B keeps zeroing only-gas poly events."""
        phase = self._phase_with_pp()
        d = _spc("C=CC", "propene(20)")
        assert phase.get_condensed_edge_daughter_bases(
            self._combined(d)) == set()

    def test_smaller_volatile_fragment_radical_does_not_qualify(self):
        """A smaller scission/volatile fragment radical (C3H7*: radical and
        neutral but NOT the proxy's heavy-atom skeleton) must not qualify --
        the composition lock is what rejects MW-window-style qualification."""
        phase = self._phase_with_pp()
        d = _spc("[CH2]CC", "C3H7(21)")
        assert phase.get_condensed_edge_daughter_bases(
            self._combined(d)) == set()

    def test_closed_shell_isomer_does_not_qualify(self):
        """A closed-shell C9H20 isomer (same composition as the proxy, zero
        radicals) is not an H-loss daughter -- (ii)/(iv) both fail."""
        phase = self._phase_with_pp()
        d = _spc("CCCCCC(C)CC", "C9H20iso(22)")
        assert phase.get_condensed_edge_daughter_bases(
            self._combined(d)) == set()

    def test_gas_vetoed_daughter_does_not_qualify(self):
        """A durably gas-vetoed species must stay gas even when it matches
        the H-loss shape (condition (v))."""
        from rmgpy.polymer import set_polymer_gas_veto
        phase = self._phase_with_pp()
        d = _spc(self.C9H19_SMILES, "C9H19v(23)")
        set_polymer_gas_veto(d)
        assert phase.get_condensed_edge_daughter_bases(
            self._combined(d)) == set()

    def test_spawned_proxy_without_triplet_still_unqualified(self):
        """The spawned-proxy branch is untouched: a Polymer daughter without
        its mu-triplet stays unqualified exactly as before (it must NOT leak
        through the new plain-Species branch)."""
        from rmgpy.polymer import Polymer
        phase = self._phase_with_pp()
        d = Polymer(label="PP_scission_tail", monomer="[CH2][CH]C",
                    Mn=100.0, Mw=200.0)
        d.is_polymer_proxy = True
        bases = phase.get_condensed_edge_daughter_bases(self._combined(d))
        assert "PP_scission_tail" not in bases

    def test_sticky_proxy_stamp_alone_still_never_qualifies(self):
        """Safety property preserved: a random ordinary species carrying a
        sticky is_polymer_proxy stamp (family.py over-stamping) does not
        qualify unless it independently passes the narrow H-loss rule."""
        phase = self._phase_with_pp()
        sticky = _spc("CCO", "OTH")
        sticky.is_polymer_proxy = True
        assert phase.get_condensed_edge_daughter_bases(
            self._combined(sticky)) == set()


class TestGateBHLossRadicalDaughterFlux:
    """Solver-level PP-shaped pin: with the REAL bound classifier, the proxy
    H-abstraction edge reaction (H + PP <=> H2 + C9H19*) must pass Gate B
    (daughter prospectively condensed) and carry nonzero real flux, while an
    only-gas poly event (propene product) stays Gate-B zeroed -- the hole
    Gate B closes must NOT reopen."""

    PP_PROXY_SMILES = "CCCC(C)CC(C)C"
    C9H19_SMILES = "CCCC(C)CC(C)[CH2]"

    def _build(self, daughter):
        phase = TestHLossRadicalDaughterQualifier._phase_with_pp()
        proxy = _spc(self.PP_PROXY_SMILES, "PP")
        mu0 = _spc("CO", "PP_mu0")
        mu1 = _spc("C=O", "PP_mu1")
        mu2 = _spc("C#N", "PP_mu2")
        h = _spc("[H]", "H")
        h2 = _spc("[H][H]", "H2")
        core = [proxy, mu0, mu1, mu2, h, h2]
        mask = np.array([False, False, False, False, True, True], dtype=bool)
        rxn = Reaction(reactants=[h, proxy], products=[h2, daughter], **_KIN)
        cfg = PolymerPoolConfig(label="PP", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None,
                                k_unzip=0.0)
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
            T=800.0, P=1.0e5, initial_mole_fractions={h: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"PP": (1.0, 5.0, 30.0)}, termination=[],
            prospective_condensed_edge_daughter_classifier=(
                phase.get_condensed_edge_daughter_bases),
            allow_default_prospective_edge=True)
        rs.initialize_model(core, [], [daughter], [rxn])
        return rs, core

    def test_h_abstraction_daughter_passes_gate_b_flux_nonzero(self):
        """RED pin (PP incident): pre-fix this poly event is Gate-B zeroed
        (gate_code 2, flux identically 0); post-fix the C9H19 daughter is
        prospectively condensed, the gate opens and the REAL flux flows."""
        d = _spc(self.C9H19_SMILES, "C9H19a(15)")
        rs, core = self._build(d)
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert bool(pm[len(core)]) is False    # edge daughter CONDENSED
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert int(np.asarray(rs.edge_reaction_gate_code)[0]) == 0
        assert float(np.asarray(rs.edge_reaction_rates)[0]) > 0.0

    def test_only_gas_poly_event_remains_gate_b_zeroed(self):
        """Converse pin: a poly event whose only products are genuinely gas
        (H2 + propene) must remain Gate-B zeroed with a live counterfactual
        -- the fix must not reopen the phase-gate hole."""
        d = _spc("C=CC", "propene(20)")
        rs, core = self._build(d)
        pm = np.asarray(rs.prospective_gas_mask, dtype=bool)
        assert bool(pm[len(core)]) is True     # propene stays GAS
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert int(np.asarray(rs.edge_reaction_gate_code)[0]) == 2
        assert float(np.asarray(rs.edge_reaction_rates)[0]) == 0.0
        assert float(np.asarray(rs.edge_reaction_rates_ungated)[0]) > 0.0


class TestHLossRadicalDaughterCoreMask:
    """P1-A (ratified core-mask extension of the 2026-07-03 H-loss
    qualifier): the a10acaed4 qualifier only affected EDGE slots. Once a
    qualifying H-loss radical daughter (e.g. C9H19 from the PP proxy) is
    PROMOTED TO CORE, get_gas_mask stage-1/2 still classified it GAS and
    static Gate B zeroed its core production rows -- PP regeneration would
    stall one step after promotion. The SAME narrow qualifier (plain
    Species, radical-bearing, neutral, same non-H element composition as a
    configured pool proxy, H-loss-only, not durably gas-vetoed) must feed
    the CORE mask construction path (PolymerPhase.get_gas_mask), so a
    qualifying daughter is condensed-classified consistently in core and
    edge. Negative pins: gas products, volatile fragments, spawned
    proxies, and gas-vetoed species classify exactly as today, in core
    too."""

    PP_PROXY_SMILES = "CCCC(C)CC(C)C"        # C9H20 PP trimer proxy
    C9H19_SMILES = "CCCC(C)CC(C)[CH2]"       # H-loss radical daughter

    @staticmethod
    def _phase():
        return TestHLossRadicalDaughterQualifier._phase_with_pp()

    # ---- predicate level: get_gas_mask over a CORE list ----

    def test_promoted_daughter_condensed_by_get_gas_mask(self):
        """THE core red pin: the qualifying C9H19 daughter, now a CORE
        species, must be condensed-classified by get_gas_mask itself (the
        core mask construction path), not just by the edge-slot flip."""
        phase = self._phase()
        core = [_spc(self.PP_PROXY_SMILES, "PP(2)"),
                _spc(self.C9H19_SMILES, "C9H19a(15)"),
                _spc("[H]", "H"), _spc("[H][H]", "H2")]
        mask = phase.get_gas_mask(core)
        assert bool(mask[1]) is False          # daughter CONDENSED in core
        assert bool(mask[2]) is True           # H stays gas
        assert bool(mask[3]) is True           # H2 stays gas

    def test_negative_pins_stay_gas_in_core_mask(self):
        """Gas products, smaller volatile fragments, closed-shell isomers
        and gas-vetoed species must classify exactly as today in core."""
        from rmgpy.polymer import set_polymer_gas_veto
        phase = self._phase()
        vetoed = _spc(self.C9H19_SMILES, "C9H19v(23)")
        set_polymer_gas_veto(vetoed)
        core = [_spc(self.PP_PROXY_SMILES, "PP(2)"),
                _spc("C=CC", "propene(20)"),
                _spc("[CH2]CC", "C3H7(21)"),
                _spc("CCCCCC(C)CC", "C9H20iso(22)"),
                vetoed]
        mask = phase.get_gas_mask(core)
        assert [bool(v) for v in mask[1:]] == [True, True, True, True]

    def test_spawned_proxy_classification_unchanged_in_core(self):
        """A spawned Polymer proxy (no configured registration) keeps its
        current get_gas_mask classification (GAS at this stage; stage-2 /
        spawned-pool machinery owns its condensation) -- the H-loss branch
        must not leak onto Polymer instances."""
        from rmgpy.polymer import Polymer
        phase = self._phase()
        d = Polymer(label="PP_scission_tail", monomer="[CH2][CH]C",
                    Mn=100.0, Mw=200.0)
        d.is_polymer_proxy = True
        core = [_spc(self.PP_PROXY_SMILES, "PP(2)"), d]
        mask = phase.get_gas_mask(core)
        assert bool(mask[1]) is True

    # ---- solver level: static Gate B on core rows ----

    def _build_core_system(self, daughter_smiles, daughter_label):
        phase = self._phase()
        proxy = _spc(self.PP_PROXY_SMILES, "PP")
        mu0 = _spc("CO", "PP_mu0")
        mu1 = _spc("C=O", "PP_mu1")
        mu2 = _spc("C#N", "PP_mu2")
        h = _spc("[H]", "H")
        h2 = _spc("[H][H]", "H2")
        d = _spc(daughter_smiles, daughter_label)
        core = [proxy, mu0, mu1, mu2, h, h2, d]
        mask = np.asarray(phase.get_gas_mask(core), dtype=bool)
        rxn = Reaction(reactants=[h, proxy], products=[h2, d], **_KIN)
        cfg = PolymerPoolConfig(label="PP", xs=2,
                                explicit_dp_to_species_index={},
                                mu_indices=(1, 2, 3), monomer_poly_index=None,
                                k_unzip=0.0)
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
            T=800.0, P=1.0e5, initial_mole_fractions={h: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"PP": (1.0, 5.0, 30.0)}, termination=[],
            prospective_classifier=phase.get_gas_mask,
            prospective_condensed_edge_daughter_classifier=(
                phase.get_condensed_edge_daughter_bases))
        rs.initialize_model(core, [rxn], [], [])
        return rs, core

    def test_promoted_daughter_production_row_not_gate_b_zeroed(self):
        """RED pin (P1-A): with the C9H19 daughter IN THE CORE, the proxy
        H-abstraction core row (H + PP <=> H2 + C9H19) was Gate-B zeroed
        (bare continue, production rate identically 0) because the core
        mask still said GAS. Post-fix the daughter is condensed in the
        core mask, the poly event has a condensed product, and its REAL
        production flux flows -- PP regeneration survives promotion."""
        rs, core = self._build_core_system(self.C9H19_SMILES, "C9H19a(15)")
        assert bool(np.asarray(rs.gas_species_mask)[6]) is False
        # R1 core-prefix parity held (initialize_model completed): the
        # prospective mask agrees with the core mask on the daughter.
        assert bool(np.asarray(rs.prospective_gas_mask)[6]) is False
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert float(rs.core_species_production_rates[6]) > 0.0

    def test_only_gas_core_poly_event_stays_gate_b_zeroed(self):
        """Converse pin: a CORE poly event whose only products are
        genuinely gas (H2 + propene) must remain Gate-B zeroed -- the
        core-mask extension must not reopen the phase-gate hole."""
        rs, core = self._build_core_system("C=CC", "propene(20)")
        assert bool(np.asarray(rs.gas_species_mask)[6]) is True
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert float(rs.core_species_production_rates[6]) == 0.0


class TestReversibleChainScissionScopePin:
    """ITEM C -- honest scope pin for the reference-state guard (no
    behavior change; pin + wording only). The 2026-07-03 monomer-gas fix
    removed the monomer-PRODUCT conflation (the routed release target
    being force-condensed). It does NOT make reversible chain-scission
    decks initialize: with the monomer now genuinely gas, a REVERSIBLE
    polymer-proxy scission `proxy <=> gas monomer + gas fragment` (the
    real PS mechanism's `polystyrene(2) <=> C8H8 + C16H18` shape, which
    read U = 0.72 decades pre-fix ONLY because the monomer was wrongly
    condensed) carries exactly one net melt participant and LEGITIMATELY
    refuses via the thermo reference-state tripwire. The resolution is
    deck/prep-level irreversibility -- open-crucible physics (the
    volatiles leave; the reverse reaction has no meaning in the melt) --
    not a code defect."""

    # polystyrene(2) trimer-scale proxy: C24H26 = C8H8 + C16H18
    PROXY_SMILES = "CC(c1ccccc1)CC(c1ccccc1)CCc1ccccc1"   # C24H26
    STYRENE_SMILES = "C=Cc1ccccc1"                        # C8H8
    FRAGMENT_SMILES = "CC(c1ccccc1)CCc1ccccc1"            # C16H18

    def _build(self, **kwargs):
        sp = {
            "PS": _spc(self.PROXY_SMILES, "polystyrene(2)"),
            "PS_mu0": _spc("CO", "polystyrene_mu0"),
            "PS_mu1": _spc("C=O", "polystyrene_mu1"),
            "PS_mu2": _spc("C#N", "polystyrene_mu2"),
            "STY": _spc(self.STYRENE_SMILES, "C8H8"),      # gas monomer
            "FRAG": _spc(self.FRAGMENT_SMILES, "C16H18"),  # gas fragment
        }
        for s in sp.values():
            s.thermo = _trivial_nasa(_LIB_COMMENT)
        core = [sp["PS"], sp["PS_mu0"], sp["PS_mu1"], sp["PS_mu2"],
                sp["STY"], sp["FRAG"]]
        mask = np.array([False, False, False, False, True, True],
                        dtype=bool)
        rxn = Reaction(reactants=[sp["PS"]],
                       products=[sp["STY"], sp["FRAG"]], **_REV_KIN)
        pool = PolymerPoolConfig(
            label="polystyrene", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            monomer_mw_g_mol=104.15, k_scission=0.0, k_unzip=0.0)
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
            T=800.0, P=1.0e5, initial_mole_fractions={sp["STY"]: 0.0},
            V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"polystyrene": (1.0, 5.0, 30.0)},
            termination=[], **kwargs)
        return rs, core, rxn, sp

    def test_reversible_scission_to_gas_products_refuses(self):
        """THE scope pin: one condensed proxy reactant, gas monomer + gas
        fragment products, reversible -- exactly one net melt participant
        -- must REFUSE via the tripwire. This is the guard working as
        specified, not a defect for a future 'fix' to erode: any change
        that lets this shape initialize without a deck-level thermo/
        irreversibility resolution reintroduces the unpaired
        reference-state error the tripwire exists to catch."""
        rs, core, rxn, sp = self._build()

        # LIVENESS: independently recomputed U for the single net melt
        # participant (the proxy) is far above the 3.0-decade refuse
        # bound -- the fixture genuinely carries the pathological shape,
        # so the raises-check below cannot pass vacuously.
        mw_proxy = sp["PS"].molecule[0].get_molecular_weight()
        u_expected = (_sackur_tetrode_decades(mw_proxy, 800.0)
                      + math.log10(1.0e5 / (constants.R * 800.0 * 1.0)))
        assert u_expected > 3.0, (
            "FIXTURE BROKEN: recomputed U for the melt proxy "
            f"({u_expected:.2f}) is not above the refuse bound")
        assert rxn.reversible

        with pytest.raises(ValueError,
                           match=r"THERMO REFERENCE-STATE TRIPWIRE"):
            rs.initialize_model(core, [rxn], [], [])

    def test_census_attributes_refusal_to_the_scission_reaction(self):
        """Documentation pin: under the deck author's explicit
        allow_unpaired_reference_state bypass the same deck initializes
        and the census names THIS reaction with the chain-scale U -- the
        refusal is attributable, physical, and resolved at deck/prep
        level (irreversible scission in an open crucible), never by
        weakening the guard."""
        rs, core, rxn, sp = self._build(allow_unpaired_reference_state=True)
        rs.initialize_model(core, [rxn], [], [])
        assert rs.reference_state_max_decades > 3.0
        assert any(str(rxn) == name
                   for name, _u in rs.reference_state_census)


class TestEdgeReactantReverseFlux:
    """Review round 51 (beta-scission-flux-diagnosis.md): RMG canonically
    stores decomposition channels in the ASSOCIATION direction -- the
    decomposition daughter is an edge REACTANT and the discovery flux is the
    reverse rate rr = kb*[core products]. simple.pyx (:443-459, :509-523)
    evaluates such rows (edge concentration = 0 -> rf vanishes; rr computable
    when ALL products are core) and accumulates edge-reactant flux with
    edge_species_rates[reactant] -= rate, so a net rate of -rr registers
    POSITIVE daughter production. polymer.pyx skipped every edge-reactant row
    before rate evaluation and accumulated the product side only, so
    beta-scission, H-loss and proxy C-C homolysis discovery flux was
    identically zero by construction (PP run 3: every core promotion in 39
    iterations came from all-core-reactant rows)."""

    @staticmethod
    def _association_fixture():
        """Run-3 shape: edge C6H13 + core propene <=> core C9H19."""
        propene = _spc("C=CC", "propene")
        c9 = _spc("C[CH]CC(C)CC(C)C", "C9H19")
        c6 = _spc("C[CH]CC(C)C", "C6H13")
        for s in (propene, c9, c6):
            s.thermo = _trivial_nasa()
        kin = Arrhenius(A=(1.0, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"),
                        T0=(1.0, "K"))
        rxn = Reaction(reactants=[c6, propene], products=[c9],
                       kinetics=kin, reversible=True)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={c9: 1.0}, V_poly=1.0,
            polymer_pools=[], mass_transfer=[],
            gas_species_mask=np.array([True, True], dtype=bool),
            constant_gas_volume=False, termination=[],
            # direct-test build, no blueprint phase: permit the default
            # (GAS) prospective edge fill so R1-EDGE does not raise.
            allow_default_prospective_edge=True)
        rs.initialize_model([propene, c9], [], [c6], [rxn])
        return rs, propene, c9, c6

    def test_association_row_registers_positive_daughter_edge_flux(self):
        """PIN (a): a core-product association row with an edge REACTANT
        (the canonical reverse-stored beta-scission channel) must show
        POSITIVE daughter edge flux equal to rr = kb*[C9H19] -- today it is
        exactly 0 because the row is skipped before rate evaluation."""
        rs, propene, c9, c6 = self._association_fixture()
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        assert rs.kb[0] > 0.0            # fixture liveness: reverse is live
        c9_conc = max(0.0, rs.y[1]) / rs.V_gas
        expected_rr = rs.kb[0] * c9_conc
        assert expected_rr > 0.0

        # rf vanishes (edge reactant has no state => zero concentration),
        # so the net rate is exactly -rr ...
        assert rs.edge_reaction_rates[0] == pytest.approx(-expected_rr)
        # ... and the reactant-side accumulator (-= rate, simple.pyx
        # :509-515 sign parity) registers POSITIVE daughter production.
        assert rs.edge_species_rates[0] > 0.0
        assert rs.edge_species_rates[0] == pytest.approx(expected_rr)
        assert rs.edge_species_rates_ungated[0] == pytest.approx(expected_rr)

    def test_association_row_stays_diagnostic_only(self):
        """PIN (a) rider: the newly-live edge row must stay DIAGNOSTIC-ONLY --
        no dn_dt writes, no core production/consumption side effects."""
        rs, propene, c9, c6 = self._association_fixture()
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        assert np.all(dn_dt == 0.0)
        assert np.all(rs.core_species_consumption_rates == 0.0)
        assert np.all(rs.core_species_production_rates == 0.0)

    def test_get_edge_reaction_rates_handles_edge_reactant_rows(self):
        """Bullet 5: the diagnostic path (get_reaction_rates feeding
        get_edge_reaction_rates) indexes masks/concentrations with raw
        participant slots as if every participant were core -- an edge
        reactant must not crash it, and it must report the same -rr net rate
        as the residual."""
        rs, propene, c9, c6 = self._association_fixture()
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))

        edge_rates = np.zeros(1, float)
        rs.get_edge_reaction_rates(rs.y, edge_rates)
        c9_conc = max(0.0, rs.y[1]) / rs.V_gas
        assert edge_rates[0] == pytest.approx(-rs.kb[0] * c9_conc)
        assert edge_rates[0] == pytest.approx(rs.edge_reaction_rates[0])

    @staticmethod
    def _homolysis_fixture():
        """Reverse-stored proxy recombination: edge C3 + edge C6 <=> core PP
        proxy (the canonical pyrolysis-initiation discovery pattern). The
        daughters are prospectively CONDENSED via the edge-daughter
        condensed-mask classifier (spec 2026-06-29, production shape for
        pool-scission daughters), so the row is a poly event with a condensed
        product and neither phase gate zeroes it."""
        a = _spc("CCC(C)CC(C)C", "A")
        mu0 = _spc("CO", "A_mu0")
        mu1 = _spc("C=O", "A_mu1")
        mu2 = _spc("C#N", "A_mu2")
        g = _spc("[CH3]", "G")
        c3 = _spc("[CH2]CC", "C3(7)")
        c6 = _spc("C[CH]CC(C)C", "C6(6)")
        for s in (a, mu0, mu1, mu2, g, c3, c6):
            s.thermo = _trivial_nasa()
        kin = Arrhenius(A=(1.0, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"),
                        T0=(1.0, "K"))
        rxn = Reaction(reactants=[c3, c6], products=[a],
                       kinetics=kin, reversible=True)
        pool = PolymerPoolConfig(label="A", xs=2,
                                 explicit_dp_to_species_index={},
                                 mu_indices=(1, 2, 3),
                                 monomer_poly_index=None, k_unzip=0.0)
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
            T=800.0, P=1.0e5, initial_mole_fractions={g: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=np.array([False] * 4 + [True], dtype=bool),
            constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
            prospective_condensed_edge_daughter_classifier=
                _condensed_edge_stub({"C3", "C6"}),
            allow_default_prospective_edge=True)
        rs.initialize_model([a, mu0, mu1, mu2, g], [], [c3, c6], [rxn])
        return rs

    def test_edge_edge_reverse_homolysis_scaled_by_product_site(self):
        """PIN (b): edge + edge <=> core PP proxy gives positive edge flux
        ONLY through rr, and rr takes its site scaling from the PRODUCT-side
        proxy pool (default mu1 site density, y[mu1]/V_poly = 5.0) --
        mirroring how forward-stored proxy rows replace _C(proxy) = 1.0 with
        the pool site density."""
        rs = self._homolysis_fixture()
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        assert rs.kb[0] > 0.0            # fixture liveness: reverse is live
        site = max(0.0, rs.y[2]) / rs.V_poly     # product pool mu1 basis
        assert site == pytest.approx(5.0)
        expected_rr = rs.kb[0] * site

        # rf == 0 exactly (BOTH reactants are edge, zero concentration):
        # the flux exists only through rr.
        assert rs.edge_reaction_rates[0] == pytest.approx(-expected_rr)
        for edge_slot in (0, 1):         # C3 and C6 daughters alike
            assert rs.edge_species_rates[edge_slot] > 0.0
            assert rs.edge_species_rates[edge_slot] == \
                pytest.approx(expected_rr)

        # Diagnostic path parity (bullet 5): same product-side site scaling.
        edge_rates = np.zeros(1, float)
        rs.get_edge_reaction_rates(rs.y, edge_rates)
        assert edge_rates[0] == pytest.approx(-expected_rr)

    def test_edge_edge_reverse_homolysis_stays_diagnostic_only(self):
        """PIN (b) rider: no dn_dt writes -- proxy state and pool moments
        untouched by the newly-live edge row."""
        rs = self._homolysis_fixture()
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(dn_dt == 0.0)

    def test_edge_product_reverse_rate_stays_uncomputable(self):
        """PIN (c): a row with an EDGE PRODUCT (no state in y) still has
        rr = 0 -- the concentration-availability hole (simple.pyx :452-453,
        Z6) is the VALID half of the old skip and must be preserved even on
        the newly-live edge-reactant path."""
        a = _spc("C", "A")
        b = _spc("CC", "B")
        e1 = _spc("[CH3]", "E1")
        e2 = _spc("C[CH2]", "E2")
        for s in (a, b, e1, e2):
            s.thermo = _trivial_nasa()
        kin = Arrhenius(A=(1.0, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"),
                        T0=(1.0, "K"))
        # edge E1 + core A <=> edge E2 + core B: rf = 0 (E1 has no state)
        # AND rr = 0 (E2 has no state) => the row contributes exactly nothing.
        rxn = Reaction(reactants=[e1, a], products=[e2, b],
                       kinetics=kin, reversible=True)
        rs = HybridPolymerSystem(  # r71 FIX 4 escape: direct-test fixture
            allow_unstamped_proxy_rows=True,
            T=800.0, P=1.0e5, initial_mole_fractions={a: 0.5, b: 0.5},
            V_poly=1.0, polymer_pools=[], mass_transfer=[],
            gas_species_mask=np.array([True, True], dtype=bool),
            constant_gas_volume=False, termination=[],
            allow_default_prospective_edge=True)
        rs.initialize_model([a, b], [], [e1, e2], [rxn])
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        assert rs.edge_reaction_rates[0] == 0.0
        assert rs.edge_species_rates[0] == 0.0
        assert rs.edge_species_rates[1] == 0.0
        assert rs.edge_species_rates_ungated[0] == 0.0
        assert rs.edge_species_rates_ungated[1] == 0.0
        assert np.all(dn_dt == 0.0)


class TestInitialExplicitPlumbing:
    """Defect fix (explicit-DP stage A prerequisite): the run path passed the
    phase-level ``{Species: moles}`` dissolved-species dict straight into the
    solver ctor arg ``initial_explicit_species``, whose contract is
    ``{pool_label: {dp: moles}}`` (polymer.pyx set_initial_conditions step 2
    looks up ``pool.label in self.initial_explicit_species``). A Species key
    never matches a str pool label, so the channel was dead on arrival: a
    deck-stated explicit-oligomer loading could never reach y0 or split off
    the tail moments. to_solver_object must translate the phase's
    Species-keyed oligomer loadings (the same source calculate_volume reads,
    polymer_input.py explicit_mu1 accounting) into the solver's pool-labeled
    shape -- and must NOT forward the dissolved gases, which already receive
    y0 through initialMoles (step 1)."""

    @staticmethod
    def _reactor_with_explicit_pool(explicit_moles=0.3):
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.polymer_input import (HybridPolymerReactor, PolymerPhase,
                                             PolymerPool)
        n2 = _spc("N#N", "N2")
        p2 = _spc("CC", "P2")  # explicit DP=2 oligomer (condensed via explicit_map)
        mu0 = _spc("CO", "poly_mu0")
        mu1 = _spc("C=O", "poly_mu1")
        mu2 = _spc("C#N", "poly_mu2")
        for m in (mu0, mu1, mu2):
            m.is_moment_dummy = True
        pool = PolymerPool(label="poly", xs=2,
                           monomer=Molecule().from_smiles("C=C"),
                           explicit_map={2: p2}, mu_species=[mu0, mu1, mu2],
                           proxy_species=_spc("CCCC", "poly"))
        phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"),
                             initial_moments={"poly": (1.0, 5.0, 30.0)},
                             # Species-keyed, exactly as compile_polymer_phase
                             # builds it: oligomer loading + dissolved gas.
                             initial_explicit={p2: explicit_moles, n2: 1.0},
                             pools=[pool], mass_transfer=[])
        reactor = HybridPolymerReactor(
            temperature=(800.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles={n2: 1.0}, polymerPhase=phase,
            terminationTime=(1.0, "s"))
        core = [n2, p2, mu0, mu1, mu2]
        return reactor, core, n2, p2

    def test_run_path_delivers_pool_labeled_initial_explicit(self):
        """RED pre-fix: the solver received the raw {Species: moles} dict
        (silently ignored -- shape mismatch). GREEN: to_solver_object hands the
        solver {pool_label: {dp: moles}} holding ONLY the oligomer loadings."""
        reactor, core, _, _ = self._reactor_with_explicit_pool()
        solver = reactor.to_solver_object(core, [], [], [])
        assert solver.initial_explicit_species == {"poly": {2: 0.3}}

    def test_initial_explicit_moles_reach_y0_and_split_off_tail_moments(self):
        """RED pre-fix: y0[P2] stayed 0 and the tail moments kept the full
        (1.0, 5.0, 30.0) -- the explicit loading never subtracted (moment
        double-count the instant the species is populated). GREEN: y0[P2] =
        0.3 mol and set_initial_conditions step 6 subtracts (N, 2N, 4N)."""
        reactor, core, _, p2 = self._reactor_with_explicit_pool()
        solver = reactor.to_solver_object(core, [], [], [])
        solver.initialize_model(core, [], [], [])
        i_p2 = core.index(p2)
        assert solver.y0[i_p2] == pytest.approx(0.3)
        # Moment moles: V_poly cancels -- mu_tot - mu_exp scaled back by V_poly.
        i_mu0, i_mu1, i_mu2 = 2, 3, 4
        assert solver.y0[i_mu0] == pytest.approx(1.0 - 0.3)
        assert solver.y0[i_mu1] == pytest.approx(5.0 - 2 * 0.3)
        assert solver.y0[i_mu2] == pytest.approx(30.0 - 4 * 0.3)

    def test_dissolved_gas_y0_via_initial_moles_unregressed(self):
        """Regression pin (GREEN before AND after): the dissolved gas keeps
        getting its y0 from initialMoles (step 1); dropping the dissolved dict
        from initial_explicit_species must not zero it or crash the build."""
        reactor, core, n2, _ = self._reactor_with_explicit_pool()
        solver = reactor.to_solver_object(core, [], [], [])
        solver.initialize_model(core, [], [], [])
        assert solver.y0[core.index(n2)] == pytest.approx(1.0)


class TestExplicitDpAutoGenPath:
    """Explicit-DP handshake stage A (producer): the deck flag explicit_dp=True
    auto-generates ONE capped oligomer at DP == cutoff (xs), registers it as a
    real core species, and compile_polymer_phase wires it into the pool's
    explicit_map so PolymerPool.to_config resolves
    explicit_dp_to_species_index == {xs: core index}. These tests drive the
    compiled run path (compile_polymer_phase -> HybridPolymerReactor ->
    to_solver_object -> residual)."""

    @staticmethod
    def _compiled_setup(explicit_dp):
        """Build (reactor, core, poly, dp_spc, mu_indices) through the real
        compile path. Mean tail DP ~5 > xs=3 so the handshake is live;
        k_unzip > 0 supplies the per-chain handshake frequency."""
        from rmgpy.polymer import Polymer
        from rmgpy.quantity import Quantity  # noqa: F401 (parity with compile)
        from rmgpy.rmg.polymer_input import (HybridPolymerReactor,
                                             PolymerPhaseBlueprint,
                                             compile_polymer_phase)

        poly = Polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=520.0, Mw=624.0, initial_mass=0.0052,
                       k_unzip=0.1)
        styrene = _spc("C=Cc1ccccc1", "styrene")
        poly.monomer_product_species = styrene

        dp_spc = poly.generate_explicit_dp_species()
        if explicit_dp:
            poly.explicit_dp = True
            poly.explicit_dp_species = dp_spc

        n2 = _spc("N#N", "N2")
        mu = {}
        for k, smi in ((0, "CO"), (1, "C=O"), (2, "C#N")):
            m = _spc(smi, f"PS_mu{k}")
            m.is_moment_dummy = True
            mu[k] = m

        species_dict = {"PS": poly, "PS_mu0": mu[0], "PS_mu1": mu[1],
                        "PS_mu2": mu[2]}
        blueprint = PolymerPhaseBlueprint(label="melt", species=["PS"],
                                          solvent="PS",
                                          density=(1050.0, "kg/m^3"))
        initial_moles = {n2: 1.0, poly: 0.01}

        phase = compile_polymer_phase(blueprint, initial_moles, species_dict)

        reactor = HybridPolymerReactor(
            temperature=(800.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles=initial_moles, polymerPhase=phase,
            terminationTime=(1.0, "s"))
        core = [n2, poly, mu[0], mu[1], mu[2], styrene, dp_spc]
        mu_indices = (2, 3, 4)
        return reactor, core, poly, dp_spc, mu_indices

    def test_flag_on_config_maps_xs_to_condensed_core_index_zero_t0(self):
        """Flag ON: the compiled solver pool holds exactly {xs: index}, the
        index is condensed-masked (validate_configuration's explicit-DP check
        passes), and with no deck loading the t=0 default is zero explicit
        moles with UNSPLIT tail moments (nothing subtracted for a species
        that starts empty)."""
        reactor, core, _, dp_spc, mu_idx = self._compiled_setup(explicit_dp=True)
        solver = reactor.to_solver_object(core, [], [], [])
        solver.initialize_model(core, [], [], [])

        i_dp = core.index(dp_spc)
        pool_cfg = next(p for p in solver.polymer_pools if p.label == "PS")
        assert pool_cfg.explicit_dp_to_species_index == {3: i_dp}
        assert not solver.gas_species_mask[i_dp]
        assert solver.y0[i_dp] == 0.0
        moments = reactor.polymerPhase.initial_moments["PS"]
        for k in range(3):
            assert solver.y0[mu_idx[k]] == pytest.approx(moments[k])

    def test_flag_off_config_stays_empty(self):
        """Flag OFF (default): explicit_dp_to_species_index == {} --
        byte-identical to today's structurally-inert behavior."""
        reactor, core, _, _, _ = self._compiled_setup(explicit_dp=False)
        solver = reactor.to_solver_object(core, [], [], [])
        pool_cfg = next(p for p in solver.polymer_pools if p.label == "PS")
        assert pool_cfg.explicit_dp_to_species_index == {}

    def test_live_handshake_deposits_into_generated_species_with_matching_drain(self):
        """Live RHS on the auto-gen path: with the flag ON the hybrid
        handshake yields F > 0 deposit into the generated oligomer index and
        drains (mu0, mu1, mu2) by exactly (F, xs*F, xs^2*F) -- pinned by
        differencing the ON residual against an otherwise-identical OFF
        build (deposit and drain are gated together on the map entry)."""
        r_on, core_on, _, dp_on, mu_idx = self._compiled_setup(explicit_dp=True)
        r_off, core_off, _, dp_off, _ = self._compiled_setup(explicit_dp=False)

        s_on = r_on.to_solver_object(core_on, [], [], [])
        s_on.initialize_model(core_on, [], [], [])
        s_off = r_off.to_solver_object(core_off, [], [], [])
        s_off.initialize_model(core_off, [], [], [])

        dn_on = s_on.residual(0.0, s_on.y, np.zeros_like(s_on.y))[0]
        dn_off = s_off.residual(0.0, s_off.y, np.zeros_like(s_off.y))[0]

        i_dp = core_on.index(dp_on)
        assert core_off.index(dp_off) == i_dp  # aligned state vectors

        f_dep = dn_on[i_dp] - dn_off[i_dp]  # F * V_poly [mol/s]
        assert f_dep > 0.0
        xs = 3
        assert dn_on[mu_idx[0]] - dn_off[mu_idx[0]] == pytest.approx(-f_dep)
        assert dn_on[mu_idx[1]] - dn_off[mu_idx[1]] == pytest.approx(-xs * f_dep)
        assert dn_on[mu_idx[2]] - dn_off[mu_idx[2]] == pytest.approx(-xs * xs * f_dep)

    def test_daughter_pools_do_not_get_explicit_dp_in_v1(self):
        """v1 LIMITATION (documented, not silent): daughter pools spawned
        mid-run are derived with explicit_dp_to_species_index == {} even when
        the parent carried explicit_dp=True -- auto-generating and
        core-registering a capped oligomer at enlarge time is out of stage-A
        scope. This test is the documentation."""
        from rmgpy.polymer import Polymer
        from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

        daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                           end_groups=["[CH3]", "[H]"], cutoff=3,
                           Mn=520.0, Mw=624.0, initial_mass=0.001)
        # Simulate a parent with the flag ON whose attributes propagated.
        daughter.explicit_dp = True
        daughter.explicit_dp_species = daughter.generate_explicit_dp_species()

        mu = []
        for k, smi in ((0, "CO"), (1, "C=O"), (2, "C#N")):
            m = _spc(smi, f"PS_d1_mu{k}")
            m.is_moment_dummy = True
            mu.append(m)
        core = [daughter] + mu
        spc_map = {s: i for i, s in enumerate(core)}

        configs = derive_daughter_pool_configs(core, spc_map, {"PS"})

        assert len(configs) == 1
        assert configs[0].label == "PS_d1"
        assert configs[0].explicit_dp_to_species_index == {}


# ---------------------------------------------------------------------------
# Radical-homolysis initiation kernel, Stage 1 (adjudicated adversarial
# round 66). Event: random backbone C-C homolysis on the parent pool at
# R = k(T) * max(mu1 - mu0, 0),  k(T) = A * T^n * exp(-Ea/(R_gas*T))
# with the bond-weighted parent debit
#   B1 = (mu2 - mu1)/(mu1 - mu0),  B2 = (mu3 - mu2)/(mu1 - mu0)
#   parent: dmu0 -= R; dmu1 -= R*B1; dmu2 -= R*B2
# and EACH of the two end-radical daughter pools credited
#   dmu0 += R; dmu1 += R*B1/2; dmu2 += R*(2*B2 - B1)/6.
# mu3 closes via the established log-Lagrange idiom mu0*(mu2/mu1)^3.
# No gas species term (homolysis releases no volatiles); no reverse leg.
# ---------------------------------------------------------------------------

_KHOM_R_GAS = 8.314  # matches the solver's QSSA_R_GAS pin (SI, J/(mol K))


def _khom_triplet(A=1.0e13, n=0.5, Ea=1.2e5):
    return dict(A=A, n=n, Ea=Ea)


def _khom_arrhenius(T, trip=None):
    trip = trip or _khom_triplet()
    return trip["A"] * T ** trip["n"] * math.exp(-trip["Ea"] / (_KHOM_R_GAS * T))


def _khom_core_and_pools(k_homolysis=None, k_scission=0.0,
                         parent_moments=(1.0, 5.0, 30.0),
                         include_daughters=True, extra_species=None):
    """Parent pool 'poly' (mu slots 2-4) + end-radical daughters
    'poly_rad_primary_end' (6-8) and 'poly_rad_secondary_end' (10-12),
    all condensed; index 0 is a gas inert."""
    Inert = _spc("N#N", "N2")
    core = [Inert,
            _spc("CCCC", "poly"),
            _spc("CO", "poly_mu0"), _spc("C=O", "poly_mu1"), _spc("C#N", "poly_mu2")]
    pools = [PolymerPoolConfig(
        label="poly", xs=2, explicit_dp_to_species_index={},
        mu_indices=(2, 3, 4), monomer_poly_index=None,
        k_scission=k_scission, k_unzip=0.0,
        k_homolysis=k_homolysis, tail_kinetics=None)]
    moments = {"poly": tuple(parent_moments)}
    if include_daughters:
        core += [_spc("[CH2]CC", "poly_rad_primary_end"),
                 _spc("CCO", "poly_rad_primary_end_mu0"),
                 _spc("CC=O", "poly_rad_primary_end_mu1"),
                 _spc("CC#N", "poly_rad_primary_end_mu2"),
                 _spc("C[CH]C", "poly_rad_secondary_end"),
                 _spc("CCCO", "poly_rad_secondary_end_mu0"),
                 _spc("CCC=O", "poly_rad_secondary_end_mu1"),
                 _spc("CCC#N", "poly_rad_secondary_end_mu2")]
        pools += [
            PolymerPoolConfig(
                label="poly_rad_primary_end", xs=2,
                explicit_dp_to_species_index={}, mu_indices=(6, 7, 8),
                monomer_poly_index=None, k_scission=0.0, k_unzip=0.0,
                tail_kinetics=None),
            PolymerPoolConfig(
                label="poly_rad_secondary_end", xs=2,
                explicit_dp_to_species_index={}, mu_indices=(10, 11, 12),
                monomer_poly_index=None, k_scission=0.0, k_unzip=0.0,
                tail_kinetics=None)]
        moments["poly_rad_primary_end"] = (0.0, 0.0, 0.0)
        moments["poly_rad_secondary_end"] = (0.0, 0.0, 0.0)
    if extra_species:
        core += extra_species
    mask = np.array([s.label in ("N2",) or (extra_species and s in extra_species)
                     for s in core], dtype=bool)
    return core, mask, pools, moments


def _khom_system(core, mask, pools, moments, T=800.0, rxns=None):
    rs = HybridPolymerSystem(
        T=T, P=1.0e5, initial_mole_fractions={core[0]: 1.0},
        V_poly=1.0, polymer_pools=pools, mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments=moments, termination=[])
    rs.initialize_model(core, rxns or [], [], [])
    return rs


class TestRadicalHomolysisKernel:

    def test_homolysis_moment_derivatives_parent_and_daughters(self):
        """Pins 1+2+3 (round-66 red list): event rate, parent debit, both
        daughter credits, and the three conservation identities -- all
        hand-derived for mu=(1,5,30), V_poly=1:

          mu3 = mu0*(mu2/mu1)^3 = 1*(30/5)^3 = 216
          denom = mu1 - mu0 = 4;  R = k(800)*4
          B1 = (30-5)/4 = 6.25;  B2 = (216-30)/4 = 46.5
          parent:   dmu0 = -R; dmu1 = -6.25*R; dmu2 = -46.5*R
          daughter: dmu0 = +R; dmu1 = +3.125*R;
                    dmu2 = +(2*46.5 - 6.25)/6 * R = +14.4583333...*R
          totals:   dmu0 = +R  (one new chain per event)
                    dmu1 = 0   (monomer units conserved, machine precision)
                    dmu2 = -R*(B1+B2)/3 = -R*52.75/3
                         == k*(mu1 - mu3)/3 = k*(5-216)/3   [legacy
                         Ziff-McGrady random-scission total: R*(B1+B2)/3 =
                         k*(mu1-mu0)*((mu2-mu1)+(mu3-mu2))/((mu1-mu0)*3)
                         = k*(mu3-mu1)/3 -- EXACT algebraic identity]
        """
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=_khom_triplet())
        rs = _khom_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        k = _khom_arrhenius(800.0)
        mu0, mu1, mu2 = 1.0, 5.0, 30.0
        mu3 = mu0 * (mu2 / mu1) ** 3          # 216.0
        denom = mu1 - mu0                     # 4.0
        R = k * denom
        B1 = (mu2 - mu1) / denom              # 6.25
        B2 = (mu3 - mu2) / denom              # 46.5

        # parent debit (pin 2)
        assert np.isclose(dn_dt[2], -R, rtol=1e-12)
        assert np.isclose(dn_dt[3], -R * B1, rtol=1e-12)
        assert np.isclose(dn_dt[4], -R * B2, rtol=1e-12)
        # each daughter credit (pin 2)
        for base in (6, 10):
            assert np.isclose(dn_dt[base], R, rtol=1e-12)
            assert np.isclose(dn_dt[base + 1], R * B1 / 2.0, rtol=1e-12)
            assert np.isclose(dn_dt[base + 2], R * (2.0 * B2 - B1) / 6.0,
                              rtol=1e-12)
        # identities (pin 3)
        total_mu1 = dn_dt[3] + dn_dt[7] + dn_dt[11]
        assert total_mu1 == 0.0, (
            f"total dmu1 must vanish at machine precision, got {total_mu1!r}")
        total_mu0 = dn_dt[2] + dn_dt[6] + dn_dt[10]
        assert np.isclose(total_mu0, R, rtol=1e-12)
        total_mu2 = dn_dt[4] + dn_dt[8] + dn_dt[12]
        legacy_mu2 = k * (mu1 - mu3) / 3.0    # Ziff-McGrady total
        assert np.isclose(total_mu2, legacy_mu2, rtol=1e-12)
        assert np.isclose(total_mu2, -R * (B1 + B2) / 3.0, rtol=1e-12)
        # no gas species term (round 66: homolysis produces no volatiles)
        assert dn_dt[0] == 0.0
        # proxies carry no kernel flux
        assert dn_dt[1] == 0.0 and dn_dt[5] == 0.0 and dn_dt[9] == 0.0

    def test_homolysis_zero_flux_at_mu1_le_mu0_and_empty_pool(self):
        """Pin 1: R = k(T)*max(mu1-mu0, 0) -- zero flux at mu1 <= mu0 and on
        an empty pool; the whole residual is exactly zero (no reactions, no
        other channels)."""
        for parent_moments in ((2.0, 2.0, 2.0), (0.0, 0.0, 0.0)):
            core, mask, pools, moments = _khom_core_and_pools(
                k_homolysis=_khom_triplet(), parent_moments=parent_moments)
            rs = _khom_system(core, mask, pools, moments)
            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            assert np.all(dn_dt == 0.0), (
                f"homolysis must carry zero flux for moments "
                f"{parent_moments}, got nonzero {dn_dt.nonzero()}")

    def test_homolysis_arrhenius_evaluated_at_runtime_temperature(self):
        """Pin 7: k(T) is evaluated at the solver's runtime T (round 66: 'a
        scalar precomputed at 1100 K will poison any ramp/TA replay'). k at
        two different temperatures matches A*T^n*exp(-Ea/(R*T)) -- the n=0.5
        triplet also catches a dropped T^n factor."""
        rates = {}
        for T in (700.0, 1100.0):
            core, mask, pools, moments = _khom_core_and_pools(
                k_homolysis=_khom_triplet())
            rs = _khom_system(core, mask, pools, moments, T=T)
            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            rates[T] = dn_dt[6]     # primary-daughter mu0 credit == R
            expected = _khom_arrhenius(T) * (5.0 - 1.0)
            assert np.isclose(rates[T], expected, rtol=1e-10), (
                f"k({T}) mismatch: got {rates[T]}, expected {expected}")
        assert rates[1100.0] > rates[700.0] * 10.0  # genuinely T-dependent

    def test_homolysis_requires_resolvable_end_radical_pools(self):
        """Loud guard: k_homolysis enabled but the end-radical daughter pools
        absent from the solver's pool configs is a hard initialize error
        (the producer creates them at model setup; a missing pool means the
        kernel's credits would silently vanish)."""
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=_khom_triplet(), include_daughters=False)
        with pytest.raises(ValueError,
                           match=r"poly.*k_homolysis.*_rad_primary_end"):
            _khom_system(core, mask, pools, moments)

    def test_solver_rejects_k_homolysis_with_positive_k_scission(self):
        """Pin 5 (solver last line of defense): k_scission > 0 AND
        k_homolysis on one pool double-counts random backbone homolysis."""
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=_khom_triplet(), k_scission=0.3)
        with pytest.raises(ValueError,
                           match=r"Pool poly.*k_homolysis.*k_scission.*mutually exclusive"):
            _khom_system(core, mask, pools, moments)

    def test_solver_rejects_k_homolysis_with_qssa_channel(self):
        """Pin 5 (solver last line of defense): QSSA random initiation AND
        k_homolysis on one pool double-counts initiation."""
        qssa = dict(initiation=dict(A=1.0e15, n=0.0, Ea=3.0e5),
                    depropagation=dict(A=1.0e13, n=0.0, Ea=8.0e4),
                    termination=dict(A=1.0e8, n=0.0, Ea=1.0e4))
        Mono = _spc("C=CC", "propene_gas")
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=None, extra_species=[Mono])
        parent = dataclasses.replace(
            pools[0], k_homolysis=_khom_triplet(),
            radical_qssa_unzip=qssa, monomer_poly_index=len(core) - 1)
        pools[0] = parent
        with pytest.raises(ValueError,
                           match=r"Pool poly.*k_homolysis.*radical_qssa_unzip.*mutually exclusive"):
            _khom_system(core, mask, pools, moments)

    def test_solver_rejects_k_homolysis_with_positive_k_unzip(self):
        """P1 (adjudicated round 67, solver last line of defense): legacy
        k_unzip (phenomenological closed-chain monomer-loss channel) AND
        k_homolysis (radical-end pools feeding explicit beta-scission/unzip
        chemistry) on one pool double-carries depolymerization.
        monomer_poly_index is wired so the k_unzip routing guard cannot mask
        this exclusion."""
        Mono = _spc("C=CC", "propene_gas")
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=None, extra_species=[Mono])
        parent = dataclasses.replace(
            pools[0], k_homolysis=_khom_triplet(),
            k_unzip=0.4, monomer_poly_index=len(core) - 1)
        pools[0] = parent
        with pytest.raises(ValueError,
                           match=r"Pool poly.*k_homolysis.*k_unzip.*mutually exclusive"):
            _khom_system(core, mask, pools, moments)

    def test_homolysis_guard_rejects_negative_daughter_mu2_credit(self, caplog):
        """P1 (adjudicated round 67): B1 >= 0 and B2 >= 0 do NOT imply
        2*B2 - B1 >= 0. mu = (1, 4, 8.2) -> mu3 = mu0*(mu2/mu1)^3 = 8.615,
        so B1 = 4.2/3 = 1.4 >= 0 and B2 = 0.415/3 ~ 0.138 >= 0 pass the
        per-moment guards, but the daughter mu2 credit factor
        2*B2 - B1 = (2*mu3 - 3*mu2 + mu1)/(mu1 - mu0) = -3.37/3 < 0 would
        credit NEGATIVE mu2 to both radical daughters. Convention (mirrors
        the existing nonrealizable-state guard in this kernel): loud
        warn-once + the kernel contributes ZERO flux."""
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=_khom_triplet(), parent_moments=(1.0, 4.0, 8.2))
        rs = _khom_system(core, mask, pools, moments)
        with caplog.at_level(logging.WARNING):
            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(dn_dt == 0.0), (
            f"kernel must contribute zero flux on a moment state whose "
            f"daughter mu2 credit factor 2*B2 - B1 is negative; nonzero at "
            f"{dn_dt.nonzero()}: {dn_dt[dn_dt != 0.0]}")
        assert any("k_homolysis kernel skipped" in r.getMessage()
                   for r in caplog.records), (
            "the out-of-domain guard must warn loudly, not skip silently")

    def test_homolysis_credits_use_stable_direct_forms_near_exhaustion(self):
        """P2 (adjudicated round 67): near DP -> 1 exhaustion (mu1 barely
        above mu0) the divided forms R*Bi = k*B*((mu_{i+1}-mu_i)/B) suffer
        ratio cancellation; the credits must be computed in the
        algebraically identical DIRECT forms:
            R*B1          = k*(mu2 - mu1)
            R*B2          = k*(mu3 - mu2)
            R*(2B2-B1)/6  = k*(2*mu3 - 3*mu2 + mu1)/6
        Pinned by EXACT (bitwise) equality against the direct forms at a
        near-exhaustion state; the divide-then-remultiply forms differ in
        the last ulps here."""
        mu0, mu1, mu2 = 1.0, 1.0 + 1e-9, 1.0 + 3e-9
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=_khom_triplet(), parent_moments=(mu0, mu1, mu2))
        rs = _khom_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        k = _khom_arrhenius(800.0)
        # bitwise replica of the solver's log-Lagrange mu3 closure
        # (_safe_mu3_from_mu012) -- exact equality below depends on it
        mu3 = float(np.exp(3.0 * np.log(mu2) - 3.0 * np.log(mu1)
                           + np.log(mu0)))
        # parent debit, stable direct forms
        assert dn_dt[2] == -(k * (mu1 - mu0))
        assert dn_dt[3] == -(k * (mu2 - mu1))
        assert dn_dt[4] == -(k * (mu3 - mu2))
        # EACH daughter credit, stable direct forms
        for base in (6, 10):
            assert dn_dt[base] == k * (mu1 - mu0)
            assert dn_dt[base + 1] == k * (mu2 - mu1) / 2.0
            assert dn_dt[base + 2] == k * (2.0 * mu3 - 3.0 * mu2 + mu1) / 6.0
        # analytic invariants survive the exhaustion regime
        assert dn_dt[3] + dn_dt[7] + dn_dt[11] == 0.0
        total_mu2 = dn_dt[4] + dn_dt[8] + dn_dt[12]
        assert np.isclose(total_mu2, k * (mu1 - mu3) / 3.0, rtol=1e-12)

    def test_solver_rejects_malformed_k_homolysis_triplet(self):
        """Pin 5 (solver field validation): a directly-constructed config
        with a malformed triplet must fail loudly, not run a garbage rate."""
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=dict(A=-1.0, n=0.0, Ea=8.0e4))
        with pytest.raises(ValueError, match=r"Pool poly.*k_homolysis.*A.*> 0"):
            _khom_system(core, mask, pools, moments)

    def test_refused_homolysis_rows_zero_flux_and_census_pairs_kernel(self):
        """Pin 6: the run-5 gas-radical association shape (gas radicals ->
        condensed proxy), refused conduit-deferred, stays refused and carries
        EXACTLY zero flux while the k_homolysis kernel is live -- and the
        census pairs them: a warn-level entry names the refused row(s) as
        superseded by the kernel."""
        R1 = _spc("[CH2]C(C)C", "isobutyl")
        R2 = _spc("C[CH]CCC", "pentyl2")
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=_khom_triplet(), extra_species=[R1, R2])
        proxy = core[1]
        refused = Reaction(reactants=[R1, R2], products=[proxy],
                           kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0,
                                              Ea=(0.0, "J/mol"), T0=(1.0, "K")),
                           reversible=False)
        # run-5 stamp (the classifier itself is pinned in polymerTest.py::
        # test_gas_radical_association_into_condensed_proxy_refused_...)
        refused.polymer_refused = True
        refused.polymer_refused_accumulating = False    # conduit-deferred

        rs_with = _khom_system(core, mask, pools, moments, rxns=[refused])
        dn_with = rs_with.residual(0.0, rs_with.y, np.zeros_like(rs_with.y))[0]

        core2, mask2, pools2, moments2 = _khom_core_and_pools(
            k_homolysis=_khom_triplet(),
            extra_species=[_spc("[CH2]C(C)C", "isobutyl"),
                           _spc("C[CH]CCC", "pentyl2")])
        rs_wo = _khom_system(core2, mask2, pools2, moments2, rxns=[])
        dn_wo = rs_wo.residual(0.0, rs_wo.y, np.zeros_like(rs_wo.y))[0]

        assert np.any(dn_wo != 0.0)                 # kernel flux is live
        assert np.array_equal(dn_with, dn_wo), (
            "refused row must add exactly zero flux next to the live kernel")
        # still censused refused (not silently un-refused by the kernel)
        assert any(c["reason"] == "conduit-deferred"
                   for c in rs_with.refused_reaction_census)
        # the pairing census names the refused rows as superseded
        census = rs_with.homolysis_supersession_census
        assert census, "kernel-live census must pair the refused rows"
        entry = next(c for c in census if c["pool"] == "poly")
        assert len(entry["superseded_rows"]) == 1

    def test_homolysis_diagnostic_reports_kernel_and_radical_pools(self, capsys):
        """Round 66: radical inventories must be traceable in run logs -- the
        POLYMER SOLVER DIAGNOSTIC block lists the kernel triplet on the parent
        pool and enumerates both end-radical pools with their moment slots."""
        core, mask, pools, moments = _khom_core_and_pools(
            k_homolysis=_khom_triplet())
        _khom_system(core, mask, pools, moments)
        out = capsys.readouterr().out
        assert "POLYMER SOLVER DIAGNOSTIC" in out
        assert "k_homolysis" in out
        assert "poly_rad_primary_end" in out
        assert "poly_rad_secondary_end" in out


# ---------------------------------------------------------------------------
# End-radical DEPROPAGATION kernel (adjudicated round 74 SS2, the run-6 wall
# fix): pool-level unzip consumption channel for radical-end daughter pools.
# For a radical-end pool with one active radical end per chain and
# k_dep(T) = A*T^n*exp(-Ea/(R_gas*T)):
#   event rate  R = k_dep*mu0        (one monomer released per event)
#   gas credit  +R of the pool's monomer volatile (monomer_poly_index)
#   dmu1 = -R;  dmu2 = -k_dep*(2*mu1 - mu0);  dmu0 = -k_dep*N1
# with N1 = mu0 * P(DP=1) from the EXISTING discrete/gamma moment closure
# (_gamma_params_from_mu012 + _gamma_prob_conditional_hybrid(1, 0, k, th)),
# floored by a smooth terminal boundary law (smoothstep over mean DP in
# [1, 2]) so the pool never stalls at a one-repeat-per-chain residue.
# SMOOTH exhaustion gate (r74 SS5: no hard max(...,0) cliff): all release
# terms carry g = 1 - sp(1 - mean), sp(x) = x^3/(x^2 + W^2), which is
# EXACTLY 1 for mean >= 1 and rolls off C2-smoothly below.
# Mass invariant: gas release rate == -dmu1 identically (bitwise).
# ---------------------------------------------------------------------------


def _kdep_triplet(A=1.0e13, n=0.5, Ea=1.2e5):
    return dict(A=A, n=n, Ea=Ea)


def _kdep_arrhenius(T, trip=None):
    trip = trip or _kdep_triplet()
    return trip["A"] * T ** trip["n"] * math.exp(-trip["Ea"] / (_KHOM_R_GAS * T))


def _kdep_expected_p1(mu0, mu1, mu2):
    """Independent replica of the adjudicated N1 closure: DP=1 chain fraction
    from the gamma closure over (mu0, mu1, mu2), discretized on half-integer
    bins and conditioned on n >= 1, floored by the smoothstep terminal
    boundary law over mean DP in [1, 2]."""
    from scipy.special import gammainc

    mean = mu1 / mu0
    t = mean - 1.0
    if t <= 0.0:
        p_fb = 1.0
    elif t >= 1.0:
        p_fb = 0.0
    else:
        p_fb = 1.0 - (3.0 * t * t - 2.0 * t * t * t)
    p_gamma = 0.0
    pdi = mu2 * mu0 / (mu1 * mu1)
    if math.isfinite(pdi) and pdi > 1.0 + 1e-6:
        k = 1.0 / (pdi - 1.0)
        theta = mean / k
        if k > 0.0 and theta > 0.0:
            F = lambda x: float(gammainc(k, x / theta))
            tail = 1.0 - F(0.5)
            if tail > 1e-12:
                p_gamma = max(0.0, F(1.5) - F(0.5)) / tail
    return min(1.0, max(0.0, max(p_gamma, p_fb)))


def _kdep_core_and_pools(k_depropagation=None, moments=(1.0, 5.0, 30.0),
                         k_unzip=0.0, radical_qssa_unzip=None,
                         k_homolysis=None, monomer=True):
    """One radical-end pool 'poly_rad_primary_end' (mu slots 3-5, condensed)
    plus a gas inert (0) and the released-monomer gas volatile (1)."""
    Inert = _spc("N#N", "N2")
    Mono = _spc("C=CC", "propene_gas")
    core = [Inert, Mono,
            _spc("[CH2]CC", "poly_rad_primary_end"),
            _spc("CCO", "poly_rad_primary_end_mu0"),
            _spc("CC=O", "poly_rad_primary_end_mu1"),
            _spc("CC#N", "poly_rad_primary_end_mu2")]
    pools = [PolymerPoolConfig(
        label="poly_rad_primary_end", xs=2,
        explicit_dp_to_species_index={}, mu_indices=(3, 4, 5),
        monomer_poly_index=(1 if monomer else None),
        monomer_mw_g_mol=42.08,
        k_scission=0.0, k_unzip=k_unzip,
        radical_qssa_unzip=radical_qssa_unzip,
        k_homolysis=k_homolysis,
        k_depropagation=k_depropagation, tail_kinetics=None)]
    all_moments = {"poly_rad_primary_end": tuple(moments)}
    mask = np.array([s.label in ("N2", "propene_gas") for s in core],
                    dtype=bool)
    return core, mask, pools, all_moments


class TestEndRadicalDepropagationKernel:

    GAS, MU0, MU1, MU2 = 1, 3, 4, 5

    def test_deprop_law_moment_derivatives_and_gas_rate(self):
        """r74 SS2 law at a healthy state mu=(1,5,30), V_poly=1, T=800:
            R    = k(800)*mu0            (mean DP 5 >= 1 -> gate g == 1)
            gas  = +R at monomer_poly_index
            dmu1 = -R
            dmu2 = -k*(2*mu1 - mu0) = -9k  (up to the documented smooth
                   O(W^2/(mean-1)) deficit of the C2 positive-part)
            dmu0 = -k*mu0*p1, p1 the gamma-closure DP=1 fraction
        (pdi = 30/25 = 1.2 -> k_shape = 5, theta = 1)."""
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet())
        rs = _khom_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        k = _kdep_arrhenius(800.0)
        mu0, mu1, mu2 = 1.0, 5.0, 30.0
        R = k * mu0
        assert np.isclose(dn_dt[self.GAS], R, rtol=1e-9)
        assert np.isclose(dn_dt[self.MU1], -R, rtol=1e-9)
        assert np.isclose(dn_dt[self.MU2], -k * (2.0 * mu1 - mu0), rtol=1e-4)
        p1 = _kdep_expected_p1(mu0, mu1, mu2)
        assert 0.0 < p1 < 1.0
        assert np.isclose(dn_dt[self.MU0], -k * mu0 * p1, rtol=1e-9)
        # proxy and inert carry no kernel flux
        assert dn_dt[0] == 0.0 and dn_dt[2] == 0.0

    def test_deprop_exact_mass_invariant(self):
        """r74 SS2 mass pin: d(condensed mass) + d(propene mass) = 0 EXACTLY.
        The gas release rate and the mu1 drain are the SAME float (bitwise),
        so monomer_MW * R closes the balance at machine precision."""
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet())
        rs = _khom_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert dn_dt[self.GAS] > 0.0
        assert dn_dt[self.GAS] == -dn_dt[self.MU1], (
            f"gas release {dn_dt[self.GAS]!r} must equal the mu1 drain "
            f"{-dn_dt[self.MU1]!r} bitwise -- anything else un-conserves "
            f"mass under MW multiplication")
        mw = pools[0].monomer_mw_g_mol
        d_condensed = pools[0].condensed_mass_g(dn_dt[self.MU0],
                                                dn_dt[self.MU1])
        d_gas = mw * dn_dt[self.GAS]
        assert d_condensed + d_gas == 0.0

    def test_deprop_n1_boundary_all_dp1_residue_terminates(self):
        """N1 boundary pin (r74: 'NO permanent dmu0=0 -- that stalls at a
        one-repeat-per-chain residue'): an all-DP-1 pool (c, c, c) has
        N1 = mu0 exactly (fallback p1 = 1, gate g = 1), so all three moments
        decay together at -k*c and the chains terminate."""
        c = 0.7
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet(), moments=(c, c, c))
        rs = _khom_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        k = _kdep_arrhenius(800.0)
        assert np.isclose(dn_dt[self.MU0], -k * c, rtol=1e-9)
        assert np.isclose(dn_dt[self.MU1], -k * c, rtol=1e-9)
        assert np.isclose(dn_dt[self.MU2], -k * c, rtol=1e-9)
        assert dn_dt[self.GAS] == -dn_dt[self.MU1]

    def test_deprop_pool_drains_to_zero_no_dp1_stall(self):
        """Live-path exhaustion pin: integrating the kernel-only pool from
        mu = (1, 3, 10) must drain the CHAIN COUNT and the REPEAT UNITS to
        ~zero (not stall at the mu0 = 1 one-repeat-per-chain residue a
        permanent dmu0 = 0 law leaves), stay finite throughout, and deliver
        the released units to the gas volatile within integration error.
        (mu2 may keep a small inert orphan: the r74 law dmu2 =
        -k*(2*mu1 - mu0) is exact per event, but the CLOSURE-estimated N1
        makes the coupled trajectory drift from a realizable distribution;
        the residue is variance bookkeeping, carries no mass and no flux,
        and is pinned bounded here.)"""
        trip = dict(A=1.0, n=0.0, Ea=0.0)   # k_dep == 1 s^-1 at any T
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=trip, moments=(1.0, 3.0, 10.0))
        rs = _khom_system(core, mask, pools, moments, T=800.0)

        y = rs.y.copy().astype(float)

        def f(yv):
            return rs.residual(0.0, yv, np.zeros_like(yv))[0]

        dt = 0.02
        for _ in range(int(60.0 / dt)):     # 60 e-folds of k_dep
            k1 = f(y)
            k2 = f(y + 0.5 * dt * k1)
            k3 = f(y + 0.5 * dt * k2)
            k4 = f(y + dt * k3)
            y = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            assert np.all(np.isfinite(y)), "state must stay finite"
        assert y[self.MU0] < 1e-3, (
            f"chain count stalled at mu0={y[self.MU0]:g} -- the DP=1 "
            f"residue never terminated")
        assert y[self.MU1] < 1e-3, (
            f"repeat units stalled at mu1={y[self.MU1]:g}")
        # documented closure-orphan bound (inert, massless, fluxless)
        assert 0.0 <= y[self.MU2] < 0.5
        # released units land in the gas volatile (mass closure over the
        # run: every drained mu1 unit is a gas monomer, nothing else moves
        # mass)
        assert np.isclose(y[self.GAS] + y[self.MU1], 3.0, atol=5e-3)

    def test_deprop_smooth_exhaustion_gate_no_cliff(self):
        """r74 SS5 (designed-in from the start): the kernel rolls off
        SMOOTHLY as the pool empties -- no max(...,0) cliff at mean DP = 1.
        Pins: (a) the gate is EXACTLY 1 at mean >= 1 (healthy law exact);
        (b) crossing mean = 1 downward changes the release rate by O(eps^3)
        (C2 rolloff), not a jump; (c) at the pathological mu1 = 0, mu0 > 0
        noise state the gas release is ~W^2-suppressed (no fountain) while
        the chain count still drains at the full -k*mu0 (no DASPK grind)."""
        k = _kdep_arrhenius(800.0)

        def rate_at(mu):
            core, mask, pools, moments = _kdep_core_and_pools(
                k_depropagation=_kdep_triplet(), moments=mu)
            rs = _khom_system(core, mask, pools, moments, T=800.0)
            return rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        # (a) exact law at mean exactly 1 from above
        dn = rate_at((1.0, 1.0, 1.0))
        assert dn[self.GAS] == k * 1.0

        # (b) no cliff crossing mean = 1: relative change is O(eps^3/W^2)
        eps = 1e-3
        dn_below = rate_at((1.0, 1.0 - eps, 1.0))
        rel_jump = abs(dn_below[self.GAS] - dn[self.GAS]) / dn[self.GAS]
        assert rel_jump < 1e-4, (
            f"release rate jumped by {rel_jump:g} across the exhaustion "
            f"boundary -- the gate must roll off smoothly (C2), not cliff")

        # (c) mu1 = 0 noise state: negligible fabrication, full chain drain
        dn0 = rate_at((0.5, 0.0, 0.0))
        assert 0.0 <= dn0[self.GAS] < k * 0.5 * 1e-3, (
            f"gas fabrication {dn0[self.GAS]:g} from an empty-unit pool")
        assert np.isclose(dn0[self.MU0], -k * 0.5, rtol=1e-9), (
            "chain count must keep draining at -k*mu0 (fallback p1 = 1) so "
            "the noise state cannot become a stiff no-outlet grind")

    def test_deprop_arrhenius_evaluated_at_runtime_temperature(self):
        """k_dep(T) is evaluated at the solver's runtime T (round 66
        precedent: 'a scalar precomputed at 1100 K will poison any ramp/TA
        replay'); the n=0.5 triplet also catches a dropped T^n factor."""
        rates = {}
        for T in (700.0, 1100.0):
            core, mask, pools, moments = _kdep_core_and_pools(
                k_depropagation=_kdep_triplet())
            rs = _khom_system(core, mask, pools, moments, T=T)
            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            rates[T] = dn_dt[self.GAS]
            assert np.isclose(rates[T], _kdep_arrhenius(T) * 1.0, rtol=1e-9)
        assert rates[1100.0] > rates[700.0] * 10.0

    def test_deprop_zero_flux_on_empty_pool(self):
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet(), moments=(0.0, 0.0, 0.0))
        rs = _khom_system(core, mask, pools, moments)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(dn_dt == 0.0)

    def test_solver_rejects_k_depropagation_without_monomer_index(self):
        """LAST line of defense: the kernel releases one monomer per event;
        without a resolvable gas destination the released units leave the
        condensed phase silently un-conserved."""
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet(), monomer=False)
        with pytest.raises(ValueError,
                           match=r"poly_rad_primary_end.*k_depropagation.*monomer"):
            _khom_system(core, mask, pools, moments)

    def test_solver_rejects_k_depropagation_with_positive_k_unzip(self):
        """Mutual exclusion (probe finding: k_unzip is the legacy
        phenomenological form of the SAME chain-end monomer-release event --
        dmu1 = -k*mu0, dmu2 = -k*(2mu1-mu0), permanent dmu0 = 0): enabling
        both double-carries depropagation."""
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet(), k_unzip=0.4)
        with pytest.raises(ValueError,
                           match=r"poly_rad_primary_end.*k_depropagation.*k_unzip.*mutually exclusive"):
            _khom_system(core, mask, pools, moments)

    def test_solver_rejects_k_depropagation_with_qssa_channel(self):
        """Mutual exclusion: the QSSA channel's depropagation block IS
        lumped chain-end depropagation on the same pool."""
        qssa = dict(initiation=dict(A=1.0e15, n=0.0, Ea=3.0e5),
                    depropagation=dict(A=1.0e13, n=0.0, Ea=8.0e4),
                    termination=dict(A=1.0e8, n=0.0, Ea=1.0e4))
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet(), radical_qssa_unzip=qssa)
        with pytest.raises(ValueError,
                           match=r"poly_rad_primary_end.*k_depropagation.*radical_qssa_unzip.*mutually exclusive"):
            _khom_system(core, mask, pools, moments)

    def test_solver_rejects_k_depropagation_with_k_homolysis(self):
        """Mutual exclusion (r74 SS3: multi-generation homolysis of
        radical-ended chains is DEFERRED -- 'do not inherit the block
        blindly'): a pool cannot carry both the closed-chain initiation
        kernel and the radical-end depropagation kernel."""
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=_kdep_triplet(), k_homolysis=_khom_triplet())
        with pytest.raises(ValueError,
                           match=r"poly_rad_primary_end.*k_depropagation.*k_homolysis.*mutually exclusive"):
            _khom_system(core, mask, pools, moments)

    def test_solver_rejects_malformed_k_depropagation_triplet(self):
        core, mask, pools, moments = _kdep_core_and_pools(
            k_depropagation=dict(A=-1.0, n=0.0, Ea=8.0e4))
        with pytest.raises(ValueError,
                           match=r"poly_rad_primary_end.*k_depropagation.*A.*> 0"):
            _khom_system(core, mask, pools, moments)


class TestDepropagationDaughterWiring:
    """Deck param -> spawned end-radical daughter pools (probe finding c):
    k_depropagation is declared on the PARENT pool's k_homolysis context;
    the producer copies it (plus the released-monomer routing, held BY
    REFERENCE) onto both spawned daughters, and
    derive_daughter_pool_configs wires it into the daughters' solver
    configs."""

    @staticmethod
    def _parent_with_kdep():
        from rmgpy.polymer import Polymer
        pp = Polymer(label='PP', monomer='[CH2][CH](C)',
                     end_groups=['[H]', '[H]'], cutoff=3,
                     Mn=1500.0, Mw=1800.0, initial_mass=0.1,
                     k_homolysis=_khom_triplet(),
                     k_depropagation=_kdep_triplet())
        pp.monomer_product_species = _spc("C=CC", "propene")
        return pp

    def test_daughters_carry_k_depropagation_and_monomer_routing(self):
        pp = self._parent_with_kdep()
        prim, sec = pp.generate_end_radical_daughters()
        for d in (prim, sec):
            assert d.k_depropagation == _kdep_triplet()
            # deep copy, never aliasing the parent's dict
            assert d.k_depropagation is not pp.k_depropagation
            # routing held BY REFERENCE (identity is load-bearing for the
            # object-keyed spc_map resolution downstream)
            assert d.monomer_product_species is pp.monomer_product_species

    def test_derive_daughter_pool_configs_wires_k_depropagation(self):
        from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
        pp = self._parent_with_kdep()
        prim, sec = pp.generate_end_radical_daughters()
        mono = pp.monomer_product_species
        core = [mono]
        for d in (prim, sec):
            core += [d, _spc("CCO", f"{d.label}_mu0"),
                     _spc("CC=O", f"{d.label}_mu1"),
                     _spc("CC#N", f"{d.label}_mu2")]
        spc_map = {s: i for i, s in enumerate(core)}
        cfgs = derive_daughter_pool_configs(core, spc_map, {"PP"})
        assert len(cfgs) == 2
        for cfg in cfgs:
            assert cfg.k_depropagation == _kdep_triplet()
            assert cfg.monomer_poly_index == spc_map[mono]

    def test_derive_daughter_pool_configs_rejects_kdep_without_routing(self):
        """A daughter carrying the kernel but whose released-monomer species
        is missing from the core map must fail LOUDLY at config derivation
        (mass would silently un-conserve), never build a mute config."""
        from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
        pp = self._parent_with_kdep()
        prim, _ = pp.generate_end_radical_daughters()
        prim.monomer_product_species = None
        core = [prim, _spc("CCO", f"{prim.label}_mu0"),
                _spc("CC=O", f"{prim.label}_mu1"),
                _spc("CC#N", f"{prim.label}_mu2")]
        spc_map = {s: i for i, s in enumerate(core)}
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*k_depropagation.*monomer"):
            derive_daughter_pool_configs(core, spc_map, {"PP"})


# ---------------------------------------------------------------------------
# Side-group homolysis initiation kernel (FR1-K1) -- KERNEL-V2 (explicit Br-
# inventory depletion, TGA-faithful multi-loss): pool-level side-group X-loss
# (e.g. aliphatic C-Br -> Br.(gas) + mid-chain backbone radical, chain length
# UNCHANGED). v2 replaces the v1 one-loss feature-pool encoding with an
# explicit auxiliary extensive state Z per (pool, channel row): Z = moles of
# remaining removable X (Br) sites, seeded Z(0) = sites_per_unit*mu1(0) and
# depleting as the X atoms leave the melt. State vector order is
# [core species | U slots (QSSA) | Z slots]. Per channel with sites = s and
# k(T) = A*T^n*exp(-Ea/(R_gas*T)) per site, at runtime T and V_poly:
#   k(T)          = A*T^n*exp(-Ea/(R_gas*T))     [s^-1 per site]
#   R             = k(T) * Z / V_poly            [mol/(m^3 s), NO extra sites
#                                                factor -- sites folded into Z0]
#   inventory     dZ/dt = -R*V_poly              [mol/s] (Br drains from melt)
#   gas credit    gas X radical += R             (same small_src -> dn*V_poly)
#   parent moments mu0/mu1/mu2 are UNCHANGED (no debit, no feature-pool credit)
#   mass contract condensed carrier mass = mu1*MW - Sum_row max(0, s*mu1 - Z)*M_X
#                 so d(condensed X mass)/dt + d(gas X mass)/dt = 0 exactly.
# v2 DELETES the SGH feature/daughter pool entirely (carrier chain_mass_defect
# is 0). An SGH carrier that ALSO moves mu1 (k_homolysis / k_depropagation /
# tail_kinetics / k_unzip / radical_qssa_unzip / live resolved row flux) is
# REFUSED at setup; k_scission (verified mu1-preserving) is allowed.
# ---------------------------------------------------------------------------


def _sgh_channel(label="aliphatic_C-Br", A=1.0e13, n=0.5, Ea=1.2e5,
                 site_selector="aliphatic", sites_per_unit=2.0,
                 gas_product="[Br]"):
    return dict(label=label, A=A, n=n, Ea=Ea, site_selector=site_selector,
                sites_per_unit=sites_per_unit, gas_product=gas_product)


def _sgh_arrhenius(T, ch=None):
    ch = ch or _sgh_channel()
    return ch["A"] * T ** ch["n"] * math.exp(-ch["Ea"] / (_KHOM_R_GAS * T))


def _sgh_mx_g_mol(gas_product="[Br]"):
    return Molecule().from_smiles(gas_product).get_molecular_weight() * 1000.0


_SGH_PARENT_MW = 106.949  # ~vinyl-bromide repeat unit [g/mol]; exact value shared


def _sgh_core_and_pools(channels=None, parent_moments=(1.0, 5.0, 30.0),
                        gas_indices="auto",
                        k_homolysis=None, k_unzip=0.0, monomer_poly_index=None,
                        extra_species=None):
    """SGH kernel-v2 fixture: a single condensed carrier pool 'poly' (mu slots
    3-5); index 0 is a gas inert, index 1 the gas Br radical, so n_core = 6.

    v2 DELETES the side-group feature/daughter pool entirely -- Br leaves the
    melt through the per-channel auxiliary Z inventory (Z = moles of removable
    X sites, Z(0) = sites*mu1(0), dZ/dt = -k*Z/V_poly), NOT a moment-debited
    daughter pool. The state vector is [core | U slots | Z slots]; with no QSSA
    U slots in these fixtures, channel row r's Z slot is the absolute ODE index
    n_core + r (= 6 + r). Read the slot from rs.sgh_z_slot[r] rather than
    hard-coding it."""
    channels = channels if channels is not None else [_sgh_channel()]
    Inert = _spc("N#N", "N2")
    BrRad = _spc("[Br]", "Br_rad")
    core = [Inert, BrRad,
            _spc("CCCC", "poly"),
            _spc("CO", "poly_mu0"), _spc("C=O", "poly_mu1"),
            _spc("C#N", "poly_mu2")]
    moments = {"poly": tuple(parent_moments)}
    pools = []
    if gas_indices == "auto":
        gas_indices = tuple(1 for _ in channels)
    parent = PolymerPoolConfig(
        label="poly", xs=2, explicit_dp_to_species_index={},
        mu_indices=(3, 4, 5), monomer_poly_index=monomer_poly_index,
        k_scission=0.0, k_unzip=k_unzip,
        monomer_mw_g_mol=_SGH_PARENT_MW,
        k_homolysis=k_homolysis,
        side_group_homolysis=channels,
        side_group_gas_indices=gas_indices,
        tail_kinetics=None)
    pools.insert(0, parent)
    if k_homolysis is not None:
        # end-radical daughters for coexistence tests
        for suffix, smis in (("_rad_primary_end", ("[CH2]CC", "OCCO", "OCC=O", "OCC#N")),
                             ("_rad_secondary_end", ("C[CH]C", "OCCCO", "OCCC=O", "OCCC#N"))):
            base = len(core)
            core += [_spc(smis[0], f"poly{suffix}"),
                     _spc(smis[1], f"poly{suffix}_mu0"),
                     _spc(smis[2], f"poly{suffix}_mu1"),
                     _spc(smis[3], f"poly{suffix}_mu2")]
            pools.append(PolymerPoolConfig(
                label=f"poly{suffix}", xs=2, explicit_dp_to_species_index={},
                mu_indices=(base + 1, base + 2, base + 3),
                monomer_poly_index=None, k_scission=0.0, k_unzip=0.0,
                monomer_mw_g_mol=_SGH_PARENT_MW, tail_kinetics=None))
            moments[f"poly{suffix}"] = (0.0, 0.0, 0.0)
    if extra_species:
        core += extra_species
    mask = np.array([s.label in ("N2", "Br_rad")
                     or (extra_species and s in extra_species)
                     for s in core], dtype=bool)
    return core, mask, pools, moments


def _sgh_system(core, mask, pools, moments, T=800.0, rxns=None):
    rs = HybridPolymerSystem(
        T=T, P=1.0e5, initial_mole_fractions={core[0]: 1.0},
        V_poly=1.0, polymer_pools=pools, mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments=moments, termination=[])
    rs.initialize_model(core, rxns or [], [], [])
    return rs


class TestSideGroupHomolysisKernel:

    def test_validator_field_rules(self):
        """Strict-key validation for the side_group_homolysis channel list
        (single source of truth: rmgpy.solver.polymer.
        validate_side_group_homolysis). Every malformed shape fails loudly,
        never a laundered no-op."""
        from rmgpy.solver.polymer import validate_side_group_homolysis

        good = _sgh_channel()
        out = validate_side_group_homolysis("P", [good])
        assert isinstance(out, list) and len(out) == 1
        assert out[0]["A"] == 1.0e13 and out[0]["sites_per_unit"] == 2.0
        assert set(out[0]) == {"label", "A", "n", "Ea", "site_selector",
                               "sites_per_unit", "gas_product"}

        cases = [
            (dict(good), r"P.*side_group_homolysis.*list"),          # not a list
            ([42.0], r"P.*side_group_homolysis.*dict"),               # entry not a dict
            ([], r"P.*side_group_homolysis.*empty"),                  # empty list
            ([{**good, "extra": 1.0}], r"P.*unknown"),                # unknown key
            ([{k: v for k, v in good.items() if k != "Ea"}], r"P.*missing"),
            # round-72 P1: the structural site selector is REQUIRED and
            # vocabulary-pinned even at the shape layer (the solver config
            # carries no monomer structure; the structural layers live at
            # deck / to_config / producer).
            ([{k: v for k, v in good.items() if k != "site_selector"}],
             r"P.*missing.*site_selector"),
            ([{**good, "site_selector": "halogenated"}],
             r"P.*site_selector.*one of"),
            ([{**good, "site_selector": 1.0}],
             r"P.*site_selector.*one of"),
            ([{**good, "A": -1.0}], r"P.*A.*> 0"),
            ([{**good, "A": float("nan")}], r"P.*A.*not finite"),
            ([{**good, "Ea": -5.0}], r"P.*Ea.*>= 0"),
            ([{**good, "sites_per_unit": 0.0}], r"P.*sites_per_unit.*> 0"),
            ([{**good, "sites_per_unit": True}], r"P.*sites_per_unit.*number"),
            ([{**good, "label": ""}], r"P.*label.*non-empty"),
            ([{**good, "gas_product": "not a smiles ]["}], r"P.*gas_product.*parse"),
            ([{**good, "gas_product": "CC"}], r"P.*gas_product.*radical"),
            ([{**good, "gas_product": "[CH3]"}], r"P.*gas_product.*monoatomic"),
            ([good, dict(good)], r"P.*duplicate"),
            # sanitized-collision duplicates are duplicates too
            ([good, {**good, "label": "aliphatic_C+Br"}], r"P.*duplicate"),
        ]
        for bad, pattern in cases:
            with pytest.raises(ValueError, match=pattern):
                validate_side_group_homolysis("P", bad)

    def test_side_group_inventory_law_parent_stationary_gas(self):
        """SGH kernel-v2 inventory law, hand-computed at mu = (1, 5, 30),
        T = 800 K, s = 2, V_poly = 1, Z(0) = sites*mu1(0) = 2*5 = 10:

          k(800)   = A*800^0.5*exp(-Ea/(R*800))
          R        = k*Z / V_poly = k*10   [mol/(m^3 s); NO extra sites factor,
                                            sites already folded into Z(0)]
          dZ/dt    = -R*V_poly = -k*10     (Br inventory drains)
          gas Br   = +R*V_poly = +k*10     (one X radical per event)

        The real v2 change: the carrier's mu0/mu1/mu2 are STATIONARY (no moment
        debit, no feature-pool credit) -- Br leaves through the Z inventory, not
        the moments. The gas rate is numerically identical to v1's k*s*mu1 at
        t0 (because Z(0) = s*mu1(0)), so moment-stationarity is what pins the
        change. dZ/dt and the gas credit are bit-exact negatives (shared R)."""
        core, mask, pools, moments = _sgh_core_and_pools()
        rs = _sgh_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        sites, mu1 = 2.0, 5.0
        z0 = sites * mu1                      # 10 == Z(0)
        R = _sgh_arrhenius(800.0) * z0        # k*Z/V_poly, V_poly = 1
        z_slot = int(rs.sgh_z_slot[0])
        assert z_slot == len(core)            # n_core + 0 QSSA U slots

        # parent moments are STATIONARY (the v2 correction over v1's debit):
        assert dn_dt[3] == 0.0
        assert dn_dt[4] == 0.0
        assert dn_dt[5] == 0.0
        # Z inventory drains at R; gas X gains R:
        assert np.isclose(dn_dt[z_slot], -R, rtol=1e-12)
        assert np.isclose(dn_dt[1], R, rtol=1e-12)
        # gas gain and inventory loss are bit-exact negatives (shared R):
        assert dn_dt[1] == -dn_dt[z_slot]
        # inert gas and the poly proxy carry no kernel flux:
        assert dn_dt[0] == 0.0 and dn_dt[2] == 0.0

    def test_side_group_mass_conservation_pin(self):
        """LOAD-BEARING mass closure (SGH kernel-v2). Br leaves the melt through
        the Z inventory; get_total_polymer_condensed_mass_g books the carrier
        mass as mu1*MW - Sum_row max(0, sites*mu1 - Z)*M_X, and the gas X radical
        gains the same moles. With the carrier mu1 STATIONARY, condensed mass =
        mu1*MW - sites*mu1*M_X + Z*M_X, so d(condensed)/dt = dZ/dt * M_X and
        d(gas X mass)/dt = dn_gas * M_X. The RHS derives both from the SAME rate
        R (gas += R, dZ -= R), so d(condensed)/dt + d(gas X)/dt == 0.0 BITWISE.
        v2 carriers carry NO chain-mass defect (Br mass rides the Z path)."""
        core, mask, pools, moments = _sgh_core_and_pools()
        rs = _sgh_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        m_x = _sgh_mx_g_mol()                 # M_X of the ejected Br [g/mol]
        parent = pools[0]
        z_slot = int(rs.sgh_z_slot[0])
        # v2 carrier holds NO chain-mass defect (there is no feature pool):
        assert parent.chain_mass_defect_g_mol == 0.0
        assert rs.sgh_mw_X[0] == m_x

        # mu1 stationary (dn_dt[4] == 0) => d(condensed)/dt = dZ/dt * M_X:
        d_condensed = dn_dt[z_slot] * m_x     # = -R*M_X
        d_gas_mass = dn_dt[1] * m_x           # = +R*M_X
        # bit-exact closure (dn_dt[1] == -dn_dt[z_slot] from the shared R):
        assert d_condensed + d_gas_mass == 0.0, (
            f"kernel must conserve condensed+gas X mass exactly, got "
            f"{d_condensed + d_gas_mass!r}")
        assert d_gas_mass > 0.0 and d_condensed == -d_gas_mass

        # Cross-check the mass tally itself against the Z path: removing delta
        # moles of X from the melt (Z <- Z0 - delta, still in [0, sites*mu1])
        # drops the condensed mass by EXACTLY delta*M_X.
        mass0 = rs.get_total_polymer_condensed_mass_g(rs.y)
        y2 = rs.y.copy()
        delta = 1.0                           # mol of X removed (Z0=10 -> 9)
        y2[z_slot] = rs.y[z_slot] - delta
        mass1 = rs.get_total_polymer_condensed_mass_g(y2)
        assert np.isclose(mass0 - mass1, delta * m_x, rtol=1e-12)

    def test_side_group_arrhenius_evaluated_at_runtime_temperature(self):
        """k(T) is evaluated at the solver's runtime T (n=0.5 also catches a
        dropped T^n factor). The gas rate is k(T)*Z(0), Z(0)=sites*mu1(0)=10
        (kernel-v2: k*Z, NOT k*sites*mu1)."""
        rates = {}
        z0 = 2.0 * 5.0                        # sites*mu1(0) = Z(0)
        for T in (700.0, 1100.0):
            core, mask, pools, moments = _sgh_core_and_pools()
            rs = _sgh_system(core, mask, pools, moments, T=T)
            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            rates[T] = dn_dt[1]
            expected = _sgh_arrhenius(T) * z0        # k(T) * Z(0)
            assert np.isclose(rates[T], expected, rtol=1e-10), (
                f"k({T}) mismatch: got {rates[T]}, expected {expected}")
        assert rates[1100.0] > rates[700.0] * 10.0

    def test_side_group_zero_flux_on_empty_pool(self):
        """Empty carrier: mu1(0)=0 -> Z(0)=sites*mu1(0)=0 -> R=k*Z/V_poly=0.
        Exactly zero flux everywhere (zero inventory is in-domain, not an
        error)."""
        core, mask, pools, moments = _sgh_core_and_pools(
            parent_moments=(0.0, 0.0, 0.0))
        rs = _sgh_system(core, mask, pools, moments)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(dn_dt == 0.0)

    def test_side_group_residual_z_excursion_censuses_not_raises(self, caplog):
        """P1-A trial-state safety (corrected design): residual() runs on the
        RAW Newton trial state, so an out-of-band Z must NOT raise. Z decays
        stiffly (dZ/dt = -kZ) and a recoverable Newton trial can undershoot
        -floor as Z -> 0; hard-raising there would turn a transient into a
        fatal integration error. The residual only CENSUSES (warn-once per
        (pool, channel)) and the RHS's local max(0., .) read sanitizes a
        negative trial Z; y is never mutated by the guard."""
        core, mask, pools, moments = _sgh_core_and_pools()
        rs = _sgh_system(core, mask, pools, moments, T=800.0)
        z_slot = int(rs.sgh_z_slot[0])
        upper = 2.0 * 5.0                     # sites*mu1 = Z(0)

        with caplog.at_level(logging.WARNING):
            # (a) trial Z FAR below -floor: NO raise (pre-P1-A this was fatal).
            y_lo = rs.y.copy()
            y_lo[z_slot] = -1.0
            dn = rs.residual(0.0, y_lo, np.zeros_like(y_lo))[0]
            # the RHS sanitizes the negative Z (max(0, Z) = 0): the SGH rate is
            # exactly zero, not negative -- Z slot and gas slot stay put.
            assert dn[z_slot] == 0.0 and dn[1] == 0.0
            n_first = sum("INVENTORY CENSUS" in r.getMessage()
                          for r in caplog.records)
            # (b) trial Z well ABOVE sites*mu1 + floor: also census-only.
            y_hi = rs.y.copy()
            y_hi[z_slot] = upper + 5.0
            rs.residual(0.0, y_hi, np.zeros_like(y_hi))      # must not raise
            # (c) warn-once: re-entering the same (pool, channel) adds no line.
            rs.residual(0.0, y_lo, np.zeros_like(y_lo))
            n_after = sum("INVENTORY CENSUS" in r.getMessage()
                          for r in caplog.records)
        assert n_first == 1, "an out-of-band trial Z must census loudly, once"
        assert n_after == n_first, "warn-once: no repeat across residual calls"
        # the guard never mutates the trial state:
        assert y_lo[z_slot] == -1.0 and y_hi[z_slot] == upper + 5.0

    def test_side_group_accepted_state_z_bound_hard_raises(self):
        """P1-A accepted-state invariant (corrected design): the HARD physical-
        band check the residual guard is deliberately too trial-unsafe to make
        now lives in _assert_sgh_inventory_accepted, which reads self.y (the
        ACCEPTED integrator solution). A Z beyond the atol-derived floor =
        max(SMALL_EPS, EXHAUSTION_FLOOR_K*atol[z_slot]) (below -floor OR above
        sites*mu1+floor) is genuine integrator corruption and HARD-RAISES; a
        small in-band excursion does not."""
        core, mask, pools, moments = _sgh_core_and_pools()
        rs = _sgh_system(core, mask, pools, moments, T=800.0)
        z_slot = int(rs.sgh_z_slot[0])
        upper = 2.0 * 5.0                     # sites*mu1 = Z(0)
        floor_z = max(1e-30, 100.0 * float(rs.atol_array[z_slot]))

        # (a) accepted Z below -floor -> hard error.
        rs.y[z_slot] = -1.0
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*ACCEPTED state"):
            rs._assert_sgh_inventory_accepted()

        # (b) accepted Z above sites*mu1 + floor -> equally refused.
        rs.y[:] = rs.y0
        rs.y[z_slot] = upper + 1.0
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*ACCEPTED state"):
            rs._assert_sgh_inventory_accepted()

        # (c) a small in-band negative (within -floor) does NOT raise.
        rs.y[:] = rs.y0
        rs.y[z_slot] = -0.25 * floor_z
        rs._assert_sgh_inventory_accepted()   # must not raise

    def test_side_group_mass_closure_at_nonunit_v_poly(self):
        """Gap A: the SGH mass balance and amount/concentration algebra must
        close at V_poly != 1. Moments are MOLES (mu1 slot = 5 regardless of
        V_poly) and Z(0) = sites*mu1(0) is on the AMOUNT basis, so at V_poly=2:
          R      = k*Z/V_poly   [mol/(m^3 s)]  (concentration-basis rate)
          dZ/dt  = -R*V_poly    [mol/s]        (amount drains)
          gas X += R*V_poly     [mol/s]        (small_src -> dn_dt*V_poly)
        d(condensed X mass)/dt + d(gas X mass)/dt == 0 BITWISE, and Z(0) is
        V_poly-independent (a moles quantity)."""
        core, mask, pools, moments = _sgh_core_and_pools()
        V_poly = 2.0
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0},
            V_poly=V_poly, polymer_pools=pools, mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments=moments, termination=[])
        rs.initialize_model(core, [], [], [])
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        z_slot = int(rs.sgh_z_slot[0])
        z0 = 2.0 * 5.0                        # sites*mu1(0), MOLES (V_poly-indep)
        assert np.isclose(rs.y0[z_slot], z0, rtol=1e-12)
        R = _sgh_arrhenius(800.0) * z0 / V_poly     # k*Z/V_poly (conc basis)
        # amount algebra closes at non-unit volume: dZ = -R*V_poly, gas =
        # +R*V_poly (both equal r_sgh*V_poly = k*z0):
        assert np.isclose(dn_dt[z_slot], -R * V_poly, rtol=1e-12)
        assert np.isclose(dn_dt[1], R * V_poly, rtol=1e-12)
        assert dn_dt[1] == -dn_dt[z_slot]     # bit-exact negatives
        # moments stationary at non-unit volume too:
        assert dn_dt[3] == 0.0 and dn_dt[4] == 0.0 and dn_dt[5] == 0.0
        # mass closure (mu1 stationary => d(condensed)/dt = dZ/dt * M_X):
        m_x = _sgh_mx_g_mol()
        assert dn_dt[z_slot] * m_x + dn_dt[1] * m_x == 0.0
        # the condensed-mass tally is on the MOLE basis (V_poly-independent):
        mass0 = rs.get_total_polymer_condensed_mass_g(rs.y)
        y2 = rs.y.copy()
        y2[z_slot] = rs.y[z_slot] - 1.0       # remove 1 mol of X (Z: 10 -> 9)
        assert np.isclose(mass0 - rs.get_total_polymer_condensed_mass_g(y2),
                          1.0 * m_x, rtol=1e-12)

    def test_side_group_coexists_with_k_scission_and_mass_closes(self):
        """Gap B: k_scission is the one mu1-PRESERVING kernel allowed alongside
        SGH. A carrier with k_scission > 0 AND side_group_homolysis must build
        (no raise); the SGH Z law still fires and mass still closes; and
        k_scission still moves mu0/mu2 while leaving mu1 (hence Z and the Br
        mass) stationary."""
        k_s = 0.3
        core, mask, pools, moments = _sgh_core_and_pools()
        pools[0] = dataclasses.replace(pools[0], k_scission=k_s)
        rs = _sgh_system(core, mask, pools, moments, T=800.0)   # must NOT raise
        assert rs.sgh_enabled[0] == 1
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        z_slot = int(rs.sgh_z_slot[0])
        z0 = 2.0 * 5.0
        R = _sgh_arrhenius(800.0) * z0        # k*Z/V_poly, V_poly = 1
        # SGH Z law is still LIVE next to k_scission:
        assert np.isclose(dn_dt[z_slot], -R, rtol=1e-12)
        assert np.isclose(dn_dt[1], R, rtol=1e-12)
        # mu1 (index 4) stationary -> Z and the Br mass stay pinned:
        assert dn_dt[4] == 0.0
        # k_scission MOVES mu0 (index 3) and mu2 (index 5); mu = (1, 5, 30),
        # mu3 = mu0*(mu2/mu1)^3 (log-Lagrange closure) = 216:
        mu0, mu1, mu2 = 1.0, 5.0, 30.0
        mu3 = mu0 * (mu2 / mu1) ** 3
        assert np.isclose(dn_dt[3], k_s * (mu1 - mu0), rtol=1e-9)         # +1.2
        assert np.isclose(dn_dt[5], k_s * (mu1 - mu3) / 3.0, rtol=1e-9)   # < 0
        assert dn_dt[3] > 0.0 and dn_dt[5] < 0.0
        # mass still closes (mu1 stationary => d(condensed)/dt = dZ/dt * M_X):
        m_x = _sgh_mx_g_mol()
        assert dn_dt[z_slot] * m_x + dn_dt[1] * m_x == 0.0

    def test_side_group_rejects_live_resolved_row_on_carrier(self):
        """Gap C: the runtime-row-flux arm of the SGH mu1-mover build guard. A
        NON-refused live flux-archetype row whose src/dst is the SGH carrier
        moves the carrier's mu1 -- which inventory-depletion SGH forbids -- and
        must be REFUSED at initialize_model (build time). Here a non-refused
        proxy-consuming row resolves to FLUX_UNRESOLVED with src_pool = the SGH
        carrier (the allow_unstamped_proxy_rows escape hatch), the live analogue
        of the refused (zero-flux, allowed) supersession row."""
        core, mask, pools, moments = _sgh_core_and_pools()
        # non-refused row: SGH carrier proxy (core[2] 'poly') -> gas Br (core[1]).
        live_row = Reaction(
            reactants=[core[2]], products=[core[1]],
            kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "J/mol"),
                               T0=(1.0, "K")),
            reversible=False)
        rs = HybridPolymerSystem(
            allow_unstamped_proxy_rows=True,   # resolve the proxy row to UNRESOLVED
            T=800.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0},
            V_poly=1.0, polymer_pools=pools, mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments=moments, termination=[])
        with pytest.raises(ValueError,
                           match=r"poly.*stationary mu1.*FLUX_UNRESOLVED"):
            rs.initialize_model(core, [live_row], [], [])

    def test_side_group_rejects_explicit_dp_on_carrier(self):
        """Gap D (P1-B, solver-build): an SGH carrier that ALSO carries a
        non-empty explicit_dp_to_species_index is a silent mass-corruption hole
        (get_total_polymer_condensed_mass_g counts explicit-DP repeat-unit mass
        in mu1 while the SGH lost-X term uses tail-only mu1) and must be REFUSED
        at initialize_model. The deck refuses this pairing, but a direct
        PolymerPoolConfig (here) bypasses the deck -- the solver-build guard is
        the backstop for BOTH deck-built and directly-built configs."""
        core, mask, pools, moments = _sgh_core_and_pools()
        # index 2 (the condensed proxy) is in range and non-gas, so the early
        # explicit-DP index checks pass and the SGH build guard is what fires.
        pools[0] = dataclasses.replace(
            pools[0], explicit_dp_to_species_index={3: 2})
        with pytest.raises(ValueError,
                           match=r"SGH-v2 carrier pool 'poly'.*explicit_dp"):
            _sgh_system(core, mask, pools, moments)

    def test_side_group_degenerate_kT_raises(self):
        """A k(T) evaluating to inf/NaN/0-overflow is a poisoned kernel --
        hard error, never integrate garbage (sibling convention)."""
        core, mask, pools, moments = _sgh_core_and_pools(
            channels=[_sgh_channel(A=1.0e300, n=100.0)])
        with pytest.raises(ValueError, match=r"poly.*side_group_homolysis.*degenerate"):
            # the poisoned k(T) fires on the FIRST residual evaluation --
            # initialize_model's set_initial_derivative already runs one
            rs = _sgh_system(core, mask, pools, moments)
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))

    def test_side_group_seeds_and_flattens_inventory_state(self):
        """SGH kernel-v2 flatten/seed contract (replaces the v1 feature-pool
        resolution + mass-defect-pin guards). One Z aux-state per (pool,
        channel), placed at the absolute ODE index n_core + num_qssa_u + row,
        seeded Z(0) = sites*mu1(0) (moles, matching the mu1 amount basis), with
        the per-row M_X recorded on sgh_mw_X and neq extended by num_sgh_z. v2
        deletes the feature pool, so the carrier carries no chain-mass defect."""
        core, mask, pools, moments = _sgh_core_and_pools()
        rs = _sgh_system(core, mask, pools, moments, T=800.0)
        n_core = len(core)
        assert rs.num_sgh_z == 1
        assert int(rs.sgh_z_slot[0]) == n_core + int(rs.num_qssa_u)   # = 6
        assert rs.neq == n_core + int(rs.num_qssa_u) + int(rs.num_sgh_z)
        # Z(0) = sites*mu1(0) = 2*5 = 10 seeded on the aux slot:
        assert np.isclose(rs.y0[int(rs.sgh_z_slot[0])], 2.0 * 5.0, rtol=1e-12)
        # per-row M_X of the ejected Br radical (the value v1 discarded):
        assert rs.sgh_mw_X[0] == _sgh_mx_g_mol("[Br]")
        # no feature pool in v2 => the carrier carries no chain-mass defect:
        assert pools[0].chain_mass_defect_g_mol == 0.0

    def test_side_group_requires_resolved_gas_indices(self):
        """An enabled channel with no resolved gas-product index would mint
        condensed-mass loss with no gas destination -- hard error."""
        core, mask, pools, moments = _sgh_core_and_pools(gas_indices=None)
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*gas"):
            _sgh_system(core, mask, pools, moments)

    def test_side_group_gas_index_must_be_gas_masked(self):
        """The gas_product index must be classified GAS: emitting the X
        radical into a condensed slot would corrupt the phase bookkeeping."""
        core, mask, pools, moments = _sgh_core_and_pools(gas_indices=(2,))
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*GAS"):
            _sgh_system(core, mask, pools, moments)

    def test_solver_rejects_side_group_with_positive_k_unzip(self):
        """Exclusion (solver last line of defense): legacy k_unzip and
        side_group_homolysis double-carry degradation on one pool."""
        Mono = _spc("C=CC", "propene_gas")
        core, mask, pools, moments = _sgh_core_and_pools(extra_species=[Mono])
        pools[0] = dataclasses.replace(
            pools[0], k_unzip=0.4, monomer_poly_index=len(core) - 1)
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*k_unzip.*mutually exclusive"):
            _sgh_system(core, mask, pools, moments)

    def test_solver_rejects_side_group_with_qssa_channel(self):
        """Exclusion (solver last line of defense): QSSA random initiation and
        side_group_homolysis on one pool double-carry initiation."""
        qssa = dict(initiation=dict(A=1.0e15, n=0.0, Ea=3.0e5),
                    depropagation=dict(A=1.0e13, n=0.0, Ea=8.0e4),
                    termination=dict(A=1.0e8, n=0.0, Ea=1.0e4))
        Mono = _spc("C=CC", "propene_gas")
        core, mask, pools, moments = _sgh_core_and_pools(extra_species=[Mono])
        pools[0] = dataclasses.replace(
            pools[0], radical_qssa_unzip=qssa,
            monomer_poly_index=len(core) - 1)
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*radical_qssa_unzip.*mutually exclusive"):
            _sgh_system(core, mask, pools, moments)

    def test_side_group_rejects_coexistence_with_k_homolysis(self):
        """SGH kernel-v2 mu1-mover guard (contract 1h + P1-B): inventory-
        depletion SGH seeds Z(0) = sites*mu1(0) ONCE and depletes it
        independently, so it REQUIRES a stationary carrier mu1. A carrier that
        also declares a mu1-moving kernel -- here k_homolysis, which debits the
        backbone-break moments -- is REFUSED at setup (a moving mu1 would strand
        or double-count the Br inventory). k_scission stays allowed elsewhere;
        this is the k_homolysis arm of the guard."""
        core, mask, pools, moments = _sgh_core_and_pools(
            k_homolysis=_khom_triplet())
        with pytest.raises(ValueError,
                           match=r"poly.*stationary mu1.*k_homolysis"):
            _sgh_system(core, mask, pools, moments)

    def test_side_group_multi_channel_one_carrier_z_rows(self):
        """SGH kernel-v2 multi-channel = multiple Z inventory rows on ONE
        carrier (no separate feature pools). Aliphatic (s=2) and aryl (s=4)
        C-Br each get their own Z slot seeded Z(0)=s*mu1(0) and drain at
        R=k*Z/V_poly; the shared gas Br slot receives the sum of the per-row
        event rates, and each row's M_X lands on sgh_mw_X. The carrier's
        moments stay stationary."""
        ch_a = _sgh_channel()                                  # s = 2
        ch_b = _sgh_channel(label="aryl_C-Br", A=5.0e12, n=0.5,
                            Ea=2.0e5, site_selector="aryl",
                            sites_per_unit=4.0)                # s = 4
        core, mask, pools, moments = _sgh_core_and_pools(channels=[ch_a, ch_b])
        rs = _sgh_system(core, mask, pools, moments, T=800.0)
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        mu1 = 5.0
        n_core = len(core)
        # two Z rows on the one carrier, at consecutive absolute slots:
        assert rs.num_sgh_z == 2
        z_a = int(rs.sgh_z_slot[0])
        z_b = int(rs.sgh_z_slot[1])
        assert z_a == n_core + int(rs.num_qssa_u)              # = 6
        assert z_b == n_core + int(rs.num_qssa_u) + 1          # = 7
        # per-row rates R = k*Z(0)/V_poly, Z(0)=sites*mu1(0), V_poly=1:
        r_a = _sgh_arrhenius(800.0, ch_a) * (2.0 * mu1)        # k_a*10
        r_b = _sgh_arrhenius(800.0, ch_b) * (4.0 * mu1)        # k_b*20
        # each Z inventory drains at its own R:
        assert np.isclose(dn_dt[z_a], -r_a, rtol=1e-12)
        assert np.isclose(dn_dt[z_b], -r_b, rtol=1e-12)
        # the carrier's moments are stationary (no per-channel debit):
        assert dn_dt[3] == 0.0 and dn_dt[4] == 0.0 and dn_dt[5] == 0.0
        # shared gas slot: sum of the two event rates:
        assert np.isclose(dn_dt[1], r_a + r_b, rtol=1e-12)
        # bit-exact closure: gas gain == total inventory loss:
        assert dn_dt[1] == -(dn_dt[z_a] + dn_dt[z_b])
        # per-channel M_X recorded on sgh_mw_X (both eject Br here):
        m_x = _sgh_mx_g_mol("[Br]")
        assert rs.sgh_mw_X[0] == m_x and rs.sgh_mw_X[1] == m_x

    def test_solver_rejects_duplicate_channel_labels(self):
        """Duplicate channel labels on one pool are a hard error (two decks
        of the same channel double-carry the same bond class)."""
        core, mask, pools, moments = _sgh_core_and_pools(
            channels=[_sgh_channel(), _sgh_channel(A=2.0e13)],
            gas_indices=(1, 1))
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*duplicate"):
            _sgh_system(core, mask, pools, moments)

    def test_solver_rejects_malformed_side_group_channels(self):
        """Solver field validation: a directly-constructed config with a
        malformed channel list must fail loudly, not run a garbage rate."""
        core, mask, pools, moments = _sgh_core_and_pools(
            channels=[dict(label="x", A=1.0e13, n=0.0)], gas_indices=(1,))
        with pytest.raises(ValueError,
                           match=r"poly.*side_group_homolysis.*missing"):
            _sgh_system(core, mask, pools, moments)

    def test_refused_rows_zero_flux_and_census_pairs_side_group_kernel(self):
        """Refused explicit gas-radical<->condensed rows remain refused
        (supersession, same contract as the sibling kernel): the
        conduit-deferred row carries EXACTLY zero flux next to the live
        kernel, and the census pairs them explicitly."""
        R1 = _spc("[CH2]C(C)C", "isobutyl")
        R2 = _spc("C[CH]CCC", "pentyl2")
        core, mask, pools, moments = _sgh_core_and_pools(
            extra_species=[R1, R2])
        proxy = core[2]
        refused = Reaction(reactants=[R1, R2], products=[proxy],
                           kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0,
                                              Ea=(0.0, "J/mol"), T0=(1.0, "K")),
                           reversible=False)
        refused.polymer_refused = True
        refused.polymer_refused_accumulating = False    # conduit-deferred

        rs_with = _sgh_system(core, mask, pools, moments, rxns=[refused])
        dn_with = rs_with.residual(0.0, rs_with.y, np.zeros_like(rs_with.y))[0]

        core2, mask2, pools2, moments2 = _sgh_core_and_pools(
            extra_species=[_spc("[CH2]C(C)C", "isobutyl"),
                           _spc("C[CH]CCC", "pentyl2")])
        rs_wo = _sgh_system(core2, mask2, pools2, moments2, rxns=[])
        dn_wo = rs_wo.residual(0.0, rs_wo.y, np.zeros_like(rs_wo.y))[0]

        assert np.any(dn_wo != 0.0)                 # kernel flux is live
        assert np.array_equal(dn_with, dn_wo), (
            "refused row must add exactly zero flux next to the live kernel")
        assert any(c["reason"] == "conduit-deferred"
                   for c in rs_with.refused_reaction_census)
        census = rs_with.side_group_supersession_census
        assert census, "kernel-live census must pair the refused rows"
        entry = next(c for c in census if c["pool"] == "poly")
        assert len(entry["superseded_rows"]) == 1

    def test_side_group_diagnostic_lists_channels_and_inventory(self, capsys):
        """The POLYMER SOLVER DIAGNOSTIC block names the side_group_homolysis
        channels on the carrier pool and, for kernel-v2, reports each channel's
        seeded Z inventory (Z(0)=sites*mu1(0)) and Z slot -- NOT the dropped v1
        feature pool or the v1 saturation banner."""
        core, mask, pools, moments = _sgh_core_and_pools()
        _sgh_system(core, mask, pools, moments)
        out = capsys.readouterr().out
        assert "POLYMER SOLVER DIAGNOSTIC" in out
        assert "side_group_homolysis" in out
        assert "aliphatic_C-Br" in out
        # v2 Z-inventory disclosure (replaces the v1 feature pool + saturation):
        assert "Z inventory slot" in out
        assert "Z(0)=" in out
        # the v1 feature pool no longer exists, so its label never appears:
        assert "poly_sidegrp_aliphatic_C_Br" not in out


# Bitwise pre-change pin for the live-VE negative control below, recorded on
# b0de7dde8 (pre-fix HEAD) with the exact fixture in
# test_negative_control_live_ve_row_with_real_polymer_keeps_flux:
# dn[:5] = [proxy, mu0, mu1, mu2, volatile].
_VE_CONTROL_PIN = np.array([0.0, 0.0, -100.0, -11900.0, 100.0])


class TestRefusalAdjudicationSurvivesRebuild:
    """r71 adjudication (PP run-5 solver stall, forensics
    /home/alon/Projects/polymer/PP/rmg/run5): the polymer_refused
    adjudication was LOST between generation-time stamping and the
    post-promotion solver rebuild, and refused EDGE rows still fed
    edge_species_rates -- which is exactly HOW the dead chain-scale radicals
    were promoted to core, where the unstamped rows ran live legacy mu1-only
    flux against an exhausted pool (mu ~ 1e-10) and collapsed DASPK.

    RED-FIRST discipline: every red assertion below was quoted FAILING on
    pre-fix HEAD (b0de7dde8) for its pinned reason before the fix landed.
    Liveness pins precede each red assert so a red can only mean "guard
    absent", never "fixture dead".
    """

    @staticmethod
    def _run5_edge_fixture():
        """Association-orientation run-5 shape: gas radical (core) +
        chain-scale radical (EDGE, prospectively CONDENSED -- the edge-daughter
        condensed-mask kept run-5's gates open) <=> pool proxy (core),
        REVERSIBLE. The promotion flux is the reverse (homolysis-discovery)
        leg rr; rf vanishes on the edge reactant's zero concentration."""
        sp = {
            "P": _spc("CC(C)CC(C)CCC", "poly"),
            "mu0": _spc("CO", "poly_mu0"), "mu1": _spc("C=O", "poly_mu1"),
            "mu2": _spc("C#N", "poly_mu2"),
            "R1": _spc("[CH2]C(C)C", "isobutyl"),
        }
        edge_r = _spc("C[CH]CCC", "pentyl2_edge")
        for s in list(sp.values()) + [edge_r]:
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [sp["P"], sp["mu0"], sp["mu1"], sp["mu2"], sp["R1"]]
        mask = np.array([False] * 4 + [True], dtype=bool)
        rxn = Reaction(
            reactants=[sp["R1"], edge_r], products=[sp["P"]],
            kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0,
                               Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=True)
        # Stamped at generation (round-63 refusal), exactly as run-5's rows
        # were (RMG.log:648-660 census).
        rxn.polymer_refused = True
        rxn.polymer_refused_accumulating = False   # conduit-deferred
        return sp, edge_r, core, mask, rxn

    @staticmethod
    def _pool(label, mu_indices, monomer_mw_g_mol=42.08):
        return PolymerPoolConfig(
            label=label, xs=2, explicit_dp_to_species_index={},
            mu_indices=mu_indices, monomer_poly_index=None,
            monomer_mw_g_mol=monomer_mw_g_mol,
            k_scission=0.0, k_unzip=0.0)

    def test_red_a_refused_edge_row_is_flux_dead_for_enlargement(self):
        """RED-A edge half (FIX 3): a refused EDGE row must contribute ZERO
        to edge_species_rates / edge_reaction_rates (the enlargement inputs
        that promoted run-5's dead radicals). The per-reaction ungated
        counterfactual stays live (diagnostic only, never flux)."""
        sp, edge_r, core, mask, rxn = self._run5_edge_fixture()
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[-1]: 0.0},
            V_poly=1.0, polymer_pools=[self._pool("poly", (1, 2, 3))],
            mass_transfer=[], gas_species_mask=mask.copy(),
            constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 50.0, 3000.0)},
            termination=[],
            allow_default_prospective_edge=True,
            prospective_condensed_edge_daughter_classifier=(
                lambda spcs: {"pentyl2_edge"}),
        )
        rs.initialize_model(core, [], [edge_r], [rxn])
        assert rs.reaction_refused[0] == 1   # liveness: stamp visible to solver
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))

        # LIVENESS PIN -- BEFORE the red asserts: the row genuinely carries
        # reverse flux (the gates are open: prospectively-condensed edge
        # daughter, condensed product). The per-reaction ungated
        # counterfactual sees it, so the zeros below cannot mean "fixture
        # dead" (Gate A/B zeroing would zero this too).
        assert rs.edge_reaction_rates_ungated[0] != 0.0, (
            "FIXTURE BROKEN, not a valid red: the refused edge row carries "
            "no reverse flux at all")

        # THE RED asserts (FIX 3): flux-dead in every enlargement input.
        assert rs.edge_reaction_rates[0] == 0.0, (
            "refused edge row leaked into edge_reaction_rates")
        assert rs.edge_species_rates[0] == 0.0, (
            "refused edge row leaked promotion flux into edge_species_rates "
            "-- this is HOW run-5's dead radicals reached core")

    def test_red_a_adjudication_survives_promotion_and_rebuild(self):
        """RED-A promotion half (FIX 2): build the row through the REAL
        generation route (make_new_reaction stamps it), place it in edge,
        lose the adjudication on the canonical object (the run-5 arrival
        state, r71: 'the rows arrived unstamped'), promote the missing
        radical through the REAL path (add_species_to_core moves the row
        edge->core), rebuild the solver -- the rebuild restamp (last honest
        chokepoint) must re-derive the refusal regardless of how it was
        lost upstream."""
        import rmgpy.data.rmg as rmg_data
        from rmgpy.data.base import ForbiddenStructures
        from rmgpy.data.kinetics import TemplateReaction
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer
        from rmgpy.rmg.model import CoreEdgeReactionModel
        import rmgpy.solver.polymer as solver_mod

        old_db = rmg_data.database
        db = RMGDatabase()
        db.forbidden_structures = ForbiddenStructures()
        rmg_data.database = db
        try:
            cerm = CoreEdgeReactionModel()
            pp = Polymer(label='polypropylene', monomer='[CH2][CH]C',
                         Mn=5000.0, Mw=8000.0, initial_mass=1.0)
            cerm._register_polymer(pp, generate_thermo=False)
            forward = TemplateReaction(
                reactants=[Molecule().from_smiles('[CH2]C(C)C'),
                           Molecule().from_smiles('C[CH]CCC')],
                products=[Molecule().from_smiles('CC(C)CC(C)CCC')],
                family='R_Recombination', is_forward=True, reversible=True,
                kinetics=Arrhenius(A=(1.0e7, 'm^3/(mol*s)'), n=0.0,
                                   Ea=(0.0, 'kJ/mol')))
            out, is_new = cerm.make_new_reaction(
                forward, check_existing=False, generate_thermo=False,
                generate_kinetics=False)
            # LIVENESS PINS: real generation route, product resolved onto the
            # registered pool Polymer, refusal stamped at generation.
            assert out is not None and is_new
            assert any(isinstance(p, Polymer) for p in out.products)
            assert out.polymer_refused is True

            r1_spc = next(s for s in out.reactants
                          if s.molecule[0].get_formula() == 'C4H9')
            r2_spc = next(s for s in out.reactants
                          if s.molecule[0].get_formula() == 'C5H11')

            cerm.add_species_to_edge(pp)      # brings the mu dummies along
            cerm.add_species_to_core(pp)      # moves pp + dummies to core
            cerm.add_species_to_core(r1_spc)
            cerm.add_species_to_edge(r2_spc)
            cerm.add_reaction_to_edge(out)

            # The adjudicated run-5 arrival state (r71, binding): by rebuild
            # time the canonical row is UNSTAMPED -- the stamps are ad-hoc
            # object attributes, deliberately NOT serialized in __reduce__
            # (rmgpy/reaction.py: "the solver re-derives them at
            # initialize_model" -- the re-derivation FIX 2 forces into
            # existence), and check_for_existing_reaction can discard a
            # stamped candidate for an unstamped canonical (FIX 1's seam,
            # pinned live in modelTest). Simulate that loss on the canonical
            # edge object.
            out.polymer_refused = False
            out.polymer_refused_accumulating = False

            # REAL promotion path (run 5, RMG.log:917 'Moved 1 reactions from
            # edge to core: H(1) + CCCC(C)C[C](C)C(16) <=> polypropylene(2)').
            moved = cerm.add_species_to_core(r2_spc)
            assert moved == [out]                       # liveness
            assert out in cerm.core.reactions
            assert out not in cerm.edge.reactions

            core_species = list(cerm.core.species)
            for s in core_species + list(cerm.edge.species):
                s.thermo = _trivial_nasa(_GAV_COMMENT)
            # Polymer thermo delegates to its internal proxy Species
            # (get_free_energy -> get_proxy_species().thermo): the reversible
            # row's Keq needs it resolvable without a thermo database.
            pp.baseline_proxy.thermo = _trivial_nasa(_GAV_COMMENT)
            idx = {s.label: i for i, s in enumerate(core_species)}
            mask = np.ones(len(core_species), dtype=bool)
            for lbl in ('polypropylene', 'polypropylene_mu0',
                        'polypropylene_mu1', 'polypropylene_mu2'):
                mask[idx[lbl]] = False
            # Run-5 mirror: the promoted chain-scale radical classifies
            # CONDENSED in the core mask (melt window) -- what kept the row
            # open through the product gates and mass-paired past the thermo
            # reference-state tripwire while its legacy flux drained the pool.
            mask[idx[r2_spc.label]] = False
            pool = self._pool('polypropylene',
                              (idx['polypropylene_mu0'],
                               idx['polypropylene_mu1'],
                               idx['polypropylene_mu2']))
            rs = HybridPolymerSystem(
                T=800.0, P=1.0e5,
                initial_mole_fractions={core_species[idx[r1_spc.label]]: 0.0},
                V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
                gas_species_mask=mask, constant_gas_volume=False,
                initial_polymer_moments={'polypropylene': (1.0, 50.0, 3000.0)},
                termination=[], allow_default_prospective_edge=True)
            rs.initialize_model(core_species, list(cerm.core.reactions),
                                list(cerm.edge.species),
                                list(cerm.edge.reactions))

            i = cerm.core.reactions.index(out)
            # THE RED asserts (FIX 2 -- rebuild restamping at the last honest
            # chokepoint): the canonical core row must arrive adjudicated.
            assert getattr(out, 'polymer_refused', False) is True, (
                "adjudication lost across promotion: canonical core row "
                "arrived unstamped at the rebuild")
            assert rs.reaction_refused[i] == 1, (
                "reaction_refused missed the refused core row")
            assert rs.reaction_flux_archetype[i] == solver_mod.FLUX_NONE, (
                "refused row was demoted onto the legacy mu1-only path")
            assert any(c['reason'] == 'conduit-deferred'
                       for c in rs.refused_reaction_census)
            dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            assert np.all(dn == 0.0), (
                "refused core row applied nonzero RHS flux (the run-5 stall "
                "mechanism: live legacy mu1 flux against the pool)")
        finally:
            rmg_data.database = old_db

    def test_red_hard_fail_on_unclassified_live_proxy_row(self):
        """RED (FIX 4): a proxy-touching row that is neither refused nor
        archetype-classified after restamping must HARD-FAIL initialize_model
        with an error naming the row -- never silently run legacy mu1-only
        flux (the run-5 stall channel)."""
        sp, core, mask = _refstate_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]],
                       **_KIN)   # irreversible: archetype NONE, not refused
        # LIVENESS PIN: the row is genuinely proxy-touching and unclassified.
        assert int(getattr(rxn, "polymer_flux_archetype", 0)) == 0
        assert not getattr(rxn, "polymer_refused", False)
        # Constructed DIRECTLY (not via _refstate_rs, whose direct-test
        # posture sets the escape hatch): this pins the PRODUCTION default,
        # allow_unstamped_proxy_rows=False -> hard-fail.
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[-1]: 0.0},
            V_poly=1.0, polymer_pools=[_gate_pool_config()], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[])
        with pytest.raises(ValueError, match="polymer_flux_archetype"):
            rs.initialize_model(core, [rxn], [], [])

    def test_unclassified_proxy_row_escape_hatch_keeps_legacy_demotion(self):
        """FIX 4 escape hatch (documented test/runner-only, mirroring
        allow_default_prospective_edge): with allow_unstamped_proxy_rows=True
        the pre-fix legacy demotion (NONE -> UNRESOLVED, warn-once) is
        preserved for direct-construction fixtures."""
        import rmgpy.solver.polymer as solver_mod
        sp, core, mask = _refstate_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]],
                       **_KIN)
        rs = _refstate_rs(core, [rxn], mask, [_gate_pool_config()],
                          {"A": (1.0, 5.0, 30.0)},
                          rs_kwargs={"allow_unstamped_proxy_rows": True})
        assert rs.reaction_flux_archetype[0] == solver_mod.FLUX_UNRESOLVED

    def test_negative_control_live_ve_row_with_real_polymer_keeps_flux(self):
        """Negative control (r71 mandate): a legitimately LIVE proxy row
        carried by a REAL Polymer object (the restamp predicate's phase
        vocabulary) must NOT be refused by the rebuild restamp, and its
        residual stays bitwise at the pre-change value pinned below."""
        from rmgpy.polymer import Polymer
        import rmgpy.solver.polymer as solver_mod
        pp = Polymer(label='poly', monomer='[CH2][CH]C',
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        mu0 = _spc("CO", "poly_mu0"); mu1 = _spc("C=O", "poly_mu1")
        mu2 = _spc("C#N", "poly_mu2"); vol = _spc("C=CC", "propene")
        for s in (pp, mu0, mu1, mu2, vol):
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        core = [pp, mu0, mu1, mu2, vol]
        mask = np.array([False] * 4 + [True], dtype=bool)
        ve = Reaction(reactants=[pp], products=[pp, vol], **_KIN)
        ve.polymer_flux_archetype = 6   # VOLATILE_EJECTION
        ve.polymer_eject_units = 1.0
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={vol: 0.0}, V_poly=1.0,
            polymer_pools=[self._pool("poly", (1, 2, 3))], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 50.0, 3000.0)},
            termination=[])
        rs.initialize_model(core, [ve], [], [])
        assert rs.reaction_refused[0] == 0, (
            "rebuild restamp refused a legitimately live VE row")
        assert rs.reaction_flux_archetype[0] == solver_mod.FLUX_VOLATILE_EJECTION
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        # Bitwise pre-change pin (recorded on b0de7dde8 BEFORE the fix):
        assert dn[4] == _VE_CONTROL_PIN[4] != 0.0
        assert np.array_equal(dn[:5], _VE_CONTROL_PIN)


class TestSameProxyRefusalFluxDead:
    """r74 adjudication (PP run-5 "H fountain", forensics
    /home/alon/Projects/polymer/PP/rmg/run5 +
    /home/alon/Projects/polymer/PP/ta-diag-run5): the two degenerate core
    rows ``H(1) + rad_end <=> rad_end`` carry the SAME resolved pool proxy
    Species on BOTH sides, so Keq collapses to G(H) alone and the reverse
    leg is a unimolecular H source (kb ~ 1e16 in this fixture) that strips
    the pool's mu1 at one H-mass per event -- the TA diagnostic volatilized
    the whole 10 mg PP sample as H2 at 889 C, element-impossible. The
    same-proxy refusal must ride the FULL r71 refusal contract: re-derived
    at the rebuild restamp (FIX 2) regardless of upstream stamp loss, zero
    core RHS flux, zero edge enlargement inputs (FIX 3), archetype ignored
    while refused.

    RED-FIRST: every red assertion below was quoted FAILING on pre-fix HEAD
    (9a7530b8b) -- pre-fix the core row produces H at +8.1e17 mol/s while
    draining pool mu1/mu2, and the edge row leaks -8.1e17 into
    edge_reaction_rates."""

    @staticmethod
    def _degenerate_fixture():
        """Run-5 degenerate row built through the REAL generation route
        (make_new_reaction handshake + fold-back + species_dict resolution),
        assembled into a core species list with pool config -- the same
        recipe as TestRefusalAdjudicationSurvivesRebuild's promotion half."""
        from rmgpy.data.kinetics import TemplateReaction
        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer
        from rmgpy.rmg.model import CoreEdgeReactionModel

        cerm = CoreEdgeReactionModel()
        pp = Polymer(label='polypropylene', monomer='[CH2][CH]C',
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        cerm._register_polymer(pp, generate_thermo=False)
        prim, sec = pp.generate_end_radical_daughters()
        cerm._register_polymer(sec, generate_thermo=False)
        rad = sec.molecule[0].copy(deep=True)
        capped = rad.copy(deep=True)
        capped.saturate_radicals()
        forward = TemplateReaction(
            reactants=[Molecule().from_smiles('[H]'), rad],
            products=[capped],
            family='R_Recombination', is_forward=True, reversible=True,
            kinetics=Arrhenius(A=(1.0e7, 'm^3/(mol*s)'), n=0.0,
                               Ea=(0.0, 'kJ/mol')))
        out, is_new = cerm.make_new_reaction(
            forward, check_existing=False, generate_thermo=False,
            generate_kinetics=False)
        assert out is not None and is_new
        # LIVENESS: the run-5 laundering mechanism itself -- the H-capped
        # adduct resolved back onto the SAME rad_end pool object.
        polymer_p = [s for s in out.products if isinstance(s, Polymer)]
        assert polymer_p and polymer_p[0] is sec, (
            "FIXTURE BROKEN: product did not fold back onto the rad_end pool")
        h_spc = next(s for s in out.reactants if not isinstance(s, Polymer))

        cerm.add_species_to_edge(sec)
        cerm.add_species_to_core(sec)     # brings the mu dummies along
        cerm.add_species_to_core(h_spc)
        core_species = list(cerm.core.species)
        for s in core_species:
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        # Polymer thermo delegates to its ACTIVE proxy Species
        # (get_free_energy -> get_proxy_species().thermo -- the end-radical
        # daughter's reactive proxy, not baseline): the reversible row's Keq
        # needs it resolvable without a thermo database.
        sec.get_proxy_species().thermo = _trivial_nasa(_GAV_COMMENT)
        sec.baseline_proxy.thermo = _trivial_nasa(_GAV_COMMENT)
        idx = {s.label: i for i, s in enumerate(core_species)}
        mask = np.ones(len(core_species), dtype=bool)
        for lbl in (sec.label, sec.label + '_mu0', sec.label + '_mu1',
                    sec.label + '_mu2'):
            mask[idx[lbl]] = False
        pool = PolymerPoolConfig(
            label=sec.label, xs=2, explicit_dp_to_species_index={},
            mu_indices=(idx[sec.label + '_mu0'], idx[sec.label + '_mu1'],
                        idx[sec.label + '_mu2']),
            monomer_poly_index=None, monomer_mw_g_mol=42.08,
            k_scission=0.0, k_unzip=0.0)
        return out, h_spc, core_species, mask, pool

    def _build_rs(self, out, h_spc, core_species, mask, pool,
                  core_rxns, edge_rxns):
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={h_spc: 0.0},
            V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={pool.label: (1.0, 50.0, 3000.0)},
            termination=[], allow_default_prospective_edge=True)
        # The adjudicated run-5 arrival state (r71, binding): by rebuild
        # time the row can arrive UNSTAMPED (stamps are ad-hoc attributes,
        # losable across canonical dedup / __reduce__) -- the rebuild
        # restamp must re-derive the refusal regardless.
        out.polymer_refused = False
        out.polymer_refused_accumulating = False
        rs.initialize_model(core_species, core_rxns, [], edge_rxns)
        return rs

    def _with_db(self, body):
        import rmgpy.data.rmg as rmg_data
        from rmgpy.data.base import ForbiddenStructures
        from rmgpy.data.rmg import RMGDatabase

        old_db = rmg_data.database
        db = RMGDatabase()
        db.forbidden_structures = ForbiddenStructures()
        rmg_data.database = db
        try:
            return body()
        finally:
            rmg_data.database = old_db

    def test_red_same_proxy_core_row_restamped_refused_and_rhs_dead(self):
        """RED (restamp + core RHS): the rebuild restamp must re-derive the
        same-proxy refusal on the core row, and the refused row must apply
        ZERO RHS flux. Pre-fix: reaction_refused[0] == 0 and the residual
        produces H at +8.1e17 mol/s while draining pool mu1/mu2 -- the H
        fountain running inside the RMG solver."""
        def body():
            out, h_spc, core_species, mask, pool = self._degenerate_fixture()
            rs = self._build_rs(out, h_spc, core_species, mask, pool,
                                [out], [])
            # LIVENESS PIN -- BEFORE the red asserts: the reverse leg
            # genuinely carries rate (kb = kf/Keq with Keq = f(G_H) alone),
            # so the zeros below cannot mean "fixture dead".
            assert rs.kb[0] > 0.0, (
                "FIXTURE BROKEN, not a valid red: the degenerate row "
                "carries no reverse rate coefficient at all")
            # THE red asserts (r74 riding r71 FIX 2 + item-18 flux gate):
            assert getattr(out, 'polymer_refused', False) is True, (
                "rebuild restamp did not re-derive the same-proxy refusal "
                "(the run-5 arrival state runs live)")
            assert rs.reaction_refused[0] == 1, (
                "reaction_refused missed the same-proxy degenerate core row")
            dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            assert np.all(dn == 0.0), (
                "same-proxy-refused core row applied nonzero RHS flux -- "
                "the H fountain (unimolecular H source, Keq = f(G_H) alone) "
                "is live in the solver")
        self._with_db(body)

    def test_red_same_proxy_edge_row_flux_dead_for_enlargement(self):
        """RED (edge half, r71 FIX 3 inheritance): the same-proxy refusal
        must zero the refused EDGE row's enlargement inputs
        (edge_reaction_rates); only the ungated per-reaction counterfactual
        survives (diagnostic, never flux). Pre-fix the row leaks -8.1e17
        into edge_reaction_rates."""
        def body():
            out, h_spc, core_species, mask, pool = self._degenerate_fixture()
            rs = self._build_rs(out, h_spc, core_species, mask, pool,
                                [], [out])
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
            # LIVENESS PIN -- BEFORE the red asserts: the row genuinely
            # carries reverse flux; the ungated counterfactual sees it, so
            # the zero below cannot mean "fixture dead".
            assert rs.edge_reaction_rates_ungated[0] != 0.0, (
                "FIXTURE BROKEN, not a valid red: the degenerate edge row "
                "carries no flux at all")
            # THE red asserts:
            assert rs.reaction_refused[0] == 1, (
                "restamp missed the same-proxy degenerate EDGE row")
            assert rs.edge_reaction_rates[0] == 0.0, (
                "same-proxy-refused edge row leaked into "
                "edge_reaction_rates (the enlargement input that promoted "
                "run-5's rows to core)")
        self._with_db(body)


class TestImpostorRowRefusalFluxDead:
    """r82 adjudication (FR1 run-2 impostor rows, forensics
    /home/alon/Projects/polymer/FR1/rmg/run2/RMG.out ~line 1774): 15
    XY_Addition_MultipleBond rows shaped
    ``Br/BrBr + <C36-scale unsaturated discrete> <=> FR1`` bridge a
    POLYMER-SIZED discrete molecule (the pool proxy minus a small
    closed-shell gas partner) into the condensed proxy. The gas partner is
    closed-shell, so the r63 all-gas-radicals conjunct never fires; the
    rows arrived at the solver rebuild UNSTAMPED and killed the run at the
    r71 unstamped-proxy hard-fail. The r82 impostor refusal must ride the
    FULL r71 refusal contract: stamped at generation, re-derived at the
    rebuild restamp regardless of upstream stamp loss, zero core RHS flux,
    zero edge enlargement inputs, conduit-deferred census reason.

    RED-FIRST: pre-fix (3521caead) the fixture row is NOT refused at
    generation and ``initialize_model`` raises the r71
    ``UNSTAMPED PROXY ROW(S)`` ValueError -- the exact run-2 death."""

    @staticmethod
    def _impostor_fixture():
        """Run-2-shaped impostor row built through the REAL generation route
        (make_new_reaction + species_dict isomorphism fold-back onto the
        registered pool): H2 + <proxy-with-one-double-bond> <=> PS pool.
        Uses a POLYSTYRENE pool because its stitched trimer proxy is itself
        genuinely chain-scale (C24H26, 314.5 g/mol / 24 heavy, clearing the r95
        absolute floor) -- a PP proxy (C9, 128 g/mol) is below the floor and
        would make an honest chain-scale impostor unrepresentable. The
        unsaturated C24H24 discrete (312.5 g/mol / 24 heavy) is the
        proxy-minus-H2 impostor, mirroring FR1 run-2's C36-scale
        proxy-minus-Br2/HBr molecules that folded onto the condensed proxy."""
        from rmgpy.data.kinetics import TemplateReaction
        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer
        from rmgpy.rmg.model import CoreEdgeReactionModel

        cerm = CoreEdgeReactionModel()
        pp = Polymer(label='PS', monomer='[CH2][CH]c1ccccc1',
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        cerm._register_polymer(pp, generate_thermo=False)
        forward = TemplateReaction(
            reactants=[Molecule().from_smiles('[H][H]'),
                       Molecule().from_smiles(
                           'C=C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1')],
            products=[Molecule().from_smiles(
                'CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1')],
            family='XY_Addition_MultipleBond', is_forward=True,
            reversible=True,
            kinetics=Arrhenius(A=(1.0e7, 'm^3/(mol*s)'), n=0.0,
                               Ea=(0.0, 'kJ/mol')))
        out, is_new = cerm.make_new_reaction(
            forward, check_existing=False, generate_thermo=False,
            generate_kinetics=False)
        assert out is not None and is_new
        # LIVENESS: the run-2 impostor mechanism itself -- the saturated
        # addition product resolved onto the registered pool Polymer via
        # species_dict isomorphism, while the unsaturated impostor stayed a
        # discrete gas Species.
        polymer_p = [s for s in out.products if isinstance(s, Polymer)]
        assert polymer_p and polymer_p[0] is pp, (
            "FIXTURE BROKEN: addition product did not fold onto the pool")
        assert not any(isinstance(s, Polymer) for s in out.reactants), (
            "FIXTURE BROKEN: impostor side must carry no Polymer")
        gas_spc = next(s for s in out.reactants
                       if s.molecule[0].get_formula() == 'H2')
        imp_spc = next(s for s in out.reactants
                       if s.molecule[0].get_formula() == 'C24H24')

        cerm.add_species_to_edge(pp)
        cerm.add_species_to_core(pp)      # brings the mu dummies along
        cerm.add_species_to_core(gas_spc)
        cerm.add_species_to_core(imp_spc)
        core_species = list(cerm.core.species)
        for s in core_species:
            s.thermo = _trivial_nasa(_GAV_COMMENT)
        # The reversible row's Keq needs the pool's proxy thermo resolvable
        # without a thermo database.
        pp.baseline_proxy.thermo = _trivial_nasa(_GAV_COMMENT)
        idx = {s.label: i for i, s in enumerate(core_species)}
        mask = np.ones(len(core_species), dtype=bool)
        for lbl in ('PS', 'PS_mu0', 'PS_mu1', 'PS_mu2'):
            mask[idx[lbl]] = False
        pool = PolymerPoolConfig(
            label='PS', xs=2, explicit_dp_to_species_index={},
            mu_indices=(idx['PS_mu0'], idx['PS_mu1'], idx['PS_mu2']),
            monomer_poly_index=None, monomer_mw_g_mol=104.15,
            k_scission=0.0, k_unzip=0.0)
        return out, gas_spc, core_species, mask, pool

    def _build_rs(self, out, gas_spc, core_species, mask, pool,
                  core_rxns, edge_rxns):
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={gas_spc: 0.0},
            V_poly=1.0, polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={pool.label: (1.0, 50.0, 3000.0)},
            termination=[], allow_default_prospective_edge=True)
        # The adjudicated run-2 arrival state (r71, binding): by rebuild
        # time the row can arrive UNSTAMPED (stamps are ad-hoc attributes,
        # losable across canonical dedup / __reduce__) -- the rebuild
        # restamp must re-derive the refusal regardless.
        out.polymer_refused = False
        out.polymer_refused_accumulating = False
        rs.initialize_model(core_species, core_rxns, [], edge_rxns)
        return rs

    def _with_db(self, body):
        import rmgpy.data.rmg as rmg_data
        from rmgpy.data.base import ForbiddenStructures
        from rmgpy.data.rmg import RMGDatabase

        old_db = rmg_data.database
        db = RMGDatabase()
        db.forbidden_structures = ForbiddenStructures()
        rmg_data.database = db
        try:
            return body()
        finally:
            rmg_data.database = old_db

    def test_red_impostor_row_refused_at_generation(self):
        """RED (generation half): make_new_reaction must stamp the impostor
        row refused conduit-deferred at classification time. Pre-fix
        (3521caead): no classifier sees the shape (the gas partner is
        closed-shell, so the r63 conjunct passes it) and polymer_refused
        stays False -- run-2's 15 rows entered the model live."""
        def body():
            out, _, _, _, _ = self._impostor_fixture()
            assert getattr(out, 'polymer_refused', False) is True, (
                "impostor row not refused at generation -- the run-2 shape "
                "enters the model unadjudicated")
            assert getattr(out, 'polymer_refused_accumulating', True) is False
        self._with_db(body)

    def test_red_impostor_core_row_restamped_refused_and_rhs_dead(self):
        """RED (rebuild restamp + core RHS): the rebuild restamp must
        re-derive the impostor refusal on an UNSTAMPED-arriving core row,
        and the refused row must apply ZERO RHS flux. Pre-fix (3521caead):
        initialize_model raises the r71 ``UNSTAMPED PROXY ROW(S)``
        ValueError -- the exact FR1 run-2 death at RMG.out ~line 1774."""
        def body():
            import rmgpy.solver.polymer as solver_mod
            out, gas_spc, core_species, mask, pool = self._impostor_fixture()
            rs = self._build_rs(out, gas_spc, core_species, mask, pool,
                                [out], [])
            # LIVENESS PIN -- BEFORE the red asserts: the row genuinely
            # carries a forward rate coefficient, so the zeros below cannot
            # mean "fixture dead".
            assert rs.kf[0] > 0.0, (
                "FIXTURE BROKEN, not a valid red: the impostor row carries "
                "no forward rate coefficient at all")
            # THE red asserts (r82 riding r71 FIX 2 + item-18 flux gate):
            assert getattr(out, 'polymer_refused', False) is True, (
                "rebuild restamp did not re-derive the impostor refusal "
                "(the run-2 arrival state runs live)")
            assert rs.reaction_refused[0] == 1, (
                "reaction_refused missed the impostor core row")
            assert rs.reaction_flux_archetype[0] == solver_mod.FLUX_NONE, (
                "refused impostor row was demoted onto the legacy "
                "mu1-only path")
            assert any(c['reason'] == 'conduit-deferred'
                       for c in rs.refused_reaction_census)
            dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
            assert np.all(dn == 0.0), (
                "impostor-refused core row applied nonzero RHS flux")
        self._with_db(body)

    def test_red_impostor_edge_row_flux_dead_for_enlargement(self):
        """RED (edge half, r71 FIX 3 inheritance): the impostor refusal must
        zero the refused EDGE row's enlargement inputs
        (edge_reaction_rates); only the ungated per-reaction counterfactual
        survives (diagnostic, never flux). Pre-fix (3521caead):
        initialize_model raises the r71 hard-fail on the unstamped edge row
        (run-2's rows were edge rows 53-67)."""
        def body():
            out, gas_spc, core_species, mask, pool = self._impostor_fixture()
            rs = self._build_rs(out, gas_spc, core_species, mask, pool,
                                [], [out])
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
            # LIVENESS PIN -- BEFORE the red asserts: the row genuinely
            # carries flux; the ungated counterfactual sees it, so the zero
            # below cannot mean "fixture dead".
            assert rs.edge_reaction_rates_ungated[0] != 0.0, (
                "FIXTURE BROKEN, not a valid red: the impostor edge row "
                "carries no flux at all")
            # THE red asserts:
            assert rs.reaction_refused[0] == 1, (
                "restamp missed the impostor EDGE row")
            assert rs.edge_reaction_rates[0] == 0.0, (
                "impostor-refused edge row leaked into edge_reaction_rates "
                "(the enlargement input)")
        self._with_db(body)


class TestSpawnGateDefectAwareMass:
    """FR1-K2 mass-consumer audit (round-72 P2): the spawn-gate snapshot's
    per-chain event mass must use the EXACT condensed mass
    (condensed_mass_g(1, E[n]) = E[n]*monomer_MW - chain_mass_defect), not
    the raw E[n]*MW -- an X-loss feature pool's chains each lost one X, so
    the moment-derived per-chain mass overstates by M_X. pool_stats grows
    the pool's defect as a third element so the python-side ledger half
    (_snapshot_event_mass) stays consistent with the engine half."""

    def _defect_rs(self):
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]],
                       **_KIN)
        rxn.polymer_flux_archetype = 1
        # E[n] = 5, MW = 28, defect = 8 -> exact per-chain mass 132, not
        # 140.
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0),
                               monomer_mw_g_mol=28.0,
                               chain_mass_defect_g_mol=8.0)
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        return rs

    def test_pool_stats_carry_defect_and_proxy_total_is_exact(self):
        """RED pin: pool_stats rows are (E[n], MW, defect) triples and
        proxy_event_mass_total uses the exact per-chain mass
        E[n]*MW - defect for the defect pool."""
        rs = self._defect_rs()
        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()
        e_n, mw, defect = pool_stats["A"]
        assert e_n == pytest.approx(5.0)
        assert mw == pytest.approx(28.0)
        assert defect == pytest.approx(8.0)
        assert proxy_total == pytest.approx(
            gross["A"] * (5.0 * 28.0 - 8.0), rel=1e-12)
        assert proxy_total > 0.0

    def test_defect_free_pool_unchanged(self):
        """Negative control: an ordinary pool (defect 0) keeps the exact
        legacy E[n]*MW event mass (the accessor reduces to mu1*MW)."""
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]],
                       **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0),
                               monomer_mw_g_mol=28.0)
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()
        e_n, mw, defect = pool_stats["A"]
        assert defect == 0.0
        assert proxy_total == pytest.approx(gross["A"] * 5.0 * 28.0,
                                            rel=1e-12)


# ---------------------------------------------------------------------------
# r81 run-8 gate (adjudicated round 81, 2026-07-06): (B) exhaustion-tail
# conditioning + (C) resurrection zero-metric guard. Forensics: PP run-7
# (/home/alon/Projects/polymer/PP/rmg/run7/RMG.log) -- after full
# depropagation every pool sits at 1e-15..1e-18 mol, DASPK hits IDID=-7 on
# the ill-conditioned exhausted tail, and Model Resurrection amplifies it by
# adding species at rate ratio 0.0 (zero flux), six loops deep.
# ---------------------------------------------------------------------------


def _r81_pool_rs(moments, k_unzip=0.0, monomer_poly_index=None):
    """One-pool fixture for the exhaustion-conditioning pins: pool A
    (proxy + mu dummies), gas seed X, and one core fold-back row
    A -> A + X (the canonical unstamped legacy shape, rate = kf*mu1).
    Default solver atol = 1e-16, so the r81 exhaustion floor is
    max(SMALL_EPS, 100*atol) = 1e-14 mol per moment slot."""
    sp = _gate17_species()
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
    mask = [False, False, False, False, True]
    rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["X"]], **_KIN)
    rs = _gate17_rs(core, mask, [rxn], moments={"A": moments},
                    k_unzip=k_unzip, monomer_poly_index=monomer_poly_index)
    return rs, core, rxn


class TestR81ExhaustionTailConditioning:
    """(B) The centralized pool-exhaustion predicate, re-adjudicated r86
    (Codex-adjudicated, spar r86): the predicate is a CENSUS/TRIPWIRE
    driver ONLY. A pool is counted exhausted when |mu0|, |mu1| AND |mu2|
    all sit below per-state floors max(SMALL_EPS, 100*atol[state]) --
    census once per pool per rebuild -- but the RHS keeps computing from
    the RAW state vector (the r81 read-projection is reverted; flags have
    ZERO RHS effect). A moment more negative than -floor stays a HARD
    error, never max(..., SMALL_EPS)."""

    # all three moments below the 1e-14 mol floor -- the run-7 tail state
    _EXHAUSTED = (1.0e-16, 5.0e-16, 3.0e-15)

    def test_exhausted_pool_rhs_stays_live_census_only(self, caplog):
        """r86 re-adjudication of the r81 zero-RHS pin: with every moment
        in the exhaustion band the census fires (once, diagnostic), but
        the RHS computes LIVE from raw y -- the proxy row keeps its exact
        legacy rate kf*mu1 = 2.0*5e-16 and the residual is NOT the zero
        vector (the r81 projected-to-zero read is reverted; the run-9d
        H-early wake-up must never see a frozen pool)."""
        rs, _, _ = _r81_pool_rs(self._EXHAUSTED, k_unzip=1.0,
                                monomer_poly_index=4)
        with caplog.at_level(logging.WARNING):
            delta, _ = rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert float(np.asarray(rs.core_reaction_rates)[0]) == \
            pytest.approx(2.0 * self._EXHAUSTED[1], rel=1e-12)
        assert float(np.max(np.abs(np.asarray(delta)))) > 0.0
        assert any("POOL EXHAUSTION CENSUS" in r.getMessage()
                   for r in caplog.records)

    def test_floor_crossing_pool_wake_up_integrates_smoothly(self):
        """r86 revert demonstration (the run-9d H-early wake-up shape in
        miniature): a pool whose moments START inside the census band is
        fed by a healthy pool and must DASSL-integrate smoothly THROUGH
        the floor crossing -- finite state throughout, monotone growth of
        the waking pool, and its own mu1-scaled fold-back row live once
        awake (gas actually produced). The RED half of the revert (the
        exactly-zero in-band projection removed) is pinned at residual
        level by test_exhausted_pool_rhs_stays_live_census_only and the
        two cross-pool continuity tests; an integration-scale RED
        observable is NOT constructible here because the pre-crossing gas
        signal sits below the integrator atol (1e-16) by construction of
        the census band -- a trajectory-level pre/post difference would
        pin integrator noise, not physics."""
        sp, core, mask = _two_pool_species()
        feed = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        feed.polymer_flux_archetype = 2          # MIGRATION: A chains -> B
        foldback = Reaction(reactants=[sp["B"]],
                            products=[sp["B"], sp["G"]], **_KIN)
        pool_a = PolymerPoolConfig(label="A", xs=2,
                                   explicit_dp_to_species_index={},
                                   mu_indices=(1, 2, 3),
                                   monomer_poly_index=None,
                                   k_scission=0.0, k_unzip=0.0,
                                   tail_kinetics=None)
        pool_b = PolymerPoolConfig(label="B", xs=2,
                                   explicit_dp_to_species_index={},
                                   mu_indices=(5, 6, 7),
                                   monomer_poly_index=None,
                                   k_scission=0.0, k_unzip=0.0,
                                   tail_kinetics=None)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0},
            V_poly=1.0, polymer_pools=[pool_a, pool_b], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0),
                                     "B": (1.9e-18, 2.6e-17, 6.1e-16)},
            termination=[],
            allow_default_prospective_edge=True,
            allow_unstamped_proxy_rows=True)
        rs.initialize_model(core, [feed, foldback], [], [])
        mu1_b_traj = []
        for t in (1.0e-3, 1.0e-2, 1.0e-1, 1.0):
            rs.advance(t)
            y = np.asarray(rs.y)
            assert np.all(np.isfinite(y)), t
            mu1_b_traj.append(float(y[6]))
        # B woke up: mu1_B grew monotonically from 2.6e-17 mol (in-band)
        # to far above the 1e-14 mol census floor -- no wall, no cliff.
        assert all(b > a for a, b in zip(mu1_b_traj, mu1_b_traj[1:]))
        assert mu1_b_traj[-1] > 1.0e-3
        # ... and B's own fold-back row is LIVE post-wake: gas accumulated.
        assert float(np.asarray(rs.y)[8]) > 0.0

    def test_tiny_mu0_with_real_mu1_is_not_clamped(self, caplog):
        """Anti-overreach pin (r81 B1: never mu0 alone): tiny mu0 with a
        nontrivial mu1 is a few LONG chains, not an empty pool -- the
        mu1-scaled proxy row keeps its exact legacy rate kf*mu1 = 10.0 and
        NO exhaustion census fires."""
        rs, _, _ = _r81_pool_rs((1.0e-16, 5.0, 30.0))
        with caplog.at_level(logging.WARNING):
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert float(np.asarray(rs.core_reaction_rates)[0]) == \
            pytest.approx(10.0, rel=1e-12)
        assert not any("POOL EXHAUSTION CENSUS" in r.getMessage()
                       for r in caplog.records)

    def test_negative_moment_beyond_floor_hard_errors(self):
        """RED (r81 negative-moment rule): a moment more negative than
        -floor is integrator corruption -- the residual must raise loudly,
        never max(..., SMALL_EPS) it back to zero."""
        rs, _, _ = _r81_pool_rs((1.0, 5.0, 30.0))
        y2 = np.array(rs.y, dtype=float).copy()
        y2[rs.polymer_pools[0].mu_indices[0]] = -1.0e-12  # << -1e-14 floor
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            rs.residual(0.0, y2, np.zeros_like(y2))

    def test_negative_moment_within_band_projects_with_census(self, caplog):
        """r81 negative-moment rule, band half (kept under r86 as a
        census-only pin): a negative inside [-floor, +floor] is ANNOUNCED
        through the exhaustion census -- never silently -- and does NOT
        raise. The RHS keeps its pre-existing local max(0, .) reads; the
        census is diagnostic only (r86: no RHS projection)."""
        rs, _, _ = _r81_pool_rs((1.0, 5.0, 30.0))
        y2 = np.array(rs.y, dtype=float).copy()
        y2[rs.polymer_pools[0].mu_indices[0]] = -5.0e-15  # inside the band
        with caplog.at_level(logging.WARNING):
            rs.residual(0.0, y2, np.zeros_like(y2))
        assert any("POOL EXHAUSTION CENSUS" in r.getMessage()
                   for r in caplog.records)

    def test_exhaustion_census_once_per_pool_per_rebuild(self, caplog):
        """RED (r81 B3): the census line carries label + raw moments +
        floors, fires ONCE per pool per rebuild (two residual calls -> one
        line), and re-fires after an initialize_model rebuild."""
        rs, core, rxn = _r81_pool_rs(self._EXHAUSTED)
        with caplog.at_level(logging.WARNING):
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        lines = [r.getMessage() for r in caplog.records
                 if "POOL EXHAUSTION CENSUS" in r.getMessage()]
        assert len(lines) == 1
        assert "pool A" in lines[0]
        assert "floor" in lines[0]
        caplog.clear()
        rs.initialize_model(list(core), [rxn], [], [])
        with caplog.at_level(logging.WARNING):
            rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        assert any("POOL EXHAUSTION CENSUS" in r.getMessage()
                   for r in caplog.records)

    def test_spawn_gate_negative_beyond_floor_hard_errors(self):
        """RED (r81 B5): the spawn-gate mu0 read at polymer.pyx:3568
        (mu0 <= SMALL_EPS -> 0.0 silently) swallowed negatives; beyond
        -floor it must now raise."""
        rs, _, _ = _r81_pool_rs((1.0, 5.0, 30.0))
        rs.y[int(rs.pool_mu0_indices[0])] = -1.0e-12
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            rs.spawn_gate_flux_snapshot()

    def test_spawn_gate_in_band_zero_carries_census(self, caplog):
        """RED (r81 B5, band half): a spawn-gate mu0 inside
        [-floor, SMALL_EPS] still defers (E[n] = 0) but now announces the
        projection through the exhaustion census instead of silence."""
        rs, _, _ = _r81_pool_rs((1.0, 5.0, 30.0))
        rs.y[int(rs.pool_mu0_indices[0])] = -5.0e-15
        with caplog.at_level(logging.WARNING):
            _, pool_stats, _ = rs.spawn_gate_flux_snapshot()
        assert pool_stats["A"][0] == 0.0
        assert any("POOL EXHAUSTION CENSUS" in r.getMessage()
                   for r in caplog.records)


class _R81FailingStepSystem(HybridPolymerSystem):
    """step() reproduces the run-7 DASPK wall unconditionally: every
    attempted integration step dies with the IDID=-7 convergence failure."""

    def step(self, step_time):
        from pydas.dassl import DASSLError
        raise DASSLError(
            "DASSL returned with an IDID = -7, Repeated convergence test "
            "failures occurred on the last attempted step in DDASSL.")


def _r81_resurrection_rs(edge_reactant_key):
    """Resurrection-guard fixture: pool A, gas seed X, gas driver X -> Y
    (positive char_rate), and ONE edge species G fed by
    sp[edge_reactant_key] -> G. edge_reactant_key='Y' (zero moles) gives the
    exact run-7 shape -- edge rate ratio 0.0; 'X' (the 1-mol seed) gives a
    strictly positive edge rate."""
    sp = _gate17_species()
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"], sp["Y"]]
    mask = np.array([False, False, False, False, True, True], dtype=bool)
    driver = Reaction(
        reactants=[sp["X"]], products=[sp["Y"]],
        kinetics=Arrhenius(A=(1.0e-3, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"),
                           T0=(298.15, "K")),
        reversible=False)
    gated = Reaction(reactants=[sp[edge_reactant_key]], products=[sp["G"]],
                     **_KIN)
    cfg = PolymerPoolConfig(label="A", xs=2,
                            explicit_dp_to_species_index={},
                            mu_indices=(1, 2, 3), monomer_poly_index=None)
    rs = _R81FailingStepSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
        V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
        gas_species_mask=mask, constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
        allow_default_prospective_edge=True,
        allow_unstamped_proxy_rows=True)
    return rs, core, [driver], [sp["G"]], [gated]


def _r81_simulate(rs, core, rxns_core, edge_spcs, rxns_edge, ms):
    from rmgpy.rmg.settings import SimulatorSettings
    from rmgpy.solver.base import TerminationTime
    rs.termination.append(TerminationTime((1.0, "s")))
    return rs.simulate(list(core), list(rxns_core), list(edge_spcs),
                       list(rxns_edge), [], [],
                       model_settings=ms,
                       simulator_settings=SimulatorSettings())


class TestR81ResurrectionZeroMetricGuard:
    """(C) Model Resurrection must never promote an object whose selected
    metric is <= 0 or below the normal movement threshold; when every
    candidate is zero it logs 'no positive resurrection candidate' and stops
    resurrecting -- the DASPK failure is surfaced honestly instead of
    growing the core on zero flux (the run-7 amplifier)."""

    def test_run7_shape_zero_edge_rates_do_not_add_species(self, caplog):
        """RED (the exact run-7 shape): DASPK error + ALL edge rates zero
        must NOT add a species at rate ratio 0.0 -- resurrection stops and
        the failure surfaces as the loud resurrection-failed error."""
        from rmgpy.rmg.settings import ModelSettings
        rs, core, rxns_core, edge_spcs, rxns_edge = _r81_resurrection_rs("Y")
        ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e-3,
                           tol_interrupt_simulation=1.0e8)
        with caplog.at_level(logging.INFO):
            with pytest.raises(ValueError,
                               match="invalid_objects could not be filled"):
                _r81_simulate(rs, core, rxns_core, edge_spcs, rxns_edge, ms)
        msgs = [r.getMessage() for r in caplog.records]
        assert not any(
            "was added to model core in model resurrection process" in m
            for m in msgs)
        assert any("no positive resurrection candidate" in m for m in msgs)

    def test_positive_edge_rate_resurrection_unchanged(self, caplog):
        """No-regression pin (r81 C3): with a strictly positive edge rate
        at/above the movement threshold, resurrection still promotes the
        candidate and simulate returns resurrected=True."""
        from rmgpy.rmg.settings import ModelSettings
        rs, core, rxns_core, edge_spcs, rxns_edge = _r81_resurrection_rs("X")
        # ignore_overall_flux_criterion keeps normal enlargement from
        # collecting G first, so the DASPK failure exercises the
        # resurrection selection itself; tol_move_to_core=0 makes any
        # strictly positive ratio qualify.
        ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=0.0,
                           tol_interrupt_simulation=1.0e8,
                           ignore_overall_flux_criterion=True)
        with caplog.at_level(logging.INFO):
            terminated, resurrected, invalid_objects, _, _, _, _ = \
                _r81_simulate(rs, core, rxns_core, edge_spcs, rxns_edge, ms)
        assert resurrected is True
        assert edge_spcs[0] in invalid_objects
        assert any(
            "was added to model core in model resurrection process"
            in r.getMessage() for r in caplog.records)


class _StiffPastCrossingSystem(HybridPolymerSystem):
    """Bounded stand-in for the FR1 tf100 signature (r86): the
    post-conversion region is a stiff wall -- any DASSL continuation whose
    OUTER TARGET time lies beyond ``t_wall`` (which sits PAST the
    conversion crossing) dies inside the integrator before the post-step
    termination check can run. Integrates normally otherwise, so tests
    stay fast and fail fast instead of hanging."""

    t_wall = 6.0

    def step(self, step_time):
        if step_time > self.t_wall:
            from pydas.dassl import DASSLError
            raise DASSLError(
                "DASSL returned with an IDID = -1, 500 steps taken on this "
                "call before reaching TOUT (the tf100 wall past the "
                "polymer-conversion crossing).")
        return HybridPolymerSystem.step(self, step_time)


def _polymer_conversion_rs(cls=HybridPolymerSystem, moments=(1.0, 5.0, 30.0),
                           monomer_mw_g_mol=42.08, chain_mass_defect_g_mol=0.0,
                           termination=None):
    """Single consuming pool A: same-pool VOLATILE_EJECTION row A -> A + G
    (archetype 6, eject a) drains mu1 at exactly d(mu1)/dt = -kf*mu1
    at exactly d(mu1)/dt = -kf*a*mu1 (empirically pinned: site = mu1/a,
    a units ejected per event, quadratic in a through the site*drain
    pairing), so with _KIN's kf = 2.0 1/s and a = 1.135 the defect-free
    polymer conversion X(t) = 1 - exp(-2.27 t) crosses 0.999 at
    t = ln(1000)/2.27 = 3.043 s."""
    sp = _gate17_species()
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"]]
    mask = np.array([False, False, False, False, True], dtype=bool)
    rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["X"]], **_KIN)
    rxn.polymer_flux_archetype = 6            # VOLATILE_EJECTION, same-pool
    rxn.polymer_eject_units = 1.135
    cfg = PolymerPoolConfig(label="A", xs=2,
                            explicit_dp_to_species_index={},
                            mu_indices=(1, 2, 3), monomer_poly_index=None,
                            monomer_mw_g_mol=monomer_mw_g_mol,
                            chain_mass_defect_g_mol=chain_mass_defect_g_mol)
    rs = cls(
        T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
        V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
        gas_species_mask=mask, constant_gas_volume=False,
        initial_polymer_moments={"A": tuple(moments)},
        termination=list(termination or []),
        allow_default_prospective_edge=True,
        allow_unstamped_proxy_rows=True)
    return rs, core, [rxn]


def _polymer_conversion_simulate(rs, core, rxns_core):
    from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
    ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e8,
                       tol_interrupt_simulation=1.0e8)
    return rs.simulate(list(core), list(rxns_core), [], [], [], [],
                       model_settings=ms,
                       simulator_settings=SimulatorSettings())


class TestTerminationPolymerConversion:
    """r86 terminationPolymerConversion: terminate the simulation when the
    defect-adjusted condensed polymer mass across ALL solver pools has
    dropped by the target fraction --
    X_polymer(t) = 1 - M_poly(t)/M_poly(0), with
    M_pool = max(0, mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol)
    summed over every pool (configured, spawned, daughters, side-group
    feature pools) plus explicit-DP oligomer carriers. Reachability is
    Codex-adjudicated (r86 NO-GO on post-step-only): once X_polymer > 0
    the next solver target time is CAPPED before the predicted crossing,
    so the crossing is bracketed before the post-chemistry stiff region
    (the FR1 tf100 wall). No bare terminationTime cut anywhere."""

    def test_termination_class_validates_bounds(self):
        from rmgpy.solver.termination import TerminationPolymerConversion
        term = TerminationPolymerConversion(0.999)
        assert term.conversion == 0.999
        for bad in (0.0, 1.0, -0.1, 1.5):
            with pytest.raises(ValueError,
                               match="terminationPolymerConversion"):
                TerminationPolymerConversion(bad)

    def test_metric_sums_all_pools_defect_adjusted_and_clamped(self):
        """M_poly = sum over ALL pools of max(0, mu1*MW - mu0*defect):
        defect pools use the schema-2.7 mass-defect formula (an FR1
        Br-loss pool must not falsely retain Br mass), and only each
        pool's FINAL contribution is clamped at zero (numerical fuzz) --
        never the individual moments."""
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0),
                          (2.0, 3.0, 10.0))
        # Give the fixture's pools the r86 mass parameters directly (the
        # PolymerPoolConfig fields the deck compiler fills).
        import dataclasses
        pool_a = dataclasses.replace(rs.polymer_pools[0],
                                     monomer_mw_g_mol=42.08,
                                     chain_mass_defect_g_mol=0.0)
        pool_b = dataclasses.replace(rs.polymer_pools[1],
                                     monomer_mw_g_mol=42.08,
                                     chain_mass_defect_g_mol=79.9)
        rs.polymer_pools = [pool_a, pool_b]
        m = rs.get_total_polymer_condensed_mass_g()
        expected_a = 5.0 * 42.08 - 1.0 * 0.0
        expected_b = 3.0 * 42.08 - 2.0 * 79.9
        assert expected_b < 0.0          # defect exceeds mass: fuzz shape
        assert m == pytest.approx(expected_a + max(0.0, expected_b))
        # The clamp is per-pool on the FINAL contribution: pool B clamps
        # to zero, it does NOT subtract from pool A.
        assert m == pytest.approx(expected_a)

    def test_metric_includes_explicit_dp_oligomer_mass(self):
        """r86 probe 3: a pool with explicit_dp species carries genuine
        polymer repeat-unit mass in an ordinary condensed core species
        OUTSIDE the mu slots (the init consistency check subtracts it from
        mu1); the conversion metric must add it back via
        _explicit_moment_contributions, not silently ignore it."""
        sp = _gate17_species()
        core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["X"],
                sp["G"]]
        mask = np.array([False, False, False, False, True, False],
                        dtype=bool)
        cfg = PolymerPoolConfig(label="A", xs=3,
                                explicit_dp_to_species_index={3: 5},
                                mu_indices=(1, 2, 3),
                                monomer_poly_index=None,
                                monomer_mw_g_mol=42.08,
                                chain_mass_defect_g_mol=0.0)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={core[4]: 1.0},
            V_poly=1.0, polymer_pools=[cfg], mass_transfer=[],
            gas_species_mask=mask, constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0)},
            initial_explicit_species={"A": {3: 0.5}},
            termination=[])
        rs.initialize_model(core, [], [], [])
        # Init subtracted the oligomer share from the mu1 STATE slot; the
        # metric must report the TOTAL: (mu1_tail + 3*0.5)*MW = 5.0*MW.
        m = rs.get_total_polymer_condensed_mass_g()
        assert float(np.asarray(rs.y)[2]) == pytest.approx(5.0 - 1.5)
        assert m == pytest.approx(5.0 * 42.08)

    def test_simulate_rejects_nonpositive_initial_polymer_mass(self):
        """M_poly(0) is frozen at simulation initialization and must be
        strictly positive -- a zero-mass (or defect-swamped) initial pool
        state is a loud ValueError, never a silent 0/0."""
        from rmgpy.solver.termination import TerminationPolymerConversion
        rs, core, rxns = _polymer_conversion_rs(
            moments=(0.0, 0.0, 0.0),
            termination=[TerminationPolymerConversion(0.999)])
        with pytest.raises(ValueError, match="M_poly\\(0\\)"):
            _polymer_conversion_simulate(rs, core, rxns)

    def test_simulate_terminates_cleanly_past_stiff_wall(self, caplog):
        """Live-path reachability pin (the r86 NO-GO shape, bounded): the
        post-conversion region is a stiff wall at t = 6.0 s, PAST the
        X = 0.999 crossing at t = ln(1000)/2.27 = 3.043 s. Without
        conservative target capping the geometric step schedule requests
        a target of 10 s from t ~ 1 s and dies inside DASSL before any
        post-step check runs (RED shape); with capping the crossing is
        bracketed, termination fires, completion status is IDENTICAL to a
        normal termination (terminated=True, no resurrection), and the
        X_polymer value is reported in the log."""
        from rmgpy.solver.termination import TerminationPolymerConversion
        rs, core, rxns = _polymer_conversion_rs(
            cls=_StiffPastCrossingSystem,
            termination=[TerminationPolymerConversion(0.999)])
        with caplog.at_level(logging.INFO):
            terminated, resurrected, invalid_objects, _, _, t_final, _ = \
                _polymer_conversion_simulate(rs, core, rxns)
        assert terminated is True
        assert resurrected is False
        assert invalid_objects == []
        # Crossing bracketed: we stopped just past ln(1000)/(kf*a), far
        # from the 6.0 s wall and further still from the 10 s uncapped
        # target.
        assert t_final == pytest.approx(np.log(1000.0) / 2.27, rel=0.05)
        assert t_final < _StiffPastCrossingSystem.t_wall
        msgs = [r.getMessage() for r in caplog.records]
        assert any("termination polymer conversion" in m and "X_polymer" in m
                   for m in msgs)

    def test_simulate_without_wall_reports_conversion_value(self, caplog):
        """Plain live-path completion (no wall subclass): the criterion
        fires on its own -- no terminationTime backstop anywhere in the
        termination list -- and reports X_polymer >= the target."""
        import re
        from rmgpy.solver.termination import TerminationPolymerConversion
        rs, core, rxns = _polymer_conversion_rs(
            termination=[TerminationPolymerConversion(0.9)])
        with caplog.at_level(logging.INFO):
            terminated, _, _, _, _, t_final, _ = \
                _polymer_conversion_simulate(rs, core, rxns)
        assert terminated is True
        assert t_final == pytest.approx(np.log(10.0) / 2.27, rel=0.05)
        msgs = [r.getMessage() for r in caplog.records
                if "termination polymer conversion" in r.getMessage()]
        assert msgs
        x_val = float(re.search(r"X_polymer = ([0-9.eE+-]+)",
                                msgs[0]).group(1))
        assert x_val >= 0.9

    def test_to_solver_object_plumb_and_no_pool_rejection(self):
        """Deck surface: terminationPolymerConversion is accepted by the
        HybridPolymerReactor input object, forwarded by the input-file
        reader, materialized as a TerminationPolymerConversion solver
        criterion -- and rejected loudly when the model has NO polymer
        pools (nothing for the metric to measure)."""
        import inspect
        from rmgpy.quantity import Quantity
        from rmgpy.polymer import Polymer
        from rmgpy.rmg.polymer_input import (HybridPolymerReactor,
                                             PolymerPhase)
        from rmgpy.solver.termination import TerminationPolymerConversion

        a = _spc("CCCC", "A")
        daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                           end_groups=["[CH3]", "[H]"], cutoff=3,
                           Mn=5000.0, Mw=6000.0, initial_mass=0.001)
        mu0 = _spc("[Ne]", "PS_d1_mu0"); mu0.reactive = False
        mu0.is_moment_dummy = True
        mu1 = _spc("[Ne]", "PS_d1_mu1"); mu1.reactive = False
        mu1.is_moment_dummy = True
        mu2 = _spc("[Ne]", "PS_d1_mu2"); mu2.reactive = False
        mu2.is_moment_dummy = True
        core = [a, daughter, mu0, mu1, mu2]
        phase = PolymerPhase(density=Quantity(1050.0, "kg/m^3"),
                             initial_moments={},
                             initial_explicit={a: 1.0}, pools=[],
                             mass_transfer=[])
        reactor = HybridPolymerReactor(
            temperature=(1000.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles={a: 1.0}, polymerPhase=phase,
            terminationPolymerConversion=0.999)
        solver = reactor.to_solver_object(core, [], [], [])
        tpc = [t for t in solver.termination
               if isinstance(t, TerminationPolymerConversion)]
        assert len(tpc) == 1
        assert tpc[0].conversion == 0.999

        # NO polymer pools anywhere (no deck pools, no core daughters):
        # loud rejection, never a silently ignored criterion.
        reactor_no_pool = HybridPolymerReactor(
            temperature=(1000.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles={a: 1.0}, polymerPhase=PolymerPhase(
                density=Quantity(1050.0, "kg/m^3"), initial_moments={},
                initial_explicit={a: 1.0}, pools=[], mass_transfer=[]),
            terminationPolymerConversion=0.999)
        with pytest.raises(ValueError, match="no polymer pools"):
            reactor_no_pool.to_solver_object([a], [], [], [])

        # Out-of-range values die at the single validation chokepoint.
        reactor_bad = HybridPolymerReactor(
            temperature=(1000.0, "K"), pressure=(1.0e5, "Pa"),
            initialMoles={a: 1.0}, polymerPhase=phase,
            terminationPolymerConversion=1.5)
        with pytest.raises(ValueError,
                           match="terminationPolymerConversion"):
            reactor_bad.to_solver_object(core, [], [], [])

        # Input-file reader carries the kwarg (static pin, same pattern as
        # the allow_unpaired_reference_state knob).
        from rmgpy.rmg.input import hybrid_polymer_reactor
        sig = inspect.signature(hybrid_polymer_reactor)
        assert "terminationPolymerConversion" in sig.parameters
        assert sig.parameters["terminationPolymerConversion"].default is None
        src = inspect.getsource(hybrid_polymer_reactor)
        assert ("terminationPolymerConversion=terminationPolymerConversion"
                in src)


class TestSpawnedPoolDemotionRefusal:
    """r91 (PP run-10 death): a stamped pool-coupled row whose required pool
    endpoint resolves to -1 BECAUSE the endpoint is a mid-run-SPAWNED,
    config-less pool (a Polymer species present in the species lists but
    absent from the solver's pool configs) must be REFUSED conduit-deferred
    (reaction_refused -> zero flux, same mechanics as every other refused
    row) instead of DEMOTED to the legacy mu1-only transfer. The legacy
    demotion applied an unapportioned mu1-only drain that drove pool
    polypropylene_mod to mu0=+0.0047 / mu1=-5.6e-5 / mu2=+0.235 at
    simulate-leg 21 (the r81 negative-beyond-floor tripwire killed run-10).
    Unresolved endpoints NOT attributable to a spawned config-less pool
    (e.g. a plain-Species daughter) keep the legacy demotion, censused
    separately as a loud anomaly."""

    @staticmethod
    def _daughter(label="poly_mod"):
        from rmgpy.polymer import Polymer
        return Polymer(label=label, monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)

    @staticmethod
    def _species():
        return {
            "Proxy": _spc("CCCC", "poly"),
            "Mu0": _spc("CO", "poly_mu0"),
            "Mu1": _spc("C=O", "poly_mu1"),
            "Mu2": _spc("C#N", "poly_mu2"),
            "Vol": _spc("C=C", "vol_gas"),
            "B": _spc("[CH3]", "B_gas"),
        }

    @staticmethod
    def _build(core, mask, rxns, moments=(1.0, 5.0, 30.0)):
        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=np.array(mask, dtype=bool),
            constant_gas_volume=False,
            initial_polymer_moments={"poly": moments}, termination=[],
        )
        rs.initialize_model(core, list(rxns), [], [])
        return rs

    def test_run10_shape_migration_to_spawned_pool_refused_not_demoted(self, caplog):
        """RED-matrix 1: the run-10 shape. A stamped MIGRATION row whose dst
        is a spawned config-less Polymer core species must be refused
        (zero flux), not demoted to legacy mu1-only."""
        import rmgpy.solver.polymer as sp_mod
        s = self._species()
        d = self._daughter()
        core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], d]
        mask = [False] * 5
        rxn = Reaction(reactants=[s["Proxy"]], products=[d], **_KIN)
        rxn.polymer_flux_archetype = sp_mod.FLUX_MIGRATION
        with caplog.at_level(logging.WARNING):
            rs = self._build(core, mask, [rxn])
        assert rs.reaction_refused[0] == 1
        assert rs.reaction_flux_archetype[0] == sp_mod.FLUX_NONE
        # zero flux: elementwise identical to the row-free system
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        rs0 = self._build(core, mask, [])
        dn0 = rs0.residual(0.0, rs0.y, np.zeros_like(rs0.y))[0]
        assert np.array_equal(dn, dn0)
        assert any("SPAWNED-POOL DEMOTION REFUSAL" in r.getMessage()
                   for r in caplog.records)
        assert not any("could not resolve their solver pool(s)"
                       in r.getMessage() for r in caplog.records)

    @pytest.mark.parametrize("shape", ["scission_fragment", "discrete_chip",
                                       "volatile_ejection"])
    def test_all_archetypes_targeting_spawned_pool_are_refused(self, shape):
        """RED-matrix 1 (archetype coverage): SCISSION_FRAGMENT with a
        spawned dst, DISCRETE_CHIP with a spawned src, and VOLATILE_EJECTION
        with a spawned (unresolved) required endpoint all refuse."""
        import rmgpy.solver.polymer as sp_mod
        s = self._species()
        d = self._daughter()
        if shape == "scission_fragment":
            core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], d]
            mask = [False] * 5
            rxn = Reaction(reactants=[s["Proxy"]], products=[d], **_KIN)
            rxn.polymer_flux_archetype = sp_mod.FLUX_SCISSION_FRAGMENT
        elif shape == "discrete_chip":
            core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], d, s["B"]]
            mask = [False] * 5 + [True]
            rxn = Reaction(reactants=[d], products=[s["B"]], **_KIN)
            rxn.polymer_flux_archetype = sp_mod.FLUX_DISCRETE_CHIP
            rxn.polymer_chip_units = 2
        else:
            core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], d, s["Vol"]]
            mask = [False] * 5 + [True]
            rxn = Reaction(reactants=[s["Proxy"]], products=[d, s["Vol"]],
                           **_KIN)
            rxn.polymer_flux_archetype = sp_mod.FLUX_VOLATILE_EJECTION
            rxn.polymer_eject_units = 0.5
        rs = self._build(core, mask, [rxn])
        assert rs.reaction_refused[0] == 1
        assert rs.reaction_flux_archetype[0] == sp_mod.FLUX_NONE
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        rs0 = self._build(core, mask, [])
        dn0 = rs0.residual(0.0, rs0.y, np.zeros_like(rs0.y))[0]
        assert np.array_equal(dn, dn0)

    def test_refused_row_contributes_zero_to_all_three_moment_legs(self):
        """RED-matrix 2 (residual-level): the run-10 shape no longer drives
        mu1 negative -- the refused row contributes exactly zero to mu0, mu1
        and mu2 (baseline demotion applied dn[mu1] = -k*mu1 < 0)."""
        import rmgpy.solver.polymer as sp_mod
        s = self._species()
        d = self._daughter()
        core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], d]
        mask = [False] * 5
        rxn = Reaction(reactants=[s["Proxy"]], products=[d], **_KIN)
        rxn.polymer_flux_archetype = sp_mod.FLUX_MIGRATION
        rs = self._build(core, mask, [rxn])
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert dn[1] == 0.0    # mu0 leg
        assert dn[2] == 0.0    # mu1 leg: baseline drove this negative
        assert dn[3] == 0.0    # mu2 leg
        assert dn[4] == 0.0    # spawned daughter gains nothing

    def test_census_line_format_and_payload(self, caplog):
        """r91 census format + solver-side census attribute."""
        import re
        import rmgpy.solver.polymer as sp_mod
        s = self._species()
        d = self._daughter()
        core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], d]
        rxn = Reaction(reactants=[s["Proxy"]], products=[d], **_KIN)
        rxn.polymer_flux_archetype = sp_mod.FLUX_MIGRATION
        with caplog.at_level(logging.WARNING):
            rs = self._build(core, [False] * 5, [rxn])
        msgs = [r.getMessage() for r in caplog.records
                if "SPAWNED-POOL DEMOTION REFUSAL" in r.getMessage()]
        assert len(msgs) == 1, msgs
        assert re.search(
            r"SPAWNED-POOL DEMOTION REFUSAL: 1 stamped pool-coupled row\(s\) "
            r"targeting unconfigured spawned pools refused conduit-deferred "
            r"instead of demoted to legacy mu1-only; item-16 pending; "
            r"archetypes=MIGRATION, pools=poly_mod, first_rows=core row 0:",
            msgs[0]), msgs[0]
        census = rs.spawned_pool_refusal_census
        assert len(census) == 1
        assert census[0]["reason"] == "conduit-deferred"
        assert census[0]["archetype"] == "MIGRATION"
        assert census[0]["pools"] == ["poly_mod"]

    def test_refusal_not_sticky_object_stays_unstamped(self):
        """The refusal is re-derived per rebuild and must NOT stamp
        polymer_refused on the reaction object: once item 16 configures the
        spawned pool, the SAME row must resolve and run fully apportioned."""
        import rmgpy.solver.polymer as sp_mod
        s = self._species()
        d = self._daughter()
        core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], d]
        rxn = Reaction(reactants=[s["Proxy"]], products=[d], **_KIN)
        rxn.polymer_flux_archetype = sp_mod.FLUX_MIGRATION
        rs = self._build(core, [False] * 5, [rxn])
        assert rs.reaction_refused[0] == 1
        assert not getattr(rxn, "polymer_refused", False)
        # negative control (RED-matrix 4): the same row against a
        # configuration that RESOLVES both pools stays fully live.
        dmu0 = _spc("CCO", "poly_mod_mu0")
        dmu1 = _spc("CC=O", "poly_mod_mu1")
        dmu2 = _spc("CC#N", "poly_mod_mu2")
        core2 = core + [dmu0, dmu1, dmu2]
        parent = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
        spawned = PolymerPoolConfig(
            label="poly_mod", xs=2, explicit_dp_to_species_index={},
            mu_indices=(5, 6, 7), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
        rs2 = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={}, V_poly=1.0,
            polymer_pools=[parent, spawned], mass_transfer=[],
            gas_species_mask=np.array([False] * 8, dtype=bool),
            constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0),
                                     "poly_mod": (0.0, 0.0, 0.0)},
            termination=[])
        rs2.initialize_model(core2, [rxn], [], [])
        assert rs2.reaction_refused[0] == 0
        assert rs2.reaction_flux_archetype[0] == sp_mod.FLUX_MIGRATION

    def test_plain_species_daughter_keeps_demotion_and_censuses_anomaly(self, caplog):
        """Negative control + anomaly census: an unresolved endpoint NOT
        attributable to a spawned config-less Polymer (plain-Species
        daughter, the test_stamped_scission_without_daughter_pool shape)
        keeps the legacy mu1-only demotion, is NOT refused, and is censused
        as a loud anomaly."""
        import rmgpy.solver.polymer as sp_mod
        s = self._species()
        plain = _spc("CCC", "poly_scission_tail")   # NOT a Polymer
        core = [s["Proxy"], s["Mu0"], s["Mu1"], s["Mu2"], plain]
        rxn = Reaction(reactants=[s["Proxy"]], products=[plain], **_KIN)
        rxn.polymer_flux_archetype = sp_mod.FLUX_SCISSION_FRAGMENT
        with caplog.at_level(logging.WARNING):
            rs = self._build(core, [False] * 5, [rxn])
        assert rs.reaction_refused[0] == 0
        assert rs.reaction_flux_archetype[0] == sp_mod.FLUX_UNRESOLVED
        assert rs.spawned_pool_refusal_census == []
        # legacy mu1 drain still applies (behavior unchanged)
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        assert np.isclose(dn[2], -kf * 5.0)
        assert any("demoted to legacy" in r.getMessage()
                   for r in caplog.records)
        assert any("DEMOTION ANOMALY" in r.getMessage()
                   for r in caplog.records)
