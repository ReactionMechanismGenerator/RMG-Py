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

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

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
        continuity), and places the daughter's event-seeded moments (mu_k = N*DP^k)
        at its own state-vector slots.

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
        N, DP = 2.0e-4, 5
        intent = SpawnIntent(parent_pool=parent, monomer=parent.backbone_group,
                             end_groups=["[H]", "[H]"], triggering_dp=DP, triggering_moles=N)
        daughter = drain_spawn_intents([intent], iteration=1, existing_pools=[parent])[0]
        assert daughter.label == "PE_d1"
        assert np.allclose(daughter.moments, [N, N * DP, N * DP * DP])

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
        # Daughter's event-seeded moments (N*DP^k) land at its own slots.
        assert rs2.y[4] == pytest.approx(N)
        assert rs2.y[5] == pytest.approx(N * DP)
        assert rs2.y[6] == pytest.approx(N * DP * DP)
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
