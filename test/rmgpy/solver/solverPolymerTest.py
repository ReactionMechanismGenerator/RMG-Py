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
        import rmgpy.solver.polymer as sp
        assert sp.FLUX_NONE == int(PolymerFluxArchetype.NONE) == 0
        assert sp.FLUX_SAME_POOL == int(PolymerFluxArchetype.SAME_POOL) == 1
        assert sp.FLUX_MIGRATION == int(PolymerFluxArchetype.MIGRATION) == 2
        assert sp.FLUX_SCISSION_FRAGMENT == int(PolymerFluxArchetype.SCISSION_FRAGMENT) == 3
        assert sp.FLUX_UNRESOLVED == int(PolymerFluxArchetype.UNRESOLVED) == 4

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

        import rmgpy.solver.polymer as sp
        assert rs.reaction_flux_archetype[0] == sp.FLUX_UNRESOLVED
        assert rs.reaction_src_pool[0] == 0
        assert rs.reaction_dst_pool[0] == -1        # gas-only products
        assert rs.reaction_flux_archetype[1] == sp.FLUX_NONE  # pure-gas stays NONE
        assert rs.reaction_src_pool[1] == -1
        assert rs.reaction_dst_pool[1] == -1

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
