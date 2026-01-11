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
