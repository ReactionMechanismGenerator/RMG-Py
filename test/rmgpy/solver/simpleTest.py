#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

import os


import numpy as np

import rmgpy.constants as constants
from rmgpy.chemkin import load_chemkin_file
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import TerminationTime
from rmgpy.solver.simple import SimpleReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


class SimpleReactorCheck:
    def test_solve(self):
        """
        Test the simple batch reactor with a simple kinetic model. Here we
        choose a kinetic model consisting of the hydrogen abstraction reaction
        CH4 + C2H5 <=> CH3 + C2H6.
        """
        ch4 = Species(
            molecule=[Molecule().from_smiles("C")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=(
                    [8.615, 9.687, 10.963, 12.301, 14.841, 16.976, 20.528],
                    "cal/(mol*K)",
                ),
                H298=(-17.714, "kcal/mol"),
                S298=(44.472, "cal/(mol*K)"),
            ),
        )
        ch3 = Species(
            molecule=[Molecule().from_smiles("[CH3]")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=(
                    [9.397, 10.123, 10.856, 11.571, 12.899, 14.055, 16.195],
                    "cal/(mol*K)",
                ),
                H298=(9.357, "kcal/mol"),
                S298=(45.174, "cal/(mol*K)"),
            ),
        )
        c2h6 = Species(
            molecule=[Molecule().from_smiles("CC")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=(
                    [12.684, 15.506, 18.326, 20.971, 25.500, 29.016, 34.595],
                    "cal/(mol*K)",
                ),
                H298=(-19.521, "kcal/mol"),
                S298=(54.799, "cal/(mol*K)"),
            ),
        )
        c2h5 = Species(
            molecule=[Molecule().from_smiles("C[CH2]")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=(
                    [11.635, 13.744, 16.085, 18.246, 21.885, 24.676, 29.107],
                    "cal/(mol*K)",
                ),
                H298=(29.496, "kcal/mol"),
                S298=(56.687, "cal/(mol*K)"),
            ),
        )

        rxn1 = Reaction(
            reactants=[c2h6, ch3],
            products=[c2h5, ch4],
            kinetics=Arrhenius(
                A=(686.375 * 6, "m^3/(mol*s)"),
                n=4.40721,
                Ea=(7.82799, "kcal/mol"),
                T0=(298.15, "K"),
            ),
        )

        core_species = [ch4, ch3, c2h6, c2h5]
        edge_species = []
        core_reactions = [rxn1]
        edge_reactions = []

        T = 1000
        P = 1.0e5
        rxn_system = SimpleReactor(
            T,
            P,
            initial_mole_fractions={c2h5: 0.1, ch3: 0.1, ch4: 0.4, c2h6: 0.4},
            n_sims=1,
            termination=[],
        )

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.array([10 ** (i / 10.0) for i in range(-130, -49)], float)

        # Integrate to get the solution at each time point
        t = []
        y = []
        reaction_rates = []
        species_rates = []
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxn_system.y.copy())
            reaction_rates.append(rxn_system.core_reaction_rates.copy())
            species_rates.append(rxn_system.core_species_rates.copy())

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        V = constants.R * rxn_system.T.value_si * np.sum(y) / rxn_system.P.value_si

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - species_rates[i, 0]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - -species_rates[i, 1]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - -species_rates[i, 2]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - species_rates[i, 3]) < 1e-6 * reaction_rates[i, 0]

        # Check that we've reached equilibrium
        assert abs(reaction_rates[-1, 0] - 0.0) < 1e-2

        # Unit test for the jacobian function:
        # Solve a reaction system and check if the analytical jacobian matches the finite difference jacobian

        h2 = Species(
            molecule=[Molecule().from_smiles("[H][H]")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([6.89, 6.97, 6.99, 7.01, 7.08, 7.22, 7.72], "cal/(mol*K)"),
                H298=(0, "kcal/mol"),
                S298=(31.23, "cal/(mol*K)"),
            ),
        )

        rxn_list = [
            Reaction(
                reactants=[c2h6],
                products=[ch3, ch3],
                kinetics=Arrhenius(
                    A=(686.375 * 6, "1/s"),
                    n=4.40721,
                    Ea=(7.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[ch3, ch3],
                products=[c2h6],
                kinetics=Arrhenius(
                    A=(686.375 * 6, "m^3/(mol*s)"),
                    n=4.40721,
                    Ea=(7.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h6, ch3],
                products=[c2h5, ch4],
                kinetics=Arrhenius(
                    A=(46.375 * 6, "m^3/(mol*s)"),
                    n=3.40721,
                    Ea=(6.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h5, ch4],
                products=[c2h6, ch3],
                kinetics=Arrhenius(
                    A=(46.375 * 6, "m^3/(mol*s)"),
                    n=3.40721,
                    Ea=(6.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h5, ch4],
                products=[ch3, ch3, ch3],
                kinetics=Arrhenius(
                    A=(246.375 * 6, "m^3/(mol*s)"),
                    n=1.40721,
                    Ea=(3.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[ch3, ch3, ch3],
                products=[c2h5, ch4],
                kinetics=Arrhenius(
                    A=(246.375 * 6, "m^6/(mol^2*s)"),
                    n=1.40721,
                    Ea=(3.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h6, ch3, ch3],
                products=[c2h5, c2h5, h2],
                kinetics=Arrhenius(
                    A=(146.375 * 6, "m^6/(mol^2*s)"),
                    n=2.40721,
                    Ea=(8.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h5, c2h5, h2],
                products=[c2h6, ch3, ch3],
                kinetics=Arrhenius(
                    A=(146.375 * 6, "m^6/(mol^2*s)"),
                    n=2.40721,
                    Ea=(8.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h6, c2h6],
                products=[ch3, ch4, c2h5],
                kinetics=Arrhenius(
                    A=(1246.375 * 6, "m^3/(mol*s)"),
                    n=0.40721,
                    Ea=(8.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[ch3, ch4, c2h5],
                products=[c2h6, c2h6],
                kinetics=Arrhenius(
                    A=(46.375 * 6, "m^6/(mol^2*s)"),
                    n=0.10721,
                    Ea=(8.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
        ]

        for rxn in rxn_list:
            core_species = [ch4, ch3, c2h6, c2h5, h2]
            edge_species = []
            core_reactions = [rxn]

            rxn_system0 = SimpleReactor(
                T,
                P,
                initial_mole_fractions={
                    ch4: 0.2,
                    ch3: 0.1,
                    c2h6: 0.35,
                    c2h5: 0.15,
                    h2: 0.2,
                },
                n_sims=1,
                termination=[],
            )
            rxn_system0.initialize_model(core_species, core_reactions, edge_species, edge_reactions)
            dydt0 = rxn_system0.residual(0.0, rxn_system0.y, np.zeros(rxn_system0.y.shape))[0]
            num_core_species = len(core_species)
            dN = 0.000001 * sum(rxn_system0.y)
            dN_array = dN * np.eye(num_core_species)

            dydt = []
            for i in range(num_core_species):
                rxn_system0.y[i] += dN
                dydt.append(rxn_system0.residual(0.0, rxn_system0.y, np.zeros(rxn_system0.y.shape))[0])
                rxn_system0.y[i] -= dN  # reset y to original y0

            # Let the solver compute the jacobian
            solver_jacobian = rxn_system0.jacobian(0.0, rxn_system0.y, dydt0, 0.0)
            # Compute the jacobian using finite differences
            jacobian = np.zeros((num_core_species, num_core_species))
            for i in range(num_core_species):
                for j in range(num_core_species):
                    jacobian[i, j] = (dydt[j][i] - dydt0[i]) / dN
                    assert abs(jacobian[i, j] - solver_jacobian[i, j]) < abs(1e-4 * jacobian[i, j])

        # print 'Solver jacobian'
        # print solver_jacobian
        # print 'Numerical jacobian'
        # print jacobian

        # Unit test for the compute rate derivative
        rxn_list = [
            Reaction(
                reactants=[c2h6],
                products=[ch3, ch3],
                kinetics=Arrhenius(
                    A=(686.375e6, "1/s"),
                    n=4.40721,
                    Ea=(7.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h6, ch3],
                products=[c2h5, ch4],
                kinetics=Arrhenius(
                    A=(46.375 * 6, "m^3/(mol*s)"),
                    n=3.40721,
                    Ea=(6.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
            Reaction(
                reactants=[c2h6, ch3, ch3],
                products=[c2h5, c2h5, h2],
                kinetics=Arrhenius(
                    A=(146.375 * 6, "m^6/(mol^2*s)"),
                    n=2.40721,
                    Ea=(8.82799, "kcal/mol"),
                    T0=(298.15, "K"),
                ),
            ),
        ]

        core_species = [ch4, ch3, c2h6, c2h5, h2]
        edge_species = []
        core_reactions = rxn_list

        rxn_system0 = SimpleReactor(
            T,
            P,
            initial_mole_fractions={
                ch4: 0.2,
                ch3: 0.1,
                c2h6: 0.35,
                c2h5: 0.15,
                h2: 0.2,
            },
            n_sims=1,
            termination=[],
        )
        rxn_system0.initialize_model(core_species, core_reactions, edge_species, edge_reactions)
        dfdt0 = rxn_system0.residual(0.0, rxn_system0.y, np.zeros(rxn_system0.y.shape))[0]
        solver_dfdk = rxn_system0.compute_rate_derivative()
        # print 'Solver d(dy/dt)/dk'
        # print solver_dfdk

        integration_time = 1e-8
        rxn_system0.termination.append(TerminationTime((integration_time, "s")))
        model_settings = ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=0)
        simulator_settings = SimulatorSettings()
        rxn_system0.simulate(
            core_species,
            core_reactions,
            [],
            [],
            [],
            [],
            model_settings=model_settings,
            simulator_settings=simulator_settings,
        )

        y0 = rxn_system0.y

        dfdk = np.zeros((num_core_species, len(rxn_list)))  # d(dy/dt)/dk

        for i in range(len(rxn_list)):
            k0 = rxn_list[i].get_rate_coefficient(T, P)
            rxn_list[i].kinetics.A.value_si = rxn_list[i].kinetics.A.value_si * (1 + 1e-3)
            dk = rxn_list[i].get_rate_coefficient(T, P) - k0

            rxn_system = SimpleReactor(
                T,
                P,
                initial_mole_fractions={
                    ch4: 0.2,
                    ch3: 0.1,
                    c2h6: 0.35,
                    c2h5: 0.15,
                    h2: 0.2,
                },
                n_sims=1,
                termination=[],
            )
            rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

            dfdt = rxn_system.residual(0.0, rxn_system.y, np.zeros(rxn_system.y.shape))[0]
            dfdk[:, i] = (dfdt - dfdt0) / dk

            rxn_system.termination.append(TerminationTime((integration_time, "s")))
            model_settings = ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=0)
            simulator_settings = SimulatorSettings()

            rxn_system.simulate(
                core_species,
                core_reactions,
                [],
                [],
                [],
                [],
                model_settings=model_settings,
                simulator_settings=simulator_settings,
            )

            rxn_list[i].kinetics.A.value_si = rxn_list[i].kinetics.A.value_si / (1 + 1e-3)  # reset A factor

        for i in range(num_core_species):
            for j in range(len(rxn_list)):
                assert abs(dfdk[i, j] - solver_dfdk[i, j]) < abs(1e-3 * dfdk[i, j])

        # print 'Numerical d(dy/dt)/dk'
        # print dfdk

        # # Visualize the simulation results
        # import pylab
        # fig = pylab.figure(figsize=(6,6))
        # pylab.subplot(2,1,1)
        # pylab.semilogx(t, y)
        # pylab.ylabel('Concentration (mol/m$^\\mathdefault{3}$)')
        # pylab.legend(['CH4', 'CH3', 'C2H6', 'C2H5'], loc=4)
        # pylab.subplot(2,1,2)
        # pylab.semilogx(t, species_rates)
        # pylab.legend(['CH4', 'CH3', 'C2H6', 'C2H5'], loc=4)
        # pylab.xlabel('Time (s)')
        # pylab.ylabel('Rate (mol/m$^\\mathdefault{3}$*s)')
        # fig.subplots_adjust(left=0.12, bottom=0.10, right=0.95, top=0.95, wspace=0.20, hspace=0.35)
        # pylab.show()

    def test_collider_model(self):
        """
        Test the solver's ability to simulate a model with collision efficiencies.
        """
        chem_file = os.path.join(os.path.dirname(__file__), "files", "collider_model", "chem.inp")
        dictionary_file = os.path.join(
            os.path.dirname(__file__),
            "files",
            "collider_model",
            "species_dictionary.txt",
        )
        species_list, reaction_list = load_chemkin_file(chem_file, dictionary_file)

        smiles_dict = {
            "H": "[H]",
            "HO2": "[O]O",
            "O2": "[O][O]",
            "Ar": "[Ar]",
            "N2": "N#N",
            "CO2": "O=C=O",
            "CH3": "[CH3]",
            "CH4": "C",
        }
        species_dict = {}
        for name, smiles in smiles_dict.items():
            mol = Molecule(smiles=smiles)
            for species in species_list:
                if species.is_isomorphic(mol):
                    species_dict[name] = species
                    break

        T = 1000  # K
        P = 10  # Pa
        initial_mole_fractions = {
            species_dict["O2"]: 0.5,
            species_dict["H"]: 0.5,
            species_dict["CO2"]: 1.0,
            species_dict["Ar"]: 4.0,
        }

        # Initialize the model
        rxn_system = SimpleReactor(
            T,
            P,
            initial_mole_fractions=initial_mole_fractions,
            n_sims=1,
            termination=None,
        )
        rxn_system.initialize_model(species_list, reaction_list, [], [])

        # Advance to time = 0.1 s
        rxn_system.advance(0.1)
        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        simulated_mole_fracs = rxn_system.y / np.sum(rxn_system.y)
        expected_mole_fracs = np.array(
            [
                0.6666667,
                0,
                0,
                0,
                0.1666667,
                0,
                0.08333333,
                0.08333333,
                2.466066000000000e-10,
                0,
                0,
                0,
                0,
                0,
            ]
        )
        for i in range(len(simulated_mole_fracs)):
            assert round(abs(simulated_mole_fracs[i] - expected_mole_fracs[i]), 7) == 0

        # Advance to time = 5 s
        rxn_system.advance(5)
        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        expected_mole_fracs = np.array(
            [
                0.6666667,
                0,
                0,
                0,
                0.1666667,
                0,
                0.08333332,
                0.08333332,
                1.233033000000000e-08,
                0,
                0,
                0,
                0,
                0,
            ]
        )
        for i in range(len(simulated_mole_fracs)):
            assert round(abs(simulated_mole_fracs[i] - expected_mole_fracs[i]), 7) == 0

        # Try a new set of conditions

        T = 850  # K
        P = 200  # Pa
        initial_mole_fractions = {
            species_dict["O2"]: 0.5,
            species_dict["H"]: 1,
            species_dict["CO2"]: 1,
            species_dict["N2"]: 4,
            species_dict["CH3"]: 1,
        }

        # Initialize the model
        rxn_system = SimpleReactor(
            T,
            P,
            initial_mole_fractions=initial_mole_fractions,
            n_sims=1,
            termination=None,
        )
        rxn_system.initialize_model(species_list, reaction_list, [], [])

        # Advance to time = 5 s
        rxn_system.advance(5)

        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        simulated_mole_fracs = rxn_system.y / np.sum(rxn_system.y)
        expected_mole_fracs = np.array(
            [
                0,
                0,
                0,
                0.5487241,
                0.137181,
                0,
                0.1083234,
                0.0685777,
                1.280687000000000e-05,
                0,
                0,
                0,
                0.1083362,
                0.02884481,
            ]
        )
        for i in range(len(simulated_mole_fracs)):
            assert round(abs(simulated_mole_fracs[i] - expected_mole_fracs[i]), 7) == 0

    def test_specific_collider_model(self):
        """
        Test the solver's ability to simulate a model with specific third body species collision efficiencies.
        """
        chem_file = os.path.join(os.path.dirname(__file__), "files", "specific_collider_model", "chem.inp")
        dictionary_file = os.path.join(
            os.path.dirname(__file__),
            "files",
            "specific_collider_model",
            "species_dictionary.txt",
        )
        species_list, reaction_list = load_chemkin_file(chem_file, dictionary_file)

        # Note that H2O is in the corresponding Chemkin file, but intentionally not defined here to check the solver.
        smiles_dict = {
            "Ar": "[Ar]",
            "N2(1)": "N#N",
            "O2": "[O][O]",
            "H": "[H]",
            "CH3": "[CH3]",
            "CH4": "C",
        }
        species_dict = {}
        for name, smiles in smiles_dict.items():
            mol = Molecule(smiles=smiles)
            for species in species_list:
                if species.is_isomorphic(mol):
                    species_dict[name] = species
                    break

        T = 1000  # K
        P = 10  # Pa
        initial_mole_fractions = {
            species_dict["Ar"]: 2.0,
            species_dict["N2(1)"]: 1.0,
            species_dict["O2"]: 0.5,
            species_dict["H"]: 0.1,
            species_dict["CH3"]: 0.1,
            species_dict["CH4"]: 0.001,
        }

        # Initialize the model
        rxn_system = SimpleReactor(
            T,
            P,
            initial_mole_fractions=initial_mole_fractions,
            n_sims=1,
            termination=None,
        )
        rxn_system.initialize_model(species_list, reaction_list, [], [])

        # Advance to time = 0.1 s
        rxn_system.advance(0.1)
        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        simulated_mole_fracs = rxn_system.y / np.sum(rxn_system.y)
        expected_mole_fracs = np.array(
            [
                0.540394532,
                0.270197216,
                0.135098608,
                0.027019722,
                0.027019722,
                0.000270202,
                0.0,
            ]
        )
        # order: Ar, N2, O2, H, CH3, CH4
        for i in range(len(simulated_mole_fracs)):
            assert round(abs(simulated_mole_fracs[i] - expected_mole_fracs[i]), 6) == 0

        # Advance to time = 5 s
        rxn_system.advance(5)
        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        expected_mole_fracs = np.array(
            [
                0.540394573,
                0.270197287,
                0.135098693,
                0.027019519,
                0.027019519,
                0.00027041,
                0.0,
            ]
        )
        # order: Ar, N2, O2, H, CH3, CH4
        for i in range(len(simulated_mole_fracs)):
            assert round(abs(simulated_mole_fracs[i] - expected_mole_fracs[i]), 6) == 0

        # Try a new set of conditions
        T = 850  # K
        P = 200  # Pa
        initial_mole_fractions = {
            species_dict["Ar"]: 1.0,
            species_dict["N2(1)"]: 0.5,
            species_dict["O2"]: 0.5,
            species_dict["H"]: 0.001,
            species_dict["CH3"]: 0.01,
            species_dict["CH4"]: 0.5,
        }

        # Initialize the model
        rxn_system = SimpleReactor(
            T,
            P,
            initial_mole_fractions=initial_mole_fractions,
            n_sims=1,
            termination=None,
        )
        rxn_system.initialize_model(species_list, reaction_list, [], [])

        # Advance to time = 5 s
        rxn_system.advance(5)
        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        simulated_mole_fracs = rxn_system.y / np.sum(rxn_system.y)
        expected_mole_fracs = np.array(
            [
                0.398247713,
                0.199123907,
                0.199123907,
                0.000398169,
                0.003982398,
                0.199123907,
                0.0,
            ]
        )
        # order: Ar, N2, O2, H, CH3, CH4
        for i in range(len(simulated_mole_fracs)):
            assert round(abs(simulated_mole_fracs[i] - expected_mole_fracs[i]), 6) == 0
