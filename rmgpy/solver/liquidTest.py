#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
import unittest

import numpy as np

import rmgpy
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.main import RMG
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import TerminationTime
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


################################################################################

class LiquidReactorCheck(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Here we choose a kinetic model consisting of the hydrogen abstraction
        reaction CH4 + C2H5 <=> CH3 + C2H6.

        Reset the loaded database
        """

        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

        Tlist = [300, 400, 500, 600, 800, 1000, 1500]
        cls.CH4 = Species(
            molecule=[Molecule().from_smiles("C")],
            thermo=ThermoData(
                Tdata=(Tlist, "K"),
                Cpdata=([8.615, 9.687, 10.963, 12.301, 14.841, 16.976, 20.528], "cal/(mol*K)"),
                H298=(-17.714, "kcal/mol"),
                S298=(44.472, "cal/(mol*K)"))
        )
        cls.CH3 = Species(
            molecule=[Molecule().from_smiles("[CH3]")],
            thermo=ThermoData(
                Tdata=(Tlist, "K"),
                Cpdata=([9.397, 10.123, 10.856, 11.571, 12.899, 14.055, 16.195], "cal/(mol*K)"),
                H298=(9.357, "kcal/mol"),
                S298=(45.174, "cal/(mol*K)"))
        )
        cls.C2H6 = Species(
            molecule=[Molecule().from_smiles("CC")],
            thermo=ThermoData(
                Tdata=(Tlist, "K"),
                Cpdata=([12.684, 15.506, 18.326, 20.971, 25.500, 29.016, 34.595], "cal/(mol*K)"),
                H298=(-19.521, "kcal/mol"),
                S298=(54.799, "cal/(mol*K)"))
        )
        cls.C2H5 = Species(
            molecule=[Molecule().from_smiles("C[CH2]")],
            thermo=ThermoData(
                Tdata=(Tlist, "K"),
                Cpdata=([11.635, 13.744, 16.085, 18.246, 21.885, 24.676, 29.107], "cal/(mol*K)"),
                H298=(29.496, "kcal/mol"),
                S298=(56.687, "cal/(mol*K)"))
        )

        cls.H2 = Species(
            molecule=[Molecule().from_smiles("[H][H]")],
            thermo=ThermoData(
                Tdata=(Tlist, "K"),
                Cpdata=([6.89, 6.97, 6.99, 7.01, 7.08, 7.22, 7.72], "cal/(mol*K)"),
                H298=(0, "kcal/mol"),
                S298=(31.23, "cal/(mol*K)"))
        )

        cls.T = 1000

        cls.file_dir = os.path.join(os.path.dirname(rmgpy.__file__), 'solver', 'files', 'liquid_phase_constSPC')

    def test_compute_flux(self):
        """
        Test the liquid batch reactor with a simple kinetic model. 
        """

        rxn1 = Reaction(
            reactants=[self.C2H6, self.CH3],
            products=[self.C2H5, self.CH4],
            kinetics=Arrhenius(A=(686.375 * 6, 'm^3/(mol*s)'), n=4.40721, Ea=(7.82799, 'kcal/mol'), T0=(298.15, 'K'))
        )

        core_species = [self.CH4, self.CH3, self.C2H6, self.C2H5]
        edge_species = []
        core_reactions = [rxn1]
        edge_reactions = []

        c0 = {self.C2H5: 0.1, self.CH3: 0.1, self.CH4: 0.4, self.C2H6: 0.4}

        rxn_system = LiquidReactor(self.T, c0, 1, termination=[])

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.array([10 ** (i / 10.0) for i in range(-130, -49)], np.float64)

        # Integrate to get the solution at each time point
        t, y, reaction_rates, species_rates = [], [], [], []
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxn_system.y.copy())
            reaction_rates.append(rxn_system.core_reaction_rates.copy())
            species_rates.append(rxn_system.core_species_rates.copy())

        # Convert the solution vectors to np arrays
        t = np.array(t, np.float64)
        reaction_rates = np.array(reaction_rates, np.float64)
        species_rates = np.array(species_rates, np.float64)

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            self.assertAlmostEqual(reaction_rates[i, 0], species_rates[i, 0], delta=1e-6 * reaction_rates[i, 0])
            self.assertAlmostEqual(reaction_rates[i, 0], -species_rates[i, 1], delta=1e-6 * reaction_rates[i, 0])
            self.assertAlmostEqual(reaction_rates[i, 0], -species_rates[i, 2], delta=1e-6 * reaction_rates[i, 0])
            self.assertAlmostEqual(reaction_rates[i, 0], species_rates[i, 3], delta=1e-6 * reaction_rates[i, 0])

        # Check that we've reached equilibrium 
        self.assertAlmostEqual(reaction_rates[-1, 0], 0.0, delta=1e-2)

    def test_jacobian(self):
        """
        Unit test for the jacobian function:
        Solve a reaction system and check if the analytical jacobian matches
        the finite difference jacobian.
        """

        core_species = [self.CH4, self.CH3, self.C2H6, self.C2H5, self.H2]
        edge_species = []
        num_core_species = len(core_species)
        c0 = {self.CH4: 0.2, self.CH3: 0.1, self.C2H6: 0.35, self.C2H5: 0.15, self.H2: 0.2}
        edge_reactions = []

        rxn_list = [
            Reaction(
                reactants=[self.C2H6],
                products=[self.CH3, self.CH3],
                kinetics=Arrhenius(
                    A=(686.375 * 6, '1/s'), n=4.40721, Ea=(7.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.CH3, self.CH3],
                products=[self.C2H6],
                kinetics=Arrhenius(
                    A=(686.375 * 6, 'm^3/(mol*s)'), n=4.40721, Ea=(7.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H6, self.CH3],
                products=[self.C2H5, self.CH4],
                kinetics=Arrhenius(
                    A=(46.375 * 6, 'm^3/(mol*s)'), n=3.40721, Ea=(6.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H5, self.CH4],
                products=[self.C2H6, self.CH3],
                kinetics=Arrhenius(
                    A=(46.375 * 6, 'm^3/(mol*s)'), n=3.40721, Ea=(6.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H5, self.CH4],
                products=[self.CH3, self.CH3, self.CH3],
                kinetics=Arrhenius(
                    A=(246.375 * 6, 'm^3/(mol*s)'), n=1.40721, Ea=(3.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.CH3, self.CH3, self.CH3],
                products=[self.C2H5, self.CH4],
                kinetics=Arrhenius(
                    A=(246.375 * 6, 'm^6/(mol^2*s)'), n=1.40721, Ea=(3.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H6, self.CH3, self.CH3],
                products=[self.C2H5, self.C2H5, self.H2],
                kinetics=Arrhenius(
                    A=(146.375 * 6, 'm^6/(mol^2*s)'), n=2.40721, Ea=(8.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H5, self.C2H5, self.H2],
                products=[self.C2H6, self.CH3, self.CH3],
                kinetics=Arrhenius(
                    A=(146.375 * 6, 'm^6/(mol^2*s)'), n=2.40721, Ea=(8.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H6, self.C2H6],
                products=[self.CH3, self.CH4, self.C2H5],
                kinetics=Arrhenius(
                    A=(1246.375 * 6, 'm^3/(mol*s)'), n=0.40721, Ea=(8.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.CH3, self.CH4, self.C2H5],
                products=[self.C2H6, self.C2H6],
                kinetics=Arrhenius(
                    A=(46.375 * 6, 'm^6/(mol^2*s)'), n=0.10721, Ea=(8.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
        ]

        # Analytical Jacobian for reaction 6
        def jacobian_rxn6(c, kf, kr, s):
            c1, c2, c3, c4 = c[s[1]], c[s[2]], c[s[3]], c[s[4]]
            jaco = np.zeros((5, 5))

            jaco[1, 1] = -4 * kf * c1 * c2
            jaco[1, 2] = -2 * kf * c1 * c1
            jaco[1, 3] = 4 * kr * c3 * c4
            jaco[1, 4] = 2 * kr * c3 * c3
            jaco[2, 1:] = 0.5 * jaco[1, 1:]
            jaco[3, 1:] = -jaco[1, 1:]
            jaco[4, 1:] = -0.5 * jaco[1, 1:]
            return jaco

        # Analytical Jacobian for reaction 7
        def jacobian_rxn7(c, kf, kr, s):
            c1, c2, c3, c4 = c[s[1]], c[s[2]], c[s[3]], c[s[4]]
            jaco = np.zeros((5, 5))

            jaco[1, 1] = -4 * kr * c1 * c2
            jaco[1, 2] = -2 * kr * c1 * c1
            jaco[1, 3] = 4 * kf * c3 * c4
            jaco[1, 4] = 2 * kf * c3 * c3
            jaco[2, 1:] = 0.5 * jaco[1, 1:]
            jaco[3, 1:] = -jaco[1, 1:]
            jaco[4, 1:] = -0.5 * jaco[1, 1:]
            return jaco

        for rxn_num, rxn in enumerate(rxn_list):
            core_reactions = [rxn]

            rxn_system0 = LiquidReactor(self.T, c0, 1, termination=[])
            rxn_system0.initialize_model(core_species, core_reactions, edge_species, edge_reactions)
            dydt0 = rxn_system0.residual(0.0, rxn_system0.y, np.zeros(rxn_system0.y.shape))[0]

            dN = .000001 * sum(rxn_system0.y)

            # Let the solver compute the jacobian
            solver_jacobian = rxn_system0.jacobian(0.0, rxn_system0.y, dydt0, 0.0)

            if rxn_num not in (6, 7):
                dydt = []
                for i in range(num_core_species):
                    rxn_system0.y[i] += dN
                    dydt.append(rxn_system0.residual(0.0, rxn_system0.y, np.zeros(rxn_system0.y.shape))[0])
                    rxn_system0.y[i] -= dN  # reset y

                # Compute the jacobian using finite differences
                jacobian = np.zeros((num_core_species, num_core_species))
                for i in range(num_core_species):
                    for j in range(num_core_species):
                        jacobian[i, j] = (dydt[j][i] - dydt0[i]) / dN
                        self.assertAlmostEqual(jacobian[i, j], solver_jacobian[i, j], delta=abs(1e-4 * jacobian[i, j]))
            # The forward finite difference is very unstable for reactions
            # 6 and 7. Use Jacobians calculated by hand instead.
            elif rxn_num == 6:
                kforward = rxn.get_rate_coefficient(self.T)
                kreverse = kforward / rxn.get_equilibrium_constant(self.T)
                jacobian = jacobian_rxn6(c0, kforward, kreverse, core_species)
                for i in range(num_core_species):
                    for j in range(num_core_species):
                        self.assertAlmostEqual(jacobian[i, j], solver_jacobian[i, j], delta=abs(1e-4 * jacobian[i, j]))
            elif rxn_num == 7:
                kforward = rxn.get_rate_coefficient(self.T)
                kreverse = kforward / rxn.get_equilibrium_constant(self.T)
                jacobian = jacobian_rxn7(c0, kforward, kreverse, core_species)
                for i in range(num_core_species):
                    for j in range(num_core_species):
                        self.assertAlmostEqual(jacobian[i, j], solver_jacobian[i, j], delta=abs(1e-4 * jacobian[i, j]))

    def test_compute_derivative(self):
        rxn_list = [
            Reaction(
                reactants=[self.C2H6],
                products=[self.CH3, self.CH3],
                kinetics=Arrhenius(
                    A=(686.375e6, '1/s'), n=4.40721, Ea=(7.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H6, self.CH3],
                products=[self.C2H5, self.CH4],
                kinetics=Arrhenius(
                    A=(46.375 * 6, 'm^3/(mol*s)'), n=3.40721, Ea=(6.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
            Reaction(
                reactants=[self.C2H6, self.CH3, self.CH3],
                products=[self.C2H5, self.C2H5, self.H2],
                kinetics=Arrhenius(
                    A=(146.375 * 6, 'm^6/(mol^2*s)'), n=2.40721, Ea=(8.82799, 'kcal/mol'), T0=(298.15, 'K'))
            ),
        ]

        core_species = [self.CH4, self.CH3, self.C2H6, self.C2H5, self.H2]
        edge_species = []
        core_reactions = rxn_list
        edge_reactions = []
        num_core_species = len(core_species)

        c0 = {self.CH4: 0.2, self.CH3: 0.1, self.C2H6: 0.35, self.C2H5: 0.15, self.H2: 0.2}

        rxn_system0 = LiquidReactor(self.T, c0, 1, termination=[])
        rxn_system0.initialize_model(core_species, core_reactions, edge_species, edge_reactions)
        dfdt0 = rxn_system0.residual(0.0, rxn_system0.y, np.zeros(rxn_system0.y.shape))[0]
        solver_dfdk = rxn_system0.compute_rate_derivative()
        # print 'Solver d(dy/dt)/dk'
        # print solver_dfdk

        integration_time = 1e-8

        model_settings = ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=0)
        simulator_settings = SimulatorSettings()

        rxn_system0.termination.append(TerminationTime((integration_time, 's')))

        rxn_system0.simulate(core_species, core_reactions, [], [], [], [],
                             model_settings=model_settings, simulator_settings=simulator_settings)

        y0 = rxn_system0.y

        dfdk = np.zeros((num_core_species, len(rxn_list)))  # d(dy/dt)/dk

        c0 = {self.CH4: 0.2, self.CH3: 0.1, self.C2H6: 0.35, self.C2H5: 0.15, self.H2: 0.2}

        for i in range(len(rxn_list)):
            k0 = rxn_list[i].get_rate_coefficient(self.T)
            rxn_list[i].kinetics.A.value_si = rxn_list[i].kinetics.A.value_si * (1 + 1e-3)
            dk = rxn_list[i].get_rate_coefficient(self.T) - k0

            rxn_system = LiquidReactor(self.T, c0, 1, termination=[])
            rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

            dfdt = rxn_system.residual(0.0, rxn_system.y, np.zeros(rxn_system.y.shape))[0]
            dfdk[:, i] = (dfdt - dfdt0) / dk

            rxn_system.termination.append(TerminationTime((integration_time, 's')))
            model_settings = ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=0)
            simulator_settings = SimulatorSettings()
            rxn_system.simulate(core_species, core_reactions, [], [], [], [],
                                model_settings=model_settings, simulator_settings=simulator_settings)

            rxn_list[i].kinetics.A.value_si = rxn_list[i].kinetics.A.value_si / (1 + 1e-3)  # reset A factor

        for i in range(num_core_species):
            for j in range(len(rxn_list)):
                self.assertAlmostEqual(dfdk[i, j], solver_dfdk[i, j], delta=abs(1e-3 * dfdk[i, j]))

    def test_store_constant_species_names(self):
        """
        Test if (i) constant species names are stored in reactor attributes and
        (ii) if attributes are not mix/equal for multiple conditions generation.
        """

        c0 = {self.C2H5: 0.1, self.CH3: 0.1, self.CH4: 0.4, self.C2H6: 0.4}
        temp = 1000

        # set up the liquid phase reactor 1
        termination_conversion = []
        termination_time = None
        sensitivity = []
        sensitivity_threshold = 0.001
        constant_species = ["CH4", "C2H6"]
        sens_conds = None
        rxn_system1 = LiquidReactor(temp, c0, 4, termination_conversion, sensitivity, sensitivity_threshold, sens_conds,
                                    constant_species)

        # set up the liquid phase reactor 2
        constant_species = ["O2", "H2O"]
        rxn_system2 = LiquidReactor(temp, c0, 4, termination_conversion, sensitivity, sensitivity_threshold, sens_conds,
                                    constant_species)
        for reactor in [rxn_system1, rxn_system2]:
            self.assertIsNotNone(reactor.const_spc_names)

        # check if Constant species are different in each liquid system
        for spc in rxn_system1.const_spc_names:
            for spc2 in rxn_system2.const_spc_names:
                self.assertIsNot(spc, spc2, 'Constant species declared in two different reactors seem mixed. '
                                            'Species "{0}" appears in both systems and should be.'.format(spc))

    def test_liquid_input_reading(self):
        """
        Check if constant concentration condition is well handled. 
        From input file reading to information storage in liquid reactor object.
        """
        rmg = RMG()
        rmg.input_file = os.path.join(self.file_dir, 'input.py')
        rmg.initialize()

        for index, reactionSystem in enumerate(rmg.reaction_systems):
            self.assertIsNotNone(reactionSystem.const_spc_names,
                                 'Reactor should contain constant species name and indices after few steps')
            self.assertIsNotNone(reactionSystem.const_spc_indices,
                                 'Reactor should contain constant species indices in the core species array')
            self.assertIs(reactionSystem.const_spc_names[0],
                          rmg.reaction_model.core.species[reactionSystem.const_spc_indices[0]].label,
                          'The constant species name from the reaction model and constantSPCnames should be equal')

    def test_corespecies_rate(self):
        """
        Test if a specific core species rate is equal to 0 over time.
        """

        c0 = {self.C2H5: 0.1, self.CH3: 0.1, self.CH4: 0.4, self.C2H6: 0.4}
        rxn1 = Reaction(
            reactants=[self.C2H6, self.CH3],
            products=[self.C2H5, self.CH4],
            kinetics=Arrhenius(A=(686.375 * 6, 'm^3/(mol*s)'), n=4.40721, Ea=(7.82799, 'kcal/mol'), T0=(298.15, 'K'))
        )

        core_species = [self.CH4, self.CH3, self.C2H6, self.C2H5]
        edge_species = []
        core_reactions = [rxn1]
        edge_reactions = []
        sensitivity = []
        termination_conversion = []
        sensitivity_threshold = 0.001
        const_species = ["CH4"]
        sens_conds = {self.C2H5: 0.1, self.CH3: 0.1, self.CH4: 0.4, self.C2H6: 0.4, 'T': self.T}

        rxn_system = LiquidReactor(self.T, c0, 1, termination_conversion, sensitivity, sensitivity_threshold,
                                   const_spc_names=const_species, sens_conditions=sens_conds)
        # The test regarding the writing of constantSPCindices from input file is check with the previous test.
        rxn_system.const_spc_indices = [0]

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.array([10 ** (i / 10.0) for i in range(-130, -49)], np.float64)

        # Integrate to get the solution at each time point
        t, y, reaction_rates, species_rates = [], [], [], []
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            self.assertEqual(rxn_system.core_species_rates[0], 0,
                             "Core species rate has to be equal to 0 for species hold constant. "
                             "Here it is equal to {0}".format(rxn_system.core_species_rates[0]))

    @classmethod
    def tearDownClass(cls):
        """
        Reset the database & liquid parameters for solution
        """
        global diffusion_limiter
        from rmgpy.kinetics.diffusionLimited import diffusion_limiter
        diffusion_limiter.enabled = False

        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

        os.remove(os.path.join(cls.file_dir, 'restart_from_seed.py'))
