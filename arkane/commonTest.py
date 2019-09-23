#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

"""
This module contains unit tests of the :mod:`arkane.common` module.
"""

import logging
import os
import shutil
import unittest

import numpy as np

import rmgpy
import rmgpy.constants as constants
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.quantity import ScalarQuantity
from rmgpy.species import Species, TransitionState
from rmgpy.thermo import NASA, ThermoData

from arkane import Arkane, input
from arkane.common import ArkaneSpecies, get_element_mass
from arkane.input import job_list
from arkane.statmech import InputError, StatMechJob

################################################################################


class CommonTest(unittest.TestCase):
    """
    Contains unit tests of Arkane's common functions.
    """

    def test_check_conformer_energy(self):
        """
        test the check_conformer_energy function with an list of energies.
        """
        v_list = [-272.2779012225, -272.2774933703, -272.2768397635, -272.2778432059, -272.278645477, -272.2789602654,
                  -272.2788749196, -272.278496709, -272.2779350675, -272.2777008843, -272.2777167286, -272.2780937643,
                  -272.2784838846, -272.2788050464, -272.2787865352, -272.2785091607, -272.2779977452, -272.2777957743,
                  -272.2779134906, -272.2781827547, -272.278443339, -272.2788244214, -272.2787748749]
        v_list = np.array(v_list, np.float64)
        v_diff = (v_list[0] - np.min(v_list)) * constants.E_h * constants.Na / 1000
        self.assertAlmostEqual(v_diff / 2.7805169838282797, 1, 5)


class TestArkaneJob(unittest.TestCase):
    """
    Contains unit tests of the Arkane module and its interactions with other RMG modules.
    """

    @classmethod
    def setUp(cls):
        """A method that is run before each unit test in this class"""
        arkane = Arkane()
        job_list = arkane.load_input_file(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                       'data', 'methoxy.py'))
        pdepjob = job_list[-1]
        cls.kineticsjob = job_list[0]
        pdepjob.active_j_rotor = True
        network = pdepjob.network
        cls.Nisom = len(network.isomers)
        cls.Nreac = len(network.reactants)
        cls.Nprod = len(network.products)
        cls.Npath = len(network.path_reactions)
        cls.PathReaction2 = network.path_reactions[2]
        cls.TminValue = pdepjob.Tmin.value
        cls.Tmaxvalue = pdepjob.Tmax.value
        cls.TmaxUnits = pdepjob.Tmax.units
        cls.TlistValue = pdepjob.Tlist.value
        cls.PminValue = pdepjob.Pmin.value
        cls.Pcount = pdepjob.Pcount
        cls.Tcount = pdepjob.Tcount
        cls.GenTlist = pdepjob.generate_T_list()
        cls.PlistValue = pdepjob.Plist.value
        cls.maximum_grain_size_value = pdepjob.maximum_grain_size.value
        cls.method = pdepjob.method
        cls.rmgmode = pdepjob.rmgmode

    # test Arkane's interactions with the network module
    def test_num_isom(self):
        """
        Test the number of isomers identified.
        """
        self.assertEqual(self.Nisom, 2)

    def test_num_reac(self):
        """
        Test the number of reactants identified.
        """
        self.assertEqual(self.Nreac, 1)

    def test_num_prod(self):
        """
        Test the number of products identified.
        """
        self.assertEqual(self.Nprod, 1)

    def test_n_path_reactions(self):
        """
        Test the whether or not RMG mode is turned on.
        """
        self.assertEqual(self.Npath, 3)

    def test_path_reactions(self):
        """
        Test a path reaction label
        """
        self.assertEqual(str(self.PathReaction2), 'CH2OH <=> methoxy')

    # test Arkane's interactions with the pdep module
    def test_temperatures_units(self):
        """
        Test the Temperature Units.
        """
        self.assertEqual(str(self.TmaxUnits), 'K')

    def test_temperatures_value(self):
        """
        Test the temperature value.
        """
        self.assertEqual(self.TminValue, 450.0)

    def test_temperatures_list(self):
        """
        Test the temperature list.
        """
        self.assertTrue(np.array_equal(self.TlistValue, np.array([450, 500, 678, 700])))

    def test_min_pressure_value(self):
        """
        Test the minimum pressure value.
        """
        self.assertEqual("%0.7f" % self.PminValue, str(0.0101325))

    def test_pressure_count(self):
        """
        Test the number pressures specified.
        """
        self.assertEqual(self.Pcount, 7)

    def test_temperature_count(self):
        """
        Test the number temperatures specified.
        """
        self.assertEqual(self.Tcount, 4)

    def test_pressure_list(self):
        """
        Test the pressure list.
        """
        self.assertTrue(np.array_equal(self.PlistValue, np.array([0.01, 0.1, 1, 3, 10, 100, 1000])))

    def test_generate_temperature_list(self):
        """
        Test the generated temperature list.
        """
        self.assertEqual(list(self.GenTlist), [450.0, 500.0, 678.0, 700.0])

    def test_maximum_grain_size_value(self):
        """
        Test the max grain size value.
        """
        self.assertEqual(self.maximum_grain_size_value, 0.5)

    def test_method(self):
        """
        Test the master equation solution method chosen.
        """
        self.assertEqual(self.method, 'modified strong collision')

    def test_rmg_mode(self):
        """
        Test the whether or not RMG mode is turned on.
        """
        self.assertEqual(self.rmgmode, False)

    # Test Arkane's interactions with the kinetics module
    def test_calculate_tst_rate_coefficient(self):
        """
        Test the calculation of the high-pressure limit rate coef for one of the kinetics jobs at Tmin and Tmax.
        """
        self.assertEqual("%0.7f" % self.kineticsjob.reaction.calculate_tst_rate_coefficient(self.TminValue),
                         str(46608.5904933))
        self.assertEqual("%0.5f" % self.kineticsjob.reaction.calculate_tst_rate_coefficient(self.Tmaxvalue),
                         str(498796.64535))

    def test_tunneling(self):
        """
        Test the whether or not tunneling has been included in a specific kinetics job.
        """
        self.assertEqual(self.kineticsjob.reaction.transition_state.tunneling, None)


class TestArkaneInput(unittest.TestCase):
    """
    Contains unit tests for loading and processing Arkane input files.
    """

    @classmethod
    def setUp(cls):
        """Preparation for all unit tests in this class."""
        cls.directory = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)), 'examples', 'arkane')
        cls.modelChemistry = "cbs-qb3"
        cls.frequencyScaleFactor = 0.99
        cls.useHinderedRotors = False
        cls.useBondCorrections = True

    def test_species(self):
        """Test loading of species input file."""
        spec = input.species('C2H4', os.path.join(self.directory, 'species', 'C2H4', 'ethene.py'))
        self.assertTrue(isinstance(spec, Species))
        self.assertEqual(len(spec.molecule), 0)

    def test_species_statmech(self):
        """Test loading of statmech job from species input file."""
        job = job_list[-1]
        self.assertTrue(isinstance(job, StatMechJob))
        job.modelChemistry = self.modelChemistry
        job.frequencyScaleFactor = self.frequencyScaleFactor
        job.includeHinderedRotors = self.useHinderedRotors
        job.applyBondEnergyCorrections = self.useBondCorrections
        job.load()
        self.assertTrue(isinstance(job.species.props['element_counts'], dict))
        self.assertEqual(job.species.props['element_counts']['C'], 2)
        self.assertEqual(job.species.props['element_counts']['H'], 4)

    def test_species_thermo(self):
        """Test thermo job execution for species from separate input file."""
        input.thermo('C2H4', 'NASA')
        job = job_list[-1]
        filepath = os.path.join(self.directory, 'reactions', 'H+C2H4=C2H5')
        job.execute(output_directory=filepath)
        self.assertTrue(os.path.isfile(os.path.join(filepath, 'output.py')))
        self.assertTrue(os.path.isfile(os.path.join(filepath, 'chem.inp')))
        os.remove(os.path.join(filepath, 'output.py'))
        os.remove(os.path.join(filepath, 'chem.inp'))

    def test_transition_state(self):
        """Test loading of transition state input file."""
        ts = input.transitionState('TS', os.path.join(self.directory, 'reactions', 'H+C2H4=C2H5', 'TS.py'))
        self.assertTrue(isinstance(ts, TransitionState))

    def test_transition_state_statmech(self):
        """Test loading of statmech job from transition state input file."""
        job = job_list[-1]
        self.assertTrue(isinstance(job, StatMechJob))
        job.modelChemistry = self.modelChemistry
        job.frequencyScaleFactor = self.frequencyScaleFactor
        job.includeHinderedRotors = self.useHinderedRotors
        job.applyBondEnergyCorrections = self.useBondCorrections
        job.load()


class TestStatmech(unittest.TestCase):
    """
    Contains unit tests of statmech.py
    """

    @classmethod
    def setUp(cls):
        """A method that is run before each unit test in this class"""
        arkane = Arkane()
        cls.job_list = arkane.load_input_file(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                           'data', 'Benzyl', 'input.py'))

    def test_gaussian_log_file_error(self):
        """Test that the proper error is raised if gaussian geometry and frequency file paths are the same"""
        job = self.job_list[-2]
        self.assertTrue(isinstance(job, StatMechJob))
        with self.assertRaises(InputError):
            job.load()


class TestArkaneSpecies(unittest.TestCase):
    """
    Contains YAML dump and load unit tests for :class:ArkaneSpecies
    """

    @classmethod
    def setUpClass(cls):
        """
        A method that is run ONCE before all unit tests in this class.
        """
        cls.arkane = Arkane()
        path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)),
                            'examples', 'arkane', 'species')
        cls.dump_path = os.path.join(path, 'C2H6')
        cls.dump_input_path = os.path.join(cls.dump_path, 'input.py')
        cls.dump_output_file = os.path.join(cls.dump_path, 'output.py')
        cls.dump_yaml_file = os.path.join(cls.dump_path, 'species', 'C2H6.yml')

        cls.load_path = os.path.join(path, 'C2H6_from_yaml')
        cls.load_input_path = os.path.join(cls.load_path, 'input.py')
        cls.load_output_file = os.path.join(cls.load_path, 'output.py')

        cls.data_path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)), 'arkane', 'data')

        if os.path.exists(cls.dump_yaml_file):
            logging.debug('removing existing yaml file {0} before running tests'.format(cls.dump_yaml_file))
            os.remove(cls.dump_yaml_file)

    def test_dump_yaml(self):
        """
        Test properly dumping the ArkaneSpecies object and respective sub-objects
        """
        job_list = self.arkane.load_input_file(self.dump_input_path)
        for job in job_list:
            job.execute(output_directory=self.dump_path)
        self.assertTrue(os.path.isfile(self.dump_output_file))

    def test_create_and_load_yaml(self):
        """
        Test properly loading the ArkaneSpecies object and respective sub-objects
        """
        # Create YAML file by running Arkane
        job_list = self.arkane.load_input_file(self.dump_input_path)
        for job in job_list:
            job.execute(output_directory=self.dump_path)

        # Load in newly created YAML file
        arkane_spc_old = job_list[0].arkane_species
        arkane_spc = ArkaneSpecies.__new__(ArkaneSpecies)
        arkane_spc.load_yaml(path=os.path.join(self.dump_path, 'species', arkane_spc_old.label + '.yml'))

        self.assertIsInstance(arkane_spc, ArkaneSpecies)  # checks make_object
        self.assertIsInstance(arkane_spc.molecular_weight, ScalarQuantity)
        self.assertIsInstance(arkane_spc.thermo, NASA)
        self.assertNotEqual(arkane_spc.author, '')
        self.assertEqual(arkane_spc.inchi, 'InChI=1S/C2H6/c1-2/h1-2H3')
        self.assertEqual(arkane_spc.inchi_key, 'OTMSDBZUPAUEDD-UHFFFAOYSA-N')
        self.assertEqual(arkane_spc.smiles, 'CC')
        self.assertTrue('8 H u0 p0 c0 {2,S}' in arkane_spc.adjacency_list)
        self.assertEqual(arkane_spc.label, 'C2H6')
        self.assertEqual(arkane_spc.frequency_scale_factor, 0.99 * 1.014)  # checks float conversion
        self.assertFalse(arkane_spc.use_bond_corrections)
        self.assertAlmostEqual(arkane_spc.conformer.modes[2].frequencies.value_si[0], 830.38202, 4)  # HarmonicOsc.
        self.assertIsInstance(arkane_spc.energy_transfer_model, SingleExponentialDown)
        self.assertFalse(arkane_spc.is_ts)
        self.assertEqual(arkane_spc.level_of_theory, 'cbs-qb3')
        self.assertIsInstance(arkane_spc.thermo_data, ThermoData)
        self.assertTrue(arkane_spc.use_hindered_rotors)
        self.assertIsInstance(arkane_spc.chemkin_thermo_string, str)
        expected_xyz = """8
C2H6
C       0.00075400    0.00119300    0.00055200
H       0.00074000    0.00117100    1.09413800
H       1.04376600    0.00117100   -0.32820200
H      -0.44760300    0.94289500   -0.32825300
C      -0.76014200   -1.20389600   -0.55748300
H      -0.76012800   -1.20387400   -1.65106900
H      -0.31178500   -2.14559800   -0.22867800
H      -1.80315400   -1.20387400   -0.22872900"""
        self.assertEqual(arkane_spc.xyz, expected_xyz)

    def test_load_existing_yaml(self):
        """
        Test that existing Arkane YAML files can still be loaded
        """
        # Load in YAML file
        arkane_spc = ArkaneSpecies.__new__(ArkaneSpecies)
        arkane_spc.load_yaml(path=os.path.join(self.load_path, 'C2H6.yml'))

        self.assertIsInstance(arkane_spc, ArkaneSpecies)  # checks make_object
        self.assertIsInstance(arkane_spc.molecular_weight, ScalarQuantity)
        self.assertIsInstance(arkane_spc.thermo, NASA)
        self.assertNotEqual(arkane_spc.author, '')
        self.assertEqual(arkane_spc.inchi, 'InChI=1S/C2H6/c1-2/h1-2H3')
        self.assertEqual(arkane_spc.inchi_key, 'OTMSDBZUPAUEDD-UHFFFAOYSA-N')
        self.assertEqual(arkane_spc.smiles, 'CC')
        self.assertTrue('8 H u0 p0 c0 {2,S}' in arkane_spc.adjacency_list)
        self.assertEqual(arkane_spc.label, 'C2H6')
        self.assertEqual(arkane_spc.frequency_scale_factor, 0.99)  # checks float conversion
        self.assertFalse(arkane_spc.use_bond_corrections)
        self.assertAlmostEqual(arkane_spc.conformer.modes[2].frequencies.value_si[0], 818.91718, 4)  # HarmonicOsc.
        self.assertIsInstance(arkane_spc.energy_transfer_model, SingleExponentialDown)
        self.assertFalse(arkane_spc.is_ts)
        self.assertTrue(arkane_spc.use_hindered_rotors)
        self.assertTrue('C 7.54e-14 1.193e-13 5.52e-14' in arkane_spc.xyz)
        self.assertIsInstance(arkane_spc.chemkin_thermo_string, str)

    def test_loading_different_versions_of_yaml(self):
        """Test loading a YAML file generated by RMG v 2.4.1 and by a more recent version"""
        arkane_spc_v_241 = ArkaneSpecies.__new__(ArkaneSpecies)
        arkane_spc_v_241.load_yaml(path=os.path.join(self.data_path, 'vinoxy_v_2.4.1.yml'))
        self.assertIsInstance(arkane_spc_v_241, ArkaneSpecies)  # checks make_object
        self.assertEqual(arkane_spc_v_241.conformer.spin_multiplicity, 2)

        arkane_current = ArkaneSpecies.__new__(ArkaneSpecies)
        arkane_current.load_yaml(path=os.path.join(self.data_path, 'vinoxy_current.yml'))
        self.assertIsInstance(arkane_current, ArkaneSpecies)  # checks make_object
        self.assertEqual(arkane_current.conformer.spin_multiplicity, 2)

    @classmethod
    def tearDownClass(cls):
        """
        A method that is run ONCE after all unit tests in this class.
        """
        path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)),
                            'examples', 'arkane', 'species')
        cls.dump_path = os.path.join(path, 'C2H6')
        cls.load_path = os.path.join(path, 'C2H6_from_yaml')
        cls.extensions_to_delete = ['pdf', 'txt', 'inp', 'csv']
        cls.files_to_delete = ['arkane.log', 'output.py']
        cls.files_to_keep = ['C2H6.yml']
        for path in [cls.dump_path, cls.load_path]:
            for name in os.listdir(path):
                item_path = os.path.join(path, name)
                if os.path.isfile(item_path):
                    extension = name.split('.')[-1]
                    if name in cls.files_to_delete or \
                            (extension in cls.extensions_to_delete and name not in cls.files_to_keep):
                        os.remove(item_path)
                else:
                    # This is a sub-directory. remove.
                    shutil.rmtree(item_path)


class TestGetMass(unittest.TestCase):
    """
    Contains unit tests of common.py
    """

    def test_get_mass(self):
        """Test that the correct mass/number/isotope is returned from get_element_mass"""
        self.assertEquals(get_element_mass(1), (1.00782503224, 1))  # test input by integer
        self.assertEquals(get_element_mass('Si'), (27.97692653465, 14))  # test string input and most common isotope
        self.assertEquals(get_element_mass('C', 13), (13.00335483507, 6))  # test specific isotope
        self.assertEquals(get_element_mass('Bk'), (247.0703073, 97))  # test a two-element array (no isotope data)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
