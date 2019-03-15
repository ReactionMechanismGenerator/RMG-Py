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
This script contains unit tests of the :mod:`arkane.common` module.
"""

import unittest
import numpy
import os
import shutil
import logging

import rmgpy
import rmgpy.constants as constants
from rmgpy.species import Species, TransitionState
from rmgpy.quantity import ScalarQuantity
from rmgpy.thermo import NASA

from arkane import Arkane, input
from arkane.common import ArkaneSpecies, get_element_mass
from arkane.statmech import InputError, StatMechJob
from arkane.input import jobList

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
        v_list = numpy.array(v_list, numpy.float64)
        v_diff = (v_list[0] - numpy.min(v_list)) * constants.E_h * constants.Na / 1000
        self.assertAlmostEqual(v_diff / 2.7805169838282797, 1, 5)


class TestArkaneJob(unittest.TestCase):
    """
    Contains unit tests of the Arkane module and its interactions with other RMG modules.
    """
    @classmethod
    def setUp(self):
        arkane = Arkane()
        jobList = arkane.loadInputFile(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                          'data', 'methoxy.py'))
        pdepjob = jobList[-1]
        self.kineticsjob = jobList[0]
        pdepjob.activeJRotor = True
        network = pdepjob.network
        self.Nisom = len(network.isomers)
        self.Nreac = len(network.reactants)
        self.Nprod = len(network.products)
        self.Npath = len(network.pathReactions)
        self.PathReaction2 = network.pathReactions[2]
        self.TminValue = pdepjob.Tmin.value
        self.Tmaxvalue = pdepjob.Tmax.value
        self.TmaxUnits = pdepjob.Tmax.units
        self.TlistValue = pdepjob.Tlist.value
        self.PminValue = pdepjob.Pmin.value
        self.Pcount = pdepjob.Pcount
        self.Tcount = pdepjob.Tcount
        self.GenTlist = pdepjob.generateTemperatureList()
        self.PlistValue = pdepjob.Plist.value
        self.maximumGrainSizeValue = pdepjob.maximumGrainSize.value
        self.method = pdepjob.method
        self.rmgmode = pdepjob.rmgmode

    # test Arkane's interactions with the network module
    def testNisom(self):
        """
        Test the number of isomers identified.
        """
        self.assertEqual(self.Nisom, 2, msg=None)

    def testNreac(self):
        """
        Test the number of reactants identified.
        """
        self.assertEqual(self.Nreac, 1, msg=None)

    def testNprod(self):
        """
        Test the number of products identified.
        """
        self.assertEqual(self.Nprod, 1, msg=None)

    def testNpathReactions(self):
        """
        Test the whether or not RMG mode is turned on.
        """
        self.assertEqual(self.Npath, 3, msg=None)

    def testPathReactions(self):
        """
        Test a path reaction label
        """
        self.assertEqual(str(self.PathReaction2), 'CH2OH <=> methoxy', msg=None)

    # test Arkane's interactions with the pdep module
    def testTemperaturesUnits(self):
        """
        Test the Temperature Units.
        """
        self.assertEqual(str(self.TmaxUnits), 'K', msg=None)

    def testTemperaturesValue(self):
        """
        Test the temperature value.
        """
        self.assertEqual(self.TminValue, 450.0, msg=None)

    def testTemperaturesList(self):
        """
        Test the temperature list.
        """
        self.assertEqual(numpy.array_equal(self.TlistValue, numpy.array([450, 500, 678, 700])), True, msg=None)

    def testPminValue(self):
        """
        Test the minimum pressure value.
        """
        self.assertEqual("%0.7f" % self.PminValue, str(0.0101325), msg=None)

    def testPcount(self):
        """
        Test the number pressures specified.
        """
        self.assertEqual(self.Pcount, 7, msg=None)

    def testTcount(self):
        """
        Test the number temperatures specified.
        """
        self.assertEqual(self.Tcount, 4, msg=None)

    def testPressureList(self):
        """
        Test the pressure list.
        """
        self.assertEqual(numpy.array_equal(self.PlistValue, numpy.array([0.01, 0.1, 1, 3, 10, 100, 1000])), True, msg=None)

    def testGenerateTemperatureList(self):
        """
        Test the generated temperature list.
        """
        self.assertEqual(list(self.GenTlist), [450.0, 500.0, 678.0, 700.0], msg=None)

    def testmaximumGrainSizeValue(self):
        """
        Test the max grain size value.
        """
        self.assertEqual(self.maximumGrainSizeValue, 0.5, msg=None)

    def testMethod(self):
        """
        Test the master equation solution method chosen.
        """
        self.assertEqual(self.method, 'modified strong collision', msg=None)

    def testRmgmode(self):
        """
        Test the whether or not RMG mode is turned on.
        """
        self.assertEqual(self.rmgmode, False, msg=None)

    # Test Arkane's interactions with the kinetics module
    def testCalculateTSTRateCoefficient(self):
        """
        Test the calculation of the high-pressure limit rate coef for one of the kinetics jobs at Tmin and Tmax.
        """
        self.assertEqual("%0.7f" % self.kineticsjob.reaction.calculateTSTRateCoefficient(self.TminValue),
                         str(46608.5904933), msg=None)
        self.assertEqual("%0.5f" % self.kineticsjob.reaction.calculateTSTRateCoefficient(self.Tmaxvalue),
                         str(498796.64535), msg=None)

    def testTunneling(self):
        """
        Test the whether or not tunneling has been included in a specific kinetics job.
        """
        self.assertEqual(self.kineticsjob.reaction.transitionState.tunneling, None, msg=None)


class TestArkaneInput(unittest.TestCase):
    """
    Contains unit tests for loading and processing Arkane input files.
    """
    @classmethod
    def setUp(self):
        """Preparation for all unit tests in this class."""
        self.directory = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)), 'examples', 'arkane')
        self.modelChemistry = "cbs-qb3"
        self.frequencyScaleFactor = 0.99
        self.useHinderedRotors = False
        self.useBondCorrections = True

    def testSpecies(self):
        """Test loading of species input file."""
        spec = input.species('C2H4', os.path.join(self.directory, 'species', 'C2H4', 'ethene.py'))
        self.assertTrue(isinstance(spec, Species))
        self.assertEqual(len(spec.molecule), 0)

    def testSpeciesStatmech(self):
        """Test loading of statmech job from species input file."""
        job = jobList[-1]
        self.assertTrue(isinstance(job, StatMechJob))
        job.modelChemistry = self.modelChemistry
        job.frequencyScaleFactor = self.frequencyScaleFactor
        job.includeHinderedRotors = self.useHinderedRotors
        job.applyBondEnergyCorrections = self.useBondCorrections
        job.load()
        self.assertTrue(isinstance(job.species.props['elementCounts'], dict))
        self.assertEqual(job.species.props['elementCounts']['C'], 2)
        self.assertEqual(job.species.props['elementCounts']['H'], 4)

    def testSpeciesThermo(self):
        """Test thermo job execution for species from separate input file."""
        input.thermo('C2H4', 'NASA')
        job = jobList[-1]
        filepath = os.path.join(self.directory, 'reactions', 'H+C2H4=C2H5', 'output.py')
        job.execute(outputFile=filepath)
        self.assertTrue(os.path.isfile(os.path.join(os.path.dirname(filepath), 'output.py')))
        self.assertTrue(os.path.isfile(os.path.join(os.path.dirname(filepath), 'chem.inp')))
        os.remove(os.path.join(os.path.dirname(filepath), 'output.py'))
        os.remove(os.path.join(os.path.dirname(filepath), 'chem.inp'))

    def testTransitionState(self):
        """Test loading of transition state input file."""
        ts = input.transitionState('TS', os.path.join(self.directory, 'reactions', 'H+C2H4=C2H5', 'TS.py'))
        self.assertTrue(isinstance(ts, TransitionState))
    
    def testTransitionStateStatmech(self):
        """Test loading of statmech job from transition state input file."""
        job = jobList[-1]
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
    def setUp(self):
        arkane = Arkane()
        self.job_list = arkane.loadInputFile(os.path.join(os.path.dirname(os.path.abspath(__file__)),
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

        if os.path.exists(cls.dump_yaml_file):
            logging.debug('removing existing yaml file {0} before running tests'.format(cls.dump_yaml_file))
            os.remove(cls.dump_yaml_file)

    def test_dump_yaml(self):
        """
        Test properly dumping the ArkaneSpecies object and respective sub-objects
        """
        jobList = self.arkane.loadInputFile(self.dump_input_path)
        for job in jobList:
            job.execute(outputFile=self.dump_output_file)
        self.assertTrue(os.path.isfile(self.dump_yaml_file))

    def test_load_yaml(self):
        """
        Test properly loading the ArkaneSpecies object and respective sub-objects
        """
        jobList = self.arkane.loadInputFile(self.load_input_path)
        for job in jobList:
            job.execute(outputFile=self.load_output_file)
        arkane_spc = jobList[0].arkane_species
        self.assertIsInstance(arkane_spc, ArkaneSpecies)  # checks make_object
        self.assertIsInstance(arkane_spc.molecular_weight, ScalarQuantity)
        self.assertIsInstance(arkane_spc.thermo, NASA)
        self.assertNotEqual(arkane_spc.author, '')
        self.assertEqual(arkane_spc.inchi, 'InChI=1S/C2H6/c1-2/h1-2H3')
        self.assertEqual(arkane_spc.smiles, 'CC')
        self.assertTrue('8 H u0 p0 c0 {2,S}' in arkane_spc.adjacency_list)
        self.assertEqual(arkane_spc.label, 'C2H6')
        self.assertEqual(arkane_spc.frequency_scale_factor, 0.99)  # checks float conversion
        self.assertFalse(arkane_spc.use_bond_corrections)
        self.assertAlmostEqual(arkane_spc.conformer.modes[2].frequencies.value_si[0], 818.91718, 4)  # HarmonicOsc.

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
                    if name in cls.files_to_delete or\
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
        """Test that the correct mass/number/isotop is returned from get_element_mass"""
        self.assertEquals(get_element_mass(1), (1.00782503224, 1))  # test input by integer
        self.assertEquals(get_element_mass('Si'), (27.97692653465, 14))  # test string input and most common isotope
        self.assertEquals(get_element_mass('C', 13), (13.00335483507, 6))  # test specific isotope
        self.assertEquals(get_element_mass('Bk'), (247.0703073, 97))  # test a two-element array (no isotope data)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
