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

import logging
import os
import shutil
import unittest
from unittest.mock import patch

from nose.plugins.attrib import attr

import pandas as pd

from rmgpy.rmg.main import RMG, initialize_log, make_profile_graph
from rmgpy.rmg.main import RMG_Memory
from rmgpy import get_path
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.model import CoreEdgeReactionModel

###################################################

originalPath = get_path()


@attr('functional')
class TestMain(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, 'rmg', 'test_data', 'mainTest')
        cls.outputDir = 'output'
        cls.databaseDirectory = settings['database.directory']

        cls.seedKinetics = os.path.join(cls.databaseDirectory, 'kinetics', 'libraries', 'testSeed')
        cls.seedKineticsEdge = os.path.join(cls.databaseDirectory, 'kinetics', 'libraries', 'testSeed_edge')

        os.mkdir(os.path.join(cls.testDir, cls.outputDir))

        cls.rmg = RMG(input_file=os.path.join(cls.testDir, 'input.py'),
                      output_directory=os.path.join(cls.testDir, cls.outputDir))

        cls.rmg.execute()

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(os.path.join(cls.testDir, cls.outputDir))

        # Delete the seed libraries created in database
        shutil.rmtree(cls.seedKinetics)
        shutil.rmtree(cls.seedKineticsEdge)

    def test_rmg_execute(self):
        """Test that RMG.execute completed successfully."""
        self.assertIsInstance(self.rmg.database, RMGDatabase)
        self.assertTrue(self.rmg.done)

    def test_rmg_increases_reactions(self):
        """Test that RMG.execute increases reactions and species."""
        self.assertTrue(len(self.rmg.reaction_model.core.reactions) > 0)
        self.assertTrue(len(self.rmg.reaction_model.core.species) > 1)
        self.assertTrue(len(self.rmg.reaction_model.edge.reactions) > 0)
        self.assertTrue(len(self.rmg.reaction_model.edge.species) > 0)

    def test_rmg_seed_mechanism_creation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seed_dir = os.path.join(self.testDir, self.outputDir, 'seed')
        self.assertTrue(os.path.exists)

        self.assertTrue(os.path.exists(os.path.join(seed_dir, 'seed')))  # kinetics library folder made

        self.assertTrue(os.path.exists(os.path.join(seed_dir, 'seed', 'dictionary.txt')))  # dictionary file made
        self.assertTrue(os.path.exists(os.path.join(seed_dir, 'seed', 'reactions.py')))  # reactions file made

    def test_rmg_seed_edge_mechanism_creation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seed_dir = os.path.join(self.testDir, self.outputDir, 'seed')
        self.assertTrue(os.path.exists)

        self.assertTrue(os.path.exists(os.path.join(seed_dir, 'seed_edge')))  # kinetics library folder made

        self.assertTrue(os.path.exists(os.path.join(seed_dir, 'seed_edge', 'dictionary.txt')))  # dictionary file made
        self.assertTrue(os.path.exists(os.path.join(seed_dir, 'seed_edge', 'reactions.py')))  # reactions file made

    def test_rmg_seed_library_creation(self):
        """Test that seed mechanisms are created in the correct database locations."""
        self.assertTrue(os.path.exists(self.seedKinetics))

    def test_rmg_seed_edge_library_creation(self):
        """Test that edge seed mechanisms are created in the correct database locations."""
        self.assertTrue(os.path.exists(self.seedKinetics))

    def test_rmg_seed_works(self):
        """Test that the created seed libraries work.

        Note: Since this test modifies the class level RMG instance,
        it can cause other tests to fail if run out of order."""
        # Load the seed libraries into the database
        self.rmg.database.load(
            path=self.databaseDirectory,
            thermo_libraries=[],
            reaction_libraries=['testSeed', 'testSeed_edge'],
            seed_mechanisms=['testSeed', 'testSeed_edge'],
            kinetics_families='default',
            kinetics_depositories=[],
            depository=False
        )

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_reaction_library_to_edge('testSeed')  # try adding seed as library
        self.assertTrue(len(self.rmg.reaction_model.edge.species) > 0)
        self.assertTrue(len(self.rmg.reaction_model.edge.reactions) > 0)

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_seed_mechanism_to_core('testSeed')  # try adding seed as seed mech
        self.assertTrue(len(self.rmg.reaction_model.core.species) > 0)
        self.assertTrue(len(self.rmg.reaction_model.core.reactions) > 0)

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_reaction_library_to_edge('testSeed_edge')  # try adding seed as library
        self.assertTrue(len(self.rmg.reaction_model.edge.species) > 0)
        self.assertTrue(len(self.rmg.reaction_model.edge.reactions) > 0)

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_seed_mechanism_to_core('testSeed_edge')  # try adding seed as seed mech
        self.assertTrue(len(self.rmg.reaction_model.core.species) > 0)
        self.assertTrue(len(self.rmg.reaction_model.core.reactions) > 0)

    def test_rmg_memory(self):
        """
        test that RMG Memory objects function properly
        """
        for rxnsys in self.rmg.reaction_systems:
            Rmem = RMG_Memory(rxnsys, None)
            Rmem.generate_cond()
            Rmem.get_cond()
            Rmem.add_t_conv_N(1.0, .2, 2)
            Rmem.generate_cond()
            Rmem.get_cond()

    def test_make_cantera_input_file(self):
        """
        This tests to ensure that a usable Cantera input file is created.
        """
        import cantera as ct

        outName = os.path.join(self.rmg.output_directory, 'cantera')
        files = os.listdir(outName)
        for f in files:
            if '.yaml' in f:
                try:
                    ct.Solution(os.path.join(outName, f))
                except:
                    self.fail('The output Cantera file is not loadable in Cantera.')


@attr('functional')
class TestRestartWithFilters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, 'rmg', 'test_data', 'restartTest')
        cls.outputDir = os.path.join(cls.testDir, 'output_w_filters')
        cls.databaseDirectory = settings['database.directory']

        os.mkdir(cls.outputDir)
        initialize_log(logging.INFO, os.path.join(cls.outputDir, 'RMG.log'))

        cls.rmg = RMG(input_file=os.path.join(cls.testDir, 'restart_w_filters.py'),
                      output_directory=os.path.join(cls.outputDir))

    def test_restart_with_filters(self):
        """
        Test that the RMG restart job with filters included completed without problems
        """
        self.rmg.execute()
        with open(os.path.join(self.outputDir, 'RMG.log'), 'r') as f:
            self.assertIn('MODEL GENERATION COMPLETED', f.read())

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


@attr('functional')
class TestRestartNoFilters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, 'rmg', 'test_data', 'restartTest')
        cls.outputDir = os.path.join(cls.testDir, 'output_no_filters')
        cls.databaseDirectory = settings['database.directory']

        os.mkdir(cls.outputDir)
        initialize_log(logging.INFO, os.path.join(cls.outputDir, 'RMG.log'))

        cls.rmg = RMG(input_file=os.path.join(cls.testDir, 'restart_no_filters.py'),
                      output_directory=os.path.join(cls.outputDir))

    def test_restart_no_filters(self):
        """
        Test that the RMG restart job with no filters included completed without problems
        """
        self.rmg.execute()
        with open(os.path.join(self.outputDir, 'RMG.log'), 'r') as f:
            self.assertIn('MODEL GENERATION COMPLETED', f.read())

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


@attr('functional')
class TestMainFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, 'rmg', 'test_data', 'mainTest')
        cls.outputDir = os.path.join(cls.testDir, 'output')
        cls.databaseDirectory = settings['database.directory']

        os.mkdir(os.path.join(cls.testDir, cls.outputDir))

        cls.max_iter = 10

        cls.rmg = RMG(input_file=os.path.join(cls.testDir, 'superminimal_input.py'),
                      output_directory=cls.outputDir)

        cls.rmg.execute(max_iterations=cls.max_iter)

    def test_save_seed_modulus(self):
        """
        Test that saveSeedModulus argument from superminimal_input.py saved the correct number of seeds
        """
        path = os.path.join(self.outputDir, 'previous_seeds')
        num_dir_actual = sum(os.path.isdir(os.path.join(path, i)) for i in os.listdir(path))
        num_dir_expected = self.max_iter//2 + 1   # +1 is for saving iteration 0
        self.assertEqual(num_dir_actual, num_dir_expected)

    def test_max_iter(self):
        """
        Test the command line argument of -i
        """
        df = pd.read_excel(os.path.join(self.outputDir, 'statistics.xls'))
        num_rows = df.shape[0]

        num_iter_actual = num_rows
        num_iter_expected = self.max_iter + 1  # +1 is for saving iteration 0
        self.assertEqual(num_iter_actual, num_iter_expected)

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


class TestProfiling(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        # Making the profile graph requires a display. See if one is available first
        cls.display_found = False

        try:
            cls.display_found = bool(os.environ['DISPLAY'])
        except KeyError:  # This means that no display was found
            pass
        cls.test_dir = os.path.join(originalPath, 'rmg', 'test_data', 'mainTest')

    @patch('rmgpy.rmg.main.logging')
    def test_make_profile_graph(self, mock_logging):
        """
        Test that `make_profile_graph` function behaves properly given the current display state
        """
        profile_file = os.path.join(self.test_dir, 'RMG.profile')
        make_profile_graph(profile_file)
        if self.display_found:  # Check that the profile graph was made
            self.assertTrue(os.path.exists(os.path.join(self.test_dir, 'RMG.profile.dot.pdf')))
        else:  # We can't test making a profile graph on this system, but at least test that this was recognized
            mock_logging.warning.assert_called_with(
                'Could not find a display, which is required in order to generate '
                'the profile graph. This '
                'is likely due to this job being run on a remote server without performing X11 forwarding '
                'or running the job through a job manager like SLURM.\n\n The graph can be generated later '
                'by running with the postprocessing flag `rmg.py -P input.py` from any directory/computer '
                'where both the input file and RMG.profile file are located and a display is available.\n\n'
                'Note that if the postprocessing flag is specified, this will force the graph generation '
                'regardless of if a display was found, which could cause this program to crash or freeze.'
            )

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""

        if cls.display_found:  # Remove output PDF
            os.remove(os.path.join(cls.test_dir, 'RMG.profile.dot.pdf'))
            os.remove(os.path.join(cls.test_dir, 'RMG.profile.dot'))
            os.remove(os.path.join(cls.test_dir, 'RMG.profile.dot.ps2'))


class TestCanteraOutput(unittest.TestCase):

    def setUp(self):
        self.chemkin_files = {"""ELEMENTS
	H
	D /2.014/
	T /3.016/
	C
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
    300.000  1000.000  5000.000

ethane(1)               H 6  C 2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

CH3(4)                  H 3  C 1            G100.000   5000.000  1337.62       1
 3.54144859E+00 4.76788187E-03-1.82149144E-06 3.28878182E-10-2.22546856E-14    2
 1.62239622E+04 1.66040083E+00 3.91546822E+00 1.84153688E-03 3.48743616E-06    3
-3.32749553E-09 8.49963443E-13 1.62856393E+04 3.51739246E-01                   4

END



REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": True,
"""ELEMENTS
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
    300.000  1000.000  5000.000

ethane(1)               H 6  C  2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

CH3(4)                  H 3  C 1            G100.000   5000.000  1337.62       1
 3.54144859E+00 4.76788187E-03-1.82149144E-06 3.28878182E-10-2.22546856E-14    2
 1.62239622E+04 1.66040083E+00 3.91546822E+00 1.84153688E-03 3.48743616E-06    3
-3.32749553E-09 8.49963443E-13 1.62856393E+04 3.51739246E-01                   4

END



REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": False,
"""ELEMENTS
	H
	D /2.014/
	T /3.016/
	C
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
    300.000  1000.000  5000.000

ethane(1)               H 6  C 2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

END

REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": False,
        }
        self.rmg = RMG()
        self.dir_name = 'temp_dir_for_testing'
        self.rmg.output_directory = os.path.join(originalPath, self.dir_name)

        self.tran_dat = '''
! Species         Shape    LJ-depth  LJ-diam   DiplMom   Polzblty  RotRelaxNum Data     
! Name            Index    epsilon/k_B sigma     mu        alpha     Zrot      Source   
ethane(1)           2     252.301     4.302     0.000     0.000     1.500    ! GRI-Mech
CH3(4)              2     144.001     3.800     0.000     0.000     0.000    ! GRI-Mech
        '''

    def tearDown(self):
        os.chdir(originalPath)
        # try to remove the tree. If testChemkinToCanteraConversion properly
        # ran, the files should already be removed.
        try:
            shutil.rmtree(self.dir_name)
        except OSError:
            pass
        # go back to the main RMG-Py directory
        os.chdir('..')

    def test_chemkin_to_cantera_conversion(self):
        """
        Tests that good and bad chemkin files raise proper exceptions
        """

        from cantera.ck2yaml import InputError

        for ck_input, works in self.chemkin_files.items():
            os.chdir(originalPath)
            os.mkdir(self.dir_name)
            os.chdir(self.dir_name)

            f = open('chem001.inp', 'w')
            f.write(ck_input)
            f.close()

            f = open('tran.dat', 'w')
            f.write(self.tran_dat)
            f.close()

            if works:
                self.rmg.generate_cantera_files(os.path.join(os.getcwd(), 'chem001.inp'))
            else:
                with self.assertRaises(InputError):
                    self.rmg.generate_cantera_files(os.path.join(os.getcwd(), 'chem001.inp'))

            # clean up
            os.chdir(originalPath)
            shutil.rmtree(self.dir_name)
