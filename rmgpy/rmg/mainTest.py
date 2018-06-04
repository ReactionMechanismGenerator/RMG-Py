#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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
import shutil 
from nose.plugins.attrib import attr
from main import RMG
from main import RMG_Memory
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy import getPath
from rmgpy.rmg.model import CoreEdgeReactionModel
###################################################

originalPath = getPath()
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

        cls.rmg = RMG(inputFile=os.path.join(cls.testDir, 'input.py'),
                      outputDirectory=os.path.join(cls.testDir, cls.outputDir))

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

    def testRMGExecute(self):
        """Test that RMG.execute completed successfully."""
        self.assertIsInstance(self.rmg.database, RMGDatabase)
        self.assertTrue(self.rmg.done)

    def testRMGIncreasesReactions(self):
        """Test that RMG.execute increases reactions and species."""
        self.assertTrue(len(self.rmg.reactionModel.core.reactions) > 0)
        self.assertTrue(len(self.rmg.reactionModel.core.species) > 1)
        self.assertTrue(len(self.rmg.reactionModel.edge.reactions) > 0)
        self.assertTrue(len(self.rmg.reactionModel.edge.species) > 0)

    def testRMGSeedMechanismCreation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seedDir = os.path.join(self.testDir, self.outputDir, 'seed')
        self.assertTrue(os.path.exists)

        self.assertTrue(os.path.exists(os.path.join(seedDir, self.rmg.name)))  # kinetics library folder made

        self.assertTrue(os.path.exists(os.path.join(seedDir, self.rmg.name, 'dictionary.txt')))  # dictionary file made
        self.assertTrue(os.path.exists(os.path.join(seedDir, self.rmg.name, 'reactions.py')))  # reactions file made

    def testRMGSeedEdgeMechanismCreation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seedDir = os.path.join(self.testDir, self.outputDir, 'seed')
        self.assertTrue(os.path.exists)

        name = self.rmg.name + '_edge'

        self.assertTrue(os.path.exists(os.path.join(seedDir, name)))  # kinetics library folder made

        self.assertTrue(os.path.exists(os.path.join(seedDir, name, 'dictionary.txt')))  # dictionary file made
        self.assertTrue(os.path.exists(os.path.join(seedDir, name, 'reactions.py')))  # reactions file made

    def testRMGSeedLibraryCreation(self):
        """Test that seed mechanisms are created in the correct database locations."""
        self.assertTrue(os.path.exists(self.seedKinetics))

    def testRMGSeedEdgeLibraryCreation(self):
        """Test that edge seed mechanisms are created in the correct database locations."""
        self.assertTrue(os.path.exists(self.seedKinetics))

    def testRMGSeedWorks(self):
        """Test that the created seed libraries work.

        Note: Since this test modifies the class level RMG instance,
        it can cause other tests to fail if run out of order."""
        # Load the seed libraries into the database
        self.rmg.database.load(
            path=self.databaseDirectory,
            thermoLibraries=[],
            reactionLibraries=['testSeed', 'testSeed_edge'],
            seedMechanisms=['testSeed', 'testSeed_edge'],
            kineticsFamilies='default',
            kineticsDepositories=[],
            depository=False
        )
        
        self.rmg.reactionModel = CoreEdgeReactionModel()
        self.rmg.reactionModel.addReactionLibraryToEdge('testSeed')  # try adding seed as library
        self.assertTrue(len(self.rmg.reactionModel.edge.species) > 0)
        self.assertTrue(len(self.rmg.reactionModel.edge.reactions) > 0)
        
        self.rmg.reactionModel = CoreEdgeReactionModel()
        self.rmg.reactionModel.addSeedMechanismToCore('testSeed')  # try adding seed as seed mech
        self.assertTrue(len(self.rmg.reactionModel.core.species) > 0)
        self.assertTrue(len(self.rmg.reactionModel.core.reactions) > 0)

        self.rmg.reactionModel = CoreEdgeReactionModel()
        self.rmg.reactionModel.addReactionLibraryToEdge('testSeed_edge')  # try adding seed as library
        self.assertTrue(len(self.rmg.reactionModel.edge.species) > 0)
        self.assertTrue(len(self.rmg.reactionModel.edge.reactions) > 0)

        self.rmg.reactionModel = CoreEdgeReactionModel()
        self.rmg.reactionModel.addSeedMechanismToCore('testSeed_edge')  # try adding seed as seed mech
        self.assertTrue(len(self.rmg.reactionModel.core.species) > 0)
        self.assertTrue(len(self.rmg.reactionModel.core.reactions) > 0)
    
    def testRMGMemory(self):
        """
        test that RMG Memory objects function properly
        """
        for rxnsys in self.rmg.reactionSystems:
            Rmem = RMG_Memory(rxnsys,None)
            Rmem.generate_cond()
            Rmem.get_cond()
            Rmem.add_t_conv_N(1.0,.2,2)
            Rmem.generate_cond()
            Rmem.get_cond()
        
    def testMakeCanteraInputFile(self):
        """
        This tests to ensure that a usable Cantera input file is created.
        """
        import cantera as ct
        
        outName = os.path.join(self.rmg.outputDirectory, 'cantera')
        files = os.listdir(outName)
        for f in files:
            if '.cti' in f:
                try:
                    ct.Solution(os.path.join(outName, f))
                except:
                    self.fail('The output Cantera file is not loadable in Cantera.')


class TestCanteraOutput(unittest.TestCase):
    
    def setUp(self):
        self.chemkin_files={"""ELEMENTS
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
        self.rmg.outputDirectory = os.path.join(originalPath, self.dir_name)

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

    def testChemkinToCanteraConversion(self):
        """
        Tests that good and bad chemkin files raise proper exceptions
        """
        
        from cantera.ck2cti import InputParseError
        
        for ck_input, works in self.chemkin_files.items():
            os.chdir(originalPath)
            os.mkdir(self.dir_name)
            os.chdir(self.dir_name)
            
            f = open('chem001.inp','w')
            f.write(ck_input)
            f.close()
            
            f = open('tran.dat','w')
            f.write(self.tran_dat)
            f.close()
            
            if works:
                self.rmg.generateCanteraFiles(os.path.join(os.getcwd(),'chem001.inp'))
            else:
                with self.assertRaises(InputParseError):
                    self.rmg.generateCanteraFiles(os.path.join(os.getcwd(),'chem001.inp'))
            
            # clean up
            os.chdir(originalPath)
            shutil.rmtree(self.dir_name)
        
