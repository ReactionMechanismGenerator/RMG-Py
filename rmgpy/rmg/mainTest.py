#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os
import unittest
import shutil 

from rmgpy.chemkin import loadChemkinFile

from main import RMG
from rmgpy.data.rmg import RMGDatabase
from rmgpy import getPath
###################################################

originalPath = getPath()

class TestMain(unittest.TestCase):

    def setUp(self):
        self.dir_name = 'temp_dir_for_testing'
        os.chdir(originalPath)
        os.mkdir(self.dir_name)
        os.chdir(self.dir_name)
        inputFile = """
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['R_Recombination'],
    kineticsEstimator = 'rate rules',
)
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
    terminationConversion={
        'ethane': 0.000000000001,
    },
    terminationTime=(1e6,'s'),
)
simulator(
atol=1e-16,
rtol=1e-8
)
model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.2,
    toleranceInterruptSimulation=0.2,
)
options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
)
        """

        f = open('input.py','w')
        f.write(inputFile)
        f.close()

        self.rmg = RMG(inputFile=os.path.join(os.getcwd(), 'input.py'), outputDirectory=os.getcwd())

    def tearDown(self):
        os.chdir(originalPath)
        shutil.rmtree(self.dir_name)
        # go back to the main RMG-Py directory
        os.chdir('..')
        # remove modular level database
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    def testRMGExecute(self):
        """
        This example is to test if RMG.execute increases the core reactions
        """
        
        self.rmg.execute()
        self.assertIsInstance(self.rmg.database, RMGDatabase)
        self.assertTrue(self.rmg.done)
        self.assertTrue(len(self.rmg.reactionModel.core.reactions) > 0)
        self.assertTrue(len(self.rmg.reactionModel.core.species) > 1)
        self.assertTrue(len(self.rmg.reactionModel.edge.reactions) > 0)
        self.assertTrue(len(self.rmg.reactionModel.edge.species) > 0)
        
    def testMakeCanteraInputFile(self):
        """
        This tests to ensure that a usable cantera input file is created
        """
        
        self.rmg.execute()
        
        import cantera as ct
        
        outName = os.path.join(self.rmg.outputDirectory, 'cantera')
        files = os.listdir(outName)
        for f in files:
            if '.cti' in f:
                try:
                    ct.Solution(os.path.join(outName, f))
                except:
                    self.assertTrue(False, 'The output cantera file is not loadable in cantera')
                    
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
        
