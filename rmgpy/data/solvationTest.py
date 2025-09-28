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

import os
from unittest import TestCase, TestLoader, TextTestRunner
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import Species
from rmgpy.data.solvation import DatabaseError, SoluteData, SolvationDatabase, SolventLibrary
from rmgpy.rmg.main import RMG

###################################################

class TestSoluteDatabase(TestCase):
    
    def setUp(self):
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))

    def tearDown(self):
        """
        Reset the database & liquid parameters for solution
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
        
    def runTest(self):
        pass
    
    def testSoluteLibrary(self):
        "Test we can obtain solute parameters from a library"
        species = Species(molecule=[Molecule(SMILES='COC=O')]) #methyl formate - we know this is in the solute library
        
        libraryData = self.database.getSoluteDataFromLibrary(species, self.database.libraries['solute'])
        self.assertEqual(len(libraryData), 3)
        
        soluteData = self.database.getSoluteData(species)
        self.assertTrue(isinstance(soluteData, SoluteData))
        
        S = soluteData.S
        self.assertEqual(S, 0.68)
        self.assertTrue(soluteData.V is not None)
     
    def testMcGowan(self):
        "Test we can calculate and set the McGowan volume for species containing H,C,O,N or S"
        self.testCases = [
                          ['CCCCCCCC', 1.2358], #n-octane, in library
                          ['C(CO)O', 0.5078], #ethylene glycol
                          ['CC#N', 0.4042], #acetonitrile
                          ['CCS', 0.5539] #ethanethiol
                           ]
        
        for smiles, volume in self.testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            soluteData = self.database.getSoluteData(species)
            soluteData.setMcGowanVolume(species) # even if it was found in library, recalculate
            self.assertTrue(soluteData.V is not None) # so if it wasn't found in library, we should have calculated it
            self.assertAlmostEqual(soluteData.V, volume) # the volume is what we expect given the atoms and bonds 
            
    
    def testDiffusivity(self):
        "Test that for a given solvent viscosity and temperature we can calculate a solute's diffusivity"
        species = Species(molecule=[Molecule(SMILES='O')])  # water
        soluteData = self.database.getSoluteData(species)
        T = 298.
        solventViscosity = 0.00089  # water is about 8.9e-4 Pa.s
        D = soluteData.getStokesDiffusivity(T, solventViscosity)  # m2/s
        self.assertAlmostEqual((D * 1e9), 1.3, 1)
        # self-diffusivity of water is about 2e-9 m2/s
        
    def testSolventLibrary(self):
        "Test we can obtain solvent parameters from a library"
        solventData = self.database.getSolventData('water')
        self.assertTrue(solventData is not None)
        self.assertEqual(solventData.s_h, 2.836)
        self.assertRaises(DatabaseError, self.database.getSolventData, 'orange_juice')
        
    def testViscosity(self):
        "Test we can calculate the solvent viscosity given a temperature and its A-E correlation parameters"
        solventData = self.database.getSolventData('water')  
        self.assertAlmostEqual(solventData.getSolventViscosity(298), 0.0009155)

    def testSoluteGeneration(self):
        "Test we can estimate Abraham solute parameters correctly using group contributions"
        
        self.testCases = [
        ['1,2-ethanediol', 'C(CO)O', 0.823, 0.685, 0.327, 2.572, 0.693, None],
        ]
        
        for name, smiles, S, B, E, L, A, V in self.testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            soluteData = self.database.getSoluteDataFromGroups(Species(molecule=[species.molecule[0]]))
            self.assertAlmostEqual(soluteData.S, S, places=2)
            self.assertAlmostEqual(soluteData.B, B, places=2)
            self.assertAlmostEqual(soluteData.E, E, places=2)
            self.assertAlmostEqual(soluteData.L, L, places=2)
            self.assertAlmostEqual(soluteData.A, A, places=2)

    def testLonePairSoluteGeneration(self):
        "Test we can obtain solute parameters via group additivity for a molecule with lone pairs"
        molecule=Molecule().fromAdjacencyList(
"""
CH2_singlet
multiplicity 1
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""")
        species = Species(molecule=[molecule])
        soluteData = self.database.getSoluteDataFromGroups(species)
        self.assertTrue(soluteData is not None)
        
    def testSoluteDataGenerationAmmonia(self):
        "Test we can obtain solute parameters via group additivity for ammonia"
        molecule=Molecule().fromAdjacencyList(
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""")
        species = Species(molecule=[molecule])
        soluteData = self.database.getSoluteDataFromGroups(species)
        self.assertTrue(soluteData is not None)
        
    def testSoluteDataGenerationAmide(self):
        "Test that we can obtain solute parameters via group additivity for an amide"
        molecule=Molecule().fromAdjacencyList(
"""
1 N u0 p1 {2,S} {3,S} {4,S}
2 H u0 {1,S}
3 C u0 {1,S} {6,S} {7,S} {8,S}
4 C u0 {1,S} {5,D} {9,S}
5 O u0 p2 {4,D}
6 H u0 {3,S}
7 H u0 {3,S}
8 H u0 {3,S}
9 H u0 {4,S}
""")
        species = Species(molecule=[molecule])
        soluteData = self.database.getSoluteDataFromGroups(species)
        self.assertTrue(soluteData is not None)

    def testSoluteDataGenerationCO(self):
        "Test that we can obtain solute parameters via group additivity for CO."
        molecule=Molecule().fromAdjacencyList(
"""
1  C u0 p1 c-1 {2,T}
2  O u0 p1 c+1 {1,T}
""")
        species = Species(molecule=[molecule])
        soluteData = self.database.getSoluteDataFromGroups(species)
        self.assertTrue(soluteData is not None)
    
    def testRadicalandLonePairGeneration(self):
        """
        Test we can obtain solute parameters via group additivity for a molecule with both lone 
        pairs and a radical
        """
        molecule=Molecule().fromAdjacencyList(
"""
[C]OH
multiplicity 2
1 C u1 p1 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
""")
        species = Species(molecule=[molecule])
        soluteData = self.database.getSoluteDataFromGroups(species)
        self.assertTrue(soluteData is not None)

    def testCorrectionGeneration(self):
        "Test we can estimate solvation thermochemistry."
        self.testCases = [
        # solventName, soluteName, soluteSMILES, Hsolv, Gsolv
        ['water', 'acetic acid', 'C(C)(=O)O', -56500, -6700*4.184],
        ['water', 'naphthalene', 'C1=CC=CC2=CC=CC=C12', -42800, -2390*4.184],
        ['1-octanol', 'octane', 'CCCCCCCC', -40080, -4180*4.184],
        ['1-octanol', 'tetrahydrofuran', 'C1CCOC1', -28320, -3930*4.184],
        ['benzene', 'toluene', 'C1(=CC=CC=C1)C', -37660, -5320*4.184],
        ['benzene', '1,4-dioxane', 'C1COCCO1', -39030, -5210*4.184]
        ]
        
        for solventName, soluteName, smiles, H, G in self.testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            soluteData = self.database.getSoluteData(species)
            solventData = self.database.getSolventData(solventName)
            solvationCorrection = self.database.getSolvationCorrection(soluteData, solventData)
            self.assertAlmostEqual(solvationCorrection.enthalpy / 10000., H / 10000., 0, msg="Solvation enthalpy discrepancy ({2:.0f}!={3:.0f}) for {0} in {1}".format(soluteName, solventName, solvationCorrection.enthalpy, H))  #0 decimal place, in 10kJ.
            self.assertAlmostEqual(solvationCorrection.gibbs / 10000., G / 10000., 0, msg="Solvation Gibbs free energy discrepancy ({2:.0f}!={3:.0f}) for {0} in {1}".format(soluteName, solventName, solvationCorrection.gibbs, G))

    def testInitialSpecies(self):
        " Test we can check whether the solvent is listed as one of the initial species in various scenarios "

        # Case 1. when SMILES for solvent is available, the molecular structures of the initial species and the solvent
        # are compared to check whether the solvent is in the initial species list

        # Case 1-1: the solvent water is not in the initialSpecies list, so it raises Exception
        rmg=RMG()
        rmg.initialSpecies = []
        solute = Species(label='n-octane', molecule=[Molecule().fromSMILES('C(CCCCC)CC')])
        rmg.initialSpecies.append(solute)
        rmg.solvent = 'water'
        solventStructure = Species().fromSMILES('O')
        self.assertRaises(Exception, self.database.checkSolventinInitialSpecies, rmg, solventStructure)

        # Case 1-2: the solvent is now octane and it is listed as the initialSpecies. Although the string
        # names of the solute and the solvent are different, because the solvent SMILES is provided,
        # it can identify the 'n-octane' as the solvent
        rmg.solvent = 'octane'
        solventStructure = Species().fromSMILES('CCCCCCCC')
        self.database.checkSolventinInitialSpecies(rmg, solventStructure)
        self.assertTrue(rmg.initialSpecies[0].isSolvent)

        # Case 2: the solvent SMILES is not provided. In this case, it can identify the species as the
        # solvent by looking at the string name.

        # Case 2-1: Since 'n-octane and 'octane' are not equal, it raises Exception
        solventStructure = None
        self.assertRaises(Exception, self.database.checkSolventinInitialSpecies, rmg, solventStructure)

        # Case 2-2: The label 'n-ocatne' is corrected to 'octane', so it is identified as the solvent
        rmg.initialSpecies[0].label = 'octane'
        self.database.checkSolventinInitialSpecies(rmg, solventStructure)
        self.assertTrue(rmg.initialSpecies[0].isSolvent)

    def testSolventMolecule(self):
        """Test that we can assign a proper solvent molecular structure when different formats are given"""

        # solventlibrary.entries['solvent_label'].item should be the instance of Species with the solvent's molecular structure
        # if the solvent database contains the solvent SMILES or adjacency list. If not, then item is None

        # Case 1: When the solventDatabase does not contain the solvent SMILES, the item attribute is None
        solventlibrary = SolventLibrary()
        solventlibrary.loadEntry(index=1, label='water', solvent=None)
        self.assertTrue(solventlibrary.entries['water'].item is None)

        # Case 2: When the solventDatabase contains the correct solvent SMILES, the item attribute is the instance of
        # Species with the correct solvent molecular structure
        solventlibrary.loadEntry(index=2, label='octane', solvent=None, molecule='CCCCCCCC')
        solvent_species = Species().fromSMILES('C(CCCCC)CC')
        self.assertTrue(solvent_species.isIsomorphic(solventlibrary.entries['octane'].item[0]))

        # Case 3: When the solventDatabase contains the correct solvent adjacency list, the item attribute is the instance of
        # the species with the correct solvent molecular structure.
        # This will display the SMILES Parse Error message from the external function, but ignore it.
        solventlibrary.loadEntry(index=3, label='ethanol', solvent=None, molecule=
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3 O u0 p2 c0 {2,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
""")
        solvent_species = Species().fromSMILES('CCO')
        self.assertTrue(solvent_species.isIsomorphic(solventlibrary.entries['ethanol'].item[0]))

        # Case 4: when the solventDatabase contains incorrect values for the molecule attribute, it raises Exception
        # This will display the SMILES Parse Error message from the external function, but ignore it.
        self.assertRaises(Exception, solventlibrary.loadEntry, index=4, label='benzene', solvent=None, molecule='ring')

        # Case 5: when the solventDatabase contains data for co-solvents.
        solventlibrary.loadEntry(index=5, label='methanol_50_water_50', solvent=None, molecule=['CO', 'O'])
        solvent_species_list = [Species().fromSMILES('CO'), Species().fromSMILES('O')]
        self.assertEqual(len(solventlibrary.entries['methanol_50_water_50'].item), 2)
        for spc1 in solventlibrary.entries['methanol_50_water_50'].item:
            self.assertTrue(any([spc1.isIsomorphic(spc2) for spc2 in solvent_species_list]))

#####################################################


if __name__ == '__main__':
    suite = TestLoader().loadTestsFromTestCase(TestSoluteDatabase)
    TextTestRunner(verbosity=2).run(suite)
