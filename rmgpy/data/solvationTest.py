#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from unittest import TestCase, TestLoader, TextTestRunner
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.solvation import DatabaseError, SoluteData, SolvationDatabase

###################################################

class TestSoluteDatabase(TestCase):
    
    def setUp(self):
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))
    
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
1  N u0 p1 {2,S} {3,S} {4,S}
2   H   u0 {1,S}
3   C u0 {1,S} {6,S} {7,S} {8,S}
4   C  u0 {1,S} {5,D} {9,S}
5   O  u0 p2 {4,D}
6   H   u0 {3,S}
7   H   u0 {3,S}
8   H   u0 {3,S}
9   H   u0 {4,S}
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
#####################################################

if __name__ == '__main__':
    suite = TestLoader().loadTestsFromTestCase(TestSoluteDatabase)
    TextTestRunner(verbosity=2).run(suite)
