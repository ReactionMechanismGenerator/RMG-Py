#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from unittest import TestCase, TestLoader, TextTestRunner
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.species import Species
from rmgpy.data.solvation import *
from rmgpy.molecule.molecule import Molecule
from rmgpy.data.base import DatabaseError

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
        species = Species(molecule=[Molecule(SMILES='COC=O')])
        soluteData = self.database.getSoluteData(species)
        T = 298
        solventViscosity = 0.001
        D = soluteData.getStokesDiffusivity(T, solventViscosity)
        self.assertAlmostEqual((D*1E12), 0.00000979)
        
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
        # from RMG-Java test runs by RWest (mostly in agreement with Jalan et. al. supplementary data)
        ['1,2-ethanediol', 'C(CO)O', 0.823, 0.685, 0.327, 2.572, 0.693, None],
        # a nitrogen case
        #['acetonitrile', 'CC#N', 0.9, 0.33, 0.237, 1.739, 0.04, None],
        # a sulfur case
        #['ethanethiol', 'CCS', 0.35, 0.24, 0.392, 2.173, 0, None]
        ]
        
        for name, smiles, S, B, E, L, A, V in self.testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            soluteData = self.database.getSoluteDataFromGroups(Species(molecule=[species.molecule[0]]))
            self.assertAlmostEqual(soluteData.S, S, places=2)
            self.assertAlmostEqual(soluteData.B, B, places=2)
            self.assertAlmostEqual(soluteData.E, E, places=2)
            self.assertAlmostEqual(soluteData.L, L, places=2)
            self.assertAlmostEqual(soluteData.A, A, places=2)

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
