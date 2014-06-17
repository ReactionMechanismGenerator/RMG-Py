#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from unittest import TestCase
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.species import Species
from rmgpy.data.solvation import *
from rmgpy.molecule.molecule import Molecule

###################################################

class TestSoluteDatabase(TestCase):
    
    def runTest(self):
        pass

    def testSoluteGeneration(self):
        "Test we can estimate Abraham solute paramaters correctly"
        
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))
        
        self.testCases = [
        # from RMG-Java test runs by RWest (mostly in agreement with Jalan et. al. supplementary data)
        # ['octane', 'CCCCCCCC', 0.127,   0.085, 0.04,    3.766,  0.0030  , 1.2358],  
#         ['water', 'O', 0.524,   0.378,  0.309,  0.802,  0.348,  0.1673],
#         ['chrysene', 'C1=CC=CC3=C1C=CC4=C2C=CC=CC2=CC=C34', 1.603, 0.317, 2.864, 10.222, 0.003, None],
#         ['anthracene', 'C2=CC=CC3=CC1=CC=CC=C1C=C23', 1.261, 0.257, 2.128, 7.796, 0.003, None],
        ['1,2-ethanediol', 'C(CO)O', 0.823, 0.685, 0.327, 2.572, 0.693, None],
        # ['1-decanol', 'C(CCCCCCCCC)O', 0.449, 0.385, 0.205, 5.614, 0.348, None],
#         ['acetylaldehyde', 'CC=O', 0.622,   0.423,  0.171,  1.415,  0.0030, 0.4061],
#         ['acetic acid', 'C(C)(=O)O', 0.508, 0.411, 0.152, 1.873, 0.591, None],
#         ['ethyloctadecanoate', 'C(C)OC(CCCCCCCCCCCCCCCCC)=O', 0.558, 0.424, 0.08, 10.344, 0.003, None],
#         ['naphthalene', 'C1=CC=CC2=CC=CC=C12', 0.919, 0.197, 1.392, 5.37, 0.003, None],
#         ['phenol', 'C1(=CC=CC=C1)O', 0.875, 0.433, 0.829, 3.771, 0.546, None],
#         ['m-cresol', 'C1=C(C=CC=C1O)C', 0.851, 0.429, 0.837, 4.247, 0.546, None],
#         ['m-hydroxybenzaldehyde', 'OC1=CC=CC(=C1)C=O', 1.346, 0.767, 0.968, 4.89, 0.546, None],
#         ['ethylbenzene', 'C(C)C1=CC=CC=C1', 0.553, 0.133, 0.664, 3.919, 0.003, None],
#         ['toluene', 'C1(=CC=CC=C1)C', 0.553, 0.133, 0.664, 3.42, 0.003, None],
#         ['1-butene', 'C=CCC', 0.167, 0.108, 0.167, 1.663, 0.003, None],
#         ['g-decalactone', 'C(C1OC(=O)CC1)CCCCC', 0.669, 0.428, 0.273, 5.482, 0.003, None],

        #['dimethyl ether']
        ]
        
        for name, smiles, S, B, E, L, A, V in self.testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            soluteData = self.database.getSoluteDataFromGroups(Species(molecule=[species.molecule[0]]))
            self.assertAlmostEqual(soluteData.S, S)
            self.assertAlmostEqual(soluteData.B, B)
            self.assertAlmostEqual(soluteData.E, E)
            self.assertAlmostEqual(soluteData.L, L)
            self.assertAlmostEqual(soluteData.A, A)

    @work_in_progress
    def testCorrectionGeneration(self):
        "Test we can estimate solvation thermochemistry."
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))
        self.testCases = [
        # solventName, soluteName, soluteSMILES, Hsolv, Gsolv, Ssolv
        ['water', 'methane', 'C', -12000, 2000*4.184, None],
        ['water', 'octane', 'CCCCCCCC', -36000, 2890*4.184, None],
        ['water', '1,2-ethanediol', 'C(CO)O', -77300, -9300*4.184, None],
        ['water', 'acetic acid', 'C(C)(=O)O', -56500, -6700*4.184, None],
        ['water', 'naphthalene', 'C1=CC=CC2=CC=CC=C12', -42800, -2390*4.184, None],
        ['water', 'm-hydroxybenzaldehyde', 'OC1=CC=CC(=C1)C=O', -70700, -9510*4.184, None],
        ['water', 'ethylbenzene', 'C(C)C1=CC=CC=C1', -39400, -800*4.184, None],
        ['water', 'toluene', 'C1(=CC=CC=C1)C', -32400, -890*4.184, None],
        ['water', 'ethane', 'CC', -17900, 1830*4.184, None],
        ['water', 'propane', 'CCC', -20400, 1960*4.184, None],
        ['water', 'ethene', 'C=C', -13700, 1270*4.184, None],
        ['water', 'propene', 'CC=C', -21600, 1270*4.184, None],
        ['water', 'dimethyl ether', 'COC', -34000, -1920*4.184, None],
        ['water', 'diethyl ether', 'C(C)OCC', -45300, -1760*4.184, None],
        ['water', 'tetrahydrofuran', 'C1CCOC1', -47300, -3470*4.184, None],
        ['water', '1,4-dioxane', 'C1COCCO1', -48400, -5050*4.184, None],
        ['water', 'methanol', 'CO', -52000, -5110*4.184, None],
        ['water', 'ethanol', 'C(C)O', -50600, -5010*4.184, None],
        ['water', '1,2 propanediol', 'C(CO)CO', -81100, None, None],
        ['1-octanol', 'methane', 'C', -3900, 510*4.184, None],
        ['1-octanol', 'octane', 'CCCCCCCC', -40080, -4180*4.184, None],
        ['1-octanol', 'ethylbenzene', 'C(C)C1=CC=CC=C1', -39910, -5080*4.184, None],
        ['1-octanol', 'toluene', 'C1(=CC=CC=C1)C', -35990, -4550*4.184, None],
        ['1-octanol', 'tetrahydrofuran', 'C1CCOC1', -28320, -3930*4.184, None],
        ['benzene', 'octane', 'CCCCCCCC', -35480, -5350*4.184, None],
        ['benzene', 'toluene', 'C1(=CC=CC=C1)C', -37660, -5320*4.184, None],
        ['benzene', '1,4-dioxane', 'C1COCCO1', -39030, -5210*4.184, None],
        ['benzene', 'methanol', 'CO', None, -2580*4.184, None],
        ['benzene', 'ethanol', 'C(CO)CO', -29620, -3420*4.184, None],
        
        ]
        
        for solventName, soluteName, smiles, H, G, S in self.testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            soluteData = self.database.getSoluteDataFromGroups(Species(molecule=[species.molecule[0]]))
            solventData = self.database.getSolventData(solventName)
            solvationCorrection = self.database.getSolvationCorrection(soluteData, solventData)
            #print("Enthalpy of solvation for {0} in {1} is {2} J/mol".format(soluteName, solventName, solvationCorrection.enthalpy))
            #print("Enthalpy: {0} {1} {2}".format(soluteName, H, solvationCorrection.enthalpy))
            #print("Gibbs: {0} {1} {2}".format(soluteName, G, solvationCorrection.gibbs))
            self.assertAlmostEqual(solvationCorrection.enthalpy / 10000., H / 10000., 0, msg="Solvation enthalpy discrepancy ({2:.0f}!={3:.0f}) for {0} in {1}".format(soluteName, solventName, solvationCorrection.enthalpy, H))  #0 decimal place, in 10kJ.

#####################################################

if __name__ == '__main__':
    myTest = TestSoluteDatabase()
    #myTest.testSoluteGeneration()
    myTest.testCorrectionGeneration()
