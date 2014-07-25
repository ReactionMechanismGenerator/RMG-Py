#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import unittest
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.species import Species
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.molecule.molecule import Molecule

################################################################################

class TestThermoDatabase(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """
    # Only load these once to save time
    database = ThermoDatabase()
    database.load(os.path.join(settings['database.directory'], 'thermo'))
#    oldDatabase = ThermoDatabase()
#    oldDatabase.loadOld(os.path.join(settings['database.directory'], '../output/RMG_database'))
    
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        
        self.database = self.__class__.database
#        self.oldDatabase = self.__class__.oldDatabase

        self.Tlist = [300, 400, 500, 600, 800, 1000, 1500]
        
        self.testCases = [
            # SMILES            symm  H298     S298     Cp300  Cp400  Cp500  Cp600  Cp800  Cp1000 Cp1500
            
            # 1,3-hexadiene decomposition products
            ['C=CC=CCC',        3,    13.45, 86.37, 29.49, 37.67, 44.54, 50.12, 58.66, 64.95, 74.71],
            ['[CH]=CC=CCC',     3,    72.55, 87.76, 29.30, 36.92, 43.18, 48.20, 55.84, 61.46, 70.18],
            ['C=[C]C=CCC',      3,    61.15, 87.08, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=C[C]=CCC',      3,    61.15, 87.08, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=CC=[C]CC',      3,    70.35, 88.18, 29.15, 36.46, 42.6, 47.6, 55.32, 61.04, 69.95],
            ['C=CC=C[CH]C',     6,    38.24, 84.41, 27.79, 35.46, 41.94, 47.43, 55.74, 61.92, 71.86],
            ['C=CC=CC[CH2]',    2,    62.45, 89.78, 28.72, 36.31, 42.63, 47.72, 55.50, 61.21, 70.05],
            ['[CH3]',           6,    34.81, 46.37, 9.14, 10.18, 10.81, 11.34, 12.57, 13.71, 15.2],
            ['C=CC=C[CH2]',     2,    46.11, 75.82, 22.54, 28.95, 34.24, 38.64, 45.14, 49.97, 57.85],
            ['[CH2]C',          6,    28.6, 59.87, 11.73, 14.47, 17.05, 19.34, 23.02, 25.91, 31.53],
            ['C=CC=[CH]',       1,    85.18, 69.37, 18.93, 23.55, 27.16, 29.92, 34.02, 37.03, 41.81],
            ['C=[CH]',          1,    71.62, 56.61, 10.01, 11.97, 13.66, 15.08, 17.32, 19.05, 21.85],
            ['[CH]=CCC',        3,    58.99, 75.0, 20.38, 25.34, 29.68, 33.36, 39.14, 43.48, 50.22],
            
            # Cyclic Structures
            ['C1CCCCC1',        12,   -29.45, 69.71, 27.20, 37.60, 46.60, 54.80, 67.50, 76.20, 88.50],
            ['C1CCC1',          8,     6.51, 63.35, 17.39, 23.91, 29.86, 34.76, 42.40, 47.98, 56.33],
            ['C1C=CC=C1',       2,    32.5, 65.5, 18.16, 24.71, 30.25, 34.7, 41.25, 45.83, 52.61],
        ]

    @work_in_progress
    def testNewThermoGeneration(self):
        """
        Test that the new ThermoDatabase generates appropriate thermo data.
        """
        
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            Cplist = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
            molecule=Molecule(SMILES=smiles)
            species = Species(molecule=molecule)
            species.generateResonanceIsomers()
            species.molecule[0]
            thermoData = self.database.getThermoDataFromGroups(species)
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.database.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.database.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            self.assertAlmostEqual(H298, thermoData.getEnthalpy(298) / 4184, places=1, msg="H298 error for {0}".format(smiles))
            self.assertAlmostEqual(S298, thermoData.getEntropy(298) / 4.184, places=1, msg="S298 error for {0}".format(smiles))
            for T, Cp in zip(self.Tlist, Cplist):
                self.assertAlmostEqual(Cp, thermoData.getHeatCapacity(T) / 4.184, places=1, msg="Cp{1} error for {0}".format(smiles,T))

    @work_in_progress
    def testSymmetryNumberGeneration(self):
        """
        Test we generate symmetry numbers correctly.
        
        This uses the new thermo database to generate the H298, used 
        to select the stablest resonance isomer.
        """
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            molecule=Molecule(SMILES=smiles)
            species = Species(molecule=molecule)
            species.generateResonanceIsomers()
            thermoData = self.database.getThermoDataFromGroups(Species(molecule=[species.molecule[0]]))
            # pick the molecule with lowest H298
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.database.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.database.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            self.assertEqual(molecule.calculateSymmetryNumber(), symm, msg="Symmetry number error for {0}".format(smiles))

#    @work_in_progress
#    def testOldThermoGeneration(self):
#        """
#        Test that the old ThermoDatabase generates relatively accurate thermo data.
#        """
#        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
#            Cplist = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
#            species = Species(molecule=[Molecule(SMILES=smiles)])
#            species.generateResonanceIsomers()
#            thermoData = self.oldDatabase.getThermoData(Species(molecule=[species.molecule[0]]))
#            molecule = species.molecule[0]
#            for mol in species.molecule[1:]:
#                thermoData0 = self.oldDatabase.getAllThermoData(Species(molecule=[mol]))[0][0]
#                for data in self.oldDatabase.getAllThermoData(Species(molecule=[mol]))[1:]:
#                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
#                        thermoData0 = data
#                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
#                    thermoData = thermoData0
#                    molecule = mol
#            
#            self.assertAlmostEqual(H298, thermoData.getEnthalpy(298) / 4184, places=1, msg="H298 error for {0}".format(smiles))
#            self.assertAlmostEqual(S298, thermoData.getEntropy(298) / 4.184, places=1, msg="S298 error for {0}".format(smiles))
#            for T, Cp in zip(self.Tlist, Cplist):
#                self.assertAlmostEqual(Cp, thermoData.getHeatCapacity(T) / 4.184, places=1, msg="Cp{1} error for {0}".format(smiles, T))

class TestThermoDatabaseAromatics(TestThermoDatabase):
    """
    Test only Aromatic species.
    
    A copy of the above class, but with different test compounds
    """
    def setUp(self):
        TestThermoDatabase.setUp(self)
        self.testCases = [
            # SMILES            symm  H298     S298     Cp300  Cp400  Cp500  Cp600  Cp800  Cp1000 Cp1500
            ['c1ccccc1', 12, 19.80, 64.24, 19.44, 26.64, 32.76, 37.80, 45.24, 50.46, 58.38],
            ['c1ccc2ccccc2c1', 4, 36.0, 79.49, 31.94, 42.88, 52.08, 59.62, 70.72, 78.68, 90.24],
        ]
    def __init__(self, *args, **kwargs):
        super(TestThermoDatabaseAromatics, self).__init__(*args, **kwargs)
        self._testMethodDoc = self._testMethodDoc.strip().split('\n')[0] + " for Aromatics.\n"

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

