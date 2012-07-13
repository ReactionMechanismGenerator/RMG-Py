#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import unittest

from rmgpy import settings
from rmgpy.species import Species
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.molecule.molecule import Molecule

################################################################################

class TestThermoDatabase(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        
        self.database = ThermoDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'thermo'))
        
        self.oldDatabase = ThermoDatabase()
        self.oldDatabase.loadOld(os.path.join(settings['database.directory'], '../output/RMG_database'))
        
        self.Tlist = [300, 400, 500, 600, 800, 1000, 1500]
        
        self.testCases = [
            # SMILES            symm  H298     S298     Cp300  Cp400  Cp500  Cp600  Cp800  Cp1000 Cp1500
            
            # 1,3-hexadiene decomposition products
            ['C=CC=CCC',        3,    13.5090, 86.5641, 29.49, 37.67, 44.54, 50.12, 58.66, 64.95, 74.71],
            ['[CH]=CC=CCC',     3,    72.6056, 87.9528, 29.30, 36.92, 43.18, 48.20, 55.84, 61.46, 70.18],
            ['C=[C]C=CCC',      3,    61.2064, 87.2754, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=C[C]=CCC',      3,    61.2064, 87.2754, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=CC=[C]CC',      3,    70.4053, 88.3718, 29.15, 36.46, 42.60, 47.60, 55.32, 61.04, 69.95],
            ['C=CC=C[CH]C',     6,    38.2926, 84.5953, 27.79, 35.46, 41.94, 47.43, 55.74, 61.92, 71.86],
            ['C=CC=CC[CH2]',    2,    62.5044, 89.9747, 28.72, 36.31, 42.63, 47.72, 55.50, 61.21, 70.05],
            ['[CH3]',           6,    35.1084, 46.3644,  9.20,  9.98, 10.75, 11.50, 12.86, 14.08, 16.29],
            ['C=CC=C[CH2]',     2,    46.1521, 75.9733, 22.54, 28.95, 34.24, 38.64, 45.14, 49.97, 57.85],
            ['[CH2]C',          6,    28.3580, 59.0565, 12.11, 14.59, 17.08, 19.35, 22.93, 25.78, 30.30],
            ['C=CC=[CH]',       1,    85.2149, 69.4966, 18.93, 23.55, 27.16, 29.92, 34.02, 37.03, 41.81],
            ['C=[CH]',          1,    71.6377, 55.8964, 10.24, 12.03, 13.71, 15.17, 17.35, 19.07, 21.82],
            ['[CH]=CCC',        3,    59.0278, 75.1332, 20.38, 25.34, 29.68, 33.36, 39.14, 43.48, 50.22],
            
            # Cyclic structures
            ['c1ccccc1',        1,    19.8389, 69.3100, 19.44, 26.64, 32.76, 37.80, 45.24, 50.46, 58.38],
            ['C1CCCCC1',        1,   -29.4456, 74.8296, 27.20, 37.60, 46.60, 54.80, 67.50, 76.20, 88.50],
            ['c1ccc2ccccc2c1',  1,    36.0639, 82.4536, 31.94, 42.88, 52.08, 59.62, 70.72, 78.68, 90.24],
            ['C1CCC1',          1,     6.5148, 67.5963, 17.39, 23.91, 29.86, 34.76, 42.40, 47.98, 56.33],
            ['C1C=CC=C1',       1,    32.5363, 67.0035, 18.16, 24.71, 30.25, 34.70, 41.25, 45.83, 52.61],
        ]

    def testNewThermoGeneration(self):
        """
        Test that the new ThermoDatabase generates appropriate thermo data.
        """
        
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            Cplist = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
            species = Species(molecule=[Molecule(SMILES=smiles)])
            species.generateResonanceIsomers()
            thermoData = self.database.getThermoData(Species(molecule=[species.molecule[0]]))
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.database.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.database.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            
            self.assertEqual(molecule.calculateSymmetryNumber(), symm)
            self.assertTrue(1 - thermoData.getEnthalpy(298) / 4184 / H298 < 0.001)
            self.assertTrue(1 - thermoData.getEntropy(298) / 4.184 / S298 < 0.001)
            for T, Cp in zip(self.Tlist, Cplist):
                self.assertTrue(1 - thermoData.getHeatCapacity(T) / 4.184 / Cp < 0.001)

    def testOldThermoGeneration(self):
        """
        Test that the old ThermoDatabase generates relatively accurate thermo data.
        """
        
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            Cplist = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
            species = Species(molecule=[Molecule(SMILES=smiles)])
            species.generateResonanceIsomers()
            thermoData = self.oldDatabase.getThermoData(Species(molecule=[species.molecule[0]]))
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.oldDatabase.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.oldDatabase.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            
            self.assertEqual(molecule.calculateSymmetryNumber(), symm)
            self.assertTrue(1 - thermoData.getEnthalpy(298) / 4184 / H298 < 0.01)
            self.assertTrue(1 - thermoData.getEntropy(298) / 4.184 / S298 < 0.01)
            for T, Cp in zip(self.Tlist, Cplist):
                self.assertTrue(1 - thermoData.getHeatCapacity(T) / 4.184 / Cp < 0.1)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
