#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

import rmgpy.data.thermo

from rmgpy.data.base import LogicNode
from rmgpy.chem.molecule import Molecule
from rmgpy.chem.pattern import MoleculePattern
from rmgpy.chem.thermo import ThermoData
from rmgpy.chem.species import Species

from rmgpy import settings

################################################################################

class ThermoDatabaseCheck(unittest.TestCase):

    def testOldThermoDatabase(self):
        """
        Check the database load functions.
        """
        
        thermoDatabase = loadThermoDatabase('database/output/RMG_Database/thermo_groups', group=True, old=True)

        for database in [thermoDatabase.groupDatabase,
            thermoDatabase.int15Database,
            thermoDatabase.gaucheDatabase,
            thermoDatabase.otherDatabase,
            thermoDatabase.radicalDatabase,
            thermoDatabase.ringDatabase]:

            # All nodes in library should be in tree and dictionary
            # All nodes in tree should be in dictionary
            for node in database.library:
                self.assertTrue(node in database.tree.parent)
                self.assertTrue(node in database.tree.children)
                self.assertTrue(node in database.dictionary)
            for node in database.tree.parent:
                self.assertTrue(node in database.tree.children)
                self.assertTrue(node in database.dictionary)
            for node in database.tree.children:
                self.assertTrue(node in database.tree.parent)
                self.assertTrue(node in database.dictionary)

            # All parents in tree should be in tree
            for node in database.tree.parent:
                parentNode = database.tree.parent[node]
                if parentNode is not None:
                    self.assertTrue(parentNode in database.tree.parent)
                    self.assertTrue(parentNode in database.tree.children)

            # All children in tree should be in tree
            for node in database.tree.children:
                for childNode in database.tree.children[node]:
                    if childNode is not None:
                        self.assertTrue(childNode in database.tree.parent)
                        self.assertTrue(childNode in database.tree.children)

            # All values in dictionary should be chemical structures
            for node in database.dictionary:
                self.assertTrue(isinstance(database.dictionary[node], MoleculePattern) or isinstance(database.dictionary[node], LogicNode))

            # All values in library should be ThermoEntry objects
            for node in database.library:
                self.assertTrue(isinstance(database.library[node], ThermoEntry), '"%s" is of unexpected type "%s".' % (node, database.library[node].__class__))
                self.assertTrue(database.library[node].model is not None or database.library[node].node != '')

        molecule = Molecule(SMILES='CC')
        print
        print thermoDatabase.generateThermoData(molecule)

    def testThermoDatabase(self):
        """
        Check the database load functions.
        """
        thermoDatabase = rmgpy.data.thermo.ThermoDatabase()
        thermoDatabase.load(os.path.join(settings['database.path'],'RMG_Database'))

        for database in [thermoDatabase.groupDatabase,
                thermoDatabase.int15Database,
                thermoDatabase.gaucheDatabase,
                thermoDatabase.otherDatabase,
                thermoDatabase.radicalDatabase,
                thermoDatabase.ringDatabase]:

            # All nodes in library should be in tree and dictionary
            # All nodes in tree should be in dictionary
            for node in database.library:
                self.assertTrue(node in database.tree.parent, 'Expected node "%s" from library to be in parent attribute of tree.' % node)
                self.assertTrue(node in database.tree.children, 'Expected node "%s" from library from library to be in children attribute of tree.' % node)
                self.assertTrue(node in database.dictionary, 'Expected node "%s" to be in dictionary.' % node)
            for node in database.tree.parent:
                self.assertTrue(node in database.tree.children, 'Expected node "%s" from parent attribute of tree to be in children attribute of tree.' % node)
                self.assertTrue(node in database.dictionary, 'Expected node "%s" from parent attribute of tree to be in dictionary.' % node)
            for node in database.tree.children:
                self.assertTrue(node in database.tree.parent, 'Expected node "%s" from children attribute of tree to be in parent attribute of tree.' % node)
                self.assertTrue(node in database.dictionary, 'Expected node "%s" from children attribute of tree to be in dictionary.' % node)

            # All parents in tree should be in tree
            for node in database.tree.parent:
                parentNode = database.tree.parent[node]
                if parentNode is not None:
                    self.assertTrue(parentNode in database.tree.parent)
                    self.assertTrue(parentNode in database.tree.children)

            # All children in tree should be in tree
            for node in database.tree.children:
                for childNode in database.tree.children[node]:
                    if childNode is not None:
                        self.assertTrue(childNode in database.tree.parent)
                        self.assertTrue(childNode in database.tree.children)

            # All values in dictionary should be chemical structures
            for node in database.dictionary:
                self.assertTrue(isinstance(database.dictionary[node], MoleculePattern) or isinstance(database.dictionary[node], LogicNode))

            # All values in library should be ThermoEntry objects
            for node in database.library:
                self.assertTrue(isinstance(database.library[node], ThermoEntry), '"%s" is of unexpected type "%s".' % (node, database.library[node].__class__))
                self.assertTrue(database.library[node].model is not None or database.library[node].node != '')

        molecule = Molecule(SMILES='CC')
        print
        print thermoDatabase.generateThermoData(molecule)

    def testThermoGeneration(self):
        
        # The test cases represent 1,3-hexadiene and all of its possible 
        # unimolecular bond fission products (i.e. products resulting from the
        # breaking of any single bond in 1,3-hexadiene)
        testCases = [
            # SMILES         symm  H298     S298       Cp300  Cp400  Cp500  Cp600  Cp800  Cp1000 Cp1500
            ['C=CC=CCC',     3,    13.45,   86.3671,   29.49, 37.67, 44.54, 50.12, 58.66, 64.95, 74.71],
            ['[CH]=CC=CCC',  3,    72.55,   87.7571,   29.3 , 36.92, 43.18, 48.20, 55.84, 61.46, 70.18],
            ['C=[C]C=CCC',   3,    61.15,   87.0771,   29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=C[C]=CCC',   3,    61.15,   87.0771,   29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=CC=[C]CC',   3,    70.35,   88.1771,   29.15, 36.46, 42.6 , 47.60, 55.32, 61.04, 69.95],
            ['C=CC=C[CH]C',  6,    38.24,   84.4098,   27.79, 35.46, 41.94, 47.43, 55.74, 61.92, 71.86],
            ['C=CC=CC[CH2]', 2,    62.45,   89.7827,   28.72, 36.31, 42.63, 47.72, 55.5 , 61.21, 70.05],
            #['[H]',          1,    49.00,    2.61  ,   -0.77, -1.36, -1.91, -2.4 , -3.16, -3.74, -4.66],
            ['[CH3]',        6,    34.81,   46.3698,    9.14, 10.18, 10.81, 11.34, 12.57, 13.71, 15.2 ],
            ['C=CC=C[CH2]',  2,    46.11,   75.8227,   22.54, 28.95, 34.24, 38.64, 45.14, 49.97, 57.85],
            ['[CH2]C',       6,    28.60,   59.8698,   11.73, 14.46, 17.05, 19.34, 23.02, 25.91, 31.53],
            ['C=CC=[CH]',    1,    85.18,   69.37  ,   18.93, 23.55, 27.16, 29.92, 34.02, 37.03, 41.81],
            ['C=[CH]',       1,    71.62,   56.61  ,   10.01, 11.97, 13.66, 15.08, 17.32, 19.05, 21.85],
            ['[CH]=CCC',     3,    58.99,   74.9971,   20.38, 25.34, 29.68, 33.36, 39.14, 43.48, 50.22],
        ]
               
        thermoDatabase = loadThermoDatabase('database/output/RMG_Database/thermo_groups', group=True, old=True)
        
        failMessage = '\n'
        
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            species.generateResonanceIsomers()
            thermoData = generateThermoData(species.molecule[0])
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = generateThermoData(mol)
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            
            if molecule.calculateSymmetryNumber() != symm:
                failMessage += "For %s, expected symmetry number of %g, got %g\n" % (smiles, symm, molecule.calculateSymmetryNumber())
            
            if abs(1 - thermoData.getEnthalpy(298) / 4184 / H298) > 0.001:
                failMessage += "For %s, expected H298 of %g kcal/mol, got %g kcal/mol\n" % (smiles, H298, thermoData.getEnthalpy(298) / 4184)
            if abs(1 - thermoData.getEntropy(298) / 4.184 / S298) > 0.001:
                failMessage += "For %s, expected S298 of %g cal/mol*K, got %g cal/mol*K\n" % (smiles, S298, thermoData.getEntropy(298) / 4.184)
            for T, Cp in zip([300, 400, 500, 600, 800, 1000, 1500], [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]):
                if abs(1 - thermoData.getHeatCapacity(T) / 4.184 / Cp) > 0.001:
                    failMessage += "For %s, expected Cp at %g K of %g cal/mol*K, got %g cal/mol*K\n" % (smiles, T, Cp, thermoData.getHeatCapacity(T) / 4.184)
        
        
        self.assertEqual(failMessage, '\n', failMessage)
            
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )

