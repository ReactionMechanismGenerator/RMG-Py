#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

from rmgdata.thermo import loadThermoDatabase, ThermoEntry
from rmgdata.base import LogicNode
from chempy.molecule import Molecule
from chempy.pattern import MoleculePattern
from chempy.thermo import ThermoGAModel

################################################################################

class ThermoDatabaseCheck(unittest.TestCase):

    def testOldThermoDatabase(self):
        """
        Check the database load functions.
        """

        thermoDatabase = loadThermoDatabase('output/RMG_Database', old=True)

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

        thermoDatabase = loadThermoDatabase('input/RMG_database', old=False)

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

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )

