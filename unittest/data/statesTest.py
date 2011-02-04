#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

from rmgpy.data.states import loadFrequencyDatabase, GroupFrequency
from rmgpy.data.base import LogicNode
from rmgpy.chem.molecule import Molecule
from rmgpy.chem.pattern import MoleculePattern

################################################################################

class ThermoDatabaseCheck(unittest.TestCase):

    def testThermoDatabase(self):
        """
        Check the database load functions.
        """

        frequencyDatabase = loadFrequencyDatabase('database/output/RMG_Database')
        
        # All nodes in library should be in tree and dictionary
        # All nodes in tree should be in dictionary
        database = frequencyDatabase.groupDatabase
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

        # All values in library should be ThermoGAModel objects or lists of length 2
        for node in database.library:
            self.assertTrue(isinstance(database.library[node], GroupFrequency) or (isinstance(database.library[node], list) and len(database.library[node]) == 2), '"%s" is of unexpected type "%s".' % (node, database.library[node].__class__))

        molecule = Molecule(SMILES='CC')
        print
        print frequencyDatabase.getFrequencies(molecule)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )

