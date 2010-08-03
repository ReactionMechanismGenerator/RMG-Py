#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

import logging

from rmgdata.kinetics import loadKineticsDatabase, KineticsEntry
from rmgdata.base import LogicNode
from chempy.molecule import Molecule
from chempy.pattern import MoleculePattern
from chempy.reaction import Reaction
from chempy.kinetics import ArrheniusEPModel

################################################################################

class KineticsDatabaseCheck(unittest.TestCase):

    def testKineticsDatabase(self):
        """
        Check the database load functions.
        """

        # Create logger
        logger = logging.getLogger()
        logger.setLevel(10)
        # Create console handler and set level to debug; send everything to stdout
        # rather than stderr
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(10)
        # Create formatter and add to console handler
        formatter = logging.Formatter('%(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        # remove old handlers!
        while logger.handlers: logger.removeHandler(logger.handlers[0])
        # Add ch to logger
        logger.addHandler(ch)
	
        kineticsDatabase = loadKineticsDatabase('output/RMG_Database')

        for label, family in kineticsDatabase.families.iteritems():

            # All nodes in library should be in tree and dictionary
            # All nodes in tree should be in dictionary
#            for nodes in family.library:
#                for node in nodes.split(';'):
#                    self.assertTrue(node in family.tree.parent)
#                    self.assertTrue(node in family.tree.children)
#                    self.assertTrue(node in family.dictionary)
            for node in family.tree.parent:
                self.assertTrue(node in family.tree.children)
                self.assertTrue(node in family.dictionary)
            for node in family.tree.children:
                self.assertTrue(node in family.tree.parent)
                self.assertTrue(node in family.dictionary)

            # All parents in tree should be in tree
            for node in family.tree.parent:
                parentNode = family.tree.parent[node]
                if parentNode is not None:
                    self.assertTrue(parentNode in family.tree.parent)
                    self.assertTrue(parentNode in family.tree.children)

            # All children in tree should be in tree
            for node in family.tree.children:
                for childNode in family.tree.children[node]:
                    if childNode is not None:
                        self.assertTrue(childNode in family.tree.parent)
                        self.assertTrue(childNode in family.tree.children)

            # All values in dictionary should be chemical structures
            for node in family.dictionary:
                self.assertTrue(isinstance(family.dictionary[node], MoleculePattern) or isinstance(family.dictionary[node], LogicNode))

            # All values in library should be ArrheniusEPModel objects or lists of length 2
            for node in family.library:
                self.assertTrue(isinstance(family.library[node], KineticsEntry), '"%s" is of unexpected type "%s".' % (node, family.library[node].__class__))
                self.assertTrue(family.library[node].model is not None or family.library[node].node != '')
                
        C2H4 = Molecule().fromAdjacencyList("""
        1 *1 C 0 {2,D} {3,S} {4,S}
        2 *2 C 0 {1,D} {5,S} {6,S}
        3    H 0 {1,S}
        4    H 0 {1,S}
        5    H 0 {2,S}
        6    H 0 {2,S}
        """)
        H = Molecule().fromAdjacencyList("""
        1 *3 H 1
        """)
        C2H5 = Molecule().fromAdjacencyList("""
        1 *1 C 0 {2,S} {3,S} {4,S} {7,S}
        2 *2 C 1 {1,S} {5,S} {6,S}
        3    H 0 {1,S}
        4    H 0 {1,S}
        5    H 0 {2,S}
        6    H 0 {2,S}
        7 *3 H 0 {1,S}
        """)
        rxn = Reaction(reactants=[C2H4, H], products=[C2H5])
        print 'GENERATING KINETICS DATA:'
        print kineticsDatabase.generateKineticsData(rxn, family='R_Addition_MultipleBond')

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )

