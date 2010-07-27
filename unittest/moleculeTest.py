#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

from chempy.molecule import *

################################################################################

class MoleculeCheck(unittest.TestCase):

    def testIsomorphism(self):
        """
        Check the graph isomorphism functions.
        """
        molecule1 = Molecule().fromSMILES('C=CC=C[CH]C')
        molecule2 = Molecule().fromSMILES('C[CH]C=CC=C')
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))

    def testSubgraphIsomorphism(self):
        """
        Check the graph isomorphism functions.
        """
        molecule = Molecule().fromSMILES('C=CC=C[CH]C')
        pattern = MoleculePattern().fromAdjacencyList("""
        1 Cd 0 {2,D}
        2 Cd 0 {1,D}
        """)
        self.assertTrue(len(molecule.resonanceForms[0].findSubgraphIsomorphisms(pattern)[1]) == 4)

    def testAdjacencyList(self):
        """
        Check the adjacency list read/write functions for a full molecule.
        """
        molecule1 = Molecule().fromAdjacencyList("""
        1 C 0       {2,D}
        2 C 0 {1,D} {3,S}
        3 C 0 {2,S} {4,D}
        4 C 0 {3,D} {5,S}
        5 C 1 {4,S} {6,S}
        6 C 0 {5,S}
        """)
        molecule2 = Molecule().fromSMILES('C=CC=C[CH]C')
        molecule1.toAdjacencyList()
        molecule2.toAdjacencyList()
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))

    def testAdjacencyListPattern(self):
        """
        Check the adjacency list read/write functions for a molecular
        substructure.
        """
        pattern1 = MoleculePattern().fromAdjacencyList("""
        1 {Cs,Os} 0  {2,S}
        2 R!H 0 {1,S}
        """)
        pattern1.toAdjacencyList()

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )