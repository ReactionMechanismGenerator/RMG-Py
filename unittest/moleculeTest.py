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

        self.assertTrue(molecule.isSubgraphIsomorphic(pattern))
        match, mapping = molecule.findSubgraphIsomorphisms(pattern)
        self.assertTrue(match)
        self.assertTrue(len(mapping) == 4, "len(mapping) = %d, should be = 4" % (len(mapping)))
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(pattern.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in pattern.atoms)

    def testSubgraphIsomorphismAgain(self):
        molecule = Molecule()
        molecule.fromAdjacencyList("""
        1 * C 0 {2,D} {7,S} {8,S}
        2 C 0 {1,D} {3,S} {9,S}
        3 C 0 {2,S} {4,D} {10,S}
        4 C 0 {3,D} {5,S} {11,S}
        5 C 0 {4,S} {6,S} {12,S} {13,S}
        6 C 0 {5,S} {14,S} {15,S} {16,S}
        7 H 0 {1,S}
        8 H 0 {1,S}
        9 H 0 {2,S}
        10 H 0 {3,S}
        11 H 0 {4,S}
        12 H 0 {5,S}
        13 H 0 {5,S}
        14 H 0 {6,S}
        15 H 0 {6,S}
        16 H 0 {6,S}
        """)
        
        pattern = MoleculePattern()
        pattern.fromAdjacencyList("""
        1 * C 0 {2,D} {3,S} {4,S}
        2   C 0 {1,D}
        3   H 0 {1,S}
        4   H 0 {1,S}
        """)

        molecule.makeHydrogensExplicit()
        
        labeled1 = molecule.getLabeledAtoms().values()[0]
        labeled2 = pattern.getLabeledAtoms().values()[0]
        
        initialMap = {labeled1: labeled2}
        self.assertTrue(molecule.isSubgraphIsomorphic(pattern, initialMap))

        initialMap = {labeled1: labeled2}
        match, mapping = molecule.findSubgraphIsomorphisms(pattern, initialMap)
        self.assertTrue(match)
        self.assertTrue(len(mapping) == 2,  "len(mapping) = %d, should be = 2" % (len(mapping)))
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(pattern.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in pattern.atoms)

    def testSubgraphIsomorphismManyLabels(self):
        molecule = Molecule() # specific case (species)
        molecule.fromAdjacencyList("""
1 *1 C  1 {2,S} {3,S}
2    C  0 {1,S} {3,S}
3    C  0 {1,S} {2,S}
        """)
        
        pattern = MoleculePattern() # general case (functional group)
        pattern.fromAdjacencyList("""
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
        """)

        labeled1 = molecule.getLabeledAtoms()
        labeled2 = pattern.getLabeledAtoms()
        initialMap = {}
        for label,atom1 in labeled1.iteritems():
            initialMap[atom1] = labeled2[label]
        self.assertTrue(molecule.isSubgraphIsomorphic(pattern, initialMap))

        match, mapping = molecule.findSubgraphIsomorphisms(pattern, initialMap)
        self.assertTrue(match)
        self.assertTrue(len(mapping) == 1)
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(pattern.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in pattern.atoms)

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

    def testSSSR(self):
        """
        Check the graph's Smallest Set of Smallest Rings function
        """
        molecule = Molecule()
        molecule.fromSMILES('C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC')
        #http://cactus.nci.nih.gov/chemical/structure/C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC/image
        sssr = molecule.getSmallestSetOfSmallestRings()
        self.assertEqual( len(sssr), 3)

    def testIsInCycle(self):

        # ethane
        molecule = Molecule().fromSMILES('CC')
        for atom in molecule.atoms:
            self.assertFalse(molecule.isAtomInCycle(atom))
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                self.assertFalse(molecule.isBondInCycle(atom1, atom2))

        # cyclohexane
        molecule = Molecule().fromInChI('InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2')
        for atom in molecule.atoms:
            if atom.isHydrogen():
                self.assertFalse(molecule.isAtomInCycle(atom))
            elif atom.isCarbon():
                self.assertTrue(molecule.isAtomInCycle(atom))
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                if atom1.isCarbon() and atom2.isCarbon():
                    self.assertTrue(molecule.isBondInCycle(atom1, atom2))
                else:
                    self.assertFalse(molecule.isBondInCycle(atom1, atom2))

    def testRotorNumber(self):
        """Count the number of internal rotors"""
        # http://cactus.nci.nih.gov/chemical/structure/C1CCCC1C/image
        test_set = [('CC', 1),
                    ('CCC', 2),
                    ('CC(C)(C)C', 4),
                    ('C1CCCC1C',1),
                    ('C=C',0)
                    ]
        fail_message = ''
        for smile,should_be in test_set:
            molecule = Molecule(SMILES=smile)
            molecule.makeHydrogensExplicit()
            rotorNumber = molecule.countInternalRotors()
            if rotorNumber!=should_be:
                fail_message+="Got rotor number of %s for %s (expected %s)\n"%(rotorNumber,smile,should_be)
        self.assertEqual(fail_message,'',fail_message)

    def testRotorNumberHard(self):
        """Count the number of internal rotors in a tricky case"""
        test_set = [('CC', 1),   # start with something simple:    H3C---CH3
                    ('CC#CC', 1) # now lengthen that middle bond: H3C-C#C-CH3
                    ]
        fail_message = ''
        for smile,should_be in test_set:
            molecule = Molecule(SMILES=smile)
            molecule.makeHydrogensExplicit()
            rotorNumber = molecule.countInternalRotors()
            if rotorNumber!=should_be:
                fail_message+="Got rotor number of %s for %s (expected %s)\n"%(rotorNumber,smile,should_be)
        self.assertEqual(fail_message,'',fail_message)

    def testLinear(self):
        """Identify linear molecules"""
        # http://cactus.nci.nih.gov/chemical/structure/C1CCCC1C/image
        test_set = [('CC', False),
                    ('CCC', False),
                    ('CC(C)(C)C', False),
                    ('C',False),
                    ('[H]',False),
                    ('O=O',True),
                    #('O=S',True),
                    ('O=C=O',True),
                    ('C#C', True),
                    ('C#CC#CC#C', True)
                    ]
        fail_message = ''
        for smile,should_be in test_set:
            molecule = Molecule(SMILES=smile)
            molecule.makeHydrogensExplicit()
            symmetryNumber = molecule.isLinear()
            if symmetryNumber!=should_be:
                fail_message+="Got linearity %s for %s (expected %s)\n"%(symmetryNumber,smile,should_be)
        self.assertEqual(fail_message,'',fail_message)

    def testH(self):
        """
        Make sure that H radicals are produced properly from various shorthands.
        """

        # InChI
        molecule = Molecule(InChI='InChI=1/H')
        self.assertTrue(len(molecule.atoms) == 1)
        H = molecule.atoms[0]
        self.assertTrue(H.isHydrogen())
        self.assertTrue(H.radicalElectrons == 1)

        # SMILES
        molecule = Molecule(SMILES='[H]')
        self.assertTrue(len(molecule.atoms) == 1)
        H = molecule.atoms[0]
        print repr(H)
        self.assertTrue(H.isHydrogen())
        self.assertTrue(H.radicalElectrons == 1)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )