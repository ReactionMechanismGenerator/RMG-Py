import unittest
from external.wip import work_in_progress

from rmgpy.molecule import Molecule
from rmgpy.molecule.kekulize import *

class KekulizeTest(unittest.TestCase):

    def setUp(self):
        """To be run before each test."""
        molecule = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {6,B} {7,S}
2  C u0 p0 c0 {1,B} {3,B} {8,S}
3  C u0 p0 c0 {2,B} {4,B} {9,S}
4  C u0 p0 c0 {3,B} {5,B} {10,S}
5  C u0 p0 c0 {4,B} {6,B} {11,S}
6  C u0 p0 c0 {1,B} {5,B} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""")
        bonds = set()
        for atom in molecule.atoms:
            bonds.update(atom.bonds.values())

        ringAtoms, ringBonds = molecule.getAromaticRings()

        self.aromaticRing = AromaticRing(ringAtoms[0], set(ringBonds[0]), bonds - set(ringBonds[0]))

    def testAromaticRing(self):
        self.aromaticRing.update()

        self.assertEqual(self.aromaticRing.endoDOF, 6)
        self.assertEqual(self.aromaticRing.exoDOF, 0)

        result = self.aromaticRing.kekulize()

        self.assertTrue(result)

    def testAromaticBond(self):
        resolved, unresolved = self.aromaticRing.processBonds()

        self.assertEqual(len(resolved), 0)
        self.assertEqual(len(unresolved), 6)

        for bond in unresolved:
            bond.update()
            self.assertEqual(bond.endoDOF, 2)
            self.assertEqual(bond.exoDOF, 0)
            self.assertTrue(bond.doublePossible)
            self.assertFalse(bond.doubleRequired)
