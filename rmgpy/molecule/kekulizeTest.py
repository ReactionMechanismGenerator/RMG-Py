################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

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
