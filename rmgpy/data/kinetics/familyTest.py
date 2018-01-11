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
import os.path
import shutil

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import ReactionRecipe
from rmgpy.molecule import Molecule
import filecmp
###################################################

class TestFamily(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        # Set up a dummy database
        cls.database = KineticsDatabase()
        cls.database.loadFamilies(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families'),
                                  families=['intra_H_migration', 'R_Addition_MultipleBond'])
        cls.family = cls.database.families['intra_H_migration']

    def testGetBackboneRoots(self):
        """
        Test the getBackboneRoots() function
        """
        backbones = self.family.getBackboneRoots()
        self.assertEquals(backbones[0].label, "RnH")

    def testGetEndRoots(self):
        """
        Test the getEndRoots() function
        """
        ends = self.family.getEndRoots()
        self.assertEquals(len(ends), 2)
        self.assertIn(self.family.groups.entries["Y_rad_out"], ends)
        self.assertIn(self.family.groups.entries["XH_out"], ends)

    def testGetTopLevelGroups(self):
        """
        Test the getTopLevelGroups() function
        """
        topGroups = self.family.getTopLevelGroups(self.family.groups.entries["RnH"])
        self.assertEquals(len(topGroups), 2)
        self.assertIn(self.family.groups.entries["R5Hall"], topGroups)
        self.assertIn(self.family.groups.entries["R6Hall"], topGroups)

    def testReactBenzeneBond(self):
        """
        Test that hydrogen addition to benzene (w/ benzene bonds) returns kekulized product.
        """
        family = self.database.families['R_Addition_MultipleBond']
        reactants = [Molecule().fromAdjacencyList("""
1  *1 C u0 p0 c0 {2,B} {6,B} {7,S}
2  *2 C u0 p0 c0 {1,B} {3,B} {8,S}
3     C u0 p0 c0 {2,B} {4,B} {9,S}
4     C u0 p0 c0 {3,B} {5,B} {10,S}
5     C u0 p0 c0 {4,B} {6,B} {11,S}
6     C u0 p0 c0 {1,B} {5,B} {12,S}
7     H u0 p0 c0 {1,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {4,S}
11    H u0 p0 c0 {5,S}
12    H u0 p0 c0 {6,S}
"""),
                     Molecule().fromAdjacencyList("1 *3 H u1 p0 c0")]
        expectedProduct = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
2  C u1 p0 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {5,D} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)
        self.assertTrue(expectedProduct.isIsomorphic(products[0]))

    def testReactBenzeneBond2(self):
        """
        Test that hydrogen addition to phenanthrene (w/ benzene bonds) returns kekulized product.
        """
        family = self.database.families['R_Addition_MultipleBond']
        reactants = [Molecule().fromAdjacencyList("""
1  *1 C u0 p0 c0 {2,B} {3,B} {6,B}
2  *2 C u0 p0 c0 {1,B} {4,B} {9,B}
3     C u0 p0 c0 {1,B} {5,B} {7,B}
4     C u0 p0 c0 {2,B} {8,B} {10,B}
5     C u0 p0 c0 {3,B} {11,B} {17,S}
6     C u0 p0 c0 {1,B} {12,B} {18,S}
7     C u0 p0 c0 {3,B} {8,B} {19,S}
8     C u0 p0 c0 {4,B} {7,B} {20,S}
9     C u0 p0 c0 {2,B} {13,B} {21,S}
10    C u0 p0 c0 {4,B} {14,B} {23,S}
11    C u0 p0 c0 {5,B} {12,B} {15,S}
12    C u0 p0 c0 {6,B} {11,B} {16,S}
13    C u0 p0 c0 {9,B} {14,B} {22,S}
14    C u0 p0 c0 {10,B} {13,B} {24,S}
15    H u0 p0 c0 {11,S}
16    H u0 p0 c0 {12,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {7,S}
20    H u0 p0 c0 {8,S}
21    H u0 p0 c0 {9,S}
22    H u0 p0 c0 {13,S}
23    H u0 p0 c0 {10,S}
24    H u0 p0 c0 {14,S}
"""),
                     Molecule().fromAdjacencyList("1 *3 H u1 p0 c0")]
        expectedProduct = Molecule().fromAdjacencyList("""
multiplicity 2
1  *1 C u0 p0 c0 {2,S} {3,S} {5,S} {15,S}
2  *2 C u1 p0 c0 {1,S} {4,S} {8,S}
3     C u0 p0 c0 {1,S} {6,S} {7,D}
4     C u0 p0 c0 {2,S} {9,D} {10,S}
5     C u0 p0 c0 {1,S} {11,D} {16,S}
6     C u0 p0 c0 {3,S} {12,D} {19,S}
7     C u0 p0 c0 {3,D} {9,S} {20,S}
8     C u0 p0 c0 {2,S} {13,D} {22,S}
9     C u0 p0 c0 {4,D} {7,S} {21,S}
10    C u0 p0 c0 {4,S} {14,D} {24,S}
11    C u0 p0 c0 {5,D} {12,S} {18,S}
12    C u0 p0 c0 {6,D} {11,S} {17,S}
13    C u0 p0 c0 {8,D} {14,S} {23,S}
14    C u0 p0 c0 {10,D} {13,S} {25,S}
15 *3 H u0 p0 c0 {1,S}
16    H u0 p0 c0 {5,S}
17    H u0 p0 c0 {12,S}
18    H u0 p0 c0 {11,S}
19    H u0 p0 c0 {6,S}
20    H u0 p0 c0 {7,S}
21    H u0 p0 c0 {9,S}
22    H u0 p0 c0 {8,S}
23    H u0 p0 c0 {13,S}
24    H u0 p0 c0 {10,S}
25    H u0 p0 c0 {14,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)
        self.assertTrue(expectedProduct.isIsomorphic(products[0]))

    def testSaveFamily(self):
        """

        This tests the the family.save method by writing a new temporary file and
        comparing it to the original source.

        """
        os.makedirs(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_copy'))
        self.family.save(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_copy'))
        try:
            self.assertTrue(filecmp.cmp(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_migration/groups.py'),
                                 os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_copy/groups.py')))
            self.assertTrue(filecmp.cmp(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_migration/rules.py'),
                                 os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_copy/rules.py')))
            self.assertTrue(filecmp.cmp(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_migration/training/reactions.py'),
                                 os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_copy/training/reactions.py')))
            self.assertTrue(filecmp.cmp(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_migration/training/dictionary.txt'),
                                 os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_copy/training/dictionary.txt')))
        finally:
            shutil.rmtree(os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families/intra_H_copy'))