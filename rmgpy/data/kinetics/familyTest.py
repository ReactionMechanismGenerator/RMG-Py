#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import filecmp
import mock
import os.path
import shutil
import unittest
import numpy as np
import logging

from rmgpy import settings
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.molecule.group import Group


###################################################

class TestFamily(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        # Set up a dummy database
        cls.database = KineticsDatabase()
        cls.database.loadFamilies(
            path=os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families'),
            families=[
                'intra_H_migration',
                'R_Addition_MultipleBond',
                'H_Abstraction',
                'Intra_ene_reaction',
                '6_membered_central_C-C_shift',
                '1,2_shiftC',
                'Intra_R_Add_Exo_scission',
                'intra_substitutionS_isomerization',
                'R_Addition_COm'
            ],
        )
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

    def test_intra_H_migration(self):
        """
        Test that the intra_H_migration family, which is its own reverse, returns a properly re-labeled product structure.
        """
        family = self.database.families['intra_H_migration']
        reactants = [Molecule().fromAdjacencyList("""
multiplicity 2
1  *2 C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
2  *4 C u0 p0 c0 {4,S} {5,S} {6,D}
3  *5 C u0 p0 c0 {1,S} {7,D} {14,S}
4  *1 C u1 p0 c0 {2,S} {8,S} {15,S}
5     C u0 p0 c0 {2,S} {10,D} {17,S}
6  *6 C u0 p0 c0 {2,D} {7,S} {19,S}
7  *7 C u0 p0 c0 {3,D} {6,S} {21,S}
8     C u0 p0 c0 {4,S} {9,D} {16,S}
9     C u0 p0 c0 {8,D} {10,S} {20,S}
10    C u0 p0 c0 {5,D} {9,S} {18,S}
11 *3 H u0 p0 c0 {1,S}
12    H u0 p0 c0 {1,S}
13    H u0 p0 c0 {1,S}
14    H u0 p0 c0 {3,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {8,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {10,S}
19    H u0 p0 c0 {6,S}
20    H u0 p0 c0 {9,S}
21    H u0 p0 c0 {7,S}
""")]
        expectedProduct = Molecule().fromAdjacencyList("""
multiplicity 2
1  *1 C u1 p0 c0 {3,S} {12,S} {13,S}
2  *5 C u0 p0 c0 {4,S} {5,S} {6,D}
3  *4 C u0 p0 c0 {1,S} {7,D} {14,S}
4  *2 C u0 p0 c0 {2,S} {11,S} {8,S} {15,S}
5     C u0 p0 c0 {2,S} {10,D} {17,S}
6  *7 C u0 p0 c0 {2,D} {7,S} {19,S}
7  *6 C u0 p0 c0 {3,D} {6,S} {21,S}
8     C u0 p0 c0 {4,S} {9,D} {16,S}
9     C u0 p0 c0 {8,D} {10,S} {20,S}
10    C u0 p0 c0 {5,D} {9,S} {18,S}
11 *3 H u0 p0 c0 {4,S}
12    H u0 p0 c0 {1,S}
13    H u0 p0 c0 {1,S}
14    H u0 p0 c0 {3,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {8,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {10,S}
19    H u0 p0 c0 {6,S}
20    H u0 p0 c0 {9,S}
21    H u0 p0 c0 {7,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)

        mapping = {}
        for label, atom in expectedProduct.getLabeledAtoms().iteritems():
            mapping[atom] = products[0].getLabeledAtom(label)

        self.assertTrue(expectedProduct.isIsomorphic(products[0], mapping))

    def test_H_Abstraction(self):
        """
        Test that the H_Abstraction family, which is its own reverse, returns a properly re-labeled product structure.
        """
        family = self.database.families['H_Abstraction']
        reactants = [Molecule().fromAdjacencyList("""
1 *1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2    C u0 p0 c0 {1,S} {3,D} {7,S}
3    C u0 p0 c0 {2,D} {8,S} {9,S}
4 *2 H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {3,S}
        """),
                     Molecule().fromAdjacencyList("1 *3 H u1 p0 c0")]
        expectedProducts = [Molecule().fromAdjacencyList("""
1 *1 H u0 p0 c0 {2,S}
2 *2 H u0 p0 c0 {1,S}
        """),
                    Molecule().fromAdjacencyList("""
1 *3 C u1 p0 c0 {2,S} {5,S} {6,S}
2    C u0 p0 c0 {1,S} {3,D} {7,S}
3    C u0 p0 c0 {2,D} {8,S} {9,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {3,S}
        """)]
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 2)

        mapping1 = {}
        for label, atom in expectedProducts[0].getLabeledAtoms().iteritems():
            mapping1[atom] = products[0].getLabeledAtom(label)

        self.assertTrue(expectedProducts[0].isIsomorphic(products[0], mapping1))

        mapping2 = {}
        for label, atom in expectedProducts[1].getLabeledAtoms().iteritems():
            mapping2[atom] = products[1].getLabeledAtom(label)

        self.assertTrue(expectedProducts[1].isIsomorphic(products[1], mapping2))

    def test_Intra_ene_reaction(self):
        """
        Test that the Intra_ene_reaction family, which is its own reverse, returns a properly re-labeled product structure.
        """
        family = self.database.families['Intra_ene_reaction']
        reactants = [Molecule().fromAdjacencyList("""
1  *1 C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  *5 C u0 p0 c0 {1,S} {5,D} {6,S}
3  *2 C u0 p0 c0 {1,S} {7,D} {11,S}
4     C u0 p0 c0 {1,S} {8,D} {12,S}
5  *4 C u0 p0 c0 {2,D} {7,S} {13,S}
6     C u0 p0 c0 {2,S} {9,D} {15,S}
7  *3 C u0 p0 c0 {3,D} {5,S} {14,S}
8     C u0 p0 c0 {4,D} {9,S} {17,S}
9     C u0 p0 c0 {6,D} {8,S} {16,S}
10 *6 H u0 p0 c0 {1,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
14    H u0 p0 c0 {7,S}
15    H u0 p0 c0 {6,S}
16    H u0 p0 c0 {9,S}
17    H u0 p0 c0 {8,S}
""")]
        expectedProduct = Molecule().fromAdjacencyList("""
1  *2 C u0 p0 c0 {2,D} {3,S} {4,S} 
2  *3 C u0 p0 c0 {1,D} {5,S} {6,S}
3  *1 C u0 p0 c0 {1,S} {7,S} {11,S} {10,S}
4     C u0 p0 c0 {1,S} {8,D} {12,S}
5  *4 C u0 p0 c0 {2,S} {7,D} {13,S}
6     C u0 p0 c0 {2,S} {9,D} {15,S}
7  *5 C u0 p0 c0 {3,S} {5,D} {14,S}
8     C u0 p0 c0 {4,D} {9,S} {17,S}
9     C u0 p0 c0 {6,D} {8,S} {16,S}
10 *6 H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
14    H u0 p0 c0 {7,S}
15    H u0 p0 c0 {6,S}
16    H u0 p0 c0 {9,S}
17    H u0 p0 c0 {8,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)

        mapping = {}
        for label, atom in expectedProduct.getLabeledAtoms().iteritems():
            mapping[atom] = products[0].getLabeledAtom(label)

        self.assertTrue(expectedProduct.isIsomorphic(products[0], mapping))

    def test_6_membered_central_CC_shift(self):
        """
        Test that the 6_membered_central_C-C_shift family, which is its own reverse, returns a properly re-labeled product structure.
        """
        family = self.database.families['6_membered_central_C-C_shift']
        reactants = [Molecule().fromAdjacencyList("""
1  *3 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  *4 C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  *2 C u0 p0 c0 {1,S} {5,T}
4  *5 C u0 p0 c0 {2,S} {6,T}
5  *1 C u0 p0 c0 {3,T} {11,S}
6  *6 C u0 p0 c0 {4,T} {12,S}
7     H u0 p0 c0 {1,S}
8     H u0 p0 c0 {1,S}
9     H u0 p0 c0 {2,S}
10    H u0 p0 c0 {2,S}
11    H u0 p0 c0 {5,S}
12    H u0 p0 c0 {6,S}
""")]
        expectedProduct = Molecule().fromAdjacencyList("""
1  *3 C u0 p0 c0 {2,S} {5,D} {7,S}
2  *4 C u0 p0 c0 {1,S} {6,D} {8,S}
3  *1 C u0 p0 c0 {5,D} {9,S} {10,S}
4  *6 C u0 p0 c0 {6,D} {11,S} {12,S}
5  *2 C u0 p0 c0 {1,D} {3,D}
6  *5 C u0 p0 c0 {2,D} {4,D}
7     H u0 p0 c0 {1,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {4,S}
12    H u0 p0 c0 {4,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)

        mapping = {}
        for label, atom in expectedProduct.getLabeledAtoms().iteritems():
            mapping[atom] = products[0].getLabeledAtom(label)

        self.assertTrue(expectedProduct.isIsomorphic(products[0], mapping))

    def test_12_shiftC(self):
        """
        Test that the 1,2_shiftC family, which is its own reverse, returns a properly re-labeled product structure.
        """
        family = self.database.families['1,2_shiftC']
        reactants = [Molecule().fromAdjacencyList("""
multiplicity 2
1  *2 C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  *1 C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
3  *3 C u1 p0 c0 {1,S} {4,S} {5,S}
4     C u0 p0 c0 {3,S} {6,D} {13,S}
5     C u0 p0 c0 {3,S} {7,D} {14,S}
6     C u0 p0 c0 {4,D} {7,S} {15,S}
7     C u0 p0 c0 {5,D} {6,S} {16,S}
8     H u0 p0 c0 {1,S}
9     H u0 p0 c0 {1,S}
10    H u0 p0 c0 {2,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {2,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {5,S}
15    H u0 p0 c0 {6,S}
16    H u0 p0 c0 {7,S}
""")]
        expectedProduct = Molecule().fromAdjacencyList("""
multiplicity 2
1  *2 C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  *1 C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3     C u0 p0 c0 {1,S} {5,D} {11,S}
4     C u0 p0 c0 {1,S} {6,D} {12,S}
5     C u0 p0 c0 {3,D} {6,S} {13,S}
6     C u0 p0 c0 {4,D} {5,S} {14,S}
7  *3 C u1 p0 c0 {1,S} {15,S} {16,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {2,S}
10    H u0 p0 c0 {2,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
14    H u0 p0 c0 {6,S}
15    H u0 p0 c0 {7,S}
16    H u0 p0 c0 {7,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)

        mapping = {}
        for label, atom in expectedProduct.getLabeledAtoms().iteritems():
            mapping[atom] = products[0].getLabeledAtom(label)

        self.assertTrue(expectedProduct.isIsomorphic(products[0], mapping))

    def test_Intra_R_Add_Exo_scission(self):
        """
        Test that the Intra_R_Add_Exo_scission family, which is its own reverse, returns a properly re-labeled product structure.
        """
        family = self.database.families['Intra_R_Add_Exo_scission']
        reactants = [Molecule().fromAdjacencyList("""
multiplicity 2
1  *3 C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
2  *2 C u0 p0 c0 {1,S} {3,B} {4,B}
3     C u0 p0 c0 {2,B} {5,B} {13,S}
4     C u0 p0 c0 {2,B} {7,B} {17,S}
5     C u0 p0 c0 {3,B} {6,B} {14,S}
6     C u0 p0 c0 {5,B} {7,B} {15,S}
7     C u0 p0 c0 {4,B} {6,B} {16,S}
8  *1 C u1 p0 c0 {1,S} {9,S} {18,S}
9     C u0 p0 c0 {8,S} {10,T}
10    C u0 p0 c0 {9,T} {19,S}
11    H u0 p0 c0 {1,S}
12    H u0 p0 c0 {1,S}
13    H u0 p0 c0 {3,S}
14    H u0 p0 c0 {5,S}
15    H u0 p0 c0 {6,S}
16    H u0 p0 c0 {7,S}
17    H u0 p0 c0 {4,S}
18    H u0 p0 c0 {8,S}
19    H u0 p0 c0 {10,S}
""")]
        expectedProduct = Molecule().fromAdjacencyList("""
multiplicity 2
1  *3 C u0 p0 c0 {2,S} {8,S} {9,S} {11,S}
2  *2 C u0 p0 c0 {1,S} {3,B} {4,B}
3     C u0 p0 c0 {2,B} {5,B} {12,S}
4     C u0 p0 c0 {2,B} {7,B} {16,S}
5     C u0 p0 c0 {3,B} {6,B} {13,S}
6     C u0 p0 c0 {5,B} {7,B} {14,S}
7     C u0 p0 c0 {4,B} {6,B} {15,S}
8  *1 C u1 p0 c0 {1,S} {17,S} {18,S}
9     C u0 p0 c0 {1,S} {10,T}
10    C u0 p0 c0 {9,T} {19,S}
11    H u0 p0 c0 {1,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {5,S}
14    H u0 p0 c0 {6,S}
15    H u0 p0 c0 {7,S}
16    H u0 p0 c0 {4,S}
17    H u0 p0 c0 {8,S}
18    H u0 p0 c0 {8,S}
19    H u0 p0 c0 {10,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)

        mapping = {}
        for label, atom in expectedProduct.getLabeledAtoms().iteritems():
            mapping[atom] = products[0].getLabeledAtom(label)

        self.assertTrue(expectedProduct.isIsomorphic(products[0], mapping))

    def test_intra_substitutionS_isomerization(self):
        """
        Test that the intra_substitutionS_isomerization family, which is its own reverse, returns a properly re-labeled product structure.
        """
        family = self.database.families['intra_substitutionS_isomerization']
        reactants = [Molecule().fromAdjacencyList("""
multiplicity 2
1  *2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2     C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  *3 C u1 p0 c0 {1,S} {2,S} {10,S}
4  *1 S u0 p2 c0 {1,S} {11,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {1,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {2,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {4,S}
""")]
        expectedProduct = Molecule().fromAdjacencyList("""
multiplicity 2
1  *2 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  *3 C u1 p0 c0 {1,S} {9,S} {10,S}
4  *1 S u0 p2 c0 {1,S} {11,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {4,S}
""")
        products = family.applyRecipe(reactants)

        self.assertEqual(len(products), 1)

        mapping = {}
        for label, atom in expectedProduct.getLabeledAtoms().iteritems():
            mapping[atom] = products[0].getLabeledAtom(label)

        self.assertTrue(expectedProduct.isIsomorphic(products[0], mapping))

    def test_r_addition_com(self):
            """
            Test that the R_addition_COm family, whose product template is generated by
            charged groups, can successfully match the reaction and returns properly product structures.
            """
            family = self.database.families['R_Addition_COm']
            reactants = [Molecule().fromAdjacencyList("""
1 *1  C u0 p1 c-1 {2,T}
2 *3  O u0 p1 c+1 {1,T}
"""),
                         Molecule().fromAdjacencyList("""
multiplicity 2
1      C u0 p0 c0 {2,D} {7,S} {8,S}
2      C u0 p0 c0 {1,D} {3,S} {9,S}
3      C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  *2  C u1 p0 c0 {3,S} {5,S} {6,S}
5      H u0 p0 c0 {4,S}
6      H u0 p0 c0 {4,S}
7      H u0 p0 c0 {1,S}
8      H u0 p0 c0 {1,S}
9      H u0 p0 c0 {2,S}
10     H u0 p0 c0 {3,S}
11     H u0 p0 c0 {3,S}
"""),
]

            expected_products = [Molecule().fromAdjacencyList("""
multiplicity 2
1      C u0 p0 c0 {2,D} {7,S} {8,S}
2      C u0 p0 c0 {1,D} {3,S} {9,S}
3      C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  *2  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  *1  C u1 p0 c0 {4,S} {6,D}
6  *3  O u0 p2 c0 {5,D}
7      H u0 p0 c0 {1,S}
8      H u0 p0 c0 {1,S}
9      H u0 p0 c0 {2,S}
10     H u0 p0 c0 {3,S}
11     H u0 p0 c0 {3,S}
12     H u0 p0 c0 {4,S}
13     H u0 p0 c0 {4,S}
"""),
]

            products = family.applyRecipe(reactants)

            self.assertEqual(len(products), 1)
            self.assertTrue(expected_products[0].isIsomorphic(products[0]))


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


class TestTreeGeneration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function run ONCE before all unit tests in this class."""
        # Set up a dummy database
        cls.database = RMGDatabase()
        cls.database.load(
            path=os.path.join(settings['test_data.directory'], 'testing_database'),
            thermoLibraries=[],
            reactionLibraries=[],
            kineticsFamilies=[],
            depository=False,
            solvation=False,
            testing=True,
        )
        cls.database.loadForbiddenStructures()
        
        cls.thermoDatabase = ThermoDatabase() #the real full Thermo Database
        cls.thermoDatabase.load(path=os.path.join(settings['database.directory'],'thermo'),libraries=['primaryThermoLibrary'])
        
        cls.kineticsDatabase = KineticsDatabase()
        cls.kineticsDatabase.loadFamilies(
            path=os.path.join(settings['test_data.directory'], 'testing_database/kinetics/families'),
            families=[
                'Singlet_Carbene_Intra_Disproportionation',
            ],
        )
        cls.family = cls.kineticsDatabase.families['Singlet_Carbene_Intra_Disproportionation']

    @classmethod
    def tearDownClass(cls):
        """A function run ONCE after all unit tests in this class."""
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    
    def test_AClearTree(self):
        """
        Test that the tree was properly cleared before generation
        """
        self.family.prepareTreeForGeneration(self.thermoDatabase)
        ents = [ent for ent in self.family.groups.entries.itervalues() if ent.index != -1]
        self.assertEquals(len(ents),1,'more than one relevant group left in groups after preparing tree for generation')
        self.assertEquals(len(self.family.rules.entries),1,'more than one group in rules.entries after preparing tree for generation' )
        root = self.family.groups.entries[self.family.rules.entries.keys()[0]]
        self.assertEquals([root],self.family.forwardTemplate.reactants)
        self.assertEquals([root],self.family.groups.top)
        
    def test_BGenerateTree(self):
        """
        test tree generation process
        """
        def objective(k1s,k2s):
            return len(k1s)*np.std(k1s)+len(k2s)*np.std(k2s)
        
        self.family.generateTree(thermoDatabase=self.thermoDatabase,obj=objective) #test input objective function
        
        self.family.prepareTreeForGeneration(self.thermoDatabase) #reclear
        
        self.family.generateTree(thermoDatabase=self.thermoDatabase) #test that default objective works
        
    def test_CParentChild(self):
        """
        test that the tree is structured properly
        """
        for entry in self.family.groups.entries.itervalues():
            for entry2 in entry.children:
                self.assertTrue(entry2 in self.family.groups.entries.itervalues())
            if entry.parent:
                self.assertTrue(entry.parent in self.family.groups.entries.itervalues())
        
        self.assertTrue(self.family.groups.entries['Root'].parent is None)
                
    def test_DRules(self):
        """
        test that there are four rules and each is under a different group
        """
        c = 0
        for rs in self.family.rules.entries.itervalues():
            self.assertLess(len(rs),2,'more than one training reaction at a node')
            if len(rs) == 1:
                c += 1
        
        self.assertEquals(c,4,'incorrect number of kinetics information, expected 4 found {0}'.format(c))
    
    def test_ERegularizationDims(self):
        """
        test that appropriate regularization dimensions have been identified
        """
        
        for entry in self.family.groups.entries.itervalues():
            if entry.children == []:
                continue
            vioObj = set() #set of violations, one atom or one bond is allowed to be in violation (if it was just created)
            pgrp = entry.item
            exts = pgrp.getExtensions()
            for grp,grpc,name,typ,indc in exts:
                if typ == 'intNewBondExt' or typ =='extNewBondExt':
                    continue
                else:
                    val,boo = self.family.evalExt(entry,grp,name)
                    if val != np.inf:
                        continue
                    atms = grp.atoms
                    if typ == 'bondExt':
                        bd = grp.getBond(atms[indc[0]],atms[indc[1]])
                        bds = bd.reg_dim
                        if boo and bds != [] and not (set(bd.order) <= set(bds)):
                            logging.error('bond regularization dimension missed')
                            vioObj.add((tuple(indc),tuple(bds),tuple(bd.order),typ))
                    elif typ == 'atomExt':
                        atypes = atms[indc[0]].reg_dim_atm
                        atype = atms[indc[0]].atomType
                        if boo and atypes != [] and not (set(atype) <= set(atypes)):
                            logging.error('atomtype regularization dimension missed')
                            vioObj.add((tuple(indc),tuple(atypes),tuple(atype),typ))
                    elif typ == 'elExt':
                        us = atms[indc[0]].reg_dim_u
                        u = atms[indc[0]].radicalElectrons
                        if boo and us != [] and not (set(u) <= set(us)):
                            logging.error('unpaired electron regularization dimension missed')
                            vioObj.add((tuple(indc),tuple(us),tuple(u),typ))
            self.assertTrue(len(vioObj) <= 1,'there were {0} regularization violations at, {1}'.format(len(vioObj),vioObj))
    
    def test_FRegularizationStructure(self):
        """
        test that the tree is structured properly after regularization
        """
        self.family.regularize()
        for entry in self.family.groups.entries.itervalues():
            if isinstance(entry.item,Group):
                for child in entry.children:
                    if isinstance(child.item,Group):
                        self.assertTrue(child.item.isSubgraphIsomorphic(entry.item,generateInitialMap=True,saveOrder=True))
    
class TestGenerateReactions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function run ONCE before all unit tests in this class."""
        # Set up a dummy database
        cls.database = RMGDatabase()
        cls.database.load(
            path=os.path.join(settings['test_data.directory'], 'testing_database'),
            thermoLibraries=[],
            reactionLibraries=[],
            kineticsFamilies=['H_Abstraction', 'R_Addition_MultipleBond','Singlet_Val6_to_triplet', 'R_Recombination'],
            depository=False,
            solvation=False,
            testing=True,
        )
        cls.database.loadForbiddenStructures()

    @classmethod
    def tearDownClass(cls):
        """A function run ONCE after all unit tests in this class."""
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    @mock.patch('rmgpy.data.kinetics.family.logging')
    def test_debug_forbidden_reverse_rxn(self, mock_logging):
        """Test that we can automatically debug when a reverse reaction is forbidden."""
        reactants = [Species().fromSMILES('CC'), Species().fromSMILES('[CH2]C=C[CH2]')]
        products = [Species().fromSMILES('C[CH2]'), Species().fromSMILES('[CH2]C=CC')]

        reaction = TemplateReaction(reactants=reactants, products=products)

        successful = self.database.kinetics.families['H_Abstraction'].addReverseAttribute(reaction)

        self.assertFalse(successful)

        mock_logging.error.assert_has_calls([
            mock.call('Expecting one matching reverse reaction, not zero in reaction family H_Abstraction for forward reaction CC + [CH2]C=C[CH2] <=> C[CH2] + [CH2]C=CC.\n'),
        ])

        mock_logging.error.assert_has_calls([
            mock.call('Error was fixed, the product is a forbidden structure when used as a reactant in the reverse direction.'),
        ])

    def test_addAtomLabelsForReaction(self):
        """Test that we can add atom labels to an existing reaction"""
        reactants = [Species().fromSMILES('C=C'), Species().fromSMILES('[OH]')]
        products = [Species().fromSMILES('[CH2]CO')]

        reaction = TemplateReaction(reactants=reactants, products=products)

        self.database.kinetics.families['R_Addition_MultipleBond'].addAtomLabelsForReaction(reaction)

        expected_reactants = [
            Molecule().fromAdjacencyList("""
1 *1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 *2 C u0 p0 c0 {1,D} {5,S} {6,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
"""),
            Molecule().fromAdjacencyList("""
multiplicity 2
1 *3 O u1 p2 c0 {2,S}
2    H u0 p0 c0 {1,S}
""")]

        expected_products = [
            Molecule().fromAdjacencyList("""
multiplicity 2
1 *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 *3 O u0 p2 c0 {1,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
""")]

        for i, reactant in enumerate(reaction.reactants):
            mapping = {}
            for label, atom in expected_reactants[i].getLabeledAtoms().iteritems():
                mapping[atom] = reactant.molecule[0].getLabeledAtom(label)

            self.assertTrue(expected_reactants[i].isIsomorphic(reactant.molecule[0], mapping))

        for i, product in enumerate(reaction.products):
            mapping = {}
            for label, atom in expected_products[i].getLabeledAtoms().iteritems():
                mapping[atom] = product.molecule[0].getLabeledAtom(label)

            self.assertTrue(expected_products[i].isIsomorphic(product.molecule[0], mapping))

    def test_irreversible_reaction(self):
        """Test that the Singlet_Val6_to_triplet and 1,2-Birad_to_alkene families generate irreversible reactions."""

        reactant = [Molecule(SMILES='O=O')]
        reactionList = self.database.kinetics.families['Singlet_Val6_to_triplet'].generateReactions(reactant)
        self.assertFalse(reactionList[0].reversible)

    def test_net_charge_of_products(self):
        """Test that __generateProductStructures() does not generate charged products"""

        reactant = [Molecule(SMILES='[NH-][NH2+]')]
        reactionList = self.database.kinetics.families['R_Recombination'].generateReactions(reactant)
        for rxn in reactionList:
            for product in rxn.products:
                self.assertEquals(product.getNetCharge(), 0)

        reactant = [Molecule(SMILES='[O-][N+]#N')]
        reactionList = self.database.kinetics.families['R_Recombination'].generateReactions(reactant)
        self.assertEquals(len(reactionList), 0)
