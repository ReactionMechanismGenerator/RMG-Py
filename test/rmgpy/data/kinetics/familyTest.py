#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
import logging
import os.path
import shutil

from unittest import mock

import numpy as np

from rmgpy import settings
import rmgpy.data.kinetics.family
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.kinetics import Arrhenius

import pytest


class TestFamily:
    @classmethod
    def setup_class(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        # Set up a dummy database
        cls.database = KineticsDatabase()
        cls.database.load_families(
            path=os.path.join(settings["test_data.directory"], "testing_database/kinetics/families"),
            families=[
                "intra_H_migration",
                "R_Addition_MultipleBond",
                "H_Abstraction",
                "Intra_ene_reaction",
                "6_membered_central_C-C_shift",
                "1,2_shiftC",
                "Intra_R_Add_Exo_scission",
                "intra_substitutionS_isomerization",
                "R_Addition_COm",
                "R_Recombination",
            ],
        )
        cls.family = cls.database.families["intra_H_migration"]

    def test_get_backbone_roots(self):
        """
        Test the get_backbone_roots() function
        """
        backbones = self.family.get_backbone_roots()
        assert backbones[0].label == "RnH"

    def test_get_end_roots(self):
        """
        Test the get_end_roots() function
        """
        ends = self.family.get_end_roots()
        assert len(ends) == 2
        assert self.family.groups.entries["Y_rad_out"] in ends
        assert self.family.groups.entries["XH_out"] in ends

    def test_get_top_level_groups(self):
        """
        Test the get_top_level_groups() function
        """
        top_groups = self.family.get_top_level_groups(self.family.groups.entries["RnH"])
        assert len(top_groups) == 4
        assert self.family.groups.entries["R5Hall"] in top_groups
        assert self.family.groups.entries["R6Hall"] in top_groups
        assert self.family.groups.entries["R2Hall"] in top_groups
        assert self.family.groups.entries["R3Hall"] in top_groups

    def test_react_benzene_bond(self):
        """
        Test that hydrogen addition to benzene (w/ benzene bonds) returns kekulized product.
        """
        family = self.database.families["R_Addition_MultipleBond"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            ),
            Molecule().from_adjacency_list("1 *3 H u1 p0 c0"),
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1
        assert expected_product.is_isomorphic(products[0])

    def test_react_benzene_bond2(self):
        """
        Test that hydrogen addition to phenanthrene (w/ benzene bonds) returns kekulized product.
        """
        family = self.database.families["R_Addition_MultipleBond"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            ),
            Molecule().from_adjacency_list("1 *3 H u1 p0 c0"),
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1
        assert expected_product.is_isomorphic(products[0])

    def test_intra__h_migration(self):
        """
        Test that the intra_H_migration family returns a properly re-labeled product structure.
        This family is its own reverse.
        """
        family = self.database.families["intra_H_migration"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            )
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1

        mapping = {}
        for label, atom in expected_product.get_all_labeled_atoms().items():
            mapping[atom] = products[0].get_labeled_atoms(label)[0]

        assert expected_product.is_isomorphic(products[0], mapping)

    def test_h_abstraction(self):
        """
        Test that the H_Abstraction family returns a properly re-labeled product structure.
        This family is its own reverse.
        """
        family = self.database.families["H_Abstraction"]
        reactants = [
            Molecule().from_adjacency_list(
                """
1 *1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2    C u0 p0 c0 {1,S} {3,D} {7,S}
3    C u0 p0 c0 {2,D} {8,S} {9,S}
4 *2 H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {3,S}
        """
            ),
            Molecule().from_adjacency_list("1 *3 H u1 p0 c0"),
        ]
        expected_products = [
            Molecule().from_adjacency_list(
                """
1 *1 H u0 p0 c0 {2,S}
2 *2 H u0 p0 c0 {1,S}
        """
            ),
            Molecule().from_adjacency_list(
                """
1 *3 C u1 p0 c0 {2,S} {5,S} {6,S}
2    C u0 p0 c0 {1,S} {3,D} {7,S}
3    C u0 p0 c0 {2,D} {8,S} {9,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {3,S}
        """
            ),
        ]
        products = family.apply_recipe(reactants)

        assert len(products) == 2

        mapping1 = {}
        for label, atom in expected_products[0].get_all_labeled_atoms().items():
            mapping1[atom] = products[0].get_labeled_atoms(label)[0]

        assert expected_products[0].is_isomorphic(products[0], mapping1)

        mapping2 = {}
        for label, atom in expected_products[1].get_all_labeled_atoms().items():
            mapping2[atom] = products[1].get_labeled_atoms(label)[0]

        assert expected_products[1].is_isomorphic(products[1], mapping2)

    def test_intra_ene_reaction(self):
        """
        Test that the Intra_ene_reaction family returns a properly re-labeled product structure.
        This family is its own reverse.
        """
        family = self.database.families["Intra_ene_reaction"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            )
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1

        mapping = {}
        for label, atom in expected_product.get_all_labeled_atoms().items():
            mapping[atom] = products[0].get_labeled_atoms(label)[0]

        assert expected_product.is_isomorphic(products[0], mapping)

    def test_6_membered_central_cc_shift(self):
        """
        Test that the 6_membered_central_C-C_shift family returns a properly re-labeled product structure.
        This family is its own reverse.
        """
        family = self.database.families["6_membered_central_C-C_shift"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            )
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1

        mapping = {}
        for label, atom in expected_product.get_all_labeled_atoms().items():
            mapping[atom] = products[0].get_labeled_atoms(label)[0]

        assert expected_product.is_isomorphic(products[0], mapping)

    def test_12_shift_c(self):
        """
        Test that the 1,2_shiftC family returns a properly re-labeled product structure.
        This family is its own reverse.
        """
        family = self.database.families["1,2_shiftC"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            )
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1

        mapping = {}
        for label, atom in expected_product.get_all_labeled_atoms().items():
            mapping[atom] = products[0].get_labeled_atoms(label)[0]

        assert expected_product.is_isomorphic(products[0], mapping)

    def test_intra_r_add_exo_scission(self):
        """
        Test that the Intra_R_Add_Exo_scission family returns a properly re-labeled product structure.
        This family is its own reverse.
        """
        family = self.database.families["Intra_R_Add_Exo_scission"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            )
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1

        mapping = {}
        for label, atom in expected_product.get_all_labeled_atoms().items():
            mapping[atom] = products[0].get_labeled_atoms(label)[0]

        assert expected_product.is_isomorphic(products[0], mapping)

    def test_intra_substitution_s_isomerization(self):
        """
        Test that the intra_substitutionS_isomerization family returns a properly re-labeled product structure.
        This family is its own reverse.
        """
        family = self.database.families["intra_substitutionS_isomerization"]
        reactants = [
            Molecule().from_adjacency_list(
                """
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
"""
            )
        ]
        expected_product = Molecule().from_adjacency_list(
            """
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
"""
        )
        products = family.apply_recipe(reactants)

        assert len(products) == 1

        mapping = {}
        for label, atom in expected_product.get_all_labeled_atoms().items():
            mapping[atom] = products[0].get_labeled_atoms(label)[0]

        assert expected_product.is_isomorphic(products[0], mapping)

    def test_r_addition_com(self):
        """
        Test that the R_Addition_COm family can successfully match the reaction and returns properly product structures.
        This family's product template is generated by charged groups.
        """
        family = self.database.families["R_Addition_COm"]
        reactants = [
            Molecule().from_adjacency_list(
                """
1 *1  C u0 p1 c-1 {2,T}
2 *3  O u0 p1 c+1 {1,T}
"""
            ),
            Molecule().from_adjacency_list(
                """
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
"""
            ),
        ]

        expected_products = [
            Molecule().from_adjacency_list(
                """
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
"""
            ),
        ]

        products = family.apply_recipe(reactants)

        assert len(products) == 1
        assert expected_products[0].is_isomorphic(products[0])

    @pytest.mark.skip(reason="WIP")
    def test_save_family(self):
        """
        This method formerly loaded a family and then wrote it back to disk for comparison
        However it did not account for formatting differences that do not affect the
        actual content of the file.

        It should be rewritten to actually re-load the family after writing it and ensure
        that the objects are the same in memory.
        """
        pass

    def test_reactant_num_id(self):
        """
        Tests that templates aren't applied to the incorrect
        number of reactants
        """
        family = self.database.families["R_Recombination"]
        spc = Molecule().from_smiles("[CH2]CC[CH2]")
        out = family._generate_reactions(reactants=[spc], forward=True)
        assert out == []


class TestTreeGeneration:
    @classmethod
    def setup_class(cls):
        """A function run ONCE before all unit tests in this class."""
        # Set up a dummy database
        cls.database = RMGDatabase()
        cls.database.load(
            path=os.path.join(settings["test_data.directory"], "testing_database"),
            thermo_libraries=[],
            reaction_libraries=[],
            kinetics_families=[],
            depository=False,
            solvation=False,
            surface=False,
            testing=True,
        )
        cls.database.load_forbidden_structures()

        cls.thermoDatabase = ThermoDatabase()  # the real full Thermo Database
        cls.thermoDatabase.load(
            path=os.path.join(settings["database.directory"], "thermo"),
            libraries=["primaryThermoLibrary"],
        )

        cls.kineticsDatabase = KineticsDatabase()
        cls.kineticsDatabase.load_families(
            path=os.path.join(settings["test_data.directory"], "testing_database/kinetics/families"),
            families=[
                "Singlet_Carbene_Intra_Disproportionation",
            ],
        )
        cls.family = cls.kineticsDatabase.families["Singlet_Carbene_Intra_Disproportionation"]
        cls.treerxns = cls.family.get_training_set(
            thermo_database=cls.thermoDatabase,
            remove_degeneracy=True,
            estimate_thermo=True,
            fix_labels=True,
            get_reverse=True,
        )

    @classmethod
    def teardown_class(cls):
        """A function run ONCE after all unit tests in this class."""
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

    def test_a_clear_tree(self):
        """
        Test that the tree was properly cleared before generation
        """
        self.family.clean_tree()
        ents = [ent for ent in self.family.groups.entries.values() if ent.index != -1]
        assert len(ents) == 1, "more than one relevant group left in groups after preparing tree for generation"
        assert len(self.family.rules.entries) == 1, "more than one group in rules.entries after preparing tree for generation"
        root = self.family.groups.entries[list(self.family.rules.entries.keys())[0]]
        assert [root] == self.family.forward_template.reactants
        assert [root] == self.family.groups.top

    def test_b_generate_tree(self):
        """
        test tree generation process
        """

        def objective(k1s, k2s):
            return len(k1s) * np.std(k1s) + len(k2s) * np.std(k2s)

        self.family.generate_tree(thermo_database=self.thermoDatabase, rxns=self.treerxns, obj=objective)  # test input objective function

        self.family.clean_tree()  # reclear

        self.family.generate_tree(thermo_database=self.thermoDatabase, rxns=self.treerxns)  # test that default objective works

    def test_c_parent_child(self):
        """
        test that the tree is structured properly
        """
        for entry in self.family.groups.entries.values():
            for entry2 in entry.children:
                assert entry2 in list(self.family.groups.entries.values())
            if entry.parent:
                assert entry.parent in list(self.family.groups.entries.values())

        assert self.family.groups.entries["Root"].parent is None

    def test_f_rules(self):
        """
        test that there are six rules and each is under a different group
        """
        template_rxn_map = self.family.get_reaction_matches(thermo_database=self.thermoDatabase, remove_degeneracy=True)
        self.family.make_bm_rules_from_template_rxn_map(template_rxn_map)

        c = 0
        for rs in self.family.rules.entries.values():
            assert len(rs) < 2, "more than one training reaction at a node"
            if len(rs) == 1:
                c += 1

        assert c == 6, "incorrect number of kinetics information, expected 6 found {0}".format(c)

    def test_d_regularization_dims(self):
        """
        test that appropriate regularization dimensions have been identified
        """
        template_rxn_map = self.family.get_reaction_matches(thermo_database=self.database.thermo, estimate_thermo=False)

        for entry in self.family.groups.entries.values():
            if entry.children == []:
                continue
            # set of violations, one atom or one bond is allowed to be in violation (if it was just created)
            vio_obj = set()
            pgrp = entry.item
            exts = pgrp.get_extensions()
            for grp, grpc, name, typ, indc in exts:
                if typ == "intNewBondExt" or typ == "extNewBondExt":
                    continue
                else:
                    val, boo = self.family.eval_ext(entry, grp, name, template_rxn_map)
                    if val != np.inf:
                        continue
                    atms = grp.atoms
                    if typ == "bondExt":
                        bd = grp.get_bond(atms[indc[0]], atms[indc[1]])
                        bds = bd.reg_dim[1]
                        if boo and bds != [] and not (set(bd.order) <= set(bds)):
                            logging.error("bond regularization dimension missed")
                            vio_obj.add((tuple(indc), tuple(bds), tuple(bd.order), typ))
                    elif typ == "atomExt":
                        atypes = atms[indc[0]].reg_dim_atm[1]
                        atype = atms[indc[0]].atomtype
                        if boo and atypes != [] and not (set(atype) <= set(atypes)):
                            logging.error("atomtype regularization dimension missed")
                            vio_obj.add((tuple(indc), tuple(atypes), tuple(atype), typ))
                    elif typ == "elExt":
                        us = atms[indc[0]].reg_dim_u[1]
                        u = atms[indc[0]].radical_electrons
                        if boo and us != [] and not (set(u) <= set(us)):
                            logging.error("unpaired electron regularization dimension missed")
                            vio_obj.add((tuple(indc), tuple(us), tuple(u), typ))
                    elif typ == "ringExt":
                        rs = atms[indc[0]].reg_dim_r[1]
                        if "inRing" in atms[indc[0]].props.keys():
                            r = atms[indc[0]].props["inRing"]
                        else:
                            r = [True, False]
                        if boo and rs != [] and not (set(r) <= set(rs)):
                            logging.error("in ring regularization dimension missed")
                            vio_obj.add((tuple(indc), tuple(rs), tuple(r), typ))
                    else:
                        raise ValueError("extension type {0} not identified within test".format(typ))

            assert len(vio_obj) <= 1, "there were {0} regularization violations at, {1}".format(len(vio_obj), vio_obj)

    def test_e_regularization_structure(self):
        """
        test that the tree is structured properly after regularization
        """
        self.family.clean_tree()
        self.family.generate_tree(thermo_database=self.thermoDatabase, rxns=self.treerxns)
        self.family.check_tree()
        self.family.regularize(thermo_database=self.thermoDatabase, rxns=self.treerxns)
        self.family.check_tree()


class TestGenerateReactions:
    @classmethod
    def setup_class(cls):
        """A function run ONCE before all unit tests in this class."""
        # Set up a dummy database
        cls.database = RMGDatabase()
        cls.database.load(
            path=os.path.join(settings["test_data.directory"], "testing_database"),
            thermo_libraries=[],
            reaction_libraries=[],
            kinetics_families=[
                "H_Abstraction",
                "R_Addition_MultipleBond",
                "Singlet_Val6_to_triplet",
                "R_Recombination",
                "Baeyer-Villiger_step1_cat",
                "Surface_Adsorption_Dissociative",
                "Surface_Abstraction_vdW",
                "Surface_Dissociation_vdW",
                "intra_H_migration",
            ],
            depository=False,
            solvation=False,
            surface=False,
            testing=True,
        )
        cls.database.load_forbidden_structures()

    @classmethod
    def teardown_class(cls):
        """A function run ONCE after all unit tests in this class."""
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

    @mock.patch("rmgpy.data.kinetics.family.logging")
    def test_debug_forbidden_reverse_rxn(self, mock_logging):
        """Test that we can automatically debug when a reverse reaction is forbidden."""
        reactants = [
            Species().from_smiles("CC"),
            Species().from_smiles("[CH2]C=C[CH2]"),
        ]
        products = [Species().from_smiles("C[CH2]"), Species().from_smiles("[CH2]C=CC")]

        reaction = TemplateReaction(reactants=reactants, products=products)

        successful = self.database.kinetics.families["H_Abstraction"].add_reverse_attribute(reaction)

        assert not successful

        mock_logging.error.assert_has_calls(
            [
                mock.call(
                    "Expecting one matching reverse reaction, not zero in reaction family H_Abstraction "
                    "for forward reaction CC + [CH2]C=C[CH2] <=> C[CH2] + [CH2]C=CC.\n"
                ),
            ]
        )

        mock_logging.error.assert_has_calls(
            [
                mock.call("Error was fixed, the product is a forbidden structure when " "used as a reactant in the reverse direction."),
            ]
        )

    def test_molecule_forbidden(self):
        forbidden_mol = Molecule(smiles="*CC.[*]")  # vdw bidentate

        mol1 = Molecule(smiles="*CC*")  # bidentate
        mol2 = Molecule(smiles="C.*")  # vdw
        mol3 = Molecule(smiles="CC*")  # chemisorbed

        fam = self.database.kinetics.families["Surface_Dissociation_vdW"]
        assert fam.is_molecule_forbidden(forbidden_mol)
        for allowed_mol in (mol1, mol2, mol3):
            assert not fam.is_molecule_forbidden(allowed_mol)

    def test_add_atom_labels_for_reaction(self):
        """Test that we can add atom labels to an existing reaction"""
        reactants = [Species().from_smiles("C=C"), Species().from_smiles("[OH]")]
        products = [Species().from_smiles("[CH2]CO")]

        reaction = TemplateReaction(reactants=reactants, products=products)

        self.database.kinetics.families["R_Addition_MultipleBond"].add_atom_labels_for_reaction(reaction)

        expected_reactants = [
            Molecule().from_adjacency_list(
                """
1 *1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 *2 C u0 p0 c0 {1,D} {5,S} {6,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
"""
            ),
            Molecule().from_adjacency_list(
                """
multiplicity 2
1 *3 O u1 p2 c0 {2,S}
2    H u0 p0 c0 {1,S}
"""
            ),
        ]

        expected_products = [
            Molecule().from_adjacency_list(
                """
multiplicity 2
1 *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 *3 O u0 p2 c0 {1,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
"""
            )
        ]

        for i, reactant in enumerate(reaction.reactants):
            mapping = {}
            for label, atom in expected_reactants[i].get_all_labeled_atoms().items():
                mapping[atom] = reactant.molecule[0].get_labeled_atoms(label)[0]

            assert expected_reactants[i].is_isomorphic(reactant.molecule[0], mapping)

        for i, product in enumerate(reaction.products):
            mapping = {}
            for label, atom in expected_products[i].get_all_labeled_atoms().items():
                mapping[atom] = product.molecule[0].get_labeled_atoms(label)[0]

            assert expected_products[i].is_isomorphic(product.molecule[0], mapping)

    def test_add_atom_labels_for_reaction_r_recombination(self):
        """Test that we can add atom labels to an existing R_Recombination reaction"""
        reactants = [Species().from_smiles("C[CH2]"), Species().from_smiles("[CH3]")]
        products = [Species().from_smiles("CCC")]

        reaction = TemplateReaction(reactants=reactants, products=products)

        self.database.kinetics.families["R_Recombination"].add_atom_labels_for_reaction(reaction)

        expected_reactants = [
            Molecule().from_adjacency_list(
                """
multiplicity 2
1   C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 * C u1 p0 c0 {1,S} {6,S} {7,S}
3   H u0 p0 c0 {1,S}
4   H u0 p0 c0 {1,S}
5   H u0 p0 c0 {1,S}
6   H u0 p0 c0 {2,S}
7   H u0 p0 c0 {2,S}
"""
            ),
            Molecule().from_adjacency_list(
                """
multiplicity 2
1 * C u1 p0 c0 {2,S} {3,S} {4,S}
2   H u0 p0 c0 {1,S}
3   H u0 p0 c0 {1,S}
4   H u0 p0 c0 {1,S}
"""
            ),
        ]

        expected_products = [
            Molecule().from_adjacency_list(
                """
1  * C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  * C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3    C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {2,S}
9    H u0 p0 c0 {3,S}
10   H u0 p0 c0 {3,S}
11   H u0 p0 c0 {3,S}
"""
            )
        ]

        for i, reactant in enumerate(reaction.reactants):
            mapping = {}
            for label, atom in expected_reactants[i].get_all_labeled_atoms().items():
                mapping[atom] = reactant.molecule[0].get_labeled_atoms(label)[0]

            assert expected_reactants[i].is_isomorphic(reactant.molecule[0], mapping)

        for i, product in enumerate(reaction.products):
            # There are two identical labels in the product, so we need to check both mappings
            # Only one of the mappings will result in isomorphic structures though
            atoms_a = expected_products[i].get_labeled_atoms("*")
            atoms_b = product.molecule[0].get_labeled_atoms("*")
            mapping1 = {atoms_a[0]: atoms_b[0], atoms_a[1]: atoms_b[1]}
            mapping2 = {atoms_a[0]: atoms_b[1], atoms_a[1]: atoms_b[0]}

            results = [
                expected_products[i].is_isomorphic(product.molecule[0], mapping1),
                expected_products[i].is_isomorphic(product.molecule[0], mapping2),
            ]

            assert any(results)
            assert not all(results)

    def test_irreversible_reaction(self):
        """Test that the Singlet_Val6_to_triplet and 1,2-Birad_to_alkene families generate irreversible reactions."""

        reactant = [Molecule(smiles="O=O")]
        reaction_list = self.database.kinetics.families["Singlet_Val6_to_triplet"].generate_reactions(reactant)
        assert not reaction_list[0].reversible

    def test_net_charge_of_products(self):
        """Test that _generate_product_structures() does not generate charged products"""

        reactant = [Molecule(smiles="[NH-][NH2+]")]
        reaction_list = self.database.kinetics.families["R_Recombination"].generate_reactions(reactant)
        for rxn in reaction_list:
            for product in rxn.products:
                assert product.get_net_charge() == 0

        reactant = [Molecule(smiles="[O-][N+]#N")]
        reaction_list = self.database.kinetics.families["R_Recombination"].generate_reactions(reactant)
        assert len(reaction_list) == 0

    def test_reactant_num_mismatch(self):
        """Test that we get no reactions for reactant/template size mismatch

        This happens often because we test every combo of molecules against all families.
        """
        reactants = [Molecule(smiles="C"), Molecule(smiles="[OH]")]
        reaction_list = self.database.kinetics.families["Singlet_Val6_to_triplet"].generate_reactions(reactants)
        assert len(reaction_list) == 0
        reaction_list = self.database.kinetics.families["Baeyer-Villiger_step1_cat"].generate_reactions(reactants)
        assert len(reaction_list) == 0
        reaction_list = self.database.kinetics.families["Surface_Adsorption_Dissociative"].generate_reactions(reactants)
        assert len(reaction_list) == 0

    def test_match_reactant_to_template_surface_site(self):
        """
        Test that an empty surface site template group matches an empty surface site Molecule and does not match
        a vdW adsorbate
        """
        family = self.database.kinetics.families["Surface_Adsorption_Dissociative"]
        empty_surface_site_template_group = [r.item for r in family.forward_template.reactants if r.item.is_surface_site()][0]

        empty_surface_site_mol = Molecule().from_adjacency_list("1 X u0")
        vdW_adsorbate = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 X u0 p0 c0
"""
        )
        empty_surface_site_matches = family._match_reactant_to_template(empty_surface_site_mol, empty_surface_site_template_group)
        vdW_matches = family._match_reactant_to_template(vdW_adsorbate, empty_surface_site_template_group)
        assert len(empty_surface_site_matches) == 1
        assert len(vdW_matches) == 0

    def test_reactant_num_mismatch_2(self):
        """Test that we get no reactions for reactant/template size mismatch

        This happens often because we test every combo of molecules against all families.
        """
        reactants = [
            Molecule().from_smiles("CC"),
            Molecule().from_adjacency_list("1 X u0"),
            Molecule().from_adjacency_list("1 X u0"),
        ]
        reaction_list = self.database.kinetics.families["Surface_Dissociation_vdW"].generate_reactions(reactants)
        assert len(reaction_list) == 0

    def test_apply_recipe_multiplicity_check(self):
        """
        Test that the multiplicity check is working correctly in the apply_recipe function
        """
        family = self.database.kinetics.families["Surface_Abstraction_vdW"]
        reacts = [Molecule(smiles="*[CH2]"), Molecule(smiles="*[CH2]")]
        reaction_list = family.generate_reactions(reacts)
        assert len(reaction_list) == 0

    def test_retaining_atom_labels_in_template_reaction(self):
        """
        Test that atom labels are not deleted from a TemplateReaction if so requested.
        """
        family = self.database.kinetics.families["intra_H_migration"]
        reacts = [Molecule(smiles="C[CH]C")]
        reaction_list_1 = family.generate_reactions(reacts)
        assert not hasattr(reaction_list_1[0], "labeled_atoms")
        reaction_list_2 = family.generate_reactions(reacts, delete_labels=False, relabel_atoms=False)
        assert hasattr(reaction_list_2[0], "labeled_atoms")
        assert [(label, str(atom)) for label, atom in reaction_list_2[0].labeled_atoms["reactants"].items()] == [
            ("*2", "C"),
            ("*1", "C."),
            ("*3", "H"),
        ]
        assert [(label, str(atom)) for label, atom in reaction_list_2[0].labeled_atoms["products"].items()] == [
            ("*1", "C"),
            ("*2", "C."),
            ("*3", "H"),
        ]

    def test_average_kinetics(self):
        """
        Test that the average kinetics are calculated correctly
        """
        k1 = Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"), T0=(1, "K"))
        k2 = Arrhenius(A=(4e13, "cm^3/(mol*s)"), n=1, Ea=(10, "kJ/mol"), T0=(1, "K"))
        kav = rmgpy.data.kinetics.family.average_kinetics([k1, k2])
        assert round(abs(kav.A.value_si - 2.0e7), 2) == 0  # m3/mol/s
        assert round(abs(kav.n.value_si - 0.5), 6) == 0
        assert round(abs(kav.Ea.value_si - 5.0e3), 2) == 0
        assert round(abs(kav.T0.value_si - 1.0), 6) == 0
        assert kav.A.units == "m^3/(mol*s)"
        assert (
            round(
                abs(
                    np.log(kav.get_rate_coefficient(300))
                    - np.average(
                        [
                            np.log(k1.get_rate_coefficient(300)),
                            np.log(k2.get_rate_coefficient(300)),
                        ]
                    )
                ),
                6,
            )
            == 0
        )
