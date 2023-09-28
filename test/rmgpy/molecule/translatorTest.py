#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This module contains unit test for the translator module.
"""

import re
from unittest.mock import patch

from rmgpy.molecule.adjlist import ConsistencyChecker
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.inchi import (
    compose_aug_inchi,
    P_LAYER_PREFIX,
    P_LAYER_SEPARATOR,
    U_LAYER_PREFIX,
    U_LAYER_SEPARATOR,
)
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.translator import *
from rmgpy.species import Species
import pytest


class TranslatorTest:
    def test_empty_molecule(self):
        """Test that we can safely return a blank identifier for an empty molecule."""
        mol = Molecule()

        assert mol.to_smiles() == ""
        assert mol.to_inchi() == ""

    @patch("rmgpy.molecule.translator.logging")
    def test_failure_message(self, mock_logging):
        """Test that we log the molecule adjlist upon failure."""
        mol = Molecule(smiles="[CH2-][N+]#N")

        with pytest.raises(ValueError, match="Unable to generate identifier type"):
            to_inchi(mol, backend="rdkit")

        mock_logging.error.assert_called_with("Unable to generate identifier for this molecule:\n{0}".format(mol.to_adjacency_list()))


class InChIGenerationTest:
    def compare(self, adjlist, aug_inchi):
        spc = Species(molecule=[Molecule().from_adjacency_list(adjlist)])
        spc.generate_resonance_structures()

        ignore_prefix = r"(InChI=1+)(S*)/"

        exp = re.split(ignore_prefix, aug_inchi)[-1]
        comp = re.split(ignore_prefix, spc.get_augmented_inchi())[-1]
        assert exp == comp

    def test_c5h5(self):
        """
        Test that the unpaired electron of 1,3-cyclopentadienyl radical always
        ends up on the 1-carbon atom.
        """

        adjlist = """
1 C 0 {2,D} {5,S}
2 C 0 {1,D} {3,S} 
3 C 0 {2,S} {4,D} 
4 C 0 {3,D} {5,S} 
5 C 1 {4,S} {1,S}
        """

        aug_inchi = "InChI=1S/C5H5/c1-2-4-5-3-1/h1-5H/u1"
        self.compare(adjlist, aug_inchi)

    def test_c7h8(self):
        """Looks a lot like toluene but with 1 double bond replaced by a biradical.

        unpaired electrons on tertiary carbon, and on carbon in para position."""
        adjlist = """
1  C u1 p0 c0 {2,S} {3,S} {4,S}
2  C u0 p0 c0 {1,S} {7,D} {10,S}
3  C u0 p0 c0 {1,S} {6,D} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {6,S} {7,S} {15,S}
6  C u0 p0 c0 {3,D} {5,S} {8,S}
7  C u0 p0 c0 {2,D} {5,S} {9,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
        """

        aug_inchi = "InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3/u2,3"
        self.compare(adjlist, aug_inchi)

    def test_c8h8(self):
        """Looks a lot like cycloctene but with 1 double bond replaced by a biradical."""

        adjlist = """
1  C u0 p0 c0 {2,S} {5,D} {9,S}
2  C u0 p0 c0 {1,S} {3,D} {10,S}
3  C u0 p0 c0 {2,D} {4,S} {11,S}
4  C u1 p0 c0 {3,S} {6,S} {12,S}
5  C u0 p0 c0 {1,D} {8,S} {14,S}
6  C u1 p0 c0 {4,S} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {13,S}
8  C u0 p0 c0 {5,S} {7,D} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
        """

        aug_inchi = "InChI=1S/C8H8/c1-2-4-6-8-7-5-3-1/h1-8H/u1,2"
        self.compare(adjlist, aug_inchi)

    def test_benzyne(self):
        adjlist = """
1  C u0 p0 c0 {2,T} {6,S}
2  C u0 p0 c0 {1,T} {3,S}
3  C u0 p0 c0 {2,S} {4,D} {7,S}
4  C u0 p0 c0 {3,D} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,D} {9,S}
6  C u0 p0 c0 {1,S} {5,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
        """
        aug_inchi = "InChI=1S/C6H4/c1-2-4-6-5-3-1/h1-4H"
        self.compare(adjlist, aug_inchi)

    def test_h(self):
        adjlist = """
multiplicity 2
1 H u1 p0 c0
"""
        aug_inchi = "InChI=1S/H/u1"
        self.compare(adjlist, aug_inchi)

    def test_c6h8(self):
        """
        Test that the 2 unpaired electrons of .CC(=C)C(C.)=C
        do not end up at the same side of the central C-C bond.
        """
        adjlist = """
1 C 0 {2,D}
2 C 0 {1,D} {3,S} {4,S}
3 C 1 {2,S}
4 C 0 {2,S} {5,S} {6,D}
5 C 1 {4,S}
6 C 0 {4,D}
        """

        aug_inchi = "InChI=1S/C6H8/c1-5(2)6(3)4/h1-4H2/u1,3"
        self.compare(adjlist, aug_inchi)

    def test_c6h10_tetrarad(self):
        adjlist = """
1  C u1 p0 c0 {3,S} {7,S} {8,S}
2  C u1 p0 c0 {4,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {5,S} {11,S}
4  C u1 p0 c0 {2,S} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {5,S} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
        """

        aug_inchi = "InChI=1S/C6H10/c1-3-5-6-4-2/h3-4H,1-2,5-6H2/u1,2,3,4"
        self.compare(adjlist, aug_inchi)

    def test_buta13diyl_triplet(self):
        """
        C=CC.C.
        """
        adjlist = """
        multiplicity 3
1  C u1 p0 c0 {2,S} {5,S} {6,S}
2  C u1 p0 c0 {1,S} {3,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""

        aug_inchi = "InChI=1S/C4H6/c1-3-4-2/h3-4H,1-2H2/u1,2"
        self.compare(adjlist, aug_inchi)

    def test_ch2o2(self):
        adjlist = """
1 C 1 {2,S} {3,S}
2 O 0 {1,S}
3 O 1 {1,S}
"""

        aug_inchi = "InChI=1/CH2O2/c2-1-3/h1H,(H,2,3)/u1,2"
        self.compare(adjlist, aug_inchi)

    def test_c7h10(self):
        adjlist = """

        1 C 1 {2,S}
2 C 0 {1,S} {3,D} {4,S}
3 C 0 {2,D}
4 C 0 {2,S} {5,S}
5 C 1 {4,S} {6,S} {7,S}
6 C 1 {5,S}
7 C 1 {5,S}
"""

        aug_inchi = "InChI=1S/C7H10/c1-6(2)5-7(3)4/h1-5H2/u1,2,3,6"
        self.compare(adjlist, aug_inchi)

    def test_c5h6o(self):
        adjlist = """
1 C 1 {2,S}
2 C 0 {1,S} {3,D}
3 C 0 {2,D} {4,S} {5,S}
4 O 1 {3,S}
5 C 0 {3,S} {6,D}
6 C 0 {5,D}
"""

        aug_inchi = "InChI=1S/C5H6O/c1-3-5(6)4-2/h3-4H,1-2H2/u1,3"
        self.compare(adjlist, aug_inchi)

    def test_c7h9(self):
        adjlist = """
1 C 0 {4,D} 
2 C 0 {5,D}
3 C 1 {6,S}
4 C 0 {1,D} {7,S}
5 C 0 {2,D} {7,S}
6 C 1 {3,S} {7,S}
7 C 1 {4,S} {5,S} {6,S}
"""

        aug_inchi = "InChI=1S/C7H9/c1-4-7(5-2)6-3/h4-6H,1-3H2/u1,2,4"
        self.compare(adjlist, aug_inchi)

    def test_c11h16(self):
        adjlist = """
1 C 0 {5,D}
2 C 1 {6,S}
3 C 1 {7,S}
4 C 0 {8,D}
5 C 0 {1,D} {9,S}
6 C 1 {2,S} {10,S}
7 C 1 {3,S} {11,S}
8 C 0 {4,D} {11,S}
9 C 0 {5,S} {11,S}
10 C 0 {6,S} {11,S}
11 C 0 {7,S} {8,S} {9,S} {10,S}
"""

        aug_inchi = "InChI=1S/C11H16/c1-5-9-11(7-3,8-4)10-6-2/h5-8H,1-4,9-10H2/u1,3,5,7"
        self.compare(adjlist, aug_inchi)

    def test_singlet_vs_triplet(self):
        adjlist_singlet = """
        1 C u0 p1 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """

        adjlist_triplet = """
        multiplicity 3
        1 C u2 p0 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """

        singlet = Species(molecule=[Molecule().from_adjacency_list(adjlist_singlet)])
        triplet = Species(molecule=[Molecule().from_adjacency_list(adjlist_triplet)])
        singlet_aug_inchi = singlet.get_augmented_inchi()
        triplet_aug_inchi = triplet.get_augmented_inchi()
        assert singlet_aug_inchi != triplet_aug_inchi

    #     def test_c6h5(self):
    #         """Test that the u-layer of phenyl shows atom 1."""
    #         adjlist = """
    # multiplicity 2
    # 1  C u0 p0 c0 {2,D} {3,S} {10,S}
    # 2  C u0 p0 c0 {1,D} {5,S} {7,S}
    # 3  C u0 p0 c0 {1,S} {6,D} {8,S}
    # 4  C u0 p0 c0 {5,D} {6,S} {11,S}
    # 5  C u0 p0 c0 {2,S} {4,D} {9,S}
    # 6  C u1 p0 c0 {3,D} {4,S}
    # 7  H u0 p0 c0 {2,S}
    # 8  H u0 p0 c0 {3,S}
    # 9  H u0 p0 c0 {5,S}
    # 10 H u0 p0 c0 {1,S}
    # 11 H u0 p0 c0 {4,S}
    # """

    #         aug_inchi = 'InChI=1S/C6H5/c1-2-4-6-5-3-1/h1-5H/u1'
    #         self.compare(adjlist, aug_inchi)

    def test_c5h6_singlet(self):
        """
        n-C5 chain with 1 lone pair at the central carbon atom
        """
        adjlist = """
        1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
        2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
        3  C u0 p1 c0 {2,S} {4,S}
        4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
        5  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
        6  H u0 p0 c0 {1,S}
        7  H u0 p0 c0 {1,S}
        8  H u0 p0 c0 {1,S}
        9  H u0 p0 c0 {2,S}
        10 H u0 p0 c0 {2,S}
        11 H u0 p0 c0 {4,S}
        12 H u0 p0 c0 {4,S}
        13 H u0 p0 c0 {5,S}
        14 H u0 p0 c0 {5,S}
        15 H u0 p0 c0 {5,S}
        """
        aug_inchi = "C5H10/c1-3-5-4-2/h3-4H2,1-2H3/lp5"
        self.compare(adjlist, aug_inchi)

    def test_aromatic_resonance_structures(self):
        """Test that different resonance structures give identical InChIs."""
        mol = Molecule().from_adjacency_list(
            """
multiplicity 2
1  C u0 p0 c0 {2,D} {14,S} {18,S}
2  C u0 p0 c0 {1,D} {3,S} {19,S}
3  C u0 p0 c0 {2,S} {4,D} {20,S}
4  C u0 p0 c0 {3,D} {5,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,D}
6  C u0 p0 c0 {5,S} {7,D} {21,S}
7  C u0 p0 c0 {6,D} {8,S} {22,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {8,D} {10,S} {23,S}
10 C u0 p0 c0 {9,S} {11,D} {24,S}
11 C u0 p0 c0 {10,D} {12,S} {25,S}
12 C u0 p0 c0 {11,S} {13,D} {26,S}
13 C u0 p0 c0 {4,S} {8,S} {12,D}
14 C u0 p0 c0 {1,S} {5,D} {15,S}
15 C u1 p0 c0 {14,S} {16,S} {17,S}
16 H u0 p0 c0 {15,S}
17 H u0 p0 c0 {15,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {12,S}
"""
        )
        res = mol.generate_resonance_structures()

        inchi_list = [struct.to_inchi() for struct in res]

        expected_inchi = "InChI=1S/C15H11/c1-11-5-4-8-15-13(11)10-9-12-6-2-3-7-14(12)15/h2-10H,1H2"

        for inchi in inchi_list:
            assert inchi == expected_inchi

    def test_disconnected_molecule(self):
        """Test that we can generate an InChI for a disconnected molecule."""
        mol = Molecule().from_smiles("CCCCO.C=O")

        inchi = "InChI=1S/C4H10O.CH2O/c1-2-3-4-5;1-2/h5H,2-4H2,1H3;1H2"

        assert mol.to_inchi() == inchi

    def test_isotopic_molecule_1(self):
        """Test that we can generate an InChI for an isotopic molecule."""
        mol = Molecule().from_smiles("[13CH4]")

        inchi = "InChI=1S/CH4/h1H4/i1+1"

        assert mol.to_inchi() == inchi

    def test_isotopic_molecule_2(self):
        """Test that we can generate an InChI for an isotopic molecule."""
        mol = Molecule().from_smiles("[13CH3]C")

        inchi = "InChI=1S/C2H6/c1-2/h1-2H3/i1+1"

        assert mol.to_inchi() == inchi

    def test_surface_molecule_rdkit(self):
        """Test InChI generation for a surface molecule using RDKit"""
        mol = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}
"""
        )
        inchi = "InChI=1S/CH3.Pt/h1H3;"

        assert to_inchi(mol, backend="rdkit") == inchi

    def test_surface_molecule_ob(self):
        """Test InChI generation for a surface molecule using OpenBabel"""
        mol = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}
"""
        )
        inchi = "InChI=1S/CH3.Pt/h1H3;"

        assert to_inchi(mol, backend="openbabel") == inchi


class SMILESGenerationTest:
    def compare(self, adjlist, smiles):
        mol = Molecule().from_adjacency_list(adjlist)
        assert smiles == mol.to_smiles()

    def test_ch4(self):
        "Test the SMILES generation for methane"

        adjlist = """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
        """
        smiles = "C"
        self.compare(adjlist, smiles)

    def test_c(self):
        "Test the SMILES generation for atomic carbon mult=(1,3,5)"
        adjlist = "1 C u0 p2 c0"
        smiles = "[C]"
        self.compare(adjlist, smiles)

        adjlist = "multiplicity 3\n1 C u2 p1 c0"
        smiles = "[C]"
        self.compare(adjlist, smiles)

        adjlist = "multiplicity 5\n1 C u4 p0 c0"
        smiles = "[C]"
        self.compare(adjlist, smiles)

    def test_various(self):
        "Test the SMILES generation for various molecules and radicals"

        # Test N2
        adjlist = """
        1 N u0 p1 c0 {2,T}
        2 N u0 p1 c0 {1,T}
        """
        smiles = "N#N"
        self.compare(adjlist, smiles)

        # Test CH4
        adjlist = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        """
        smiles = "C"
        self.compare(adjlist, smiles)

        # Test H2O
        adjlist = """
        1 O u0 p2 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """
        smiles = "O"
        self.compare(adjlist, smiles)

        # Test C2H6
        adjlist = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """
        smiles = "CC"
        self.compare(adjlist, smiles)

        # Test H2
        adjlist = """
        1 H u0 p0 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        """
        smiles = "[H][H]"
        self.compare(adjlist, smiles)

        # Test H2O2
        adjlist = """
        1 O u0 p2 c0 {2,S} {3,S}
        2 O u0 p2 c0 {1,S} {4,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {2,S}
        """
        smiles = "OO"
        self.compare(adjlist, smiles)

        # Test C3H8
        adjlist = """
        1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
        2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
        3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
        4  H u0 p0 c0 {1,S}
        5  H u0 p0 c0 {1,S}
        6  H u0 p0 c0 {1,S}
        7  H u0 p0 c0 {2,S}
        8  H u0 p0 c0 {2,S}
        9  H u0 p0 c0 {3,S}
        10 H u0 p0 c0 {3,S}
        11 H u0 p0 c0 {3,S}
        """
        smiles = "CCC"
        self.compare(adjlist, smiles)

        # Test Ar
        adjlist = """
        1 Ar u0 p4 c0
        """
        smiles = "[Ar]"
        self.compare(adjlist, smiles)

        # Test He
        adjlist = """
        1 He u0 p1 c0
        """
        smiles = "[He]"
        self.compare(adjlist, smiles)

        # Test CH4O
        adjlist = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 O u0 p2 c0 {1,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        """
        smiles = "CO"
        self.compare(adjlist, smiles)

        # Test CO2
        adjlist = """
        1 O u0 p2 c0 {2,D}
        2 C u0 p0 c0 {1,D} {3,D}
        3 O u0 p2 c0 {2,D}
        """
        smiles = "O=C=O"
        self.compare(adjlist, smiles)

        # Test CO
        adjlist = """
        1 C u0 p1 c-1 {2,T}
        2 O u0 p1 c+1 {1,T}
        """
        smiles = "[C-]#[O+]"
        self.compare(adjlist, smiles)

        # Test C2H4
        adjlist = """
        1 C u0 p0 c0 {2,D} {3,S} {4,S}
        2 C u0 p0 c0 {1,D} {5,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {2,S}
        6 H u0 p0 c0 {2,S}
        """
        smiles = "C=C"
        self.compare(adjlist, smiles)

        # Test O2
        adjlist = """
        1 O u0 p2 c0 {2,D}
        2 O u0 p2 c0 {1,D}
        """
        smiles = "O=O"
        self.compare(adjlist, smiles)

        # Test CH3
        adjlist = """
        multiplicity 2
        1 C u1 p0 c0 {2,S} {3,S} {4,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        """
        smiles = "[CH3]"
        self.compare(adjlist, smiles)

        # Test HO
        adjlist = """
        multiplicity 2
        1 O u1 p2 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        """
        smiles = "[OH]"
        self.compare(adjlist, smiles)

        # Test C2H5
        adjlist = """
        multiplicity 2
        1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
        2 C u1 p0 c0 {1,S} {3,S} {4,S}
        3 H u0 p0 c0 {2,S}
        4 H u0 p0 c0 {2,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {1,S}
        7 H u0 p0 c0 {1,S}
        """
        smiles = "C[CH2]"
        self.compare(adjlist, smiles)

        # Test O
        adjlist = """
        multiplicity 3
        1 O u2 p2 c0
        """
        smiles = "[O]"
        self.compare(adjlist, smiles)

        # Test HO2
        adjlist = """
        multiplicity 2
        1 O u1 p2 c0 {2,S}
        2 O u0 p2 c0 {1,S} {3,S}
        3 H u0 p0 c0 {2,S}
        """
        smiles = "[O]O"
        self.compare(adjlist, smiles)

        # Test CH
        adjlist = """
        multiplicity 4
        1 C u3 p0 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        """
        smiles = "[CH]"
        self.compare(adjlist, smiles)

        # Test H
        adjlist = """
        multiplicity 2
        1 H u1 p0 c0
        """
        smiles = "[H]"
        self.compare(adjlist, smiles)

        # Test C
        adjlist = """
        multiplicity 5
        1 C u4 p0 c0
        """
        smiles = "[C]"
        self.compare(adjlist, smiles)

        # Test O2
        adjlist = """
        multiplicity 3
        1 O u1 p2 c0 {2,S}
        2 O u1 p2 c0 {1,S}
        """
        smiles = "[O][O]"
        self.compare(adjlist, smiles)

        # Test CF
        adjlist = """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 F u0 p3 c0 {1,S}
        """
        smiles = "[C]F"
        self.compare(adjlist, smiles)

        # Test CCl
        adjlist = """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 Cl u0 p3 c0 {1,S}
        """
        smiles = "[C]Cl"
        self.compare(adjlist, smiles)

        # Test CBr
        adjlist = """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 Br u0 p3 c0 {1,S}
        """
        smiles = "[C]Br"
        self.compare(adjlist, smiles)

        # Test CC-triplet-1
        adjlist = '''
        multiplicity 3
        1 *3 C u1 p1 c0 {2,S}
        2 *3 C u1 p1 c0 {1,S}
        '''
        smiles = '[C][C]'
        self.compare(adjlist, smiles)

        # Test CC-triplet-2
        adjlist = '''
        multiplicity 3
        1 *3 C u1 p0 c0 {2,T}
        2 *3 C u1 p0 c0 {1,T}
        '''
        smiles = '[C]#[C]'
        self.compare(adjlist, smiles)

        # todo: Test CC-singlet-1
        # We couldn't test the case where C forms quadruple bond with C
        # because `$ `is not added to rdkit until Jan 2022,
        # and RMG environment usually has rdkit <= 2022.03
        # as of Sep 2023.
        # This test should be added after we update the rdkit dependency.

        # adjlist = '''
        # 1 *3 C u0 p0 c0 {2,Q}
        # 2 *3 C u0 p0 c0 {1,Q}
        # '''
        # smiles = '[C]$[C]'
        # self.compare(adjlist, smiles)

        # Test CC-singlet-2
        # We couldn't test the case where C forms quadruple bond with C
        adjlist = '''
        1 *3 C u0 p1 c0 {2,D}
        2 *3 C u0 p1 c0 {1,D}
        '''
        smiles = '[C]=[C]'
        self.compare(adjlist, smiles)

    def test_aromatics(self):
        """Test that different aromatics representations returns different SMILES."""
        mol1 = Molecule().from_adjacency_list(
            """
1  O u0 p2 c0 {6,S} {9,S}
2  C u0 p0 c0 {3,D} {5,S} {11,S}
3  C u0 p0 c0 {2,D} {4,S} {12,S}
4  C u0 p0 c0 {3,S} {6,D} {13,S}
5  C u0 p0 c0 {2,S} {7,D} {10,S}
6  C u0 p0 c0 {1,S} {4,D} {7,S}
7  C u0 p0 c0 {5,D} {6,S} {8,S}
8  C u0 p0 c0 {7,S} {14,S} {15,S} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
"""
        )
        mol2 = Molecule().from_adjacency_list(
            """
1  O u0 p2 c0 {6,S} {9,S}
2  C u0 p0 c0 {3,S} {5,D} {11,S}
3  C u0 p0 c0 {2,S} {4,D} {12,S}
4  C u0 p0 c0 {3,D} {6,S} {13,S}
5  C u0 p0 c0 {2,D} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {8,S}
8  C u0 p0 c0 {7,S} {14,S} {15,S} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
"""
        )
        mol3 = Molecule().from_adjacency_list(
            """
1  O u0 p2 c0 {6,S} {9,S}
2  C u0 p0 c0 {3,B} {5,B} {11,S}
3  C u0 p0 c0 {2,B} {4,B} {12,S}
4  C u0 p0 c0 {3,B} {6,B} {13,S}
5  C u0 p0 c0 {2,B} {7,B} {10,S}
6  C u0 p0 c0 {1,S} {4,B} {7,B}
7  C u0 p0 c0 {5,B} {6,B} {8,S}
8  C u0 p0 c0 {7,S} {14,S} {15,S} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
"""
        )

        smiles1 = mol1.to_smiles()
        smiles2 = mol2.to_smiles()
        smiles3 = mol3.to_smiles()

        assert smiles1 != smiles2
        assert smiles2 != smiles3
        assert smiles1 != smiles3

    def test_surface_molecule_rdkit(self):
        """Test InChI generation for a surface molecule using RDKit"""
        mol = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}
"""
        )
        smiles = "C[Pt]"

        assert to_smiles(mol, backend="rdkit") == smiles

    def test_surface_molecule_ob(self):
        """Test InChI generation for a surface molecule using OpenBabel"""
        mol = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}
"""
        )
        smiles = "C[Pt]"

        assert to_smiles(mol, backend="openbabel") == smiles


class ParsingTest:
    def setup_class(self):
        self.methane = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""
        )
        self.methylamine = Molecule().from_adjacency_list(
            """
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""
        )

    def test_from_augmented_inchi(self):
        aug_inchi = "InChI=1S/CH4/h1H4"
        mol = from_augmented_inchi(Molecule(), aug_inchi)
        assert not mol.inchi == ""
        assert mol.is_isomorphic(self.methane)

        aug_inchi = "InChI=1/CH4/h1H4"
        mol = from_augmented_inchi(Molecule(), aug_inchi)
        assert not mol.inchi == ""
        assert mol.is_isomorphic(self.methane)

    def compare(self, adjlist, smiles):
        """
        Compare result of parsing an adjacency list and a SMILES string.

        The adjacency list is presumed correct and this is to test the SMILES parser.
        """
        mol1 = Molecule().from_adjacency_list(adjlist)
        mol2 = Molecule(smiles=smiles)
        assert mol1.is_isomorphic(mol2), "Parsing SMILES={!r} gave unexpected molecule\n{}".format(smiles, mol2.to_adjacency_list())

    def test_from_smiles(self):
        smiles = "C"
        mol = from_smiles(Molecule(), smiles)
        assert mol.is_isomorphic(self.methane)

        # Test that atomtypes that rely on lone pairs for identity are typed correctly
        smiles = "CN"
        mol = from_smiles(Molecule(), smiles)
        assert mol.atoms[1].atomtype == ATOMTYPES["N3s"]

        # Test N2
        adjlist = """
        1 N u0 p1 c0 {2,T}
        2 N u0 p1 c0 {1,T}
        """
        smiles = "N#N"
        self.compare(adjlist, smiles)

        # Test CH4
        adjlist = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        """
        smiles = "C"
        self.compare(adjlist, smiles)

        # Test H2O
        adjlist = """
        1 O u0 p2 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """
        smiles = "O"
        self.compare(adjlist, smiles)

        # Test C2H6
        adjlist = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """
        smiles = "CC"
        self.compare(adjlist, smiles)

        # Test H2
        adjlist = """
        1 H u0 p0 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        """
        smiles = "[H][H]"
        self.compare(adjlist, smiles)

        # Test H2O2
        adjlist = """
        1 O u0 p2 c0 {2,S} {3,S}
        2 O u0 p2 c0 {1,S} {4,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {2,S}
        """
        smiles = "OO"
        self.compare(adjlist, smiles)

        # Test C3H8
        adjlist = """
        1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
        2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
        3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
        4  H u0 p0 c0 {1,S}
        5  H u0 p0 c0 {1,S}
        6  H u0 p0 c0 {1,S}
        7  H u0 p0 c0 {2,S}
        8  H u0 p0 c0 {2,S}
        9  H u0 p0 c0 {3,S}
        10 H u0 p0 c0 {3,S}
        11 H u0 p0 c0 {3,S}
        """
        smiles = "CCC"
        self.compare(adjlist, smiles)

        # Test Ar
        adjlist = """
        1 Ar u0 p4 c0
        """
        smiles = "[Ar]"
        self.compare(adjlist, smiles)

        # Test He
        adjlist = """
        1 He u0 p1 c0
        """
        smiles = "[He]"
        self.compare(adjlist, smiles)

        # Test CH4O
        adjlist = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 O u0 p2 c0 {1,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        """
        smiles = "CO"
        self.compare(adjlist, smiles)

        # Test CO2
        adjlist = """
        1 O u0 p2 c0 {2,D}
        2 C u0 p0 c0 {1,D} {3,D}
        3 O u0 p2 c0 {2,D}
        """
        smiles = "O=C=O"
        self.compare(adjlist, smiles)

        # Test CO
        adjlist = """
        1 C u0 p1 c-1 {2,T}
        2 O u0 p1 c+1 {1,T}
        """
        smiles = "[C-]#[O+]"
        self.compare(adjlist, smiles)

        # Test C2H4
        adjlist = """
        1 C u0 p0 c0 {2,D} {3,S} {4,S}
        2 C u0 p0 c0 {1,D} {5,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {2,S}
        6 H u0 p0 c0 {2,S}
        """
        smiles = "C=C"
        self.compare(adjlist, smiles)

        # Test O2
        adjlist = """
        1 O u0 p2 c0 {2,D}
        2 O u0 p2 c0 {1,D}
        """
        smiles = "O=O"
        self.compare(adjlist, smiles)

        # Test CH3
        adjlist = """
        multiplicity 2
        1 C u1 p0 c0 {2,S} {3,S} {4,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        """
        smiles = "[CH3]"
        self.compare(adjlist, smiles)

        # Test HO
        adjlist = """
        multiplicity 2
        1 O u1 p2 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        """
        smiles = "[OH]"
        self.compare(adjlist, smiles)

        # Test C2H5
        adjlist = """
        multiplicity 2
        1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
        2 C u1 p0 c0 {1,S} {3,S} {4,S}
        3 H u0 p0 c0 {2,S}
        4 H u0 p0 c0 {2,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {1,S}
        7 H u0 p0 c0 {1,S}
        """
        smiles = "C[CH2]"
        self.compare(adjlist, smiles)

        # Test O
        adjlist = """
        multiplicity 3
        1 O u2 p2 c0
        """
        smiles = "[O]"
        self.compare(adjlist, smiles)

        # Test HO2
        adjlist = """
        multiplicity 2
        1 O u1 p2 c0 {2,S}
        2 O u0 p2 c0 {1,S} {3,S}
        3 H u0 p0 c0 {2,S}
        """
        smiles = "[O]O"
        self.compare(adjlist, smiles)

        # Test CH, methylidyne.
        # Wikipedia reports:
        # The ground state is a doublet radical with one unpaired electron,
        # and the first two excited states are a quartet radical with three
        # unpaired electrons and a doublet radical with one unpaired electron.
        # With the quartet radical only 71 kJ above the ground state, a sample
        # of methylidyne exists as a mixture of electronic states even at
        # room temperature, giving rise to complex reactions.
        #
        adjlist = """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        """
        smiles = "[CH]"
        self.compare(adjlist, smiles)

        # Test H
        adjlist = """
        multiplicity 2
        1 H u1 p0 c0
        """
        smiles = "[H]"
        self.compare(adjlist, smiles)

        # Test atomic C, which is triplet in ground state
        adjlist = """
        multiplicity 3
        1 C u2 p1 c0
        """
        smiles = "[C]"
        self.compare(adjlist, smiles)

        # Test O2
        adjlist = """
        multiplicity 3
        1 O u1 p2 c0 {2,S}
        2 O u1 p2 c0 {1,S}
        """
        smiles = "[O][O]"
        self.compare(adjlist, smiles)

    def test_from_inchi(self):
        inchi = "InChI=1S/CH4/h1H4"
        mol = from_inchi(Molecule(), inchi)
        assert mol.is_isomorphic(self.methane)
        # Test that atomtypes that rely on lone pairs for identity are typed correctly
        inchi = "InChI=1S/CH5N/c1-2/h2H2,1H3"
        mol = from_inchi(Molecule(), inchi)
        assert mol.atoms[1].atomtype == ATOMTYPES["N3s"]

    # current implementation of SMARTS is broken
    def test_from_smarts(self):
        smarts = "[CH4]"
        mol = from_smarts(Molecule(), smarts)
        assert mol.is_isomorphic(self.methane)

    def test_incorrect_identifier_type(self):
        """Test that the appropriate error is raised for identifier/type mismatch."""
        with pytest.raises(ValueError) as cm:
            Molecule().from_smiles("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")

        assert "Improper identifier type" in str(cm.exconly())

    def test_read_inchikey_error(self):
        """Test that the correct error is raised when reading an InChIKey"""
        with pytest.raises(ValueError) as cm:
            Molecule().from_inchi("InChIKey=UHOVQNZJYSORNB-UHFFFAOYSA-N")

        assert "InChIKey is a write-only format" in str(cm.exconly())


class InChIParsingTest:
    def compare(self, inchi, u_indices=None, p_indices=None):
        u_layer = U_LAYER_PREFIX + U_LAYER_SEPARATOR.join(map(str, u_indices)) if u_indices else None
        p_layer = P_LAYER_PREFIX + P_LAYER_SEPARATOR.join(map(str, p_indices)) if p_indices else None

        aug_inchi = compose_aug_inchi(inchi, u_layer, p_layer)

        mol = from_augmented_inchi(Molecule(), aug_inchi)
        ConsistencyChecker.check_multiplicity(mol.get_radical_count(), mol.multiplicity)

        for at in mol.atoms:
            ConsistencyChecker.check_partial_charge(at)

        spc = Species(molecule=[mol])
        spc.generate_resonance_structures()

        ignore_prefix = r"(InChI=1+)(S*)/"
        aug_inchi_expected = re.split(ignore_prefix, aug_inchi)[-1]
        aug_inchi_computed = re.split(ignore_prefix, spc.get_augmented_inchi())[-1]
        assert aug_inchi_expected == aug_inchi_computed

        return mol

    def test_ethane_parsing(self):
        inchi = "C2H6/c1-2/h1-2H3"
        self.compare(inchi)

    def test_ethyl_parsing(self):
        inchi = "C2H5/c1-2/h1H2,2H3"
        u_indices = [1]
        self.compare(inchi, u_indices)

    def test_ch3_parsing(self):
        inchi = "CH3/h1H3"
        u_indices = [1]
        self.compare(inchi, u_indices)

    def test_h2_parsing(self):
        inchi = "H2/h1H"
        self.compare(inchi)

    def test_c2h4_biradical_parsing(self):
        inchi = "C2H4/c1-2/h1-2H2"
        u_indices = [1, 2]
        self.compare(inchi, u_indices)

    def test_c2h3_triradical_parsing(self):
        inchi = "C2H3/c1-2/h1H,2H2"
        u_indices = [1, 1, 2]
        self.compare(inchi, u_indices)

    def test_c3h6_biradical_parsing(self):
        inchi = "C3H6/c1-3-2/h1-3H2"
        u_indices = [1, 2]
        self.compare(inchi, u_indices)

    def test_c2h3o3(self):
        inchi = "C2H3O3/c1-2(3)5-4/h4H,1H2"
        u_indices = [1]
        self.compare(inchi, u_indices)

    def test_c2h2(self):
        inchi = "C2H2/c1-2/h1-2H"
        u_indices = [1, 2]
        self.compare(inchi, u_indices)

    def test_o2(self):
        inchi = "O2/c1-2"
        u_indices = [1, 2]
        self.compare(inchi, u_indices)

    def test_tri_radical_zwitter_mult4(self):
        inchi = "C6H11/c1-3-5-6-4-2/h5H,1-4,6H2"
        u_indices = [1, 2, 5]
        self.compare(inchi, u_indices)

    def test_tri_radical_double_bond_mult4(self):
        inchi = "C4H7/c1-3-4-2/h3H,1-2,4H2"
        u_indices = [1, 2, 3]
        self.compare(inchi, u_indices)

    def test_tri_radical2_double_bond_mult4(self):
        inchi = "C6H9/c1-4-6(3)5-2/h1,4-6H,2H2,3H3"
        u_indices = [1, 2, 5]
        self.compare(inchi, u_indices)

    def test_quadri_radical_double_bond_zwitter_mult5(self):
        inchi = "C8H14/c1-4-6-7-8(3)5-2/h5-6,8H,1-2,4,7H2,3H3"
        u_indices = [1, 2, 5, 6]
        self.compare(inchi, u_indices)

    def test_quadri2_double_bond_mult5(self):
        inchi = "C8H14/c1-5-7(3)8(4)6-2/h5-8H,1-2H2,3-4H3"
        u_indices = [1, 2, 5, 6]
        self.compare(inchi, u_indices)

    def test_c5h6o(self):
        inchi = "C5H6O/c6-5-3-1-2-4-5/h1-3,5H,4H2"
        u_indices = [2, 6]
        self.compare(inchi, u_indices)

    def test_c5h6o2(self):
        inchi = "C5H6O/c1-5-3-2-4-6-5/h2-5H,1H2"
        u_indices = [1, 3]
        self.compare(inchi, u_indices)

    def test_c5h6o3(self):
        inchi = "C5H6O/c1-5-3-2-4-6-5/h2-5H,1H2"
        u_indices = [1, 2, 3, 4]
        self.compare(inchi, u_indices)

    @pytest.mark.skip(reason="WIP")
    def test_co(self):
        inchi = "CO/c1-2"
        p_indices = [1, 2]
        mol = self.compare(inchi, [], p_indices)

        assert mol.atoms[1].lone_pairs == 1  # Oxygen

        assert mol.atoms[0].charge == -1
        assert mol.atoms[1].charge == 1

    def test_triplet_methylene(self):
        inchi = "CH2/h1H2"

        u_indices = [1, 1]
        self.compare(inchi, u_indices)

    def test_singlet_methylene(self):
        inchi = "CH2/h1H2"

        p_indices = [1]
        self.compare(inchi, u_indices=[], p_indices=p_indices)

    def test_c4h6o(self):
        inchi = "C4H6O/c1-2-3-4-5/h2H,3H2,1H3"
        u_indices = [2, 4]
        mol = self.compare(inchi, u_indices)
        for at in mol.atoms:
            if at.is_oxygen():
                assert at.lone_pairs == 2

    def test_c6h6(self):
        inchi = "C6H6/c1-3-5-6-4-2/h1,6H,2,5H2"
        u_indices = [1, 3]
        self.compare(inchi, u_indices)

    def test_c4h6o2(self):
        inchi = "C4H6O/c1-2-3-4-5/h2,4H,1,3H2"
        u_indices = [4, 5]
        self.compare(inchi, u_indices)

    def test_co_triplet(self):
        adjlist = """
        multiplicity 3
        1 C u2 p0 c0 {2,D}
        2 O u0 p2 c0 {1,D}

        """
        spc = Species(molecule=[Molecule().from_adjacency_list(adjlist)])
        aug_inchi = spc.get_augmented_inchi()

        assert Species(molecule=[Molecule().from_augmented_inchi(aug_inchi)]).is_isomorphic(spc) == True

    def test_ccco_triplet(self):
        adjlist = """
        multiplicity 3
1 C u0 p0 c0 {2,D} {5,S} {6,S}
2 C u0 p0 c0 {1,D} {3,S} {7,S}
3 C u1 p0 c0 {2,S} {4,S} {8,S}
4 O u1 p2 c0 {3,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
        """
        mol = Molecule().from_adjacency_list(adjlist)

        spc = Species(molecule=[mol])
        spc.generate_resonance_structures()
        aug_inchi = spc.get_augmented_inchi()

        assert Species(molecule=[Molecule().from_augmented_inchi(aug_inchi)]).is_isomorphic(spc) == True

    def test_c3h4(self):
        inchi = "C3H4/c1-3-2/h1,3H,2H2"
        u_indices = [1, 1]
        self.compare(inchi, u_indices)

    def test_c6h8o2(self):
        inchi = "C6H8O2/c1-3-5(7)6(8)4-2/h3-6H,1-2H2"
        u_indices = [7, 8]
        self.compare(inchi, u_indices)

    def test_c3h3o3(self):
        inchi = "C3H3O3/c1-2-5-3-6-4/h1-3H"
        u_indices = [1, 3, 4]
        self.compare(inchi, u_indices)

    def test_ch2o2(self):
        inchi = "CH2O2/c2-1-3/h1H,(H,2,3)"
        u_indices = [1, 2]
        self.compare(inchi, u_indices)

    def test_c2h2o3(self):
        inchi = "C2H2O3/c1-5-2(3)4/h1H2"
        u_indices = [1, 3]
        self.compare(inchi, u_indices)

    def test_c3h4o4(self):
        inchi = "C3H4O4/c4-3(5)1-2-7-6/h1-3,6H"
        u_indices = [4, 5]
        self.compare(inchi, u_indices)

    @pytest.mark.skip(reason="WIP")
    def test_c6h6o4(self):
        """
        This test used to pass with OpenBabel < 3.0, but I think the inchi is invalid?
        or at least not standard.
        OpenBabel reports:
            Problems/mismatches: Mobile-H( Hydrogens: Locations or number, Number; Charge(s): Do not match)
        and cactus.nci.nih.gov converts it to InChI=1S/C6H7O4/c1-2-4-9-6(7)3-5-10-8/h2-3,8H,1,5H2/q+1
        which at least doesn't make OpenBabel complain. However, both have a net charge
        and cause RMG to crash. I'm not sure what the molecule was ever supposed to represent.
        """
        inchi = "InChI=1S/C6H6O4/c1-2-4-9-6(7)3-5-10-8/h2-3H,1,5H2"
        u_indices = [1, 3, 4, 8]
        self.compare(inchi, u_indices)

    def test_c3h2o3(self):
        inchi = "InChI=1S/C3H2O3/c1-2-3(4)6-5/h1H2"
        u_indices = [2, 5]

        self.compare(inchi, u_indices)

    def test_c6h6o6(self):
        inchi = "C6H6O6/c7-6(2-5-12-9)10-3-1-4-11-8/h1,7H,4-5H2"
        u_indices = [2, 3, 8, 9]
        self.compare(inchi, u_indices)

    def test_c3h2(self):
        inchi = "C3H2/c1-3-2/h1-2H"
        u_indices = [1, 1]
        self.compare(inchi, u_indices)

    def test_c3h4(self):
        inchi = "InChI=1S/C3H4/c1-3-2/h1,3H,2H2"
        u_indices = [1, 1]
        self.compare(inchi, u_indices)

    def test_c6h8(self):
        inchi = "InChI=1S/C6H8/c1-3-5-6-4-2/h1,4H,2,5-6H2"
        u_indices = [1, 1, 3, 3]
        self.compare(inchi, u_indices)

    def test_c6h10(self):
        inchi = "InChI=1S/C6H10/c1-3-5-6-4-2/h3-4H,1-2,5-6H2"
        u_indices = [1, 3]
        self.compare(inchi, u_indices)

    def test_ammonia(self):
        inchi = "InChI=1S/H3N/h1H3"
        self.compare(inchi)

    @pytest.mark.skip(reason="WIP")
    def test_ammonium(self):
        """
        has same inchi as ammonia but gets a proton layer: /p+1
        """
        inchi = "InChI=1S/H3N/h1H3/p+1"
        self.compare(inchi)

    def test_h2s(self):
        inchi = "InChI=1S/H2S/h1H2"
        self.compare(inchi)

    def test_pyridine(self):
        inchi = "InChI=1S/C5H5N/c1-2-4-6-5-3-1/h1-5H"
        self.compare(inchi)

    def test_pyrimidine(self):
        inchi = "InChI=1S/C4H4N2/c1-2-5-4-6-3-1/h1-4H"
        self.compare(inchi)

    @pytest.mark.skip(reason="WIP")
    def test_nitrate(self):
        """
        - Mobile H spread over oxygen 2, 3, 4
        - Negative charge (3 lone pairs) spread out over oxygen 2, 3, 4
        - Nitrogen 1 positively charged

        """
        inchi = "InChI=1S/HNO3/c2-1(3)4/h(H,2,3,4)"
        p_indices = [-1, 3, 3, 3]  # ???
        self.compare(inchi, [], p_indices)

    def test_no(self):
        inchi = "InChI=1S/NO/c1-2"
        u_indices = [1]
        self.compare(inchi, u_indices)

    def test_isotopic_molecule_1(self):
        """Test that we can parse an InChI for an isotopic molecule."""
        mol = Molecule().from_inchi("InChI=1S/CH4/h1H4/i1+1")

        assert len(mol.atoms), 4
        assert [atom.element.isotope for atom in mol.atoms].count(13) == 1

    def test_isotopic_molecule_2(self):
        """Test that we can parse an InChI for an isotopic molecule."""
        mol = Molecule().from_inchi("InChI=1S/C2H6/c1-2/h1-2H3/i1+1")

        assert len(mol.atoms), 6
        assert [atom.element.isotope for atom in mol.atoms].count(13) == 1
