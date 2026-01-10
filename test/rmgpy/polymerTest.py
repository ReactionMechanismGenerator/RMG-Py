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
This module contains unit tests of the rmgpy.polymer module.
"""

import numpy as np
import pytest

from rmgpy.exceptions import InputError
from rmgpy.molecule import Molecule
from rmgpy.polymer import (
    Polymer,
    stitch_molecules_by_labeled_atoms,
    find_labeled_atom,
    get_label_1_label_2_atoms,
    LABELS_1,
    LABELS_2,
)
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData


class TestPolymer:
    """
    Contains unit tests for the Polymer class.
    """
    @pytest.fixture(autouse=True)
    def setup_species(self):
        """
        A method that is run before each unit test in this class.
        """
        ps_adj = """multiplicity 3
                    1 *1 C u1 p0 c0 {2,S} {9,S} {10,S}
                    2 *2 C u1 p0 c0 {1,S} {3,S} {11,S}
                    3    C u0 p0 c0 {2,S} {4,S} {8,D}
                    4    C u0 p0 c0 {3,S} {5,D} {12,S}
                    5    C u0 p0 c0 {4,D} {6,S} {13,S}
                    6    C u0 p0 c0 {5,S} {7,D} {14,S}
                    7    C u0 p0 c0 {6,D} {8,S} {15,S}
                    8    C u0 p0 c0 {3,D} {7,S} {16,S}
                    9    H u0 p0 c0 {1,S}
                    10   H u0 p0 c0 {1,S}
                    11   H u0 p0 c0 {2,S}
                    12   H u0 p0 c0 {4,S}
                    13   H u0 p0 c0 {5,S}
                    14   H u0 p0 c0 {6,S}
                    15   H u0 p0 c0 {7,S}
                    16   H u0 p0 c0 {8,S}"""
        ps_feature_adj = """multiplicity 4
                            1  *1 C u1 p0 c0 {2,S} {9,S} {10,S}
                            2  *2 C u1 p0 c0 {1,S} {3,S} {11,S}
                            3     C u0 p0 c0 {2,S} {4,S} {8,D}
                            4     C u0 p0 c0 {3,S} {5,D} {12,S}
                            5     C u0 p0 c0 {4,D} {6,S} {13,S}
                            6     C u1 p0 c0 {5,S} {7,D}
                            7     C u0 p0 c0 {6,D} {8,S} {14,S}
                            8     C u0 p0 c0 {3,D} {7,S} {15,S}
                            9     H u0 p0 c0 {1,S}
                            10    H u0 p0 c0 {1,S}
                            11    H u0 p0 c0 {2,S}
                            12    H u0 p0 c0 {4,S}
                            13    H u0 p0 c0 {5,S}
                            14    H u0 p0 c0 {7,S}
                            15    H u0 p0 c0 {8,S}"""
        ps_smiles = '[CH2][CH]c1ccccc1'
        pe_smiles = '[CH2][CH2]'
        self.polymer_1 = Polymer(
                 label='PS_1',
                 monomer=ps_adj,
                 end_groups=['[CH3]', '[H]'],
                 cutoff=3,
                 Mn=5000.0,
                 Mw=6000.0,
                 initial_mass=1.0,
        )

        self.polymer_2 = Polymer(
                 label='PS_2',
                 monomer=ps_smiles,
                 end_groups=['[CH3]', '[H]'],
                 cutoff=5,
                 Mn=3000.0,
                 Mw=10000.0,
                 initial_mass=1.0,
        )

        self.polymer_3 = Polymer(
                 label='PE_1',
                 monomer=pe_smiles,
                 end_groups=['[H]', '[H]'],
                 cutoff=10,
                 Mn=1000.0,
                 Mw=2500.0,
                 initial_mass=1.0,
        )

        self.polymer_4 = Polymer(
                 label='PS_3',
                 monomer=ps_smiles,
                 feature_monomer=ps_feature_adj,
                 end_groups=['[CH3]', '[H]'],
                 cutoff=5,
                 Mn=3000.0,
                 Mw=10000.0,
                 initial_mass=1.0,
        )

        self.ethylene_diradical_labeled_adj = """multiplicity 3
                                                 1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
                                                 2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}
                                                 3    H u0 p0 c0 {1,S}
                                                 4    H u0 p0 c0 {1,S}
                                                 5    H u0 p0 c0 {2,S}
                                                 6    H u0 p0 c0 {2,S}"""

        yield
        # teardown here if necessary

    def test_repr(self):
        """
        Test Polymer representation.
        """
        expected_repr_1 = "<Polymer 'PS_1' Mn=5000.0 Mw=6000.0 Cutoff=3>"
        expected_repr_2 = "<Polymer 'PS_2' Mn=3000.0 Mw=10000.0 Cutoff=5>"
        expected_repr_3 = "<Polymer 'PE_1' Mn=1000.0 Mw=2500.0 Cutoff=10>"
        repr_1 = repr(self.polymer_1)
        repr_2 = repr(self.polymer_2)
        repr_3 = repr(self.polymer_3)
        assert repr_1 == expected_repr_1
        assert repr_2 == expected_repr_2
        assert repr_3 == expected_repr_3

    def test_equality(self):
        """Test that we can perform equality comparison with Species objects"""
        assert self.polymer_1 != self.polymer_2
        assert self.polymer_1 == self.polymer_1
        assert self.polymer_2 == self.polymer_2

    def test_to_adjacency_list(self):
        """
        Test that to_adjacency_list() works as expected.
        """
        adj = self.polymer_1.copy().to_adjacency_list()
        partial_expected_adj = """PS_1
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  H u0 p0 c0 {1,S}
3  H u0 p0 c0 {1,S}
4  H u0 p0 c0 {1,S}
5  C u0 p0 c0 {6,S} {13,S} {14,S} {22,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,B} {12,B}
8  C u0 p0 c0 {7,B} {9,B} {16,S}
9  C u0 p0 c0 {8,B} {10,B} {17,S}
10 C u0 p0 c0 {9,B} {11,B} {18,S}
11 C u0 p0 c0 {10,B} {12,B} {19,S}
12 C u0 p0 c0 {7,B} {11,B} {20,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {12,S}
21 C u0 p0 c0 {22,S} {29,S} {30,S} {38,S}
22 C u0 p0 c0 {5,S} {21,S} {23,S} {31,S}
23 C u0 p0 c0 {22,S} {24,B} {28,B}
24 C u0 p0 c0 {23,B} {25,B} {32,S}
25 C u0 p0 c0 {24,B} {26,B} {33,S}
26 C u0 p0 c0 {25,B} {27,B} {34,S}
27 C u0 p0 c0 {26,B} {28,B} {35,S}
28 C u0 p0 c0 {23,B} {27,B} {36,S}
29 H u0 p0 c0 {21,S}
30 H u0 p0 c0 {21,S}
31 H u0 p0 c0 {22,S}
32 H u0 p0 c0 {24,S}
33 H u0 p0 c0 {25,S}
34 H u0 p0 c0 {26,S}
35 H u0 p0 c0 {27,S}
36 H u0 p0 c0 {28,S}
37 C u0 p0 c0 {38,S} {45,S} {46,S} {53,S}
38 C u0 p0 c0 {21,S} {37,S} {39,S} {47,S}
39 C u0 p0 c0 {38,S} {40,B} {44,B}
40 C u0 p0 c0 {39,B} {41,B} {48,S}
41 C u0 p0 c0 {40,B} {42,B} {49,S}
42 C u0 p0 c0 {41,B} {43,B} {50,S}
43 C u0 p0 c0 {42,B} {44,B} {51,S}
44 C u0 p0 c0 {39,B} {43,B} {52,S}
45 H u0 p0 c0 {37,S}
46 H u0 p0 c0 {37,S}
47 H u0 p0 c0 {38,S}
48 H u0 p0 c0 {40,S}
49 H u0 p0 c0 {41,S}
50 H u0 p0 c0 {42,S}
51 H u0 p0 c0 {43,S}
52 H u0 p0 c0 {44,S}
53 H u0 p0 c0 {37,S}


PS_1
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  H u0 p0 c0 {1,S}
3  H u0 p0 c0 {1,S}
4  H u0 p0 c0 {1,S}
5  C u0 p0 c0 {6,S} {13,S} {14,S} {22,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,S} {12,D}
8  C u0 p0 c0 {7,S} {9,D} {16,S}
9  C u0 p0 c0 {8,D} {10,S} {17,S}
10 C u0 p0 c0 {9,S} {11,D} {18,S}
11 C u0 p0 c0 {10,D} {12,S} {19,S}
12 C u0 p0 c0 {7,D} {11,S} {20,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {12,S}
21 C u0 p0 c0 {22,S} {29,S} {30,S} {38,S}
22 C u0 p0 c0 {5,S} {21,S} {23,S} {31,S}
23 C u0 p0 c0 {22,S} {24,S} {28,D}
24 C u0 p0 c0 {23,S} {25,D} {32,S}
25 C u0 p0 c0 {24,D} {26,S} {33,S}
26 C u0 p0 c0 {25,S} {27,D} {34,S}
27 C u0 p0 c0 {26,D} {28,S} {35,S}
28 C u0 p0 c0 {23,D} {27,S} {36,S}
29 H u0 p0 c0 {21,S}
30 H u0 p0 c0 {21,S}
31 H u0 p0 c0 {22,S}
32 H u0 p0 c0 {24,S}
33 H u0 p0 c0 {25,S}
34 H u0 p0 c0 {26,S}
35 H u0 p0 c0 {27,S}
36 H u0 p0 c0 {28,S}
37 C u0 p0 c0 {38,S} {45,S} {46,S} {53,S}
38 C u0 p0 c0 {21,S} {37,S} {39,S} {47,S}
39 C u0 p0 c0 {38,S} {40,S} {44,D}
40 C u0 p0 c0 {39,S} {41,D} {48,S}
41 C u0 p0 c0 {40,D} {42,S} {49,S}
42 C u0 p0 c0 {41,S} {43,D} {50,S}
43 C u0 p0 c0 {42,D} {44,S} {51,S}
44 C u0 p0 c0 {39,D} {43,S} {52,S}
45 H u0 p0 c0 {37,S}
46 H u0 p0 c0 {37,S}
47 H u0 p0 c0 {38,S}
48 H u0 p0 c0 {40,S}
49 H u0 p0 c0 {41,S}
50 H u0 p0 c0 {42,S}
51 H u0 p0 c0 {43,S}
52 H u0 p0 c0 {44,S}
53 H u0 p0 c0 {37,S}"""
        assert partial_expected_adj in adj

    def test_to_smiles(self):
        """
        Test that to_smiles() works as expected.
        """
        smiles = self.polymer_1.monomer.copy(deep=True).to_smiles()
        expected_smiles = ['[CH2][CH]c1ccccc1', '[CH2][CH]C1=CC=CC=C1']
        assert smiles in expected_smiles

    def test_copy(self):
        """Test that we can make a copy of a Polymer object."""
        poly_cp = self.polymer_1.copy()
        assert id(self.polymer_1) != id(poly_cp)
        assert self.polymer_1.is_isomorphic(poly_cp)
        assert self.polymer_1.label == poly_cp.label
        assert self.polymer_1.index == poly_cp.index

    def test_fingerprint_property(self):
        """Test that the fingerprint property works"""
        assert self.polymer_1.fingerprint == "Polymer_C08H08N00O00S00_3"
        assert self.polymer_2.fingerprint == "Polymer_C08H08N00O00S00_5"
        assert self.polymer_3.fingerprint == "Polymer_C02H04N00O00S00_10"

    def test_baseline_proxy(self):
        """Test that the baseline_proxy property works"""
        assert len(self.polymer_1.baseline_proxy.molecule[0].atoms) == 53
        assert len(self.polymer_2.baseline_proxy.molecule[0].atoms) == 53
        assert len(self.polymer_3.baseline_proxy.molecule[0].atoms) == 20

    def test_feature_proxy(self):
        """Test that the feature_proxy property works"""
        assert self.polymer_1.feature_proxy is None
        assert self.polymer_2.feature_proxy is None
        assert self.polymer_3.feature_proxy is None
        assert len(self.polymer_4.feature_monomer.atoms) == 15
        assert len(self.polymer_4.feature_proxy.molecule[0].atoms) == 52

    def test_is_isomorphic(self):
        """Test that the Polymer.is_isomorphic works"""
        spc_1 = Species(smiles='CC(CC(CC(C)c1ccccc1)c1ccccc1)c1ccccc1')
        mol_1 = Molecule(smiles='[CH2][CH]c1ccccc1')
        pol_1 = Polymer(label='PS_10',
                        monomer='[CH2][CH]c1ccccc1',
                        end_groups=['[CH3]', '[H]'],
                        cutoff=30,
                        Mn=7000.0,
                        Mw=8000.0,
                        initial_mass=10.0)
        assert not spc_1.is_isomorphic(pol_1)
        assert not pol_1.is_isomorphic(spc_1)
        assert not pol_1.is_isomorphic(mol_1)
        assert self.polymer_1.is_isomorphic(self.polymer_1)
        assert self.polymer_1.is_isomorphic(pol_1)
        assert any(spc_1.is_isomorphic(monomer) for monomer in pol_1.molecule)

    def test_polymer_label(self):
        """Test that the polymer label"""
        assert self.polymer_1.label == "PS_1"
        assert self.polymer_2.label == "PS_2"
        assert self.polymer_3.label == "PE_1"

    def test_calculate_moments_from_distribution(self):
        """Test that the moments are calculated correctly from the distribution"""
        expected_moments = [2.00000000e-01, 9.60163842e+00, 5.53148762e+02]
        assert all(np.isclose(v1, v2) for v1, v2 in zip(self.polymer_1.moments, expected_moments))

        expected_moments = [3.33333333e-01, 9.60163842e+00, 9.21914603e+02]
        assert all(np.isclose(v1, v2) for v1, v2 in zip(self.polymer_2.moments, expected_moments))

        expected_moments = [1.00000000e+00, 3.56466018e+01, 3.17670055e+03]
        assert all(np.isclose(v1, v2) for v1, v2 in zip(self.polymer_3.moments, expected_moments))

    def test_get_closing_moment(self):
        """Test that the closing moment is calculated correctly"""
        expected_closing_moment_1 = 3.8240e+4
        assert np.isclose(self.polymer_1.get_closing_moment(), expected_closing_moment_1)

        expected_closing_moment_2 = 2.95063022e+5
        assert np.isclose(self.polymer_2.get_closing_moment(), expected_closing_moment_2)

        expected_closing_moment_3 = 7.07741122e+5
        assert np.isclose(self.polymer_3.get_closing_moment(), expected_closing_moment_3)

        assert np.isclose(self.polymer_3.get_closing_moment([2.00000000e-01, 9.60163842e+00, 5.53148762e+02]),
                          expected_closing_moment_1)

        assert self.polymer_3.get_closing_moment([1, 2, -5]) == 0.0
        assert self.polymer_3.get_closing_moment([1, 2, 1e-30]) == 0.0

    def test_validate_monomer_raises_when_not_labeled_and_not_diradical(self):
        """CC should be closed shell, no labels, not a diradical => error"""
        with pytest.raises(InputError):
            Polymer(label="bad",
                    monomer="CC",
                    end_groups=["[H]", "[H]"],
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0)

    def test_validate_monomer_auto_labels_when_diradical(self):
        """Use PE example from your tests: should have 2 radicals, get *1/*2 assigned"""
        p = Polymer(label="auto_label",
                    monomer="[CH2][CH2]",
                    end_groups=["[H]", "[H]"],
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0, )
        assert find_labeled_atom(p.monomer, LABELS_1) is not None
        assert find_labeled_atom(p.monomer, LABELS_2) is not None

    def test_validate_monomer_rejects_invalid_type(self):
        """monomer must be str or Molecule"""
        with pytest.raises(InputError):
            Polymer(label="bad_type",
                    monomer=123,  # not str/Molecule
                    end_groups=["[H]", "[H]"],
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0)

    def test_validate_end_groups_default_assigns_labels(self):
        """When end_groups=None, should assign [H] with correct labels"""
        p = Polymer(label="default_ends",
                    monomer=self.ethylene_diradical_labeled_adj,
                    end_groups=None,  # defaults to [H], [H]
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0)
        assert len(p.end_groups) == 2
        head, tail = p.end_groups
        assert find_labeled_atom(head, LABELS_1) is not None
        assert find_labeled_atom(tail, LABELS_2) is not None

    def test_validate_end_groups_wrong_length_raises(self):
        """end_groups must be length 2"""
        with pytest.raises(InputError):
            Polymer(label="bad_ends_len",
                    monomer=self.ethylene_diradical_labeled_adj,
                    end_groups=["[H]"],  # wrong length
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0)

    def test_validate_end_groups_non_radical_raises(self):
        """C is not a radical end-group"""
        with pytest.raises(InputError):
            Polymer(label="bad_end_rad",
                    monomer=self.ethylene_diradical_labeled_adj,
                    end_groups=["C", "[H]"],
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0)

    def test_validate_end_groups_invalid_labels_raises(self):
        """end-groups must have correct *1/*2 labels"""
        bad = Molecule().from_adjacency_list(_methyl_radical_adj("*3"))  # invalid label
        ok = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        with pytest.raises(InputError):
            Polymer(label="bad_end_label",
                    monomer=self.ethylene_diradical_labeled_adj,
                    end_groups=[bad, ok],
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0)

    def test_find_label_1_2_atoms_helper(self):
        """Test the get_label_1_label_2_atoms helper function."""
        mol = Molecule().from_adjacency_list(self.ethylene_diradical_labeled_adj)
        pair = get_label_1_label_2_atoms(mol)
        assert pair is not None
        i1, i2 = pair
        assert i1 != i2
        assert mol.atoms[i1].label in LABELS_1
        assert mol.atoms[i2].label in LABELS_2

    def test_stitch_returns_none_if_any_input_none(self):
        """Test that stitch_molecules_by_labeled_atoms returns None if any input is None."""
        mol = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        assert stitch_molecules_by_labeled_atoms(None, mol) is None
        assert stitch_molecules_by_labeled_atoms(mol, None) is None

    def test_stitch_for_p1(self):
        """Test stitch_molecules_by_labeled_atoms for polymer_1 head wing + methyl radical."""
        p = self.polymer_1.copy()
        head_wing = p._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        scission_fragment = stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert scission_fragment.to_smiles() in ['CCC(C)C1=CC=CC=C1']

    def test_stitch_raises_when_missing_labels(self):
        """Test that stitch_molecules_by_labeled_atoms raises ValueError when labels missing."""
        left = Molecule().from_adjacency_list(_methyl_radical_adj("*1"))
        right = Molecule(smiles="[CH3]")  # radical but no *2 label
        with pytest.raises(ValueError):
            stitch_molecules_by_labeled_atoms(left, right)

    def test_stitch_raises_when_stitch_sites_not_mono_radicals(self):
        """Test that stitch_molecules_by_labeled_atoms raises ValueError when stitch sites not mono-radicals."""
        left = Molecule().from_adjacency_list(_methyl_closed_shell_labeled_adj("*1"))  # labeled but u0
        right = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        with pytest.raises(ValueError):
            stitch_molecules_by_labeled_atoms(left, right)

    def test_stitch_does_not_mutate_inputs_and_clears_labels_in_product(self):
        """Test that stitch_molecules_by_labeled_atoms does not mutate inputs and clears labels in product."""
        left = Molecule().from_adjacency_list(_methyl_radical_adj("*1"))
        right = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))

        left_before_labels = [a.label for a in left.atoms]
        right_before_labels = [a.label for a in right.atoms]
        left_before_rad = left.get_radical_count()
        right_before_rad = right.get_radical_count()

        merged = stitch_molecules_by_labeled_atoms(left, right)
        assert merged is not None

        # inputs unchanged (function deep-copies)
        assert [a.label for a in left.atoms] == left_before_labels
        assert [a.label for a in right.atoms] == right_before_labels
        assert left.get_radical_count() == left_before_rad
        assert right.get_radical_count() == right_before_rad

        # product has no stitch labels remaining at the join sites
        assert find_labeled_atom(merged, LABELS_1) is None
        assert find_labeled_atom(merged, LABELS_2) is None
        assert merged.get_radical_count() == 0

    def test_proxy_species_cached_and_has_no_remaining_labels(self):
        """Test that Polymer.baseline_proxy is cached and has no remaining labels."""
        p = Polymer(label="proxy_label_cleanup",
                    monomer=self.ethylene_diradical_labeled_adj,
                    end_groups=["[H]", "[H]"],
                    cutoff=3,
                    Mn=1000.0,
                    Mw=2000.0,
                    initial_mass=1.0)
        spc1 = p.baseline_proxy
        spc2 = p.baseline_proxy
        assert spc1 is spc2  # caching

        mol = spc1.molecule[0]
        assert find_labeled_atom(mol, LABELS_1) is None
        assert find_labeled_atom(mol, LABELS_2) is None

    def test_is_isomorphic_feature_mismatch_false(self):
        """Test that is_isomorphic returns False when feature_monomer differs from monomer."""
        base = Polymer(label="base",
                       monomer=self.ethylene_diradical_labeled_adj,
                       end_groups=["[H]", "[H]"],
                       cutoff=3,
                       Mn=1000.0,
                       Mw=2000.0,
                       initial_mass=1.0)
        feat = Polymer(label="feat",
                       monomer=self.ethylene_diradical_labeled_adj,
                       feature_monomer='[CH2][CH]',
                       end_groups=["[H]", "[H]"],
                       cutoff=3,
                       Mn=1000.0,
                       Mw=2000.0,
                       initial_mass=1.0)
        assert base.is_isomorphic(feat) is False
        assert feat.is_isomorphic(base) is False

    def test_init_from_moments_recovers_distribution(self):
        """Test that initializing a Polymer from moments recovers the same Mn/Mw."""
        p1 = Polymer(label="p1",
                     monomer=self.ethylene_diradical_labeled_adj,
                     end_groups=["[H]", "[H]"],
                     cutoff=3,
                     Mn=1000.0,
                     Mw=2000.0,
                     initial_mass=1.0)
        p2 = Polymer(label="p2",
                     monomer=self.ethylene_diradical_labeled_adj,
                     end_groups=["[H]", "[H]"],
                     cutoff=3,
                     moments=p1.moments.tolist(),
                     initial_mass=1.0)
        assert np.isclose(p2.Mn, p1.Mn)
        assert np.isclose(p2.Mw, p1.Mw)

    def test_wing_pattern(self):
        """Test that _stitch_wing produces expected patterns."""
        head_wing = self.polymer_1._wing_pattern('head')
        tail_wing = self.polymer_1._wing_pattern('tail')

        atom_1 = find_labeled_atom(head_wing, LABELS_1)
        assert atom_1 is None
        atom_2 = find_labeled_atom(head_wing, LABELS_2)
        assert atom_2 is None

        for atom in head_wing.atoms:
            assert atom.radical_electrons == 0
            assert atom.label == ''
        for atom in tail_wing.atoms:
            assert atom.radical_electrons == 0
            assert atom.label == ''

    def test_assert_end_group(self):
        """Test _assert_end_group method."""
        tail = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        self.polymer_1._assert_end_group(tail, want_label="*2")

        tail = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        with pytest.raises(ValueError):
            self.polymer_1._assert_end_group(tail, want_label="*1")

        bad = Molecule().from_adjacency_list("""multiplicity 2
                                               1 *2 C u1 p0 c0 {2,S} {3,S} {4,S}
                                               2 *1 H u0 p0 c0 {1,S}
                                               3 H u0 p0 c0 {1,S}
                                               4 H u0 p0 c0 {1,S}""")
        with pytest.raises(ValueError):
            self.polymer_1._assert_end_group(bad, want_label="*2")

        closed = Molecule().from_adjacency_list(_methyl_closed_shell_labeled_adj("*2"))
        with pytest.raises(ValueError):
            self.polymer_1._assert_end_group(closed, want_label="*2")

        bad = Molecule().from_adjacency_list("""1 *2 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
                                               2 H u0 p0 c0 {1,S}
                                               3 H u0 p0 c0 {1,S}
                                               4 H u0 p0 c0 {1,S}
                                               5 H u0 p0 c0 {1,S}""")
        with pytest.raises(ValueError):
            self.polymer_1._assert_end_group(bad, want_label="*2")

    def test_assert_feature_unit(self):
        """Test _assert_feature_unit method."""

        feat = Molecule().from_adjacency_list(self.ethylene_diradical_labeled_adj)
        self.polymer_1._assert_feature_unit(feat)

        bad = Molecule().from_adjacency_list("""multiplicity 2
                                               1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
                                               2    C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
                                               3    H u0 p0 c0 {1,S}
                                               4    H u0 p0 c0 {1,S}
                                               5    H u0 p0 c0 {2,S}
                                               6    H u0 p0 c0 {2,S}
                                               7    H u0 p0 c0 {2,S}""")
        with pytest.raises(ValueError):
            self.polymer_1._assert_feature_unit(bad)

        bad = Molecule().from_adjacency_list("""multiplicity 3
                                               1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
                                               2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}
                                               3 *3 H u0 p0 c0 {1,S}
                                               4    H u0 p0 c0 {1,S}
                                               5    H u0 p0 c0 {2,S}
                                               6    H u0 p0 c0 {2,S}""")
        with pytest.raises(ValueError):
            self.polymer_1._assert_feature_unit(bad)

        bad = Molecule().from_adjacency_list("""multiplicity 2
                                               1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
                                               2 *2 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
                                               3    H u0 p0 c0 {1,S}
                                               4    H u0 p0 c0 {1,S}
                                               5    H u0 p0 c0 {2,S}
                                               6    H u0 p0 c0 {2,S}
                                               7    H u0 p0 c0 {2,S}""")
        with pytest.raises(ValueError):
            self.polymer_1._assert_feature_unit(bad)

        bad = Molecule().from_adjacency_list("""1 *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
                                               2 *2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
                                               3    H u0 p0 c0 {1,S}
                                               4    H u0 p0 c0 {1,S}
                                               5    H u0 p0 c0 {1,S}
                                               6    H u0 p0 c0 {2,S}
                                               7    H u0 p0 c0 {2,S}
                                               8    H u0 p0 c0 {2,S}""")
        with pytest.raises(ValueError):
            self.polymer_1._assert_feature_unit(bad)

    def test_extract_remainder_removes_atoms_preserves_bonds_and_strips_labels(self):
        """
        _extract_remainder() should:
          - remove the requested atoms (including attached H's if they are part of the removed subgraph)
          - preserve bonds between remaining atoms
          - strip labels on the copied atoms (for deterministic downstream relabeling)
          - not mutate the input molecule
        """
        # Build a simple 3-carbon chain with explicit hydrogens (RMG expands SMILES to explicit H)
        complex_mol = Molecule(smiles="CCC")
        complex_mol.atoms[0].label = "*1"
        complex_mol.atoms[2].label = "*2"

        # Remove the terminal carbon AND its attached hydrogens (as a subgraph match would)
        tail_c = complex_mol.atoms[2]
        atoms_to_remove = {tail_c}

        for nbr in tail_c.bonds:
            # RMG Atom has is_hydrogen(); this mirrors "remove the whole terminal group"
            if nbr.is_hydrogen():
                atoms_to_remove.add(nbr)

        remainder = self.polymer_1._extract_remainder(complex_mol, atoms_to_remove)

        # 1) No dangling atoms left in the remainder
        assert all(len(a.bonds) > 0 for a in remainder.atoms), "Remainder contains disconnected atoms."

        # 2) Heavy atoms: should now be an ethane-like fragment (2 carbons)
        remainder_c = [a for a in remainder.atoms if a.is_carbon()]
        assert len(remainder_c) == 2

        # 3) Bond preservation among the remaining carbons
        c0, c1 = remainder_c
        assert c1 in c0.bonds
        assert float(c0.bonds[c1].order) == 1.0

        # 4) Label stripping on copied atoms
        assert all(a.label == "" for a in remainder.atoms)

        # 5) Input molecule not mutated
        assert complex_mol.atoms[0].label == "*1"
        assert complex_mol.atoms[2].label == "*2"
        assert len(complex_mol.atoms) > len(remainder.atoms)


    def test_create_reacted_copy_modification_from_baseline_proxy(self):
        """
        If the reacted product contains both head and tail wings, create_reacted_copy()
        should return a Polymer with a non-None feature_monomer that has exactly one *1 and one *2.
        """
        p = self.polymer_1.copy()

        reacted_proxy = p.baseline_proxy.molecule[0].copy(deep=True)
        new_p = p.create_reacted_copy(reacted_proxy)

        assert new_p is not None
        assert isinstance(new_p, Polymer)
        assert new_p.feature_monomer is not None

        # Feature unit should contain exactly one *1 and one *2 label (per _assert_feature_unit contract)
        labels = [a.label for a in new_p.feature_monomer.atoms if a.label]
        assert labels.count("*1") == 1
        assert labels.count("*2") == 1

        # If original polymer had moments, new polymer should preserve distribution
        assert np.isclose(new_p.Mn, p.Mn)
        assert np.isclose(new_p.Mw, p.Mw)

    def test_create_reacted_copy_head_scission_returns_scission_tail_polymer(self):
        """
        Construct a fragment with ONLY the head wing present (no tail wing),
        so create_reacted_copy() should classify it as head-side scission and return a Polymer
        with a NEW tail end-group labeled *2 and mono-radical.
        """
        p = self.polymer_1.copy()

        # Build: [HeadWing] -- [methyl(*2)]
        # Head wing has an open *1 site (from baseline's remaining reactive site) that can stitch to *2.
        head_wing = p._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))

        scission_fragment = stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert scission_fragment is not None
        new_p = p.create_reacted_copy(scission_fragment)

        assert new_p is not None
        assert isinstance(new_p, Polymer)
        assert new_p.feature_monomer is None
        assert new_p.label.endswith("_scission_tail")

        # New tail end-group should be labeled *2 and mono-radical
        new_tail = new_p.end_groups[1]
        labels = {a.label for a in new_tail.atoms if a.label}
        assert labels == {"*2"}
        assert new_tail.get_radical_count() == 1

        # Transport heuristic: halve Mn/Mw
        assert np.isclose(new_p.Mn, p.Mn / 2.0)
        assert np.isclose(new_p.Mw, p.Mw / 2.0)

    def test_create_reacted_copy_tail_scission_returns_scission_head_polymer(self):
        """
        Construct a fragment with ONLY the tail wing present (no head wing),
        so create_reacted_copy() should classify it as tail-side scission and return a Polymer
        with a NEW head end-group labeled *1 and mono-radical.
        """
        p = self.polymer_1.copy()

        # Build: [methyl(*1)] -- [TailWing]
        # Tail wing has an open *2 site (from baseline's remaining reactive site) that can stitch to *1.
        tail_wing = p._stitch_wing("tail")
        methyl_star1 = Molecule().from_adjacency_list(_methyl_radical_adj("*1"))

        scission_fragment = stitch_molecules_by_labeled_atoms(methyl_star1, tail_wing)
        assert scission_fragment is not None
        # scission_fragment.atoms[0].increment_radical()
        new_p = p.create_reacted_copy(scission_fragment)

        assert new_p is not None
        assert isinstance(new_p, Polymer)
        assert new_p.feature_monomer is None
        assert new_p.label.endswith("_scission_tail")

        # New head end-group should be labeled *1 and mono-radical
        new_head = new_p.end_groups[0]
        labels = {a.label for a in new_head.atoms if a.label}
        assert labels == {"*1"}
        assert new_head.get_radical_count() == 1

        assert np.isclose(new_p.Mn, p.Mn / 2.0)
        assert np.isclose(new_p.Mw, p.Mw / 2.0)

    def test_create_reacted_copy_returns_none_for_small_molecule(self):
        """
        If reacted_proxy has no recognizable wings, create_reacted_copy() should return None.
        """
        p = self.polymer_1.copy()
        small = Molecule(smiles="CC")  # no head/tail wing subgraphs
        assert p.create_reacted_copy(small) is None


class TestPolymerThermo:
    """
    Contains unit tests for Polymer thermodynamic and property delegation.
    """

    @pytest.fixture(autouse=True)
    def setup_polymer(self):
        """
        A method that is run before each unit test in this class.
        """
        pe_adj = """
        multiplicity 3
        1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
        2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}
        3    H u0 p0 c0 {1,S}
        4    H u0 p0 c0 {1,S}
        5    H u0 p0 c0 {2,S}
        6    H u0 p0 c0 {2,S}
        """
        self.pe_polymer = Polymer(
            label='PE_Test',
            monomer=pe_adj,
            end_groups=['[CH3]', '[CH3]'],
            cutoff=4,
            Mn=1000.0,
            Mw=2000.0,
            initial_mass=1.0
        )
        self.proxy = self.pe_polymer.get_proxy_species()

        self.dummy_thermo = NASA(
            polynomials=[
                NASAPolynomial(coeffs=[1, 1, 1, 1, 1, 1, 1], Tmin=(298, 'K'), Tmax=(1000, 'K')),
                NASAPolynomial(coeffs=[2, 2, 2, 2, 2, 2, 2], Tmin=(1000, 'K'), Tmax=(3000, 'K'))
            ],
            Tmin=(298, 'K'), Tmax=(3000, 'K'), Cp0=(30, 'J/(mol*K)'), CpInf=(100, 'J/(mol*K)')
        )
        self.proxy.thermo = self.dummy_thermo

    def test_get_thermo_data_delegates_to_proxy(self):
        """Test that get_thermo_data returns the proxy's thermo object."""
        # Act
        thermo = self.pe_polymer.get_thermo_data()

        # Assert
        assert thermo is self.dummy_thermo
        assert self.pe_polymer.thermo is self.dummy_thermo  # Check sync behavior

    def test_thermo_properties_delegate_correctly(self):
        """Test get_enthalpy, entropy, heat_capacity, etc. return proxy values."""
        T = 500.0

        # Direct calls to polymer methods
        H_pol = self.pe_polymer.get_enthalpy(T)
        S_pol = self.pe_polymer.get_entropy(T)
        Cp_pol = self.pe_polymer.get_heat_capacity(T)

        # Expected calls to proxy thermo
        H_exp = self.dummy_thermo.get_enthalpy(T)
        S_exp = self.dummy_thermo.get_entropy(T)
        Cp_exp = self.dummy_thermo.get_heat_capacity(T)

        assert H_pol == H_exp
        assert S_pol == S_exp
        assert Cp_pol == Cp_exp

    def test_get_bulk_heat_capacity_scales_by_dp(self):
        """Test get_bulk_heat_capacity scales per-site Cp by DP."""
        T = 500.0
        DP = 100.0

        Cp_site = self.dummy_thermo.get_heat_capacity(T)
        Cp_bulk = self.pe_polymer.get_bulk_heat_capacity(T, DP)

        assert np.isclose(Cp_bulk, Cp_site * DP)

    def test_calculate_cp0_cpinf_delegate(self):
        """Test calculate_cp0 and calculate_cpinf don't crash and return floats."""
        # These methods typically require molecule structure analysis.
        # Since we haven't loaded the full RMG database, exact values aren't guaranteed,
        # but we verify the delegation path exists.
        cp0 = self.pe_polymer.calculate_cp0()
        cpinf = self.pe_polymer.calculate_cpinf()

        assert isinstance(cp0, float)
        assert isinstance(cpinf, float)

    def test_multiplicity_delegation(self):
        """Test multiplicity property."""
        # PE proxy (trimer) is constructed from radical stitching.
        # Since end_groups are radicals and monomer is diradical, proper stitching
        # results in a closed-shell alkane chain (multiplicity 1).
        mult = self.pe_polymer.multiplicity
        assert mult == 1

    def test_molecular_weight_delegation(self):
        """Test molecular_weight property returns proxy MW (per-site), not Mn (bulk)."""
        # Mn is 1000 g/mol
        # Proxy is trimer: CH3-(CH2CH2)3-CH3 (approx C8H18, MW ~114 g/mol)
        mw = self.pe_polymer.molecular_weight.value_si  # kg/mol
        mn = self.pe_polymer.Mn / 1000.0  # convert Mn to kg/mol
        assert mw < mn
        assert mw > 0.0

    def test_is_identical_delegation(self):
        """Test is_identical compares proxies."""
        # 1. Identity
        p2 = self.pe_polymer.copy()
        assert self.pe_polymer.is_identical(p2)

        # 2. Difference (Changing end groups changes the proxy)
        pe_adj = """
        multiplicity 3
        1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
        2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}
        3    H u0 p0 c0 {1,S}
        4    H u0 p0 c0 {1,S}
        5    H u0 p0 c0 {2,S}
        6    H u0 p0 c0 {2,S}
        """
        p3 = Polymer(
            label='PE_Diff',
            monomer=pe_adj,
            end_groups=['[H]', '[H]'],  # Hydrogen ends vs Methyl ends
            cutoff=4,Mn=1000.0,
            Mw=2000.0,
        )
        assert not self.pe_polymer.is_identical(p3)

    def test_transport_delegation(self):
        """Test generate_transport_data delegates to proxy."""
        dummy_trans = TransportData(sigma=(3.0, 'angstrom'), epsilon=(100.0, 'K'))
        self.proxy.transport_data = dummy_trans
        trans = self.pe_polymer.generate_transport_data()
        assert trans is dummy_trans
        assert self.pe_polymer.transport_data is dummy_trans


def _methyl_radical_adj(label: str) -> str:
    """CH3 rad with a label on the radical carbon"""
    return f"""multiplicity 2
               1 {label} C u1 p0 c0 {{2,S}} {{3,S}} {{4,S}}
               2 H u0 p0 c0 {{1,S}}
               3 H u0 p0 c0 {{1,S}}
               4 H u0 p0 c0 {{1,S}}"""


def _methyl_closed_shell_labeled_adj(label: str) -> str:
    """CH4 (closed shell) but labeled (bad for stitching: radical_electrons == 0)"""
    return f"""multiplicity 1
               1 {label} C u0 p0 c0 {{2,S}} {{3,S}} {{4,S}} {{5,S}}
               2 H u0 p0 c0 {{1,S}}
               3 H u0 p0 c0 {{1,S}}
               4 H u0 p0 c0 {{1,S}}
               5 H u0 p0 c0 {{1,S}}"""
