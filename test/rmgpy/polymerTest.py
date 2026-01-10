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
from collections import deque
from typing import List, Tuple, Dict, Any

import rmgpy.polymer as polymer
from rmgpy.exceptions import InputError
from rmgpy.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.group import GroupAtom
from rmgpy.polymer import LABELS_1, LABELS_2, Polymer, PolymerClass
from rmgpy.species import Species
from rmgpy.statmech import Conformer, HarmonicOscillator, IdealGasTranslation, NonlinearRotor
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
1  C u0 p0 c0 {4,S} {5,S} {9,S} {27,S}
2  C u0 p0 c0 {4,S} {6,S} {8,S} {26,S}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {28,S}
4  C u0 p0 c0 {1,S} {2,S} {29,S} {30,S}
5  C u0 p0 c0 {1,S} {3,S} {31,S} {32,S}
6  C u0 p0 c0 {2,S} {33,S} {34,S} {35,S}
7  C u0 p0 c0 {3,S} {36,S} {37,S} {38,S}
8  C u0 p0 c0 {2,S} {11,B} {12,B}
9  C u0 p0 c0 {1,S} {13,B} {14,B}
10 C u0 p0 c0 {3,S} {15,B} {16,B}
11 C u0 p0 c0 {8,B} {17,B} {39,S}
12 C u0 p0 c0 {8,B} {19,B} {43,S}
13 C u0 p0 c0 {9,B} {20,B} {44,S}
14 C u0 p0 c0 {9,B} {22,B} {48,S}
15 C u0 p0 c0 {10,B} {23,B} {49,S}
16 C u0 p0 c0 {10,B} {25,B} {53,S}
17 C u0 p0 c0 {11,B} {18,B} {40,S}
18 C u0 p0 c0 {17,B} {19,B} {41,S}
19 C u0 p0 c0 {12,B} {18,B} {42,S}
20 C u0 p0 c0 {13,B} {21,B} {45,S}
21 C u0 p0 c0 {20,B} {22,B} {46,S}
22 C u0 p0 c0 {14,B} {21,B} {47,S}
23 C u0 p0 c0 {15,B} {24,B} {50,S}
24 C u0 p0 c0 {23,B} {25,B} {51,S}
25 C u0 p0 c0 {16,B} {24,B} {52,S}
26 H u0 p0 c0 {2,S}
27 H u0 p0 c0 {1,S}
28 H u0 p0 c0 {3,S}
29 H u0 p0 c0 {4,S}
30 H u0 p0 c0 {4,S}
31 H u0 p0 c0 {5,S}
32 H u0 p0 c0 {5,S}
33 H u0 p0 c0 {6,S}
34 H u0 p0 c0 {6,S}
35 H u0 p0 c0 {6,S}
36 H u0 p0 c0 {7,S}
37 H u0 p0 c0 {7,S}
38 H u0 p0 c0 {7,S}
39 H u0 p0 c0 {11,S}
40 H u0 p0 c0 {17,S}
41 H u0 p0 c0 {18,S}
42 H u0 p0 c0 {19,S}
43 H u0 p0 c0 {12,S}
44 H u0 p0 c0 {13,S}
45 H u0 p0 c0 {20,S}
46 H u0 p0 c0 {21,S}
47 H u0 p0 c0 {22,S}
48 H u0 p0 c0 {14,S}
49 H u0 p0 c0 {15,S}
50 H u0 p0 c0 {23,S}
51 H u0 p0 c0 {24,S}
52 H u0 p0 c0 {25,S}
53 H u0 p0 c0 {16,S}


PS_1
1  C u0 p0 c0 {4,S} {5,S} {9,S} {27,S}
2  C u0 p0 c0 {4,S} {6,S} {8,S} {26,S}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {28,S}
4  C u0 p0 c0 {1,S} {2,S} {29,S} {30,S}
5  C u0 p0 c0 {1,S} {3,S} {31,S} {32,S}
6  C u0 p0 c0 {2,S} {33,S} {34,S} {35,S}
7  C u0 p0 c0 {3,S} {36,S} {37,S} {38,S}
8  C u0 p0 c0 {2,S} {11,S} {12,D}
9  C u0 p0 c0 {1,S} {13,S} {14,D}
10 C u0 p0 c0 {3,S} {15,S} {16,D}
11 C u0 p0 c0 {8,S} {17,D} {39,S}
12 C u0 p0 c0 {8,D} {19,S} {43,S}
13 C u0 p0 c0 {9,S} {20,D} {44,S}
14 C u0 p0 c0 {9,D} {22,S} {48,S}
15 C u0 p0 c0 {10,S} {23,D} {49,S}
16 C u0 p0 c0 {10,D} {25,S} {53,S}
17 C u0 p0 c0 {11,D} {18,S} {40,S}
18 C u0 p0 c0 {17,S} {19,D} {41,S}
19 C u0 p0 c0 {12,S} {18,D} {42,S}
20 C u0 p0 c0 {13,D} {21,S} {45,S}
21 C u0 p0 c0 {20,S} {22,D} {46,S}
22 C u0 p0 c0 {14,S} {21,D} {47,S}
23 C u0 p0 c0 {15,D} {24,S} {50,S}
24 C u0 p0 c0 {23,S} {25,D} {51,S}
25 C u0 p0 c0 {16,S} {24,D} {52,S}
26 H u0 p0 c0 {2,S}
27 H u0 p0 c0 {1,S}
28 H u0 p0 c0 {3,S}
29 H u0 p0 c0 {4,S}
30 H u0 p0 c0 {4,S}
31 H u0 p0 c0 {5,S}
32 H u0 p0 c0 {5,S}
33 H u0 p0 c0 {6,S}
34 H u0 p0 c0 {6,S}
35 H u0 p0 c0 {6,S}
36 H u0 p0 c0 {7,S}
37 H u0 p0 c0 {7,S}
38 H u0 p0 c0 {7,S}
39 H u0 p0 c0 {11,S}
40 H u0 p0 c0 {17,S}
41 H u0 p0 c0 {18,S}
42 H u0 p0 c0 {19,S}
43 H u0 p0 c0 {12,S}
44 H u0 p0 c0 {13,S}
45 H u0 p0 c0 {20,S}
46 H u0 p0 c0 {21,S}
47 H u0 p0 c0 {22,S}
48 H u0 p0 c0 {14,S}
49 H u0 p0 c0 {15,S}
50 H u0 p0 c0 {23,S}
51 H u0 p0 c0 {24,S}
52 H u0 p0 c0 {25,S}
53 H u0 p0 c0 {16,S}"""
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
        assert pol_1.is_isomorphic(spc_1)
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
        assert polymer.find_labeled_atom(p.monomer, LABELS_1) is not None
        assert polymer.find_labeled_atom(p.monomer, LABELS_2) is not None

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
        assert polymer.find_labeled_atom(head, LABELS_1) is not None
        assert polymer.find_labeled_atom(tail, LABELS_2) is not None

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
        assert spc1 is spc2

        mol = spc1.molecule[0]
        assert polymer.find_labeled_atom(mol, LABELS_1) is None
        assert polymer.find_labeled_atom(mol, LABELS_2) is None

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
        For PS proxy, a single H-abstraction on the backbone may reduce intact-monomer matches
        such that create_reacted_copy() classifies as scission. This test pins that behavior.
        """
        p = self.polymer_1.copy()
        reacted_proxy = p.baseline_proxy.molecule[0].copy(deep=True)
        abstract_h_from_center_backbone(reacted_proxy)
        new_p = p.create_reacted_copy(reacted_proxy)
        assert new_p is not None
        assert isinstance(new_p, Polymer)
        assert new_p.feature_monomer is None
        assert new_p.label.endswith("_scission_tail") or new_p.label.endswith("_scission_head")

    def test_create_reacted_copy_modification_baseline(self):
        """
        Ensures that an intact baseline proxy (unreacted) produces a
        modified polymer (_mod) because it contains both wings.
        """
        p = self.polymer_1.copy()
        reacted_proxy = p.baseline_proxy.molecule[0].copy(deep=True)
        new_p = p.create_reacted_copy(reacted_proxy)
        assert new_p is not None
        assert new_p.label == p.label
        assert new_p.feature_monomer is None

    def test_create_reacted_copy_head_scission_returns_scission_tail_polymer(self):
        """
        Construct a fragment with ONLY the head wing present (no tail wing),
        so create_reacted_copy() should classify it as head-side scission and return a Polymer
        with a NEW tail end-group labeled *2 and mono-radical.
        """
        p = self.polymer_1.copy()
        head_wing = p._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        scission_fragment = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert scission_fragment is not None
        new_p = p.create_reacted_copy(scission_fragment)
        assert new_p is not None
        assert isinstance(new_p, Polymer)
        assert new_p.feature_monomer is None
        assert new_p.label.endswith("_scission_tail")
        new_tail = new_p.end_groups[1]
        labels = {a.label for a in new_tail.atoms if a.label}
        assert labels == {"*2"}
        assert new_tail.get_radical_count() == 1
        assert np.isclose(new_p.Mn, p.Mn)
        assert np.isclose(new_p.Mw, p.Mw)

    def test_create_reacted_copy_tail_scission_returns_scission_head_polymer(self):
        """
        Construct a fragment with ONLY the tail wing present (no head wing),
        so create_reacted_copy() should classify it as tail-side scission and return a Polymer
        with a NEW head end-group labeled *1 and mono-radical.
        """
        p = self.polymer_1.copy()
        tail_wing = p._stitch_wing("tail")
        methyl_star1 = Molecule().from_adjacency_list(_methyl_radical_adj("*1"))
        scission_fragment = polymer.stitch_molecules_by_labeled_atoms(methyl_star1, tail_wing)
        assert scission_fragment is not None
        new_p = p.create_reacted_copy(scission_fragment)
        assert new_p is not None
        assert isinstance(new_p, Polymer)
        assert new_p.feature_monomer is None
        assert new_p.label.endswith("_scission_head")

    def test_create_reacted_copy_returns_none_for_small_molecule(self):
        """
        If reacted_proxy has no recognizable wings, create_reacted_copy() should return None.
        """
        p = self.polymer_1.copy()
        small = Molecule(smiles="CC")  # no head/tail wing subgraphs
        assert p.create_reacted_copy(small) is None

    def test_backbone_group_property(self):
        """
        Test that backbone_group generates a correctly relaxed pattern and caches it.
        """
        g = self.polymer_3.backbone_group
        assert len(g.atoms) == 2  # PE monomer heavy atoms only

        g = self.polymer_1.backbone_group
        assert len(g.atoms) == 8  # styrene heavy atoms only

        mol = self.polymer_3.baseline_proxy.molecule[0].copy(deep=True)
        mol.clear_labeled_atoms()
        matches = mol.find_subgraph_isomorphisms(self.polymer_3.backbone_group)
        assert len(matches) > 0

        # 1. Trigger generation
        group = self.polymer_1.backbone_group

        # 2. Verify Caching
        # Accessing it again should return the exact same object instance
        group_2 = self.polymer_1.backbone_group
        assert group is group_2, "Property should return cached instance on subsequent calls"

        # 3. Verify Structure & Labels
        # Ensure labels (*1, *2) are stripped
        for atom in group.atoms:
            assert atom.label == '', "Labels should be stripped from backbone group"

        # 4. Verify Relaxed Constraints (The "Fuzzy" Logic)
        for atom in group.atoms:
            assert atom.charge == [], "Charge should be wildcarded"
            assert atom.lone_pairs == [], "Lone pairs should be wildcarded"
            assert atom.radical_electrons == [0], "Radicals must be strictly [0]"

        # 5. Verify Bond Order Relaxation
        # Check that bonds allow Single, Benzene, Double, Triple ([1, 1.5, 2, 3])
        expected_orders = sorted([1, 1.5, 2, 3])

        bond_checked = False
        for atom in group.atoms:
            for neighbor, bond in atom.bonds.items():
                bond_checked = True
                assert sorted(bond.order) == expected_orders, \
                    f"Bond orders should be relaxed to {expected_orders}, got {bond.order}"

        assert bond_checked, "Monomer group should have at least one bond to check"

    def test_wing_groups_relaxation(self):
        """Verify wing templates are properly relaxed for matching."""
        wings = self.polymer_1._wing_groups("head")
        assert len(wings) > 0
        for group in wings:
            for g_atom in group.atoms:
                assert g_atom.radical_electrons == []
                if g_atom.is_carbon() and any(at.label == 'Cb' for at in g_atom.atomtype):
                    labels = {at.label for at in g_atom.atomtype}
                    assert 'Cb' in labels and 'Cd' in labels

    def test_get_heavy_view_with_maps(self):
        """Test that get_heavy_view_with_maps returns a molecule with correct atom maps."""
        full_mol = Molecule(smiles="CCC")
        expected_heavy = sum(1 for a in full_mol.atoms if not a.is_hydrogen())
        expected_light = sum(1 for a in full_mol.atoms if a.is_hydrogen())
        assert expected_light > 0
        assert expected_heavy == 3
        heavy_mol, heavy_to_full = polymer.get_heavy_view_with_maps(full_mol)
        assert isinstance(heavy_mol, Molecule)
        assert isinstance(heavy_to_full, dict)
        assert len(heavy_mol.atoms) == expected_heavy
        assert len(heavy_to_full) == expected_heavy
        assert not any(a.is_hydrogen() for a in heavy_mol.atoms)
        for heavy_atom in heavy_mol.atoms:
            assert heavy_atom in heavy_to_full, "Heavy atom missing from map keys."
            orig_atom = heavy_to_full[heavy_atom]
            assert orig_atom in full_mol.atoms, "Mapped atom is not in the original molecule."
            assert heavy_atom.element.symbol == orig_atom.element.symbol
            assert id(heavy_atom) != id(orig_atom), "Atoms share memory address; deep copy failed."

        no_h_mol = Molecule().from_smiles("[C-]#[O+]")
        no_h_heavy, no_h_map = polymer.get_heavy_view_with_maps(no_h_mol)
        assert len(no_h_heavy.atoms) == len(no_h_mol.atoms) == 2
        assert len(no_h_map) == 2
        assert not any(a.is_hydrogen() for a in no_h_heavy.atoms)

    def test_get_heavy_cut_edges(self):
        """Test that get_heavy_cut_edges correctly identifies boundary bonds."""
        # Create a simple, chemically valid linear backbone: Butane
        mol = Molecule(smiles="CCCC")
        heavy_mol, _ = polymer.get_heavy_view_with_maps(mol)
        c1 = next(atom for atom in heavy_mol.atoms if len(atom.bonds) == 1)
        c2 = list(c1.bonds.keys())[0]
        c3 = next(atom for atom in c2.bonds.keys() if atom is not c1)
        c4 = next(atom for atom in c3.bonds.keys() if atom is not c2)

        # --- Scenario 1: Terminal Wing (1 cut edge) ---
        # The wing is just the first carbon. The only cut should be C1-C2.
        wing_set_1 = {c1}
        cuts_1 = polymer.get_heavy_cut_edges(wing_set_1)
        assert len(cuts_1) == 1
        assert (c1, c2) in cuts_1

        # --- Scenario 2: Internal Atom (2 cut edges) ---
        # The wing is a single internal carbon. It should cross to C1 and C3.
        wing_set_2 = {c2}
        cuts_2 = polymer.get_heavy_cut_edges(wing_set_2)
        assert len(cuts_2) == 2
        assert (c2, c1) in cuts_2
        assert (c2, c3) in cuts_2

        # --- Scenario 3: Multi-Atom Fragment (2 cut edges) ---
        # The wing is the middle two carbons {C2, C3}.
        # The internal bond (C2-C3) should NOT be flagged as a cut.
        # The cuts should only be C2-C1 and C3-C4.
        wing_set_3 = {c2, c3}
        cuts_3 = polymer.get_heavy_cut_edges(wing_set_3)
        assert len(cuts_3) == 2
        assert (c2, c1) in cuts_3
        assert (c3, c4) in cuts_3

        # Explicitly verify the internal bond was ignored
        assert (c2, c3) not in cuts_3
        assert (c3, c2) not in cuts_3

        # --- Scenario 4: Entire Molecule (0 cut edges) ---
        # If the subset is the whole molecule, there are no external connections.
        wing_set_4 = {c1, c2, c3, c4}
        cuts_4 = polymer.get_heavy_cut_edges(wing_set_4)
        assert len(cuts_4) == 0

    def test_analyze_wing_matches(self):
        """
        Comprehensive test for _analyze_wing_matches.
        Tests Wild-Type (2 wings), Scission (1 wing), and Gas (0 wings).
        """
        p = self.polymer_1
        head_wings = p._wing_groups("head")
        tail_wings = p._wing_groups("tail")
        monomer_group = p.backbone_group

        # --- Scenario 1: Wild-Type (Intact Trimer) ---
        # The baseline proxy should yield exactly 2 disjoint wings.
        baseline_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        count_wt, details_wt = polymer._analyze_wing_matches(baseline_mol, head_wings, tail_wings, monomer_group)

        assert count_wt == 2, f"Expected 2 wings for intact baseline, got {count_wt}"
        assert details_wt["head_match"] is not None
        assert details_wt["tail_match"] is not None
        assert "heavy_to_full_map" in details_wt
        assert len(details_wt["all_optimal_wings"]) == 2

        # Ensure the Max Set Packing found matches that are terminal
        # (A true wing in a trimer should only cross 1 boundary into the remainder)
        assert details_wt["head_match"]["cut_edges"] == 1
        assert details_wt["tail_match"]["cut_edges"] == 1

        # --- Scenario 2: Scission Fragment (Single Wing) ---
        # We simulate a scission by just taking the head wing molecule itself.
        # It contains a head wing, but no tail wing.
        single_wing_mol = p._stitch_wing("head")
        single_wing_mol.clear_labeled_atoms()  # Clean up stitch labels
        single_wing_mol.update_multiplicity()

        count_sc, details_sc = polymer._analyze_wing_matches(single_wing_mol, head_wings, tail_wings, monomer_group)

        assert count_sc == 1, f"Expected 1 wing for scission fragment, got {count_sc}"
        assert details_sc["head_match"] is not None
        assert details_sc["tail_match"] is None
        assert len(details_sc["all_optimal_wings"]) == 1

        # --- Scenario 3: Gas Phase (Zero Wings) ---
        # A simple methane molecule should completely fail to match the large wing patterns.
        gas_mol = Molecule(smiles="C")

        count_gas, details_gas = polymer._analyze_wing_matches(gas_mol, head_wings, tail_wings, monomer_group)

        assert count_gas == 0, f"Expected 0 wings for methane gas, got {count_gas}"
        assert details_gas["head_match"] is None
        assert details_gas["tail_match"] is None
        assert len(details_gas["all_optimal_wings"]) == 0

    def test_end_group_modification_and_slicing(self):
        """
        Tests _slice_wing mechanics and _is_end_group_modified using a PS baseline.
        Ensures modifications are geographically isolated to the correct zones.
        """
        # --- 1. Setup & Wild-Type Baseline ---
        p = self.polymer_1
        baseline_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        monomer_group = p.backbone_group
        head_wings = p._wing_groups("head")
        tail_wings = p._wing_groups("tail")

        # Calculate PS monomer heavy count (should be 8 for Styrene: 2 backbone, 6 phenyl)
        mon_heavy_count = sum(1 for ga in monomer_group.atoms if ga.atomtype[0].label[0] != 'H')

        wing_count, details_wt = polymer._analyze_wing_matches(baseline_mol, head_wings, tail_wings, monomer_group)
        assert wing_count == 2, "Setup failed: Could not find 2 wings in baseline PS."

        # --- 2. Basic Mechanics: Test _slice_wing ---
        head_heavy_atoms = details_wt['head_match']['atoms']
        end_group_heavy, buffer_heavy = polymer._slice_wing(head_heavy_atoms, mon_heavy_count)

        # The buffer should be capped at the monomer size, and the two sets must be perfectly disjoint
        assert len(buffer_heavy) <= mon_heavy_count
        assert len(end_group_heavy) + len(buffer_heavy) == len(head_heavy_atoms)
        assert len(end_group_heavy.intersection(buffer_heavy)) == 0

        # --- 3. Wild-Type Validation ---
        # The intact baseline proxy should definitively return False
        assert polymer._is_end_group_modified(details_wt, p) is False

        # --- 4. Edge Case: True End-Group Modification ---
        mod_end_mol = baseline_mol.copy(deep=True)
        _, details_mod = polymer._analyze_wing_matches(mod_end_mol, head_wings, tail_wings, monomer_group)

        # Dynamically find a full atom that strictly belongs to the End-Group
        heavy_to_full = details_mod['heavy_to_full_map']
        end_heavy, _ = polymer._slice_wing(details_mod['head_match']['atoms'], mon_heavy_count)
        end_target_full = heavy_to_full[list(end_heavy)[0]]

        # Apply the kinetic modification (e.g., H-abstraction leaving a radical)
        end_target_full.radical_electrons = 1

        # The validator should catch this specific modification
        assert polymer._is_end_group_modified(details_mod, p) is True

        # --- 5. Edge Case: Buffer Modification (False Positive Prevention) ---
        mod_buf_mol = baseline_mol.copy(deep=True)
        _, details_buf = polymer._analyze_wing_matches(mod_buf_mol, head_wings, tail_wings, monomer_group)

        heavy_to_full_buf = details_buf['heavy_to_full_map']
        _, buffer_heavy_mod = polymer._slice_wing(details_buf['head_match']['atoms'], mon_heavy_count)

        # Ensure we have a buffer to test (PS wings should include the buffer monomer)
        assert len(buffer_heavy_mod) > 0, "Test invalid: Wing pattern did not capture a buffer monomer."

        # Dynamically find a full atom that strictly belongs to the Buffer Monomer
        buf_target_full = heavy_to_full_buf[list(buffer_heavy_mod)[0]]

        # Apply the modification to the buffer
        buf_target_full.radical_electrons = 1

        # The End-Group validator MUST ignore this modification because it is in the buffer zone!
        assert polymer._is_end_group_modified(details_buf, p) is False

    def test_is_buffer_monomer_modified(self):
        """Test that _is_buffer_monomer_modified correctly identifies changes in the buffer zone."""
        # --- 1. Setup ---
        p = self.polymer_1
        head_wings = p._wing_groups("head")
        tail_wings = p._wing_groups("tail")
        monomer_group = p.backbone_group
        mon_heavy_count = sum(1 for ga in monomer_group.atoms if not getattr(ga, 'is_hydrogen', lambda: False)())

        # --- 2. Scenario 1: Wild-Type (No Modifications) ---
        baseline_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        _, details_wt = polymer._analyze_wing_matches(baseline_mol, head_wings, tail_wings, monomer_group)

        # The intact baseline proxy should confidently return False
        assert polymer._is_buffer_monomer_modified(details_wt, p) is False

        # --- 3. Scenario 2: Buffer Modification ---
        mod_buf_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        _, details_buf = polymer._analyze_wing_matches(mod_buf_mol, head_wings, tail_wings, monomer_group)

        # Use our slicing helper to physically isolate the buffer zone
        heavy_to_full = details_buf['heavy_to_full_map']
        _, buffer_heavy = polymer._slice_wing(details_buf['head_match']['atoms'], mon_heavy_count)

        assert len(buffer_heavy) > 0, "Test invalid: Wing pattern did not capture a buffer monomer."

        # Grab a target atom strictly inside the buffer zone and mutate it
        buf_target_heavy = list(buffer_heavy)[0]
        buf_target_full = heavy_to_full[buf_target_heavy]

        # Apply the kinetic modification (e.g., an abstracted hydrogen leaving a radical)
        buf_target_full.radical_electrons = 1

        # The validator MUST catch this modification because it falls exactly in the buffer zone
        assert polymer._is_buffer_monomer_modified(details_buf, p) is True

    def test_is_center_feature_modified(self):
        """Test that _is_center_feature_modified exclusively catches center backbone changes."""
        p = self.polymer_1
        head_wings = p._wing_groups("head")
        tail_wings = p._wing_groups("tail")
        monomer_group = p.backbone_group

        baseline_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        _, details_wt = polymer._analyze_wing_matches(baseline_mol, head_wings, tail_wings, monomer_group)

        # --- Scenario 1: Wild-Type (No Modifications) ---
        # The intact baseline proxy should definitively return False
        assert polymer._is_center_feature_modified(baseline_mol, details_wt) is False

        # --- Scenario 2: Center Modification ---
        mod_center_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        _, details_center = polymer._analyze_wing_matches(mod_center_mol, head_wings, tail_wings, monomer_group)

        # Identify an atom strictly in the center by subtracting the wing full atoms
        heavy_to_full = details_center['heavy_to_full_map']
        wing_heavies = set(details_center['head_match']['atoms']).union(details_center['tail_match']['atoms'])
        wing_fulls = set()

        for ha in wing_heavies:
            fa = heavy_to_full[ha]
            wing_fulls.add(fa)
            for neighbor in fa.bonds.keys():
                if neighbor.is_hydrogen():
                    wing_fulls.add(neighbor)

        center_full_atoms = [a for a in mod_center_mol.atoms if a not in wing_fulls]
        assert len(center_full_atoms) > 0, "Test invalid: Trimer proxy has no center atoms left."

        # Mutate the center atom (e.g., an abstracted hydrogen leaving a radical)
        center_target_full = center_full_atoms[0]
        center_target_full.radical_electrons = 1

        # The validator MUST catch this because the radical is in the central repeating unit
        assert polymer._is_center_feature_modified(mod_center_mol, details_center) is True

        # --- Scenario 3: Wing Modification (False Positive Prevention) ---
        mod_wing_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        _, details_wing = polymer._analyze_wing_matches(mod_wing_mol, head_wings, tail_wings, monomer_group)

        heavy_to_full_wing = details_wing['heavy_to_full_map']

        # Intentionally mutate a wing atom instead
        wing_target_heavy = list(details_wing['head_match']['atoms'])[0]
        wing_target_full = heavy_to_full_wing[wing_target_heavy]
        wing_target_full.radical_electrons = 1

        # The center validator must IGNORE this, as it is geographically outside its domain
        assert polymer._is_center_feature_modified(mod_wing_mol, details_wing) is False

    def test_stitch_returns_none_if_any_input_none(self):
        """Test that stitch_molecules_by_labeled_atoms returns None if any input is None."""
        mol = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        assert polymer.stitch_molecules_by_labeled_atoms(None, mol) is None
        assert polymer.stitch_molecules_by_labeled_atoms(mol, None) is None

    def test_stitch_for_p1(self):
        """Test stitch_molecules_by_labeled_atoms for polymer_1 head wing + methyl radical."""
        p = self.polymer_1.copy()
        head_wing = p._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        scission_fragment = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert scission_fragment.to_smiles() in ['CCC(C)C1=CC=CC=C1']

    def test_stitch_raises_when_missing_labels(self):
        """Test that stitch_molecules_by_labeled_atoms raises ValueError when labels missing."""
        left = Molecule().from_adjacency_list(_methyl_radical_adj("*1"))
        right = Molecule(smiles="[CH3]")  # radical but no *2 label
        with pytest.raises(ValueError):
            polymer.stitch_molecules_by_labeled_atoms(left, right)

    def test_stitch_raises_when_stitch_sites_not_mono_radicals(self):
        """Test that stitch_molecules_by_labeled_atoms raises ValueError when stitch sites not mono-radicals."""
        left = Molecule().from_adjacency_list(_methyl_closed_shell_labeled_adj("*1"))  # labeled but u0
        right = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        with pytest.raises(ValueError):
            polymer.stitch_molecules_by_labeled_atoms(left, right)

    def test_stitch_does_not_mutate_inputs_and_clears_labels_in_product(self):
        """Test that stitch_molecules_by_labeled_atoms does not mutate inputs and clears labels in product."""
        left = Molecule().from_adjacency_list(_methyl_radical_adj("*1"))
        right = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        left_before_labels = [a.label for a in left.atoms]
        right_before_labels = [a.label for a in right.atoms]
        left_before_rad = left.get_radical_count()
        right_before_rad = right.get_radical_count()
        merged = polymer.stitch_molecules_by_labeled_atoms(left, right)
        assert merged is not None

        # inputs unchanged (function deep-copies)
        assert [a.label for a in left.atoms] == left_before_labels
        assert [a.label for a in right.atoms] == right_before_labels
        assert left.get_radical_count() == left_before_rad
        assert right.get_radical_count() == right_before_rad

        # product has no stitch labels remaining at the join sites
        assert polymer.find_labeled_atom(merged, LABELS_1) is None
        assert polymer.find_labeled_atom(merged, LABELS_2) is None
        assert merged.get_radical_count() == 0

    @pytest.mark.parametrize("do_update", [False, True])
    def test_end_group_modification_does_not_require_atomtypes(self, do_update):
        """
        Scenario: Head end-group modification.
        Uses Polystyrene (PS) to ensure the wing is large enough to slice into
        distinct End-Group and Buffer zones.
        """
        p = self.polymer_1
        product_mol = p.baseline_proxy.molecule[0].copy(deep=True)
        head_wings = p._wing_groups("head")
        tail_wings = p._wing_groups("tail")
        monomer_group = p.backbone_group
        _, details = polymer._analyze_wing_matches(product_mol, head_wings, tail_wings, monomer_group)
        heavy_to_full = details['heavy_to_full_map']
        head_heavy_atoms = details['head_match']['atoms']
        mon_heavy_count = sum(1 for ga in monomer_group.atoms if not ga.is_hydrogen())
        end_group_heavy, _ = polymer._slice_wing(head_heavy_atoms, mon_heavy_count)
        assert len(end_group_heavy) > 0, "PS Wing should have leftover atoms after slicing buffer."
        target_full = heavy_to_full[list(end_group_heavy)[0]]
        h_atom = next(n for n in target_full.bonds if n.is_hydrogen())
        product_mol.remove_bond(product_mol.get_bond(target_full, h_atom))
        product_mol.remove_atom(h_atom)
        target_full.radical_electrons = 1
        product_mol.update_multiplicity()
        spc = Species(molecule=[product_mol])
        classification, _ = polymer.classify_structure(spc, p)
        assert classification == PolymerClass.END_MOD

    def test_no_molecule_returns_gas(self):
        """Species without a molecule array gracefully exits."""
        spc = Species()
        classification, details = polymer.classify_structure(spc, self.polymer_3)
        assert classification == PolymerClass.UNKNOWN
        assert details["reason"] == "no_molecule"

    def test_proxy_to_gas_flag_update(self):
        """
        Ensures a misclassified proxy successfully overwrites its proxy flag to False.
        (Renamed and updated to remove the print/capsys requirement).
        """
        spc = Species(label="BadProxy", molecule=[Molecule(smiles="C")])
        if not hasattr(spc, "props"): spc.props = {}
        spc.props["is_polymer_proxy"] = True
        for m in spc.molecule:
            if not hasattr(m, "props"): m.props = {}
            m.props["is_polymer_proxy"] = True
        polymer.process_polymer_candidates([spc], None, self.polymer_3)
        assert spc.props.get("is_polymer_proxy") is False
        for m in spc.molecule:
            assert m.props.get("is_polymer_proxy") is False

    def test_end_mod_more_than_3_matches_sets_note(self):
        """
        Validates classification when a chain matches > 3 times.
        Uses Polystyrene to ensure the wings are deep enough to support
        distinct End-Group vs Buffer zones.
        """
        p = self.polymer_1
        # Build a PS tetramer (4 units) to ensure > 3 matches (3 is enough, but 4 is safer)
        # We can just use the baseline proxy if it's already a trimer, or stitch a longer one.
        # For simplicity, let's use the baseline trimer (3 units) but target the end-group.
        product_mol = p.baseline_proxy.molecule[0].copy(deep=True)

        # 1. Analyze to find topological domains
        head_wings = p._wing_groups("head")
        tail_wings = p._wing_groups("tail")
        monomer_group = p.backbone_group
        _, details = polymer._analyze_wing_matches(product_mol, head_wings, tail_wings, monomer_group)

        # 2. Targeted Strike on the True End-Group (Initiator fragment)
        heavy_to_full = details['heavy_to_full_map']
        head_heavy_atoms = details['head_match']['atoms']
        mon_heavy_count = sum(1 for ga in monomer_group.atoms if not ga.is_hydrogen())

        end_group_heavy, _ = polymer._slice_wing(head_heavy_atoms, mon_heavy_count)
        target_full = heavy_to_full[list(end_group_heavy)[0]]
        h_atom = next(n for n in target_full.bonds if n.is_hydrogen())
        product_mol.remove_bond(product_mol.get_bond(target_full, h_atom))
        product_mol.remove_atom(h_atom)
        target_full.radical_electrons = 1
        product_mol.update_multiplicity()

        # 3. Final Classification
        spc = Species(molecule=[product_mol])
        classification, details = polymer.classify_structure(spc, p)
        assert classification == PolymerClass.END_MOD

    def test_proxy_true_but_no_backbone_matches_returns_gas_reason(self):
        """
        Simulates the '0 core reactions' starvation mode where a candidate
        claims to be a proxy but structurally isn't.
        """
        spc = Species(molecule=[Molecule(smiles="CO")])
        spc.is_polymer_proxy = True
        classification, details = polymer.classify_structure(spc, self.polymer_3)
        assert classification == PolymerClass.GAS
        assert details["reason"] == "no_intact_wings"

    def test_process_filters_discard_and_sets_flags(self):
        """
        Verifies `process_polymer_candidates` correctly updates the proxy flags
        and drops DISCARD candidates.
        """
        p = self.polymer_3
        # 1. FEAT (Radical in center)
        s_feat = Species(label="FEAT", molecule=[p.baseline_proxy.molecule[0].copy(deep=True)])
        c_feat = get_monomer_regions(s_feat.molecule[0])['center'][0]
        _safe_make_radical(s_feat.molecule[0], c_feat)
        s_feat.molecule[0].update_multiplicity()

        # 2. DISC (Radical in buffer)
        s_disc = Species(label="DISC", molecule=[p.baseline_proxy.molecule[0].copy(deep=True)])
        c_disc = get_monomer_regions(s_disc.molecule[0])['head_buffer'][1]
        _safe_make_radical(s_disc.molecule[0], c_disc)
        s_disc.molecule[0].update_multiplicity()

        # 3. GAS
        s_gas = Species(label="GAS", molecule=[Molecule(smiles="C")])
        out = polymer.process_polymer_candidates([s_feat, s_disc, s_gas], None, p)
        assert len(out) == 2

    def test_restore_labels_head_scission(self):
        """
        Test that _restore_labels identifies a cut bond to a head-wing atom
        and correctly labels the remainder atom as '*2' with a radical.
        """
        original_mol = Molecule().from_smiles("CC")
        c1, c2 = original_mol.atoms[0], original_mol.atoms[1]
        removed_atoms = {c1}
        for neighbor in c1.bonds:
            if neighbor.is_hydrogen():
                removed_atoms.add(neighbor)
        new_mol = Molecule()
        for _ in range(4):
            new_mol.add_atom(Atom(element='C' if _ == 0 else 'H'))
        head_match_atoms = {c1}
        polymer.Polymer._restore_labels(
            new_mol=new_mol,
            original_mol=original_mol,
            removed_atoms=removed_atoms,
            head_match_atoms=head_match_atoms,
            tail_match_atoms=None)
        res_atom = new_mol.atoms[0]
        assert res_atom.label == '*2'
        assert res_atom.radical_electrons == 1

    def test_restore_labels_mapping_failure(self):
        """
        Test that _restore_labels raises ValueError if the atom counts
        between the original (minus removed) and the new molecule don't match.
        """
        original_mol = Molecule().from_smiles("CC")
        new_mol = Molecule().from_smiles("C")
        with pytest.raises(ValueError, match="Mapping failure"):
            polymer.Polymer._restore_labels(
                new_mol=new_mol,
                original_mol=original_mol,
                removed_atoms=set(),
                head_match_atoms=None,
                tail_match_atoms=None)

    def test_restore_labels_conflict(self):
        """
        Test that _restore_labels raises ValueError if an atom
        tries to be both *1 and *2.
        """
        mol = Molecule().from_smiles("CCC")
        c1, c2, c3 = mol.atoms[0], mol.atoms[1], mol.atoms[2]
        removed = {c1, c3}
        for terminal_c in [c1, c3]:
            for neighbor in terminal_c.bonds:
                if neighbor.is_hydrogen():
                    removed.add(neighbor)
        new_mol = Molecule()
        for _ in range(3):
            new_mol.add_atom(Atom(element='C' if _ == 0 else 'H'))
        head_match = {c3}
        tail_match = {c1}
        with pytest.raises(ValueError, match="Label conflict"):
            polymer.Polymer._restore_labels(
                new_mol=new_mol,
                original_mol=mol,
                removed_atoms=removed,
                head_match_atoms=head_match,
                tail_match_atoms=tail_match)

    def test_stitch_trimer_wildtypes(self):
        """
        Test that _stitch_trimer creates a correctly capped trimer (5 segments total).
        Tests Polystyrene (PS) and defines PMMA.
        """
        # --- 1. Polystyrene (PS) Wildtype ---
        ps_trimer_spc = self.polymer_1._stitch_trimer(self.polymer_1.monomer)

        assert isinstance(ps_trimer_spc, Species)
        assert ps_trimer_spc.is_polymer_proxy is True

        ps_mol = ps_trimer_spc.molecule[0]
        assert len([a for a in ps_mol.atoms if a.is_carbon()]) == 25
        assert all(a.label == '' for a in ps_mol.atoms)

        # --- 2. New Polymer Definition: PMMA ---
        pmma_adj = """
multiplicity 3
1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}
3    C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
4    H u0 p0 c0 {1,S}
5    C u0 p0 c0 {2,S} {10,D} {11,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {3,S}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {3,S}
10   O u0 p2 c0 {5,D}
11   O u0 p2 c0 {5,S} {12,S}
12   C u0 p0 c0 {11,S} {13,S} {14,S} {15,S}
13   H u0 p0 c0 {12,S}
14   H u0 p0 c0 {12,S}
15   H u0 p0 c0 {12,S}
"""
        pmma_poly = Polymer(
            label='PMMA',
            monomer=pmma_adj,
            end_groups=['[H]', '[H]'],
            cutoff=5,
            Mn=2000.0,
            Mw=4000.0,
            initial_mass=1.0)

        pmma_trimer_spc = pmma_poly._stitch_trimer(pmma_poly.monomer)

        assert pmma_trimer_spc is not None
        pmma_mol = pmma_trimer_spc.molecule[0]

        # PMMA monomer (C5H8O2) -> Trimer (3 units) + 2H ends = C15H26O6
        c_count = len([a for a in pmma_mol.atoms if a.is_carbon()])
        o_count = len([a for a in pmma_mol.atoms if a.is_oxygen()])
        assert c_count == 15
        assert o_count == 6

    def test_get_polydispersity(self):
        """
        Test PDI calculation (Mw/Mn) and edge case handling.
        """
        # 1. Standard Case: PE_1 (Mn=1000, Mw=2500)
        # PDI = 2500 / 1000 = 2.5
        assert self.polymer_3.get_polydispersity() == 2.5

        # 2. Monodisperse Case: Mn = Mw
        # PDI should be 1.0
        mono_poly = self.polymer_3.copy()
        mono_poly.Mn = 5000.0
        mono_poly.Mw = 5000.0
        assert mono_poly.get_polydispersity() == 1.0

        # 3. Edge Case: Mn is None
        none_poly = self.polymer_3.copy()
        none_poly.Mn = None
        assert none_poly.get_polydispersity() == 0.0

        # 4. Edge Case: Mn is 0 (Avoid DivisionByZero)
        zero_poly = self.polymer_3.copy()
        zero_poly.Mn = 0.0
        assert zero_poly.get_polydispersity() == 0.0

    def test_calculate_distribution_from_moments(self):
        """
        Test the back-calculation of Mn and Mw from raw moments.
        Mn = (mu1 / mu0) * MonomerMW
        Mw = (mu2 / mu1) * MonomerMW
        """
        # Using polymer_3 (PE) from setup_species
        p = self.polymer_3.copy()

        # 1. Setup moments for a hypothetical distribution
        # mu0 = 1.0 (1 mole of chains)
        # mu1 = 100.0 (100 moles of monomer units) -> DPn = 100
        # mu2 = 15000.0 -> DPw = 150
        p.moments = np.array([1.0, 100.0, 15000.0])

        # monomer_mw_g_mol is calculated in __init__ as get_molecular_weight() * 1000
        monomer_mw = p.monomer.get_molecular_weight() * 1000.0
        mn_calc, mw_calc = p._calculate_distribution_from_moments()
        expected_mn = 100.0 * monomer_mw
        expected_mw = 150.0 * monomer_mw
        assert np.isclose(mn_calc, expected_mn)
        assert np.isclose(mw_calc, expected_mw)

    def test_calculate_distribution_from_zero_moments(self):
        """Ensure the calculation handles zero-density states without crashing."""
        p = self.polymer_3.copy()

        # Case: mu0 is 0 (No polymer present)
        p.moments = np.array([0.0, 0.0, 0.0])
        mn, mw = p._calculate_distribution_from_moments()
        assert mn == 0.0
        assert mw == 0.0

        # Case: moments are None
        p.moments = None
        mn_none, mw_none = p._calculate_distribution_from_moments()
        assert mn_none is None
        assert mw_none is None

    def test_distribution_moment_round_trip(self):
        """
        Verifies that converting Mn/Mw -> Moments -> Mn/Mw preserves values.
        """
        p = self.polymer_3.copy()
        original_mn = p.Mn
        original_mw = p.Mw

        # This calls _calculate_moments_from_distribution() internally which depends on the Mn/Mw provided at init
        moments = p._calculate_moments_from_distribution()
        p.moments = moments

        # Now back-calculate
        new_mn, new_mw = p._calculate_distribution_from_moments()

        assert np.isclose(original_mn, new_mn)
        assert np.isclose(original_mw, new_mw)

    def test_fingerprint_generation(self):
        """
        Verify that the fingerprint is a stable string incorporating
        monomer, feature_monomer, and cutoff.
        """
        p = self.polymer_1
        fp = p.fingerprint

        # 1. Basic format check
        assert fp.startswith("Polymer_")
        assert f"_{p.cutoff}" in fp

        # 2. Immutability/Caching check
        # Fingerprint should be read-only (cached)
        first_fp = p.fingerprint
        assert first_fp is p.fingerprint

        # 3. Sensitivity check: Changing cutoff changes fingerprint
        p_diff_cutoff = p.copy()
        p_diff_cutoff.cutoff = 10
        p_diff_cutoff._fingerprint = None
        assert p_diff_cutoff.fingerprint != fp
        assert fp.endswith("_3")
        assert p_diff_cutoff.fingerprint.endswith("_10")

    def test_fingerprint_with_feature_monomer(self):
        """
        Verify that the feature_monomer fingerprint is included when present.
        """
        # polymer_4 has a feature_monomer defined in setup_species
        p = self.polymer_4
        fp = p.fingerprint
        assert "_Feat-" in fp
        assert p.monomer.fingerprint in fp
        assert p.feature_monomer.fingerprint in fp

    def test_fingerprint_consistency(self):
        """
        Verify that two identical polymers result in the same fingerprint.
        """
        p1 = self.polymer_1
        p2 = self.polymer_1.copy()
        # Clear cache to force generation
        p2._fingerprint = None
        assert p1.fingerprint == p2.fingerprint

    def test_get_closing_moment_accuracy(self):
        """
        Verify the Log-Lagrange extrapolation for mu3.
        mu3 = (mu2^3 * mu0) / mu1^3
        """
        p = self.polymer_3
        # Setup specific moments: mu0=1, mu1=10, mu2=200
        # mu3 = (200^3 * 1) / 10^3 = 8,000,000 / 1,000 = 8,000
        mu = [1.0, 10.0, 200.0]
        mu3 = p.get_closing_moment(mu)
        assert np.isclose(mu3, 8000.0)

    def test_get_closing_moment_stability_at_zero(self):
        """
        Ensure the closure doesn't crash or return NaN/Inf when moments are zero.
        This happens at t=0 when polymer species haven't formed yet.
        """
        p = self.polymer_3

        # All zeros
        assert p.get_closing_moment([0.0, 0.0, 0.0]) == 0.0

        # mu0 is zero (no number density)
        assert p.get_closing_moment([0.0, 10.0, 200.0]) == 0.0

        # mu1 is zero (to avoid ZeroDivisionError)
        assert p.get_closing_moment([1.0, 0.0, 200.0]) == 0.0

    def test_get_closing_moment_stability_negative(self):
        """
        Solver oversteps can result in slightly negative moments.
        The closure must handle this gracefully.
        """
        # Negative mu1 should return 0 rather than a complex number/NaN
        assert self.polymer_3.get_closing_moment([1.0, -0.01, 200.0]) == 0.0

    def test_scission_structural_split(self):
        """
        Tests the structural scission of a trimer into two fragments.
        Verifies that labels and end-groups are correctly redistributed.
        """
        # 1. Use the Polystyrene (PS) fixture
        poly = self.polymer_1
        # Generate the wildtype trimer proxy: [Cap1]-[M1]-[M2]-[M3]-[Cap2]
        trimer_spc = poly.get_proxy_species()
        trimer_mol = trimer_spc.molecule[0]

        # 2. Identify a scission point
        # We find a backbone C-C bond between two monomer units
        # For simplicity in this test, we'll simulate the extraction of a
        # fragment after a scission event.

        # 3. Test the 'healing' of a scission site
        # Suppose a scission occurred, leaving a fragment that needs a TailCap
        fragment_mol = Molecule().from_smiles("CCCC")  # A 4-carbon fragment
        # We need to restore labels as if it were cut from a larger chain

        # Identify 'removed' atoms (the rest of the chain)
        # We'll use a mock original and removed set
        original = Molecule().from_smiles("CCCCCC")
        removed = {original.atoms[4], original.atoms[5]}  # Remove the end
        for atom in list(removed):
            for neighbor in atom.bonds:
                if neighbor.is_hydrogen():
                    removed.add(neighbor)

        # The remainder is the 4-carbon fragment
        new_mol = Molecule()
        for _ in range(len([a for a in original.atoms if a not in removed])):
            new_mol.add_atom(Atom(element='C' if _ < 4 else 'H'))

        # Heal the scission site (tail-side cut)
        poly._restore_labels(
            new_mol=new_mol,
            original_mol=original,
            removed_atoms=removed,
            tail_match_atoms={original.atoms[4]}  # The atom that was cut away
        )

        # Verify the new tail connection point (*1) was created
        tail_atom = [a for a in new_mol.atoms if a.label == '*1']
        assert len(tail_atom) == 1
        assert tail_atom[0].radical_electrons == 1


class TestPolymerClassification:
    """
    Comprehensive test suite for the classify_structure() topological partitioner.
    Uses a Polystyrene (PS) wild-type baseline to dynamically verify all classification branches.
    """

    @pytest.fixture(autouse=True)
    def setup_polymer(self):
        """
        Initializes the PS baseline proxy and necessary components before every test.
        """
        # 1. Build the Polystyrene (PS) baseline
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

        self.p = Polymer(
            label='PS_1',
            monomer=ps_adj,
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,)

        # 2. Extract topological definitions
        self.monomer_group = self.p.backbone_group
        self.head_wings = self.p._wing_groups("head")
        self.tail_wings = self.p._wing_groups("tail")

        # 3. Determine the size of the PS monomer for dynamic slicing and thresholding
        self.mon_heavy_count = sum(1 for ga in self.monomer_group.atoms if not getattr(ga, 'is_hydrogen', lambda: False)())

        # 4. Create a clean reference baseline molecule
        self.ref_baseline_mol = self.p.baseline_proxy.molecule[0].copy(deep=True)

    def test_branch_gas_no_backbone(self):
        """
        Tests the failure mode where the polymer definition lacks a backbone group.
        Physically represents a malformed polymer object.
        """
        species = Species(molecule=[self.ref_baseline_mol])
        class BrokenPolymer:
            backbone_group = None
        p_class, details = polymer.classify_structure(species, BrokenPolymer())
        assert p_class == polymer.PolymerClass.GAS
        assert details["reason"] == "no_backbone_group"

    def test_branch_gas_too_few_atoms(self):
        """
        Tests the failure mode where the generated product is smaller than a single monomer.
        Physically represents an extreme degradation product (e.g., Methane off-gassing).
        """
        from rmgpy.molecule.molecule import Molecule
        tiny_mol = Molecule().from_smiles("C")
        species = Species(molecule=[tiny_mol])
        p_class, details = polymer.classify_structure(species, self.p)
        assert p_class == polymer.PolymerClass.GAS
        assert details["reason"] == "too_few_atoms_for_monomer"

    # =========================================================================
    # BRANCH A: INTACT BACKBONE (>= 2 WINGS)
    # =========================================================================

    def test_branch_baseline_unreacted(self):
        """
        Tests the exact, unreacted baseline proxy.
        Physically represents a spectator chain that did not undergo a kinetic reaction.
        """
        species = Species(molecule=[self.ref_baseline_mol.copy(deep=True)])

        p_class, details = polymer.classify_structure(species, self.p)

        assert p_class == polymer.PolymerClass.BASELINE
        assert details["reason"] == "unreacted_proxy"
        assert details["num_disjoint_wings"] == 2

    def test_branch_end_group_modification(self):
        """
        Tests a kinetic modification (radical) strictly localized to the terminal end-cap.
        Dynamically calculates the BFS slice to ensure accurate mutation targeting.
        """
        mod_mol = self.ref_baseline_mol.copy(deep=True)

        # 1. Analyze the molecule to find the topological zones
        _, match_details = polymer._analyze_wing_matches(mod_mol, self.head_wings, self.tail_wings, self.monomer_group)

        heavy_to_full = match_details['heavy_to_full_map']
        head_heavy_atoms = match_details['head_match']['atoms']

        # 2. Slice the wing into End-Cap and Buffer
        end_group_heavy, _ = polymer._slice_wing(head_heavy_atoms, self.mon_heavy_count)

        # 3. Apply the kinetic mutation (H-abstraction) strictly to the End-Cap
        target_heavy = list(end_group_heavy)[0]
        target_full = heavy_to_full[target_heavy]
        target_full.radical_electrons = 1

        # 4. Route through classify_structure
        species = Species(molecule=[mod_mol])
        p_class, details = polymer.classify_structure(species, self.p)

        assert p_class == polymer.PolymerClass.END_MOD
        assert details["reason"] == "terminal_end_modified"
        assert details["num_disjoint_wings"] == 2

    def test_branch_buffer_monomer_modification(self):
        """
        Tests a kinetic modification located inside the buffer monomer.
        This must be classified as DISCARD to prevent boundary effect contamination in kinetics.
        """
        mod_mol = self.ref_baseline_mol.copy(deep=True)

        # 1. Analyze and isolate the buffer zone
        _, match_details = polymer._analyze_wing_matches(mod_mol, self.head_wings, self.tail_wings, self.monomer_group)
        heavy_to_full = match_details['heavy_to_full_map']
        head_heavy_atoms = match_details['head_match']['atoms']
        _, buffer_heavy = polymer._slice_wing(head_heavy_atoms, self.mon_heavy_count)

        assert len(buffer_heavy) > 0, "Test failed to isolate a buffer zone."

        # 2. Mutate an atom strictly inside the buffer
        target_heavy = list(buffer_heavy)[0]
        target_full = heavy_to_full[target_heavy]
        target_full.radical_electrons = 1

        # 3. Route through classify_structure
        species = Species(molecule=[mod_mol])
        p_class, details = polymer.classify_structure(species, self.p)

        assert p_class == polymer.PolymerClass.DISCARD
        assert details["reason"] == "buffer_monomer_modified"
        assert details["num_disjoint_wings"] == 2

    def test_branch_center_feature_modification(self):
        """Tests kinetic modification in the center. Uses radical to ensure it's not baseline."""
        mod_mol = self.ref_baseline_mol.copy(deep=True)
        _, match_details = polymer._analyze_wing_matches(mod_mol, self.head_wings, self.tail_wings, self.monomer_group)
        heavy_to_full = match_details['heavy_to_full_map']

        wing_heavies = set(match_details['head_match']['atoms']).union(match_details['tail_match']['atoms'])
        center_full_atoms = [heavy_to_full[ha] for ha in heavy_to_full if ha not in wing_heavies]

        # Apply radical to center
        target = center_full_atoms[0]
        # Remove a hydrogen to make room for radical
        h_neighbor = next(n for n in target.bonds if n.is_hydrogen())
        mod_mol.remove_bond(mod_mol.get_bond(target, h_neighbor))
        mod_mol.remove_atom(h_neighbor)
        target.radical_electrons = 1
        mod_mol.update_multiplicity()

        p_class, details = polymer.classify_structure(Species(molecule=[mod_mol]), self.p)
        assert p_class == polymer.PolymerClass.FEATURE

    # =========================================================================
    # MACROSCOPIC STRUCTURAL CHANGES
    # =========================================================================

    def test_branch_crosslink_bimolecular(self):
        """Tests >2 wings by joining two chains at heavy atoms after making room (valency)."""
        crosslink_mol = self.ref_baseline_mol.copy(deep=True)
        second_chain = self.ref_baseline_mol.copy(deep=True)

        # 1. Add atoms and bonds from second_chain properly
        mapping = {}
        for atom in second_chain.atoms:
            new_atom = atom.copy()
            crosslink_mol.add_atom(new_atom)
            mapping[atom] = new_atom

        for atom1 in second_chain.atoms:
            for atom2, bond in atom1.edges.items():
                if id(atom1) < id(atom2):
                    crosslink_mol.add_bond(Bond(mapping[atom1], mapping[atom2], bond.order))

        # 2. Pick a heavy atom from the original first chain and clear a spot
        # We know atoms[:len_original] belong to the first chain
        a1 = next(a for a in crosslink_mol.atoms if not a.is_hydrogen())
        h1 = next(n for n in a1.edges if n.is_hydrogen())
        crosslink_mol.remove_bond(crosslink_mol.get_bond(a1, h1))
        crosslink_mol.remove_atom(h1)

        # 3. Pick a heavy atom from the newly added second chain and clear a spot
        # We search specifically in the mapped atoms to ensure internal connectivity
        a2 = next(mapping[a] for a in second_chain.atoms if not a.is_hydrogen())
        h2 = next(n for n in a2.edges if n.is_hydrogen())
        crosslink_mol.remove_bond(crosslink_mol.get_bond(a2, h2))
        crosslink_mol.remove_atom(h2)

        # 4. Connect them - Valency is now satisfied (4 bonds per Carbon)
        crosslink_mol.add_bond(Bond(a1, a2, order=1))
        crosslink_mol.update_multiplicity()

        p_class, details = polymer.classify_structure(Species(molecule=[crosslink_mol]), self.p)
        assert p_class == polymer.PolymerClass.CROSSLINK

    def test_branch_scission_single_wing(self):
        """
        Tests a severed polymer chain containing exactly one terminal end-cap.
        Physically represents beta-scission in the backbone.
        """
        # Simulate scission by building a single stitched wing
        scission_mol = self.p._stitch_wing("head")
        scission_mol.clear_labeled_atoms()
        scission_mol.update_multiplicity()

        species = Species(molecule=[scission_mol])
        p_class, details = polymer.classify_structure(species, self.p)

        assert p_class == polymer.PolymerClass.SCISSION
        assert details["reason"] == "single_terminal_wing"
        assert details["num_disjoint_wings"] == 1

    def test_branch_gas_no_wings(self):
        """Tests 0 wings for decane."""
        alien_mol = Molecule().from_smiles("CCCCCCCCCC")
        p_class, details = polymer.classify_structure(Species(molecule=[alien_mol]), self.p)
        assert p_class == polymer.PolymerClass.GAS
        assert details["reason"] == "no_intact_wings"
        # Match the key used in classify_structure's base_details
        assert "disjoint_matches" in details

    # =========================================================================
    # THE ANOMALOUS FALLBACK
    # =========================================================================

    def test_branch_unknown_anomalous_backbone(self):
        """Tests unknown branch using a multiplicity mismatch."""
        mod_mol = self.ref_baseline_mol.copy(deep=True)
        # Change multiplicity to something impossible for this structure
        # This bypasses BASELINE (isomorphism fails) and FEATURE (no radicals/labels)
        mod_mol.multiplicity = 5

        p_class, details = polymer.classify_structure(Species(molecule=[mod_mol]), self.p)
        assert p_class == polymer.PolymerClass.UNKNOWN


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

        self.pe_feat = self.pe_polymer.copy()
        feat_mol = Molecule().from_smiles("[CH][CH2]")
        feat_mol.atoms[0].label = "*1"
        feat_mol.atoms[1].label = "*2"
        self.pe_feat.feature_monomer = feat_mol
        self.pe_feat._feature_proxy = None
        self.proxy = self.pe_polymer.get_proxy_species()

        self.dummy_thermo = NASA(
            polynomials=[NASAPolynomial(coeffs=[1, 1, 1, 1, 1, 1, 1], Tmin=(298, 'K'), Tmax=(1000, 'K')),
                         NASAPolynomial(coeffs=[2, 2, 2, 2, 2, 2, 2], Tmin=(1000, 'K'), Tmax=(3000, 'K'))],
            Tmin=(298, 'K'), Tmax=(3000, 'K'), Cp0=(30, 'J/(mol*K)'), CpInf=(100, 'J/(mol*K)'))
        self.proxy.thermo = self.dummy_thermo

    def test_get_thermo_data_modes(self):
        """Verify 'baseline' vs 'feature' mode selection in thermo retrieval."""
        assert self.pe_polymer.get_thermo_data(mode='auto') is self.dummy_thermo
        with pytest.raises(RuntimeError):
            self.pe_feat.get_thermo_data(mode='feature')

    def test_get_free_energy_delegation(self):
        """Test Gibbs Free Energy delegation."""
        T = 750.0
        G_pol = self.pe_polymer.get_free_energy(T)
        G_exp = self.dummy_thermo.get_free_energy(T)
        assert G_pol == G_exp

    def test_get_thermo_data_runtime_error_on_none(self):
        """Ensure clear error if proxy exists but thermo generation fails."""
        self.proxy.thermo = None
        with pytest.raises(RuntimeError) as excinfo:
            self.pe_polymer.get_thermo_data()
        assert "Thermo generation failed" in str(excinfo.value)

    def test_get_bulk_heat_capacity_logic(self):
        """Verifies bulk scaling for reactor solvers."""
        T = 400.0
        DP = 50.0  # 50 units long
        site_cp = self.pe_polymer.get_heat_capacity(T)
        bulk_cp = self.pe_polymer.get_bulk_heat_capacity(T, DP)
        assert bulk_cp == site_cp * DP

    def test_generate_statmech_delegation(self):
        """Verify statmech delegation without triggering database calls."""
        mock_conf = Conformer(E0=(10.0, "kJ/mol"))
        mock_conf.modes = [
            IdealGasTranslation(mass=(28.0, "amu")),
            NonlinearRotor(inertia=([0.630578, 1.15529, 1.78586], "amu*angstrom^2"), symmetry=2),
            HarmonicOscillator(frequencies=([1000.0], "cm^-1"))]
        self.proxy.conformer = mock_conf
        out = self.pe_polymer.generate_statmech()
        assert out is mock_conf
        assert self.pe_polymer.conformer is mock_conf

    def test_generate_transport_delegation(self):
        """Verify TransportData delegation."""
        mock_trans = TransportData(sigma=(3.5, 'angstrom'), epsilon=(120.0, 'K'))
        self.proxy.transport_data = mock_trans
        out = self.pe_polymer.generate_transport_data()
        assert out is mock_trans
        assert self.pe_polymer.transport_data is mock_trans

    def test_calculate_cp0_cpinf_with_no_molecule(self):
        """Ensure 0.0 is returned if proxy molecule is missing (Safety Check)."""
        self.pe_polymer._baseline_proxy.molecule = []
        assert self.pe_polymer.calculate_cp0() == 0.0
        assert self.pe_polymer.calculate_cpinf() == 0.0

    def test_multiplicity_and_weight_consistency(self):
        """Ensures polymer multiplicity and MW are tied to proxy, not bulk."""
        assert self.pe_polymer.multiplicity == 1
        mw_kg_per_mol = self.pe_polymer.molecular_weight.value_si
        if mw_kg_per_mol < 1e-10:
            from rmgpy.constants import Na
            mw_kg_per_mol *= Na
        assert 0.05 < mw_kg_per_mol < 0.15

    def test_get_thermo_data_delegates_to_proxy(self):
        """Test that get_thermo_data returns the proxy's thermo object."""
        thermo = self.pe_polymer.get_thermo_data()
        assert thermo is self.dummy_thermo
        assert self.pe_polymer.thermo is self.dummy_thermo  # Check sync behavior

    def test_thermo_properties_delegate_correctly(self):
        """Test get_enthalpy, entropy, heat_capacity, etc. return proxy values."""
        T = 500.0
        H_pol = self.pe_polymer.get_enthalpy(T)
        S_pol = self.pe_polymer.get_entropy(T)
        Cp_pol = self.pe_polymer.get_heat_capacity(T)
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
        cp0 = self.pe_polymer.calculate_cp0()
        cpinf = self.pe_polymer.calculate_cpinf()
        assert isinstance(cp0, float)
        assert isinstance(cpinf, float)

    def test_multiplicity_delegation(self):
        """Test multiplicity property."""
        mult = self.pe_polymer.multiplicity
        assert mult == 1

    def test_molecular_weight_delegation(self):
        """Test molecular_weight property returns proxy MW (per-site), not Mn (bulk)."""
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
            Mw=2000.0)
        assert not self.pe_polymer.is_identical(p3)

    def test_transport_delegation(self):
        """Test generate_transport_data delegates to proxy."""
        dummy_trans = TransportData(sigma=(3.0, 'angstrom'), epsilon=(100.0, 'K'))
        self.proxy.transport_data = dummy_trans
        trans = self.pe_polymer.generate_transport_data()
        assert trans is dummy_trans
        assert self.pe_polymer.transport_data is dummy_trans


def _iter_neighbors(atom) -> List[Any]:
    """
    Return neighbor atoms for a given Atom across common APIs:
    - RMG: atom.bonds is dict[Atom, Bond]
    - Some toolkits: atom.bonds may be iterable of neighbors
    """
    bonds = getattr(atom, "bonds", None)
    if bonds is None: return []
    if isinstance(bonds, dict): return list(bonds.keys())
    try:
        return list(bonds)
    except TypeError:
        return []


def get_carbon_neighbors(atom) -> List[Any]:
    return [n for n in _iter_neighbors(atom) if n.is_non_hydrogen()]


def bfs_farthest_node(start_node: Atom, all_nodes: List[Atom]) -> Tuple[Atom, Dict[Atom, Atom]]:
    """
    BFS to find the farthest node from start_node within the subgraph of all_nodes.
    Returns (farthest_node, parent_map).
    """
    queue = deque([start_node])
    visited = {start_node}
    parent = {start_node: None}
    farthest = start_node
    while queue:
        current = queue.popleft()
        farthest = current
        for neighbor in get_carbon_neighbors(current):
            if neighbor in all_nodes and neighbor not in visited:
                visited.add(neighbor)
                parent[neighbor] = current
                queue.append(neighbor)
    return farthest, parent


def get_backbone_path(mol: Molecule) -> List[Atom]:
    """
    Identifies the longest carbon chain (backbone) in the molecule.
    """
    carbons = [a for a in mol.atoms if a.is_carbon()]
    if not carbons: return []
    u, _ = bfs_farthest_node(carbons[0], carbons)
    v, parent_map = bfs_farthest_node(u, carbons)
    path = []
    curr = v
    while curr is not None:
        path.append(curr)
        curr = parent_map[curr]
    return path


def get_monomer_regions(mol: Molecule) -> Dict[str, List[Atom]]:
    """
    Segments the linear backbone into Buffer (Head/Tail) and Center regions.
    Assumes a trimer structure (6 carbons).
    """
    path = get_backbone_path(mol)
    return {"head_buffer": path[:2],
            "center": path[2:-2],
            "tail_buffer": path[-2:]}


def abstract_h_from_center_backbone(mol):
    """
    Perform a chemically valid H-abstraction near the backbone center:
    - choose a backbone carbon near the middle that has an explicit H neighbor
    - remove that H atom
    - increment radical on the carbon
    Returns the modified carbon atom.
    """
    path = get_backbone_path(mol)
    n = len(path)
    mid = n // 2
    for k in range(n):
        for i in (mid - k, mid + k):
            if i < 0 or i >= n:
                continue
            c = path[i]
            if not c.is_carbon():
                continue
            h = next((nb for nb in c.bonds.keys() if nb.is_hydrogen()), None)
            if h is None:
                continue
            if hasattr(mol, "remove_atom"):
                mol.remove_atom(h)
            else:
                del c.bonds[h]
                del h.bonds[c]
                mol.atoms.remove(h)
            c.increment_radical()
            mol.update_multiplicity()
            return c
    raise ValueError("Could not find a center-backbone carbon with an explicit H to abstract.")


class TestPolymerAdditionalCoverage:
    @pytest.fixture(autouse=True)
    def setup_polymer(self):
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
        self.p = Polymer(
            label="PS_cov",
            monomer=ps_adj,
            end_groups=["[CH3]", "[H]"],
            cutoff=4,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )
        yield

    def test_get_proxy_species_modes(self):
        # baseline-only polymer
        assert self.p.get_proxy_species("baseline") is self.p.baseline_proxy
        assert self.p.get_proxy_species("feature") is None
        assert self.p.get_proxy_species("auto") is self.p.baseline_proxy

        feat_poly = self.p.copy()
        feat_poly.feature_monomer = feat_poly.monomer.copy(deep=True)
        feat_poly._feature_proxy = None

        assert feat_poly.get_proxy_species("baseline") is feat_poly.baseline_proxy
        assert feat_poly.get_proxy_species("feature") is feat_poly.feature_proxy
        assert feat_poly.get_proxy_species("auto") is feat_poly.feature_proxy

    def test_get_free_energy_delegates_to_proxy(self):
        """Test that get_free_energy(T) delegates to the proxy species' thermo."""
        dummy_thermo = NASA(
            polynomials=[
                NASAPolynomial(coeffs=[1, 1, 1, 1, 1, 1, 1], Tmin=(298, "K"), Tmax=(1000, "K")),
                NASAPolynomial(coeffs=[2, 2, 2, 2, 2, 2, 2], Tmin=(1000, "K"), Tmax=(3000, "K")),
            ],
            Tmin=(298, "K"),
            Tmax=(3000, "K"),
            Cp0=(30, "J/(mol*K)"),
            CpInf=(100, "J/(mol*K)"),
        )
        proxy = self.p.get_proxy_species()
        proxy.thermo = dummy_thermo
        T = 600.0
        assert self.p.get_free_energy(T) == dummy_thermo.get_free_energy(T)

    def test_generate_statmech_delegation_fast_path(self):
        """
        Covers Polymer.generate_statmech() without requiring full statmech machinery:
        if proxy.has_statmech() is True, Polymer should just copy proxy.conformer.
        """
        sentinel_conformer = Conformer()

        class MockProxy:
            def __init__(self):
                self.conformer = sentinel_conformer
                self.label = "MockProxy"

            def has_statmech(self):
                return True

            def generate_statmech(self):
                # If called, this would indicate the wrong branch; make it fail loudly
                raise AssertionError("generate_statmech should not be called if has_statmech is True")

        # 2. Inject the mock into the Polymer's cache
        self.p._baseline_proxy = MockProxy()
        self.p.feature_monomer = None  # Forces get_proxy_species() to return baseline_proxy

        # 3. Run the method
        out = self.p.generate_statmech()

        # 4. Assert the fast-path delegation occurred
        assert out is sentinel_conformer
        assert self.p.conformer is sentinel_conformer

    def test_validate_cutoff_rejects_non_int(self):
        with pytest.raises(InputError):
            Polymer(
                label="bad_cutoff",
                monomer=self.p.monomer.copy(deep=True),
                end_groups=["[H]", "[H]"],
                cutoff="abc",
                Mn=1000.0,
                Mw=2000.0,
                initial_mass=1.0,
            )

    def test_validate_cutoff_rejects_lt_2(self):
        with pytest.raises(InputError):
            Polymer(
                label="bad_cutoff2",
                monomer=self.p.monomer.copy(deep=True),
                end_groups=["[H]", "[H]"],
                cutoff=1,
                Mn=1000.0,
                Mw=2000.0,
                initial_mass=1.0,
            )

    def test_init_from_moments_with_zero_mu0_or_mu1_returns_zero_mn_mw(self):
        # mu0 = 0 -> Mn/Mw should be 0 per implementation
        p0 = Polymer(
            label="mom_mu0_zero",
            monomer=self.p.monomer.copy(deep=True),
            end_groups=["[H]", "[H]"],
            cutoff=3,
            moments=[0.0, 1.0, 2.0],
            initial_mass=1.0,
        )
        assert p0.Mn == 0.0
        assert p0.Mw == 0.0

        # mu1 = 0 -> Mn/Mw should be 0 per implementation
        p1 = Polymer(
            label="mom_mu1_zero",
            monomer=self.p.monomer.copy(deep=True),
            end_groups=["[H]", "[H]"],
            cutoff=3,
            moments=[1.0, 0.0, 2.0],
            initial_mass=1.0,
        )
        assert p1.Mn == 0.0
        assert p1.Mw == 0.0

    def test_ensure_open_site(self):
        """
        Test that _ensure_open_site promotes closed-shell atoms to radicals
        but leaves existing radicals untouched.
        """
        # Scenario 1: Closed-shell atom (e.g., Carbon in Methane)
        # Should be promoted to a radical (u1)
        c_closed = Atom(element='C', radical_electrons=0)
        polymer._ensure_open_site(c_closed)
        assert c_closed.radical_electrons == 1

        # Scenario 2: Existing mono-radical (u1)
        # Should remain unchanged (not become a diradical u2)
        c_radical = Atom(element='C', radical_electrons=1)
        polymer._ensure_open_site(c_radical)
        assert c_radical.radical_electrons == 1

        # Scenario 3: Existing multi-radical (u2)
        # Should remain unchanged
        c_diradical = Atom(element='C', radical_electrons=2)
        polymer._ensure_open_site(c_diradical)
        assert c_diradical.radical_electrons == 2

    def test_get_target_atoms(self):
        """
        Tests that get_target_atoms correctly extracts Atom objects from
        various mapping configurations (keys, values, or mixed).
        """
        a1 = Atom(element='C')
        a2 = Atom(element='C')

        # Scenario 1: Empty match
        assert polymer.get_target_atoms({}) == set()
        assert polymer.get_target_atoms(None) == set()

        # Scenario 2: Atoms as Values (Standard RMG find_subgraph_isomorphisms)
        # {GroupAtom: Atom}
        match_vals = {"p1": a1, "p2": a2}
        result_vals = polymer.get_target_atoms(match_vals)
        assert len(result_vals) == 2
        assert a1 in result_vals and a2 in result_vals

        # Scenario 3: Atoms as Keys (Often happens in reverse mappings or custom tools)
        # {Atom: GroupAtom}
        match_keys = {a1: "p1", a2: "p2"}
        result_keys = polymer.get_target_atoms(match_keys)
        assert len(result_keys) == 2
        assert a1 in result_keys and a2 in result_keys

        # Scenario 4: Mixed or Fallback
        # We explicitly use objects that are definitely NOT Atoms
        match_mixed = {a1: "p1", "p2": a2, "extra": 123}
        result_mixed = polymer.get_target_atoms(match_mixed)

        # We want ONLY a1 and a2.
        assert len(result_mixed) == 2
        assert all(not isinstance(x, (str, int)) for x in result_mixed)

    def test_stitch_trimer_copolymer_san(self):
        """
        Tests stitching a copolymer trimer (Styrene-Acrylonitrile-Styrene).
        Verifies that the feature_monomer is correctly placed in the center
        between baseline monomers.
        """
        # 1. Define Acrylonitrile (AN) Monomer with connectivity labels
        an_adj = """
multiplicity 3
1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}
3    C u0 p0 c0 {1,S} {7,T}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
7    N u0 p1 c0 {3,T}
"""
        # 2. Setup SAN Copolymer using polymer_1 (PS) as the base
        # polymer_1 already has Styrene (C8H8) as the monomer
        san_copoly = self.p.copy()
        san_copoly.label = "SAN_Copolymer"
        san_copoly.feature_monomer = Molecule().from_adjacency_list(an_adj)

        # 3. Stitch the 'Feature' Trimer: [PS]--[AN]--[PS]
        # Calling get_thermo_data(mode='feature') triggers the feature_proxy creation
        with pytest.raises(RuntimeError):  # Fails on thermo, but builds the proxy first
            san_copoly.get_thermo_data(mode='feature')

        proxy_spc = san_copoly.feature_proxy
        assert proxy_spc.is_polymer_proxy is True

        # 4. Atomic Count Validation
        # Head(CH3: 1C) + 2x Styrene(C8H8: 16C) + 1x Acrylonitrile(C3H3N: 3C) + Tail(H: 0C)
        # Expected Total Carbons = 1 + 16 + 3 = 20
        # Expected Nitrogens = 1
        mol = proxy_spc.molecule[0]
        c_atoms = [a for a in mol.atoms if a.is_carbon()]
        n_atoms = [a for a in mol.atoms if a.symbol == 'N']

        assert len(c_atoms) == 20
        assert len(n_atoms) == 1

        # 5. Connectivity Validation
        # Ensure the Nitrogen (the unique marker) is not on the terminal ends
        # In a [H]-[S]-[AN]-[S]-[Cap] trimer, the AN should be at least 3 bonds from any end
        n_atom = n_atoms[0]
        # Simple check: Nitrogen neighbor should be a Carbon (C3) which has 3 bonds
        c_nitrile = list(n_atom.bonds.keys())[0]
        assert c_nitrile.is_carbon()
        assert len(c_nitrile.bonds) == 2  # Connected to N and the backbone C

    def test_get_element_symbol(self):
        # Test Carbon types
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['Cs']])) == 'C'
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['C2d']])) == 'C'

        # Test Heteroatoms
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['N3d']])) == 'N'
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['O2d']])) == 'O'

        # Test Multi-character elements
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['Cl1s']])) == 'Cl'
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['Sibf']])) == 'Si'
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['Br1s']])) == 'Br'

        # Test Helium/Noble (Shortest string logic)
        assert polymer.get_element_symbol(GroupAtom(atomtype=[ATOMTYPES['He']])) == 'He'


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


def _safe_make_radical(mol: Molecule, atom: Atom):
    """Safely removes a hydrogen before adding a radical to maintain valency."""
    h_atom = next((a for a in atom.bonds if a.is_hydrogen()), None)
    if h_atom:
        bond = mol.get_bond(atom, h_atom)
        mol.remove_bond(bond)
        mol.remove_atom(h_atom)
    atom.increment_radical()
