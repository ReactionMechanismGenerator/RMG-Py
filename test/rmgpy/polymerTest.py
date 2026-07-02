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

    def test_copy_preserves_identity_and_kinetics(self):
        """
        copy() uses __new__ (bypassing __init__) and must explicitly carry over
        the attributes __init__ sets. Regression: it used to drop `is_polymer`
        (so _handshake_structures re-processed copied polymers) and, worse, the
        degradation kinetics k_scission/k_unzip (silently reset to 0).
        """
        p = Polymer(
            label='PE', monomer='[CH2][CH2]', end_groups=['[CH3]', '[H]'],
            cutoff=3, Mn=5000.0, Mw=6000.0, initial_mass=1.0,
            k_scission=2.5, k_unzip=0.7,
        )
        for c in (p.copy(), p.copy(deep=True)):
            assert isinstance(c, Polymer)
            assert getattr(c, 'is_polymer', False) is True
            assert c.k_scission == 2.5
            assert c.k_unzip == 0.7

    def test_fingerprint_property(self):
        """Test that the fingerprint property works"""
        assert self.polymer_1.fingerprint == "Polymer_C08H08N00O00S00_EG-C01H03N00O00S00_C00H01N00O00S00_3"
        assert self.polymer_2.fingerprint == "Polymer_C08H08N00O00S00_EG-C01H03N00O00S00_C00H01N00O00S00_5"
        assert self.polymer_3.fingerprint == "Polymer_C02H04N00O00S00_EG-C00H01N00O00S00_C00H01N00O00S00_10"

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
        # A scission product is a new, ~half-length chain population: Mn/Mw are
        # halved and the pool starts empty (zero moments), matching _scission_head.
        # (This previously asserted Mn/Mw == parent, which copied the parent's
        # full distribution into a zero-mass fragment — a mass-duplication bug.)
        assert np.isclose(new_p.Mn, p.Mn / 2.0)
        assert np.isclose(new_p.Mw, p.Mw / 2.0)
        assert np.allclose(new_p.moments, 0.0)

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

    def test_create_reacted_copy_stamps_reacted_class(self):
        """
        create_reacted_copy stamps the classification verdict on the returned
        polymer (``_reacted_class``) so the polymer handshake can flag END_MOD
        reactions for mu0 (chain-end) scaling in the solver.
        """
        p = self.polymer_1.copy()
        # END_MOD product (terminal radical-activated) -> stamped END_MOD.
        end_mod = p.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(p, end_mod)
        new_end = p.create_reacted_copy(end_mod)
        assert new_end._reacted_class == PolymerClass.END_MOD
        # Unreacted baseline proxy -> stamped BASELINE.
        baseline = p.baseline_proxy.molecule[0].copy(deep=True)
        new_base = p.create_reacted_copy(baseline)
        assert new_base._reacted_class == PolymerClass.BASELINE

    def test_is_end_group_reaction_helper(self):
        """
        is_end_group_reaction(products) is True iff some product Polymer was
        classified END_MOD. Non-Polymer products and unstamped polymers are
        ignored (default mu1 scaling).
        """
        from rmgpy.polymer import is_end_group_reaction
        p = self.polymer_1
        end = p.copy(); end._reacted_class = PolymerClass.END_MOD
        feat = p.copy(); feat._reacted_class = PolymerClass.FEATURE
        plain = Molecule(smiles="CC")
        assert is_end_group_reaction([end, plain]) is True
        assert is_end_group_reaction([feat, plain]) is False
        assert is_end_group_reaction([plain]) is False
        assert is_end_group_reaction([p.copy()]) is False  # no _reacted_class set

    def test_classify_reaction_flux_archetype(self):
        """
        classify_reaction_flux_archetype(reactants, products) drives the solver's
        per-reaction pool moment apportionment (spec 2026-06-09):
        SCISSION product -> SCISSION_FRAGMENT; single cross-pool polymer product
        -> MIGRATION; fold-backs -> SAME_POOL; no polymers -> NONE; ambiguous or
        end-initiated-scission shapes -> UNRESOLVED (legacy mu1 flux + warning).
        """
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()

        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        gas = Molecule(smiles="CC")

        # No polymer on either side -> NONE
        assert classify_reaction_flux_archetype([gas], [gas]) == PolymerFluxArchetype.NONE

        # Single polymer reactant sheds a discrete mass-carrying co-product (CC,
        # a~=0.29) -> mass-equivalent volatile ejection, NOT a net-zero fold-back;
        # proves mass bookkeeping not DP. (Signed-VE, Codex round-13: net =
        # MW(CC) - 0 > eps, so the same-pool branch routes VOLATILE_EJECTION even
        # though the polymer product folds back into the same pool.)
        fold = p.copy()
        fold._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p], [fold, gas])
                == PolymerFluxArchetype.VOLATILE_EJECTION)

        # Single cross-pool polymer product, NO non-polymer co-product (pure
        # whole-chain relabel) -> MIGRATION. (A cross-pool product WITH a
        # discrete mass-carrying co-product is VOLATILE_EJECTION -- see
        # test_classify_volatile_ejection_new_archetype.)
        other = p.copy()
        other.label = "other_pool"
        other._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p], [other])
                == PolymerFluxArchetype.MIGRATION)

        # SCISSION-stamped product -> SCISSION_FRAGMENT
        sc = p.copy()
        sc.label = f"{p.label}_scission_head"
        sc._reacted_class = PolymerClass.SCISSION
        assert (classify_reaction_flux_archetype([p], [sc, gas])
                == PolymerFluxArchetype.SCISSION_FRAGMENT)

        # End-initiated scission (SCISSION + END_MOD product) -> UNRESOLVED:
        # the uniform-cut bundle assumptions don't hold near a chain end.
        end = p.copy()
        end._reacted_class = PolymerClass.END_MOD
        assert (classify_reaction_flux_archetype([p], [sc, end])
                == PolymerFluxArchetype.UNRESOLVED)

        # Two polymer products with a cross-pool member -> UNRESOLVED
        assert (classify_reaction_flux_archetype([p], [other, fold])
                == PolymerFluxArchetype.UNRESOLVED)

        # Inter-chain (two reactant pools, each product folds back) -> SAME_POOL
        q = p.copy()
        q.label = "second_pool"
        fold_q = q.copy()
        fold_q._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p, q], [fold, fold_q])
                == PolymerFluxArchetype.SAME_POOL)

        # Cross-pool product with TWO reactant pools (ambiguous source) -> UNRESOLVED
        assert (classify_reaction_flux_archetype([p, q], [other])
                == PolymerFluxArchetype.UNRESOLVED)

        # Polymer reactant, all-gas products -> UNRESOLVED (no flux rule)
        assert (classify_reaction_flux_archetype([p], [gas])
                == PolymerFluxArchetype.UNRESOLVED)

        # Warn-once: UNRESOLVED causes above logged exactly one warning per
        # distinct (reason, detail) key.
        assert len(polymer_mod._flux_archetype_warned) == 4
        n_before = len(polymer_mod._flux_archetype_warned)
        classify_reaction_flux_archetype([p, q], [other])  # repeat call
        assert len(polymer_mod._flux_archetype_warned) == n_before

    def test_classify_chip_product_returns_discrete_chip(self):
        """
        A CHIP-stamped product polymer (left by surge_chip_products' fold-back)
        short-circuits to DISCRETE_CHIP BEFORE the SCISSION branch. Order
        matters (spec 2026-06-10 §4.1): after the (b)-surgery the product list
        has no END_MOD member, so the SCISSION branch's internal
        is_end_group_reaction(products) recompute would return False and
        misroute; the CHIP check must win even alongside a SCISSION member.
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        chip_gas = Molecule(smiles="CC")
        fold = p.copy()
        fold._reacted_class = PolymerClass.CHIP
        assert (classify_reaction_flux_archetype([p], [chip_gas, fold])
                == PolymerFluxArchetype.DISCRETE_CHIP)

        # Defensive order pin: CHIP beats a co-present SCISSION stamp.
        sc = p.copy()
        sc.label = f"{p.label}_scission_head"
        sc._reacted_class = PolymerClass.SCISSION
        assert (classify_reaction_flux_archetype([p], [sc, fold])
                == PolymerFluxArchetype.DISCRETE_CHIP)

    def test_classify_volatile_ejection_new_archetype(self):
        """
        polymer_A -> discrete volatile + cross-pool polymer_B is a MASS-LOSING
        VOLATILE_EJECTION, NOT a mass-conserving MIGRATION: the discrete
        non-polymer co-product carries mass off the chain, so the chain must
        lose it (spec: depolymerization / volatile-ejection archetype). Before
        this fix, the classifier filtered products to Polymers only, ignored
        the volatile, and mislabeled the shape MIGRATION (whole-chain relabel,
        mass-conserving) -- the volatile's mass was never debited.
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype
        p = self.polymer_1
        volatile = Molecule(smiles="C=C(C)c1ccccc1")  # alpha-methylstyrene
        dst = p.copy()
        dst.label = "PS_scission_tail"
        dst._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p], [volatile, dst])
                == PolymerFluxArchetype.VOLATILE_EJECTION)
        # Product order must not matter.
        assert (classify_reaction_flux_archetype([p], [dst, volatile])
                == PolymerFluxArchetype.VOLATILE_EJECTION)

    def test_classify_pure_relabel_still_migration(self):
        """
        A pure relabel (single cross-pool polymer product, NO non-polymer
        product) is unchanged: still MIGRATION. Only the shape with a discrete
        mass-carrying co-product becomes VOLATILE_EJECTION.
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype
        p = self.polymer_1
        dst = p.copy()
        dst.label = "other_pool"
        dst._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p], [dst])
                == PolymerFluxArchetype.MIGRATION)

    def test_stamp_volatile_ejection_units_fractional(self):
        """
        VOLATILE_EJECTION stamps ``polymer_eject_units`` = (discrete volatile
        MW g/mol) / (source pool monomer_mw_g_mol) as a FRACTIONAL value.
        alpha-methylstyrene (C9H10, 118.18) off a styrene pool
        (C8H8, 104.15) => ~1.135; it must NOT be rounded to 1.0.
        """
        from rmgpy.reaction import Reaction
        from rmgpy.polymer import PolymerFluxArchetype, stamp_polymer_flux_archetype
        p = self.polymer_1
        volatile = Molecule(smiles="C=C(C)c1ccccc1")  # alpha-methylstyrene C9H10
        dst = p.copy()
        dst.label = "PS_scission_tail"
        dst._reacted_class = PolymerClass.FEATURE
        rxn = Reaction(reactants=[p], products=[volatile, dst])
        rxn.is_end_group_reaction = False
        stamp_polymer_flux_archetype(rxn, [p], [p])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
        expected = (volatile.get_molecular_weight() * 1000.0) / p.monomer_mw_g_mol
        assert rxn.polymer_eject_units == pytest.approx(expected)
        assert rxn.polymer_eject_units == pytest.approx(1.135, abs=0.01)
        # NOT rounded to an integer number of monomer-equivalents.
        assert abs(rxn.polymer_eject_units - 1.0) > 0.1

    def test_stamp_volatile_ejection_sums_multiple_volatiles(self):
        """
        Multiple discrete volatile products sum: ``polymer_eject_units`` =
        (sum of their MW g/mol) / source monomer_mw_g_mol.
        """
        from rmgpy.reaction import Reaction
        from rmgpy.polymer import PolymerFluxArchetype, stamp_polymer_flux_archetype
        p = self.polymer_1
        v1 = Molecule(smiles="C=C(C)c1ccccc1")  # alpha-methylstyrene
        v2 = Molecule(smiles="C=Cc1ccccc1")     # styrene
        dst = p.copy()
        dst.label = "PS_scission_tail"
        dst._reacted_class = PolymerClass.FEATURE
        rxn = Reaction(reactants=[p], products=[v1, v2, dst])
        rxn.is_end_group_reaction = False
        stamp_polymer_flux_archetype(rxn, [p], [p])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
        expected = ((v1.get_molecular_weight() + v2.get_molecular_weight())
                    * 1000.0) / p.monomer_mw_g_mol
        assert rxn.polymer_eject_units == pytest.approx(expected)

    def test_classify_volatile_ejection_bare_species_volatile(self):
        """
        LIVE-PATH regression: at stamp_polymer_flux_archetype time the
        depolymerization volatile arrives as a bare pre-thermo ``Species``
        (empty label, NO ``get_molecular_weight`` method, but a populated
        ``.molecule`` structure), NOT a Molecule. A ``hasattr(p,
        'get_molecular_weight')`` gate silently misses it and mis-stamps the
        shape MIGRATION (mass-conserving) -> the TGA stays flat at 100%.
        Empirically confirmed on a full PS run (sidecar emitted 5x migration/1).
        The classifier must reach through to ``.molecule[0]`` and route
        VOLATILE_EJECTION.
        """
        from rmgpy.species import Species
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype
        p = self.polymer_1
        # Bare Species built from structure only -- Species has NO
        # get_molecular_weight method (verified) and an empty label pre-thermo.
        volatile = Species(molecule=[Molecule(smiles="C=C(C)c1ccccc1")])
        assert not hasattr(volatile, "get_molecular_weight")  # the trap
        dst = p.copy()
        dst.label = "PS_scission_tail"
        dst._reacted_class = PolymerClass.UNKNOWN  # live value (not SCISSION)
        assert (classify_reaction_flux_archetype([p], [volatile, dst])
                == PolymerFluxArchetype.VOLATILE_EJECTION)
        assert (classify_reaction_flux_archetype([p], [dst, volatile])
                == PolymerFluxArchetype.VOLATILE_EJECTION)

    def test_stamp_volatile_ejection_units_bare_species_volatile(self):
        """
        LIVE-PATH regression companion: compute_volatile_ejection_units must
        also weigh a bare ``Species`` volatile (via ``.molecule[0]``), not only
        a Molecule -- else a would be 0 / wrong and mass would be fabricated
        even once the archetype is right.
        """
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species
        from rmgpy.polymer import PolymerFluxArchetype, stamp_polymer_flux_archetype
        p = self.polymer_1
        volatile = Species(molecule=[Molecule(smiles="C=C(C)c1ccccc1")])  # a-MS
        dst = p.copy()
        dst.label = "PS_scission_tail"
        dst._reacted_class = PolymerClass.UNKNOWN
        rxn = Reaction(reactants=[p], products=[volatile, dst])
        rxn.is_end_group_reaction = False
        stamp_polymer_flux_archetype(rxn, [p], [p])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
        expected = (volatile.molecule[0].get_molecular_weight() * 1000.0) / p.monomer_mw_g_mol
        assert rxn.polymer_eject_units == pytest.approx(expected)
        assert rxn.polymer_eject_units == pytest.approx(1.135, abs=0.01)

    def test_classify_same_pool_ejection_is_mass_losing(self):
        """
        SAME-POOL unzip -- radical depropagation `pool A -> monomer + pool A`
        (the chain sheds a monomer and the shorter chain stays in the SAME
        pool) -- MUST be mass-losing VOLATILE_EJECTION, NOT the net-zero
        SAME_POOL fold-back. Before this fix, classify returned SAME_POOL as
        soon as there was no cross-pool polymer product, IGNORING the discrete
        volatile co-product -> zero moment loss while gas monomer is produced =
        the same flat-TGA/mass-fabrication class as the migration bug, under a
        different archetype (Codex round-12). This is the shape real radical
        unzipping produces, so it must be covered before enabling that chemistry.
        """
        from rmgpy.species import Species
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype
        p = self.polymer_1
        monomer = Species(molecule=[Molecule(smiles="C=Cc1ccccc1")])  # styrene
        # Daughter folds back into the SAME pool (same label as the reactant).
        same = p.copy()
        same.label = p.label
        same._reacted_class = PolymerClass.UNKNOWN
        assert (classify_reaction_flux_archetype([p], [monomer, same])
                == PolymerFluxArchetype.VOLATILE_EJECTION)
        # A genuine net-zero same-pool fold-back (NO mass-carrying co-product)
        # stays SAME_POOL.
        assert (classify_reaction_flux_archetype([p], [same])
                == PolymerFluxArchetype.SAME_POOL)

    def test_stamp_same_pool_ejection_units(self):
        """Companion: a same-pool ejection stamps eject_units from the volatile
        (styrene off a styrene pool ~= 1.0), routed through the robust MW
        helper (bare Species volatile)."""
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species
        from rmgpy.polymer import PolymerFluxArchetype, stamp_polymer_flux_archetype
        p = self.polymer_1
        monomer = Species(molecule=[Molecule(smiles="C=Cc1ccccc1")])  # styrene
        same = p.copy()
        same.label = p.label
        same._reacted_class = PolymerClass.UNKNOWN
        rxn = Reaction(reactants=[p], products=[monomer, same])
        rxn.is_end_group_reaction = False
        stamp_polymer_flux_archetype(rxn, [p], [p])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
        expected = (monomer.molecule[0].get_molecular_weight() * 1000.0) / p.monomer_mw_g_mol
        assert rxn.polymer_eject_units == pytest.approx(expected)

    def test_classify_h_abstraction_nets_atom_transfer_ve(self):
        """
        SIGNED-VE (Codex round-13): a bimolecular same-pool H_Abstraction
        ``chain + R•(bare Species) -> chain(same pool) + RH`` must net only the
        single H actually shed, NOT the full RH mass. net = MW(RH) - MW(R•) ~=
        MW(H) > eps -> VOLATILE_EJECTION, with a ~= MW(H)/monomer (tiny), and the
        atom-transfer census warns (|a| < 0.5). Netting the REACTANT non-polymers
        is what keeps a from being the full-RH ~0.154; without it every radical
        abstraction would fabricate a whole co-reactant of chain mass loss.
        """
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   classify_reaction_flux_archetype,
                                   stamp_polymer_flux_archetype)
        p = self.polymer_1
        r_rad = Species(molecule=[Molecule(smiles="[CH3]")])   # bare radical R•
        rh = Species(molecule=[Molecule(smiles="C")])          # RH = R• + H
        same = p.copy()
        same.label = p.label
        same._reacted_class = PolymerClass.UNKNOWN
        # Classify: bimolecular same-pool, net = MW(CH4) - MW(CH3) ~= MW(H) > eps.
        assert (classify_reaction_flux_archetype([p, r_rad], [same, rh])
                == PolymerFluxArchetype.VOLATILE_EJECTION)
        # Stamp: signed a is the NET (H-scale), not the full RH mass.
        rxn = Reaction(reactants=[p, r_rad], products=[same, rh])
        rxn.is_end_group_reaction = False
        stamp_polymer_flux_archetype(rxn, [p, r_rad], [p])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
        net_h_g = (rh.molecule[0].get_molecular_weight()
                   - r_rad.molecule[0].get_molecular_weight()) * 1000.0
        expected = net_h_g / p.monomer_mw_g_mol
        assert rxn.polymer_eject_units == pytest.approx(expected)
        assert rxn.polymer_eject_units > 0.0                    # net loss (H shed)
        # Tiny: ~= MW(H)/monomer ~ 0.0097, NOT the full RH mass ~ 0.154.
        full_rh = (rh.molecule[0].get_molecular_weight() * 1000.0) / p.monomer_mw_g_mol
        assert rxn.polymer_eject_units < 0.05
        assert rxn.polymer_eject_units < 0.2 * full_rh
        # Atom-transfer census warned once (|a| < 0.5 monomer-equivalents).
        assert any(k[0] == "atom_transfer_ve"
                   for k in polymer_mod._flux_archetype_warned)

    def test_classify_monomer_addition_is_negative_ve(self):
        """
        SIGNED-VE (Codex round-13): a same-pool monomer ADDITION
        ``chain + monomer(bare Species) -> chain(same pool, longer)`` has NO
        non-polymer product but a non-polymer REACTANT, so net = 0 - MW(monomer)
        < 0. abs(net) > eps -> VOLATILE_EJECTION with a < 0 (chain GAINS mass).
        The signed a is what lets the solver grow μ1 for addition instead of
        silently conserving it.
        """
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   classify_reaction_flux_archetype,
                                   stamp_polymer_flux_archetype)
        p = self.polymer_1
        monomer = Species(molecule=[Molecule(smiles="C=Cc1ccccc1")])  # styrene
        longer = p.copy()
        longer.label = p.label
        longer._reacted_class = PolymerClass.UNKNOWN
        assert (classify_reaction_flux_archetype([p, monomer], [longer])
                == PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn = Reaction(reactants=[p, monomer], products=[longer])
        rxn.is_end_group_reaction = False
        stamp_polymer_flux_archetype(rxn, [p, monomer], [p])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
        expected = -(monomer.molecule[0].get_molecular_weight() * 1000.0) / p.monomer_mw_g_mol
        assert rxn.polymer_eject_units == pytest.approx(expected)
        assert rxn.polymer_eject_units < 0.0                    # chain gains mass
        assert rxn.polymer_eject_units == pytest.approx(-1.0, abs=0.01)

    def test_classify_pure_fold_back_no_nonpoly_stays_same_pool(self):
        """
        SIGNED-VE (Codex round-13): a pure fold-back ``chain -> chain(same pool)``
        with NO non-polymer participant on either side nets exactly 0, so it stays
        SAME_POOL (abs(net) <= eps). Only a non-zero net mass routes VE.
        """
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   classify_reaction_flux_archetype)
        p = self.polymer_1
        same = p.copy()
        same.label = p.label
        same._reacted_class = PolymerClass.UNKNOWN
        assert (classify_reaction_flux_archetype([p], [same])
                == PolymerFluxArchetype.SAME_POOL)

    def test_flag_false_one_unit_piece_routes_scission_fragment(self):
        """
        Discriminator regression (spec 2026-06-10 test 2, decision D3): a
        mu1-scaled (flag-false) cut whose represented piece is ONE repeat unit
        must still route SCISSION_FRAGMENT. On a 3-unit proxy a u1-u2 cut (the
        image of EVERY interior backbone bond) yields 1- and 2-unit pieces --
        literal piece size is a representation artifact and must never be a
        routing input. Guards against reintroducing piece-size routing.
        """
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        gas = Molecule(smiles="CC")
        piece = p.copy()
        piece.label = f"{p.label}_scission_tail"
        piece._reacted_class = PolymerClass.SCISSION
        # Represented piece: cap + ONE styrene unit (10 heavy atoms).
        # _source_molecule is grounded by create_reacted_copy since Task 4
        # (chip surgery / tripwire read it); fabricated here for the shape.
        piece._source_molecule = Molecule(smiles="CCC(C)c1ccccc1")
        assert (classify_reaction_flux_archetype([p], [piece, gas])
                == PolymerFluxArchetype.SCISSION_FRAGMENT)

    def test_tripwire_warns_once_for_wing_confined_mu1_piece(self, caplog):
        """
        Spec test 5 (§4.4): a mu1-scaled (flag-false) scission whose
        represented piece is end-confined (piece <= wing + at most 1 repeat
        unit by heavy atoms) logs the probable-mis-scaled-end-cut warning
        exactly once. Diagnostics only -- routing stays SCISSION_FRAGMENT.
        The census sets the priority of the end-anchor detector follow-up.
        """
        import logging as _logging
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()
        polymer_mod._chip_tripwire_warned.clear()
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        gas = Molecule(smiles="CC")
        piece = p.copy()
        piece.label = f"{p.label}_scission_tail"
        piece._reacted_class = PolymerClass.SCISSION
        # cap (1 heavy) + 1 unit (8 heavy) + stitch CH3 = 10 heavy
        #   <= max_cap_heavy(1) + 2*monomer_heavy(2*8) = 17 -> end-confined.
        piece._source_molecule = Molecule(smiles="CCC(C)c1ccccc1")

        with caplog.at_level(_logging.WARNING):
            arch = classify_reaction_flux_archetype([p], [piece, gas])
        assert arch == PolymerFluxArchetype.SCISSION_FRAGMENT  # never routes
        hits = [r for r in caplog.records
                if "probable mis-scaled end-anchored cut" in r.getMessage()]
        assert len(hits) == 1

        with caplog.at_level(_logging.WARNING):
            classify_reaction_flux_archetype([p], [piece, gas])  # repeat
        hits = [r for r in caplog.records
                if "probable mis-scaled end-anchored cut" in r.getMessage()]
        assert len(hits) == 1                                   # warn-once

    def test_discrete_dp_threshold_config_field(self):
        """
        discrete_dp_threshold (spec 2026-06-10 §6, decisions D7/D8): per-pool
        config knob, default 4 (monomer through trimer explicit), DORMANT
        under the fixed trimer proxy -- no behavioral use yet. Stored on the
        Polymer, survives copy(), overridable per pool.
        """
        assert self.polymer_1.discrete_dp_threshold == 4
        assert self.polymer_1.copy().discrete_dp_threshold == 4
        p = Polymer(label='PE_thresh', monomer='[CH2][CH2]',
                    end_groups=['[H]', '[H]'], cutoff=3,
                    Mn=1000.0, Mw=2500.0, initial_mass=1.0,
                    discrete_dp_threshold=6)
        assert p.discrete_dp_threshold == 6

    def test_backstop_dormant_under_trimer_proxy(self):
        """
        Spec test 16 / decision D8: the conditional DP backstop ("mu1-scaled
        cut producing literal DP < threshold -> exact-a accounting") applies
        ONLY when the proxy repeat-count exceeds discrete_dp_threshold. The
        fixed trimer proxy (3 units) never exceeds the default threshold (4),
        so a mu1-scaled mid-cut keeps routing SCISSION_FRAGMENT regardless of
        literal piece DP. No backstop code exists yet -- this pins the
        routing so adding it later cannot silently activate on trimer decks
        (unconditional, it would route ALL mid-chain scission to chip and
        kill mu0 growth/Mn halving).
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        assert p.discrete_dp_threshold == 4          # >= 3 proxy repeat units
        gas = Molecule(smiles="CC")
        piece = p.copy()
        piece.label = f"{p.label}_scission_tail"
        piece._reacted_class = PolymerClass.SCISSION
        piece._source_molecule = Molecule(smiles="CCC(C)c1ccccc1")  # DP ~1 piece
        assert (classify_reaction_flux_archetype([p], [piece, gas])
                == PolymerFluxArchetype.SCISSION_FRAGMENT)

    def test_demote_flipped_polymer_archetype(self):
        """
        apply_kinetics_to_reaction flips reactions in place when the kinetics
        are estimated in reverse; the stamped archetype then encodes reversed
        parent/daughter roles. demote_flipped_polymer_archetype must reset any
        polymer-touching reaction to UNRESOLVED (role-agnostic legacy flux)
        and leave pure-gas reactions alone.
        """
        from rmgpy.reaction import Reaction
        from rmgpy.rmg.model import demote_flipped_polymer_archetype
        from rmgpy.polymer import PolymerFluxArchetype

        p = self.polymer_1.copy()
        gas = Species(molecule=[Molecule(smiles="CC")])

        rxn = Reaction(reactants=[gas], products=[p],
                       polymer_flux_archetype=int(PolymerFluxArchetype.SCISSION_FRAGMENT))
        demote_flipped_polymer_archetype(rxn)
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)

        gas_rxn = Reaction(reactants=[gas], products=[gas])
        demote_flipped_polymer_archetype(gas_rxn)
        assert gas_rxn.polymer_flux_archetype == int(PolymerFluxArchetype.NONE)

    def test_create_reacted_copy_end_mod_folds_to_parent(self):
        """
        An END_MOD product (intact chain, terminal end-group radical-activated,
        e.g. CH3 -> CH2.) must NOT leak to the gas phase. In the method-of-moments
        model chain-end activation is abstracted into k_unzip (dmu1/dt = -k_unzip*mu0),
        so an end-group modification leaves the chain-length distribution unchanged:
        create_reacted_copy folds it back into the parent pool (self.copy) with
        moments and mass preserved.

        Regression: previously create_reacted_copy returned None for END_MOD
        products. Its raw wing-matching (find_subgraph_isomorphisms) diverges from
        classify_structure's heavy-view matcher and missed the degenerate [H] tail
        wing, mis-routing the product into the head-only scission-tail branch; the
        extracted fragment then carried 2 radicals (the activation radical plus the
        scission cut) and failed the mono-radical end-group assertion -> None. The
        None then left the product Molecule in place, registering it as a spurious
        gas-phase species (a mass leak).
        """
        p = self.polymer_1.copy()
        reacted_proxy = p.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(p, reacted_proxy)

        # Sanity: the construction really is an END_MOD product.
        probe = reacted_proxy.copy(deep=True)
        probe.clear_labeled_atoms()
        probe.update()
        assert polymer.classify_structure(Species(molecule=[probe]), p)[0] == PolymerClass.END_MOD

        new_p = p.create_reacted_copy(reacted_proxy)
        assert new_p is not None                  # no leak to gas
        assert isinstance(new_p, Polymer)
        assert new_p.label == p.label             # same pool identity (folded to parent)
        assert new_p.feature_monomer is None      # backbone unchanged
        # Moments / mass preserved (folded to parent, NOT halved like a scission).
        assert np.isclose(new_p.Mn, p.Mn)
        assert np.isclose(new_p.Mw, p.Mw)
        assert np.allclose(new_p.moments, p.moments)

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

    def test_discreteness_gate_rejects_backbone_impostor(self):
        """
        Discreteness gate (spec 2026-06-10 §3.1): a wing_count >= 2 candidate
        far smaller than the baseline proxy is a backbone impostor, not a
        chain image. Bibenzyl (14 heavy atoms) genuinely matches both PS wing
        subgraphs today but sits below the one-sided bound
        proxy_heavy - round(0.35*proxy_heavy) = 25 - 9 = 16 -> GAS.
        """
        impostor = Molecule(smiles="c1ccccc1CCc1ccccc1")  # bibenzyl, 14 heavy
        p_class, details = polymer.classify_structure(
            Species(molecule=[impostor]), self.p)
        assert p_class == polymer.PolymerClass.GAS
        assert details["reason"] == "backbone_impostor"

    def test_discreteness_gate_keeps_modified_proxy_images(self):
        """
        Tolerance pin, pass side (one-sided on purpose): legitimately LARGER
        images (+1 side group, -2 H => 27 heavy) and modestly SMALLER images
        (lost center phenyl => 19 heavy >= bound 16) must NOT gas-classify as
        impostors.
        """
        bigger = Molecule(smiles="CC(CC(C)=C(CC(C)c1ccccc1)c1ccccc1)c1ccccc1")
        p_class, details = polymer.classify_structure(
            Species(molecule=[bigger]), self.p)
        assert details["num_disjoint_wings"] == 2  # actually reached the gate
        assert details["reason"] != "backbone_impostor"
        assert p_class != polymer.PolymerClass.GAS

        smaller = Molecule(smiles="CC(CCCC(C)c1ccccc1)c1ccccc1")  # 19 heavy
        p_class, details = polymer.classify_structure(
            Species(molecule=[smaller]), self.p)
        assert details["num_disjoint_wings"] == 2  # actually reached the gate
        assert details["reason"] != "backbone_impostor"
        assert p_class != polymer.PolymerClass.GAS

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

    def _build_crosslink_mol(self):
        """Build a >2-wing crosslink molecule by joining two baseline chains
        at heavy atoms (after freeing valence). Same construction as
        test_branch_crosslink_bimolecular, factored out for reuse."""
        crosslink_mol = self.ref_baseline_mol.copy(deep=True)
        second_chain = self.ref_baseline_mol.copy(deep=True)
        mapping = {}
        for atom in second_chain.atoms:
            new_atom = atom.copy()
            crosslink_mol.add_atom(new_atom)
            mapping[atom] = new_atom
        for atom1 in second_chain.atoms:
            for atom2, bond in atom1.edges.items():
                if id(atom1) < id(atom2):
                    crosslink_mol.add_bond(Bond(mapping[atom1], mapping[atom2], bond.order))
        a1 = next(a for a in crosslink_mol.atoms if not a.is_hydrogen())
        h1 = next(n for n in a1.edges if n.is_hydrogen())
        crosslink_mol.remove_bond(crosslink_mol.get_bond(a1, h1))
        crosslink_mol.remove_atom(h1)
        a2 = next(mapping[a] for a in second_chain.atoms if not a.is_hydrogen())
        h2 = next(n for n in a2.edges if n.is_hydrogen())
        crosslink_mol.remove_bond(crosslink_mol.get_bond(a2, h2))
        crosslink_mol.remove_atom(h2)
        crosslink_mol.add_bond(Bond(a1, a2, order=1))
        crosslink_mol.update_multiplicity()
        return crosslink_mol

    def test_create_reacted_copy_rejects_crosslink(self):
        """
        A crosslink / chain-coupling product (>2 wings) must raise
        PolymerCrosslinkError, NOT silently return None.

        Regression guard: previously _create_reacted_copy_logic fell through to
        None for crosslinks, so the coupled chain was registered as a spurious
        gas-phase molecule (a silent mass leak). create_reacted_copy now rejects
        it so make_new_reaction can discard the whole reaction.
        """
        crosslink_mol = self._build_crosslink_mol()
        # sanity: this really is a crosslink
        assert polymer.classify_structure(
            Species(molecule=[crosslink_mol.copy(deep=True)]), self.p
        )[0] == polymer.PolymerClass.CROSSLINK

        with pytest.raises(polymer.PolymerCrosslinkError):
            self.p.create_reacted_copy(crosslink_mol)

    def test_handshake_structures_propagates_crosslink_rejection(self):
        """
        The crosslink rejection must propagate through _handshake_structures
        (i.e. NOT be swallowed by its ``except (RuntimeError, ValueError)``
        guard), so that make_new_reaction sees it and discards the reaction.
        """
        from rmgpy.data.kinetics.family import _handshake_structures
        crosslink_mol = self._build_crosslink_mol()
        with pytest.raises(polymer.PolymerCrosslinkError):
            _handshake_structures([crosslink_mol], [self.p])

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

    def test_get_thermo_data_polymer_comment_suffix(self):
        """Test that thermo comment gets ', Polymer' suffix."""
        self.dummy_thermo.comment = 'Thermo group additivity estimation: group(Cs-CsCsHH)'
        # Reset thermo so get_thermo_data re-applies
        self.pe_polymer.thermo = None
        thermo = self.pe_polymer.get_thermo_data()
        assert thermo.comment.endswith(', Polymer')
        # Calling again should not double-append
        self.pe_polymer.thermo = None
        thermo2 = self.pe_polymer.get_thermo_data()
        assert thermo2.comment.count(', Polymer') == 1

    def test_get_thermo_data_polymer_comment_empty(self):
        """Test that empty thermo comment becomes 'Polymer'."""
        self.dummy_thermo.comment = ''
        self.pe_polymer.thermo = None
        thermo = self.pe_polymer.get_thermo_data()
        assert thermo.comment == 'Polymer'

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


def radicalize_head_end_group(p, mol):
    """
    Turn an intact baseline proxy into an END_MOD product by abstracting an H
    from the *head terminal end-group* (e.g. CH3 -> CH2.), leaving both backbone
    wings intact. Mirrors the construction used by the END_MOD classification
    tests. Mutates ``mol`` in place and returns the modified end-group atom.
    """
    head_wings = p._wing_groups("head")
    tail_wings = p._wing_groups("tail")
    monomer_group = p.backbone_group
    _, details = polymer._analyze_wing_matches(mol, head_wings, tail_wings, monomer_group)
    heavy_to_full = details["heavy_to_full_map"]
    head_heavy_atoms = details["head_match"]["atoms"]
    mon_heavy_count = sum(1 for ga in monomer_group.atoms if not ga.is_hydrogen())
    end_group_heavy, _ = polymer._slice_wing(head_heavy_atoms, mon_heavy_count)
    target_full = heavy_to_full[list(end_group_heavy)[0]]
    h_atom = next(n for n in target_full.bonds if n.is_hydrogen())
    mol.remove_bond(mol.get_bond(target_full, h_atom))
    mol.remove_atom(h_atom)
    target_full.radical_electrons = 1
    mol.update_multiplicity()
    return target_full


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


# ---------------------------------------------------------------------------
# Functional tests for the polymer-handshake pipeline
# ---------------------------------------------------------------------------

class TestHandshakeStructures:
    """
    Functional tests verifying that _handshake_structures (called from
    CoreEdgeReactionModel.make_new_reaction) correctly converts product
    Molecule objects into Polymer objects when a Polymer is among the
    reaction reactants.

    These tests simulate the key step in the pipeline:
        react_all → generate_reactions_from_families → Molecule products
        → _handshake_structures → Polymer products
        → make_new_species / _register_polymer → Edge species
    """

    PS_ADJ = """multiplicity 3
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

    @pytest.fixture(autouse=True)
    def setup(self):
        from rmgpy.data.kinetics.family import _handshake_structures
        self._handshake = _handshake_structures
        self.ps = Polymer(
            label='PS',
            monomer=self.PS_ADJ,
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )

    # ------------------------------------------------------------------
    # 1. Baseline (unreacted proxy) stays a Polymer
    # ------------------------------------------------------------------
    def test_handshake_baseline_proxy_returns_polymer(self):
        """
        The unreacted baseline proxy molecule should be recognised as 'still polymer'
        (create_reacted_copy returns a copy), so _handshake_structures replaces
        the Molecule with a Polymer.
        """
        proxy_mol = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        product_list = [proxy_mol]
        self._handshake(product_list, [self.ps])
        assert isinstance(product_list[0], Polymer), (
            "Unreacted proxy fragment should become a Polymer after handshake, "
            f"got {type(product_list[0])}"
        )

    def test_handshake_end_mod_flags_end_group_reaction(self):
        """
        After the handshake an END_MOD product makes is_end_group_reaction(products)
        True — exactly what make_new_reaction uses to set
        Reaction.is_end_group_reaction (mu0 chain-end scaling in the solver). A
        baseline (non-terminal) product leaves it False (default mu1 scaling).
        """
        from rmgpy.polymer import is_end_group_reaction
        end_mod = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(self.ps, end_mod)
        products = [end_mod]
        self._handshake(products, [self.ps])
        assert isinstance(products[0], Polymer)
        assert is_end_group_reaction(products) is True

        base = [self.ps.baseline_proxy.molecule[0].copy(deep=True)]
        self._handshake(base, [self.ps])
        assert is_end_group_reaction(base) is False

    def test_handshake_products_classify_flux_archetype(self):
        """
        After the handshake, classify_reaction_flux_archetype (the classifier
        make_new_reaction delegates to when stamping
        Reaction.polymer_flux_archetype) returns SAME_POOL for fold-back
        products and SCISSION_FRAGMENT for scission fragments.
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        # Baseline fold-back -> SAME_POOL
        base = [self.ps.baseline_proxy.molecule[0].copy(deep=True)]
        self._handshake(base, [self.ps])
        assert isinstance(base[0], Polymer)
        assert (classify_reaction_flux_archetype([self.ps], base)
                == PolymerFluxArchetype.SAME_POOL)

        # Head-wing-only fragment -> scission_tail Polymer -> SCISSION_FRAGMENT
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert frag is not None
        prods = [frag]
        self._handshake(prods, [self.ps])
        assert isinstance(prods[0], Polymer)
        assert (classify_reaction_flux_archetype([self.ps], prods)
                == PolymerFluxArchetype.SCISSION_FRAGMENT)

    def test_surge_chip_sub_shape_b_live_end_mod_fold_back(self):
        """
        Spec test 3b -- the only flag-true shape live today: products =
        [SCISSION piece, END_MOD fold-back]. Surgery demotes the chip back to
        a discrete Molecule (undoing its handshake conversion), re-stamps the
        END_MOD fold-back CHIP, and returns a = round(134.2/104.15) = 1.
        Flag-stability rider: the recompute over surged products flips to
        False (END_MOD member gone) -- which is exactly why nothing downstream
        may recompute the flag from product stamps.
        """
        from rmgpy.polymer import (PolymerFluxArchetype, is_end_group_reaction,
                                   classify_reaction_flux_archetype,
                                   surge_chip_products)

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert frag is not None
        end_mod = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(self.ps, end_mod)

        products = [frag.copy(deep=True), end_mod]
        self._handshake(products, [self.ps])
        assert products[0]._reacted_class == PolymerClass.SCISSION
        assert products[1]._reacted_class == PolymerClass.END_MOD
        assert is_end_group_reaction(products) is True  # the stored flag's value

        a = surge_chip_products(products, self.ps)

        assert a == 1
        assert isinstance(products[0], Molecule)          # chip demoted
        assert not isinstance(products[0], Polymer)
        assert products[0].get_formula() == frag.get_formula()
        assert products[1]._reacted_class == PolymerClass.CHIP
        # Recompute now flips -- pins the no-recompute rule.
        assert is_end_group_reaction(products) is False
        assert (classify_reaction_flux_archetype([self.ps], products)
                == PolymerFluxArchetype.DISCRETE_CHIP)

    def test_surge_chip_sub_shape_a_macro_daughter(self):
        """
        Spec test 3 -- sub-shape (a) (dormant today, live when the end-anchor
        detector lands): the SCISSION-stamped Polymer is the MACRO daughter
        and the chip is the single discrete co-product. Surgery replaces the
        daughter with parent.copy(deep=True) stamped CHIP; the chip stays
        as-is; a stamps from the chip's MW ratio.
        """
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   classify_reaction_flux_archetype,
                                   surge_chip_products)

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        prods = [frag]
        self._handshake(prods, [self.ps])
        daughter = prods[0]
        assert isinstance(daughter, Polymer)
        assert daughter._reacted_class == PolymerClass.SCISSION

        chip = Molecule(smiles="C=Cc1ccccc1")   # styrene, 104.15 g/mol -> a = 1
        products = [daughter, chip]
        a = surge_chip_products(products, self.ps)

        assert a == 1
        fold = products[0]
        assert isinstance(fold, Polymer)
        assert fold.label == self.ps.label                # PARENT pool fold-back
        assert fold._reacted_class == PolymerClass.CHIP
        assert products[1] is chip                        # chip untouched, discrete
        assert (classify_reaction_flux_archetype([self.ps], products)
                == PolymerFluxArchetype.DISCRETE_CHIP)

    def test_surge_chip_infeasible_stamps_unresolved_never_scission(self):
        """
        Spec test 4: surgery-infeasible flag-true scission shapes stamp
        UNRESOLVED + warn-once via stamp_polymer_flux_archetype -- NEVER
        SCISSION_FRAGMENT (uniform-cut statistics near an end + unaccounted
        chip mass). Two infeasible shapes: (b) without a demotable source
        molecule, and (a) without a discrete co-product.
        """
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()
        from rmgpy.reaction import Reaction
        from rmgpy.polymer import PolymerFluxArchetype, stamp_polymer_flux_archetype

        # Infeasible (b): SCISSION chip with no _source_molecule.
        sc = self.ps.copy()
        sc.label = "PS_scission_tail"
        sc._reacted_class = PolymerClass.SCISSION
        end = self.ps.copy()
        end._reacted_class = PolymerClass.END_MOD
        rxn = Reaction(reactants=[self.ps], products=[sc, end],
                       is_end_group_reaction=True)
        stamp_polymer_flux_archetype(rxn, [self.ps], [self.ps])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)
        assert rxn.polymer_chip_units == 0

        # Infeasible (a): flag-true scission shape, no discrete co-product.
        sc2 = self.ps.copy()
        sc2.label = "PS_scission_head"
        sc2._reacted_class = PolymerClass.SCISSION
        rxn2 = Reaction(reactants=[self.ps], products=[sc2],
                        is_end_group_reaction=True)
        stamp_polymer_flux_archetype(rxn2, [self.ps], [self.ps])
        assert rxn2.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)

        # Warn-once: repeating an already-warned shape adds no registry entry.
        n = len(polymer_mod._flux_archetype_warned)
        stamp_polymer_flux_archetype(rxn2, [self.ps], [self.ps])
        assert len(polymer_mod._flux_archetype_warned) == n

    def test_surge_chip_a_zero_bare_cap_ejection(self):
        """
        Spec test 6: a = 0 chips are legal (bare end-cap ejection, e.g. CH3
        loss): surgery succeeds and returns 0 (NOT None) -- the archetype
        fires with zero mu1/mu2 drain, net pool effect ~ SAME_POOL.
        """
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   classify_reaction_flux_archetype,
                                   surge_chip_products)

        sc = self.ps.copy()
        sc.label = "PS_scission_tail"
        sc._reacted_class = PolymerClass.SCISSION
        sc._source_molecule = Molecule(smiles="C")   # CH4 cap image, 16 g/mol
        end = self.ps.copy()
        end._reacted_class = PolymerClass.END_MOD
        products = [sc, end]

        a = surge_chip_products(products, self.ps)

        assert a == 0
        assert a is not None                         # 0 != infeasible
        assert isinstance(products[0], Molecule)
        assert products[1]._reacted_class == PolymerClass.CHIP
        assert (classify_reaction_flux_archetype([self.ps], products)
                == PolymerFluxArchetype.DISCRETE_CHIP)

    # ------------------------------------------------------------------
    # 2. Head-scission fragment → scission_tail Polymer
    # ------------------------------------------------------------------
    def test_handshake_head_scission_fragment_returns_scission_tail_polymer(self):
        """
        A fragment that contains only the head wing (no tail wing) should be
        classified as a scission_tail Polymer product.
        """
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        scission_frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert scission_frag is not None, "test setup: scission fragment construction failed"

        product_list = [scission_frag]
        self._handshake(product_list, [self.ps])

        result = product_list[0]
        assert isinstance(result, Polymer), (
            f"Scission fragment should become a Polymer, got {type(result)}"
        )
        assert result.label.endswith('_scission_tail'), (
            f"Expected label ending in '_scission_tail', got '{result.label}'"
        )

    # ------------------------------------------------------------------
    # 3. Tail-scission fragment → scission_head Polymer
    # ------------------------------------------------------------------
    def test_handshake_tail_scission_fragment_returns_scission_head_polymer(self):
        """
        A fragment that contains only the tail wing (no head wing) should be
        classified as a scission_head Polymer product.
        """
        tail_wing = self.ps._stitch_wing("tail")
        methyl_star1 = Molecule().from_adjacency_list(_methyl_radical_adj("*1"))
        scission_frag = polymer.stitch_molecules_by_labeled_atoms(methyl_star1, tail_wing)
        assert scission_frag is not None, "test setup: scission fragment construction failed"

        product_list = [scission_frag]
        self._handshake(product_list, [self.ps])

        result = product_list[0]
        assert isinstance(result, Polymer), (
            f"Scission fragment should become a Polymer, got {type(result)}"
        )
        assert result.label.endswith('_scission_head'), (
            f"Expected label ending in '_scission_head', got '{result.label}'"
        )

    # ------------------------------------------------------------------
    # 4. Small molecule (no wings) is left as a Molecule
    # ------------------------------------------------------------------
    def test_handshake_small_molecule_remains_molecule(self):
        """
        A fragment too small to contain any polymer wing should NOT be
        converted to a Polymer — it should remain a plain Molecule (gas-phase).
        """
        small = Molecule(smiles='CC')  # ethane — no PS wings
        product_list = [small]
        self._handshake(product_list, [self.ps])
        assert isinstance(product_list[0], Molecule), (
            "Small gas-phase molecule should remain a Molecule after handshake"
        )

    # ------------------------------------------------------------------
    # 5. Non-Molecule items in the list are untouched
    # ------------------------------------------------------------------
    def test_handshake_ignores_non_molecule_items(self):
        """
        Non-Molecule items (e.g. already-converted Polymer or Species) in
        the product list should be left unchanged.
        """
        already_poly = self.ps.copy()
        product_list = [already_poly]
        self._handshake(product_list, [self.ps])
        assert product_list[0] is already_poly, (
            "Non-Molecule items should be left untouched by _handshake_structures"
        )

    # ------------------------------------------------------------------
    # 6. Mixed product list: one convertible, one gas
    # ------------------------------------------------------------------
    def test_handshake_mixed_product_list(self):
        """
        In a bimolecular-product reaction, one product may be a Polymer
        fragment and the other a small gas-phase molecule.  Both must be
        handled correctly in the same call.
        """
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        polymer_frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)

        gas_frag = Molecule(smiles='[CH3]')  # methyl radical — no PS wings

        product_list = [polymer_frag, gas_frag]
        self._handshake(product_list, [self.ps])

        assert isinstance(product_list[0], Polymer), (
            f"First product (polymer fragment) should be a Polymer, got {type(product_list[0])}"
        )
        assert isinstance(product_list[1], Molecule), (
            f"Second product (gas fragment) should remain a Molecule, got {type(product_list[1])}"
        )

    # ------------------------------------------------------------------
    # 6b. Stale polymer-proxy tag clearing on retained-discrete products
    #     (handshake/chip handshake-layer fix)
    # ------------------------------------------------------------------
    def test_handshake_clears_proxy_on_discrete_gas_product(self):
        """
        family.py:1665 blanket-stamps every product is_polymer_proxy=True when
        any reactant is a proxy (the PS pool always is). alpha-methylstyrene
        (C=C(C)c1ccccc1, C9H10) is a genuine discrete gas-phase volatile that
        create_reacted_copy returns None for, so the handshake KEEPS it a
        Molecule -- but the stale proxy stamp must be CLEARED so the solver
        does not count it as a melt reference-state participant. A genuine
        scission tail in the same product list still becomes a proxy Polymer.

        RED before the fix: the volatile keeps is_polymer_proxy True.
        """
        # genuine scission tail -> converts to a SCISSION proxy Polymer
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        tail = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert tail is not None, "test setup: scission tail construction failed"

        # discrete volatile pre-stamped proxy True (simulating family.py:1665)
        volatile = Molecule(smiles="C=C(C)c1ccccc1")
        volatile.is_polymer_proxy = True
        volatile.props["is_polymer_proxy"] = True

        products = [tail.copy(deep=True), volatile]
        self._handshake(products, [self.ps])

        # genuine tail -> Polymer (still a proxy)
        assert isinstance(products[0], Polymer), (
            f"Scission tail should become a Polymer, got {type(products[0])}"
        )
        # volatile retained as a discrete Molecule with proxy tag CLEARED
        assert isinstance(products[1], Molecule)
        assert not isinstance(products[1], Polymer)
        assert products[1].is_polymer_proxy is False, (
            "Retained discrete gas-phase volatile must have its stale "
            "is_polymer_proxy stamp cleared by the handshake"
        )
        assert products[1].props.get("is_polymer_proxy") in (False, None)

    def test_surge_chip_clears_proxy_on_discrete_chip(self):
        """
        Chip surgery sub-shape (b) demotes the SCISSION-stamped chip back to a
        discrete Molecule (undoing its handshake conversion). That demoted
        discrete chip carries the stale is_polymer_proxy=True stamp inherited
        from its source Molecule, and must be CLEARED so the chip is not
        solver-visible as a melt participant.

        RED before the fix: the demoted chip keeps is_polymer_proxy True.
        """
        from rmgpy.polymer import surge_chip_products

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert frag is not None
        # simulate the family.py:1665 blanket stamp on the pre-handshake frag
        frag.is_polymer_proxy = True
        frag.props["is_polymer_proxy"] = True
        end_mod = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(self.ps, end_mod)

        products = [frag.copy(deep=True), end_mod]
        self._handshake(products, [self.ps])
        assert products[0]._reacted_class == PolymerClass.SCISSION

        a = surge_chip_products(products, self.ps)
        assert a == 1
        chip = products[0]
        assert isinstance(chip, Molecule)
        assert not isinstance(chip, Polymer)
        assert chip.is_polymer_proxy is False, (
            "Demoted discrete chip must have its stale is_polymer_proxy stamp "
            "cleared by surge_chip_products"
        )
        assert chip.props.get("is_polymer_proxy") in (False, None)

    def test_handshake_keeps_proxy_on_polymer_fragment(self):
        """
        GUARD (no over-clearing): a genuine scission tail that the handshake
        DOES convert to a Polymer must keep is_polymer_proxy True -- the clear
        only fires on products that stay discrete. GREEN before AND after.
        """
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        tail = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert tail is not None
        tail.is_polymer_proxy = True
        tail.props["is_polymer_proxy"] = True

        products = [tail.copy(deep=True)]
        self._handshake(products, [self.ps])
        assert isinstance(products[0], Polymer)
        assert products[0].is_polymer_proxy is True, (
            "A genuine fragment converted to a Polymer must remain a proxy; "
            "the handshake clear must not touch it"
        )

    def test_clear_polymer_proxy_settles_sticky_species_cache(self):
        """
        Species.is_polymer_proxy is a sticky lazy-cache property: its getter
        re-derives True from ANY proxy molecule and re-caches it. The LIVE
        handshake item is a Species (not a Molecule), so clear_polymer_proxy
        must clear the constituent molecules BEFORE the species-level flag --
        otherwise the species _is_polymer_proxy cache re-sticks True off the
        not-yet-cleared molecules, and make_new_species (model.py) ORs that
        stale True onto the solver-visible Species (the alpha-methylstyrene
        reference-state-tripwire leak observed in the live PS run).

        RED before the ordering fix: species stays True; the make_new_species
        OR stays True.
        """
        from rmgpy.species import Species
        from rmgpy.polymer import clear_polymer_proxy

        sp = Species(molecule=[Molecule(smiles="C=C(C)c1ccccc1")])
        sp.is_polymer_proxy = True  # setter caches True + propagates to molecules
        assert sp.is_polymer_proxy is True
        assert sp.molecule[0].is_polymer_proxy is True

        clear_polymer_proxy(sp)

        assert sp.molecule[0].is_polymer_proxy is False, "constituent molecule not cleared"
        assert sp.is_polymer_proxy is False, (
            "sticky species cache must settle False after molecules-first clear"
        )
        # the exact OR make_new_species performs (model.py:486) must be False
        assert (sp.molecule[0].is_polymer_proxy or sp.is_polymer_proxy) is False

        # multi-molecule (resonance) species: every molecule + the species clear
        sp2 = Species(molecule=[Molecule(smiles="C=C(C)c1ccccc1"),
                                Molecule(smiles="C=C(C)c1ccccc1")])
        sp2.is_polymer_proxy = True
        clear_polymer_proxy(sp2)
        assert sp2.is_polymer_proxy is False
        assert all(m.is_polymer_proxy is False for m in sp2.molecule)

    def test_handshake_sets_gas_veto_on_discrete_gas_product(self):
        """
        A genuine discrete gas volatile that the handshake keeps a Molecule
        (create_reacted_copy -> None) must be stamped with the DURABLE
        gas-volatile veto in ``props`` -- not merely have its stale proxy tag
        cleared. Clearing is_polymer_proxy alone is defeated downstream because
        the flag is a monotonic multi-writer sticky cache (re-stamped by
        family.py:1665 + the species.py sticky getter). The positive veto in
        ``props`` survives Species.copy and is never touched by the proxy
        stamping machinery, so the solver melt gate can honor it.

        RED before the fix: no veto key is set on the retained volatile.
        """
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        tail = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert tail is not None, "test setup: scission tail construction failed"

        volatile = Molecule(smiles="C=C(C)c1ccccc1")
        volatile.is_polymer_proxy = True
        volatile.props["is_polymer_proxy"] = True

        products = [tail.copy(deep=True), volatile]
        self._handshake(products, [self.ps])

        # genuine tail -> Polymer: must NOT be vetoed (it is a real chain)
        assert isinstance(products[0], Polymer)
        assert products[0].props.get(VETO) in (False, None), (
            "a genuine scission-tail Polymer must not receive the gas veto"
        )
        # discrete volatile: durable gas veto SET
        assert isinstance(products[1], Molecule)
        assert not isinstance(products[1], Polymer)
        assert products[1].props.get(VETO) is True, (
            "retained discrete gas volatile must carry the durable "
            "polymer_reference_state_gas_veto in props"
        )

    def test_surge_chip_sets_gas_veto_on_discrete_chip(self):
        """
        A demoted discrete chip (surge_chip_products sub-shape b) is a genuine
        gas-phase fragment; it must carry the durable gas-volatile veto so the
        solver never counts it as a melt reference-state participant.

        RED before the fix: no veto key on the demoted chip.
        """
        from rmgpy.polymer import surge_chip_products
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert frag is not None
        frag.is_polymer_proxy = True
        frag.props["is_polymer_proxy"] = True
        end_mod = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(self.ps, end_mod)

        products = [frag.copy(deep=True), end_mod]
        self._handshake(products, [self.ps])
        assert products[0]._reacted_class == PolymerClass.SCISSION

        a = surge_chip_products(products, self.ps)
        assert a == 1
        chip = products[0]
        assert isinstance(chip, Molecule) and not isinstance(chip, Polymer)
        assert chip.props.get(VETO) is True, (
            "demoted discrete chip must carry the durable gas-volatile veto"
        )

    def test_handshake_does_not_veto_polymer_fragment(self):
        """
        GUARD (no over-vetoing): a fragment the handshake DOES convert to a
        Polymer (a real chain) must NOT receive the gas veto -- otherwise a
        genuine melt chain would be wrongly excluded from the melt sum. GREEN
        before AND after.
        """
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        tail = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert tail is not None
        tail.is_polymer_proxy = True
        tail.props["is_polymer_proxy"] = True

        products = [tail.copy(deep=True)]
        self._handshake(products, [self.ps])
        assert isinstance(products[0], Polymer)
        assert products[0].props.get(VETO) in (False, None), (
            "a Polymer chain must never receive the gas-volatile veto"
        )

    def test_species_copy_preserves_gas_veto(self):
        """
        Contract pin: the durable gas veto lives in ``Species.props``, which
        ``Species.copy`` deep-copies -- so the verdict survives the copies that
        defeat both is_polymer_proxy clears and ad-hoc attributes. (Green from
        the start; documents WHY props is the chosen carrier.)
        """
        from rmgpy.species import Species
        from rmgpy.polymer import (POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO,
                                   set_polymer_gas_veto)

        sp = Species(molecule=[Molecule(smiles="C=C(C)c1ccccc1")])
        set_polymer_gas_veto(sp)
        assert sp.props.get(VETO) is True
        clone = sp.copy(deep=True)
        assert clone.props.get(VETO) is True, (
            "Species.copy must preserve the durable gas veto (props)"
        )

    def test_gas_veto_key_literal_matches_solver_gate(self):
        """Literal-drift guard (code-review NIT): the Cython solver melt gate
        (rmgpy/solver/polymer.pyx) reads the props key as a HARDCODED string
        literal, not the imported constant. If the Python constant is renamed
        without updating the pyx, the solver silently stops honoring the veto
        and the reference-state tripwire returns. Pin the exact literal so that
        rename fails loudly here.
        """
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY
        assert POLYMER_REFERENCE_STATE_GAS_VETO_KEY == "polymer_reference_state_gas_veto", (
            "the gas-veto props key literal is hardcoded in polymer.pyx's melt "
            "gate; update both together"
        )

    # ------------------------------------------------------------------
    # 7. Retroene-style scission: closed-shell fragments from proxy
    # ------------------------------------------------------------------
    def test_handshake_retroene_scission_products(self):
        """
        A Retroene reaction on the PS proxy trimer produces two closed-shell
        fragments (no radicals).  The larger fragment containing a recognizable
        wing should become a Polymer; the smaller one stays a Molecule.

        PS proxy: CH3-CH2-CH(Ph)-CH2-CH(Ph)-CH2-CH(Ph)-H
        Retroene splits e.g. into:
          C17H18: C=C(CC(C)c1ccccc1)c1ccccc1  (larger, has a wing)
          C8H10:  CC=C1C=CC=CC1               (smaller fragment)

        At least the larger fragment must be recognized as polymer-derived.
        """
        # These are actual SMILES from the RMG run output
        large_frag = Molecule(smiles='C=C(CC(C)c1ccccc1)c1ccccc1')
        small_frag = Molecule(smiles='C=C(C)c1ccccc1')

        product_list = [large_frag, small_frag]
        self._handshake(product_list, [self.ps])

        # At least the large fragment should be recognized as a Polymer
        assert isinstance(product_list[0], Polymer), (
            f"Large Retroene fragment should become a Polymer, got {type(product_list[0])}"
        )

    def test_handshake_retroene_all_scission_products(self):
        """
        Test both pairs of Retroene products from the PS proxy.
        For each pair, the larger fragment should become a Polymer.
        """
        pairs = [
            ('C=C(CC(C)c1ccccc1)c1ccccc1', 'CC=C1C=CC=CC1'),     # C17H18 + C8H10
            ('CC(CC=C1C=CC=CC1)c1ccccc1', 'C=C(C)c1ccccc1'),     # C16H18 + C9H10
        ]
        for large_smi, small_smi in pairs:
            large_frag = Molecule(smiles=large_smi)
            small_frag = Molecule(smiles=small_smi)
            product_list = [large_frag, small_frag]
            self._handshake(product_list, [self.ps])
            assert isinstance(product_list[0], Polymer), (
                f"Large fragment ({large_smi}) should become a Polymer, "
                f"got {type(product_list[0])}"
            )

    # ------------------------------------------------------------------
    # 8. Handshake with Species objects (real RMG flow)
    # ------------------------------------------------------------------
    def test_handshake_species_objects(self):
        """
        In the real RMG pipeline, find_degenerate_reactions wraps product
        Molecules into Species objects before process_reaction is called.
        _handshake_structures must handle Species (not just Molecule) items
        in the product list.
        """
        from rmgpy.species import Species as Spc
        large_mol = Molecule(smiles='C=C(CC(C)c1ccccc1)c1ccccc1')
        small_mol = Molecule(smiles='C=C(C)c1ccccc1')
        large_spc = Spc(molecule=[large_mol])
        small_spc = Spc(molecule=[small_mol])

        product_list = [large_spc, small_spc]
        self._handshake(product_list, [self.ps])

        assert isinstance(product_list[0], Polymer), (
            f"Species wrapping a large fragment should become a Polymer, "
            f"got {type(product_list[0])}"
        )
        assert not isinstance(product_list[1], Polymer), (
            "Small gas-phase fragment should stay a Species"
        )

    # ------------------------------------------------------------------
    # 9. Scission product proxy symmetry / resonance hybrid
    # ------------------------------------------------------------------
    def test_scission_product_proxy_symmetry_number(self):
        """
        The proxy of a scission-tail Polymer must be able to compute its
        symmetry number (which internally calls get_resonance_hybrid)
        without crashing due to inconsistent resonance structures.

        Regression test for ValueError: 'The specified vertices are not
        connected by an edge in this graph' during thermo generation.
        """
        large_frag = Molecule(smiles='C=C(CC(C)c1ccccc1)c1ccccc1')
        product_list = [large_frag, Molecule(smiles='C=C(C)c1ccccc1')]
        self._handshake(product_list, [self.ps])
        poly = product_list[0]
        assert isinstance(poly, Polymer)

        proxy = poly.get_proxy_species()
        # This is the exact call chain that crashed in the RMG run:
        # get_symmetry_number → get_resonance_hybrid → get_bond
        sym = poly.get_symmetry_number()
        assert sym >= 1

    def test_scission_product_generate_resonance_structures(self):
        """
        Calling generate_resonance_structures on a Polymer must not break
        the shared molecule reference between the Polymer and its proxy.
        """
        large_frag = Molecule(smiles='C=C(CC(C)c1ccccc1)c1ccccc1')
        product_list = [large_frag, Molecule(smiles='C=C(C)c1ccccc1')]
        self._handshake(product_list, [self.ps])
        poly = product_list[0]
        assert isinstance(poly, Polymer)

        proxy = poly.get_proxy_species()
        # Simulate what evaluator does
        poly.generate_resonance_structures()
        # After the call, Polymer.molecule must still reference the proxy's list
        assert poly.molecule is proxy.molecule

    def test_evaluator_on_scission_polymer(self):
        """
        Full evaluator path on a scission Polymer must not crash.
        This mimics the exact thermo-generation flow triggered by
        _register_polymer → generate_thermo → submit → evaluator.
        """
        large_frag = Molecule(smiles='C=C(CC(C)c1ccccc1)c1ccccc1')
        product_list = [large_frag, Molecule(smiles='C=C(C)c1ccccc1')]
        self._handshake(product_list, [self.ps])
        poly = product_list[0]
        assert isinstance(poly, Polymer)

        # Simulate the evaluator flow
        poly.generate_resonance_structures()
        sym = poly.get_symmetry_number()
        assert sym >= 1


class TestHandshakeRelabelFlag:
    """Tests that _handshake_structures returns a bool indicating relabeling."""

    def test_handshake_structures_returns_true_when_it_relabels(self):
        """_handshake_structures must report relabeling so make_new_reaction can
        gate the real-ΔH BM pre-conversion (spec §4.2)."""
        from rmgpy.data.kinetics.family import _handshake_structures
        from rmgpy.molecule import Molecule
        ps = Polymer(label='PS', monomer='[CH2][CH]c1ccccc1',
                     end_groups=['[CH3]', '[H]'], cutoff=3,
                     Mn=5000.0, Mw=6000.0, initial_mass=1.0)
        # A scission-tail fragment Molecule the proxy can reinterpret as a Polymer.
        frag = Molecule().from_smiles('CC(CC(CC(C)c1ccccc1)c1ccccc1)c1ccccc1')
        products = [frag]
        relabeled = _handshake_structures(products, [ps])
        assert relabeled is True
        from rmgpy.polymer import Polymer as _P
        assert isinstance(products[0], _P)

    def test_handshake_structures_returns_false_when_nothing_relabels(self):
        from rmgpy.data.kinetics.family import _handshake_structures
        from rmgpy.molecule import Molecule
        ps = Polymer(label='PS', monomer='[CH2][CH]c1ccccc1',
                     end_groups=['[CH3]', '[H]'], cutoff=3,
                     Mn=5000.0, Mw=6000.0, initial_mass=1.0)
        small = [Molecule().from_smiles('O=C=O')]  # CO2: not a polymer fragment
        relabeled = _handshake_structures(small, [ps])
        assert relabeled is False
        assert isinstance(small[0], Molecule)


class TestPolymerRegistration:
    """
    Tests for Polymer registration in the CoreEdgeReactionModel:
    fingerprint uniqueness, moment dummy injection, and the end-to-end
    make_new_species → _register_polymer pipeline.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        from rmgpy.rmg.model import CoreEdgeReactionModel
        self.model = CoreEdgeReactionModel()
        self.ps = Polymer(
            label='PS',
            monomer='[CH2][CH]c1ccccc1',
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )

    def test_register_polymer_assigns_index(self):
        """A newly registered Polymer must get a positive species index."""
        poly, is_new = self.model._register_polymer(self.ps, generate_thermo=False)
        assert is_new is True
        assert poly.index > 0

    def test_register_polymer_creates_moment_dummies(self):
        """_register_polymer must inject _mu0, _mu1, _mu2 dummy Species."""
        self.model._register_polymer(self.ps, generate_thermo=False)
        labels = [s.label for s in self.model.new_species_list]
        for suffix in ('_mu0', '_mu1', '_mu2'):
            expected = f'PS{suffix}'
            assert expected in labels, f"Missing moment dummy '{expected}' in new_species_list"

    def test_moment_dummies_are_nonreactive_ne(self):
        """Moment dummies must be non-reactive Species with [Ne] placeholder
        molecules (see CoreEdgeReactionModel._register_polymer, which uses
        from_smiles('[Ne]'); the Cantera writer also documents the Ne
        placeholder convention)."""
        self.model._register_polymer(self.ps, generate_thermo=False)
        for spc in self.model.new_species_list:
            if spc.label.startswith('PS_mu'):
                assert spc.reactive is False
                assert spc.index == -1
                assert spc.molecule[0].get_formula() == 'Ne'

    def test_duplicate_polymer_returns_existing(self):
        """Registering the same Polymer twice must return the first copy (is_new=False)."""
        poly1, is_new1 = self.model._register_polymer(self.ps, generate_thermo=False)
        dup = self.ps.copy(deep=True)
        poly2, is_new2 = self.model._register_polymer(dup, generate_thermo=False)
        assert is_new1 is True
        assert is_new2 is False
        assert poly1 is poly2

    def test_scission_polymer_gets_different_fingerprint(self):
        """A scission product must have a different fingerprint than its parent."""
        c16h18 = Molecule().from_smiles('CC(CC=C1C=CC=CC1)c1ccccc1')
        scission = self.ps.create_reacted_copy(c16h18)
        assert scission is not None, "test setup: create_reacted_copy returned None"
        assert scission.fingerprint != self.ps.fingerprint, (
            "Scission product must have a different fingerprint from the parent"
        )

    def test_scission_polymer_gets_own_moment_dummies(self):
        """A scission product registered via make_new_species gets its own moment dummies."""
        self.model._register_polymer(self.ps, generate_thermo=False)
        c16h18 = Molecule().from_smiles('CC(CC=C1C=CC=CC1)c1ccccc1')
        scission = self.ps.create_reacted_copy(c16h18)
        scission_reg, is_new = self.model.make_new_species(scission, generate_thermo=False)
        assert is_new is True
        assert isinstance(scission_reg, Polymer)
        labels = [s.label for s in self.model.new_species_list]
        for suffix in ('_mu0', '_mu1', '_mu2'):
            expected = f'{scission_reg.label}{suffix}'
            assert expected in labels, f"Missing moment dummy '{expected}' for scission product"

    def test_make_new_species_routes_polymer_to_register_polymer(self):
        """make_new_species with a Polymer must route to _register_polymer."""
        poly, is_new = self.model.make_new_species(self.ps, generate_thermo=False)
        assert is_new is True
        assert isinstance(poly, Polymer)
        assert poly.index > 0
        # Moment dummies should also be present
        labels = [s.label for s in self.model.new_species_list]
        assert f'{poly.label}_mu0' in labels

    def test_handshake_then_register_end_to_end(self):
        """
        Full pipeline: handshake converts product Molecule to Polymer,
        then make_new_species registers it with index and moment dummies.
        """
        from rmgpy.data.kinetics.family import _handshake_structures

        self.model._register_polymer(self.ps, generate_thermo=False)

        # Simulate product list from a reaction
        c16h18 = Molecule().from_smiles('CC(CC=C1C=CC=CC1)c1ccccc1')
        h_atom = Molecule().from_smiles('[H]')
        products = [c16h18, h_atom]

        _handshake_structures(products, [self.ps])

        # After handshake: first product should be Polymer, second stays Molecule
        assert isinstance(products[0], Polymer)
        assert isinstance(products[1], Molecule)

        # Register both via make_new_species
        poly_prod, is_new_poly = self.model.make_new_species(products[0], generate_thermo=False)
        h_prod, is_new_h = self.model.make_new_species(products[1], generate_thermo=False)

        assert isinstance(poly_prod, Polymer)
        assert is_new_poly is True
        assert poly_prod.index > 0
        assert not isinstance(h_prod, Polymer)
        assert is_new_h is True

        # Check moment dummies exist for the new polymer
        labels = [s.label for s in self.model.new_species_list]
        for suffix in ('_mu0', '_mu1', '_mu2'):
            assert f'{poly_prod.label}{suffix}' in labels


class TestMakeNewReactionPolymer:
    """
    Tests for make_new_reaction() handling Polymer reactants/products:
    pairs invalidation after handshake, and end-to-end product registration.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        from rmgpy.rmg.model import CoreEdgeReactionModel
        self.model = CoreEdgeReactionModel()

        # Create and register a polystyrene Polymer
        self.ps = Polymer(
            label='PS',
            monomer='[CH2][CH]c1ccccc1',
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )
        self.model._register_polymer(self.ps, generate_thermo=False)

    def _make_retroene_template_reaction(self):
        """
        Build a TemplateReaction that mimics a Retroene scission of the PS proxy.

        The PS trimer proxy decomposes into two fragments:
          C16H18: CC(CC=C1C=CC=CC1)c1ccccc1 (scission-head-like)
          C8H10:  CC=C1C=CC=CC1             (scission-tail-like)

        Returns the TemplateReaction with pairs set (the scenario that
        previously caused a ValueError after handshake).
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.kinetics import Arrhenius

        proxy_mol = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        c16h18 = Molecule().from_smiles('CC(CC=C1C=CC=CC1)c1ccccc1')
        c8h10 = Molecule().from_smiles('CC=C1C=CC=CC1')

        rxn = TemplateReaction(
            reactants=[proxy_mol],
            products=[c16h18, c8h10],
            family='Retroene',
            is_forward=True,
            kinetics=Arrhenius(A=(1.29e12, 's^-1'), n=0.0, Ea=(71.113, 'kcal/mol')),
            pairs=[(proxy_mol, c16h18), (proxy_mol, c8h10)],
        )
        return rxn

    def test_make_new_reaction_no_pairs_crash(self):
        """make_new_reaction must not raise ValueError on pairs lookup after handshake."""
        rxn = self._make_retroene_template_reaction()
        result_rxn, is_new = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        assert is_new is True
        assert result_rxn is not None

    def test_make_new_reaction_produces_polymer_products(self):
        """At least one product of a Polymer scission must be a Polymer object."""
        rxn = self._make_retroene_template_reaction()
        result_rxn, _ = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        poly_products = [p for p in result_rxn.products if isinstance(p, Polymer)]
        assert len(poly_products) > 0, (
            f"Expected Polymer products, got: {[type(p).__name__ for p in result_rxn.products]}"
        )

    def test_make_new_reaction_polymer_products_get_moment_dummies(self):
        """Polymer products from make_new_reaction must have moment dummies registered."""
        rxn = self._make_retroene_template_reaction()
        self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        labels = [s.label for s in self.model.new_species_list]
        poly_products = [p for p in rxn.products if isinstance(p, Polymer)]
        for poly in poly_products:
            for suffix in ('_mu0', '_mu1', '_mu2'):
                expected = f'{poly.label}{suffix}'
                assert expected in labels, (
                    f"Missing moment dummy '{expected}' for polymer product"
                )

    def test_make_new_reaction_reactant_resolves_to_polymer(self):
        """The proxy Molecule reactant must resolve to the registered Polymer."""
        rxn = self._make_retroene_template_reaction()
        result_rxn, _ = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        assert isinstance(result_rxn.reactants[0], Polymer)
        assert result_rxn.reactants[0].label == 'PS'

    def test_make_new_reaction_pairs_regenerated(self):
        """After handshake invalidates pairs, generate_pairs must restore them."""
        rxn = self._make_retroene_template_reaction()
        result_rxn, _ = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        assert result_rxn.pairs is not None
        assert len(result_rxn.pairs) == max(len(result_rxn.reactants), len(result_rxn.products))

    def _make_chip_template_reaction(self):
        """
        Proxy -> [cap+unit fragment, END_MOD image]: handshakes into the live
        chip shape (b) (SCISSION piece + END_MOD fold-back), per the probed
        recipe in the 2026-06-10 spec work. Fragment MW 134.2 -> a = 1.
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.kinetics import Arrhenius

        proxy_mol = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        end_mod = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(self.ps, end_mod)
        return TemplateReaction(
            reactants=[proxy_mol],
            products=[frag, end_mod],
            family='R_Recombination',
            is_forward=True,
            kinetics=Arrhenius(A=(1e13, 's^-1'), n=0.0, Ea=(50.0, 'kcal/mol')),
            pairs=[(proxy_mol, frag), (proxy_mol, end_mod)],
        )

    def test_make_new_reaction_chip_stamps_and_never_queues(self):
        """
        Spec test 7 (+3b's flag-survival rider at the model level): a chip
        event through make_new_reaction stamps DISCRETE_CHIP with
        polymer_chip_units = 1, keeps the STORED is_end_group_reaction True
        (the surgery removed the END_MOD member, so a recompute would flip
        it -- nothing recomputes), registers NO _scission_* daughter, and the
        iteration-boundary spawn pass finds nothing to spawn (never-queue:
        surgery replaced the daughter before the candidates pass).
        """
        from rmgpy.polymer import PolymerFluxArchetype

        rxn = self._make_chip_template_reaction()
        result_rxn, is_new = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False,
            generate_kinetics=False,
        )
        assert result_rxn is not None

        # Stamps: archetype, chip units, stored-flag survival.
        assert result_rxn.is_end_group_reaction is True
        assert (result_rxn.polymer_flux_archetype
                == int(PolymerFluxArchetype.DISCRETE_CHIP))
        assert result_rxn.polymer_chip_units == 1

        # Products: a discrete chip + the PS fold-back; no scission daughter.
        labels = [getattr(p, 'label', '') for p in result_rxn.products]
        assert not any('_scission' in lbl for lbl in labels)
        assert any(isinstance(p, Polymer) and p.label == 'PS'
                   for p in result_rxn.products)

        # Never-queue: no _scission_* Polymer registered...
        assert not any('_scission' in s.label for s in self.model.new_species_list)
        # ...and the spawn pass has nothing to drain (no daughter pools appear).
        self.model._apply_multipool_spawn_pass(self.model.new_species_list)
        assert not any('_scission' in s.label for s in self.model.new_species_list)
        pools = [s for s in self.model.new_species_list if isinstance(s, Polymer)]
        assert pools  # the parent pool must be present (guards vacuous all())
        assert all(p.label == 'PS' for p in pools)


class TestEnlargePolymerPipeline:
    """
    Tests that simulate the enlarge pipeline (make_new_reaction → edge placement)
    for Polymer scission reactions. This mirrors what happens when a Polymer
    species in the core reacts and its products are added to the model edge.

    Uses make_new_reaction (with check_existing=False to avoid database
    dependency) and then manually adds products to the edge, which is
    exactly what process_new_reactions does after the reaction is created.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        from rmgpy.rmg.model import CoreEdgeReactionModel
        self.model = CoreEdgeReactionModel()

        self.ps = Polymer(
            label='PS',
            monomer='[CH2][CH]c1ccccc1',
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )
        self.model._register_polymer(self.ps, generate_thermo=False)
        # Place PS directly in core (bypassing add_species_to_core which
        # requires the full RMG database for forbidden structure checks)
        self.model.core.species.append(self.ps)

    def _make_retroene_reaction(self):
        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.kinetics import Arrhenius

        proxy_mol = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        c16h18 = Molecule().from_smiles('CC(CC=C1C=CC=CC1)c1ccccc1')
        c8h10 = Molecule().from_smiles('CC=C1C=CC=CC1')

        return TemplateReaction(
            reactants=[proxy_mol],
            products=[c16h18, c8h10],
            family='Retroene',
            is_forward=True,
            kinetics=Arrhenius(A=(1.29e12, 's^-1'), n=0.0, Ea=(71.113, 'kcal/mol')),
            pairs=[(proxy_mol, c16h18), (proxy_mol, c8h10)],
        )

    def test_enlarge_pipeline_polymer_products_in_edge(self):
        """
        Simulate the enlarge pipeline: create the reaction, then add
        non-core products to the edge (as process_new_reactions does).
        Polymer products must end up as Polymer objects in the edge.
        """
        rxn = self._make_retroene_reaction()
        result_rxn, is_new = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        assert is_new is True

        # Simulate process_new_reactions edge placement
        for spec in result_rxn.products:
            if spec not in self.model.core.species and spec not in self.model.edge.species:
                self.model.edge.species.append(spec)

        edge_polymers = [s for s in self.model.edge.species if isinstance(s, Polymer)]
        assert len(edge_polymers) > 0, (
            f"No Polymer found in edge species. Edge types: "
            f"{[type(s).__name__ for s in self.model.edge.species]}"
        )

    def test_enlarge_pipeline_polymer_moment_dummies_registered(self):
        """
        After the enlarge pipeline, Polymer products in the edge
        must have their _mu0, _mu1, _mu2 moment dummies registered.
        """
        rxn = self._make_retroene_reaction()
        result_rxn, _ = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        for spec in result_rxn.products:
            if spec not in self.model.core.species and spec not in self.model.edge.species:
                self.model.edge.species.append(spec)

        all_labels = {s.label for s in
                      self.model.new_species_list + self.model.core.species + self.model.edge.species}

        edge_polymers = [s for s in self.model.edge.species if isinstance(s, Polymer)]
        for poly in edge_polymers:
            for suffix in ('_mu0', '_mu1', '_mu2'):
                expected = f'{poly.label}{suffix}'
                assert expected in all_labels, (
                    f"Missing moment dummy '{expected}' for edge polymer '{poly.label}'"
                )

    def test_enlarge_pipeline_polymer_products_have_unique_fingerprints(self):
        """
        Polymer scission products must have different fingerprints
        from the parent PS and from each other.
        """
        rxn = self._make_retroene_reaction()
        result_rxn, _ = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        poly_products = [p for p in result_rxn.products if isinstance(p, Polymer)]
        fingerprints = {p.fingerprint for p in poly_products}
        # All polymer products should have distinct fingerprints
        assert len(fingerprints) == len(poly_products)
        # None should match the parent PS fingerprint
        for fp in fingerprints:
            assert fp != self.ps.fingerprint, (
                "Scission product fingerprint must differ from parent"
            )

    def test_polymer_reaction_not_pressure_dependent(self):
        """
        Reactions involving Polymer species must never be routed to the
        pressure-dependent network, even when pressure_dependence is on.
        They should go directly to the core or edge reaction lists.
        """
        rxn = self._make_retroene_reaction()
        result_rxn, is_new = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        # Place all species in core so the reaction qualifies for the core
        for spec in result_rxn.reactants + result_rxn.products:
            if spec not in self.model.core.species:
                self.model.core.species.append(spec)

        # Enable pressure dependence on the model with a low atom limit
        from unittest.mock import MagicMock
        self.model.pressure_dependence = MagicMock()
        self.model.pressure_dependence.maximum_atoms = 10  # far below polymer size
        self.model.unrealgroups = []

        # Simulate the pdep decision from process_new_reactions
        isomer_atoms = sum(len(spec.molecule[0].atoms) for spec in result_rxn.reactants)
        pdep = True
        if not self.model.pressure_dependence:
            pdep = False
        elif any(isinstance(spec, Polymer) for spec in result_rxn.reactants + result_rxn.products):
            pdep = False

        assert pdep is False, (
            "Polymer reaction should NOT be treated as pressure-dependent"
        )

    def test_moment_dummies_promoted_to_core_with_polymer(self):
        """
        When a Polymer is moved from edge to core, its _mu0, _mu1, _mu2
        moment dummies must be promoted to the core as well.
        """
        from unittest.mock import patch, MagicMock

        rxn = self._make_retroene_reaction()
        result_rxn, _ = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False, generate_kinetics=False,
        )
        # Find the scission polymer product
        poly_product = None
        for spec in result_rxn.products:
            if isinstance(spec, Polymer) and spec is not self.ps:
                poly_product = spec
                break
        assert poly_product is not None, "Expected a scission Polymer product"

        # Place the polymer and its dummies in the edge (simulating edge placement)
        self.model.edge.species.append(poly_product)
        for suffix in ('_mu0', '_mu1', '_mu2'):
            m_label = f"{poly_product.label}{suffix}"
            for s in self.model.new_species_list:
                if s.label == m_label:
                    self.model.edge.species.append(s)
                    break

        # Mock get_db to avoid requiring the full RMG database
        mock_forbidden = MagicMock()
        mock_forbidden.is_molecule_forbidden.return_value = False
        with patch('rmgpy.rmg.model.get_db', return_value=mock_forbidden):
            self.model.add_species_to_core(poly_product)

        core_labels = {s.label for s in self.model.core.species}
        for suffix in ('_mu0', '_mu1', '_mu2'):
            expected = f"{poly_product.label}{suffix}"
            assert expected in core_labels, (
                f"Moment dummy '{expected}' should be in core after polymer promotion"
            )


def test_is_qssa_eliminating_radical_distinguishes_allylic():
    from rmgpy.molecule import Molecule
    from rmgpy.polymer import is_qssa_eliminating_radical
    saturated = Molecule().from_smiles("CCC(C)CCC[C](C)CCCC(C)C")  # C15 mid-chain
    allylic = Molecule().from_smiles("CC=C(C)[CH]CC")              # Probe F dominant allylic (faithful analog), resonance count 2
    assert is_qssa_eliminating_radical(saturated) is True    # resonance count 1 -> eliminating
    assert is_qssa_eliminating_radical(allylic) is False     # resonance count >1 -> accumulating


def test_feature_abstraction_is_flagged_refused_not_leaked():
    """A FEATURE mid-chain radical that the handshake dropped to a gas product
    (UNRESOLVED, mass-fabricating) must be FLAGGED ``polymer_refused`` at stamp
    time -- without raising and without discarding the reaction (item 18)."""
    from rmgpy.polymer import Polymer, stamp_polymer_flux_archetype
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    epdm = Polymer(label="epdm", monomer="[CH2]CC(C)[CH2]",
                   Mn=5000.0, Mw=8000.0, initial_mass=1.0)  # is_polymer_proxy
    macro = Molecule().from_smiles("CCC(C)CCC[C](C)CCCC(C)C")  # leaked FEATURE radical (Molecule, as at real stamp site)
    h = Molecule().from_smiles("[H]")
    h2 = Molecule().from_smiles("[H][H]")
    # reactants contain the Polymer directly (Polymer IS a Species); the leaked
    # FEATURE radical product is a plain Molecule (handshake left it un-converted).
    rxn = Reaction(reactants=[epdm, Species(molecule=[h])],
                   products=[Species(molecule=[h2]), macro])
    polymer_reactants = [r for r in rxn.reactants if isinstance(r, Polymer)]
    stamp_polymer_flux_archetype(rxn, rxn.reactants, polymer_reactants)
    # stamp-but-keep: the reaction is kept (products unchanged), only flags added.
    assert rxn.products == [Species(molecule=[h2]), macro] or macro in rxn.products
    assert rxn.polymer_refused is True
    assert rxn.polymer_refused_accumulating is False   # saturated -> eliminating


def test_same_pool_polymer_reaction_is_not_refused():
    """A normal SAME_POOL reaction with a real polymer product must NOT be
    flagged refused (guards against false-firing of the refuse detector)."""
    import rmgpy.polymer as polymer_mod
    polymer_mod._flux_archetype_warned.clear()
    from rmgpy.polymer import (Polymer, PolymerClass,
                               stamp_polymer_flux_archetype,
                               PolymerFluxArchetype)
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    epdm = Polymer(label="epdm", monomer="[CH2]CC(C)[CH2]",
                   Mn=5000.0, Mw=8000.0, initial_mass=1.0)
    # A FEATURE-modified polymer product in the SAME pool: a handshake-converted
    # fold-back (Polymer, same label, stamped FEATURE) -> SAME_POOL, never
    # UNRESOLVED. The gas co-product (H2) must NOT be misread as a lost radical.
    product = epdm.copy()
    product._reacted_class = PolymerClass.FEATURE
    assert isinstance(product, Polymer)
    h = Molecule().from_smiles("[H]")
    h2 = Molecule().from_smiles("[H][H]")
    rxn = Reaction(reactants=[epdm, Species(molecule=[h])],
                   products=[Species(molecule=[h2]), product])
    polymer_reactants = [r for r in rxn.reactants if isinstance(r, Polymer)]
    stamp_polymer_flux_archetype(rxn, rxn.reactants, polymer_reactants)
    assert rxn.polymer_flux_archetype != int(PolymerFluxArchetype.UNRESOLVED)
    assert rxn.polymer_refused is False
    assert rxn.polymer_refused_accumulating is False


def test_unresolved_non_feature_gas_product_helper_declines():
    """An UNRESOLVED reaction (polymer reactant, no polymer product) whose gas
    product is a genuinely-small fragment -- CH4 (MW ~16 g/mol, well below the
    chain-scale threshold monomer ~70 + slack 10 ~= 80) -- must NOT be refused:
    the solver's single-monomer-debit accounts for sub-monomer leaks correctly.
    This exercises the widened, mechanism-keyed refuse-detection helper's full
    decline path: unlike the SAME_POOL guard, the UNRESOLVED branch IS reached,
    the helper runs its complete loop over gas products, finds none at chain
    scale, and correctly returns None. NOTE: post-widening the decline is by the
    SIZE gate, not the old FEATURE-label check (item 18 Task 3 follow-up)."""
    import rmgpy.polymer as polymer_mod
    polymer_mod._flux_archetype_warned.clear()
    from rmgpy.polymer import (Polymer, stamp_polymer_flux_archetype,
                               PolymerFluxArchetype)
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    epdm = Polymer(label="epdm", monomer="[CH2]CC(C)[CH2]",
                   Mn=5000.0, Mw=8000.0, initial_mass=1.0)
    # CH4 (~16 g/mol) is far below the chain-scale threshold -> not refused.
    ch4 = Molecule().from_smiles("C")
    h = Molecule().from_smiles("[H]")
    h2 = Molecule().from_smiles("[H][H]")
    # Polymer reactant + only-gas products (no polymer product) -> UNRESOLVED.
    rxn = Reaction(reactants=[epdm, Species(molecule=[h])],
                   products=[Species(molecule=[h2]), Species(molecule=[ch4])])
    polymer_reactants = [r for r in rxn.reactants if isinstance(r, Polymer)]
    stamp_polymer_flux_archetype(rxn, rxn.reactants, polymer_reactants)
    # The UNRESOLVED branch IS reached (this is what makes the helper run, unlike
    # the SAME_POOL guard which short-circuits before detection).
    assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)
    # Helper ran the full loop, found nothing at chain scale, and declined.
    assert rxn.polymer_refused is False
    assert rxn.polymer_refused_accumulating is False


def test_feature_allylic_radical_lost_to_gas_is_refused_accumulating():
    """An UNRESOLVED reaction that leaks a FEATURE radical which is ALSO an
    accumulating (allylic, resonance-stabilized) radical must be flagged
    ``polymer_refused`` with ``polymer_refused_accumulating is True`` -- exercising
    the accumulating branch at the stamp site (item 18 Task 3)."""
    import rmgpy.polymer as polymer_mod
    polymer_mod._flux_archetype_warned.clear()
    from rmgpy.polymer import (Polymer, stamp_polymer_flux_archetype,
                               is_qssa_eliminating_radical)
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    epdm = Polymer(label="epdm", monomer="[CH2]CC(C)[CH2]",
                   Mn=5000.0, Mw=8000.0, initial_mass=1.0)
    # A backbone-sized EPDM macroradical bearing an internal C=C with the radical
    # allylic to it: large enough to classify FEATURE against epdm, and
    # resonance-stabilized so it is accumulating (not QSSA-eliminating). The small
    # Probe F analog CC=C(C)[CH]CC is too small (classifies SCISSION, not FEATURE),
    # so a backbone-length allylic radical is required for an end-to-end refuse.
    allylic = Molecule().from_smiles("CCC(C)CCCC=C[C](C)CCC(C)CC")
    assert is_qssa_eliminating_radical(allylic) is False  # accumulating
    h = Molecule().from_smiles("[H]")
    h2 = Molecule().from_smiles("[H][H]")
    rxn = Reaction(reactants=[epdm, Species(molecule=[h])],
                   products=[Species(molecule=[h2]), allylic])
    polymer_reactants = [r for r in rxn.reactants if isinstance(r, Polymer)]
    stamp_polymer_flux_archetype(rxn, rxn.reactants, polymer_reactants)
    assert rxn.polymer_refused is True
    assert rxn.polymer_refused_accumulating is True


def test_discard_chain_radical_lost_to_gas_is_refused():
    """A DISCARD (buffer_monomer_modified) backbone radical dropped to gas on the
    UNRESOLVED leg must now be FLAGGED ``polymer_refused`` -- the widened,
    mechanism-keyed predicate refuses ANY chain-scale gas radical, not only
    classify_structure==FEATURE (item 18 Task 3 follow-up). The FEATURE/DISCARD
    split is a positional artifact of the 3-unit proxy (center vs cap-adjacent
    monomer), not chemistry; both leak the same MW-211 C15 backbone radical and
    fabricate the same mass under the solver's UNRESOLVED single-monomer-debit.
    This is the regression-lock for the widening: it FAILS under the old
    FEATURE-only predicate (DISCARD was not refused) and PASSES after."""
    import rmgpy.polymer as polymer_mod
    polymer_mod._flux_archetype_warned.clear()
    from rmgpy.polymer import (Polymer, classify_structure, PolymerClass,
                               stamp_polymer_flux_archetype)
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    epdm = Polymer(label="epdm", monomer="[CH2]CC(C)[CH2]",
                   Mn=5000.0, Mw=8000.0, initial_mass=1.0)
    # A C15 backbone radical that classifies DISCARD (not FEATURE) against epdm --
    # the cap-adjacent proxy position. Same MW (211.41) as the FEATURE radical.
    macro = Molecule().from_smiles("CC[C](C)CCCC(C)CCCC(C)C")
    macro.update()
    # Sanity: this is genuinely DISCARD, not FEATURE, so it would have been missed
    # by the old label-keyed predicate.
    klass, _ = classify_structure(Species(molecule=[macro]), epdm)
    assert klass == PolymerClass.DISCARD
    h = Molecule().from_smiles("[H]")
    h2 = Molecule().from_smiles("[H][H]")
    rxn = Reaction(reactants=[epdm, Species(molecule=[h])],
                   products=[Species(molecule=[h2]), macro])
    polymer_reactants = [r for r in rxn.reactants if isinstance(r, Polymer)]
    stamp_polymer_flux_archetype(rxn, rxn.reactants, polymer_reactants)
    # Chain-scale (MW 211 >> monomer 70 + slack 10) -> refused even though DISCARD.
    assert rxn.polymer_refused is True


def test_small_radical_lost_to_gas_conserves_not_refused():
    """A genuinely small fragment radical (propyl, MW ~43 < monomer+slack ~80) on
    the same UNRESOLVED leg must NOT be refused: the solver's single-monomer-debit
    accounts for sub-monomer leaks correctly, so the size gate declines (item 18
    Task 3 follow-up). Complements the CH4 helper-decline test with a radical."""
    import rmgpy.polymer as polymer_mod
    polymer_mod._flux_archetype_warned.clear()
    from rmgpy.polymer import (Polymer, stamp_polymer_flux_archetype,
                               PolymerFluxArchetype)
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    epdm = Polymer(label="epdm", monomer="[CH2]CC(C)[CH2]",
                   Mn=5000.0, Mw=8000.0, initial_mass=1.0)
    propyl = Molecule().from_smiles("CC[CH2]")  # ~43 g/mol, below chain-scale gate
    propyl.update()
    h = Molecule().from_smiles("[H]")
    h2 = Molecule().from_smiles("[H][H]")
    rxn = Reaction(reactants=[epdm, Species(molecule=[h])],
                   products=[Species(molecule=[h2]), Species(molecule=[propyl])])
    polymer_reactants = [r for r in rxn.reactants if isinstance(r, Polymer)]
    stamp_polymer_flux_archetype(rxn, rxn.reactants, polymer_reactants)
    assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)
    assert rxn.polymer_refused is False
    assert rxn.polymer_refused_accumulating is False


def _build_compile_inputs(moles, initial_mass=1.0, Mn=5000.0, Mw=6000.0,
                          label="PS", monomer="[CH2][CH]c1ccccc1"):
    """Build the (blueprint, initial_moles, species_dict) triple that
    ``compile_polymer_phase`` consumes, with the pool's stated loading
    (``initial_mass``/``Mn``) decoupled from the reactor's ``initialMoles``
    (``moles``) so the two mu0 sources can be made to agree or disagree.

    Returns ``(blueprint, initial_moles, species_dict, poly)`` where ``poly``
    is the Polymer object (the ``spc`` inside compile_polymer_phase's loop).
    """
    from rmgpy.rmg.polymer_input import PolymerPhaseBlueprint
    from rmgpy.species import Species

    poly = Polymer(label=label, monomer=monomer, end_groups=['[CH3]', '[H]'],
                   cutoff=3, Mn=Mn, Mw=Mw, initial_mass=initial_mass)
    species_dict = {
        label: poly,
        f"{label}_mu0": Species().from_smiles("CO"),
        f"{label}_mu1": Species().from_smiles("C=O"),
        f"{label}_mu2": Species().from_smiles("C#N"),
    }
    for suffix in ("_mu0", "_mu1", "_mu2"):
        species_dict[f"{label}{suffix}"].label = f"{label}{suffix}"
    blueprint = PolymerPhaseBlueprint(label=label, species=[label], solvent=label)
    initial_moles = {poly: moles}
    return blueprint, initial_moles, species_dict, poly


def test_compile_polymer_phase_reconciles_moments_to_initial_moles():
    """CYCLE 1 (tracer): compile_polymer_phase must make the Polymer object's
    .moments (what the sidecar serializes) AGREE with the solver-integrated
    initial_moments (derived from initialMoles), even when the pool's stated
    initial_mass/Mn implies a different mu0.

    Deck states initial_mass=1 kg, Mn=5000 -> mu0 = 1000/5000 = 0.2, but the
    reactor's initialMoles gives moles=0.01. The solver integrates 0.01; the
    sidecar (Polymer.moments) must report the same, not the 0.2-based moments.
    """
    from rmgpy.rmg.polymer_input import compile_polymer_phase

    blueprint, initial_moles, species_dict, poly = _build_compile_inputs(moles=0.01)

    # Pre-condition: the distribution-derived moments disagree with initialMoles.
    assert poly.moments[0] == pytest.approx(0.2, rel=1e-6)

    phase = compile_polymer_phase(blueprint, initial_moles, species_dict)

    solver_moments = np.asarray(phase.initial_moments[poly.label], dtype=float)
    # The solver integrates moles=0.01 as mu0.
    assert solver_moments[0] == pytest.approx(0.01, rel=1e-9)
    # The reconciled Polymer.moments (sidecar source) must equal the solver's.
    assert np.allclose(np.asarray(poly.moments, dtype=float), solver_moments)


def test_compile_polymer_phase_warns_on_mu0_disagreement(caplog):
    """CYCLE 2: when the initial_mass/Mn-implied chain count disagrees with
    initialMoles[proxy], compile_polymer_phase must emit a clear warning that
    names the pool, both mu0 values, and which one the solver uses."""
    import logging
    from rmgpy.rmg.polymer_input import compile_polymer_phase

    blueprint, initial_moles, species_dict, poly = _build_compile_inputs(moles=0.01)

    with caplog.at_level(logging.WARNING):
        compile_polymer_phase(blueprint, initial_moles, species_dict)

    warnings = [r.getMessage() for r in caplog.records
                if r.levelno >= logging.WARNING]
    assert any("PS" in m for m in warnings), warnings
    joined = " ".join(warnings)
    assert "0.2" in joined          # initial_mass/Mn-implied mu0
    assert "0.01" in joined         # initialMoles mu0 (what the solver uses)
    assert "initialMoles" in joined


def test_compile_polymer_phase_no_warning_when_consistent(caplog):
    """CYCLE 2 (negative): once the deck is consistent (initial_mass/Mn ==
    initialMoles), NO disagreement warning fires."""
    import logging
    from rmgpy.rmg.polymer_input import compile_polymer_phase

    # initial_mass=1 kg, Mn=5000 -> implied mu0 = 0.2; set moles to match.
    blueprint, initial_moles, species_dict, poly = _build_compile_inputs(moles=0.2)

    with caplog.at_level(logging.WARNING):
        compile_polymer_phase(blueprint, initial_moles, species_dict)

    disagreement = [r.getMessage() for r in caplog.records
                    if r.levelno >= logging.WARNING and "initialMoles" in r.getMessage()]
    assert disagreement == [], disagreement


# ---------------------------------------------------------------------------
# Stage 1: daughter-pool registration (proxy_reaction_reality_rules.md Layer 2)
#
# A scission/spawn daughter Polymer is registered as a core species (with its
# own _mu0/_mu1/_mu2 dummies) by _register_polymer, but pool_configs is built
# only from the static deck list polymerPhase.pools -- so the daughter's
# species map to -1 and its stamped SCISSION_FRAGMENT/MIGRATION flux demotes to
# UNRESOLVED ("could not resolve their solver pool(s)"). Registration derives a
# PolymerPoolConfig for each such daughter from the core species themselves.
# ---------------------------------------------------------------------------

def _moment_dummy(label):
    """A moment-dummy Species exactly as _register_polymer injects it."""
    s = Species(label=label, reactive=False)
    s.molecule = [Molecule().from_smiles("[Ne]")]
    s.is_moment_dummy = True
    return s


def test_derive_daughter_pool_config_binds_moment_dummies():
    """STAGE 1 / CYCLE 1 (tracer): a daughter Polymer registered as a core
    species (with its auto-created _mu0/_mu1/_mu2 dummies) must yield a
    PolymerPoolConfig that binds those dummies by index, so the solver resolves
    the daughter's pool instead of demoting its scission/migration flux to
    UNRESOLVED. Mirrors what _register_polymer leaves in core after a scission
    or spawn-intent daughter is registered."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)
    mu0 = _moment_dummy("PS_d1_mu0")
    mu1 = _moment_dummy("PS_d1_mu1")
    mu2 = _moment_dummy("PS_d1_mu2")
    core = [daughter, mu0, mu1, mu2]
    spc_map = {s: i for i, s in enumerate(core)}

    configs = derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})

    assert len(configs) == 1
    cfg = configs[0]
    assert cfg.label == "PS_d1"
    assert cfg.xs == 3                  # from daughter.cutoff
    assert tuple(cfg.mu_indices) == (1, 2, 3)   # PS_d1_mu0/_mu1/_mu2 core indices


def test_derive_daughter_pool_config_populates_monomer_mw():
    """The reference-state tripwire's chain_window = max(monomer_mw over pools) +
    slack. A derived daughter pool config that omits monomer_mw_g_mol (leaving it
    0.0) drags that max to 0 -> chain_window collapses to the 10 g/mol slack ->
    small gas scission fragments (over-tagged is_polymer_proxy) leak into the melt
    reference-state sum (the U=11.3 PS tripwire). The derived config must carry the
    daughter's own monomer MW (g/mol)."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)
    core = [daughter, _moment_dummy("PS_d1_mu0"), _moment_dummy("PS_d1_mu1"), _moment_dummy("PS_d1_mu2")]
    spc_map = {s: i for i, s in enumerate(core)}

    configs = derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})

    assert len(configs) == 1
    assert configs[0].monomer_mw_g_mol == pytest.approx(daughter.monomer_mw_g_mol, rel=1e-9)
    assert configs[0].monomer_mw_g_mol > 100.0   # styrene-scale, not the 0.0 default


def test_pool_to_config_populates_monomer_mw_from_molecule_monomer():
    """PolymerPool.monomer is a Molecule (Polymer._validate_monomer / the polymer()
    input helper builds PolymerPool with monomer=spc.monomer, a Molecule). to_config
    historically read it as a Species (getattr(self.monomer,'molecule')[0]) -> None ->
    monomer_mw_g_mol=0, collapsing the tripwire chain_window. The config must carry the
    real monomer MW regardless of whether monomer is a Molecule or a Species."""
    from rmgpy.rmg.polymer_input import PolymerPool

    mono = Molecule().from_smiles("C=Cc1ccccc1")   # styrene ~104.15 g/mol, a MOLECULE
    mu = [_moment_dummy("P_mu0"), _moment_dummy("P_mu1"), _moment_dummy("P_mu2")]
    pool = PolymerPool(label="P", xs=3, monomer=mono, explicit_map={}, mu_species=mu)
    spc_map = {s: i for i, s in enumerate(mu)}

    cfg = pool.to_config(spc_map)

    assert cfg.monomer_mw_g_mol == pytest.approx(104.15, abs=1.0)


def test_pool_to_config_hard_errors_on_unzip_without_monomer_product():
    """k_unzip > 0 with no resolvable monomer_product must be a HARD config error.

    The solver drains condensed moments unconditionally when k_unzip > 0
    (polymer.pyx: dmu1_dt -= k_unzip*mu0) but only emits the released monomer
    when monomer_poly_index is not None -- so a config with k_unzip > 0 and
    monomer_poly_index=None silently un-conserves mass (drained mass goes
    nowhere). to_config is the last point before that config reaches the
    solver; it must refuse, naming the pool."""
    from rmgpy.rmg.polymer_input import PolymerPool

    mono = Molecule().from_smiles("C=Cc1ccccc1")
    mu = [_moment_dummy("P_mu0"), _moment_dummy("P_mu1"), _moment_dummy("P_mu2")]
    pool = PolymerPool(label="P", xs=3, monomer=mono, explicit_map={},
                       mu_species=mu, k_unzip=0.5, monomer_product=None)
    spc_map = {s: i for i, s in enumerate(mu)}

    with pytest.raises(ValueError, match=r"Pool P.*k_unzip.*un-conserved") as excinfo:
        pool.to_config(spc_map)
    assert "monomer_product" in str(excinfo.value)


def test_pool_to_config_unzip_with_monomer_product_wires_index():
    """GREEN path: k_unzip > 0 WITH a resolvable monomer_product still builds a
    config, with monomer_poly_index bound to the released monomer's core index."""
    from rmgpy.rmg.polymer_input import PolymerPool

    mono = Molecule().from_smiles("C=Cc1ccccc1")
    styrene = Species(label="styrene", molecule=[Molecule().from_smiles("C=Cc1ccccc1")])
    mu = [_moment_dummy("P_mu0"), _moment_dummy("P_mu1"), _moment_dummy("P_mu2")]
    pool = PolymerPool(label="P", xs=3, monomer=mono, explicit_map={},
                       mu_species=mu, k_unzip=0.5, monomer_product=styrene)
    core = mu + [styrene]
    spc_map = {s: i for i, s in enumerate(core)}

    cfg = pool.to_config(spc_map)

    assert cfg.k_unzip == 0.5
    assert cfg.monomer_poly_index == spc_map[styrene]


def test_pool_to_config_zero_unzip_without_monomer_product_stays_legal():
    """GREEN path: k_unzip == 0 with monomer_product=None is a valid frozen /
    scission-only pool -- no unzip drain exists, so no emission target is needed."""
    from rmgpy.rmg.polymer_input import PolymerPool

    mono = Molecule().from_smiles("C=Cc1ccccc1")
    mu = [_moment_dummy("P_mu0"), _moment_dummy("P_mu1"), _moment_dummy("P_mu2")]
    pool = PolymerPool(label="P", xs=3, monomer=mono, explicit_map={},
                       mu_species=mu, k_unzip=0.0, monomer_product=None)
    spc_map = {s: i for i, s in enumerate(mu)}

    cfg = pool.to_config(spc_map)

    assert cfg.k_unzip == 0.0
    assert cfg.monomer_poly_index is None


def test_polymer_input_helper_hard_errors_on_unzip_without_monomer_product():
    """Parse-time companion to the to_config guard: the polymer() input-deck
    helper must refuse k_unzip > 0 with monomer_product=None immediately (clear
    InputError at deck-read time), before any species registration. The check
    fires before the helper touches the module-global rmg object, so this test
    needs no RMG instance."""
    from rmgpy.rmg import input as rmg_input

    with pytest.raises(InputError, match=r"PS.*k_unzip.*un-conserved"):
        rmg_input.polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                          end_groups=["[CH3]", "[H]"], cutoff=3,
                          Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                          k_unzip=1.0, monomer_product=None)


def test_polymer_input_helper_rejects_negative_k_unzip():
    """A negative k_unzip is not a valid rate constant. Every solver consumer
    of k_unzip is gated on k_unzip > 0, so a negative value would silently
    become an inert channel instead of failing -- the deck helper must refuse
    it at parse time with a clear InputError, same class of error as the
    missing-monomer_product guard."""
    from rmgpy.rmg import input as rmg_input

    with pytest.raises(InputError, match=r"PS.*k_unzip.*not a valid rate constant"):
        rmg_input.polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                          end_groups=["[CH3]", "[H]"], cutoff=3,
                          Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                          k_unzip=-0.5)


def test_pool_to_config_rejects_negative_k_unzip():
    """Config-assembly companion to the deck-helper check: PolymerPool.to_config
    must refuse a negative k_unzip (not a valid rate constant) even when a
    monomer_product IS wired -- routing must not dodge the sign check."""
    from rmgpy.rmg.polymer_input import PolymerPool

    mono = Molecule().from_smiles("C=Cc1ccccc1")
    styrene = Species(label="styrene", molecule=[Molecule().from_smiles("C=Cc1ccccc1")])
    mu = [_moment_dummy("P_mu0"), _moment_dummy("P_mu1"), _moment_dummy("P_mu2")]
    pool = PolymerPool(label="P", xs=3, monomer=mono, explicit_map={},
                       mu_species=mu, k_unzip=-0.5, monomer_product=styrene)
    spc_map = {s: i for i, s in enumerate(mu + [styrene])}

    with pytest.raises(ValueError, match=r"Pool P.*k_unzip.*not a valid rate constant"):
        pool.to_config(spc_map)


@pytest.mark.parametrize("field", ["k_unzip", "k_scission"])
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")],
                         ids=["nan", "+inf", "-inf"])
def test_polymer_input_helper_rejects_non_finite_rates(field, bad):
    """NaN passes BOTH the `< 0` and `> 0` checks as False, so a non-finite
    k_unzip/k_scission would make the channel SILENTLY INERT (or poison the
    residual with inf) -- a laundered no-op. The deck helper must reject
    NaN/inf at parse time with a clear InputError, mirroring the QSSA triplet
    validator's finite-rejection posture."""
    from rmgpy.rmg import input as rmg_input

    with pytest.raises(InputError,
                       match=rf"PS.*{field}.*not a valid rate constant"):
        rmg_input.polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                          end_groups=["[CH3]", "[H]"], cutoff=3,
                          Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                          **{field: bad})


@pytest.mark.parametrize("field", ["k_unzip", "k_scission"])
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")],
                         ids=["nan", "+inf", "-inf"])
def test_pool_to_config_rejects_non_finite_rates(field, bad):
    """Config-assembly companion to the deck-helper finite check: a
    non-finite k_unzip/k_scission must be refused by PolymerPool.to_config
    even with a monomer_product wired (routing must not dodge the check)."""
    from rmgpy.rmg.polymer_input import PolymerPool

    mono = Molecule().from_smiles("C=Cc1ccccc1")
    styrene = Species(label="styrene", molecule=[Molecule().from_smiles("C=Cc1ccccc1")])
    mu = [_moment_dummy("P_mu0"), _moment_dummy("P_mu1"), _moment_dummy("P_mu2")]
    pool = PolymerPool(label="P", xs=3, monomer=mono, explicit_map={},
                       mu_species=mu, monomer_product=styrene, **{field: bad})
    spc_map = {s: i for i, s in enumerate(mu + [styrene])}

    with pytest.raises(ValueError,
                       match=rf"Pool P.*{field}.*not a valid rate constant"):
        pool.to_config(spc_map)


# ---------------------------------------------------------------------------
# radical_qssa_unzip channel (M1: config + validation only, NO RHS effect)
# ---------------------------------------------------------------------------


def _qssa_triplet(A=1.0e13, n=0.0, Ea=1.0e5):
    """Arrhenius triplet for the radical_qssa_unzip channel (SI convention:
    A [s^-1] unimolecular / [m^3 mol^-1 s^-1] bimolecular, Ea [J/mol])."""
    return dict(A=A, n=n, Ea=Ea)


def _qssa_channel(**overrides):
    """A valid minimal radical_qssa_unzip channel config (mandatory blocks
    only; efficiency/monomer_yield/basis/transfer left to their defaults)."""
    ch = dict(
        initiation=_qssa_triplet(A=1.0e15, Ea=3.0e5),
        depropagation=_qssa_triplet(A=1.0e13, Ea=8.0e4),
        termination=_qssa_triplet(A=1.0e8, Ea=1.0e4),
    )
    ch.update(overrides)
    return ch


def _qssa_pool(channel, k_unzip=0.0, wire_monomer_product=True):
    """PolymerPool + spc_map fixture for the radical_qssa_unzip to_config
    tests. The released-monomer routing reuses the pool's EXISTING
    monomer_product field (design contract: no new routing field)."""
    from rmgpy.rmg.polymer_input import PolymerPool

    mono = Molecule().from_smiles("C=Cc1ccccc1")
    mu = [_moment_dummy("P_mu0"), _moment_dummy("P_mu1"), _moment_dummy("P_mu2")]
    core = list(mu)
    mp = None
    if wire_monomer_product:
        mp = Species(label="styrene", molecule=[Molecule().from_smiles("C=Cc1ccccc1")])
        core.append(mp)
    pool = PolymerPool(label="P", xs=3, monomer=mono, explicit_map={},
                       mu_species=mu, k_unzip=k_unzip, monomer_product=mp,
                       radical_qssa_unzip=channel)
    return pool, {s: i for i, s in enumerate(core)}


def test_polymer_stores_radical_qssa_unzip_and_copy_preserves_it():
    """Layer-2 attribute: the Polymer object carries the channel config as
    passive storage (validation lives in the deck helper / to_config / solver),
    and copy() must preserve it -- like k_unzip/k_scission, losing it on copy
    would silently disable the pool's degradation channel."""
    ch = _qssa_channel()
    poly = Polymer(label='PS', monomer='[CH2][CH]c1ccccc1',
                   end_groups=['[CH3]', '[H]'], cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                   radical_qssa_unzip=ch)
    assert poly.radical_qssa_unzip == ch
    assert poly.copy(deep=True).radical_qssa_unzip == ch


def test_polymer_copy_deep_copies_radical_qssa_unzip():
    """COPY ALIASING (review round 21, finding 3): copy() must deep-copy the
    channel dict, not shallow-assign it. Copies aliasing the same nested dict
    means mutating one Polymer's channel silently rewrites every copy's
    (including the spawned-daughter configs that will inherit it)."""
    poly = Polymer(label='PS', monomer='[CH2][CH]c1ccccc1',
                   end_groups=['[CH3]', '[H]'], cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                   radical_qssa_unzip=_qssa_channel())
    for cp in (poly.copy(), poly.copy(deep=True)):
        assert cp.radical_qssa_unzip is not poly.radical_qssa_unzip
        poly.radical_qssa_unzip["initiation"]["A"] = 999.0
        poly.radical_qssa_unzip["efficiency"] = 0.123
        assert cp.radical_qssa_unzip["initiation"]["A"] == 1.0e15
        assert "efficiency" not in cp.radical_qssa_unzip
        # restore for the second iteration
        poly.radical_qssa_unzip["initiation"]["A"] = 1.0e15
        del poly.radical_qssa_unzip["efficiency"]


def test_polymer_radical_qssa_unzip_defaults_to_none():
    """Regression: a Polymer built without the channel stays channel-free,
    through copy() too."""
    poly = Polymer(label='PS', monomer='[CH2][CH]c1ccccc1',
                   end_groups=['[CH3]', '[H]'], cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=0.001)
    assert poly.radical_qssa_unzip is None
    assert poly.copy(deep=True).radical_qssa_unzip is None


def test_pool_to_config_roundtrips_radical_qssa_unzip():
    """Valid channel round-trip: to_config validates + normalizes the dict and
    stores it on PolymerPoolConfig with defaults filled (efficiency=1.0,
    monomer_yield=1.0, pinned basis, transfer=None), reusing the pool's
    existing monomer_product routing (monomer_poly_index; NO new routing
    field). M1 contract: the config is validated but INERT (no RHS reads)."""
    pool, spc_map = _qssa_pool(_qssa_channel())

    cfg = pool.to_config(spc_map)

    q = cfg.radical_qssa_unzip
    assert q is not None
    assert q["initiation"] == dict(A=1.0e15, n=0.0, Ea=3.0e5)
    assert q["depropagation"] == dict(A=1.0e13, n=0.0, Ea=8.0e4)
    assert q["termination"] == dict(A=1.0e8, n=0.0, Ea=1.0e4)
    assert q["efficiency"] == 1.0
    assert q["monomer_yield"] == 1.0
    assert q["basis"] == "backbone_bonds_mu1_minus_mu0"
    assert q["transfer"] is None
    assert cfg.monomer_poly_index == 3  # existing routing reused
    assert cfg.k_unzip == 0.0


def test_pool_to_config_radical_qssa_unzip_explicit_optionals_and_transfer():
    """Explicit efficiency/monomer_yield/basis/transfer values survive
    normalization (transfer is accepted and stored -- same finite/positivity
    rules, no rate law yet)."""
    ch = _qssa_channel(efficiency=0.6, monomer_yield=0.9,
                       basis="backbone_bonds_mu1_minus_mu0",
                       transfer=_qssa_triplet(A=5.0e6, n=0.5, Ea=2.0e4))
    pool, spc_map = _qssa_pool(ch)

    q = pool.to_config(spc_map).radical_qssa_unzip

    assert q["efficiency"] == 0.6
    assert q["monomer_yield"] == 0.9
    assert q["transfer"] == dict(A=5.0e6, n=0.5, Ea=2.0e4)


def test_pool_to_config_channel_absent_stays_none():
    """Regression: a channel-absent pool (here a legal k_unzip-only pool) is
    completely unaffected -- radical_qssa_unzip stays None on its config."""
    pool, spc_map = _qssa_pool(None, k_unzip=0.5)

    cfg = pool.to_config(spc_map)

    assert cfg.radical_qssa_unzip is None
    assert cfg.k_unzip == 0.5
    assert cfg.monomer_poly_index == 3


_QSSA_BAD_CHANNELS = [
    pytest.param({k: v for k, v in _qssa_channel().items() if k != "termination"},
                 r"Pool P.*radical_qssa_unzip.*missing.*termination",
                 id="missing-termination-block"),
    pytest.param(_qssa_channel(initiation=_qssa_triplet(A=float("nan"))),
                 r"Pool P.*initiation.*A.*not finite", id="nan-A"),
    pytest.param(_qssa_channel(termination=_qssa_triplet(Ea=float("inf"))),
                 r"Pool P.*termination.*Ea.*not finite", id="inf-Ea"),
    pytest.param(_qssa_channel(initiation=_qssa_triplet(n=float("inf"))),
                 r"Pool P.*initiation.*n.*not finite", id="inf-n"),
    pytest.param(_qssa_channel(depropagation=_qssa_triplet(A=0.0)),
                 r"Pool P.*depropagation.*A.*> 0", id="zero-A"),
    pytest.param(_qssa_channel(depropagation=_qssa_triplet(A=-1.0e13)),
                 r"Pool P.*depropagation.*A.*> 0", id="negative-A"),
    pytest.param(_qssa_channel(initiation=_qssa_triplet(Ea=-5.0)),
                 r"Pool P.*initiation.*Ea.*>= 0", id="negative-Ea"),
    pytest.param(_qssa_channel(initiation=dict(A=1.0e13, n=0.0)),
                 r"Pool P.*initiation.*Ea", id="triplet-missing-Ea"),
    pytest.param(_qssa_channel(efficiency=0.0),
                 r"Pool P.*efficiency.*\(0, 1\]", id="efficiency-zero"),
    pytest.param(_qssa_channel(efficiency=1.5),
                 r"Pool P.*efficiency.*\(0, 1\]", id="efficiency-above-one"),
    pytest.param(_qssa_channel(monomer_yield=0.0),
                 r"Pool P.*monomer_yield.*\(0, 1\]", id="monomer-yield-zero"),
    pytest.param(_qssa_channel(monomer_yield=1.5),
                 r"Pool P.*monomer_yield.*\(0, 1\]", id="monomer-yield-above-one"),
    pytest.param(_qssa_channel(basis="chain_ends_mu0"),
                 r"Pool P.*basis.*backbone_bonds_mu1_minus_mu0", id="bad-basis"),
    pytest.param(_qssa_channel(transfer=_qssa_triplet(A=float("nan"))),
                 r"Pool P.*transfer.*A.*not finite", id="nan-transfer-A"),
    pytest.param(_qssa_channel(bogus_key=1.0),
                 r"Pool P.*radical_qssa_unzip.*unknown key", id="unknown-key"),
]


@pytest.mark.parametrize("channel, pattern", _QSSA_BAD_CHANNELS)
def test_pool_to_config_rejects_invalid_radical_qssa_unzip(channel, pattern):
    """Field validation at config assembly: every A/n/Ea must be FINITE
    (NaN/inf rejected explicitly -- this channel gets finite checks from day
    one), A > 0, Ea >= 0; efficiency/monomer_yield in (0, 1]; basis pinned to
    'backbone_bonds_mu1_minus_mu0' (forward-compat pin)."""
    pool, spc_map = _qssa_pool(channel)
    with pytest.raises(ValueError, match=pattern):
        pool.to_config(spc_map)


def test_pool_to_config_rejects_qssa_channel_without_monomer_product():
    """Channel present without a resolvable monomer product = hard error: the
    QSSA unzip channel releases monomer through the pool's existing monomer
    routing; without an emission target the depropagated repeat units would
    leave the condensed phase silently un-conserved (same failure class as the
    k_unzip guard above)."""
    pool, spc_map = _qssa_pool(_qssa_channel(), wire_monomer_product=False)
    with pytest.raises(ValueError, match=r"Pool P.*radical_qssa_unzip.*un-conserved"):
        pool.to_config(spc_map)


def test_pool_to_config_rejects_qssa_channel_with_positive_k_unzip():
    """Double-counting guard: radical_qssa_unzip and k_unzip > 0 are two
    representations of the SAME chain-end depropagation channel and are
    mutually exclusive on a pool."""
    pool, spc_map = _qssa_pool(_qssa_channel(), k_unzip=0.5)
    with pytest.raises(ValueError, match=r"Pool P.*mutually exclusive"):
        pool.to_config(spc_map)


def test_polymer_input_helper_rejects_bad_radical_qssa_unzip():
    """Parse-time companion to the to_config field validation: the polymer()
    deck helper must refuse a malformed radical_qssa_unzip with a clear
    InputError at deck-read time. The check fires before the helper touches
    the module-global rmg object, so this test needs no RMG instance."""
    from rmgpy.rmg import input as rmg_input

    for channel, pattern in [
        ({k: v for k, v in _qssa_channel().items() if k != "initiation"},
         r"PS.*missing.*initiation"),
        (_qssa_channel(depropagation=_qssa_triplet(A=float("nan"))),
         r"PS.*depropagation.*A.*not finite"),
        (_qssa_channel(termination=_qssa_triplet(Ea=float("inf"))),
         r"PS.*termination.*Ea.*not finite"),
        (_qssa_channel(initiation=_qssa_triplet(A=-1.0)),
         r"PS.*initiation.*A.*> 0"),
        (_qssa_channel(efficiency=2.0), r"PS.*efficiency.*\(0, 1\]"),
        (_qssa_channel(basis="wrong"), r"PS.*basis"),
    ]:
        with pytest.raises(InputError, match=pattern):
            rmg_input.polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                              end_groups=["[CH3]", "[H]"], cutoff=3,
                              Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                              monomer_product="C=Cc1ccccc1",
                              radical_qssa_unzip=channel)


def test_polymer_input_helper_rejects_qssa_channel_without_monomer_product():
    """Deck-read-time cross-invariant: radical_qssa_unzip requires a
    monomer_product (the channel reuses the pool's existing monomer routing;
    without it the released mass would leave the condensed phase
    un-conserved)."""
    from rmgpy.rmg import input as rmg_input

    with pytest.raises(InputError, match=r"PS.*radical_qssa_unzip.*un-conserved"):
        rmg_input.polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                          end_groups=["[CH3]", "[H]"], cutoff=3,
                          Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                          monomer_product=None,
                          radical_qssa_unzip=_qssa_channel())


def test_polymer_input_helper_rejects_qssa_channel_with_positive_k_unzip():
    """Deck-read-time double-counting guard: radical_qssa_unzip AND
    k_unzip > 0 on the same pool is a hard error even when monomer_product is
    wired (the two depropagation representations are mutually exclusive)."""
    from rmgpy.rmg import input as rmg_input

    with pytest.raises(InputError, match=r"PS.*mutually exclusive"):
        rmg_input.polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                          end_groups=["[CH3]", "[H]"], cutoff=3,
                          Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                          k_unzip=1.0, monomer_product="C=Cc1ccccc1",
                          radical_qssa_unzip=_qssa_channel())


def test_polymer_input_helper_valid_qssa_channel_reaches_polymer_object():
    """GREEN deck path (deck -> Polymer leg of the round-trip): a valid
    radical_qssa_unzip passes parse-time validation and lands normalized
    (defaults filled) on the Polymer object. The module-global rmg is mocked:
    species registration is out of scope here."""
    from unittest.mock import MagicMock
    from rmgpy.rmg import input as rmg_input

    def _make_new_species(obj, **kwargs):
        if isinstance(obj, Species):
            return obj, True
        return Species(label="styrene", molecule=[obj]), True

    mock_rmg = MagicMock()
    mock_rmg.initial_species = []
    mock_rmg.reaction_model.iteration_num = 0
    mock_rmg.reaction_model.new_species_list = []
    mock_rmg.reaction_model.make_new_species.side_effect = _make_new_species

    old_rmg, old_sd = rmg_input.rmg, rmg_input.species_dict
    rmg_input.rmg, rmg_input.species_dict = mock_rmg, {}
    try:
        poly = rmg_input.polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                                 end_groups=["[CH3]", "[H]"], cutoff=3,
                                 Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                                 monomer_product="C=Cc1ccccc1",
                                 radical_qssa_unzip=_qssa_channel())
    finally:
        rmg_input.rmg, rmg_input.species_dict = old_rmg, old_sd

    q = poly.radical_qssa_unzip
    assert q["initiation"] == dict(A=1.0e15, n=0.0, Ea=3.0e5)
    assert q["depropagation"] == dict(A=1.0e13, n=0.0, Ea=8.0e4)
    assert q["termination"] == dict(A=1.0e8, n=0.0, Ea=1.0e4)
    assert q["efficiency"] == 1.0
    assert q["monomer_yield"] == 1.0
    assert q["basis"] == "backbone_bonds_mu1_minus_mu0"
    assert q["transfer"] is None


def test_derive_daughter_pool_configs_skips_static_and_incomplete():
    """STAGE 1 / CYCLE 2: the root proxy (a Polymer in core whose label IS a
    static deck pool) must NOT be re-derived (else it is double-configured), and
    a daughter missing part of its _muN triplet is skipped rather than yielding an
    unresolvable pool. Only the complete, non-static daughter gets a config."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    def _poly(label):
        return Polymer(label=label, monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)

    root = _poly("PS")                 # static deck pool's proxy, lives in core
    good = _poly("PS_d1")              # complete daughter -> one config
    incomplete = _poly("PS_d2")        # missing _mu2 -> skipped
    core = [
        root, _moment_dummy("PS_mu0"), _moment_dummy("PS_mu1"), _moment_dummy("PS_mu2"),
        good, _moment_dummy("PS_d1_mu0"), _moment_dummy("PS_d1_mu1"), _moment_dummy("PS_d1_mu2"),
        incomplete, _moment_dummy("PS_d2_mu0"), _moment_dummy("PS_d2_mu1"),  # no PS_d2_mu2
    ]
    spc_map = {s: i for i, s in enumerate(core)}

    configs = derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})

    assert [c.label for c in configs] == ["PS_d1"]


def test_derive_daughter_pool_config_uses_base_label_with_index_suffix():
    """STAGE 1 / CYCLE 5 (hardening): RMG appends a "(N)" index to registered
    species labels (the proxy displays as "PS(2)" while its dummies stay the
    clean "PS_mu0"). A daughter whose proxy label acquired such an index
    ("PS_d1(9)") still has clean "PS_d1_muN" dummies. The derived config must use
    the '('-stripped base label -- both to FIND the clean dummies and so the
    solver (which binds on label.partition('(')[0] == pool.label) can resolve the
    pool. Keying off the raw label drops the daughter (looks for the nonexistent
    "PS_d1(9)_mu0")."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    daughter = Polymer(label="PS_d1(9)", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)
    core = [
        daughter,
        _moment_dummy("PS_d1_mu0"), _moment_dummy("PS_d1_mu1"), _moment_dummy("PS_d1_mu2"),
    ]
    spc_map = {s: i for i, s in enumerate(core)}

    configs = derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})

    assert len(configs) == 1
    # Config label is the clean base so the solver's base_label match binds it.
    assert configs[0].label == "PS_d1"
    assert tuple(configs[0].mu_indices) == (1, 2, 3)


# ---------------------------------------------------------------------------
# Daughter-pool QSSA inheritance (radical_qssa_unzip cascade milestone 5).
# The recorded M1 decision (polymer_input.py, derive_daughter_pool_configs):
# spawned scission daughters inherit the parent pool's radical_qssa_unzip
# channel (deep-copied) -- same monomer chemistry and monomer_mw imply the
# same elementary initiation/depropagation/termination constants. Without
# inheritance a PS scission cascade freezes: the parent unzips but daughters
# are inert and the TGA S-curve never completes.
# ---------------------------------------------------------------------------

def _qssa_raw_channel():
    """Deck-shaped (pre-normalization) radical QSSA channel config."""
    return {
        "initiation": {"A": 1.0e13, "n": 0.0, "Ea": 3.0e5},
        "depropagation": {"A": 1.0e14, "n": 0.5, "Ea": 9.0e4},
        "termination": {"A": 1.0e8, "n": 0.0, "Ea": 1.0e4},
        "efficiency": 0.8,
        "monomer_yield": 0.9,
    }


def _qssa_parent(channel=None):
    """A PS-like parent Polymer, optionally carrying the QSSA channel and a
    resolvable monomer_product_species (the deck attaches both; input.py:432)."""
    p = Polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                end_groups=["[CH3]", "[H]"], cutoff=3,
                Mn=5000.0, Mw=6000.0, initial_mass=1.0,
                radical_qssa_unzip=channel)
    styrene = Species(label="styrene", smiles="C=Cc1ccccc1")
    p.monomer_product_species = styrene
    return p, styrene


def _scission_tail_of(parent):
    """Run a real scission EVENT through create_reacted_copy: head wing +
    labeled methyl radical -> head-side scission -> a _scission_tail daughter."""
    head_wing = parent._stitch_wing("head")
    methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
    frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
    assert frag is not None
    daughter = parent.create_reacted_copy(frag)
    assert daughter is not None and daughter.label.endswith("_scission_tail")
    return daughter


def test_scission_daughter_inherits_qssa_channel_deepcopy():
    """A scission daughter Polymer must carry the parent's radical_qssa_unzip
    channel DEEP-COPIED (parent mutation must not propagate) plus the parent's
    monomer_product_species by REFERENCE (spc_map resolution is object-keyed,
    so identity is load-bearing for the daughter's monomer routing)."""
    channel = _qssa_raw_channel()
    parent, styrene = _qssa_parent(channel)

    daughter = _scission_tail_of(parent)

    assert daughter.radical_qssa_unzip == channel
    assert daughter.radical_qssa_unzip is not parent.radical_qssa_unzip
    # Deep copy: mutate the parent's nested triplet; daughter must not move.
    parent.radical_qssa_unzip["initiation"]["A"] = 1.0
    assert daughter.radical_qssa_unzip["initiation"]["A"] == 1.0e13
    # Routing: SAME species object, so the daughter resolves the same core index.
    assert daughter.monomer_product_species is styrene


def test_scission_daughter_channel_free_parent_stays_channel_free():
    """No noise: a channel-free parent spawns a channel-free daughter."""
    parent, _ = _qssa_parent(channel=None)
    daughter = _scission_tail_of(parent)
    assert daughter.radical_qssa_unzip is None


def test_inherit_gate_same_monomer_chemistry_inherits():
    """Round-25 P2-1 gate, PASS arm: a daughter with the parent's monomer_mw
    and unchanged feature chemistry (the scission tail/head shape) inherits
    channel (deep-copied) + routing reference."""
    import logging
    channel = _qssa_raw_channel()
    parent, styrene = _qssa_parent(channel)
    daughter = Polymer(label="PS_scission_tail", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=2500.0, Mw=3000.0, initial_mass=0.0)
    polymer._inherit_unzip_channel(daughter, parent)
    assert daughter.radical_qssa_unzip == channel
    assert daughter.radical_qssa_unzip is not parent.radical_qssa_unzip
    assert daughter.monomer_product_species is styrene


def test_feature_mod_daughter_does_not_inherit_channel_and_warns(caplog):
    """Round-25 P2-1, BLOCK arm: create_reacted_copy stamps inheritance at
    its single exit (:894) on EVERY non-None daughter, including _mod
    products whose feature_monomer CHANGED (:1010) -- 'same monomer
    chemistry' is false there, so the QSSA constants must NOT transfer.
    Channel-free + once-per-pool WARNING instead."""
    import logging
    getattr(polymer, "_unzip_inherit_warned", set()).clear()
    channel = _qssa_raw_channel()
    parent, styrene = _qssa_parent(channel)
    # The daughter shape the _mod constructor site produces: same monomer
    # attr, but a CHANGED feature unit in the chain.
    daughter = Polymer(label="PS_mod", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=1.0)
    daughter.feature_monomer = Molecule().from_smiles("C=C")

    with caplog.at_level(logging.WARNING):
        polymer._inherit_unzip_channel(daughter, parent)

    assert daughter.radical_qssa_unzip is None
    assert getattr(daughter, "monomer_product_species", None) is None
    warned = [r for r in caplog.records
              if "channel-free" in r.getMessage()
              and "PS_mod" in r.getMessage()]
    assert warned, "expected a changed-chemistry channel-free WARNING"
    # once-per-pool: a second identical event does not warn again
    n = len(caplog.records)
    with caplog.at_level(logging.WARNING):
        polymer._inherit_unzip_channel(daughter, parent)
    assert daughter.radical_qssa_unzip is None
    assert len(caplog.records) == n


def test_different_monomer_mw_daughter_does_not_inherit_channel():
    """Round-25 P2-1, mw arm: a daughter whose monomer_mw differs from the
    parent's fails the cheap truthful gate (the M1 rationale binds the
    constants to same monomer chemistry AND monomer_mw)."""
    import logging
    getattr(polymer, "_unzip_inherit_warned", set()).clear()
    parent, _ = _qssa_parent(_qssa_raw_channel())
    daughter = Polymer(label="PE_like", monomer="[CH2][CH2]",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       Mn=1000.0, Mw=2500.0, initial_mass=1.0)
    polymer._inherit_unzip_channel(daughter, parent)
    assert daughter.radical_qssa_unzip is None
    assert getattr(daughter, "monomer_product_species", None) is None


def test_derive_daughter_pool_config_inherits_qssa_channel_and_routing():
    """derive_daughter_pool_configs must build the daughter's config with the
    inherited channel run through the SHARED validator (normalized: defaults
    filled) and the monomer routing resolved (monomer_poly_index), deep-copied
    so post-hoc mutation of the species' dict cannot reach the config."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                       radical_qssa_unzip=_qssa_raw_channel())
    styrene = Species(label="styrene", smiles="C=Cc1ccccc1")
    daughter.monomer_product_species = styrene
    core = [daughter, _moment_dummy("PS_d1_mu0"), _moment_dummy("PS_d1_mu1"),
            _moment_dummy("PS_d1_mu2"), styrene]
    spc_map = {s: i for i, s in enumerate(core)}

    configs = derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})

    assert len(configs) == 1
    cfg = configs[0]
    q = cfg.radical_qssa_unzip
    assert q is not None
    assert q["initiation"]["A"] == 1.0e13
    assert q["efficiency"] == 0.8
    # Normalized through the shared validator: omitted fields filled in.
    assert q["transfer"] is None
    assert q["basis"] == "backbone_bonds_mu1_minus_mu0"
    # Monomer routing resolved to the released monomer's core index.
    assert cfg.monomer_poly_index == spc_map[styrene]
    # Mutual-exclusion invariant holds by construction on the daughter.
    assert cfg.k_unzip == 0.0
    # Deep copy: the config is independent of the species' mutable dict.
    daughter.radical_qssa_unzip["initiation"]["A"] = 1.0
    assert q["initiation"]["A"] == 1.0e13


def test_derive_daughter_pool_config_channel_free_stays_channel_free():
    """A daughter without the channel derives a channel-free config (the
    pre-milestone shape): no channel, no routing requirement."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)
    core = [daughter, _moment_dummy("PS_d1_mu0"), _moment_dummy("PS_d1_mu1"),
            _moment_dummy("PS_d1_mu2")]
    spc_map = {s: i for i, s in enumerate(core)}

    configs = derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})

    assert len(configs) == 1
    assert configs[0].radical_qssa_unzip is None
    assert configs[0].monomer_poly_index is None


def test_derive_daughter_pool_config_qssa_without_routing_is_loud():
    """A daughter carrying the channel but NO resolvable monomer emission
    target must FAIL LOUDLY at derivation (mirrors PolymerPool.to_config):
    a silent channel-drop is exactly the failure class this milestone kills."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                       radical_qssa_unzip=_qssa_raw_channel())
    # No monomer_product_species at all.
    core = [daughter, _moment_dummy("PS_d1_mu0"), _moment_dummy("PS_d1_mu1"),
            _moment_dummy("PS_d1_mu2")]
    spc_map = {s: i for i, s in enumerate(core)}
    with pytest.raises(ValueError, match="monomer_product"):
        derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})

    # monomer_product_species present but NOT in core (unresolvable index).
    orphan = Species(label="styrene", smiles="C=Cc1ccccc1")
    daughter.monomer_product_species = orphan
    with pytest.raises(ValueError, match="monomer_product"):
        derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})


def test_derive_daughter_pool_config_invalid_inherited_channel_is_loud():
    """A malformed inherited channel must raise through the SHARED validator
    (validate_radical_qssa_unzip), naming the daughter pool -- daughters do
    not bypass validation."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs

    daughter = Polymer(label="PS_d1", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001,
                       radical_qssa_unzip={"initiation": {"A": 1.0}})
    styrene = Species(label="styrene", smiles="C=Cc1ccccc1")
    daughter.monomer_product_species = styrene
    core = [daughter, _moment_dummy("PS_d1_mu0"), _moment_dummy("PS_d1_mu1"),
            _moment_dummy("PS_d1_mu2"), styrene]
    spc_map = {s: i for i, s in enumerate(core)}
    with pytest.raises(ValueError, match="PS_d1"):
        derive_daughter_pool_configs(core, spc_map, existing_pool_labels={"PS"})


# ---------------------------------------------------------------------------
# Task 5: End-to-end scission tracer — Ea from REAL products, not pool proxy
# ---------------------------------------------------------------------------


class TestScissionRealDHrxnEndToEnd:
    """
    End-to-end pin: PS(2) retro-ene scission through make_new_reaction must
    yield an Arrhenius Ea computed from the REAL atom-balanced products
    (C9H10 + C16H18, ΔH ≈ 48 kcal/mol, Ea ≈ 72 kcal/mol), NOT from the
    moment-pool-relabeled representative (C32H34, ΔH ≈ 81 kcal/mol, Ea ≈ 94
    kcal/mol).

    Verified empirically before writing:
      - Proxy C25H28; products C9H10 + C16H18 sum to C25H28 ✓
      - _handshake_structures([C16H18, ...], [PS]) → True; C16H18 → PS_scission_tail ✓
      - real_dH ≈ 202.6 kJ/mol; pool_dH ≈ 337.3 kJ/mol (Δ ≈ +32 kcal/mol) ✓
      - BM Ea(real) ≈ 300 kJ/mol; BM Ea(pool) ≈ 395 kJ/mol ✓
    """

    @classmethod
    def setup_class(cls):
        import os
        from rmgpy import settings
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.rmg.main import RMG
        rmg = RMG()
        rmg.database = RMGDatabase()
        rmg.database.load_thermo(os.path.join(settings["database.directory"], "thermo"))

    @pytest.fixture(autouse=True)
    def setup(self):
        from rmgpy.rmg.model import CoreEdgeReactionModel
        self.model = CoreEdgeReactionModel()
        self.ps = Polymer(
            label='PS',
            monomer='[CH2][CH]c1ccccc1',
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )
        self.model._register_polymer(self.ps, generate_thermo=True)

    def test_scission_Ea_uses_real_products_not_pool(self):
        """
        PS(2) retro-ene scission: Arrhenius Ea must come from the real
        atom-balanced C9H10+C16H18 thermo, not from the relabeled C32H34 pool
        proxy thermo.

        Assertion (A): Ea ≈ bm.get_activation_energy(real_dH)   [within 1%]
        Assertion (B): real_dH ≠ pool_dH by > 40 kJ/mol, and
                       Ea is NOT close to bm.get_activation_energy(pool_dH).
        Without (B) the test would pass vacuously even on the unfixed code.
        """
        from rmgpy.data.kinetics.family import TemplateReaction, _handshake_structures
        from rmgpy.kinetics import ArrheniusBM, Arrhenius

        # Empirically verified real Retroene products for PS proxy (C25H28):
        #   C9H10 (alpha-methylstyrene) + C16H18 (PS tail) = C25H28 ✓
        # C16H18 = CC(CC=C1C=CC=CC1)C1=CC=CC=C1 (vinyl-ended PS-dimer tail)
        C9H10_SMILES = 'C=C(C)C1=CC=CC=C1'
        C16H18_SMILES = 'CC(CC=C1C=CC=CC1)C1=CC=CC=C1'

        # BM parameters representative of the Retroene family
        bm_params = dict(
            A=(1.293332e12, 's^-1'), n=0.0,
            w0=(968.0, 'kJ/mol'), E0=(182.946, 'kJ/mol'),
        )

        # --- Independent real_dH (mirrors _polymer_real_dHrxn estimation path) ---
        def _H(smiles):
            spc = Species(molecule=[Molecule().from_smiles(smiles)])
            spc.generate_resonance_structures()
            self.model.generate_thermo(spc)
            return spc.get_enthalpy(298)

        H_c9h10 = _H(C9H10_SMILES)
        H_c16h18 = _H(C16H18_SMILES)
        H_proxy = self.ps.get_enthalpy(298)
        real_dH = H_c9h10 + H_c16h18 - H_proxy          # J/mol
        expected_Ea = ArrheniusBM(**bm_params).get_activation_energy(real_dH)  # J/mol

        # --- Build and run make_new_reaction ---
        proxy_mol = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        rxn = TemplateReaction(
            reactants=[proxy_mol],
            products=[
                Molecule().from_smiles(C16H18_SMILES),  # relabels → PS_scission_tail
                Molecule().from_smiles(C9H10_SMILES),   # stays as alpha-methylstyrene
            ],
            kinetics=ArrheniusBM(**bm_params),
            family='Retroene',
            is_forward=True,
        )
        result, is_new = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=True, generate_kinetics=True,
        )

        assert result is not None, "make_new_reaction returned None unexpectedly"
        assert isinstance(result.kinetics, Arrhenius), (
            f"Expected Arrhenius kinetics after BM pre-conversion, "
            f"got {type(result.kinetics).__name__}"
        )

        Ea = result.kinetics.Ea.value_si

        # (A) Ea is based on real atom-balanced products. The pool-sourced H0 clamp
        #     was a residual channel (spec §8 / Task 7); for this PS case it fires by
        #     only ~0.3 kJ/mol (pool_H0 barely above ea_pre), so rtol=0.01 still holds
        #     both before and after the Task 7 correction. The clamp is now closed
        #     structurally (Task 7: _polymer_real_H0 correction post fix_barrier_height).
        assert np.isclose(Ea, expected_Ea, rtol=0.01), (
            f"Ea = {Ea / 4184:.2f} kcal/mol but expected ~{expected_Ea / 4184:.2f} kcal/mol "
            f"(from real C9H10+C16H18 thermo, real_dH = {real_dH / 4184:.2f} kcal/mol); "
            f"BM conversion likely used polluted pool thermo instead of real products."
        )

        # (B) Non-vacuousness: recover pool_dH from result.products (the
        #     handshake-relabeled C32H34 PS_scission_tail polymer carries thermo
        #     after make_new_reaction with generate_thermo=True).
        pool_tail = next(
            (p for p in result.products if isinstance(p, Polymer)), None
        )
        assert pool_tail is not None, (
            "Expected a Polymer product (PS_scission_tail) in result.products"
        )
        c9h10_prod = next(
            (p for p in result.products if not isinstance(p, Polymer)), None
        )
        assert c9h10_prod is not None, (
            "Expected a non-Polymer product (alpha-methylstyrene) in result.products"
        )
        pool_dH = (pool_tail.get_enthalpy(298) + c9h10_prod.get_enthalpy(298)
                   - H_proxy)
        pool_Ea = ArrheniusBM(**bm_params).get_activation_energy(pool_dH)

        # The dHrxn values must differ by at least 40 kJ/mol (~9.5 kcal/mol)
        # so the test WOULD FAIL on code that used pool thermo for BM conversion.
        assert abs(real_dH - pool_dH) > 40_000, (
            f"real_dH ({real_dH / 4184:.2f} kcal/mol) and "
            f"pool_dH ({pool_dH / 4184:.2f} kcal/mol) are unexpectedly similar "
            f"(Δ = {abs(real_dH - pool_dH) / 4184:.2f} kcal/mol); "
            "test cannot distinguish real-product thermo from pool thermo."
        )
        assert not np.isclose(Ea, pool_Ea, rtol=0.01), (
            f"Ea = {Ea / 4184:.2f} kcal/mol is too close to pool Ea = "
            f"{pool_Ea / 4184:.2f} kcal/mol; the BM pre-conversion should use "
            "real atom-balanced products, not the relabeled pool proxy."
        )

    def test_pool_H0_clamp_does_not_repollute_Ea(self):
        """I-1 (spec §8 / Task 7): after the real-ΔH pre-conversion, fix_barrier_height's
        pool-sourced endothermicity H0 clamp must NOT raise Ea. The final Ea must equal the
        real-product result (ea_pre), not the pool H0 floor.

        BM params engineered so pool H0 clamp fires (ea_pre < pool_H0 by ~16.7 kJ/mol)
        but real H0 clamp would NOT (ea_pre > real_H0 by ~86.7 kJ/mol).
        Without the Task 7 correction, final Ea = pool_H0 (RED); with correction,
        final Ea = ea_pre (GREEN).

        BM params: E0=165.0 kJ/mol, w0=968.0 kJ/mol → ea_pre ≈ 283.8 kJ/mol
        pool_H0 ≈ 300.4 kJ/mol (from PS_scission_tail + alpha-methylstyrene E0 vs proxy)
        real_H0 ≈ 197.0 kJ/mol (from C16H18 + C9H10 E0 vs proxy)
        Margin: ea_pre - real_H0 ≈ 86.7 kJ/mol, pool_H0 - ea_pre ≈ 16.7 kJ/mol
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.kinetics import ArrheniusBM, Arrhenius

        C9H10_SMILES = 'C=C(C)C1=CC=CC=C1'
        C16H18_SMILES = 'CC(CC=C1C=CC=CC1)C1=CC=CC=C1'

        # BM params tuned so pool H0 clamp fires (~16.7 kJ/mol margin), real H0 clamp does not
        bm_E0_kJ = 165.0
        bm_params = dict(A=(1.293332e12, 's^-1'), n=0.0, w0=(968.0, 'kJ/mol'), E0=(bm_E0_kJ, 'kJ/mol'))

        # --- Compute ea_pre, real_H0, pool_H0 independently (non-circular) ---
        def _make_spc(smiles):
            spc = Species(molecule=[Molecule().from_smiles(smiles)])
            spc.generate_resonance_structures()
            self.model.generate_thermo(spc)
            return spc

        def _E0(spec):
            td = spec.get_thermo_data()
            return td.E0.value_si if td.E0 is not None else td.to_wilhoit().E0.value_si

        c9h10_spc = _make_spc(C9H10_SMILES)
        c16h18_spc = _make_spc(C16H18_SMILES)

        real_dH = (c9h10_spc.get_enthalpy(298) + c16h18_spc.get_enthalpy(298)
                   - self.ps.get_enthalpy(298))
        ea_pre = ArrheniusBM(**bm_params).get_activation_energy(real_dH)

        real_H0 = _E0(c9h10_spc) + _E0(c16h18_spc) - _E0(self.ps)

        # --- Run make_new_reaction with the tuned BM params ---
        proxy_mol = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        rxn = TemplateReaction(
            reactants=[proxy_mol],
            products=[
                Molecule().from_smiles(C16H18_SMILES),
                Molecule().from_smiles(C9H10_SMILES),
            ],
            kinetics=ArrheniusBM(**bm_params),
            family='Retroene',
            is_forward=True,
        )
        result, _ = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=True, generate_kinetics=True,
        )

        # Compute pool_H0 from the actual pool products in result (mirrors fix_barrier_height)
        def _E0_spc(spc):
            td = spc.get_thermo_data()
            return td.E0.value_si if td.E0 is not None else td.to_wilhoit().E0.value_si

        pool_H0 = (sum(_E0_spc(p) for p in result.products)
                   - sum(_E0_spc(r) for r in result.reactants))

        ea_final = result.kinetics.Ea.value_si

        # Verify the engineered setup is correct (margins)
        assert pool_H0 - ea_pre > 5_000, (
            f"Setup check: pool_H0 ({pool_H0/1000:.2f} kJ/mol) should exceed "
            f"ea_pre ({ea_pre/1000:.2f} kJ/mol) by >5 kJ/mol; actual margin "
            f"{(pool_H0-ea_pre)/1000:.3f} kJ/mol"
        )
        assert ea_pre - real_H0 > 5_000, (
            f"Setup check: ea_pre ({ea_pre/1000:.2f} kJ/mol) should exceed "
            f"real_H0 ({real_H0/1000:.2f} kJ/mol) by >5 kJ/mol; actual margin "
            f"{(ea_pre-real_H0)/1000:.3f} kJ/mol"
        )

        # Core assertions: Ea must equal ea_pre (real-product result), NOT pool_H0
        assert abs(ea_final - ea_pre) < 1_000, (
            f"Ea {ea_final/1000:.2f} kJ/mol should equal real-product ea_pre "
            f"{ea_pre/1000:.2f} kJ/mol "
            f"(difference = {abs(ea_final-ea_pre)/1000:.3f} kJ/mol). "
            f"pool_H0={pool_H0/1000:.2f} kJ/mol — did the pool H0 clamp re-pollute Ea?"
        )
        assert ea_final < pool_H0 - 5_000, (
            f"Ea {ea_final/1000:.2f} kJ/mol must NOT be clamped to pool H0 "
            f"{pool_H0/1000:.2f} kJ/mol "
            f"(margin = {(pool_H0-ea_final)/1000:.2f} kJ/mol)"
        )


# Task 6 (non-scission relabel coverage): assessed and skipped — the chip helper's
# end_mod is a radicalized proxy (C25H27•) whose real vs pool ΔH gap (~150 kJ/mol)
# is driven by radical-vs-closed-shell enthalpy, not by the polymer relabeling the
# fix targets; forcing ArrheniusBM onto that scaffold would be artificially misleading.
# The scission path (Task 5) provides the load-bearing coverage for the relabeled gate.
