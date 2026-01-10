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

from rmgpy.polymer import Polymer
from rmgpy.molecule import Molecule
from rmgpy.species import Species


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
        ps_smiles = '[CH2][CH]c1ccccc1'
        pe_smiles = '[CH2][CH2]'
        self.polymer_1 = Polymer(
                 label='PS_1',
                 monomer=ps_adj,
                 end_groups=['C[C](C)c1ccccc1', '[H]'],
                 cutoff=3,
                 Mn=5000.0,
                 Mw=6000.0,
                 initial_mass=1.0,
        )

        self.polymer_2 = Polymer(
                 label='PS_2',
                 monomer=ps_smiles,
                 end_groups=['C[C](C)c1ccccc1', '[H]'],
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
        adj = self.polymer_1.to_adjacency_list()
        partial_expected_adj = """PS_1
multiplicity 3
1  *1 C u1 p0 c0 {2,S} {9,S} {10,S}
2  *2 C u1 p0 c0 {1,S} {3,S} {11,S}
3     C u0 p0 c0 {2,S} {4,B} {8,B}
4     C u0 p0 c0 {3,B} {5,B} {12,S}
5     C u0 p0 c0 {4,B} {6,B} {13,S}
6     C u0 p0 c0 {5,B} {7,B} {14,S}
7     C u0 p0 c0 {6,B} {8,B} {15,S}
8     C u0 p0 c0 {3,B} {7,B} {16,S}
9     H u0 p0 c0 {1,S}
10    H u0 p0 c0 {1,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
14    H u0 p0 c0 {6,S}
15    H u0 p0 c0 {7,S}
16    H u0 p0 c0 {8,S}


PS_1
multiplicity 3
1  *1 C u1 p0 c0 {2,S} {9,S} {10,S}
2  *2 C u0 p0 c0 {1,S} {3,D} {11,S}
3     C u0 p0 c0 {2,D} {4,S} {8,S}
4     C u0 p0 c0 {3,S} {5,D} {12,S}
5     C u0 p0 c0 {4,D} {6,S} {13,S}
6     C u0 p0 c0 {5,S} {7,D} {14,S}
7     C u0 p0 c0 {6,D} {8,S} {15,S}
8     C u1 p0 c0 {3,S} {7,S} {16,S}
9     H u0 p0 c0 {1,S}
10    H u0 p0 c0 {1,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
14    H u0 p0 c0 {6,S}
15    H u0 p0 c0 {7,S}
16    H u0 p0 c0 {8,S}"""
        assert partial_expected_adj in adj

    def test_to_smiles(self):
        """
        Test that to_smiles() works as expected.
        """
        smiles = self.polymer_1.monomers[0].to_smiles()
        expected_smiles = '[CH2][CH]c1ccccc1'
        assert smiles == expected_smiles

    def test_copy(self):
        """Test that we can make a copy of a Polymer object."""
        poly_cp = self.polymer_1.copy()
        assert id(self.polymer_1) != id(poly_cp)
        assert self.polymer_1.is_isomorphic(poly_cp)
        assert self.polymer_1.label == poly_cp.label
        assert self.polymer_1.index == poly_cp.index

    def test_fingerprint_property(self):
        """Test that the fingerprint property works"""
        assert self.polymer_1.fingerprint == "Polymer_5000.0_6000.0_3_C08H08N00O00S00"
        assert self.polymer_2.fingerprint == "Polymer_3000.0_10000.0_5_C08H08N00O00S00"
        assert self.polymer_3.fingerprint == "Polymer_1000.0_2500.0_10_C02H04N00O00S00"

    def test_is_isomorphic(self):
        """Test that the Polymer.is_isomorphic works"""
        spc_1 = Species(smiles='[CH2][CH]c1ccccc1')
        mol_1 = Molecule(smiles='[CH2][CH]c1ccccc1')
        pol_1 = Polymer(
                 label='PS_10',
                 monomer='[CH2][CH]c1ccccc1',
                 end_groups=['C[C](C)c1ccccc1', '[H]'],
                 cutoff=30,
                 Mn=7000.0,
                 Mw=8000.0,
                 initial_mass=10.0,
        )
        assert not spc_1.is_isomorphic(pol_1)
        assert not pol_1.is_isomorphic(spc_1)
        assert not pol_1.is_isomorphic(mol_1)
        assert self.polymer_1.is_isomorphic(self.polymer_1)
        assert self.polymer_1.is_isomorphic(pol_1)
        assert any(spc_1.is_isomorphic(monomer) for monomer in pol_1.monomers)

    def test_polymer_label(self):
        """Test that the polymer label"""
        assert self.polymer_1.label == "PS_1"
        assert self.polymer_2.label == "PS_2"
        assert self.polymer_3.label == "PE_1"

    def test_monomer_atom_labels(self):
        """Test that the monomer atom labels are assigned correctly"""
        for polymer in [self.polymer_1, self.polymer_2, self.polymer_3]:
            assert polymer.monomer_head is not None
            assert polymer.monomer_tail is not None
            assert polymer.monomer_head != polymer.monomer_tail
            head_atom = polymer.monomers[0].atoms[polymer.monomer_head]
            tail_atom = polymer.monomers[0].atoms[polymer.monomer_tail]
            assert head_atom.is_carbon()
            assert tail_atom.is_carbon()
            assert head_atom.radical_electrons == 1
            assert tail_atom.radical_electrons == 1
            assert head_atom is not tail_atom

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

