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

import unittest
import molecule

from .resonance import generate_resonance_structures
from .filtration import get_octet_deviation_list, get_octet_deviation, filter_structures,charge_filtration, get_charge_span_list
################################################################################


class FiltrationTest(unittest.TestCase):

    def basic_filtration_test(self):
        """Test that structures with higher octet deviation get filtered out"""
        adj1 = """
        multiplicity 2
        1 N u0 p1 c0 {2,D} {3,S}
        2 O u0 p2 c0 {1,D}
        3 O u1 p2 c0 {1,S}
        """
        adj2 = """
        multiplicity 2
        1 N u1 p1 c0 {2,S} {3,S}
        2 O u0 p2 c+1 {1,S}
        3 O u0 p3 c-1 {1,S}
        """
        adj3 = """
        multiplicity 2
        1 N u1 p0 c0 {2,T} {3,S}
        2 O u0 p1 c+1 {1,T}
        3 O u0 p3 c-1 {1,S}
        """

        mol1 = molecule.Molecule().fromAdjacencyList(adj1)
        mol2 = molecule.Molecule().fromAdjacencyList(adj2)
        mol3 = molecule.Molecule().fromAdjacencyList(adj3)

        mol_list = [mol1,mol2,mol3]
        octet_deviation_list = get_octet_deviation_list(mol_list)
        filtered_list = filter_structures(mol_list)

        self.assertEqual(octet_deviation_list,[1,3,3])
        self.assertEqual(len(filtered_list), 1)
        self.assertTrue(all([atom.charge == 0 for atom in filtered_list[0].vertices]))

    def penalty_for_N_val_9_test(self):
        """Test that N atoms with valance 9 get penalized in the octet deviation score"""
        adj = """
        multiplicity 2
        1 N u1 p0 c0 {2,S} {3,T}
        2 O u0 p2 c0 {1,S} {4,S}
        3 N u0 p1 c0 {1,T}
        4 H u0 p0 c0 {2,S}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        octet_deviation = get_octet_deviation(mol)
        self.assertEqual(octet_deviation, 2)

    def penalty_for_O4tc_test(self):
        """Test that an O4tc atomType with octet 8 gets penalized in the octet deviation score"""
        adj = """
        1 S u0 p1 c0 {2,S} {3,T}
        2 O u0 p3 c-1 {1,S}
        3 O u0 p1 c+1 {1,T}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        octet_deviation = get_octet_deviation(mol)
        self.assertEqual(octet_deviation, 1)
        self.assertEqual(mol.vertices[2].atomType.label, 'O4tc')

    def penalty_birads_replacing_lone_pairs_test(self):
        """Test that birads on `S u2 p0` are penalized"""
        adj = """
        multiplicity 3
        1 S u2 p0 c0 {2,D} {3,D}
        2 O u0 p2 c0 {1,D}
        3 O u0 p2 c0 {1,D}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        mol.update()
        mol_list = generate_resonance_structures(mol, keep_isomorphic=False, filter_structures=True)
        for mol in mol_list:
            if mol.reactive:
                for atom in mol.vertices:
                    if atom.isSulfur():
                        self.assertNotEquals(atom.radicalElectrons, 2)
        self.assertEqual(len(mol_list), 3)
        self.assertEqual(sum([1 for mol in mol_list if mol.reactive]), 2)

    def penalty_for_s_triple_s_test(self):
        """Test that an S#S substructure in a molecule gets penalized in the octet deviation score"""
        adj = """
        1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
        2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
        3  S u0 p0 c0 {1,S} {4,T} {11,D}
        4  S u0 p1 c0 {2,S} {3,T}
        5  H u0 p0 c0 {1,S}
        6  H u0 p0 c0 {1,S}
        7  H u0 p0 c0 {1,S}
        8  H u0 p0 c0 {2,S}
        9  H u0 p0 c0 {2,S}
        10 H u0 p0 c0 {2,S}
        11 O u0 p2 c0 {3,D}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        octet_deviation = get_octet_deviation(mol)
        self.assertEqual(octet_deviation, 1)

    def electronegativity_test(self):
        """Test that structures with charge separation are only kept if they obey the electronegativity rule

        (If a structure must have charge separation, negative charges will be assigned to more electronegative atoms,
        whereas positive charges will be assigned to less electronegative atoms)

        In this test, only the three structures with no charge separation and the structure where both partial charges
        are on the nitrogen atoms should be kept."""
        adj1 = """
        multiplicity 2
        1 N u0 p1 c0 {2,S} {3,D}
        2 N u1 p1 c0 {1,S} {4,S}
        3 S u0 p1 c0 {1,D} {5,D}
        4 H u0 p0 c0 {2,S}
        5 O u0 p2 c0 {3,D}
        """
        adj2 = """
        multiplicity 2
        1 N u0 p1 c0 {2,D} {3,S}
        2 N u0 p1 c0 {1,D} {4,S}
        3 S u1 p1 c0 {1,S} {5,D}
        4 H u0 p0 c0 {2,S}
        5 O u0 p2 c0 {3,D}
        """
        adj3 = """
        multiplicity 2
        1 N u0 p1 c0 {2,D} {3,S}
        2 N u0 p1 c0 {1,D} {4,S}
        3 S u0 p2 c0 {1,S} {5,S}
        4 H u0 p0 c0 {2,S}
        5 O u1 p2 c0 {3,S}
        """
        adj4 = """
        multiplicity 2
        1 N u1 p0 c+1 {2,S} {3,D}
        2 N u0 p2 c-1 {1,S} {4,S}
        3 S u0 p1 c0 {1,D} {5,D}
        4 H u0 p0 c0 {2,S}
        5 O u0 p2 c0 {3,D}
        """
        adj5 = """
        multiplicity 2
        1 N u1 p0 c+1 {2,D} {3,S}
        2 N u0 p1 c0 {1,D} {4,S}
        3 S u0 p2 c-1 {1,S} {5,D}
        4 H u0 p0 c0 {2,S}
        5 O u0 p2 c0 {3,D}
        """
        adj6 = """
        multiplicity 2
        1 N u1 p1 c0 {2,S} {3,S}
        2 N u0 p2 c-1 {1,S} {4,S}
        3 S u0 p1 c+1 {1,S} {5,D}
        4 H u0 p0 c0 {2,S}
        5 O u0 p2 c0 {3,D}
        """
        adj7 = """
        multiplicity 2
        1 N u1 p0 c+1 {2,D} {3,S}
        2 N u0 p1 c0 {1,D} {4,S}
        3 S u0 p2 c0 {1,S} {5,S}
        4 H u0 p0 c0 {2,S}
        5 O u0 p3 c-1 {3,S}
        """

        mol_list = [molecule.Molecule().fromAdjacencyList(adj1),
        molecule.Molecule().fromAdjacencyList(adj2),
        molecule.Molecule().fromAdjacencyList(adj3),
        molecule.Molecule().fromAdjacencyList(adj4),
        molecule.Molecule().fromAdjacencyList(adj5),
        molecule.Molecule().fromAdjacencyList(adj6),
        molecule.Molecule().fromAdjacencyList(adj7)]

        for mol in mol_list:
            mol.update()  # the charge_filtration uses the atom.sortingLabel attribute

        filtered_list = charge_filtration(mol_list, get_charge_span_list(mol_list))
        self.assertEqual(len(filtered_list), 4)
        self.assertTrue(any([mol.getChargeSpan() == 1 for mol in filtered_list]))
        for mol in filtered_list:
            if mol.getChargeSpan() == 1:
                for atom in mol.vertices:
                    if abs(atom.charge) == 1:
                        self.assertTrue(atom.isNitrogen())
