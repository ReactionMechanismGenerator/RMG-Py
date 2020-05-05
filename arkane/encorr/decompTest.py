#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
This script contains unit tests for the :mod:`arkane.encorr.decomp` module.
"""

import unittest
from collections import Counter

from arkane.encorr.decomp import substruct_decomp, get_substructs, _smi_to_mol


class TestBAC(unittest.TestCase):
    """
    A class for testing that the decomp module functions properly.
    """

    @classmethod
    def setUpClass(cls):
        cls.smi = ('[CH3:1][O:2][C:3](=[C:4]1[CH:5]2[O:6][CH:7]3[N:8]1[CH:9]=[CH:10][CH:11]23)'
                   '[S:12](=[O:13])(=[O:14])[Cl:15]')
        cls.mol = _smi_to_mol(cls.smi)
        cls.idx_map = {atom.GetIdx(): atom.GetAtomMapNum() for atom in cls.mol.GetAtoms()}

    def test_substruct_decomp(self):
        """
        Test that an RDKit molecule can be correctly decomposed into a
        list of lists of atom indices where each sublist contains a
        substructure.
        """
        map_num_list_sorted = [
            [1, 2],                      # C-O bond
            [2],                         # O center with 2 neighbors
            [2, 3],                      # O-C bond
            [3],                         # C center with 3 neigbors
            [3, 4],                      # C=C bond
            [3, 12],                     # C-S bond
            [4, 5, 6, 7, 8, 9, 10, 11],  # Bridged 3-ring complex
            [12],                        # S center with 4 neighbors
            [12, 13],                    # S=O bond
            [12, 14],                    # S=O bond
            [12, 15],                    # S-Cl bond
        ]

        substruct_idxs = substruct_decomp(self.mol)
        substruct_map_nums_sorted = sorted(
            sorted(self.idx_map[idx] for idx in substruct) for substruct in substruct_idxs
        )

        self.assertListEqual(substruct_map_nums_sorted, map_num_list_sorted)

    def test_get_substructs(self):
        """
        Test that the correct substructure counts are returned by
        get_substructs.
        """
        substruct_counts = Counter({  # Canonical RDKit SMILES
            'O=S': 2,
            'CO': 2,
            'SCl': 1,
            'CS': 1,
            'C=C': 1,
            'C1=CN2CC3OC2C13': 1,
            'CS(=O)(=O)Cl': 1,
            'C=C(O)S': 1,
            'COC': 1
        })

        self.assertEqual(get_substructs(self.smi), substruct_counts)
