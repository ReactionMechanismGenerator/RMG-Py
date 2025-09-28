#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import unittest
from numpy import testing

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.vf2 import VF2
from rmgpy.molecule.graph import getVertexConnectivityValue
################################################################################

class TestVF2(unittest.TestCase):
    """
    Contains unit tests of the methods for computing symmetry numbers for a
    given Molecule object.
    """

    def setUp(self):
        self.vf2 = VF2()
        self.mol = Molecule().fromSMILES("CC(=O)C[CH2]")
        self.mol2 = self.mol.copy(deep=True)

    def test_import_graph(self):
        """Test that we can add graphs to the object and that they are sorted"""

        self.mol.sortVertices()
        
        ordered_original_connectivity_order = [getVertexConnectivityValue(atom) for atom in self.mol.atoms]

        self.vf2.graphA = self.mol
        self.vf2.graphB = self.mol2
        
        final_connectivity_order = [getVertexConnectivityValue(atom) for atom in self.vf2.graphA.atoms]
        final_connectivity_order2 = [getVertexConnectivityValue(atom) for atom in self.vf2.graphB.atoms]

        testing.assert_array_equal(final_connectivity_order,final_connectivity_order2)
        testing.assert_array_equal(final_connectivity_order,ordered_original_connectivity_order)

    def test_feasible(self):
        """
        Test that feasibility returns correct values on highly functional molecule
        
        `feasible` method isn't perfect in assigning values but it should do a good
        job on highly functional values
        """
        
        self.vf2.graphA = self.mol
        self.vf2.graphB = self.mol2

        for atom1 in self.vf2.graphA.atoms:
            for atom2 in self.vf2.graphB.atoms:
                # same connectivity values should result in `feasible` being true
                if getVertexConnectivityValue(atom1) == getVertexConnectivityValue(atom2):
                    self.assertTrue(self.vf2.feasible(atom1, atom2))
                else: # different connectivity values should return false
                    self.assertFalse(self.vf2.feasible(atom1, atom2))

    def test_clear_mapping(self):
        """Test that vertex mapping is cleared after isomorphism."""
        self.vf2.isIsomorphic(self.mol, self.mol2, None)

        for atom in self.mol.atoms:
            self.assertIsNone(atom.mapping)
            self.assertFalse(atom.terminal)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
