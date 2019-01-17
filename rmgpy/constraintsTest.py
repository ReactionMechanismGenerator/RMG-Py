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

"""
This script contains unit tests of the :mod:`rmgpy.constraints` module.
"""

import unittest
import mock

from rmgpy.rmg.main import RMG
from rmgpy.constraints import failsSpeciesConstraints
from rmgpy.species import Species
from rmgpy.molecule import Molecule
import rmgpy.rmg.input

################################################################################

class TestFailsSpeciesConstraints(unittest.TestCase):
    """
    Contains unit tests of the failsSpeciesConstraints function.
    """

    @classmethod
    def setUpClass(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        cls.rmg = RMG()
        rmgpy.rmg.input.rmg = cls.rmg
        rmgpy.rmg.input.generatedSpeciesConstraints(
            maximumCarbonAtoms=2,
            maximumOxygenAtoms=1,
            maximumNitrogenAtoms=1,
            maximumSiliconAtoms=1,
            maximumSulfurAtoms=1,
            maximumHeavyAtoms=3,
            maximumRadicalElectrons=2,
            maximumSingletCarbenes=1,
            maximumCarbeneRadicals=0,
        )

    @classmethod
    def tearDownClass(cls):
        """
        A function run ONCE after all unit tests in this class.
        """
        rmgpy.rmg.input.rmg = None

    @mock.patch('rmgpy.constraints.logging')
    def testConstraintsNotLoaded(self, mock_logging):
        """
        Test what happens when constraints are not loaded.
        """
        # Reset module level rmg variable in rmgpy.rmg.input
        rmgpy.rmg.input.rmg = None

        mol = Molecule(SMILES='C')

        self.assertFalse(failsSpeciesConstraints(mol))

        mock_logging.debug.assert_called_with('Species constraints could not be found.')

        # Restore module level rmg variable in rmgpy.rmg.input
        rmgpy.rmg.input.rmg = self.rmg

    def testSpeciesInput(self):
        """
        Test that failsSpeciesConstraints can handle a Species object.
        """
        spc = Species().fromSMILES('C')

        self.assertFalse(failsSpeciesConstraints(spc))

    def testExplicitlyAllowedMolecules(self):
        """
        Test that we can explicitly allow molecules in species constraints.
        """
        mol = Molecule(SMILES='CCCC')
        self.assertTrue(failsSpeciesConstraints(mol))

        self.rmg.speciesConstraints['explicitlyAllowedMolecules'] = [Molecule(SMILES='CCCC')]
        self.assertFalse(failsSpeciesConstraints(mol))

    def testCarbonConstraint(self):
        """
        Test that we can constrain the max number of carbon atoms.
        """
        mol1 = Molecule(SMILES='CC')
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule(SMILES='CCC')
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testOxygenConstraint(self):
        """
        Test that we can constrain the max number of oxygen atoms.
        """
        mol1 = Molecule(SMILES='C=O')
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule(SMILES='OC=O')
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testNitrogenConstraint(self):
        """
        Test that we can constrain the max number of nitrogen atoms.
        """
        mol1 = Molecule(SMILES='CN')
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule(SMILES='NCN')
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testSiliconConstraint(self):
        """
        Test that we can constrain the max number of silicon atoms.
        """
        mol1 = Molecule(SMILES='[SiH4]')
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule(SMILES='[SiH3][SiH3]')
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testSulfurConstraint(self):
        """
        Test that we can constrain the max number of sulfur atoms.
        """
        mol1 = Molecule(SMILES='CS')
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule(SMILES='SCS')
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testHeavyConstraint(self):
        """
        Test that we can constrain the max number of heavy atoms.
        """
        mol1 = Molecule(SMILES='CCO')
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule(SMILES='CCN=O')
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testRadicalConstraint(self):
        """
        Test that we can constrain the max number of radical electrons.
        """
        mol1 = Molecule(SMILES='[CH2][CH2]')
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule(SMILES='[CH2][CH][CH2]')
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testCarbeneConstraint(self):
        """
        Test that we can constrain the max number of singlet carbenes.
        """
        mol1 = Molecule().fromAdjacencyList("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""")
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule().fromAdjacencyList("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 C u0 p1 c0 {1,S} {4,S}
4 H u0 p0 c0 {3,S}
""")
        self.assertTrue(failsSpeciesConstraints(mol2))

    def testCarbeneRadicalConstraint(self):
        """
        Test that we can constrain the max number of radical electrons with a carbene.
        """
        mol1 = Molecule().fromAdjacencyList("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""")
        self.assertFalse(failsSpeciesConstraints(mol1))

        mol2 = Molecule().fromAdjacencyList("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
""")
        self.assertTrue(failsSpeciesConstraints(mol2))
