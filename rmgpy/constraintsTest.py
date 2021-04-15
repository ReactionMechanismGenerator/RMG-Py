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

"""
This script contains unit tests of the :mod:`rmgpy.constraints` module.
"""

import unittest
from unittest import mock

from rmgpy.rmg.main import RMG
from rmgpy.constraints import fails_species_constraints
from rmgpy.species import Species
from rmgpy.molecule import Molecule
import rmgpy.rmg.input


################################################################################

class TestFailsSpeciesConstraints(unittest.TestCase):
    """
    Contains unit tests of the fails_species_constraints function.
    """

    @classmethod
    def setUpClass(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        cls.rmg = RMG()
        rmgpy.rmg.input.rmg = cls.rmg
        rmgpy.rmg.input.generated_species_constraints(
            maximumCarbonAtoms=2,
            maximumOxygenAtoms=1,
            maximumNitrogenAtoms=1,
            maximumSiliconAtoms=1,
            maximumSulfurAtoms=1,
            maximumSurfaceSites=2,
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
    def test_constraints_not_loaded(self, mock_logging):
        """
        Test what happens when constraints are not loaded.
        """
        # Reset module level rmg variable in rmgpy.rmg.input
        rmgpy.rmg.input.rmg = None

        mol = Molecule(smiles='C')

        self.assertFalse(fails_species_constraints(mol))

        mock_logging.debug.assert_called_with('Species constraints could not be found.')

        # Restore module level rmg variable in rmgpy.rmg.input
        rmgpy.rmg.input.rmg = self.rmg

    def test_species_input(self):
        """
        Test that fails_species_constraints can handle a Species object.
        """
        spc = Species().from_smiles('C')

        self.assertFalse(fails_species_constraints(spc))

    def test_explicitly_allowed_molecules(self):
        """
        Test that we can explicitly allow molecules in species constraints.
        """
        mol = Molecule(smiles='CCCC')
        self.assertTrue(fails_species_constraints(mol))

        self.rmg.species_constraints['explicitlyAllowedMolecules'] = [Molecule(smiles='CCCC')]
        self.assertFalse(fails_species_constraints(mol))

    def test_carbon_constraint(self):
        """
        Test that we can constrain the max number of carbon atoms.
        """
        mol1 = Molecule(smiles='CC')
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule(smiles='CCC')
        self.assertTrue(fails_species_constraints(mol2))

    def test_oxygen_constraint(self):
        """
        Test that we can constrain the max number of oxygen atoms.
        """
        mol1 = Molecule(smiles='C=O')
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule(smiles='OC=O')
        self.assertTrue(fails_species_constraints(mol2))

    def test_nitrogen_constraint(self):
        """
        Test that we can constrain the max number of nitrogen atoms.
        """
        mol1 = Molecule(smiles='CN')
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule(smiles='NCN')
        self.assertTrue(fails_species_constraints(mol2))

    def test_silicon_constraint(self):
        """
        Test that we can constrain the max number of silicon atoms.
        """
        mol1 = Molecule(smiles='[SiH4]')
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule(smiles='[SiH3][SiH3]')
        self.assertTrue(fails_species_constraints(mol2))

    def test_sulfur_constraint(self):
        """
        Test that we can constrain the max number of sulfur atoms.
        """
        mol1 = Molecule(smiles='CS')
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule(smiles='SCS')
        self.assertTrue(fails_species_constraints(mol2))

    def test_surface_site_constraint(self):
        """
        Test that we can constrain the max number of surface sites.
        """

        mol_1site = Molecule().from_adjacency_list("""
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 X u0 p0 c0 {2,D}
""")
        mol_2site = Molecule().from_adjacency_list("""
1 C u0 p0 c0 {2,D} {3,D}
2 C u0 p0 c0 {1,D} {4,D}
3 X u0 p0 c0 {1,D}
4 X u0 p0 c0 {2,D}
""")

        mol_3site_vdW = Molecule().from_adjacency_list("""
1 C u0 p0 c0 {2,D} {3,D}
2 C u0 p0 c0 {1,D} {4,D}
3 X u0 p0 c0 {1,D}
4 X u0 p0 c0 {2,D}
6 X u0 p0 c0
""")

        mol_3site = Molecule().from_adjacency_list("""
1 C u0 p0 c0 {4,S} {2,D} {7,S}
2 C u0 p0 c0 {1,D} {3,S} {8,S}
3 C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 X u0 p0 c0 {1,S}
8 X u0 p0 c0 {2,S}
9 X u0 p0 c0 {3,S}
""")
        max_carbon = self.rmg.species_constraints['maximumCarbonAtoms']
        max_heavy_atoms = self.rmg.species_constraints['maximumHeavyAtoms']

        self.rmg.species_constraints['maximumCarbonAtoms'] = 3
        self.rmg.species_constraints['maximumHeavyAtoms'] = 6

        self.assertFalse(fails_species_constraints(mol_1site))
        self.assertFalse(fails_species_constraints(mol_2site))

        self.assertTrue(fails_species_constraints(mol_3site_vdW))
        self.assertTrue(fails_species_constraints(mol_3site))

        self.rmg.species_constraints['maximumCarbonAtoms'] = max_carbon
        self.rmg.species_constraints['maximumHeavyAtoms'] = max_heavy_atoms

    def test_heavy_constraint(self):
        """
        Test that we can constrain the max number of heavy atoms.
        """
        mol1 = Molecule(smiles='CCO')
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule(smiles='CCN=O')
        self.assertTrue(fails_species_constraints(mol2))

    def test_radical_constraint(self):
        """
        Test that we can constrain the max number of radical electrons.
        """
        mol1 = Molecule(smiles='[CH2][CH2]')
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule(smiles='[CH2][CH][CH2]')
        self.assertTrue(fails_species_constraints(mol2))

    def test_carbene_constraint(self):
        """
        Test that we can constrain the max number of singlet carbenes.
        """
        mol1 = Molecule().from_adjacency_list("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""")
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule().from_adjacency_list("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 C u0 p1 c0 {1,S} {4,S}
4 H u0 p0 c0 {3,S}
""")
        self.assertTrue(fails_species_constraints(mol2))

    def test_carbene_radical_constraint(self):
        """
        Test that we can constrain the max number of radical electrons with a carbene.
        """
        mol1 = Molecule().from_adjacency_list("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""")
        self.assertFalse(fails_species_constraints(mol1))

        mol2 = Molecule().from_adjacency_list("""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
""")
        self.assertTrue(fails_species_constraints(mol2))
