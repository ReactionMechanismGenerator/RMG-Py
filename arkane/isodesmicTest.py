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

"""
This script contains unit tests of the :mod:`arkane.isodesmic` module.
"""

import unittest

import numpy as np

from arkane.isodesmic import ErrorCancelingSpecies, ErrorCancelingReaction, SpeciesConstraints
from rmgpy.molecule import Molecule
from rmgpy.species import Species

################################################################################


class TestErrorCancelingReactionAndSpecies(unittest.TestCase):
    """
    Tests that ErrorCancelingReaction objects and ErrorCancelingSpecies object are properly implemented
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.molecule1 = Molecule(SMILES='CC')
        self.molecule2 = Molecule(SMILES='[CH3]')
        self.species = Species(SMILES='CC')

    def test_error_canceling_species(self):
        """
        Test that ErrorCancelingSpecies can be created properly
        """
        error_canceling_species = ErrorCancelingSpecies(self.molecule1, (123.4, 'kcal/mol'), 'test', (100.0, 'kJ/mol'))
        self.assertIsInstance(error_canceling_species, ErrorCancelingSpecies)
        self.assertAlmostEqual(error_canceling_species.low_level_hf298.value_si, 123.4*4184)

        # For target species the high level data is not given
        target_species = ErrorCancelingSpecies(self.molecule2, (10.1, 'J/mol'), 'test')
        self.assertIs(target_species.high_level_hf298, None)

    def test_molecule_input_in_error_canceling_species(self):
        """
        Test that an exception is raised if an rmgpy Molecule object is not passed to an ErrorCancelingSpecies
        """
        with self.assertRaises(ValueError):
            ErrorCancelingSpecies(self.species, (100.0, 'J/mol'), 'test')

    def test_error_canceling_reactions(self):
        """
        Test that ErrorCancelingReaction object can be created and that hf298 can be calculated for the target
        """
        # Take ethane as the target
        ethane = ErrorCancelingSpecies(self.molecule1, (100.0, 'kJ/mol'), 'test')
        methyl = ErrorCancelingSpecies(self.molecule2, (20.0, 'kcal/mol'), 'test', (21000.0, 'cal/mol'))

        # This reaction is not an isodesmic reaction, but that does not matter for the unit test
        rxn = ErrorCancelingReaction(ethane, {methyl: 2})
        self.assertAlmostEqual(rxn.calculate_target_thermo().value_si, 2*21000.0*4.184-(2*20.0*4184-100.0*1000))

    def test_model_chemistry_consistency(self):
        """
        Test that ErrorCancelingReaction objects properly check that all species use the same model chemistry
        """
        # Take ethane as the target
        ethane = ErrorCancelingSpecies(self.molecule1, (100.0, 'kJ/mol'), 'test_A')
        methyl = ErrorCancelingSpecies(self.molecule2, (20.0, 'kcal/mol'), 'test_B', (21000.0, 'cal/mol'))

        # This should throw an exception because the model chemistry is different
        with self.assertRaises(ValueError):
            _ = ErrorCancelingReaction(ethane, {methyl: 2})


class TestSpeciesConstraints(unittest.TestCase):
    """
    A class for testing that the SpeciesConstraint class functions properly
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        # Give all species a low level Hf298 of 100 J/mol--this is not important for this test
        hf = (100.0, 'J/mol')

        self.propene = ErrorCancelingSpecies(Molecule(SMILES='CC=C'), hf, 'test')
        self.butane = ErrorCancelingSpecies(Molecule(SMILES='CCCC'), hf, 'test')
        self.benzene = ErrorCancelingSpecies(Molecule(SMILES='c1ccccc1'), hf, 'test')
        self.caffeine = ErrorCancelingSpecies(Molecule(SMILES='CN1C=NC2=C1C(=O)N(C(=O)N2C)C'), hf, 'test')
        self.ethyne = ErrorCancelingSpecies(Molecule(SMILES='C#C'), hf, 'test')

    def test_initializing_constraint_map(self):
        """
        Test that the constraint map is properly initialized when a SpeciesConstraints object is initialized
        """
        caffeine_consts = SpeciesConstraints(self.caffeine, [self.butane, self.benzene])
        self.assertEqual(caffeine_consts.constraint_map, {'H': 0, 'C': 1, 'O': 2, 'N': 3,
                                                          'C=O': 4, 'C-N': 5, 'C-H': 6, 'C=C': 7, 'C=N': 8, 'C-C': 9,
                                                          '5_ring': 10, '6_ring': 11})

        no_rings = SpeciesConstraints(self.caffeine, [self.butane, self.benzene], conserve_ring_size=False)
        self.assertEqual(no_rings.constraint_map, {'H': 0, 'C': 1, 'O': 2, 'N': 3,
                                                   'C=O': 4, 'C-N': 5, 'C-H': 6, 'C=C': 7, 'C=N': 8, 'C-C': 9})

        atoms_only = SpeciesConstraints(self.caffeine, [self.butane], conserve_ring_size=False, conserve_bonds=False)
        self.assertEqual(atoms_only.constraint_map, {'H': 0, 'C': 1, 'O': 2, 'N': 3})

    def test_enumerating_constraints(self):
        """
        Test that a SpeciesConstraints object can properly enumerate the constraints of a given ErrorCancelingSpecies
        """
        spcs_consts = SpeciesConstraints(self.benzene, [])
        self.assertEqual(set(spcs_consts.constraint_map.keys()), {'C', 'H', 'C=C', 'C-C', 'C-H', '6_ring'})

        # Now that we have confirmed that the correct keys are present, overwrite the constraint map to set the order
        spcs_consts.constraint_map = {'H': 0, 'C': 1, 'C=C': 2, 'C-C': 3, 'C-H': 4, '6_ring': 5}

        self.assertTrue(np.array_equal(spcs_consts._enumerate_constraints(self.propene), np.array([6, 3, 1, 1, 6, 0])))
        self.assertTrue(np.array_equal(spcs_consts._enumerate_constraints(self.butane), np.array([10, 4, 0, 3, 10, 0])))
        self.assertTrue(np.array_equal(spcs_consts._enumerate_constraints(self.benzene), np.array([6, 6, 3, 3, 6, 1])))

        # Caffeine and ethyne should return None since they have features not found in benzene
        self.assertIs(spcs_consts._enumerate_constraints(self.caffeine), None)
        self.assertIs(spcs_consts._enumerate_constraints(self.ethyne), None)

    def test_calculating_constraints(self):
        """
        Test that a SpeciesConstraints object can properly return the target constraint vector and the constraint matrix
        """
        spcs_consts = SpeciesConstraints(self.caffeine, [self.propene, self.butane, self.benzene, self.ethyne])
        self.assertEqual(set(spcs_consts.constraint_map.keys()), {'H', 'C', 'O', 'N', 'C=O', 'C-N', 'C-H', 'C=C', 'C=N',
                                                                  'C-C', '5_ring', '6_ring'})

        # Now that we have confirmed that the correct keys are present, overwrite the constraint map to set the order
        spcs_consts.constraint_map = ({'H': 0, 'C': 1, 'O': 2, 'N': 3,
                                       'C=O': 4, 'C-N': 5, 'C-H': 6, 'C=C': 7, 'C=N': 8, 'C-C': 9,
                                       '5_ring': 10, '6_ring': 11})

        target_consts, consts_matrix = spcs_consts.calculate_constraints()

        # First, test that ethyne is not included in the reference set
        self.assertEqual(spcs_consts.reference_species, [self.propene, self.butane, self.benzene])

        # Then, test the output of the calculation
        self.assertTrue(np.array_equal(target_consts, np.array([10, 8, 2, 4, 2, 10, 10, 1, 1, 1, 1, 1])))
        self.assertTrue(np.array_equal(consts_matrix, np.array([[6, 3, 0, 0, 0, 0, 6, 1, 0, 1, 0, 0],
                                                                [10, 4, 0, 0, 0, 0, 10, 0, 0, 3, 0, 0],
                                                                [6, 6, 0, 0, 0, 0, 6, 3, 0, 3, 0, 1]])))


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
