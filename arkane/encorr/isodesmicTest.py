#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

from rmgpy.molecule import Molecule
from rmgpy.species import Species

from arkane.encorr.isodesmic import ErrorCancelingSpecies, ErrorCancelingReaction

################################################################################


class TestErrorCancelingReactionAndSpecies(unittest.TestCase):
    """
    Tests that ErrorCancelingReaction objects and ErrorCancelingSpecies object are properly implemented
    """

    @classmethod
    def setUpClass(cls):
        """
        A method called before each unit test in this class.
        """
        cls.molecule1 = Molecule(smiles='CC')
        cls.molecule2 = Molecule(smiles='[CH3]')
        cls.species = Species(smiles='CC')

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
            ErrorCancelingReaction(ethane, {methyl: 2})


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
