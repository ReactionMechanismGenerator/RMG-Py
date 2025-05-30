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
This script contains unit tests of the :mod:`arkane.isodesmic` module.
"""

import numpy as np

from rmgpy.molecule import Molecule
from rmgpy.species import Species

from arkane.encorr.isodesmic import (
    ErrorCancelingScheme,
    ErrorCancelingSpecies,
    ErrorCancelingReaction,
    IsodesmicScheme,
    SpeciesConstraints,
)
from arkane.modelchem import LevelOfTheory
import pytest


class TestErrorCancelingReactionAndSpecies:
    """
    Tests that ErrorCancelingReaction objects and ErrorCancelingSpecies object are properly implemented
    """

    @classmethod
    def setup_class(cls):
        """
        A method called before each unit test in this class.
        """
        cls.molecule1 = Molecule(smiles="CC")
        cls.molecule2 = Molecule(smiles="[CH3]")
        cls.species = Species(smiles="CC")

    def test_error_canceling_species(self):
        """
        Test that ErrorCancelingSpecies can be created properly
        """
        lot = LevelOfTheory("test")
        error_canceling_species = ErrorCancelingSpecies(self.molecule1, (123.4, "kcal/mol"), lot, (100.0, "kJ/mol"))
        assert isinstance(error_canceling_species, ErrorCancelingSpecies)
        assert round(abs(error_canceling_species.low_level_hf298.value_si - 123.4 * 4184), 7) == 0

        # For target species the high level data is not given
        target_species = ErrorCancelingSpecies(self.molecule2, (10.1, "J/mol"), lot)
        assert target_species.high_level_hf298 is None

    def test_molecule_input_in_error_canceling_species(self):
        """
        Test that an exception is raised if an rmgpy Molecule object is not passed to an ErrorCancelingSpecies
        """
        with pytest.raises(ValueError):
            ErrorCancelingSpecies(self.species, (100.0, "J/mol"), LevelOfTheory("test"))

    def test_error_canceling_reactions(self):
        """
        Test that ErrorCancelingReaction object can be created and that hf298 can be calculated for the target
        """
        # Take ethane as the target
        lot = LevelOfTheory("test")
        ethane = ErrorCancelingSpecies(self.molecule1, (100.0, "kJ/mol"), lot)
        methyl = ErrorCancelingSpecies(self.molecule2, (20.0, "kcal/mol"), lot, (21000.0, "cal/mol"))

        # This reaction is not an isodesmic reaction, but that does not matter for the unit test
        rxn = ErrorCancelingReaction(ethane, {methyl: 2})
        assert round(abs(rxn.calculate_target_thermo().value_si - (2 * 21000.0 * 4.184 - (2 * 20.0 * 4184 - 100.0 * 1000))), 7) == 0

    def test_level_of_theory_consistency(self):
        """
        Test that ErrorCancelingReaction objects properly check that all species use the same level of theory
        """
        # Take ethane as the target
        ethane = ErrorCancelingSpecies(self.molecule1, (100.0, "kJ/mol"), LevelOfTheory("test_A"))
        methyl = ErrorCancelingSpecies(
            self.molecule2,
            (20.0, "kcal/mol"),
            LevelOfTheory("test_B"),
            (21000.0, "cal/mol"),
        )

        # This should throw an exception because the model chemistry is different
        with pytest.raises(ValueError):
            ErrorCancelingReaction(ethane, {methyl: 2})


class TestSpeciesConstraints:
    """
    A class for testing that the SpeciesConstraint class functions properly
    """

    @classmethod
    def setup_class(cls):
        """
        A method called before each unit test in this class.
        """
        # Give all species a low level Hf298 of 100 J/mol--this is not important for this test
        hf = (100.0, "J/mol")

        lot = LevelOfTheory("test")
        cls.propene = ErrorCancelingSpecies(Molecule(smiles="CC=C"), hf, lot)
        cls.butane = ErrorCancelingSpecies(Molecule(smiles="CCCC"), hf, lot)
        cls.benzene = ErrorCancelingSpecies(Molecule(smiles="c1ccccc1"), hf, lot)
        cls.caffeine = ErrorCancelingSpecies(Molecule(smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C"), hf, lot)
        cls.ethyne = ErrorCancelingSpecies(Molecule(smiles="C#C"), hf, lot)

    def test_initializing_constraint_map(self):
        """
        Test that the constraint map is properly initialized when a SpeciesConstraints object is initialized
        """
        consts = SpeciesConstraints(self.caffeine, [self.butane, self.benzene])
        caffeine_features = consts._get_all_constraints(self.caffeine)
        caffeine_constraint_list = [feat.__repr__() for feat in caffeine_features]

        assert set(caffeine_constraint_list) == {
            "C=O",
            "C-N",
            "C-H",
            "C=C",
            "C=N",
            "C-C",
            "5_ring",
            "6_ring",
        }

        no_rings = SpeciesConstraints(self.caffeine, [self.butane, self.benzene], conserve_ring_size=False)
        caffeine_features = no_rings._get_all_constraints(self.caffeine)
        caffeine_constraint_list = [feat.__repr__() for feat in caffeine_features]
        assert set(caffeine_constraint_list) == {"C=O", "C-N", "C-H", "C=C", "C=N", "C-C"}

    def test_enumerating_constraints(self):
        """
        Test that a SpeciesConstraints object can properly enumerate the constraints of a given ErrorCancelingSpecies
        """
        spcs_consts = SpeciesConstraints(self.benzene, [])
        benzene_features = spcs_consts._get_all_constraints(self.benzene)
        benzene_constraint_list = [feat.__repr__() for feat in benzene_features]
        assert set(benzene_constraint_list) == {"C=C", "C-C", "C-H", "6_ring"}

        target_constraints, _ = spcs_consts._enumerate_constraints([benzene_features])
        benzene_constraints = target_constraints

        assert np.array_equal(
            benzene_constraints,
            np.array([1, 3, 6, 3]), # ['6_ring', 'C-C', 'C-H', 'C=C']
        )

        spcs_consts.all_reference_species = [self.propene]
        propene_features = spcs_consts._get_all_constraints(self.propene)
        _, reference_constraints = spcs_consts._enumerate_constraints([benzene_features, propene_features])
        propene_constraints = reference_constraints[0]
        assert np.array_equal(
            propene_constraints,
            np.array([0, 1, 6, 1]), # ['6_ring', 'C-C', 'C-H', 'C=C']
        )

        spcs_consts.all_reference_species = [self.butane]
        butane_features = spcs_consts._get_all_constraints(self.butane)
        _, reference_constraints = spcs_consts._enumerate_constraints([benzene_features, butane_features])
        butane_constraints = reference_constraints[0]
        assert np.array_equal(
            butane_constraints,
            np.array([0, 3, 10, 0]), # ['6_ring', 'C-C', 'C-H', 'C=C']
        )

        # Caffeine and ethyne should return empty list since they have features not found in benzene
        spcs_consts.all_reference_species = [self.caffeine]
        caffeine_features = spcs_consts._get_all_constraints(self.caffeine)
        _, reference_constraints = spcs_consts._enumerate_constraints([benzene_features, caffeine_features])
        assert len(reference_constraints) == 0

        spcs_consts.all_reference_species = [self.ethyne]
        ethyne_features = spcs_consts._get_all_constraints(self.ethyne)
        _, reference_constraints = spcs_consts._enumerate_constraints([benzene_features, ethyne_features])
        assert len(reference_constraints) == 0

    def test_calculating_constraints(self):
        """
        Test that a SpeciesConstraints object can properly return the target constraint vector and the constraint matrix
        """
        spcs_consts = SpeciesConstraints(self.caffeine, [self.propene, self.butane, self.benzene, self.ethyne])
        caffeine_features = spcs_consts._get_all_constraints(self.caffeine)
        propene_features = spcs_consts._get_all_constraints(self.propene)
        butane_features = spcs_consts._get_all_constraints(self.butane)
        benzene_features = spcs_consts._get_all_constraints(self.benzene)
        ethyne_features = spcs_consts._get_all_constraints(self.ethyne)

        caffeine_feature_list = [feat.__repr__() for feat in caffeine_features]
        assert set(caffeine_feature_list) == {
            "C=O",
            "C-N",
            "C-H",
            "C=C",
            "C=N",
            "C-C",
            "5_ring",
            "6_ring",
        }

        target_consts, consts_matrix = spcs_consts.calculate_constraints()

        # First, test that ethyne is not included in the reference set
        assert spcs_consts.reference_species == [self.propene, self.butane, self.benzene]

        # Then, test the output of the calculation
        assert np.array_equal(target_consts, np.array([ 1,  1,  1, 10, 10,  1,  1,  2,  0,  8, 10,  4,  2])) # ['5_ring', '6_ring', 'C-C', 'C-H', 'C-N', 'C=C', 'C=N', 'C=O']
        assert np.array_equal(
            consts_matrix,
            np.array(
                [
                    [ 0,  0,  1,  6,  0,  1,  0,  0,  0,  3,  6,  0,  0],  # ['5_ring', '6_ring', 'C-C', 'C-H', 'C-N', 'C=C', 'C=N', 'C=O']
                    [ 0,  0,  3, 10,  0,  0,  0,  0,  0,  4, 10,  0,  0],  # ['5_ring', '6_ring', 'C-C', 'C-H', 'C-N', 'C=C', 'C=N', 'C=O']
                    [ 0,  1,  3,  6,  0,  3,  0,  0,  0,  6,  6,  0,  0],  # ['5_ring', '6_ring', 'C-C', 'C-H', 'C-N', 'C=C', 'C=N', 'C=O']
                ]
            ),
        )


class TestErrorCancelingScheme:
    """
    A class for testing that the ErrorCancelingScheme class functions properly
    """

    @classmethod
    def setup_class(cls):
        lot = LevelOfTheory("test")
        cls.propene = ErrorCancelingSpecies(Molecule(smiles="CC=C"), (100, "kJ/mol"), lot, (105, "kJ/mol"))
        cls.propane = ErrorCancelingSpecies(Molecule(smiles="CCC"), (75, "kJ/mol"), lot, (80, "kJ/mol"))
        cls.butane = ErrorCancelingSpecies(Molecule(smiles="CCCC"), (150, "kJ/mol"), lot, (145, "kJ/mol"))
        cls.butene = ErrorCancelingSpecies(Molecule(smiles="C=CCC"), (175, "kJ/mol"), lot, (180, "kJ/mol"))
        cls.pentane = ErrorCancelingSpecies(Molecule(smiles="CCCCC"), (200, "kJ/mol"), lot, (190, "kJ/mol"))
        cls.pentene = ErrorCancelingSpecies(Molecule(smiles="C=CCCC"), (225, "kJ/mol"), lot, (220, "kJ/mol"))
        cls.hexane = ErrorCancelingSpecies(Molecule(smiles="CCCCCC"), (250, "kJ/mol"), lot, (260, "kJ/mol"))
        cls.hexene = ErrorCancelingSpecies(Molecule(smiles="C=CCCCC"), (275, "kJ/mol"), lot, (275, "kJ/mol"))
        cls.benzene = ErrorCancelingSpecies(Molecule(smiles="c1ccccc1"), (-50, "kJ/mol"), lot, (-80, "kJ/mol"))
        cls.caffeine = ErrorCancelingSpecies(Molecule(smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C"), (300, "kJ/mol"), lot)
        cls.ethyne = ErrorCancelingSpecies(Molecule(smiles="C#C"), (200, "kJ/mol"), lot)

    def test_creating_error_canceling_schemes(self):
        scheme = ErrorCancelingScheme(
            target=self.propene,
            reference_set=[self.butane, self.benzene, self.caffeine, self.ethyne],
            isodesmic_class="rc2",
            conserve_ring_size=True,
            limit_charges=True,
            limit_scope=True,
        )

        assert scheme.reference_species == [self.butane]

        isodesmic_scheme = IsodesmicScheme(self.propene, [self.butane, self.benzene, self.caffeine, self.ethyne])

        assert isodesmic_scheme.reference_species == [self.butane, self.benzene]

    def test_find_error_canceling_reaction(self):
        """
        Test that the MILP problem can be solved to find a single isodesmic reaction
        """
        scheme = IsodesmicScheme(
            self.propene,
            [self.propane, self.butane, self.butene, self.caffeine, self.ethyne],
        )

        # Note that caffeine and ethyne will not be allowed, so for the full set the indices are [0, 1, 2]
        rxn, _ = scheme._find_error_canceling_reaction([0, 1, 2])
        assert rxn.species[self.butane] == -1
        assert rxn.species[self.propane] == 1
        assert rxn.species[self.butene] == 1

    def test_multiple_error_canceling_reactions(self):
        """
        Test that multiple error canceling reactions can be found
        """
        scheme = IsodesmicScheme(
            self.propene,
            [
                self.propane,
                self.butane,
                self.butene,
                self.pentane,
                self.pentene,
                self.hexane,
                self.hexene,
                self.benzene,
            ],
        )

        reaction_list = scheme.multiple_error_canceling_reaction_search(n_reactions_max=20)
        assert len(reaction_list) == 6
        reaction_string = reaction_list.__repr__()
        # Consider both permutations of the products in the reaction string
        rxn_str1 = "<ErrorCancelingReaction 1*C=CC + 1*CCCC <=> 1*CCC + 1*C=CCC >"
        rxn_str2 = "<ErrorCancelingReaction 1*C=CC + 1*CCCC <=> 1*C=CCC + 1*CCC >"
        assert any(rxn_string in reaction_string for rxn_string in [rxn_str1, rxn_str2])

    def test_calculate_target_enthalpy(self):
        """
        Test that ErrorCancelingScheme is able to calculate thermochemistry for the target species
        """
        scheme = IsodesmicScheme(
            self.propene,
            [
                self.propane,
                self.butane,
                self.butene,
                self.pentane,
                self.pentene,
                self.hexane,
                self.hexene,
                self.benzene,
            ],
        )

        target_thermo, rxn_list = scheme.calculate_target_enthalpy(n_reactions_max=3)
        assert target_thermo.value_si == 110000.0
        assert isinstance(rxn_list[0], ErrorCancelingReaction)
