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

import os


import numpy as np

from rmgpy import settings
from rmgpy.chemkin import load_chemkin_file
from rmgpy.data.base import Entry, DatabaseError, ForbiddenStructures
from rmgpy.data.kinetics.common import (
    save_entry,
    find_degenerate_reactions,
    ensure_independent_atom_ids,
)
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species
import pytest


def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.load(
        path=os.path.join(settings["test_data.directory"], "testing_database"),
        thermo_libraries=["primaryThermoLibrary"],
        reaction_libraries=["GRI-Mech3.0"],
        kinetics_families=[
            "R_Recombination",
            "Disproportionation",
            "R_Addition_MultipleBond",
            "H_Abstraction",
            "intra_H_migration",
        ],
        testing=True,
        depository=False,
        solvation=False,
        surface=False,
    )
    # load empty forbidden structures to avoid any dependence on forbidden structures
    # for these tests
    for family in database.kinetics.families.values():
        family.forbidden = ForbiddenStructures()
    database.forbidden_structures = ForbiddenStructures()

    # Prepare the database by loading training reactions and averaging the rate rules
    for family in database.kinetics.families.values():
        if not family.auto_generated:
            family.add_rules_from_training(thermo_database=database.thermo)
            family.fill_rules_by_averaging_up(verbose=True)


def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg

    rmg.database = None


#####################################################


class TestKineticsDatabase:
    def test_load_families_incorrect(self):
        """Test invalid methods for loading kinetics families"""
        path = os.path.join(settings["test_data.directory"], "testing_database", "kinetics", "families")
        database = KineticsDatabase()
        database.load_recommended_families(os.path.join(path, "recommended.py"))

        with pytest.raises(DatabaseError):
            database.load_families(path, families="random")
        with pytest.raises(DatabaseError):
            database.load_families(path, families=["!H_Abstraction", "Disproportionation"])
        with pytest.raises(DatabaseError):
            database.load_families(path, families=["fake_family"])

    def test_load_families_correct(self):
        """Test valid methods for loading kinetics families."""
        path = os.path.join(settings["test_data.directory"], "testing_database", "kinetics", "families")
        database = KineticsDatabase()
        database.load_recommended_families(os.path.join(path, "recommended.py"))

        try:
            database.load_families(path, families=[])
        except DatabaseError:
            assert False, "Unable to load families using list []"

        try:
            database.load_families(path, families="none")
        except DatabaseError:
            assert False, "Unable to load families using keyword 'none'"

        try:
            database.load_families(path, families="default")
        except DatabaseError:
            assert False, "Unable to load families using keyword 'default'"

        try:
            database.load_families(path, families=["default", "pah"])
        except DatabaseError:
            assert False, "Unable to load families using list ['default', 'pah']"

        try:
            database.load_families(path, families=["R_Addition_MultipleBond"])
        except DatabaseError:
            assert False, "Unable to load families using list ['R_Addition_MultipleBond']"

        try:
            database.load_families(path, families=["!H_Abstraction", "!Disproportionation"])
        except DatabaseError:
            assert False, "Unable to load families using list ['!H_Abstraction', '!Disproportionation']"

        try:
            database.load_families(path, families="!pah")
        except DatabaseError:
            assert False, "Unable to load families using keyword '!pah'"

        try:
            database.load_families(path, families=["H_Abstraction", "pah"])
        except DatabaseError:
            assert False, "Unable to load families using list ['H_Abstraction', 'pah']"


class TestReactionDegeneracy:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        cls.database = database

    def assert_correct_reaction_degeneracy(
        self,
        reactants,
        expected_rxn_num,
        expected_degeneracy,
        family_label=None,
        products=None,
        adjlists=False,
    ):
        """
        Generates reactions for the provided species and checks the results
        against the expected values.

        Args:
            reactants: list of SMILES for the reacting species
            family_label: label of the reaction family to react in
            expected_rxn_num: number of independent reaction expected
            expected_degeneracy: set of expected degeneracy values
            products: list of SMILES for the desired products (optional)
            adjlists: bool indicating if the input format is adjacency lists (optional)
                      assumes that the input is SMILES if False or unspecified

        Returns:
            list of the generated reactions for further analysis if desired
        """
        method = Molecule.from_adjacency_list if adjlists else Molecule.from_smiles

        reactants = [method(Molecule(), identifier) for identifier in reactants]
        if products is not None:
            products = [method(Molecule(), identifier) for identifier in products]
        else:
            products = None

        families = [family_label] if family_label is not None else None

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, products, only_families=families)

        assert len(reaction_list) == expected_rxn_num, "Expected {0} reactions, not {1} for {2} in {3}.".format(
            expected_rxn_num, len(reaction_list), reactants, family_label
        )

        degeneracy = set([rxn.degeneracy for rxn in reaction_list])

        assert degeneracy == expected_degeneracy, "Expected degeneracies of {0}, not {1} for {2} in {3}.".format(
            expected_degeneracy, degeneracy, reactants, family_label
        )

        return reaction_list

    def test_r_addition_multiple_bond_benzene(self):
        """Test that the proper degeneracy is calculated for H addition to benzene"""
        family_label = "R_Addition_MultipleBond"
        reactants = ["c1ccccc1", "[H]"]

        correct_rxn_num = 1
        correct_degeneracy = {6}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

    def test_r_addition_multiple_bond_methyl_naphthalene(self):
        """Test that the proper degeneracy is calculated for H addition to methylnaphthalene"""
        family_label = "R_Addition_MultipleBond"
        reactants = ["C1=CC=C2C=CC=CC2=C1C", "[H]"]
        products = ["C[C]1CC=CC2=CC=CC=C12"]

        correct_rxn_num = 1
        correct_degeneracy = {1}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)

    def test_r_recombination_phenyl(self):
        """Test that the proper degeneracy is calculated for phenyl + H recombination"""
        family_label = "R_Recombination"
        reactants = ["[c]1ccccc1", "[H]"]

        correct_rxn_num = 1
        correct_degeneracy = {1}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

    def test_r_recombination_h(self):
        """Test that the proper degeneracy is calculated for H + H recombination"""
        family_label = "R_Recombination"
        reactants = ["[H]", "[H]"]

        correct_rxn_num = 1
        correct_degeneracy = {0.5}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

    def test_degeneracy_for_methyl_methyl_recombination(self):
        """Test that the proper degeneracy is calculated for methyl + methyl recombination"""

        family_label = "R_Recombination"
        reactants = [
            """
            multiplicity 2
            1 C u1 p0 c0 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 H u0 p0 c0 {1,S}
            """,
            """
            multiplicity 2
            1 C u1 p0 c0 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 H u0 p0 c0 {1,S}
            """,
        ]

        correct_rxn_num = 1
        correct_degeneracy = {0.5}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    def test_degeneracy_for_methyl_labeled_methyl_recombination(self):
        """Test that the proper degeneracy is calculated for methyl + labeled methyl recombination"""

        family_label = "R_Recombination"
        reactants = [
            """
            multiplicity 2
            1 C u1 p0 c0 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 H u0 p0 c0 {1,S}
            """,
            """
            multiplicity 2
            1 C u1 p0 c0 i13 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 H u0 p0 c0 {1,S}
            """,
        ]

        correct_rxn_num = 1
        correct_degeneracy = {1}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    def test_degeneracy_for_ethyl_ethyl_disproportionation(self):
        """Test that the proper degeneracy is calculated for ethyl + ethyl disproportionation"""

        family_label = "Disproportionation"
        reactants = [
            """
            multiplicity 2
            1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
            2 C u1 p0 c0 {1,S} {3,S} {4,S}
            3 H u0 p0 c0 {2,S}
            4 H u0 p0 c0 {2,S}
            5 H u0 p0 c0 {1,S}
            6 H u0 p0 c0 {1,S}
            7 H u0 p0 c0 {1,S}
            """,
            """
            multiplicity 2
            1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
            2 C u1 p0 c0 {1,S} {3,S} {4,S}
            3 H u0 p0 c0 {2,S}
            4 H u0 p0 c0 {2,S}
            5 H u0 p0 c0 {1,S}
            6 H u0 p0 c0 {1,S}
            7 H u0 p0 c0 {1,S}
            """,
        ]

        correct_rxn_num = 1
        correct_degeneracy = {3}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    def test_degeneracy_for_ethyl_labeled_ethyl_disproportionation(self):
        """Test that the proper degeneracy is calculated for ethyl + labeled ethyl disproportionation"""

        family_label = "Disproportionation"
        reactants = [
            """
            multiplicity 2
            1 C u0 p0 c0 i13 {2,S} {5,S} {6,S} {7,S}
            2 C u1 p0 c0 {1,S} {3,S} {4,S}
            3 H u0 p0 c0 {2,S}
            4 H u0 p0 c0 {2,S}
            5 H u0 p0 c0 {1,S}
            6 H u0 p0 c0 {1,S}
            7 H u0 p0 c0 {1,S}
            """,
            """
            multiplicity 2
            1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
            2 C u1 p0 c0 {1,S} {3,S} {4,S}
            3 H u0 p0 c0 {2,S}
            4 H u0 p0 c0 {2,S}
            5 H u0 p0 c0 {1,S}
            6 H u0 p0 c0 {1,S}
            7 H u0 p0 c0 {1,S}
            """,
        ]

        correct_rxn_num = 2
        correct_degeneracy = {3}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    @pytest.mark.skip(reason="WIP")
    def test_degeneracy_does_not_include_identical_atom_labels(self):
        """
        Test that rxns with identical atom ids are not counted twice for degeneracy

        Uses [H] + CC=C[CH]C -> H2 + [CH2]C=C[CH]C as an example. Since the reactant
        is symmetric, there should be a single reaction with a degeneracy of 6.

        Marked work_in_progress because the current multiple TS algorithm will
        differentiate the reactions based on template, resulting in 2 reactions
        each with a degeneracy of 6.
        """

        family_label = "H_Abstraction"
        reactants = ["[H]", "CC=C[CH]C"]
        products = ["[H][H]", "[CH2]C=C[CH]C"]

        correct_rxn_num = 1
        correct_degeneracy = {6}

        self.assert_correct_reaction_degeneracy(
            reactants,
            correct_rxn_num,
            correct_degeneracy,
            family_label,
            products=products,
        )

    def test_degeneracy_keeps_separate_transition_states_separated(self):
        """
        Test that rxns with multiple transition states are kept as separate reactions

        Uses C[C]=C + C=C[CH2] -> C=C=C + C=CC as an example. This reaction should have
        two transition states, which should occur regardless of reactant order.
        """
        family_label = "Disproportionation"
        reactants = ["C[C]=C", "C=C[CH2]"]
        products = ["C=C=C", "CC=C"]

        correct_rxn_num = 2
        correct_degeneracy = {1, 6}

        self.assert_correct_reaction_degeneracy(
            reactants,
            correct_rxn_num,
            correct_degeneracy,
            family_label,
            products=products,
        )

    def test_separate_transition_states_generated_regardless_of_reactant_order(self):
        """
        ensure rxns with multiple transition states are kept as separate reactions

        this test uses C[C]=C + C=C[CH2] -> C=C=C + C=CC as an example.
        This reaction should have two transition states, which should occur regardless
        of the order .
        """
        mol_a = Molecule().from_smiles("C=[C]C")
        mol_b = Molecule().from_smiles("C=C[CH2]")
        mol_c = Molecule().from_smiles("C=C=C")
        mol_d = Molecule().from_smiles("C=CC")

        family = database.kinetics.families["Disproportionation"]
        reaction_list = family._generate_reactions([mol_a, mol_b], products=[mol_c, mol_d])

        swapped_reaction_list = family._generate_reactions([mol_b, mol_a], products=[mol_c, mol_d])

        # eliminate rxns that do not match products
        templates = {}
        for rxn in reaction_list:
            try:
                templates[rxn.template[0]] += 1
            except KeyError:
                templates[rxn.template[0]] = 1
        reverse_templates = {}
        for rxn in swapped_reaction_list:
            try:
                reverse_templates[rxn.template[0]] += 1
            except KeyError:
                reverse_templates[rxn.template[0]] = 1

        assert reverse_templates == templates, "The reaction output did not output all the transition states in either order of reactants"

    def test_propyl_propyl_reaction_is_the_half_propyl_butyl(self):
        """
        test that propyl propyl r-recombination is the same rate as propyl butyl

        this test assures that r-recombination reactions from the same rate rule
        with identical reactants have half the reaction rate since there is a
        symmetrical transition state.
        """
        family_label = "R_Recombination"
        propyl = "CC[CH2]"
        butyl = "CCC[CH2]"

        rxn_list_pp = self.assert_correct_reaction_degeneracy([propyl, propyl], 1, {0.5}, family_label)
        rxn_list_pb = self.assert_correct_reaction_degeneracy([propyl, butyl], 1, {1}, family_label)

        family = self.database.kinetics.families[family_label]

        pp_reaction = rxn_list_pp[0]
        pb_reaction = rxn_list_pb[0]

        # get kinetics for each reaction
        pp_kinetics_list = family.get_kinetics(
            pp_reaction,
            pp_reaction.template,
            degeneracy=pp_reaction.degeneracy,
            estimator="rate rules",
        )
        assert (
            len(pp_kinetics_list) == 1
        ), "The propyl and propyl recombination should only return one reaction. It returned {0}. " "Here is the full kinetics: {1}".format(
            len(pp_kinetics_list), pp_kinetics_list
        )

        pb_kinetics_list = family.get_kinetics(
            pb_reaction,
            pb_reaction.template,
            degeneracy=pb_reaction.degeneracy,
            estimator="rate rules",
        )
        assert (
            len(pb_kinetics_list) == 1
        ), "The propyl and butyl recombination should only return one reaction. It returned {0}. " "Here is the full kinetics: {1}".format(
            len(pb_kinetics_list), pb_kinetics_list
        )

        # the same reaction group must be found or this test will not work
        assert pb_kinetics_list[0][0].comment in pp_kinetics_list[0][0].comment, (
            "This test found different kinetics for the two groups, so it will not function as expected\n"
            + str(pp_kinetics_list)
            + str(pb_kinetics_list)
        )

        # test that the kinetics are correct
        assert round(abs(pp_kinetics_list[0][0].get_rate_coefficient(300) * 2 - pb_kinetics_list[0][0].get_rate_coefficient(300)), 7) == 0

    def test_identical_reactants_have_similar_kinetics(self):
        """
        tests identical reactants have the same kinetics than different reactants.

        this test assures that r addition multiple bond reactions from the same
        rate rule have the same reaction rate if the reactants are identicaal
        since little changes in the reactant or transition state symmetry.

        This method should be more robust than just checking
        the degeneracy of reactions.
        """
        family_label = "R_Addition_MultipleBond"
        butenyl = "C=CC[CH2]"
        pentenyl = "C=CCC[CH2]"
        symmetric_product = ["[CH2]CC([CH2])CCC=C"]
        asymmetric_product = ["[CH2]CCC([CH2])CCC=C"]

        rxn_list_bb = self.assert_correct_reaction_degeneracy([butenyl, butenyl], 1, {1}, family_label, products=symmetric_product)
        rxn_list_bp = self.assert_correct_reaction_degeneracy([butenyl, pentenyl], 1, {1}, family_label, products=asymmetric_product)

        family = self.database.kinetics.families[family_label]

        bb_reaction = rxn_list_bb[0]
        bp_reaction = rxn_list_bp[0]

        bb_kinetics_list = family.get_kinetics(
            bb_reaction,
            bb_reaction.template,
            degeneracy=bb_reaction.degeneracy,
            estimator="rate rules",
        )
        assert (
            len(bb_kinetics_list) == 1
        ), "The butenyl and butenyl addition should only return one reaction. It returned {0}. " "Here is the full kinetics: {1}".format(
            len(bb_kinetics_list), bb_kinetics_list
        )

        bp_kinetics_list = family.get_kinetics(
            bp_reaction,
            bp_reaction.template,
            degeneracy=bp_reaction.degeneracy,
            estimator="rate rules",
        )
        assert (
            len(bp_kinetics_list) == 1
        ), "The butenyl and pentenyl addition should only return one reaction. It returned {0}. " "Here is the full kinetics: {1}".format(
            len(bp_kinetics_list), bp_kinetics_list
        )

        # the same reaction group must be found or this test will not work
        assert bp_kinetics_list[0][0].comment in bb_kinetics_list[0][0].comment, (
            "This test found different kinetics for the two groups, so it will not function as expected\n"
            + str(bb_kinetics_list)
            + str(bp_kinetics_list)
        )

        # test that the kinetics are correct
        assert round(abs(bb_kinetics_list[0][0].get_rate_coefficient(300) - bp_kinetics_list[0][0].get_rate_coefficient(300)), 7) == 0

    def test_reaction_degeneracy_independent_of_generatereactions_direction(self):
        """
        test_reaction_degeneracy_independent_of_generatereactions_direction

        Ensure the returned kinetics have the same degeneracy irrespective of
        whether _generate_reactions has forward = True or False
        """

        family = database.kinetics.families["Disproportionation"]

        mol_a = Molecule().from_smiles("C[CH2]")
        mol_b = Molecule().from_smiles("C[CH2]")
        mol_c = Molecule().from_smiles("C=C")
        mol_d = Molecule().from_smiles("CC")

        mol_a.assign_atom_ids()
        mol_b.assign_atom_ids()
        mol_c.assign_atom_ids()
        mol_d.assign_atom_ids()

        # generate reactions in both directions
        forward_reactions = family._generate_reactions([mol_a, mol_b], products=[mol_c, mol_d], forward=True)
        reverse_reactions = family._generate_reactions([mol_c, mol_d], products=[mol_a, mol_b], forward=False)

        forward_reactions = find_degenerate_reactions(forward_reactions)
        reverse_reactions = find_degenerate_reactions(reverse_reactions)

        assert (
            forward_reactions[0].degeneracy == reverse_reactions[0].degeneracy
        ), "the kinetics from forward and reverse directions had different degeneracies, {} and {} " "respectively".format(
            forward_reactions[0].degeneracy, reverse_reactions[0].degeneracy
        )

    def test_degeneracy_same_reactant_different_resonance_structure(self):
        """Test if degeneracy is correct when reacting different resonance structures."""
        family_label = "Disproportionation"
        reactants = ["CC=C[CH2]", "CC=C[CH2]"]
        products = ["CC=CC", "C=CC=C"]

        correct_rxn_num = 1
        correct_degeneracy = {3}

        reaction_list = self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)

        assert set(reaction_list[0].template) == {"C_rad/H2/Cd", "Cmethyl_Csrad/H/Cd"}

    def test_degeneracy_multiple_ts_different_template(self):
        """Test that reactions from different transition states are marked as duplicates."""
        family_label = "intra_H_migration"
        reactants = ["CCCC[CH]CCCCC"]
        products = ["[CH2]CCCCCCCCC"]

        correct_rxn_num = 2
        correct_degeneracy = {3}

        reaction_list = self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)

        assert reaction_list[0].duplicate
        assert reaction_list[1].duplicate

    def test_degeneracy_multiple_resonance_different_template(self):
        """Test that reactions from different resonance structures are not kept."""
        family_label = "H_Abstraction"
        reactants = ["c1ccccc1", "[CH3]"]

        correct_rxn_num = 1
        correct_degeneracy = {6}

        reaction_list = self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

        assert not reaction_list[0].duplicate

    def test_degeneracy_resonance_keep_isomorphic(self):
        """Test that we get the correct degeneracy for [CH2]C=C[CH2] + [H].

        Incorrect results would be obtained if isomorphic resonance structures are not kept.
        """
        family_label = "R_Recombination"
        reactants = ["[CH2]C=C[CH2]", "[OH]"]
        products = ["[CH2]C(O)C=C"]

        correct_rxn_num = 1
        correct_degeneracy = {2}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)


class TestKineticsCommentsParsing:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        cls.database = database

    def test_parse_kinetics(self):
        species, reactions = load_chemkin_file(
            os.path.join(settings["test_data.directory"], "parsing_data", "chem_annotated.inp"),
            os.path.join(
                settings["test_data.directory"],
                "parsing_data",
                "species_dictionary.txt",
            ),
        )

        sources = []
        for reaction in reactions:
            sources.append(self.database.kinetics.extract_source_from_comments(reaction))

        # Source 0 comes from a kinetics library
        assert "Library" in sources[0]
        assert sources[0]["Library"] == "GRI-Mech3.0"

        reconstructed_kinetics = self.database.kinetics.reconstruct_kinetics_from_source(reactions[0], sources[0], fix_barrier_height=True)
        A = reconstructed_kinetics.A.value_si
        n = reconstructed_kinetics.n.value_si
        assert round(abs(reactions[0].kinetics.A.value_si - A), 7) == 0
        assert round(abs(reactions[0].kinetics.n.value_si - n), 7) == 0

        # Source 1 comes from a single exact match to a rate rule
        assert "Rate Rules" in sources[1]
        assert sources[1]["Rate Rules"][0] == "Disproportionation"
        rules = sources[1]["Rate Rules"][1]["rules"]

        assert len(rules) == 1
        assert rules[0][0].label == "O_pri_rad;Cmethyl_Csrad"

        reconstructed_kinetics = self.database.kinetics.reconstruct_kinetics_from_source(reactions[1], sources[1], fix_barrier_height=True)
        A = reconstructed_kinetics.A.value_si
        n = reconstructed_kinetics.n.value_si
        assert round(abs(reactions[1].kinetics.A.value_si - A), 7) == 0
        assert round(abs(reactions[1].kinetics.n.value_si - n), 7) == 0

        # Source 2 comes from an averaged rate rule that even contains a rate rule from a training reaction
        assert "Rate Rules" in sources[2]
        assert sources[2]["Rate Rules"][0] == "Disproportionation"
        expected_rules = [
            "O2b;O_Csrad",
            "O_atom_triplet;O_Csrad",
            "CH2_triplet;O_Csrad",
            "O_pri_rad;O_Csrad",
            "O_rad/NonDeC;O_Csrad",
            "O_rad/NonDeO;O_Csrad",
            "Cd_pri_rad;O_Csrad",
            "CO_pri_rad;O_Csrad",
            "C_methyl;O_Csrad",
            "C_rad/H2/Cs;O_Csrad",
            "C_rad/H2/Cd;O_Csrad",
            "C_rad/H2/O;O_Csrad",
            "C_rad/H/NonDeC;O_Csrad",
            "C_rad/Cs3;O_Csrad",
            "H_rad;O_Csrad",
        ]

        rules = sources[2]["Rate Rules"][1]["rules"]
        training = sources[2]["Rate Rules"][1]["training"]

        actual_rule_labels = [rule.label for rule, weight in rules]

        assert len(rules) == len(expected_rules)
        for rule in expected_rules:
            assert rule in actual_rule_labels

        assert len(training) == 1
        assert training[0][1].index == 0  # Assert that the index of that training reaction is 1

        reconstructed_kinetics = self.database.kinetics.reconstruct_kinetics_from_source(reactions[2], sources[2], fix_barrier_height=True)
        A = reconstructed_kinetics.A.value_si
        n = reconstructed_kinetics.n.value_si
        A = round(A, -int(np.floor(np.log10(abs(A)))) + 3)  # Do some rounding since chemkin format kinetics are rounded
        n = round(n, 3)
        assert round(abs(reactions[2].kinetics.A.value_si - A), 7) == 0
        assert round(abs(reactions[2].kinetics.n.value_si - n), 7) == 0

        # Source 3 comes from a training reaction match
        assert "Training" in sources[3]
        family_label = sources[3]["Training"][0]
        training_rxn = sources[3]["Training"][1]

        assert family_label == "Disproportionation"
        assert training_rxn.label == "C2H + CH3O <=> C2H2 + CH2O"

        reconstructed_kinetics = self.database.kinetics.reconstruct_kinetics_from_source(reactions[3], sources[3], fix_barrier_height=True)
        A = reconstructed_kinetics.A.value_si
        n = reconstructed_kinetics.n.value_si
        assert round(abs(reactions[3].kinetics.A.value_si - A), 7) == 0
        assert round(abs(reactions[3].kinetics.n.value_si - n), 7) == 0

        # Source 3 comes from a pdep reaction
        assert "PDep" in sources[4]
        assert sources[4]["PDep"] == 7


class TestKinetics:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""

        global database
        cls.database = database

        cls.species, cls.reactions = load_chemkin_file(
            os.path.join(settings["test_data.directory"], "parsing_data", "chem_annotated.inp"),
            os.path.join(
                settings["test_data.directory"],
                "parsing_data",
                "species_dictionary.txt",
            ),
        )

    def test_react_molecules(self):
        """
        Test that reaction generation for Molecule objects works.
        """

        molecule_tuple = (Molecule(smiles="CC"), Molecule(smiles="[CH3]"))

        reaction_list = self.database.kinetics.react_molecules(molecule_tuple)

        assert reaction_list is not None
        assert all([isinstance(rxn, TemplateReaction) for rxn in reaction_list])

    def test_ensure_independent_atom_ids(self):
        """
        Ensure ensure_independent_atom_ids modifies atom labels
        """
        s1 = Species().from_smiles("CCC")
        s2 = Species().from_smiles("C=C[CH]C")
        assert s2.molecule[0].atoms[0].id == -1

        ensure_independent_atom_ids([s1, s2])
        # checks atom id
        assert s2.molecule[0].atoms[0].id != -1
        # checks second resonance structure id
        assert s2.molecule[1].atoms[0].id != -1

    def test_ensure_independent_atom_ids_no_resonance(self):
        """
        Ensure ensure_independent_atom_ids does not generate resonance
        """
        s1 = Species().from_smiles("CCC")
        s2 = Species().from_smiles("C=C[CH]C")
        assert s2.molecule[0].atoms[0].id == -1

        ensure_independent_atom_ids([s1, s2], resonance=False)
        # checks resonance structures
        assert len(s2.molecule) == 1
        # checks that atom ids are changed
        for atom in s2.molecule[0].atoms:
            assert atom.id != -1

    def test_save_entry(self):
        """
        tests that save entry can run
        """
        reactions = self.reactions

        fname = "testfile.txt"
        fid = open("testfile.txt", "w")

        wd = os.getcwd()
        wdir = wd + "/" + fname

        rxn = reactions[0]
        entry = Entry(
            index=1,
            label=str(rxn),
            item=rxn,
            short_desc="sdes",
            long_desc="lsdes",
            data="stuff",
            rank=0,
        )
        save_entry(fid, entry)

        fid.close()

        os.remove(wdir)

    def test_duplicates(self):
        """
        tests that kinetics libraries load properly and that
        the duplicate related routines run without error
        """
        lib = self.database.kinetics.libraries["GRI-Mech3.0"]
        out = lib.check_for_duplicates(True)
        assert out is None
        out = lib.convert_duplicates_to_multi()
        assert out is None

    def test_add_reverse_attribute(self):
        """
        tests that the add_reverse_attribute method gets the reverse degeneracy correct
        """
        from rmgpy.data.rmg import get_db
        from rmgpy.data.kinetics.family import TemplateReaction

        adjlist = [
            """
            multiplicity 2
            1 H u0 p0 c0 {7,S}
            2 H u0 p0 c0 {4,S}
            3 C u1 p0 c0 {5,S} {7,S} {8,S}
            4 C u0 p0 c0 {2,S} {6,S} {7,D}
            5 H u0 p0 c0 {3,S}
            6 H u0 p0 c0 {4,S}
            7 C u0 p0 c0 {1,S} {3,S} {4,D}
            8 H u0 p0 c0 {3,S}
            """,
            """
            1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
            2 C u0 p0 c0 i13 {1,S} {3,D} {7,S}
            3 C u0 p0 c0 {2,D} {8,S} {9,S}
            4 H u0 p0 c0 {1,S}
            5 H u0 p0 c0 {1,S}
            6 H u0 p0 c0 {1,S}
            7 H u0 p0 c0 {2,S}
            8 H u0 p0 c0 {3,S}
            9 H u0 p0 c0 {3,S}
            """,
            """
            multiplicity 2
            1 H u0 p0 c0 {7,S}
            2 H u0 p0 c0 {4,S}
            3 C u1 p0 c0 {5,S} {7,S} {8,S}
            4 C u0 p0 c0 {2,S} {6,S} {7,D}
            5 H u0 p0 c0 {3,S}
            6 H u0 p0 c0 {4,S}
            7 C u0 p0 c0 i13 {1,S} {3,S} {4,D}
            8 H u0 p0 c0 {3,S}
            """,
            """
            1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
            2 C u0 p0 c0 {1,S} {3,D} {7,S}
            3 C u0 p0 c0 {2,D} {8,S} {9,S}
            4 H u0 p0 c0 {1,S}
            5 H u0 p0 c0 {1,S}
            6 H u0 p0 c0 {1,S}
            7 H u0 p0 c0 {2,S}
            8 H u0 p0 c0 {3,S}
            9 H u0 p0 c0 {3,S}
            """,
        ]
        family = get_db("kinetics").families["H_Abstraction"]
        r1 = Species(molecule=[Molecule().from_adjacency_list(adjlist[0])])
        r2 = Species(molecule=[Molecule().from_adjacency_list(adjlist[1])])
        p1 = Species(molecule=[Molecule().from_adjacency_list(adjlist[2])])
        p2 = Species(molecule=[Molecule().from_adjacency_list(adjlist[3])])
        r1.generate_resonance_structures(keep_isomorphic=True)
        p1.generate_resonance_structures(keep_isomorphic=True)

        rxn = TemplateReaction(reactants=[r1, r2], products=[p1, p2])

        rxn.degeneracy = family.calculate_degeneracy(rxn)
        assert rxn.degeneracy == 6

        family.add_reverse_attribute(rxn)

        assert rxn.reverse.degeneracy == 6

    def test_calculate_degeneracy_for_non_reactive_molecule(self):
        """
        tests that the calculate_degeneracy method gets the degeneracy correct for unreactive molecules
        and that _generate_reactions work correctly with the react_non_reactive flag set to `True`.
        """
        from rmgpy.data.rmg import get_db
        from rmgpy.data.kinetics.family import TemplateReaction

        adjlist = [
            """
        multiplicity 2
        1 H u1 p0 c0""",
            """
        multiplicity 2
        1 O u1 p1 c+1 {2,D}
        2 N u0 p2 c-1 {1,D}""",
            """
        1 O u0 p1 c+1 {2,D} {3,S}
        2 N u0 p2 c-1 {1,D}
        3 H u0 p0 c0 {1,S}""",
        ]

        family = get_db("kinetics").families["R_Recombination"]
        r1 = Species(molecule=[Molecule().from_adjacency_list(adjlist[0])])
        r2 = Species(molecule=[Molecule().from_adjacency_list(adjlist[1])])  # r2 is not the representative structure of
        # NO, but it is the correct structure participating in this reaction
        p1 = Species(molecule=[Molecule().from_adjacency_list(adjlist[2])])
        r2.generate_resonance_structures(keep_isomorphic=True)

        rxn = TemplateReaction(reactants=[r1, r2], products=[p1])
        rxn.degeneracy = family.calculate_degeneracy(rxn)
        assert rxn.degeneracy == 1

    def test_generate_reactions_from_families_with_resonance(self):
        """Test that we can generate reactions from families with resonance structures"""
        reactants = [
            Molecule().from_smiles("CC=C[CH2]"),
            Molecule().from_smiles("[OH]"),
        ]
        expected_product_1 = Molecule().from_smiles("CC=CCO")
        expected_product_2 = Molecule().from_smiles("CC(O)C=C")

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, only_families=["R_Recombination"], resonance=True)

        assert len(reaction_list) == 2

        case_1 = reaction_list[0].products[0].is_isomorphic(expected_product_1) and reaction_list[1].products[0].is_isomorphic(expected_product_2)
        case_2 = reaction_list[0].products[0].is_isomorphic(expected_product_2) and reaction_list[1].products[0].is_isomorphic(expected_product_1)

        # Only one case should be true
        assert case_1 ^ case_2

    def test_generate_reactions_from_families_no_resonance(self):
        """Test that we can generate reactions from families without resonance structures"""
        reactants = [
            Molecule().from_smiles("CC=C[CH2]"),
            Molecule().from_smiles("[OH]"),
        ]
        expected_product = Molecule().from_smiles("CC=CCO")

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, only_families=["R_Recombination"], resonance=False)

        assert len(reaction_list) == 1

        assert reaction_list[0].products[0].is_isomorphic(expected_product)

    def test_generate_reactions_from_families_product_resonance(self):
        """Test that we can specify the product resonance structure when generating reactions"""
        reactants = [
            Molecule().from_smiles("CCC=C"),
            Molecule().from_smiles("[H]"),
        ]
        products = [
            Molecule().from_smiles("CC=C[CH2]"),
            Molecule().from_smiles("[H][H]"),
        ]

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, products, only_families=["H_Abstraction"], resonance=True)

        assert len(reaction_list) == 1
        assert reaction_list[0].degeneracy == 2

    def test_generate_reactions_from_families_product_resonance2(self):
        """Test that we can specify the no product resonance structure when generating reactions"""
        reactants = [
            Molecule().from_smiles("CCC=C"),
            Molecule().from_smiles("[H]"),
        ]
        products = [
            Molecule().from_smiles("CC=C[CH2]"),
            Molecule().from_smiles("[H][H]"),
        ]

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, products, only_families=["H_Abstraction"], resonance=False)
        assert len(reaction_list) == 0

    def test_generate_reactions_from_libraries(self):
        """Test that we can generate reactions from libraries"""
        reactants = [
            Molecule().from_smiles("CC=O"),
            Molecule().from_smiles("[H]"),
        ]

        reaction_list = self.database.kinetics.generate_reactions_from_libraries(reactants)

        assert len(reaction_list) == 3

    def test_generate_reactions_from_libraries2(self):
        """Test that we can generate reactions from libraries specifying products"""
        reactants = [
            Molecule().from_smiles("CC=O"),
            Molecule().from_smiles("[H]"),
        ]
        products = [
            Molecule().from_smiles("[CH2]C=O"),
            Molecule().from_smiles("[H][H]"),
        ]
        reaction_list_2 = self.database.kinetics.generate_reactions_from_libraries(reactants, products)

        assert len(reaction_list_2) == 1

    def test_add_atom_labels_for_reaction(self):
        """Test that add_atom_labels_for_reaction can identify reactions with resonance
        The molecule [CH]=C=C has resonance in this reaction"""
        from rmgpy.data.rmg import get_db

        reactants = [
            Molecule().from_smiles("C=C=C"),
            Molecule().from_smiles("[CH]=C=C"),
        ]
        products = [
            Molecule().from_smiles("C#C[CH2]"),
            Molecule().from_smiles("C#CC"),
        ]
        reaction = TemplateReaction(reactants=reactants, products=products, family="H_Abstraction")
        reaction.ensure_species(reactant_resonance=True, product_resonance=True)
        family = get_db("kinetics").families["H_Abstraction"]
        family.add_atom_labels_for_reaction(reaction, output_with_resonance=False)

        # test that the reaction has labels
        found_labels = []
        for species in reaction.reactants:
            for atom in species.molecule[0].atoms:
                if atom.label != "":
                    found_labels.append(atom.label)
        assert len(found_labels) == 3
        assert "*1" in found_labels
        assert "*2" in found_labels
        assert "*3" in found_labels

        # test for the products too
        found_labels = []
        for species in reaction.products:
            for atom in species.molecule[0].atoms:
                if atom.label != "":
                    found_labels.append(atom.label)
        assert len(found_labels) == 3
        assert "*1" in found_labels
        assert "*2" in found_labels
        assert "*3" in found_labels

    def test_add_atom_labels_for_reaction_2(self):
        """Test that add_atom_labels_for_reaction can identify reactions with identical references
        The molecule [CH]=C=C has resonance in this reaction"""
        from rmgpy.data.rmg import get_db

        s1 = Species().from_smiles("C=C=C")
        s2 = Species().from_smiles("C=C=[CH]")
        s3 = Species().from_smiles("C#CC")
        s2.generate_resonance_structures()
        reactants = [s1, s2]
        products = [s2, s3]
        reaction = TemplateReaction(reactants=reactants, products=products, family="H_Abstraction")
        family = get_db("kinetics").families["H_Abstraction"]
        print(reaction.reactants)
        print(reaction.products)
        family.add_atom_labels_for_reaction(reaction, output_with_resonance=False)

        # test that the reaction has labels
        found_labels = []
        for species in reaction.reactants:
            for atom in species.molecule[0].atoms:
                if atom.label != "":
                    found_labels.append(atom.label)
        assert len(found_labels) == 3, "wrong number of labels found {0}".format(found_labels)
        assert "*1" in found_labels
        assert "*2" in found_labels
        assert "*3" in found_labels

        # test for the products too
        found_labels = []
        for species in reaction.products:
            for atom in species.molecule[0].atoms:
                if atom.label != "":
                    found_labels.append(atom.label)
        assert len(found_labels) == 3
        assert "*1" in found_labels
        assert "*2" in found_labels
        assert "*3" in found_labels

    def test_add_atom_labels_for_reaction_3(self):
        """Test that add_atom_labels_for_reaction can identify reactions with resonance and isotopes"""
        from rmgpy.data.rmg import get_db

        mr0 = Molecule().from_adjacency_list(
            "1    C u0 p0 c0 i13 {3,D} {4,S} {5,S}\n2 *1 C u0 p0 c0 {3,D} {6,S} {7,S}\n3    C u0 p0 c0 {1,D} {2,D}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {2,S}\n7 *4 H u0 p0 c0 {2,S}\n"
        )
        mr1a = Molecule().from_adjacency_list(
            "multiplicity 2\n1    C u0 p0 c0 i13 {2,D} {4,S} {5,S}\n2    C u0 p0 c0 {1,D} {3,D}\n3 *1 C u1 p0 c0 {2,D} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n"
        )
        mr1b = Molecule().from_adjacency_list(
            "multiplicity 2\n1    C u1 p0 c0 i13 {2,S} {4,S} {5,S}\n2    C u0 p0 c0 {1,S} {3,T}\n3 *1 C u0 p0 c0 {2,T} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n"
        )
        mp1a = Molecule().from_adjacency_list(
            "multiplicity 2\n1    C u0 p0 c0 {2,D} {4,S} {5,S}\n2    C u0 p0 c0 {1,D} {3,D}\n3 *1 C u1 p0 c0 i13 {2,D} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n"
        )
        mp1b = Molecule().from_adjacency_list(
            "multiplicity 2\n1    C u1 p0 c0 {2,S} {4,S} {5,S}\n2    C u0 p0 c0 {1,S} {3,T}\n3 *1 C u0 p0 c0 i13 {2,T} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n"
        )
        s1 = Species(molecule=[mr0])
        s2 = Species(molecule=[mr1a, mr1b])
        s3 = Species(molecule=[mp1a, mp1b])
        reactants = [s1, s2]
        products = [s1, s3]
        reaction = TemplateReaction(reactants=reactants, products=products, family="H_Abstraction")
        family = get_db("kinetics").families["H_Abstraction"]
        print(reaction.reactants)
        print(reaction.products)
        family.add_atom_labels_for_reaction(reaction, output_with_resonance=False)

        # test that the reaction has labels
        found_labels = []
        for species in reaction.reactants:
            for atom in species.molecule[0].atoms:
                if atom.label != "":
                    found_labels.append(atom.label)
        assert len(found_labels) == 3, "wrong number of labels found {0}".format(found_labels)
        assert "*1" in found_labels
        assert "*2" in found_labels
        assert "*3" in found_labels

        # test for the products too
        found_labels = []
        for species in reaction.products:
            for atom in species.molecule[0].atoms:
                if atom.label != "":
                    found_labels.append(atom.label)
        assert len(found_labels) == 3
        assert "*1" in found_labels
        assert "*2" in found_labels
        assert "*3" in found_labels

    def test_species_preserved_after_generate_reactions(self):
        """
        Test that Species objects do not retain changes after generating reactions

        This tests a case involving identical reactants
        """
        reactant1 = Species(index=1, label="ethyl", smiles="C[CH2]")
        reactant1_copy = reactant1.copy(deep=True)  # These copies record the state of the original attributes
        expected_product_1 = Species(smiles="CC")
        expected_product_2 = Species(smiles="C=C")

        reaction_list = self.database.kinetics.generate_reactions_from_families(
            [reactant1, reactant1], only_families=["Disproportionation"], resonance=True
        )

        # First confirm that we get the expected reaction
        assert len(reaction_list) == 1
        reaction = reaction_list[0]
        case_1 = reaction.products[0].is_isomorphic(expected_product_1) and reaction.products[1].is_isomorphic(expected_product_2)
        case_2 = reaction.products[0].is_isomorphic(expected_product_2) and reaction.products[1].is_isomorphic(expected_product_1)
        # Only one case should be true
        assert case_1 or case_2

        reactant1_out, reactant2_out = reaction.reactants

        # The species in the output reaction should be new objects because we made them from Molecule objects
        assert reactant1 is not reactant1_out
        assert reactant1 is not reactant2_out

        # The molecule objects should be the same references for one output species
        assert reactant1.molecule[0] is not reactant1_out.molecule[0]
        assert reactant1.molecule[0] is reactant2_out.molecule[0]

        # They should be isomorphic
        assert reactant1.is_isomorphic(reactant1_out)
        assert reactant1.is_isomorphic(reactant2_out)
        assert reactant1_copy.is_isomorphic(reactant1_out)
        assert reactant1_copy.is_isomorphic(reactant2_out)

        # Now, we only care whether the original reactants have deviated from the copies
        # The output reactants will be replaced by the original reactants in CERM.check_for_existing_species
        assert reactant1.index == reactant1_copy.index
        assert reactant1.label == reactant1_copy.label
        assert reactant1.props == reactant1_copy.props
        assert reactant1.molecule == reactant1_copy.molecule
        assert reactant1.molecule[0].get_all_labeled_atoms() == reactant1_copy.molecule[0].get_all_labeled_atoms()
        assert reactant1.molecule[0].props == reactant1_copy.molecule[0].props

    def test_species_preserved_after_generate_reactions_2(self):
        """
        Test that Species objects do not retain changes after generating reactions

        This tests a case involving benzene bond modification
        """
        reactant1 = Species(index=1, label="methyl", smiles="[CH3]")
        reactant2 = Species(index=2, label="benzene", smiles="c1ccccc1")
        reactant2.generate_resonance_structures()  # Only benzene has resonance structures
        reactant1_copy = reactant1.copy(deep=True)  # These copies record the state of the original attributes
        reactant2_copy = reactant2.copy(deep=True)
        expected_product = Species(smiles="CC1[CH]C=CC=C1")

        reaction_list = self.database.kinetics.generate_reactions_from_families(
            [reactant1, reactant2],
            only_families=["R_Addition_MultipleBond"],
            resonance=True,
        )

        # First confirm that we get the expected reaction
        assert len(reaction_list) == 1
        reaction = reaction_list[0]
        assert reaction.products[0].is_isomorphic(expected_product)

        reactant1_out, reactant2_out = reaction.reactants

        # The species in the output reaction should be new objects because we made them from Molecule objects
        assert reactant1 is not reactant1_out
        assert reactant2 is not reactant2_out

        # However, the molecule objects should be the same references
        assert reactant1.molecule[0] is reactant1_out.molecule[0]
        assert reactant2.molecule[0] is reactant2_out.molecule[0]

        # They should be isomorphic
        assert reactant1.is_isomorphic(reactant1_out)
        assert reactant2.is_isomorphic(reactant2_out)
        assert reactant1_copy.is_isomorphic(reactant1_out)
        assert reactant2_copy.is_isomorphic(reactant2_out)

        # Now, we only care whether the original reactants have deviated from the copies
        # The output reactants will be replaced by the original reactants in CERM.check_for_existing_species
        assert reactant1.index == reactant1_copy.index
        assert reactant2.index == reactant2_copy.index
        assert reactant1.label == reactant1_copy.label
        assert reactant2.label == reactant2_copy.label
        assert reactant1.props == reactant1_copy.props
        assert reactant2.props == reactant2_copy.props
        assert reactant1.molecule == reactant1_copy.molecule
        assert reactant2.molecule == reactant2_copy.molecule
        assert reactant1.molecule[0].get_all_labeled_atoms() == reactant1_copy.molecule[0].get_all_labeled_atoms()
        assert reactant2.molecule[0].get_all_labeled_atoms() == reactant2_copy.molecule[0].get_all_labeled_atoms()
        assert reactant1.molecule[0].props == reactant1_copy.molecule[0].props
        assert reactant2.molecule[0].props == reactant2_copy.molecule[0].props
