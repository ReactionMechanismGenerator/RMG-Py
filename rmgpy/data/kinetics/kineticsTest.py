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

import os
import unittest
from external.wip import work_in_progress

import numpy

from rmgpy import settings
from rmgpy.chemkin import loadChemkinFile
from rmgpy.data.base import Entry, DatabaseError, ForbiddenStructures
from rmgpy.data.kinetics.common import saveEntry, filter_reactions, find_degenerate_reactions, ensure_independent_atom_ids
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species


###################################################

def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.load(
        path=os.path.join(settings['test_data.directory'], 'testing_database'),
        thermoLibraries=['primaryThermoLibrary'],
        reactionLibraries=['GRI-Mech3.0'],
        kineticsFamilies=[
            'R_Recombination',
            'Disproportionation',
            'R_Addition_MultipleBond',
            'H_Abstraction',
            'intra_H_migration',
        ],
        testing=True,
        depository=False,
        solvation=False,
    )
    #load empty forbidden structures to avoid any dependence on forbidden structures
    #for these tests
    for family in database.kinetics.families.values():
        family.forbidden = ForbiddenStructures()
    database.forbiddenStructures = ForbiddenStructures()

    # Prepare the database by loading training reactions and averaging the rate rules
    for family in database.kinetics.families.values():
        family.addKineticsRulesFromTrainingSet(thermoDatabase=database.thermo)
        family.fillKineticsRulesByAveragingUp(verbose=True)


def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg
    rmg.database = None

#####################################################


class TestKineticsDatabase(unittest.TestCase):

    def testLoadFamilies(self):
        """
        Test that the loadFamilies function raises the correct exceptions
        """
        path = os.path.join(settings['test_data.directory'],'testing_database','kinetics','families')
        database = KineticsDatabase()

        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families='random')
        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families=['!H_Abstraction','Disproportionation'])
        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families=['fake_family'])
        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families=[])

class TestReactionDegeneracy(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database

    def assert_correct_reaction_degeneracy(self, reactants, expected_rxn_num, expected_degeneracy,
                                           family_label=None, products=None, adjlists=False):
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
        method = Molecule.fromAdjacencyList if adjlists else Molecule.fromSMILES

        reactants = [method(Molecule(), identifier) for identifier in reactants]
        if products is not None:
            products = [method(Molecule(), identifier) for identifier in products]
        else:
            products = None

        families = [family_label] if family_label is not None else None

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, products,
                                                                                only_families=families)

        self.assertEqual(len(reaction_list), expected_rxn_num,
                         'Expected {0} reactions, not {1} for {2} in {3}.'.format(expected_rxn_num,
                                                                                  len(reaction_list),
                                                                                  reactants,
                                                                                  family_label))

        degeneracy = set([rxn.degeneracy for rxn in reaction_list])

        self.assertEqual(degeneracy, expected_degeneracy,
                         'Expected degeneracies of {0}, not {1} for {2} in {3}.'.format(expected_degeneracy,
                                                                                        degeneracy,
                                                                                        reactants,
                                                                                        family_label))

        return reaction_list

    def testR_Addition_MultipleBondBenzene(self):
        """Test that the proper degeneracy is calculated for H addition to benzene"""
        family_label = 'R_Addition_MultipleBond'
        reactants = ['c1ccccc1', '[H]']

        correct_rxn_num = 1
        correct_degeneracy = {6}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

    def testR_Addition_MultipleBondMethylNaphthalene(self):
        """Test that the proper degeneracy is calculated for H addition to methylnaphthalene"""
        family_label = 'R_Addition_MultipleBond'
        reactants = ['C1=CC=C2C=CC=CC2=C1C', '[H]']
        products = ['C[C]1CC=CC2=CC=CC=C12']

        correct_rxn_num = 1
        correct_degeneracy = {1}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)

    def testR_RecombinationPhenyl(self):
        """Test that the proper degeneracy is calculated for phenyl + H recombination"""
        family_label = 'R_Recombination'
        reactants = ['[c]1ccccc1', '[H]']

        correct_rxn_num = 1
        correct_degeneracy = {1}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

    def testR_RecombinationH(self):
        """Test that the proper degeneracy is calculated for H + H recombination"""
        family_label = 'R_Recombination'
        reactants = ['[H]', '[H]']

        correct_rxn_num = 1
        correct_degeneracy = {0.5}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

    def test_degeneracy_for_methyl_methyl_recombination(self):
        """Test that the proper degeneracy is calculated for methyl + methyl recombination"""

        family_label = 'R_Recombination'
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
            """
        ]

        correct_rxn_num = 1
        correct_degeneracy = {0.5}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    def test_degeneracy_for_methyl_labeled_methyl_recombination(self):
        """Test that the proper degeneracy is calculated for methyl + labeled methyl recombination"""

        family_label = 'R_Recombination'
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
            """
        ]

        correct_rxn_num = 1
        correct_degeneracy = {1}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    def test_degeneracy_for_ethyl_ethyl_disproportionation(self):
        """Test that the proper degeneracy is calculated for ethyl + ethyl disproportionation"""

        family_label = 'Disproportionation'
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
            """
        ]

        correct_rxn_num = 1
        correct_degeneracy = {3}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    def test_degeneracy_for_ethyl_labeled_ethyl_disproportionation(self):
        """Test that the proper degeneracy is calculated for ethyl + labeled ethyl disproportionation"""

        family_label = 'Disproportionation'
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
            """
        ]

        correct_rxn_num = 2
        correct_degeneracy = {3}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, adjlists=True)

    @work_in_progress
    def test_degeneracy_does_not_include_identical_atom_labels(self):
        """
        Test that rxns with identical atom ids are not counted twice for degeneracy

        Uses [H] + CC=C[CH]C -> H2 + [CH2]C=C[CH]C as an example. Since the reactant
        is symmetric, there should be a single reaction with a degeneracy of 6.

        Marked work_in_progress because the current multiple TS algorithm will
        differentiate the reactions based on template, resulting in 2 reactions
        each with a degeneracy of 6.
        """

        family_label = 'H_Abstraction'
        reactants = ['[H]', 'CC=C[CH]C']
        products = ['[H][H]', '[CH2]C=C[CH]C']

        correct_rxn_num = 1
        correct_degeneracy = {6}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products=products)

    def test_degeneracy_keeps_separate_transition_states_separated(self):
        """
        Test that rxns with multiple transition states are kept as separate reactions
        
        Uses C[C]=C + C=C[CH2] -> C=C=C + C=CC as an example. This reaction should have
        two transition states, which should occur regardless of reactant order.
        """
        family_label = 'Disproportionation'
        reactants = ['C[C]=C', 'C=C[CH2]']
        products = ['C=C=C', 'CC=C']

        correct_rxn_num = 2
        correct_degeneracy = {1, 6}

        self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products=products)

    def test_separate_transition_states_generated_regardless_of_reactant_order(self):
        """
        ensure rxns with multiple transition states are kept as separate reactions
        
        this test uses C[C]=C + C=C[CH2] -> C=C=C + C=CC as an example. 
        This reaction should have two transition states, which should occur regardless
        of the order .
        """
        molA = Molecule().fromSMILES('C=[C]C')
        molB = Molecule().fromSMILES('C=C[CH2]')
        molC = Molecule().fromSMILES('C=C=C')
        molD = Molecule().fromSMILES('C=CC')
        reactionList = database.kinetics.families['Disproportionation']._KineticsFamily__generateReactions([molA, molB], products=[molC,molD])
        
        swapped_reactionList = database.kinetics.families['Disproportionation']._KineticsFamily__generateReactions([molB, molA], products=[molC,molD])
        
        
        # eliminate rxns that do not match products
        templates = {}
        for rxn in reactionList:
            try:
                templates[rxn.template[0]] += 1
            except KeyError:
                templates[rxn.template[0]] = 1
        reverseTemplates = {}
        for rxn in swapped_reactionList:  
            try:
                reverseTemplates[rxn.template[0]] += 1
            except KeyError:
                reverseTemplates[rxn.template[0]] = 1

        self.assertEqual(reverseTemplates, templates,'The reaction output did not output all the transition states in either order of reactants')

    def test_propyl_propyl_reaction_is_the_half_propyl_butyl(self):
        """
        test that propyl propyl r-recombination is the same rate as propyl butyl

        this test assures that r-recombination reactions from the same rate rule
        with identical reactants have half the reaction rate since there is a 
        symmetrical transition state.
        """
        family_label = 'R_Recombination'
        propyl = 'CC[CH2]'
        butyl = 'CCC[CH2]'

        rxn_list_pp = self.assert_correct_reaction_degeneracy([propyl, propyl], 1, {0.5}, family_label)
        rxn_list_pb = self.assert_correct_reaction_degeneracy([propyl, butyl], 1, {1}, family_label)

        family = self.database.kinetics.families[family_label]

        pp_reaction = rxn_list_pp[0]
        pb_reaction = rxn_list_pb[0]

        # get kinetics for each reaction
        pp_kinetics_list = family.getKinetics(pp_reaction, pp_reaction.template,
                                              degeneracy=pp_reaction.degeneracy,
                                              estimator='rate rules')
        self.assertEqual(len(pp_kinetics_list), 1,
                         'The propyl and propyl recombination should only return one reaction. \
                         It returned {0}. Here is the full kinetics: {1}'.format(len(pp_kinetics_list), pp_kinetics_list))

        pb_kinetics_list = family.getKinetics(pb_reaction, pb_reaction.template,
                                              degeneracy=pb_reaction.degeneracy,
                                              estimator='rate rules')
        self.assertEqual(len(pb_kinetics_list), 1,
                         'The propyl and butyl recombination should only return one reaction. \
                         It returned {0}. Here is the full kinetics: {1}'.format(len(pb_kinetics_list), pb_kinetics_list))

        # the same reaction group must be found or this test will not work
        self.assertIn(pb_kinetics_list[0][0].comment, pp_kinetics_list[0][0].comment,
                      'This test found different kinetics for the two groups, so it will not function as expected\n' +
                      str(pp_kinetics_list)+str(pb_kinetics_list))

        # test that the kinetics are correct
        self.assertAlmostEqual(pp_kinetics_list[0][0].getRateCoefficient(300) * 2, pb_kinetics_list[0][0].getRateCoefficient(300))

    def test_identical_reactants_have_similar_kinetics(self):
        """
        tests identical reactants have the same kinetics than different reactants.
        
        this test assures that r addition multiple bond reactions from the same 
        rate rule have the same reaction rate if the reactants are identicaal 
        since little changes in the reactant or transition state symmetry. 
        
        This method should be more robust than just checking
        the degeneracy of reactions.
        """
        family_label = 'R_Addition_MultipleBond'
        butenyl = 'C=CC[CH2]'
        pentenyl = 'C=CCC[CH2]'
        symmetric_product = ['[CH2]CC([CH2])CCC=C']
        asymmetric_product = ['[CH2]CCC([CH2])CCC=C']

        rxn_list_bb = self.assert_correct_reaction_degeneracy([butenyl, butenyl], 1, {1}, family_label, products=symmetric_product)
        rxn_list_bp = self.assert_correct_reaction_degeneracy([butenyl, pentenyl], 1, {1}, family_label, products=asymmetric_product)

        family = self.database.kinetics.families[family_label]

        bb_reaction = rxn_list_bb[0]
        bp_reaction = rxn_list_bp[0]

        bb_kinetics_list = family.getKinetics(bb_reaction, bb_reaction.template,
                                              degeneracy=bb_reaction.degeneracy,
                                              estimator='rate rules')
        self.assertEqual(len(bb_kinetics_list), 1,
                         'The butenyl and butenyl addition should only return one reaction. \
                         It returned {0}. Here is the full kinetics: {1}'.format(len(bb_kinetics_list), bb_kinetics_list))

        bp_kinetics_list = family.getKinetics(bp_reaction, bp_reaction.template,
                                              degeneracy=bp_reaction.degeneracy,
                                              estimator='rate rules')
        self.assertEqual(len(bp_kinetics_list), 1,
                         'The butenyl and pentenyl addition should only return one reaction. \
                         It returned {0}. Here is the full kinetics: {1}'.format(len(bp_kinetics_list), bp_kinetics_list))

        # the same reaction group must be found or this test will not work
        self.assertIn(bp_kinetics_list[0][0].comment, bb_kinetics_list[0][0].comment,
                      'This test found different kinetics for the two groups, so it will not function as expected\n' +
                      str(bb_kinetics_list)+str(bp_kinetics_list))

        # test that the kinetics are correct
        self.assertAlmostEqual(bb_kinetics_list[0][0].getRateCoefficient(300), bp_kinetics_list[0][0].getRateCoefficient(300))
        
    def test_reaction_degeneracy_independent_of_generatereactions_direction(self):
        """
        test_reaction_degeneracy_independent_of_generatereactions_direction
        
        Ensure the returned kinetics have the same degeneracy irrespective of
        whether __generateReactions has forward = True or False
        """

        family = database.kinetics.families['Disproportionation']

        molA = Molecule().fromSMILES('C[CH2]')
        molB = Molecule().fromSMILES('C[CH2]')
        molC = Molecule().fromSMILES('C=C')
        molD = Molecule().fromSMILES('CC')
        
        molA.assignAtomIDs()
        molB.assignAtomIDs()
        molC.assignAtomIDs()
        molD.assignAtomIDs()

        # generate reactions in both directions
        forward_reactions = family._KineticsFamily__generateReactions([molA, molB], products=[molC, molD], forward=True)
        reverse_reactions = family._KineticsFamily__generateReactions([molC, molD], products=[molA, molB], forward=False)

        forward_reactions = find_degenerate_reactions(forward_reactions)
        reverse_reactions = find_degenerate_reactions(reverse_reactions)

        self.assertEqual(forward_reactions[0].degeneracy, reverse_reactions[0].degeneracy,
                         'the kinetics from forward and reverse directions had different degeneracies, {} and {} respectively'.format(forward_reactions[0].degeneracy, reverse_reactions[0].degeneracy))

    def test_degeneracy_same_reactant_different_resonance_structure(self):
        """Test if degeneracy is correct when reacting different resonance structures."""
        family_label = 'Disproportionation'
        reactants = ['CC=C[CH2]', 'CC=C[CH2]']
        products = ['CC=CC', 'C=CC=C']

        correct_rxn_num = 1
        correct_degeneracy = {3}

        reaction_list = self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)

        self.assertEqual(set(reaction_list[0].template), {'C_rad/H2/Cd', 'Cmethyl_Csrad/H/Cd'})

    def test_degeneracy_multiple_ts_different_template(self):
        """Test that reactions from different transition states are marked as duplicates."""
        family_label = 'intra_H_migration'
        reactants = ['CCCC[CH]CCCCC']
        products = ['[CH2]CCCCCCCCC']

        correct_rxn_num = 2
        correct_degeneracy = {3}

        reaction_list = self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)

        self.assertTrue(reaction_list[0].duplicate)
        self.assertTrue(reaction_list[1].duplicate)

    def test_degeneracy_multiple_resonance_different_template(self):
        """Test that reactions from different resonance structures are not kept."""
        family_label = 'H_Abstraction'
        reactants = ['c1ccccc1', '[CH3]']

        correct_rxn_num = 1
        correct_degeneracy = {6}

        reaction_list = self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label)

        self.assertFalse(reaction_list[0].duplicate)


class TestKineticsCommentsParsing(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database

    def testParseKinetics(self):
        species, reactions = loadChemkinFile(os.path.join(settings['test_data.directory'], 'parsing_data','chem_annotated.inp'),
                                             os.path.join(settings['test_data.directory'], 'parsing_data','species_dictionary.txt')
                                                       )
        
        sources = []
        for reaction in reactions:
            sources.append(self.database.kinetics.extractSourceFromComments(reaction))
              
        # Source 0 comes from a kinetics library
        self.assertTrue('Library' in sources[0])
        self.assertEqual(sources[0]['Library'], 'GRI-Mech3.0')
        
        reconstructedKinetics = self.database.kinetics.reconstructKineticsFromSource(reactions[0],sources[0],fixBarrierHeight=True)
        A = reconstructedKinetics.A.value_si
        n = reconstructedKinetics.n.value_si        
        self.assertAlmostEqual(reactions[0].kinetics.A.value_si,A)
        self.assertAlmostEqual(reactions[0].kinetics.n.value_si,n)
        
        
        # Source 1 comes from a single exact match to a rate rule
        self.assertTrue('Rate Rules' in sources[1])
        self.assertEqual(sources[1]['Rate Rules'][0],'Disproportionation')
        rules = sources[1]['Rate Rules'][1]['rules']

        self.assertEqual(len(rules),1)
        self.assertEqual(rules[0][0].label, 'O_pri_rad;Cmethyl_Csrad')
        
        
        reconstructedKinetics = self.database.kinetics.reconstructKineticsFromSource(reactions[1],sources[1],fixBarrierHeight=True)
        A = reconstructedKinetics.A.value_si
        n = reconstructedKinetics.n.value_si
        self.assertAlmostEqual(reactions[1].kinetics.A.value_si,A)
        self.assertAlmostEqual(reactions[1].kinetics.n.value_si,n)
        
        # Source 2 comes from an averaged rate rule that even contains a rate rule from a training reaction
        self.assertTrue('Rate Rules' in sources[2])
        self.assertEqual(sources[2]['Rate Rules'][0],'Disproportionation')
        expectedRules = ['O2b;O_Csrad', 'O_atom_triplet;O_Csrad', 'CH2_triplet;O_Csrad', 'O_pri_rad;O_Csrad', 
                        'O_rad/NonDeC;O_Csrad','O_rad/NonDeO;O_Csrad', 'Cd_pri_rad;O_Csrad', 'CO_pri_rad;O_Csrad','C_methyl;O_Csrad','C_rad/H2/Cs;O_Csrad','C_rad/H2/Cd;O_Csrad',
                        'C_rad/H2/O;O_Csrad','C_rad/H/NonDeC;O_Csrad','C_rad/Cs3;O_Csrad','H_rad;O_Csrad']
        
        rules = sources[2]['Rate Rules'][1]['rules']
        training = sources[2]['Rate Rules'][1]['training']
        

        actualRuleLabels = [rule.label for rule, weight in rules]
        
        self.assertEqual(len(rules),len(expectedRules))
        for rule in expectedRules:
            self.assertTrue(rule in actualRuleLabels)
        
                    
        self.assertEqual(len(training),1)
        self.assertEqual(training[0][1].index,0)  # Assert that the index of that training reaction is 1
        
        reconstructedKinetics = self.database.kinetics.reconstructKineticsFromSource(reactions[2],sources[2],fixBarrierHeight=True)
        A = reconstructedKinetics.A.value_si
        n = reconstructedKinetics.n.value_si
        A = round(A, -int(numpy.floor(numpy.log10(abs(A))))+3)  # Do some rounding since chemkin format kinetics are rounded
        n = round(n,3)
        self.assertAlmostEqual(reactions[2].kinetics.A.value_si,A)
        self.assertAlmostEqual(reactions[2].kinetics.n.value_si,n)

        # Source 3 comes from a training reaction match
        self.assertTrue('Training' in sources[3])
        familyLabel = sources[3]['Training'][0]
        trainingRxn = sources[3]['Training'][1]
        
        self.assertEqual(familyLabel,'Disproportionation')
        self.assertEqual(trainingRxn.label, 'C2H + CH3O <=> C2H2 + CH2O')
        
        reconstructedKinetics = self.database.kinetics.reconstructKineticsFromSource(reactions[3],sources[3], fixBarrierHeight=True)
        A = reconstructedKinetics.A.value_si
        n = reconstructedKinetics.n.value_si
        self.assertAlmostEqual(reactions[3].kinetics.A.value_si,A)
        self.assertAlmostEqual(reactions[3].kinetics.n.value_si,n)
        
        # Source 3 comes from a pdep reaction        
        self.assertTrue('PDep' in sources[4])
        self.assertEqual(sources[4]['PDep'], 7)

class TestKinetics(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""

        global database
        self.database = database
        
        self.species, self.reactions = loadChemkinFile(
            os.path.join(settings['test_data.directory'], 'parsing_data', 'chem_annotated.inp'),
            os.path.join(settings['test_data.directory'], 'parsing_data', 'species_dictionary.txt')
        )
        
    def test_filter_reactions(self):
        """
        tests that filter reactions removes reactions that are missing
        any reactants or products
        """
        
        reactions=self.reactions
        
        reactants = []
        products = []
        for x in reactions:
            reactants += x.reactants
            products += x.products
        
        lrset = set(reactants[6:])
        mlrset = {reactants[i].molecule[0] for i in range(6,len(reactants))}
        
        reactants = set(reactants)
        products = set(products)
        mreactants = {i.molecule[0] for i in reactants}
        mproducts = {i.molecule[0] for i in products}

        newmreactants = list(mreactants-mlrset)
        newmproducts = list(mproducts-mlrset)

        out = filter_reactions(newmreactants, newmproducts, reactions)
            
        rset = list(set(reactions) - set(out))

        msets = [set(i.reactants+i.products) for i in rset]
        
        for i, iset in enumerate(msets): #test that all the reactions we removed are missing a reactant or a product
            self.assertTrue(iset & lrset != set(),msg='reaction {0} removed improperly'.format(rset[i]))
        
        outsets = [set(i.reactants+i.products) for i in out]
            
        for i, iset in enumerate(outsets): #test that all the reactions left in aren't missing any reactants or products
            self.assertTrue(iset & lrset == set(),msg='reaction {0} left in improperly, should have removed in based on presence of {1}'.format(out[i],iset & lrset))

    def test_react_molecules(self):
        """
        Test that reaction generation for Molecule objects works.
        """

        moleculeTuple = (Molecule(SMILES='CC'), Molecule(SMILES='[CH3]'))

        reactionList = self.database.kinetics.react_molecules(moleculeTuple)

        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def test_ensure_independent_atom_ids(self):
        """
        Ensure ensure_independent_atom_ids modifies atomlabels
        """
        s1 = Species().fromSMILES('CCC')
        s2 = Species().fromSMILES('C=C[CH]C')
        self.assertEqual(s2.molecule[0].atoms[0].id, -1)

        ensure_independent_atom_ids([s1, s2])
        # checks atom id
        self.assertNotEqual(s2.molecule[0].atoms[0].id, -1)
        # checks second resonance structure id
        self.assertNotEqual(s2.molecule[1].atoms[0].id, -1)

    def test_ensure_independent_atom_ids_no_resonance(self):
        """
        Ensure ensure_independent_atom_ids does not generate resonance
        """
        s1 = Species().fromSMILES('CCC')
        s2 = Species().fromSMILES('C=C[CH]C')
        self.assertEqual(s2.molecule[0].atoms[0].id, -1)

        ensure_independent_atom_ids([s1, s2],resonance=False)
        # checks resonance structures
        self.assertEqual(len(s2.molecule),1)
        # checks that atom ids are changed
        for atom in s2.molecule[0].atoms:
            self.assertNotEqual(atom.id, -1)

    def testSaveEntry(self):
        """
        tests that save entry can run
        """
        reactions=self.reactions
        
        fname  = 'testfile.txt'
        fid = open('testfile.txt','w')
        
        wd = os.getcwd()
        wdir = wd+'/'+fname
        
        rxn = reactions[0]
        entry = Entry(index=1,label=str(rxn),item=rxn,shortDesc='sdes',longDesc='lsdes',data='stuff',rank=0)
        saveEntry(fid,entry)
        
        fid.close()
        
        os.remove(wdir)
        
    def testDuplicates(self):
        """
        tests that kinetics libraries load properly and that
        the duplicate related routines run without error
        """
        lib = self.database.kinetics.libraries['GRI-Mech3.0']
        out = lib.checkForDuplicates(True)
        self.assertIsNone(out)
        out = lib.convertDuplicatesToMulti()
        self.assertIsNone(out)

    def testaddReverseAttribute(self):
        """
        tests that the addReverseAttribute method gets the reverse degeneracy correct
        """
        from rmgpy.data.rmg import getDB
        from rmgpy.data.kinetics.family import TemplateReaction
        adjlist = ['''
        multiplicity 2
        1 H u0 p0 c0 {7,S}
        2 H u0 p0 c0 {4,S}
        3 C u1 p0 c0 {5,S} {7,S} {8,S}
        4 C u0 p0 c0 {2,S} {6,S} {7,D}
        5 H u0 p0 c0 {3,S}
        6 H u0 p0 c0 {4,S}
        7 C u0 p0 c0 {1,S} {3,S} {4,D}
        8 H u0 p0 c0 {3,S}
        ''',
          '''
        1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
        2 C u0 p0 c0 i13 {1,S} {3,D} {7,S}
        3 C u0 p0 c0 {2,D} {8,S} {9,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {1,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {3,S}
        9 H u0 p0 c0 {3,S}
        ''',
                '''
        multiplicity 2
        1 H u0 p0 c0 {7,S}
        2 H u0 p0 c0 {4,S}
        3 C u1 p0 c0 {5,S} {7,S} {8,S}
        4 C u0 p0 c0 {2,S} {6,S} {7,D}
        5 H u0 p0 c0 {3,S}
        6 H u0 p0 c0 {4,S}
        7 C u0 p0 c0 i13 {1,S} {3,S} {4,D}
        8 H u0 p0 c0 {3,S}
        ''',
          '''
        1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
        2 C u0 p0 c0 {1,S} {3,D} {7,S}
        3 C u0 p0 c0 {2,D} {8,S} {9,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {1,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {3,S}
        9 H u0 p0 c0 {3,S}
        '''
          ]
        family = getDB('kinetics').families['H_Abstraction']
        r1 = Species(molecule=[Molecule().fromAdjacencyList(adjlist[0])])
        r2 = Species(molecule=[Molecule().fromAdjacencyList(adjlist[1])])
        p1 = Species(molecule=[Molecule().fromAdjacencyList(adjlist[2])])
        p2 = Species(molecule=[Molecule().fromAdjacencyList(adjlist[3])])
        r1.generate_resonance_structures(keepIsomorphic=True)
        p1.generate_resonance_structures(keepIsomorphic=True)
        
        
        rxn = TemplateReaction(reactants = [r1, r2], 
                               products = [p1, p2]
)
        
        rxn.degeneracy = family.calculateDegeneracy(rxn)
        self.assertEqual(rxn.degeneracy, 6)
        
        family.addReverseAttribute(rxn)
        
        self.assertEqual(rxn.reverse.degeneracy, 6)

    def test_generate_reactions_from_families_with_resonance(self):
        """Test that we can generate reactions from families with resonance structures"""
        reactants = [
            Molecule().fromSMILES('CC=C[CH2]'),
            Molecule().fromSMILES('[OH]'),
        ]
        expected_product_1 = Molecule().fromSMILES('CC=CCO')
        expected_product_2 = Molecule().fromSMILES('CC(O)C=C')

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, only_families=['R_Recombination'], resonance=True)

        self.assertEqual(len(reaction_list), 2)

        case_1 = reaction_list[0].products[0].isIsomorphic(expected_product_1) and reaction_list[1].products[0].isIsomorphic(expected_product_2)
        case_2 = reaction_list[0].products[0].isIsomorphic(expected_product_2) and reaction_list[1].products[0].isIsomorphic(expected_product_1)

        # Only one case should be true
        self.assertTrue(case_1 ^ case_2)

    def test_generate_reactions_from_families_no_resonance(self):
        """Test that we can generate reactions from families without resonance structures"""
        reactants = [
            Molecule().fromSMILES('CC=C[CH2]'),
            Molecule().fromSMILES('[OH]'),
        ]
        expected_product = Molecule().fromSMILES('CC=CCO')

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, only_families=['R_Recombination'], resonance=False)

        self.assertEqual(len(reaction_list), 1)

        self.assertTrue(reaction_list[0].products[0].isIsomorphic(expected_product))

    def test_generate_reactions_from_families_product_resonance(self):
        """Test that we can specify the product resonance structure when generating reactions"""
        reactants = [
            Molecule().fromSMILES('CCC=C'),
            Molecule().fromSMILES('[H]'),
        ]
        products = [
            Molecule().fromSMILES('CC=C[CH2]'),
            Molecule().fromSMILES('[H][H]'),
        ]

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, products, only_families=['H_Abstraction'], resonance=True)

        self.assertEqual(len(reaction_list), 1)
        self.assertEqual(reaction_list[0].degeneracy, 2)



    def test_generate_reactions_from_families_product_resonance2(self):
        """Test that we can specify the no product resonance structure when generating reactions"""
        reactants = [
            Molecule().fromSMILES('CCC=C'),
            Molecule().fromSMILES('[H]'),
        ]
        products = [
            Molecule().fromSMILES('CC=C[CH2]'),
            Molecule().fromSMILES('[H][H]'),
        ]

        reaction_list = self.database.kinetics.generate_reactions_from_families(reactants, products, only_families=['H_Abstraction'], resonance=False)
        self.assertEqual(len(reaction_list), 0)

        self.assertTrue(isinstance(products[0],Species))
        self.assertEqual(len(products[0].molecule),1)

    def test_generate_reactions_from_libraries(self):
        """Test that we can generate reactions from libraries"""
        reactants = [
            Molecule().fromSMILES('CC=O'),
            Molecule().fromSMILES('[H]'),
        ]

        reaction_list = self.database.kinetics.generate_reactions_from_libraries(reactants)

        self.assertEqual(len(reaction_list), 3)

    def test_generate_reactions_from_libraries2(self):
        """Test that we can generate reactions from libraries specifying products"""
        reactants = [
            Molecule().fromSMILES('CC=O'),
            Molecule().fromSMILES('[H]'),
        ]
        products = [
            Molecule().fromSMILES('[CH2]C=O'),
            Molecule().fromSMILES('[H][H]'),
        ]
        reaction_list_2 = self.database.kinetics.generate_reactions_from_libraries(reactants, products)

        self.assertEqual(len(reaction_list_2), 1)

    def test_add_atom_labels_for_reaction(self):
        """Test that addAtomLabelsForReaction can identify reactions with resonance
        The molecule [CH]=C=C has resonance in this reaction"""
        from rmgpy.data.rmg import getDB
        reactants = [
            Molecule().fromSMILES('C=C=C'),
            Molecule().fromSMILES('[CH]=C=C'),
        ]
        products = [
            Molecule().fromSMILES('C#C[CH2]'),
            Molecule().fromSMILES('C#CC'),
        ]
        reaction = TemplateReaction(reactants =reactants,
                                    products = products,
                                    family = 'H_Abstraction')
        reaction.ensure_species(reactant_resonance=True, product_resonance=True)
        family = getDB('kinetics').families['H_Abstraction']
        family.addAtomLabelsForReaction(reaction, output_with_resonance=False)

        # test that the reaction has labels
        found_labels = []
        for species in reaction.reactants:
            for atom in species.molecule[0].atoms:
                if atom.label != '':
                    found_labels.append(atom.label)
        self.assertEqual(len(found_labels), 3)
        self.assertIn('*1',found_labels)
        self.assertIn('*2',found_labels)
        self.assertIn('*3',found_labels)

        # test for the products too
        found_labels = []
        for species in reaction.products:
            for atom in species.molecule[0].atoms:
                if atom.label != '':
                    found_labels.append(atom.label)
        self.assertEqual(len(found_labels), 3)
        self.assertIn('*1',found_labels)
        self.assertIn('*2',found_labels)
        self.assertIn('*3',found_labels)

    def test_add_atom_labels_for_reaction_2(self):
        """Test that addAtomLabelsForReaction can identify reactions with identical references
        The molecule [CH]=C=C has resonance in this reaction"""
        from rmgpy.data.rmg import getDB
        s1 = Species().fromSMILES('C=C=C')
        s2 = Species().fromSMILES('C=C=[CH]')
        s3 = Species().fromSMILES('C#CC')
        s2.generate_resonance_structures()
        reactants = [s1,s2]
        products = [s2,s3]
        reaction = TemplateReaction(reactants =reactants,
                                    products = products,
                                    family = 'H_Abstraction')
        family = getDB('kinetics').families['H_Abstraction']
        print reaction.reactants
        print reaction.products
        family.addAtomLabelsForReaction(reaction, output_with_resonance=False)

        # test that the reaction has labels
        found_labels = []
        for species in reaction.reactants:
            for atom in species.molecule[0].atoms:
                if atom.label != '':
                    found_labels.append(atom.label)
        self.assertEqual(len(found_labels), 3,'wrong number of labels found {0}'.format(found_labels))
        self.assertIn('*1',found_labels)
        self.assertIn('*2',found_labels)
        self.assertIn('*3',found_labels)

        # test for the products too
        found_labels = []
        for species in reaction.products:
            for atom in species.molecule[0].atoms:
                if atom.label != '':
                    found_labels.append(atom.label)
        self.assertEqual(len(found_labels), 3)
        self.assertIn('*1',found_labels)
        self.assertIn('*2',found_labels)
        self.assertIn('*3',found_labels)

    def test_add_atom_labels_for_reaction_3(self):
        """Test that addAtomLabelsForReaction can identify reactions with resonance and isotopes"""
        from rmgpy.data.rmg import getDB
        mr0 = Molecule().fromAdjacencyList('1    C u0 p0 c0 i13 {3,D} {4,S} {5,S}\n2 *1 C u0 p0 c0 {3,D} {6,S} {7,S}\n3    C u0 p0 c0 {1,D} {2,D}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {2,S}\n7 *4 H u0 p0 c0 {2,S}\n')
        mr1a = Molecule().fromAdjacencyList('multiplicity 2\n1    C u0 p0 c0 i13 {2,D} {4,S} {5,S}\n2    C u0 p0 c0 {1,D} {3,D}\n3 *1 C u1 p0 c0 {2,D} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n')
        mr1b = Molecule().fromAdjacencyList('multiplicity 2\n1    C u1 p0 c0 i13 {2,S} {4,S} {5,S}\n2    C u0 p0 c0 {1,S} {3,T}\n3 *1 C u0 p0 c0 {2,T} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n')
        mp1a = Molecule().fromAdjacencyList('multiplicity 2\n1    C u0 p0 c0 {2,D} {4,S} {5,S}\n2    C u0 p0 c0 {1,D} {3,D}\n3 *1 C u1 p0 c0 i13 {2,D} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n')
        mp1b = Molecule().fromAdjacencyList('multiplicity 2\n1    C u1 p0 c0 {2,S} {4,S} {5,S}\n2    C u0 p0 c0 {1,S} {3,T}\n3 *1 C u0 p0 c0 i13 {2,T} {6,S}\n4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {1,S}\n6    H u0 p0 c0 {3,S}\n')
        s1 = Species(molecule = [mr0])
        s2 = Species(molecule = [mr1a,mr1b])
        s3 = Species(molecule = [mp1a,mp1b])
        reactants = [s1,s2]
        products = [s1,s3]
        reaction = TemplateReaction(reactants =reactants,
                                    products = products,
                                    family = 'H_Abstraction')
        family = getDB('kinetics').families['H_Abstraction']
        print reaction.reactants
        print reaction.products
        family.addAtomLabelsForReaction(reaction, output_with_resonance=False)

        # test that the reaction has labels
        found_labels = []
        for species in reaction.reactants:
            for atom in species.molecule[0].atoms:
                if atom.label != '':
                    found_labels.append(atom.label)
        self.assertEqual(len(found_labels), 3,'wrong number of labels found {0}'.format(found_labels))
        self.assertIn('*1',found_labels)
        self.assertIn('*2',found_labels)
        self.assertIn('*3',found_labels)

        # test for the products too
        found_labels = []
        for species in reaction.products:
            for atom in species.molecule[0].atoms:
                if atom.label != '':
                    found_labels.append(atom.label)
        self.assertEqual(len(found_labels), 3)
        self.assertIn('*1',found_labels)
        self.assertIn('*2',found_labels)
        self.assertIn('*3',found_labels)

