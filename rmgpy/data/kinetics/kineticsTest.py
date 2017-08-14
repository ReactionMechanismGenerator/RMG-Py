################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os
import unittest 
from external.wip import work_in_progress
import itertools

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.base import DatabaseError
import numpy
from rmgpy.molecule.molecule import Molecule
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.react import findDegeneracies, react, reactSpecies, _labelListOfSpecies
from rmgpy.data.base import ForbiddenStructures
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
            'H_Abstraction'
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

    def testR_Addition_MultipleBondBenzene(self):
        """Test that the proper degeneracy is calculated for H addition to benzene"""
        family = 'R_Addition_MultipleBond'
        reactants = [
            Molecule().fromSMILES('c1ccccc1'),
            Molecule().fromSMILES('[H]'),
        ]
        # assign atom IDs
        for reactant in reactants: reactant.assignAtomIDs()

        reactants = [mol.generateResonanceIsomers() for mol in reactants]

        combinations = itertools.product(reactants[0], reactants[1])

        reactionList = []
        for combi in combinations:
            reactionList.extend(self.database.kinetics.families[family].generateReactions(combi))

        reactionList = findDegeneracies(reactionList)

        self.assertEqual(len(reactionList), 1)
        for rxn in reactionList:
            self.assertEqual(rxn.degeneracy, 6)

    def testR_Addition_MultipleBondMethylNaphthalene(self):
        """Test that the proper degeneracy is calculated for H addition to methylnaphthalene"""
        family = 'R_Addition_MultipleBond'
        reactants = [
            Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1C'),
            Molecule().fromSMILES('[H]'),
        ]
        # assign atom IDs
        for reactant in reactants: reactant.assignAtomIDs()
        
        reactants = [mol.generateResonanceIsomers() for mol in reactants]

        combinations = itertools.product(reactants[0], reactants[1])

        reactionList = []
        for combi in combinations:
            reactionList.extend(self.database.kinetics.families[family].generateReactions(combi))

        product = Species().fromSMILES('C[C]1CC=CC2=CC=CC=C12')
        product.generateResonanceIsomers()

        targetReactions = []
        for rxn in reactionList:
            for spc in rxn.products:
                if product.isIsomorphic(spc):
                    targetReactions.append(rxn)

        targetReactions = findDegeneracies(targetReactions)

        self.assertEqual(len(targetReactions), 1)
        for rxn in targetReactions:
            self.assertEqual(rxn.degeneracy, 1)

    def testR_RecombinationPhenyl(self):
        """Test that the proper degeneracy is calculated for phenyl + H recombination"""
        family = 'R_Recombination'
        reactants = [
            Molecule().fromSMILES('[c]1ccccc1'),
            Molecule().fromSMILES('[H]'),
        ]

        # assign atom IDs
        for reactant in reactants: reactant.assignAtomIDs()

        reactants = [mol.generateResonanceIsomers() for mol in reactants]

        combinations = itertools.product(reactants[0], reactants[1])

        reactionList = []
        for combi in combinations:
            reactionList.extend(self.database.kinetics.families[family].generateReactions(combi))

        reactionList = findDegeneracies(reactionList)

        self.assertEqual(len(reactionList), 1)
        for rxn in reactionList:
            self.assertEqual(rxn.degeneracy, 1)

    def testR_RecombinationH(self):
        """Test that the proper degeneracy is calculated for H + H recombination"""
        family = 'R_Recombination'
        reactants = [
            Molecule().fromSMILES('[H]'),
            Molecule().fromSMILES('[H]'),
        ]
        for reactant in reactants: reactant.assignAtomIDs()

        reactionList = self.database.kinetics.families[family].generateReactions(reactants)

        reactionList = findDegeneracies(reactionList)

        self.assertEqual(len(reactionList), 1)
        self.assertEqual(reactionList[0].degeneracy, 1)

    def test_degeneracy_for_methyl_methyl_recombination(self):
        """Test that the proper degeneracy is calculated for methyl + methyl recombination"""

        correct_degeneracy = 0.5
        rxn_family_str = 'R_Recombination'
        adj_lists = [
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

        self.compare_degeneracy_of_reaction(adj_lists,rxn_family_str,correct_degeneracy)

    def test_degeneracy_for_methyl_labeled_methyl_recombination(self):
        """Test that the proper degeneracy is calculated for methyl + labeled methyl recombination"""

        correct_degeneracy = 1
        rxn_family_str = 'R_Recombination'
        adj_lists = [
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

        self.compare_degeneracy_of_reaction(adj_lists,rxn_family_str,correct_degeneracy)

    def test_degeneracy_for_ethyl_ethyl_disproportionation(self):
        """Test that the proper degeneracy is calculated for ethyl + ethyl disproportionation"""

        correct_degeneracy = 3
        rxn_family_str = 'Disproportionation'
        adj_lists = [
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

        self.compare_degeneracy_of_reaction(adj_lists,rxn_family_str,correct_degeneracy)

    def test_degeneracy_for_ethyl_labeled_ethyl_disproportionation(self):
        """Test that the proper degeneracy is calculated for ethyl + labeled ethyl disproportionation"""

        correct_degeneracy = 3
        rxn_family_str = 'Disproportionation'
        adj_lists = [
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
        expected_products = 2
        self.compare_degeneracy_of_reaction(adj_lists,rxn_family_str,correct_degeneracy * expected_products,expected_products)

    def compare_degeneracy_of_reaction(self, reactants_adj_list,
                                       rxn_family_str,
                                       num_expected_degenerate_products,
                                       num_independent_reactions = 1):
        """
        given:

        `reactants_adj_list`: a list of adjacency lists (of reactants)
        `reaction_family_str`: the string representation of the reaction family
        `num_expected_degenerate_products`: the total number of degenerate reactions
                        which should be found by generateReactions.
        `num_independent_rxns`: the number of reaction objects expected from generateReactions

        performs:

        a check to ensure that the number of degenerate reactions is what is
        expected.
        """

        found_degeneracy, reaction = self.find_reaction_degeneracy(reactants_adj_list,rxn_family_str,
                                 num_independent_reactions)
        self.assertEqual(found_degeneracy, num_expected_degenerate_products,'degeneracy returned ({0}) is not the correct value ({1}) for reaction {2}'.format(found_degeneracy, num_expected_degenerate_products,reaction))

    def find_reaction_degeneracy(self, reactants_adj_list,rxn_family_str,
                                 num_independent_reactions = 1):
        """
        given:

        reactants_adj_list: a list of adjacency lists of the reactants
        `reaction_family_str`: the string representation of the reaction family
        `num_independent_rxns`: the number of reaction objects expected from generateReactions

        returns:

        a tuple with the total degeneracy and a list of reaction objects
        """
        family = self.database.kinetics.families[rxn_family_str]
        reactants = [Molecule().fromAdjacencyList(reactants_adj_list[0]),
                     Molecule().fromAdjacencyList(reactants_adj_list[1])]

        for reactant in reactants: reactant.assignAtomIDs()
        reactions = family.generateReactions(reactants)
        reactions = findDegeneracies(reactions)
        self.assertEqual(len(reactions), num_independent_reactions,'only {1} reaction(s) should be produced. Produced reactions {0}'.format(reactions,num_independent_reactions))

        return sum([reaction.degeneracy for reaction in reactions]), reactions

    def test_degeneracy_does_not_include_identical_atom_labels(self):
        """
        ensure rxns with identical atom_ids are not counted twice for degeneracy
        
        this test uses [H] + CC=C[CH]C -> H2 + [CH2]C=C[CH]C as an example. Since
        the reactant is symmetric with the middle carbon, the degeneracy should be
        6.
        """
        spcA = Species().fromSMILES('[H]')
        spcB = Species().fromSMILES('CC=C[CH]C')
        spcB.generateResonanceIsomers(keepIsomorphic=True)
        spcTuples = [(spcA,spcB)]
        
        reactionList = list(react(*spcTuples))
        
        # find reaction with a specific product
        specific_product = Species().fromSMILES('[CH2]C=C[CH]C')
        
        specific_product.generateResonanceIsomers()
        
        specific_reaction = None
        for rxn in reactionList:
            if any([specific_product.isIsomorphic(product) for product in rxn.products]):
                specific_reaction = rxn
                break
        self.assertIsNotNone(specific_reaction,'no reaction found with the specified product')
        
        self.assertEqual(specific_reaction.degeneracy, 6,'The reaction output the wrong degeneracy of {}.'.format(specific_reaction.degeneracy))
    def test_degeneracy_keeps_separate_transition_states_separated(self):
        """
        ensure rxns with multiple transition states are kept as separate reactions
        
        this test uses C[C]=C + C=C[CH2] -> C=C=C + C=CC as an example. 
        This reaction should have two transition states, which should occur regardless
        of the order .
        """
        spcA = Species().fromSMILES('C[C]=C')
        spcB = Species().fromSMILES('C=C[CH2]')
        spcTuples = [(spcA,spcB)]
        reactionList = list(react(*spcTuples))
        # find reaction with a specific product
        specific_products = [Species().fromSMILES('C=C=C'),
                             Species().fromSMILES('CC=C'),]
        
        # eliminate rxns that do not match products
        isomorphic_rxns = 0
        for rxn in reactionList:
            #  rxn contains all products
            if all([any([specific_product.isIsomorphic(product) for product in rxn.products]) for specific_product in specific_products]):
                isomorphic_rxns += 1

        self.assertEqual(isomorphic_rxns, 2,'The reaction output did not output all the transition states in either order of reactants')
     
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

    def test_degeneracy_keeps_track_of_both_rate_rules_from_resonance_isomers(self):
        """
        rxns that have multiple resonance structures hitting different rate rules should 
        be kept separate when findDegeneracy is used.

        this test uses [H] + CC=C[CH]C -> H2 + [CH2]C=C[CH]C as an example. 
        This reaction should have two transition states.
        """
        spcA = Species().fromSMILES('[H]')
        spcB = Species().fromSMILES('CC=C[CH]C')
        spcB.generateResonanceIsomers(keepIsomorphic=True)
        spcTuples = [(spcA,spcB)]
        
        reactionList = list(react(*spcTuples))
        
        # find reaction with a specific product
        specific_product = Species().fromSMILES('CC=C[CH][CH2]')
        specific_product.generateResonanceIsomers()
        
        specific_reactions_found = 0
        templates_found = []
        for rxn in reactionList:
            if any([specific_product.isIsomorphic(product) for product in rxn.products]):
                specific_reactions_found += 1
                templates_found.append(rxn.template)
        
        self.assertEqual(specific_reactions_found, 2,'The reaction output did not contain 2 transition states.')
        self.assertNotEqual(templates_found[0],templates_found[1],'The reactions should have different templates')

    def test_propyl_propyl_reaction_is_the_half_propyl_butyl(self):
        """
        test that propyl propyl r-recombination is the same rate as propyl butyl

        this test assures that r-recombination reactions from the same rate rule
        with identical reactants have half the reaction rate since there is a 
        symmetrical transition state.
        """
        rxn_family_str = 'R_Recombination'
        propyl_adj_list = """
            multiplicity 2
            1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
            2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
            3  C u1 p0 c0 {2,S} {4,S} {5,S}
            4  H u0 p0 c0 {3,S}
            5  H u0 p0 c0 {3,S}
            6  H u0 p0 c0 {1,S}
            7  H u0 p0 c0 {1,S}
            8  H u0 p0 c0 {1,S}
            9  H u0 p0 c0 {2,S}
            10 H u0 p0 c0 {2,S}

            """
        butyl_adj_list = """
            multiplicity 2
            1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
            2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
            3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
            4  C u1 p0 c0 {3,S} {5,S} {6,S}
            5  H u0 p0 c0 {4,S}
            6  H u0 p0 c0 {4,S}
            7  H u0 p0 c0 {1,S}
            8  H u0 p0 c0 {1,S}
            9  H u0 p0 c0 {1,S}
            10 H u0 p0 c0 {2,S}
            11 H u0 p0 c0 {2,S}
            12 H u0 p0 c0 {3,S}
            13 H u0 p0 c0 {3,S}
            """

        family = self.database.kinetics.families[rxn_family_str]

        # get reaction objects and their degeneracy
        pp_degeneracy, pp_reactions = self.find_reaction_degeneracy([propyl_adj_list,propyl_adj_list],rxn_family_str)
        pb_degeneracy, pb_reactions = self.find_reaction_degeneracy([propyl_adj_list,butyl_adj_list],rxn_family_str)

        # since output is a list of 1
        pp_reaction = pp_reactions[0]
        pb_reaction = pb_reactions[0]

        # get kinetics for each reaction
        pp_kinetics_list = family.getKinetics(pp_reaction, pp_reaction.template,
                                              degeneracy=pp_reaction.degeneracy,
                                              estimator = 'rate rules')
        self.assertEqual(len(pp_kinetics_list), 1, 'The propyl and propyl recombination should only return one reaction. It returned {0}. Here is the full kinetics: {1}'.format(len(pp_kinetics_list),pp_kinetics_list))

        pb_kinetics_list = family.getKinetics(pb_reaction, pb_reaction.template,
                                              degeneracy=pb_reaction.degeneracy,
                                              estimator = 'rate rules')
        self.assertEqual(len(pb_kinetics_list), 1, 'The propyl and butyl recombination should only return one reaction. It returned {0}. Here is the full kinetics: {1}'.format(len(pb_kinetics_list),pb_kinetics_list))

        # the same reaction group must be found or this test will not work
        self.assertIn(pb_kinetics_list[0][0].comment,pp_kinetics_list[0][0].comment,
                         'this test found different kinetics for the two groups, so it will not function as expected\n' +
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
        rxn_family_str = 'R_Addition_MultipleBond'
        butenyl_adj_list = """
            multiplicity 2
            1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
            2 C u0 p0 c0 {1,S} {4,D} {7,S}
            3 C u1 p0 c0 {1,S} {8,S} {9,S}
            4 C u0 p0 c0 {2,D} {10,S} {11,S}
            5 H u0 p0 c0 {1,S}
            6 H u0 p0 c0 {1,S}
            7 H u0 p0 c0 {2,S}
            8 H u0 p0 c0 {3,S}
            9 H u0 p0 c0 {3,S}
            10 H u0 p0 c0 {4,S}
            11 H u0 p0 c0 {4,S}
            """
        pentenyl_adj_list = """
            multiplicity 2
            1 C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
            2 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
            3 C u0 p0 c0 {1,S} {5,D} {10,S}
            4 C u1 p0 c0 {2,S} {11,S} {12,S}
            5 C u0 p0 c0 {3,D} {13,S} {14,S}
            6 H u0 p0 c0 {2,S}
            7 H u0 p0 c0 {2,S}
            8 H u0 p0 c0 {1,S}
            9 H u0 p0 c0 {1,S}
            10 H u0 p0 c0 {3,S}
            11 H u0 p0 c0 {4,S}
            12 H u0 p0 c0 {4,S}
            13 H u0 p0 c0 {5,S}
            14 H u0 p0 c0 {5,S}
            """

        family = self.database.kinetics.families[rxn_family_str]

        # get reaction objects and their degeneracy
        pp_degeneracy, pp_reactions = self.find_reaction_degeneracy([butenyl_adj_list,butenyl_adj_list],rxn_family_str, num_independent_reactions=2)
        pb_degeneracy, pb_reactions = self.find_reaction_degeneracy([butenyl_adj_list,pentenyl_adj_list],rxn_family_str, num_independent_reactions=4)

        # find the correct reaction from the list
        symmetric_product=Molecule().fromAdjacencyList('''
            multiplicity 3
            1 C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
            2 C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
            3 C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
            4 C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
            5 C u0 p0 c0 {4,S} {8,D} {16,S}
            6 C u1 p0 c0 {1,S} {19,S} {20,S}
            7 C u1 p0 c0 {3,S} {17,S} {18,S}
            8 C u0 p0 c0 {5,D} {21,S} {22,S}
            9 H u0 p0 c0 {1,S}
            10 H u0 p0 c0 {2,S}
            11 H u0 p0 c0 {2,S}
            12 H u0 p0 c0 {3,S}
            13 H u0 p0 c0 {3,S}
            14 H u0 p0 c0 {4,S}
            15 H u0 p0 c0 {4,S}
            16 H u0 p0 c0 {5,S}
            17 H u0 p0 c0 {7,S}
            18 H u0 p0 c0 {7,S}
            19 H u0 p0 c0 {6,S}
            20 H u0 p0 c0 {6,S}
            21 H u0 p0 c0 {8,S}
            22 H u0 p0 c0 {8,S}
            ''')
        asymmetric_product = Molecule().fromAdjacencyList('''
            multiplicity 3
            1 C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
            2 C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
            3 C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
            4 C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
            5 C u0 p0 c0 {2,S} {8,S} {15,S} {16,S}
            6 C u0 p0 c0 {4,S} {9,D} {19,S}
            7 C u1 p0 c0 {1,S} {22,S} {23,S}
            8 C u1 p0 c0 {5,S} {20,S} {21,S}
            9 C u0 p0 c0 {6,D} {24,S} {25,S}
            10 H u0 p0 c0 {1,S}
            11 H u0 p0 c0 {2,S}
            12 H u0 p0 c0 {2,S}
            13 H u0 p0 c0 {3,S}
            14 H u0 p0 c0 {3,S}
            15 H u0 p0 c0 {5,S}
            16 H u0 p0 c0 {5,S}
            17 H u0 p0 c0 {4,S}
            18 H u0 p0 c0 {4,S}
            19 H u0 p0 c0 {6,S}
            20 H u0 p0 c0 {8,S}
            21 H u0 p0 c0 {8,S}
            22 H u0 p0 c0 {7,S}
            23 H u0 p0 c0 {7,S}
            24 H u0 p0 c0 {9,S}
            25 H u0 p0 c0 {9,S}
            ''')

        pp_reaction = next((reaction for reaction in pp_reactions if reaction.products[0].isIsomorphic(symmetric_product)),None)
        pb_reaction = next((reaction for reaction in pb_reactions if reaction.products[0].isIsomorphic(asymmetric_product)),None)

        pp_kinetics_list = family.getKinetics(pp_reaction, pp_reaction.template,
                                              degeneracy=pp_reaction.degeneracy,
                                              estimator = 'rate rules')
        self.assertEqual(len(pp_kinetics_list), 1, 'The propyl and propyl recombination should only return one reaction. It returned {0}. Here is the full kinetics: {1}'.format(len(pp_kinetics_list),pp_kinetics_list))

        pb_kinetics_list = family.getKinetics(pb_reaction, pb_reaction.template,
                                              degeneracy=pb_reaction.degeneracy,
                                              estimator = 'rate rules')
        self.assertEqual(len(pb_kinetics_list), 1,  'The propyl and butyl recombination should only return one reaction. It returned {0}. Here is the full kinetics: {1}'.format(len(pb_kinetics_list),pb_kinetics_list))

        # the same reaction group must be found or this test will not work
        self.assertIn(pb_kinetics_list[0][0].comment,pp_kinetics_list[0][0].comment,
                         'this test found different kinetics for the two groups, so it will not function as expected\n' +
                         str(pp_kinetics_list)+str(pb_kinetics_list))

        # test that the kinetics are correct
        self.assertAlmostEqual(pp_kinetics_list[0][0].getRateCoefficient(300), pb_kinetics_list[0][0].getRateCoefficient(300))
        
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

        forward_reactions = findDegeneracies(forward_reactions)
        reverse_reactions = findDegeneracies(reverse_reactions)

        self.assertEqual(forward_reactions[0].degeneracy, reverse_reactions[0].degeneracy,
                         'the kinetics from forward and reverse directions had different degeneracies, {} and {} respectively'.format(forward_reactions[0].degeneracy, reverse_reactions[0].degeneracy))

class TestKineticsCommentsParsing(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database

    def testParseKinetics(self):
        from rmgpy.chemkin import loadChemkinFile
        import rmgpy
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
        from rmgpy.chemkin import loadChemkinFile
        """A function that is run ONCE before all unit tests in this class."""

        global database
        self.database = database
        
        for family in self.database.kinetics.families.values():
            family.addKineticsRulesFromTrainingSet(thermoDatabase=self.database.thermo)    
            family.fillKineticsRulesByAveragingUp(verbose=True)
    
        self.species, self.reactions = loadChemkinFile(os.path.join(settings['test_data.directory'], 'parsing_data','chem_annotated.inp'),
                                             os.path.join(settings['test_data.directory'], 'parsing_data','species_dictionary.txt')
                                                    )
        
    def testFilterReactions(self):
        """
        tests that filter reactions removes reactions that are missing
        any reactants or products
        """
        
        from rmgpy.data.kinetics.common import filterReactions

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

        out = filterReactions(newmreactants,newmproducts,reactions)
            
        rset = list(set(reactions) - set(out))

        msets = [set(i.reactants+i.products) for i in rset]
        
        for i, iset in enumerate(msets): #test that all the reactions we removed are missing a reactant or a product
            self.assertTrue(iset & lrset != set(),msg='reaction {0} removed improperly'.format(rset[i]))
        
        outsets = [set(i.reactants+i.products) for i in out]
            
        for i, iset in enumerate(outsets): #test that all the reactions left in aren't missing any reactants or products
            self.assertTrue(iset & lrset == set(),msg='reaction {0} left in improperly, should have removed in based on presence of {1}'.format(out[i],iset & lrset))
        
        
    def testSaveEntry(self):
        """
        tests that save entry can run
        """
        from rmgpy.data.kinetics.common import saveEntry
        import os
        from rmgpy.data.base import Entry
        
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
        r1.generateResonanceIsomers(keepIsomorphic=True)
        p1.generateResonanceIsomers(keepIsomorphic=True)
        
        
        rxn = TemplateReaction(reactants = [r1, r2], 
                               products = [p1, p2]
)
        
        rxn.degeneracy = family.calculateDegeneracy(rxn)
        self.assertEqual(rxn.degeneracy, 6)
        
        family.addReverseAttribute(rxn)
        
        self.assertEqual(rxn.reverse.degeneracy, 6)