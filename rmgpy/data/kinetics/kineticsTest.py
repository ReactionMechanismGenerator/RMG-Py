import os
import unittest 
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.base import DatabaseError
import numpy
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.rmg import RMGDatabase
###################################################

def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.load(
        path=settings['database.directory'],
        thermoLibraries=['primaryThermoLibrary'],
        reactionLibraries=[],
        kineticsFamilies=[
            'R_Recombination',
            'Disproportionation',
            'R_Addition_MultipleBond',
            'intra_H_migration'
        ],
    )
    for family in database.kinetics.families.values():
        family.addKineticsRulesFromTrainingSet(thermoDatabase=database.thermo)
        family.fillKineticsRulesByAveragingUp(verbose=True)

def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg
    rmg.database = None

class TestKineticsDatabase(unittest.TestCase):
    
    def testLoadFamilies(self):
        """
        Test that the loadFamilies function raises the correct exceptions
        """
        path = os.path.join(settings['database.directory'],'kinetics','families')
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

    def test_degeneracy_for_methyl_methyl_recombination(self):
        """Test that the proper degeneracy is calculated for methyl + methyl recombination"""

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

        correct_degeneracy = 6
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

        reactions = family.generateReactions(reactants)
        self.assertEqual(len(reactions), num_independent_reactions,'only {1} reaction(s) should be produced. Produced reactions {0}'.format(reactions,num_independent_reactions))

        return sum([reaction.degeneracy for reaction in reactions]), reactions

    def test_propyl_propyl_reaction_is_the_same_as_propyl_butyl(self):
        """
        test that propyl propyl r-recombination is the same rate as propyl butyl
        
        this test assures that r-recombination reactions from the same rate rule
        have the same reaction rate since they have both symmetrical transition
        states and reactants, which should cancel out in TST
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
        self.assertIn(pp_kinetics_list[0][0].comment,pb_kinetics_list[0][0].comment,
                         'this test found different kinetics for the two groups, so it will not function as expected\n' + 
                         str(pp_kinetics_list)+str(pb_kinetics_list))
        
        # test that the kinetics are correct
        self.assertAlmostEqual(pp_kinetics_list[0][0].getRateCoefficient(300), pb_kinetics_list[0][0].getRateCoefficient(300))

    def test_identical_reactants_have_faster_kinetics(self):
        """
        tests identical reactants have faster kinetics that different reactants.
        
        this test assures that r addition multiple bond reactions from the same 
        rate rule have the faster reaction rate if the reactants are identicaal 
        since they have symmetrical reactants, with little change in the 
        transition state symmetry. This should be more robust than just checking
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
        self.assertNotEqual(pp_kinetics_list[0][0].getRateCoefficient(300), pb_kinetics_list[0][0].getRateCoefficient(300))
        self.assertAlmostEqual(pp_kinetics_list[0][0].getRateCoefficient(300) / 2, pb_kinetics_list[0][0].getRateCoefficient(300))
        
class TestKineticsCommentsParsing(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database
        
    def testParseKinetics(self):
        from rmgpy.chemkin import loadChemkinFile
        import rmgpy
        species, reactions = loadChemkinFile(os.path.join(os.path.dirname(rmgpy.__file__), 'data','kinetics','parsing_data','chem_annotated.inp'),
                                             os.path.join(os.path.dirname(rmgpy.__file__), 'data','kinetics','parsing_data','species_dictionary.txt')
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
        self.assertEqual(training[0][1].index,1)  # Assert that the index of that training reaction is 1
        
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
        