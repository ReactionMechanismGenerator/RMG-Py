import os
import unittest 
from external.wip import work_in_progress
from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.base import DatabaseError
import numpy
from rmgpy.molecule.molecule import Molecule
###################################################

class TestKineticsDatabase(unittest.TestCase):
    
        
    def setUp(self):
        self.path = os.path.join(settings['database.directory'],'kinetics','families')
        self.database = KineticsDatabase()
    
    
    def test_load_families_throws_proper_exceptions(self):
        """
        Test that the loadFamilies function raises the correct exceptions
        """
        
        
        with self.assertRaises(DatabaseError):
            self.database.loadFamilies(self.path, families='random')
        with self.assertRaises(DatabaseError):
            self.database.loadFamilies(self.path, families=['!H_Abstraction','Disproportionation'])
        with self.assertRaises(DatabaseError):
            self.database.loadFamilies(self.path, families=['fake_family'])
        with self.assertRaises(DatabaseError):
            self.database.loadFamilies(self.path, families=[])
            
    def test_proper_degeneracy_calculated_for_methyl_methyl_recombindation(self):

        correct_degeneracy = 2
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
        
        
        
    def test_proper_degeneracy_calculated_for_methyl_enriched_methyl_recombindation(self):

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
        
    def test_proper_degeneracy_calculated_for_ethyl_ethyl_disproportionation(self):

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
        
    def test_proper_degeneracy_calculated_for_ethyl_labeled_ethyl_disproportionation(self):

        correct_degeneracy = 6
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
        self.compare_degeneracy_of_reaction(adj_lists,rxn_family_str,correct_degeneracy,expected_products)
        
    def compare_degeneracy_of_reaction(self, reactants_adj_list, 
                                       rxn_family_str, 
                                       correct_degeneracy_value,
                                       num_products = 1):
        """
        input a list of adjacency lists (of reactants), 
        the reaction family name,
        and the correct degeneracy value
        """
        
        degeneracy, reaction = self.find_reaction_degeneracy(reactants_adj_list,rxn_family_str,
                                 num_products)
        self.assertEqual(degeneracy, correct_degeneracy_value,'degeneracy returned ({0}) is not the correct value ({1}) for reaction {2}'.format(degeneracy, correct_degeneracy_value,reaction)) 

    def find_reaction_degeneracy(self, reactants_adj_list,rxn_family_str,
                                 num_products = 1):
        """
        given the adjacency list of reactants and the reaction family, returns
        a tuple with the degeneracy and the reaction object.
        """
        self.database.loadFamilies(self.path, families=[rxn_family_str])
        family = self.database.families[rxn_family_str]
        reactants = [Molecule().fromAdjacencyList(reactants_adj_list[0]),
                     Molecule().fromAdjacencyList(reactants_adj_list[1])] 

        reactions = family.generateReactions(reactants)
        self.assertEqual(len(reactions), num_products,'only {1} reaction(s) should be produced. Produced reactions {0}'.format(reactions,num_products))

        return sum([family.calculateDegeneracy(reaction) for reaction in reactions]), reaction

    def test_propyl_propyl_reaction_is_faster_than_propyl_butyl(self):
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
        
        self.database.loadFamilies(self.path, families=[rxn_family_str])
        family = self.database.families[rxn_family_str]
        
        # get reaction objects and their degeneracy
        pp_degeneracy, pp_reaction = self.find_reaction_degeneracy([propyl_adj_list,propyl_adj_list],rxn_family_str)
        pb_degeneracy, pb_reaction = self.find_reaction_degeneracy([propyl_adj_list,butyl_adj_list],rxn_family_str)
        
        #label reaction objects
        family.addAtomLabelsForReaction(pp_reaction)
        family.addAtomLabelsForReaction(pb_reaction)
        
        # get kinetics for each reaction
        template = family.getReactionTemplateLabels(pp_reaction)
        pp_kinetics_list = family.getKinetics(pp_reaction, template,
                                              degeneracy=pp_degeneracy,
                                              estimator = 'rate rules')
        self.assertEqual(len(pp_kinetics_list), 1, pp_kinetics_list)
        
        template = family.getReactionTemplateLabels(pb_reaction)
        pb_kinetics_list = family.getKinetics(pb_reaction, template,
                                              degeneracy=pb_degeneracy,
                                              estimator = 'rate rules')
        self.assertEqual(len(pb_kinetics_list), 1, pb_kinetics_list)
        
        # the same reaction group must be found or this test will not work
        self.assertEqual(pp_kinetics_list[0][0].comment,pb_kinetics_list[0][0].comment,
                         'this test found different kinetics for the two groups, so it will not function as expected')
        
        # test that the kinetics are correct
        self.assertNotEqual(pp_kinetics_list[0][0].getRateCoefficient(300), pb_kinetics_list[0][0].getRateCoefficient(300))
        self.assertEqual(pp_kinetics_list[0][0].getRateCoefficient(300) * 2, pb_kinetics_list[0][0].getRateCoefficient(300))
        
        
        
class TestKineticsCommentsParsing(unittest.TestCase):
    from rmgpy.data.rmg import RMGDatabase 

    database=RMGDatabase()
    database.load(settings['database.directory'], 
                      kineticsFamilies=['Disproportionation'], 
                      kineticsDepositories=[],
                      thermoLibraries=['primaryThermoLibrary'],   # Use just the primary thermo library, which contains necessary small molecular thermo
                      reactionLibraries=[],
                      )

    # Prepare the database by loading training reactions but not averaging the rate rules
    for family in database.kinetics.families.values():
        family.addKineticsRulesFromTrainingSet(thermoDatabase=database.thermo)    
        family.fillKineticsRulesByAveragingUp(verbose=True)
        
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        
        self.database = self.__class__.database
        
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
        