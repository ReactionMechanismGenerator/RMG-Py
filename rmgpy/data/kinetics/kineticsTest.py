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
        
        
    def compare_degeneracy_of_reaction(self, reactants_adj_list, rxn_family_str, correct_degeneracy_value):
        """
        input a list of adjacency lists (of reactants), 
        the reaction family name,
        and the correct degeneracy value
        """
        
        self.database.loadFamilies(self.path, families=[rxn_family_str])
        family = self.database.families[rxn_family_str]
        reactants = [Molecule().fromAdjacencyList(reactants_adj_list[0]),
                     Molecule().fromAdjacencyList(reactants_adj_list[1])] 
        reactions = family.generateReactions(reactants)
        self.assertEqual(len(reactions), 1,'only one reaction should be produced. Produced reactions {0}'.format(reactions))
        reaction = reactions[0]
        degeneracy = family.calculateDegeneracy(reaction)
        self.assertEqual(degeneracy, correct_degeneracy_value,'degeneracy returned ({0}) is not the correct value ({1}) for reaction {2}'.format(degeneracy, correct_degeneracy_value,reaction)) 

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
        