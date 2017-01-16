import os
import unittest
import logging
from external.wip import work_in_progress 

from .main import RMG, CoreEdgeReactionModel
from rmgpy.species import Species
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.molecule import Molecule
from rmgpy.rmg.react import react
###################################################

class TestRMGWorkFlow(unittest.TestCase):

    def setUp(self):
        """
        A method that is run before each unit test in this class.
        """
        # set-up RMG object
        self.rmg = RMG()
        self.rmg.reactionModel = CoreEdgeReactionModel()

        self.rmg_dummy = RMG()
        self.rmg_dummy.reactionModel = CoreEdgeReactionModel()

        # load kinetic database and forbidden structures
        self.rmg.database = RMGDatabase()
        path = os.path.join(settings['database.directory'])

        # forbidden structure loading
        self.rmg.database.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
        # kinetics family Disproportionation loading
        self.rmg.database.loadKinetics(os.path.join(path, 'kinetics'), \
                                       kineticsFamilies=['R_Addition_MultipleBond'],reactionLibraries=[])

    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
        
    @work_in_progress
    def testDeterministicReactionTemplateMatching(self):
        """
        Test RMG work flow can match reaction template for kinetics estimation 
        deterministically. 

        In this test, a change of molecules order in a reacting species should 
        not change the reaction template matched.
        """

        # react
        spc = Species().fromSMILES("O=C[C]=C")
        spc.generateResonanceIsomers()
        newReactions = []		
        newReactions.extend(react((spc,)))

        # process newly generated reactions to make sure no duplicated reactions
        self.rmg.reactionModel.processNewReactions(newReactions, spc, None)

        # try to pick out the target reaction 
        mol_H = Molecule().fromSMILES("[H]")
        mol_C3H2O = Molecule().fromSMILES("C=C=C=O")

        target_rxns = findTargetRxnsContaining(mol_H, mol_C3H2O, \
                                               self.rmg.reactionModel.edge.reactions)
        self.assertEqual(len(target_rxns), 1)

        # reverse the order of molecules in spc
        spc.molecule = list(reversed(spc.molecule))

        # react again
        newReactions_reverse = []
        newReactions_reverse.extend(react((spc,)))

        # process newly generated reactions again to make sure no duplicated reactions
        self.rmg_dummy.reactionModel.processNewReactions(newReactions_reverse, spc, None)

        # try to pick out the target reaction 
        target_rxns_reverse = findTargetRxnsContaining(mol_H, mol_C3H2O, \
                                                       self.rmg_dummy.reactionModel.edge.reactions)
        self.assertEqual(len(target_rxns_reverse), 1)

        # whatever order of molecules in spc, the reaction template matched should be same
        self.assertEqual(target_rxns[0].template, target_rxns_reverse[0].template)

    def testCheckForExistingSpeciesForBiAromatics(self):
        """
        Test RMG checkForExistingSpecies can correctly check isomorphism for biaromatics. 
        In this test, DPP is a species already stored in rmg speciesDict, mol_test is a newly
        created molecule which has one kekulized benzene ring and one double_bond-single_bond
        benzene ring.
        """

        rmg_test = RMG()
        rmg_test.reactionModel = CoreEdgeReactionModel()
        DPP = Species().fromSMILES('C1=CC=C(C=C1)CCCC1C=CC=CC=1')
        DPP.generateResonanceIsomers()
        formula = DPP.molecule[0].getFormula()
        if formula in rmg_test.reactionModel.speciesDict:
            rmg_test.reactionModel.speciesDict[formula].append(DPP)
        else:
            rmg_test.reactionModel.speciesDict[formula] = [DPP]

        mol_test = Molecule().fromAdjacencyList(
"""
1     C u0 p0 c0 {2,S} {3,S} {16,S} {17,S}
2     C u0 p0 c0 {1,S} {4,S} {18,S} {19,S}
3     C u0 p0 c0 {1,S} {5,S} {20,S} {21,S}
4     C u0 p0 c0 {2,S} {6,B} {7,B}
5     C u0 p0 c0 {3,S} {8,D} {9,S}
6     C u0 p0 c0 {4,B} {10,B} {22,S}
7     C u0 p0 c0 {4,B} {12,B} {24,S}
8     C u0 p0 c0 {5,D} {14,S} {27,S}
9     C u0 p0 c0 {5,S} {15,D} {28,S}
10    C u0 p0 c0 {6,B} {11,B} {23,S}
11    C u0 p0 c0 {10,B} {12,B} {25,S}
12    C u0 p0 c0 {7,B} {11,B} {26,S}
13    C u0 p0 c0 {14,D} {15,S} {29,S}
14    C u0 p0 c0 {8,S} {13,D} {30,S}
15    C u0 p0 c0 {9,D} {13,S} {31,S}
16    H u0 p0 c0 {1,S}
17    H u0 p0 c0 {1,S}
18    H u0 p0 c0 {2,S}
19    H u0 p0 c0 {2,S}
20    H u0 p0 c0 {3,S}
21    H u0 p0 c0 {3,S}
22    H u0 p0 c0 {6,S}
23    H u0 p0 c0 {10,S}
24    H u0 p0 c0 {7,S}
25    H u0 p0 c0 {11,S}
26    H u0 p0 c0 {12,S}
27    H u0 p0 c0 {8,S}
28    H u0 p0 c0 {9,S}
29    H u0 p0 c0 {13,S}
30    H u0 p0 c0 {14,S}
31    H u0 p0 c0 {15,S}
""")
        found, spec = rmg_test.reactionModel.checkForExistingSpecies(mol_test)
        assert found == True

        
def findTargetRxnsContaining(mol1, mol2, reactions):
    target_rxns = []
    for rxn in reactions:
        reactants = rxn.reactants
        products = rxn.products
        rxn_specs = reactants + products
        for rxn_spec in rxn_specs:
            if rxn_spec.isIsomorphic(mol1):
                for rxn_spec1 in rxn_specs:
                    if rxn_spec1.isIsomorphic(mol2):
                        target_rxns.append(rxn)
    return target_rxns
