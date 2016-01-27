import os
import unittest
import logging
from external.wip import work_in_progress 

from .main import RMG, CoreEdgeReactionModel
from .model import Species
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.molecule import Molecule
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
                                       kineticsFamilies=['R_Addition_MultipleBond'])

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
        newReactions.extend(self.rmg.reactionModel.react(self.rmg.database, spc))

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
        newReactions_reverse.extend(self.rmg_dummy.reactionModel.react(self.rmg.database, spc))

        # process newly generated reactions again to make sure no duplicated reactions
        self.rmg_dummy.reactionModel.processNewReactions(newReactions_reverse, spc, None)

        # try to pick out the target reaction 
        target_rxns_reverse = findTargetRxnsContaining(mol_H, mol_C3H2O, \
                                                       self.rmg_dummy.reactionModel.edge.reactions)
        self.assertEqual(len(target_rxns_reverse), 1)

        # whatever order of molecules in spc, the reaction template matched should be same
        self.assertEqual(target_rxns[0].template, target_rxns_reverse[0].template)

        
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
