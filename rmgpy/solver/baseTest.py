import unittest
import pickle
import os.path

from rmgpy.tools.loader import loadRMGPyJob
import rmgpy

from rmgpy.solver.base import *

class ConcentrationPrinter:
    def __init__(self):
        self.species_names = []
        self.data = []

    def update(self, subject):
        self.data.append((subject.t , subject.coreSpeciesConcentrations))

class ReactionSystemTest(unittest.TestCase):

    def setUp(self):
        self.listener = ConcentrationPrinter()

        folder = os.path.join(os.path.dirname(rmgpy.__file__),'solver/files/listener/')
        inputFile = os.path.join(folder, 'input.py')
        chemkinFile = os.path.join(folder, 'chemkin/chem.inp')
        spc_dict = os.path.join(folder, 'chemkin/species_dictionary.txt')

        self.rmg = loadRMGPyJob(inputFile, chemkinFile, spc_dict, generateImages=False)


    def testAttachDetach(self):
        """
        Test that a ReactionSystem listener can be attached/detached.
        """
        #create observable

        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        self.assertNotEqual(reactionSystem._observers, [])
        
        reactionSystem.detach(self.listener)    
        self.assertEquals(reactionSystem._observers, [])

    def testListen(self):
        """
        Test that data can be retrieved from an attached ReactionSystem listener.
        """
        #create observable
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)

        reactionModel = self.rmg.reactionModel

        self.assertEqual(self.listener.data, [])

        # run simulation:
        terminated, obj,sspcs,srxns = reactionSystem.simulate(
            coreSpecies = reactionModel.core.species,
            coreReactions = reactionModel.core.reactions,
            edgeSpecies = reactionModel.edge.species,
            edgeReactions = reactionModel.edge.reactions,
            surfaceSpecies = [],
            surfaceReactions = [],
            toleranceKeepInEdge = 0,
            toleranceMoveToCore = 1,
            toleranceInterruptSimulation = 1,
        ) 

        self.assertNotEqual(self.listener.data, [])

    def testPickle(self):
        """
        Test that a ReactionSystem object can be un/pickled.
        """
        rxnSys1 = self.rmg.reactionSystems[0]
        rxnSys = pickle.loads(pickle.dumps(rxnSys1))

        self.assertIsNotNone(rxnSys)
        self.assertTrue(isinstance(rxnSys, rmgpy.solver.simple.SimpleReactor))
        self.assertEqual(rxnSys.T.value_si, rxnSys1.T.value_si)
        self.assertEqual(rxnSys.P.value_si, rxnSys1.P.value_si)
        self.assertEqual(rxnSys.termination[0].conversion, rxnSys1.termination[0].conversion)
        self.assertEqual(rxnSys.termination[1].time.value_si, rxnSys1.termination[1].time.value_si)
        
        
if __name__ == '__main__':
    unittest.main()