
import unittest

from .model import *


class MockMolecule(object):
    """docstring for MockMolecule"""
    def __init__(self, label):
        super(MockMolecule, self).__init__()
        self.label = label
        
class ReductionReactionTest(unittest.TestCase):

    def setUp(self):
        from rmgpy.reaction import Reaction
        from .model import ReductionReaction

        mol1 = MockMolecule(label='mol1')
        mol2 = MockMolecule(label='mol2')
        mol3 = MockMolecule(label='mol3')
        mol4 = MockMolecule(label='mol4')
        
        self.rxn = Reaction(reactants=[mol1, mol2], products=[mol3, mol4])
        
        self.rrxn = ReductionReaction(self.rxn)


    def tearDown(self):
        del self.rrxn


    def testConstructor(self):
        rrxn = self.rrxn
        rxn = self.rxn

        self.assertIsNotNone(rrxn)

        # attributes
        self.assertIsNotNone(rrxn.reactants, rxn.reactants)
        self.assertIs(rrxn.products, rxn.products)
        self.assertIs(rrxn.rmgReaction, rxn)
        self.assertIsNotNone(rrxn.stoichio)
        self.assertIsNone(rrxn.kf)
        self.assertIsNone(rrxn.kb)


        # stoichio
        for k,d in self.rrxn.stoichio.iteritems():
            for k,v in d.iteritems():
                self.assertEquals(v, 1)



    def testReduce(self):
        import pickle
        reaction = pickle.loads(pickle.dumps(self.rrxn))

if __name__ == '__main__':
    unittest.main()        