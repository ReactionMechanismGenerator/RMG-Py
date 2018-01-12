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
