
import os.path
import unittest

from .reduction import *

class ReduceTest(unittest.TestCase):

    #MINIMAL
    wd = os.path.join('rmgpy/reduction/test_data/minimal/')
    inputFile = os.path.join(wd, 'input.py')
    reductionFile = os.path.join(wd, 'reduction_input.py')
    chemkinFile = os.path.join(wd, 'chemkin','chem.inp')
    spc_dict = os.path.join(wd, 'chemkin','species_dictionary.txt')


    @classmethod
    def setUpClass(cls):
        from .input import load

        super(ReduceTest, cls).setUpClass()

        rmg, targets, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spc_dict)
        cls.rmg = rmg
        cls.targets = targets
        cls.error = error

        reactionModel = rmg.reactionModel
        initialize(rmg.outputDirectory, reactionModel.core.reactions)
    

    def test_compute_conversion(self):
        rmg = ReduceTest.rmg
        target = ReduceTest.targets[0]
        reactionModel = rmg.reactionModel

        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]

        conv = compute_conversion(target, reactionModel, reactionSystem,\
         rmg.absoluteTolerance, rmg.relativeTolerance)
        self.assertIsNotNone(conv)


    def test_reduce_compute(self):
        rmg = ReduceTest.rmg
        targets = ReduceTest.targets
        reactionModel = rmg.reactionModel


        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]

        observables = compute_observables(targets, reactionModel, reactionSystem, \
            rmg.absoluteTolerance, rmg.relativeTolerance)

        tols = [0.7, 1e-3, 1e-6]
        for tol in tols:
            conv, important_rxns = reduce_model(tol, targets, reactionModel, rmg, index)
            self.assertIsNotNone(conv)

if __name__ == '__main__':
    unittest.main()