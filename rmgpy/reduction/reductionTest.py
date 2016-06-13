
import os.path
import unittest

import rmgpy

from .reduction import *

class ReduceTest(unittest.TestCase):

    #MINIMAL
    wd = os.path.join(os.path.dirname(rmgpy.__file__), 'reduction/test_data/minimal/')
    inputFile = os.path.join(wd, 'input.py')
    reductionFile = os.path.join(wd, 'reduction_input.py')
    chemkinFile = os.path.join(wd, 'chemkin','chem_annotated.inp')
    spcDict = os.path.join(wd, 'chemkin','species_dictionary.txt')


    @classmethod
    def setUpClass(cls):
        from .input import load

        super(ReduceTest, cls).setUpClass()

        rmg, targets, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spcDict)
        cls.rmg = rmg
        cls.targets = targets
        cls.error = error

        reactionModel = rmg.reactionModel
        initialize(rmg.outputDirectory, reactionModel.core.reactions)
    

    def testComputeConversion(self):
        rmg = ReduceTest.rmg
        target = ReduceTest.targets[0]
        reactionModel = rmg.reactionModel

        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]

        conv = computeConversion(target, reactionModel, reactionSystem,\
         rmg.absoluteTolerance, rmg.relativeTolerance)
        self.assertIsNotNone(conv)


    def testReduceCompute(self):
        rmg = ReduceTest.rmg
        targets = ReduceTest.targets
        reactionModel = rmg.reactionModel


        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]

        observables = computeObservables(targets, reactionModel, reactionSystem, \
            rmg.absoluteTolerance, rmg.relativeTolerance)

        tols = [0.7, 1e-3, 1e-6]
        for tol in tols:
            conv, importantRxns = reduceModel(tol, targets, reactionModel, rmg, index)
            self.assertIsNotNone(conv)

if __name__ == '__main__':
    unittest.main()