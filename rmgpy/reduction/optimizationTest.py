import os
import sys
import unittest
import numpy as np

from rmgpy.scoop_framework.framework import TestScoopCommon
from rmgpy.scoop_framework.util import logger as logging

try:
    from scoop import futures, _control, shared
except ImportError, e:
    logging.debug("Could not properly import SCOOP.")

from .input import load
from .reduction import initialize, compute_observables

from .optimization import *

def funcOptimize(rmg, targets):
    reactionModel = rmg.reactionModel

    initialize(rmg.outputDirectory, reactionModel.core.reactions)

    error = OptimizeTest.error

    atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
    index = 0
    reactionSystem = rmg.reactionSystems[index]

    #compute original target observables
    observables = compute_observables(targets, reactionModel, reactionSystem, \
     rmg.absoluteTolerance, rmg.relativeTolerance)

    # optimize reduction tolerance
    tol, important_rxns = optimize(targets, reactionModel, rmg, index, error, observables)

    try:
        assert len(important_rxns) == 34
    except AssertionError:
        return False

    return True
    

class OptimizeTest(TestScoopCommon):

    #MINIMAL
    wd = os.path.join('rmgpy/reduction/test_data/minimal/')
    inputFile = os.path.join(wd, 'input.py')
    reductionFile = os.path.join(wd, 'reduction_input.py')
    chemkinFile = os.path.join(wd, 'chemkin','chem.inp')
    spc_dict = os.path.join(wd, 'chemkin','species_dictionary.txt')

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    @classmethod
    def setUpClass(cls):
        super(OptimizeTest, cls).setUpClass()
        rmg, targets, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spc_dict)
        cls.rmg = rmg
        cls.targets = targets
        cls.error = error


    def test_optimize(self):
        rmg = OptimizeTest.rmg
        targets = OptimizeTest.targets
        
        result = futures._startup(funcOptimize, rmg, targets)
        self.assertEquals(result, True)  

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main() 
