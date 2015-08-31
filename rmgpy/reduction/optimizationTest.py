import os
import sys
import unittest
import numpy as np

from rmgpy.scoop_framework.framework import TestScoopCommon
from rmgpy.scoop_framework.util import logger as logging

try:
    from scoop import futures, _control, shared
except ImportError, e:
    import logging
    logging.debug("Could not properly import SCOOP.")

from .input import load
from .reduction import initialize, compute_conversion

from .optimization import *

def funcOptimize(rmg, target_label):
    reactionModel = rmg.reactionModel

    initialize(rmg.outputDirectory, reactionModel.core.reactions)

    error = OptimizeTest.error

    atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
    index = 0
    reactionSystem = rmg.reactionSystems[index]

    #compute original target conversion
    Xorig = compute_conversion(target_label, reactionModel, reactionSystem, index,\
     rmg.absoluteTolerance, rmg.relativeTolerance)

    # optimize reduction tolerance
    tol, important_rxns = optimize(target_label, reactionModel, rmg, index, error, Xorig)

    try:
        assert np.allclose([1e-06], [tol])
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
        rmg, target_label, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spc_dict)
        cls.rmg = rmg
        cls.target_label = target_label
        cls.error = error


    def test_optimize(self):
        rmg = OptimizeTest.rmg
        target_label = OptimizeTest.target_label
        
        result = futures._startup(funcOptimize, rmg, target_label)
        self.assertEquals(result, True)  

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main() 
