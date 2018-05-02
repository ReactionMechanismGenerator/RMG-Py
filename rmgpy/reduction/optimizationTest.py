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

import os
import unittest
from external.wip import work_in_progress

from rmgpy.scoop_framework.framework import TestScoopCommon
from rmgpy.scoop_framework.util import logger as logging

import rmgpy

try:
    from scoop import futures
except ImportError, e:
    logging.debug("Could not properly import SCOOP.")

from .input import load
from .reduction import initialize, computeObservables

from .optimization import *

def funcOptimize(rmg, targets):
    reactionModel = rmg.reactionModel

    initialize(rmg.outputDirectory, reactionModel.core.reactions)

    error = OptimizeTest.error

    index = 0
    reactionSystem = rmg.reactionSystems[index]

    #compute original target observables
    observables = computeObservables(targets, reactionModel, reactionSystem,
                                     rmg.absoluteTolerance, rmg.relativeTolerance)

    # optimize reduction tolerance
    tol, importantRxns = optimize(targets, reactionModel, rmg, index, error, observables)

    try:
        assert len(importantRxns) == 30
    except AssertionError:
        return False

    return True
    
@work_in_progress
class OptimizeTest(TestScoopCommon):

    #MINIMAL
    wd = os.path.join(os.path.dirname(rmgpy.__file__),'reduction/test_data/minimal/')
    inputFile = os.path.join(wd, 'input.py')
    reductionFile = os.path.join(wd, 'reduction_input.py')
    chemkinFile = os.path.join(wd, 'chemkin','chem_annotated.inp')
    spcDict = os.path.join(wd, 'chemkin','species_dictionary.txt')

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    @classmethod
    def setUpClass(cls):
        super(OptimizeTest, cls).setUpClass()
        rmg, targets, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spcDict)
        cls.rmg = rmg
        cls.targets = targets
        cls.error = error


    def testOptimize(self):
        rmg = OptimizeTest.rmg
        targets = OptimizeTest.targets
        
        result = futures._startup(funcOptimize, rmg, targets)
        self.assertEquals(result, True)  

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main() 
