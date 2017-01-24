
import os.path
import unittest

import rmgpy
from external.wip import work_in_progress

from .reduction import *

from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

class ReduceFunctionalTest(unittest.TestCase):

    #MINIMAL
    wd = os.path.join(os.path.dirname(rmgpy.__file__), 'reduction/test_data/minimal/')
    inputFile = os.path.join(wd, 'input.py')
    reductionFile = os.path.join(wd, 'reduction_input.py')
    chemkinFile = os.path.join(wd, 'chemkin','chem_annotated.inp')
    spcDict = os.path.join(wd, 'chemkin','species_dictionary.txt')


    @classmethod
    def setUpClass(cls):
        from .input import load

        super(ReduceFunctionalTest, cls).setUpClass()

        rmg, targets, error = load(cls.inputFile, cls.reductionFile, cls.chemkinFile, cls.spcDict)
        cls.rmg = rmg
        cls.targets = targets
        cls.error = error

        reactionModel = rmg.reactionModel
        initialize(rmg.outputDirectory, reactionModel.core.reactions)
    

    def testComputeConversion(self):
        rmg = ReduceFunctionalTest.rmg
        target = ReduceFunctionalTest.targets[0]
        reactionModel = rmg.reactionModel

        atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
        index = 0
        reactionSystem = rmg.reactionSystems[index]

        conv = computeConversion(target, reactionModel, reactionSystem,\
         rmg.absoluteTolerance, rmg.relativeTolerance)
        self.assertIsNotNone(conv)


    def testReduceCompute(self):
        rmg = ReduceFunctionalTest.rmg
        targets = ReduceFunctionalTest.targets
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

class ReduceUnitTest(unittest.TestCase):
    

    @work_in_progress
    def testAllEntriesAccessibleInSearchTargetIndex(self):
        butene1 = Species()
        butene1.fromSMILES('C=CCC')
        butene1.label = 'C4H8'
        
        butene2 = Species()
        butene2.fromSMILES('CC=CC')
        butene2.label = 'C4H8'
        
        
        species_list =[butene1,butene2]
        # make sure different species with same label
        assert not species_list[0].isIsomorphic(species_list[1])
        assert species_list[0].label == species_list[1].label
        
        # make fake reactionModel object to fit in with the unittest
        reaction_model = CoreEdgeReactionModel()
        reaction_model.core.species = species_list
        
        # ensure second species index is returned when it's label is used
        # in `searchTargetIndex`.
        input_index = 1
        output_index = searchTargetIndex(species_list[input_index].label,reaction_model)
        self.assertEqual(input_index,output_index,'searchTargetIndex will not return the second occurance of species with the same label.')
if __name__ == '__main__':
    unittest.main()