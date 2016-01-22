import unittest
import os.path
import shutil

import rmgpy
from rmgpy.rmg.main import RMG
from rmgpy.tools.generate_reactions import *

class GenerateReactionsTest(unittest.TestCase):

    def test(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/generate')
        
        inputFile = os.path.join(folder,'input.py')
        
        rmg = RMG()
        rmg = execute(rmg, inputFile, folder)

        self.assertIsNotNone(rmg)
        self.assertIsNotNone(rmg.reactionModel.outputSpeciesList)
        self.assertIsNotNone(rmg.reactionModel.outputReactionList)


        shutil.rmtree(os.path.join(folder,'pdep'))

    def testDuplicateReaction(self):
        """
        Test that the radical addition reaction

        HCJ=O + CH2O = [CH2]OC=O

        present in the reaction library "Methylformate",
        only appears once in the model.

        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/generate/duplicates')
        
        inputFile = os.path.join(folder,'input.py')
        
        rmg = RMG()
        rmg = execute(rmg, inputFile, folder)

        self.assertIsNotNone(rmg)

        """
        13 reactions will be created,12 from the seed mechanism,
        1 additional one through the reaction family.
        """
        self.assertEquals(13, len(rmg.reactionModel.core.reactions))

        shutil.rmtree(os.path.join(folder,'pdep'))

    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
