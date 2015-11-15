import unittest
import os.path
import shutil

from rmgpy.rmg.main import RMG
from rmgpy.tools.generate_reactions import *

class GenerateReactionsTest(unittest.TestCase):

    def test(self):
        folder = os.path.join(os.getcwd(),'rmgpy/tools/data/generate')
        
        inputFile = os.path.join(folder,'input.py')
        
        rmg = RMG()
        rmg = execute(rmg, inputFile, folder)

        self.assertIsNotNone(rmg)
        self.assertIsNotNone(rmg.reactionModel.outputSpeciesList)
        self.assertIsNotNone(rmg.reactionModel.outputReactionList)


        shutil.rmtree(os.path.join(folder,'pdep'))
        shutil.rmtree(os.path.join(folder,'plot'))
        shutil.rmtree(os.path.join(folder,'solver'))