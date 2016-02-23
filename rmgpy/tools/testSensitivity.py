import unittest
import os
import os.path
import shutil

from rmgpy.tools.sensitivity import *
import rmgpy

class SensitivityTest(unittest.TestCase):

    def test_minimal(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/sens')
        
        inputFile = os.path.join(folder,'input.py')
        chemkinFile = os.path.join(folder,'chem.inp')
        dictFile = os.path.join(folder,'species_dictionary.txt')
        
        runSensitivity(inputFile, chemkinFile, dictFile)

        simfile = os.path.join(folder,'solver', 'simulation_1_17.csv')
        sensfile = os.path.join(folder,'solver', 'sensitivity_1_SPC_1.csv')

        self.assertTrue(os.path.isfile(simfile))
        self.assertTrue(os.path.isfile(sensfile))
        
        os.remove(simfile)
        os.remove(sensfile)