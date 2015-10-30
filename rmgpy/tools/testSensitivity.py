import unittest
import os
import os.path
import shutil

from rmgpy.tools.sensitivity import *

class SensitivityTest(unittest.TestCase):

    def test_minimal(self):
        folder = os.path.join(os.getcwd(),'rmgpy/tools/data')
        
        inputFile = os.path.join(folder,'input.py')
        chemkinFile = os.path.join(folder,'chem.inp')
        dictFile = os.path.join(folder,'species_dictionary.txt')
        
        runSensitivity(inputFile, chemkinFile, dictFile)

        outputdir = os.path.join(folder, 'species/')
        simfile = os.path.join(folder,'simulation_1.csv')
        sensfile = os.path.join(folder,'sensitivity_1_SPC_1.csv')

        self.assertTrue(os.path.isfile(simfile))
        self.assertTrue(os.path.isfile(sensfile))
        
        shutil.rmtree(outputdir)
        os.remove(simfile)
        os.remove(sensfile)