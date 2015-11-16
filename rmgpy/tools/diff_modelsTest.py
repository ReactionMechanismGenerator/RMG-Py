import unittest
import os
import os.path
import shutil

from rmgpy.tools.diff_models import *

class DiffModelsTest(unittest.TestCase):

    def test_identical_models(self):
        folder = os.path.join(os.getcwd(),'rmgpy/tools/data/diffmodels')
        
        chemkin1 = os.path.join(folder,'chem1.inp')
        speciesDict1 = os.path.join(folder,'species_dictionary1.txt')

        chemkin2 = os.path.join(folder,'chem2.inp')
        speciesDict2 = os.path.join(folder,'species_dictionary2.txt')
        
        kwargs = {
            'wd': folder,
        }

        execute(chemkin1, speciesDict1, None, chemkin2, speciesDict2, None, **kwargs)

        shutil.rmtree(os.path.join(folder,'species1'))
        shutil.rmtree(os.path.join(folder,'species2'))
        os.remove(os.path.join(folder,'diff.html'))