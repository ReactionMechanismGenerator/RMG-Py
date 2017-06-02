import unittest
import os
import os.path
import shutil

import rmgpy
from rmgpy.tools.fluxdiagram import *

class FluxDiagramTest(unittest.TestCase):

    def test_avi(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools','data','flux')
        
        inputFile = os.path.join(folder,'input.py')
        run(inputFile)

        outputdir = os.path.join(folder,'flux')
        simfile = os.path.join(outputdir,'1','flux_diagram.avi')

        speciesdir = os.path.join(folder,'species')

        self.assertTrue(os.path.isfile(simfile))
        
        shutil.rmtree(outputdir)
        shutil.rmtree(speciesdir)

    def tearDown(self):
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None