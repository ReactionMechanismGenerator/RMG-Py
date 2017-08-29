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

import unittest
import os
import os.path
import shutil

from rmgpy.tools.sensitivity import runSensitivity
import rmgpy

class SensitivityTest(unittest.TestCase):

    def test_minimal(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'tools/data/sens/simple')
        
        inputFile = os.path.join(folder, 'input.py')
        chemkinFile = os.path.join(folder, 'chem.inp')
        dictFile = os.path.join(folder, 'species_dictionary.txt')
        
        runSensitivity(inputFile, chemkinFile, dictFile)

        simfile = os.path.join(folder, 'solver', 'simulation_1_13.csv')
        sensfile = os.path.join(folder, 'solver', 'sensitivity_1_SPC_1.csv')

        self.assertTrue(os.path.isfile(simfile))
        self.assertTrue(os.path.isfile(sensfile))
        
        shutil.rmtree(os.path.join(folder, 'solver'))

    def test_liquid(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'tools/data/sens/liquid')

        inputFile = os.path.join(folder, 'input.py')
        chemkinFile = os.path.join(folder, 'chem.inp')
        dictFile = os.path.join(folder, 'species_dictionary.txt')

        runSensitivity(inputFile, chemkinFile, dictFile, diffusionLimited=False)

        simfile = os.path.join(folder, 'solver', 'simulation_1_28.csv')
        sensfile = os.path.join(folder, 'solver', 'sensitivity_1_SPC_1.csv')

        self.assertTrue(os.path.isfile(simfile))
        self.assertTrue(os.path.isfile(sensfile))

        shutil.rmtree(os.path.join(folder, 'solver'))

    def tearDown(self):
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
