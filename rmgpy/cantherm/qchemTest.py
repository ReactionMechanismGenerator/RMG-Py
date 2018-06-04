#!/usr/bin/env python
# encoding: utf-8 -*-

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

import numpy
import unittest
import os

from rmgpy.cantherm.qchem import QchemLog
from rmgpy.statmech import Conformer, IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator, HinderedRotor
import rmgpy.constants as constants
from external.wip import work_in_progress
################################################################################

class QChemTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.io.qchem module, used for reading
    and writing Qchem files.
    """
    def testNumberOfAtomsFromQchemLog(self):
        """
        Uses a Qchem log files to test that 
        number of atoms can be properly read.
        """
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','npropyl.out'))        
        self.assertEqual(log.getNumberOfAtoms(), 10)
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','co.out'))        
        self.assertEqual(log.getNumberOfAtoms(), 2) 

    def testEnergyFromQchemLog(self):
        """
        Uses a Qchem log files to test that 
        molecular energies can be properly read.
        """        
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','npropyl.out'))      
        self.assertAlmostEqual(log.loadEnergy(), -310896203.5432524, 1e-5)
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','co.out'))        
        self.assertAlmostEqual(log.loadEnergy(), -297402545.0217114, 1e-5)   
        
    def testLoadVibrationsFromQchemLog(self):
        """
        Uses a Qchem log files to test that 
        molecular energies can be properly read.
        """        
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','npropyl.out'))    
        conformer = log.loadConformer()    
        self.assertEqual(len(conformer.modes[2]._frequencies.getValue()), 24)    
        self.assertEqual(conformer.modes[2]._frequencies.getValue()[5], 881.79)       
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','co.out'))        
        conformer = log.loadConformer() 
        self.assertEqual(len(conformer.modes[2]._frequencies.getValue()), 1)         
        self.assertEqual(conformer.modes[2]._frequencies.getValue(), 2253.16)    
                           
    def testLoadNpropylModesFromQchemLog(self):
        """
        Uses a Qchem log file for npropyl to test that its
        molecular modes can be properly read.
        """
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','npropyl.out'))
        conformer = log.loadConformer()

        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,NonlinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HinderedRotor)]) == 0)
        
    def testSpinMultiplicityFromQchemLog(self):
        """
        Uses a Qchem log file for npropyl to test that its
        molecular degrees of freedom can be properly read.
        """
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','npropyl.out'))
        conformer = log.loadConformer()
        self.assertEqual(conformer.spinMultiplicity, 2)
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','co.out'))
        conformer = log.loadConformer()
        self.assertEqual(conformer.spinMultiplicity, 1)
    
    def testLoadCOModesFromQchemLog(self):
        """
        Uses a Qchem log file for CO to test that its
        molecular degrees of freedom can be properly read.
        """
        log = QchemLog(os.path.join(os.path.dirname(__file__),'data','co.out'))
        conformer = log.loadConformer()
        E0 = log.loadEnergy()
        
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,LinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,NonlinearRotor)]) == 0)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HinderedRotor)]) == 0)

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
