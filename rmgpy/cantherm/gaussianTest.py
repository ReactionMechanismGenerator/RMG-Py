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

from rmgpy.cantherm.gaussian import GaussianLog
from rmgpy.statmech import IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator, HinderedRotor
import rmgpy.constants as constants
from external.wip import work_in_progress
################################################################################

class GaussianTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.io.gaussian module, used for reading
    and writing Gaussian files.
    """
    @work_in_progress
    def testLoadEthyleneFromGaussianLog_CBSQB3(self):
        """
        Uses a Gaussian03 log file for ethylene (C2H4) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'data','ethylene.log'))
        conformer = log.loadConformer()
        E0 = log.loadEnergy()
        
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,NonlinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HinderedRotor)]) == 0)

        trans = [mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)][0]
        rot = [mode for mode in conformer.modes if isinstance(mode,NonlinearRotor)][0]
        vib = [mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)][0]
        Tlist = numpy.array([298.15], numpy.float64)
        self.assertAlmostEqual(trans.getPartitionFunction(Tlist), 5.83338e6, delta=1e1)
        self.assertAlmostEqual(rot.getPartitionFunction(Tlist), 2.59622e3, delta=1e-2)
        self.assertAlmostEqual(vib.getPartitionFunction(Tlist), 1.0481e0, delta=1e-4)

        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -78.467452, 4)
        self.assertEqual(conformer.spinMultiplicity, 1)
        self.assertEqual(conformer.opticalIsomers, 1)

    def testLoadOxygenFromGaussianLog(self):
        """
        Uses a Gaussian03 log file for oxygen (O2) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'data','oxygen.log'))
        conformer = log.loadConformer()
        E0 = log.loadEnergy()
        
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,LinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HinderedRotor)]) == 0)

        trans = [mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)][0]
        rot = [mode for mode in conformer.modes if isinstance(mode,LinearRotor)][0]
        vib = [mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)][0]
        Tlist = numpy.array([298.15], numpy.float64)
        self.assertAlmostEqual(trans.getPartitionFunction(Tlist), 7.11169e6, delta=1e1)
        self.assertAlmostEqual(rot.getPartitionFunction(Tlist), 7.13316e1, delta=1e-4)
        self.assertAlmostEqual(vib.getPartitionFunction(Tlist), 1.00037e0, delta=1e-4)
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -150.3784877, 4)
        self.assertEqual(conformer.spinMultiplicity, 3)
        self.assertEqual(conformer.opticalIsomers, 1)

    @work_in_progress
    def testLoadEthyleneFromGaussianLog_G3(self):
        """
        Uses a Gaussian03 log file for ethylene (C2H4) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'data','ethylene_G3.log'))
        conformer = log.loadConformer()
        E0 = log.loadEnergy()
        
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,NonlinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode,HinderedRotor)]) == 0)

        trans = [mode for mode in conformer.modes if isinstance(mode,IdealGasTranslation)][0]
        rot = [mode for mode in conformer.modes if isinstance(mode,NonlinearRotor)][0]
        vib = [mode for mode in conformer.modes if isinstance(mode,HarmonicOscillator)][0]
        Tlist = numpy.array([298.15], numpy.float64)
        
        self.assertAlmostEqual(trans.getPartitionFunction(Tlist), 5.83338e6, delta=1e1)
        self.assertAlmostEqual(rot.getPartitionFunction(Tlist), 2.53410e3, delta=1e-2)
        self.assertAlmostEqual(vib.getPartitionFunction(Tlist), 1.0304e0, delta=1e-4)

        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -78.562189, 4)
        self.assertEqual(conformer.spinMultiplicity, 1)
        self.assertEqual(conformer.opticalIsomers, 1)

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
