#!/usr/bin/env python
# encoding: utf-8 -*-

import numpy
import unittest
import os

from rmgpy.cantherm.gaussian import GaussianLog
from rmgpy.statmech import Conformer, IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator, HinderedRotor
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

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'files','ethylene.log'))
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

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'files','oxygen.log'))
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

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'files','ethylene_G3.log'))
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
