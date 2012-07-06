#!/usr/bin/env python
# encoding: utf-8 -*-

import numpy
import unittest
import os

from rmgpy.cantherm.gaussian import *
from rmgpy.statmech import *

################################################################################

class GaussianTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.io.gaussian module, used for reading
    and writing Gaussian files.
    """
    
    def testLoadEthyleneFromGaussianLog(self):
        """
        Uses a Gaussian03 log file for ethylene (C2H4) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'ethylene.log'))
        s = log.loadStates()
        E0 = log.loadEnergy()
        
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,Translation)]) == 1)
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,RigidRotor)]) == 1)
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HinderedRotor)]) == 0)

        trans = [mode for mode in s.modes if isinstance(mode,Translation)][0]
        rot = [mode for mode in s.modes if isinstance(mode,RigidRotor)][0]
        vib = [mode for mode in s.modes if isinstance(mode,HarmonicOscillator)][0]
        Tlist = numpy.array([298.15], numpy.float64)
        self.assertAlmostEqual(trans.getPartitionFunction(Tlist) / 5.83338e6, 1.0, 3)
        self.assertAlmostEqual(rot.getPartitionFunction(Tlist) / 2.59622e3, 1.0, 3)
        self.assertAlmostEqual(vib.getPartitionFunction(Tlist) / 1.0481e0, 1.0, 3)

        self.assertAlmostEqual(E0 / 6.02214179e23 / 4.35974394e-18 / -78.563169, 1.0, 2)
        self.assertEqual(s.spinMultiplicity, 1)

    def testLoadOxygenFromGaussianLog(self):
        """
        Uses a Gaussian03 log file for oxygen (O2) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__),'oxygen.log'))
        s = log.loadStates()
        E0 = log.loadEnergy()
        
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,Translation)]) == 1)
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,RigidRotor)]) == 1)
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HinderedRotor)]) == 0)

        trans = [mode for mode in s.modes if isinstance(mode,Translation)][0]
        rot = [mode for mode in s.modes if isinstance(mode,RigidRotor)][0]
        vib = [mode for mode in s.modes if isinstance(mode,HarmonicOscillator)][0]
        Tlist = numpy.array([298.15], numpy.float64)
        self.assertAlmostEqual(trans.getPartitionFunction(Tlist) / 7.11169e6, 1.0, 3)
        self.assertAlmostEqual(rot.getPartitionFunction(Tlist) / 7.13316e1, 1.0, 3)
        self.assertAlmostEqual(vib.getPartitionFunction(Tlist) / 1.000037e0, 1.0, 3)

        self.assertAlmostEqual(E0 / 6.02214179e23 / 4.35974394e-18 / -150.374756, 1.0, 4)
        self.assertEqual(s.spinMultiplicity, 3)

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
