#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest
import sys
sys.path.append('.')

from chempy.states import *

################################################################################

class StatesTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.states module, used for working with
    molecular degrees of freedom.
    """
    
    def testModesForEthylene(self):
        """
        Uses data for ethylene (C2H4) to test the various modes. The data comes
        from a CBS-QB3 calculation using Gaussian03.
        """

        Tlist = numpy.array([298.15], numpy.float64)

        trans = Translation(mass=0.02803, volume=1.0, dimension=3)
        rot = RigidRotor(linear=False, inertia=[5.6952e-47, 2.7758e-46, 3.3454e-46], symmetry=1)
        vib = HarmonicOscillator(frequencies=[834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3, 3121.6, 3136.7, 3192.5, 3221.0])

        self.assertAlmostEqual(trans.getPartitionFunction(Tlist) * 1.3806504e-23 * 298.15 / 101325 / 5.83338e6, 1.0, 3)
        self.assertAlmostEqual(rot.getPartitionFunction(Tlist) / 2.59622e3, 1.0, 3)
        self.assertAlmostEqual(vib.getPartitionFunction(Tlist) / 1.0481e0, 1.0, 3)

        self.assertAlmostEqual(trans.getHeatCapacity(Tlist) / 4.184 / 2.981, 1.0, 3)
        self.assertAlmostEqual(rot.getHeatCapacity(Tlist) / 4.184 / 2.981, 1.0, 3)
        self.assertAlmostEqual(vib.getHeatCapacity(Tlist) / 4.184 / 2.133, 1.0, 3)

        self.assertAlmostEqual(trans.getEnthalpy(Tlist) / 8.314472 / Tlist / 1.5, 1.0, 3)
        self.assertAlmostEqual(rot.getEnthalpy(Tlist) / 8.314472 / Tlist / 1.5, 1.0, 3)
        self.assertAlmostEqual(vib.getEnthalpy(Tlist) / 8.314472 / Tlist / 0.221258, 1.0, 3)

        self.assertAlmostEqual(trans.getEntropy(Tlist) / 4.184 / 152.094, 1.0, 3)
        self.assertAlmostEqual(rot.getEntropy(Tlist) / 4.184 / 18.604, 1.0, 3)
        self.assertAlmostEqual(vib.getEntropy(Tlist) / 4.184 / 0.533, 1.0, 3)

        states = StatesModel(modes=[rot, vib], spinMultiplicity=1)
        
        dE = 10.0
        Elist = numpy.arange(0, 100001, dE, numpy.float64)
        rho = states.getDensityOfStates(Elist)
        self.assertAlmostEqual(numpy.sum(rho * numpy.exp(-Elist / 8.314472 / 298.15) * dE) / states.getPartitionFunction(Tlist), 1.0, 2)
        
    def testModesForOxygen(self):
        """
        Uses data for oxygen (O2) to test the various modes. The data comes
        from a CBS-QB3 calculation using Gaussian03.
        """

        Tlist = numpy.array([298.15], numpy.float64)

        trans = Translation(mass=0.03199, volume=1.0, dimension=3)
        rot = RigidRotor(linear=True, inertia=[1.9271e-46], symmetry=2)
        vib = HarmonicOscillator(frequencies=[1637.9])

        self.assertAlmostEqual(trans.getPartitionFunction(Tlist) * 1.3806504e-23 * 298.15 / 101325 / 7.11169e6, 1.0, 3)
        self.assertAlmostEqual(rot.getPartitionFunction(Tlist) / 7.13316e1, 1.0, 3)
        self.assertAlmostEqual(vib.getPartitionFunction(Tlist) / 1.000037e0, 1.0, 3)

        self.assertAlmostEqual(trans.getHeatCapacity(Tlist) / 4.184 / 2.981, 1.0, 3)
        self.assertAlmostEqual(rot.getHeatCapacity(Tlist) / 4.184 / 1.987, 1.0, 3)
        self.assertAlmostEqual(vib.getHeatCapacity(Tlist) / 4.184 / 0.046, 1.0, 2)

        self.assertAlmostEqual(trans.getEnthalpy(Tlist) / 8.314472 / Tlist / 1.5, 1.0, 3)
        self.assertAlmostEqual(rot.getEnthalpy(Tlist) / 8.314472 / Tlist / 1.0, 1.0, 3)
        self.assertAlmostEqual(vib.getEnthalpy(Tlist) / 8.314472 / Tlist / 0.0029199, 1.0, 3)

        self.assertAlmostEqual(trans.getEntropy(Tlist) / 4.184 / 152.488, 1.0, 3)
        self.assertAlmostEqual(rot.getEntropy(Tlist) / 4.184 / 10.467, 1.0, 3)
        self.assertAlmostEqual(vib.getEntropy(Tlist) / 4.184 / 0.00654, 1.0, 2)

        states = StatesModel(modes=[rot, vib], spinMultiplicity=3)
        
        dE = 10.0
        Elist = numpy.arange(0, 100001, dE, numpy.float64)
        rho = states.getDensityOfStates(Elist)
        self.assertAlmostEqual(numpy.sum(rho * numpy.exp(-Elist / 8.314472 / 298.15) * dE) / states.getPartitionFunction(Tlist), 1.0, 2)
    
    def testHinderedRotorDensityOfStates(self):
        """
        Test that the density of states and the partition function of the
        hindered rotor are self-consistent. This is turned off because the
        density of states is for the classical limit only, while the partition
        function is not.
        """

        hr = HinderedRotor(inertia=3e-46, barrier=0.5*4184, symmetry=3)
        dE = 10.0
        Elist = numpy.arange(0, 100001, dE, numpy.float64)
        rho = hr.getDensityOfStates(Elist)

#        Tlist = 1000.0 / numpy.arange(0.5, 3.5, 0.1, numpy.float64)
#        Q = numpy.zeros_like(Tlist)
#        for i in range(len(Tlist)):
#            Q[i] = numpy.sum(rho * numpy.exp(-Elist / 8.314472 / Tlist[i]) * dE)
#        import pylab
#        pylab.semilogy(1000.0 / Tlist, Q, '--k', 1000.0 / Tlist, hr.getPartitionFunction(Tlist), '-k')
#        pylab.show()

        Tlist = numpy.array([298.15], numpy.float64)
        self.assertTrue(0.9 < numpy.sum(rho * numpy.exp(-Elist / 8.314472 / Tlist[0]) * dE) / hr.getPartitionFunction(Tlist) < 1.1)
        Tlist = numpy.array([1000.0], numpy.float64)
        self.assertTrue(0.9 < numpy.sum(rho * numpy.exp(-Elist / 8.314472 / Tlist[0]) * dE) / hr.getPartitionFunction(Tlist) < 1.1)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
