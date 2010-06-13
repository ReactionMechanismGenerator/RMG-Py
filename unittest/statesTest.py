#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest
import sys
sys.path.append('.')

from chempy.states import *

################################################################################

class StatesTest(unittest.TestCase):

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

		self.assertAlmostEqual(trans.getHeatCapacity(Tlist) * 1.987 / 2.981, 1.0, 3)
		self.assertAlmostEqual(rot.getHeatCapacity(Tlist) * 1.987 / 2.981, 1.0, 3)
		self.assertAlmostEqual(vib.getHeatCapacity(Tlist) * 1.987 / 2.133, 1.0, 3)

		self.assertAlmostEqual(trans.getEnthalpy(Tlist) / 1.5, 1.0, 3)
		self.assertAlmostEqual(rot.getEnthalpy(Tlist) / 1.5, 1.0, 3)
		self.assertAlmostEqual(vib.getEnthalpy(Tlist) / 0.221258, 1.0, 3)

		self.assertAlmostEqual(trans.getEntropy(Tlist) * 1.987 / 152.094, 1.0, 3)
		self.assertAlmostEqual(rot.getEntropy(Tlist) * 1.987 / 18.604, 1.0, 3)
		self.assertAlmostEqual(vib.getEntropy(Tlist) * 1.987 / 0.533, 1.0, 3)

		states = StatesModel()
		states.modes = [rot, vib]

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

		self.assertAlmostEqual(trans.getHeatCapacity(Tlist) * 1.987 / 2.981, 1.0, 3)
		self.assertAlmostEqual(rot.getHeatCapacity(Tlist) * 1.987 / 1.987, 1.0, 3)
		self.assertAlmostEqual(vib.getHeatCapacity(Tlist) * 1.987 / 0.046, 1.0, 2)

		self.assertAlmostEqual(trans.getEnthalpy(Tlist) / 1.5, 1.0, 3)
		self.assertAlmostEqual(rot.getEnthalpy(Tlist) / 1.0, 1.0, 3)
		self.assertAlmostEqual(vib.getEnthalpy(Tlist) / 0.0029199, 1.0, 3)

		self.assertAlmostEqual(trans.getEntropy(Tlist) * 1.987 / 152.488, 1.0, 3)
		self.assertAlmostEqual(rot.getEntropy(Tlist) * 1.987 / 10.467, 1.0, 3)
		self.assertAlmostEqual(vib.getEntropy(Tlist) * 1.987 / 0.00654, 1.0, 2)

		states = StatesModel()
		states.modes = [rot, vib]

		dE = 10.0
		Elist = numpy.arange(0, 100001, dE, numpy.float64)
		rho = states.getDensityOfStates(Elist)
		self.assertAlmostEqual(numpy.sum(rho * numpy.exp(-Elist / 8.314472 / 298.15) * dE) / states.getPartitionFunction(Tlist), 1.0, 2)

	def testLoadEthyleneFromGaussianLog(self):
		"""
		Uses a Gaussian03 log file for ethylene (C2H4) to test that its
		molecular degrees of freedom can be properly read.
		"""

		s = StatesModel()
		s.loadFromGaussianLog('unittest/ethylene.log')
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,Translation)]) == 1)
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,RigidRotor)]) == 1)
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HarmonicOscillator)]) == 1)
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HinderedRotor)]) == 0)

		trans = [mode for mode in s.modes if isinstance(mode,Translation)][0]
		rot = [mode for mode in s.modes if isinstance(mode,RigidRotor)][0]
		vib = [mode for mode in s.modes if isinstance(mode,HarmonicOscillator)][0]
		Tlist = numpy.array([298.15], numpy.float64)
		self.assertAlmostEqual(trans.getPartitionFunction(Tlist) * 1.3806504e-23 * 298.15 / 101325 / 5.83338e6, 1.0, 3)
		self.assertAlmostEqual(rot.getPartitionFunction(Tlist) / 2.59622e3, 1.0, 3)
		self.assertAlmostEqual(vib.getPartitionFunction(Tlist) / 1.0481e0, 1.0, 3)

		self.assertAlmostEqual(s.E0 / 6.02214179e23 / 4.35974394e-18 / -78.563169, 1.0, 2)
		self.assertEqual(s.spinMultiplicity, 1)

	def testLoadOxygenFromGaussianLog(self):
		"""
		Uses a Gaussian03 log file for oxygen (O2) to test that its
		molecular degrees of freedom can be properly read.
		"""

		s = StatesModel()
		s.loadFromGaussianLog('unittest/oxygen.log')
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,Translation)]) == 1)
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,RigidRotor)]) == 1)
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HarmonicOscillator)]) == 1)
		self.assertTrue(len([mode for mode in s.modes if isinstance(mode,HinderedRotor)]) == 0)

		trans = [mode for mode in s.modes if isinstance(mode,Translation)][0]
		rot = [mode for mode in s.modes if isinstance(mode,RigidRotor)][0]
		vib = [mode for mode in s.modes if isinstance(mode,HarmonicOscillator)][0]
		Tlist = numpy.array([298.15], numpy.float64)
		self.assertAlmostEqual(trans.getPartitionFunction(Tlist) * 1.3806504e-23 * 298.15 / 101325 / 7.11169e6, 1.0, 3)
		self.assertAlmostEqual(rot.getPartitionFunction(Tlist) / 7.13316e1, 1.0, 3)
		self.assertAlmostEqual(vib.getPartitionFunction(Tlist) / 1.000037e0, 1.0, 3)

		self.assertAlmostEqual(s.E0 / 6.02214179e23 / 4.35974394e-18 / -150.374756, 1.0, 4)
		self.assertEqual(s.spinMultiplicity, 2)

if __name__ == '__main__':
	unittest.main()
