#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest
import sys
sys.path.append('.')

from chempy.species import Species, TransitionState
from chempy.reaction import *
from chempy.states import *
from chempy.kinetics import ArrheniusModel

################################################################################

class ReactionTest(unittest.TestCase):
	"""
	Contains unit tests for the chempy.reaction module, used for working with
	chemical reaction objects.
	"""
	
	def testTSTCalculation(self):
		"""
		A test of the transition state theory k(T) calculation function,
		using the reaction H + C2H4 -> C2H5.
		"""
		
		states = StatesModel(
			modes = [Translation(mass=0.0280313, volume=1, dimension=3), RigidRotor(linear=False, inertia=[5.69516e-47, 2.77584e-46, 3.34536e-46], symmetry=4), HarmonicOscillator(frequencies=[834.499, 973.312, 975.369, 1067.13, 1238.46, 1379.46, 1472.29, 1691.34, 3121.57, 3136.7, 3192.46, 3220.98])],
			E0=-205882860.949,
			spinMultiplicity=1,
		)
		ethylene = Species(states=states)
		
		states = StatesModel(
			modes = [Translation(mass=0.00100783, volume=1, dimension=3), HarmonicOscillator(frequencies=[])],
			E0=-1318675.56138,
			spinMultiplicity=2,
		)
		hydrogen = Species(states=states)
		
		states = StatesModel(
			modes = [Translation(mass=0.0290391, volume=1, dimension=3), RigidRotor(linear=False, inertia=[8.07491e-47, 3.69475e-46, 3.9885e-46], symmetry=1), HarmonicOscillator(frequencies=[466.816, 815.399, 974.674, 1061.98, 1190.71, 1402.03, 1467, 1472.46, 1490.98, 2972.34, 2994.88, 3089.96, 3141.01, 3241.96])],
			E0=-207340036.867,
			spinMultiplicity=2,
		)
		ethyl = Species(states=states)
		
		states = StatesModel(
			modes = [Translation(mass=0.0290391, volume=1, dimension=3), RigidRotor(linear=False, inertia=[1.2553e-46, 3.68827e-46, 3.80416e-46], symmetry=2), HarmonicOscillator(frequencies=[241.47, 272.706, 833.984, 961.614, 974.994, 1052.32, 1238.23, 1364.42, 1471.38, 1655.51, 3128.29, 3140.3, 3201.94, 3229.51])],
			E0=-207188826.467,
			spinMultiplicity=2,
		)
		TS = TransitionState(states=states, frequency=-309.3437)
		
		reaction = Reaction(reactants=[hydrogen, ethylene], products=[ethyl])
	
		import numpy
		Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
		klist = reaction.calculateTSTRateCoefficient(Tlist, TS, tunneling='')
		arrhenius = ArrheniusModel().fitToData(Tlist, klist)
		klist2 = arrhenius.getRateCoefficient(Tlist)
		
		# Check that the correct Arrhenius parameters are returned
		self.assertAlmostEqual(arrhenius.A/458.87, 1.0, 2)
		self.assertAlmostEqual(arrhenius.n/0.978, 1.0, 2)
		self.assertAlmostEqual(arrhenius.Ea/10194, 1.0, 2)
		# Check that the fit is satisfactory
		for i in range(len(Tlist)):
			self.assertTrue(abs(1 - klist2[i] / klist[i]) < 0.01)
		
		
if __name__ == '__main__':
	unittest.main()
