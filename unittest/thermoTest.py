#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest
import sys
sys.path.append('.')

import chempy.constants as constants
from chempy.thermo import *

################################################################################

class ThermoTest(unittest.TestCase):
	"""
	Contains unit tests for the chempy.thermo module, used for working with
	thermodynamics models.
	"""
	
	def testWilhoit(self):
		"""
		Tests the Wilhoit thermodynamics model functions.
		"""
		
		# CC(=O)O[O]
		wilhoit = WilhoitModel(cp0=4.0*constants.R, cpInf=21.0*constants.R, a0=-3.95, a1=9.26, a2=-15.6, a3=8.55, B=500.0, H0=-6.151e+04, S0=-790.2)
		
		Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
		Cplist0 = [ 64.398,  94.765, 116.464, 131.392, 141.658, 148.830, 153.948, 157.683, 160.469, 162.589]
		Hlist0 = [-166312., -150244., -128990., -104110., -76742.9, -47652.6, -17347.1, 13834.8, 45663.0, 77978.1]
		Slist0 = [287.421, 341.892, 384.685, 420.369, 450.861, 477.360, 500.708, 521.521, 540.262, 557.284]
		Glist0 = [-223797., -287002., -359801., -440406., -527604., -620485., -718338., -820599., -926809., -1036590.]

		Cplist = wilhoit.getHeatCapacity(Tlist)
		Hlist = wilhoit.getEnthalpy(Tlist)
		Slist = wilhoit.getEntropy(Tlist)
		Glist = wilhoit.getFreeEnergy(Tlist)

		for i in range(len(Tlist)):
			self.assertAlmostEqual(Cplist[i] / Cplist0[i], 1.0, 4)
			self.assertAlmostEqual( Hlist[i] /  Hlist0[i], 1.0, 4)
			self.assertAlmostEqual( Slist[i] /  Slist0[i], 1.0, 4)
			self.assertAlmostEqual( Glist[i] /  Glist0[i], 1.0, 4)

if __name__ == '__main__':
	unittest.main()
