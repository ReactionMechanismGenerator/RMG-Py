#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('.')

import math
		
from rmg.kinetics.model import *
import rmg.constants

################################################################################

class KineticsCheck(unittest.TestCase):                          

	def testKineticsModel(self):
		"""
		
		"""
		
		kinetics = KineticsModel(Tmin=300, Tmax=600, rank=0, comment='')
		
		for T in range(200, 1500, 10):
			self.assertEqual(kinetics.isTemperatureInRange(T), 300 <= T <= 600)
	
	def testArrheniusModel(self):
		"""
		
		"""
		
		Tlist = [T for T in range(300, 1500, 10)]
		
		A = 1.0e0; Ea = rmg.constants.R; n = 0.0
		kinetics = ArrheniusModel(A, Ea, n)
		for T in Tlist:
			self.assertAlmostEqual(kinetics.getRateConstant(T), A * math.exp(-Ea / rmg.constants.R / T), 4)
	
		A = 1.0e10; Ea = 50.0; n = 1.0
		kinetics = ArrheniusModel(A, Ea, n)
		for T in Tlist:
			self.assertAlmostEqual(kinetics.getRateConstant(T), A * T ** n * math.exp(-Ea / rmg.constants.R / T), 4)
	
	def testArrheniusEPModel(self):
		"""
		
		"""
		
		Tlist = [T for T in range(300, 1500, 10)]
		dHrxnList = [dHrxn for dHrxn in range(-500000, 500000, 10000)]
		
		A = 1.0e0; Ea = rmg.constants.R; n = 0.0; alpha = 0.0
		kinetics = ArrheniusEPModel(A, Ea, n, alpha)
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * math.exp(-Ea / rmg.constants.R / T), 4)
	
		A = 1.0e10; Ea = 50.0; n = 1.0; alpha = 0.0
		kinetics = ArrheniusEPModel(A, Ea, n, alpha)
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * T ** n * math.exp(-Ea / rmg.constants.R / T), 4)
		
		A = 3.0e5; Ea = 100.0; n = 0.0; alpha = -0.5
		kinetics = ArrheniusEPModel(A, Ea, n, alpha)
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * math.exp(-(Ea + alpha * dHrxn) / rmg.constants.R / T), 4)
		
		A = 6.6e10; Ea = 250.0; n = 1.0; alpha = 0.5
		kinetics = ArrheniusEPModel(A, Ea, n, alpha)
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * T ** n * math.exp(-(Ea + alpha * dHrxn) / rmg.constants.R / T), 4)
		
################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )