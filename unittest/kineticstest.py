#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

import math
		
from rmg.reaction import *
import rmg.constants

################################################################################

class KineticsCheck(unittest.TestCase):                          

	def testKinetics(self):
		"""
		
		"""
		
		kinetics = Kinetics([300, 600], 0, '')
		
		for T in range(200, 1500, 10):
			self.assertEqual(kinetics.isTemperatureInRange(T), 300 <= T <= 600)
	
	def testArrheniusKinetics(self):
		"""
		
		"""
		
		A = 1.0e0; Ea = rmg.constants.R; n = 0.0
		kinetics = ArrheniusKinetics(A, Ea, n)
		
		Tlist = [T for T in range(300, 1500, 10)]
		for T in Tlist:
			self.assertAlmostEqual(kinetics.getRateConstant(T), A * math.exp(-Ea / rmg.constants.R / T), 4)
	
		A = 1.0e10; Ea = 50.0; n = 1.0
		kinetics = ArrheniusKinetics(A, Ea, n)
		for T in Tlist:
			self.assertAlmostEqual(kinetics.getRateConstant(T), A * T ** n * math.exp(-Ea / rmg.constants.R / T), 4)
	
	def testArrheniusEPKinetics(self):
		"""
		
		"""
		
		A = 1.0e0; Ea = rmg.constants.R; n = 0.0; alpha = 0.0
		kinetics = ArrheniusEPKinetics(A, Ea, n, alpha)
		
		Tlist = [T for T in range(300, 1500, 10)]
		dHrxnList = [dHrxn for dHrxn in range(-500000, 500000, 10000)]
		
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * math.exp(-Ea / rmg.constants.R / T), 4)
	
		A = 1.0e10; Ea = 50.0; n = 1.0; alpha = 0.0
		kinetics = ArrheniusEPKinetics(A, Ea, n, alpha)
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * T ** n * math.exp(-Ea / rmg.constants.R / T), 4)
		
		A = 3.0e5; Ea = 100.0; n = 0.0; alpha = -0.5
		kinetics = ArrheniusEPKinetics(A, Ea, n, alpha)
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * math.exp(-(Ea + alpha * dHrxn) / rmg.constants.R / T), 4)
		
		A = 6.6e10; Ea = 250.0; n = 1.0; alpha = 0.5
		kinetics = ArrheniusEPKinetics(A, Ea, n, alpha)
		for dHrxn in dHrxnList:
			for T in Tlist:
				self.assertAlmostEqual(kinetics.getRateConstant(T, dHrxn), A * T ** n * math.exp(-(Ea + alpha * dHrxn) / rmg.constants.R / T), 4)
		
################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )