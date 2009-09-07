#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

import math
		
import rmg.main as main
import rmg.data as data
import rmg.species as species
import rmg.reaction as reaction
import rmg.thermo as thermo

################################################################################

class ThermoGACheck(unittest.TestCase):                          

	def testHeatCapacity(self):
		"""
		A test of the method of calculating heat capacity from heat capacity data.
		The Cp data is selected to follow the formula
		
		.. math:: C_p(T) = 0.1 T - 20 \hspace{20pt} 300 < T < 1500
		
		The implementation specified Cp values at 300, 400, 500, 600, 800, 1000,
		and 1500 K and interpolated when necessary. Above 1500 K the Cp value is
		assumed to be equal to Cp(1500 K).
		"""
		
		# Heat capacity: 
		#		Cp = 0.1 * T - 20.0		300.0 < T < 1500.0
		#		Cp = 130.0				T > 1500.0
		thermoData = species.ThermoGAData(0, 0, [10, 20, 30, 40, 60, 80, 130.0])
		
		Tlist = [T for T in range(300, 1500, 10)]
		for T in Tlist:
			Cp = 0.1 * T - 20.0
			self.assertAlmostEqual(thermoData.getHeatCapacity(T), Cp, 4)
		
		T = 1500; Cp = 0.1 * T - 20.0
		Tlist = [T for T in range(1600, 2000, 50)]
		for T in Tlist:
			self.assertAlmostEqual(thermoData.getHeatCapacity(T), Cp, 4)
		
	def testEnthalpy(self):
		"""
		A test of the method of calculating enthalpy from heat capacity data.
		The Cp data is selected to follow the formula
		
		.. math:: C_p(T) = 0.1 T - 20 \hspace{20pt} 300 < T < 1500
		
		The corresponding enthalpy formula is
		
		.. math:: S(T) = S_0 + 0.1 (T^2 - T_0^2) - 20 (T - T_0) \hspace{20pt} 300 < T < 1500
		
		where :math:`T_0` is taken to be 300 K.
		"""
		
		H0 = 800000.0
		
		thermoData = species.ThermoGAData(H0, 0, [10, 20, 30, 40, 60, 80, 130.0])
		
		Tlist = [T for T in range(300, 1500, 10)]
		for T in Tlist:
			H = H0 + 0.05 * (T**2 - 300.0**2) - 20 * (T - 300.0)
			self.assertAlmostEqual(thermoData.getEnthalpy(T), H, 4)
		
		T = 1500; H0 += 0.05 * (T**2 - 300.0**2) - 20 * (T - 300.0)
		Tlist = [T for T in range(1600, 2000, 50)]
		for T in Tlist:
			H = H0 + 130 * (T - 1500.0)
			self.assertAlmostEqual(thermoData.getEnthalpy(T), H, 4)
	
	def testEntropy(self):
		"""
		A test of the method of calculating entropy from heat capacity data.
		The Cp data is selected to follow the formula
		
		.. math:: C_p(T) = 0.1 T - 20 \hspace{20pt} 300 < T < 1500
		
		The corresponding entropy formula is
		
		.. math:: S(T) = S_0 + 0.1 (T - T_0) - 20 \ln \left( \frac{T}{T_0} \right) \hspace{20pt} 300 < T < 1500
		
		where :math:`T_0` is taken to be 300 K.
		"""
		
		S0 = 500.0
		thermoData = species.ThermoGAData(0, S0, [10, 20, 30, 40, 60, 80, 130.0])
		
		Tlist = [T for T in range(300, 1500, 10)]
		for T in Tlist:
			S = S0 + 0.1 * (T - 300.0) - 20 * math.log(T/300.0)
			self.assertAlmostEqual(thermoData.getEntropy(T), S, 4)
		
		T = 1500; S0 += 0.1 * (T - 300.0) - 20 * math.log(T/300.0)
		Tlist = [T for T in range(1600, 2000, 50)]
		for T in Tlist:
			S = S0 + 130 * math.log(T/1500.0)
			self.assertAlmostEqual(thermoData.getEntropy(T), S, 4)
		
################################################################################

class ThermoEstimationCheck(unittest.TestCase):                          
	
	def test1C2H6(self):
		
		C2H6 = species.Species()
		C2H6.fromSMILES('CC')
		C2H6.getThermoData()
#		print 'ethane'
#		print 'H(298 K) = %s' % (C2H6.getEnthalpy(298) / 4184)
#		print 'S(298 K) = %s' % (C2H6.getEntropy(298) / 4.184)
#		print 'G(298 K) = %s' % (C2H6.getFreeEnergy(298) / 4184)
#		print 'Cp(298 K) = %s' % (C2H6.getHeatCapacity(300) / 4.184)
		self.assertAlmostEqual(C2H6.getEnthalpy(298) / 4184, -20.4, 1)
		self.assertAlmostEqual(C2H6.getEntropy(298) / 4.184, 55.1, 1)
		self.assertAlmostEqual(C2H6.getFreeEnergy(298) / 4184, -36.8, 1)
		self.assertAlmostEqual(C2H6.getHeatCapacity(298) / 4.184, 12.5, 1)

	def test2CH3(self):
		
		CH3 = species.Species()
		CH3.fromSMILES('[CH3]')
		CH3.getThermoData()
#		print 'methyl'
#		print 'H(298 K) = %s' % (CH3.getEnthalpy(298) / 4184)
#		print 'S(298 K) = %s' % (CH3.getEntropy(298) / 4.184)
#		print 'G(298 K) = %s' % (CH3.getFreeEnergy(298) / 4184)
#		print 'Cp(298 K) = %s' % (CH3.getHeatCapacity(300) / 4.184)
		self.assertAlmostEqual(CH3.getEnthalpy(298) / 4184, 34.8, 1)
		self.assertAlmostEqual(CH3.getEntropy(298) / 4.184, 46.4, 1)
		self.assertAlmostEqual(CH3.getFreeEnergy(298) / 4184, 21.0, 1)
		self.assertAlmostEqual(CH3.getHeatCapacity(298) / 4.184, 9.54, 1)

################################################################################

if __name__ == '__main__':
	
	# Show debug messages (as databases are loading)
	main.initializeLog(10)

	# Load databases
	databasePath = '../data/RMG_database'
	
	# Create and load thermo databases
	thermo.thermoDatabase = thermo.ThermoDatabaseSet()
	thermo.thermoDatabase.load(databasePath + '/')

	# Create and load forbidden structures
	thermo.forbiddenStructures = data.Dictionary()
	thermo.forbiddenStructures.load(databasePath + '/forbiddenStructure.txt')
	thermo.forbiddenStructures.toStructure()

	unittest.main()