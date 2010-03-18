#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

import math
		
import rmg.data as data
import rmg.species as species
import rmg.structure as structure
import rmg.thermo as thermo
import rmg.log as logging

# Run this whether being run as __main__ or called by other unit test suite:
# Load databases
databasePath = '../data/RMG_database'

# Create and load thermo databases
species.thermoDatabase = species.ThermoDatabaseSet()
species.thermoDatabase.load(databasePath + '/')

# Create and load forbidden structures
thermo.forbiddenStructures = data.Dictionary()
thermo.forbiddenStructures.load(databasePath + '/forbiddenStructure.txt')
thermo.forbiddenStructures.toStructure()
	

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
		thermoData = thermo.ThermoGAData(0, 0, [10, 20, 30, 40, 60, 80, 130.0])
		
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
		
		thermoData = thermo.ThermoGAData(H0, 0, [10, 20, 30, 40, 60, 80, 130.0])
		
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
		thermoData = thermo.ThermoGAData(0, S0, [10, 20, 30, 40, 60, 80, 130.0])
		
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
		C2H6.generateThermoData(thermoClass=thermo.ThermoGAData)
#		print 'ethane'
#		print 'H(298 K) = %s' % (C2H6.getEnthalpy(298) / 4184)
#		print 'S(298 K) = %s' % (C2H6.getEntropy(298) / 4.184)
#		print 'G(298 K) = %s' % (C2H6.getFreeEnergy(298) / 4184)
#		print 'Cp(298 K) = %s' % (C2H6.getHeatCapacity(300) / 4.184)
		self.assertAlmostEqual(C2H6.getEnthalpy(298) / 4184, -20.4, 1)
		self.assertAlmostEqual(C2H6.getEntropy(298) / 4.184, 55.1, 1)
		self.assertAlmostEqual(C2H6.getFreeEnergy(298) / 4184, -36.8, 1)
		self.assertAlmostEqual(C2H6.getHeatCapacity(298) / 4.184, 12.38, 1)

	def test2CH3(self):
		
		CH3 = species.Species()
		CH3.fromSMILES('[CH3]')
		CH3.generateThermoData(thermoClass=thermo.ThermoGAData)
#		print 'methyl'
#		print 'H(298 K) = %s' % (CH3.getEnthalpy(298) / 4184)
#		print 'S(298 K) = %s' % (CH3.getEntropy(298) / 4.184)
#		print 'G(298 K) = %s' % (CH3.getFreeEnergy(298) / 4184)
#		print 'Cp(298 K) = %s' % (CH3.getHeatCapacity(300) / 4.184)
		self.assertAlmostEqual(CH3.getEnthalpy(298) / 4184, 34.8, 1)
		self.assertAlmostEqual(CH3.getEntropy(298) / 4.184, 46.4, 1)
		self.assertAlmostEqual(CH3.getFreeEnergy(298) / 4184, 21.0, 1)
		self.assertAlmostEqual(CH3.getHeatCapacity(300) / 4.184, 9.14, 1)

	def test3C6H9(self):
		"""
		This tests the effect of abstraction of a hydrogen from each of the
		six sites of a 1,3-hexadiene molecule to ensure that the appropriate
		thermo parameters are determined for each. It also checks the original
		1,3-hexadiene. The comparison is done in units of kcal/mol for enthalpy
		and free energy and cal/mol*K for entropy and heat capacity.
		"""

		C6H10 = species.Species()
		C6H10.fromSMILES('C=CC=CCC')
		C6H10.getResonanceIsomers()
		C6H10.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(C6H10.getEnthalpy(298) / 4184, 13.45, 1)
		self.assertAlmostEqual(C6H10.getEntropy(298) / 4.184, 86.37, 1)
		self.assertAlmostEqual(C6H10.getFreeEnergy(298) / 4184, -12.3, 1)
		self.assertAlmostEqual(C6H10.getHeatCapacity(298) / 4.184, 29.49, 1)

		C6H9 = species.Species()
		C6H9.fromSMILES('[CH]=CC=CCC')
		C6H9.getResonanceIsomers()
		C6H9.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(C6H9.getEnthalpy(298) / 4184, 72.55, 1)
		self.assertAlmostEqual(C6H9.getEntropy(298) / 4.184, 87.76, 1)
		self.assertAlmostEqual(C6H9.getHeatCapacity(298) / 4.184, 29.30, 1)

		C6H9 = species.Species()
		C6H9.fromSMILES('C=[C]C=CCC')
		C6H9.getResonanceIsomers()
		C6H9.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(C6H9.getEnthalpy(298) / 4184, 61.15, 1)
		self.assertAlmostEqual(C6H9.getEntropy(298) / 4.184, 87.07, 1)
		self.assertAlmostEqual(C6H9.getHeatCapacity(298) / 4.184, 29.68, 1)

		C6H9 = species.Species()
		C6H9.fromSMILES('C=C[C]=CCC')
		C6H9.getResonanceIsomers()
		C6H9.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(C6H9.getEnthalpy(298) / 4184, 61.15, 1)
		self.assertAlmostEqual(C6H9.getEntropy(298) / 4.184, 87.07, 1)
		self.assertAlmostEqual(C6H9.getHeatCapacity(298) / 4.184, 29.68, 1)

		C6H9 = species.Species()
		C6H9.fromSMILES('C=CC=[C]CC')
		C6H9.getResonanceIsomers()
		C6H9.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(C6H9.getEnthalpy(298) / 4184, 70.35, 1)
		self.assertAlmostEqual(C6H9.getEntropy(298) / 4.184, 88.18, 1)
		self.assertAlmostEqual(C6H9.getHeatCapacity(298) / 4.184, 29.15, 1)

		C6H9 = species.Species()
		C6H9.fromSMILES('C=CC=C[CH]C')
		C6H9.getResonanceIsomers()
		C6H9.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(C6H9.getEnthalpy(298) / 4184, 38.24, 1)
		self.assertAlmostEqual(C6H9.getEntropy(298) / 4.184, 84.41, 1)
		self.assertAlmostEqual(C6H9.getHeatCapacity(298) / 4.184, 27.79, 1)

		C6H9 = species.Species()
		C6H9.fromSMILES('C=CC=CC[CH2]')
		C6H9.getResonanceIsomers()
		C6H9.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(C6H9.getEnthalpy(298) / 4184, 62.44, 1)
		self.assertAlmostEqual(C6H9.getEntropy(298) / 4.184, 89.78, 1)
		self.assertAlmostEqual(C6H9.getHeatCapacity(298) / 4.184, 28.72, 1)

	def test4H(self):
		"""
		This tests the lookup of thermo data from the primary thermo library
		using H as a test species. The comparison is done in units of kcal/mol
		for enthalpy and free energy and cal/mol*K for entropy and heat
		capacity.
		"""

		H = species.Species()
		H.fromSMILES('[H]')
		H.getResonanceIsomers()
		H.generateThermoData(thermoClass=thermo.ThermoGAData)
		self.assertAlmostEqual(H.getEnthalpy(298) / 4184, 52.103, 1)
		self.assertAlmostEqual(H.getEntropy(298) / 4.184, 27.419, 1)
		self.assertAlmostEqual(H.getHeatCapacity(298) / 4.184, 4.968, 1)

################################################################################


class ThermoGAtoWilhoitCheck(unittest.TestCase):                          
	"""Test conversion of Group Additivity thermo to Wilhoit thermo"""
	def testWilhoitCreated(self):
		"""Can we make Wilhoit polynomial data"""
		GAthermoData = thermo.ThermoGAData(0, 0, [10, 20, 30, 40, 60, 80, 130.0])
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=2, rotors=0, linear=True )
		self.assertTrue(isinstance(WilhoitData,thermo.ThermoWilhoitData),"Didn't make ThermoWilhoitData instance")
		self.assertTrue(isinstance(WilhoitData,thermo.ThermoData),"Didn't make any kind of ThermoData instance")
		# well, if we didn't cause an error, I guess that's good enough for now.

	def testHeatCapacity(self):
		"""Check the Wilhoit Cp matches the GA Cp for propane.
		
		Uses Propane as a test-case. atoms=11, rotors=2, linear=False
		"""
		
		propane = structure.Structure(SMILES='CCC')
		propane.updateAtomTypes()
		GAthermoData = species.getThermoData(propane,required_class=thermo.ThermoGAData)
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=11, rotors=2, linear=False)
		
		Tlist = thermo.ThermoGAData.CpTlist # just check at defined data points
		for T in Tlist:
			ga = GAthermoData.getHeatCapacity(T)
			wil = WilhoitData.getHeatCapacity(T)
			err = abs(ga-wil)
			limit = 4.0 # J/mol/K
			self.assertTrue(err<limit,"GA (%.1f) and Wilhoit (%.1f) differ by more than %s J/mol/K at %dK"%(ga,wil,limit,T))
	
	def testEnthalpy(self):
		"""Check the Wilhoit H matches the GA H for propane.
		
		Uses Propane as a test-case. atoms=11, rotors=2, linear=False
		"""
		
		propane = structure.Structure(SMILES='CCC')
		propane.updateAtomTypes()
		GAthermoData = species.getThermoData(propane,required_class=thermo.ThermoGAData)
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=11, rotors=2, linear=False)
		
		Tlist = thermo.ThermoGAData.CpTlist # just check at defined data points
		for T in Tlist:
			ga = GAthermoData.getEnthalpy(T)
			wil = WilhoitData.getEnthalpy(T)
			err = abs(ga-wil)
			limit = 2000.0 # J/mol  # the wilhoit should be more accurate then trapezoid integration of GA, so wouldn't want them to be exactly the same
			self.assertTrue(err<limit,"GA (%.1f) and Wilhoit (%.1f) differ by more than %s J/mol at %dK"%(ga,wil,limit,T))
			
	def testEntropy(self):
		"""Check the Wilhoit S matches the GA S for propane.
		
		Uses Propane as a test-case. atoms=11, rotors=2, linear=False
		"""
		
		propane = structure.Structure(SMILES='CCC')
		propane.updateAtomTypes()
		GAthermoData = species.getThermoData(propane,required_class=thermo.ThermoGAData)
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=11, rotors=2, linear=False)
		
		Tlist = thermo.ThermoGAData.CpTlist # just check at defined data points
		for T in Tlist:
			ga = GAthermoData.getEntropy(T)
			wil = WilhoitData.getEntropy(T)
			err = abs(ga-wil)
			limit = 4.0 # J/mol/K
			self.assertTrue(err<limit,"GA (%.1f) and Wilhoit (%.1f) differ by more than %s J/mol/K at %dK"%(ga,wil,limit,T))
			
			
			
class ThermoWilhoitToNASACheck(unittest.TestCase):                          
	"""Test conversion from Wilhoit to NASA polynomials"""
	def testNASAcreated(self):
		"""Can we make NASA polynomial data"""
		cp0, cpInf, a0, a1, a2, a3, I, J = (1.0,1.0,1.0,1.0,1.0,1.0, 1.0, 1.0)
		comment = "Stupid thermo."
		WilhoitData = thermo.ThermoWilhoitData( cp0, cpInf, a0, a1, a2, a3, I, J, comment=comment)
		NASAthermoData = thermo.convertWilhoitToNASA(WilhoitData)
		# well, if we didn't cause an error, I guess that's good enough for now.

	def testHeatCapacity(self):
		"""Check the NASA Cp matches the GA Cp for propane.
		
		Uses Propane as a test-case. atoms=11, rotors=2, linear=False
		"""
		
		propane = structure.Structure(SMILES='CCC')
		propane.updateAtomTypes()
		GAthermoData = species.getThermoData(propane,required_class=thermo.ThermoGAData)
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=11, rotors=2, linear=False)
		NASAthermoData = thermo.convertWilhoitToNASA(WilhoitData)
		
		Tlist = thermo.ThermoGAData.CpTlist # just check at defined data points
		for T in Tlist:
			ga = GAthermoData.getHeatCapacity(T)
			nasa = NASAthermoData.getHeatCapacity(T)
			err = abs(ga-nasa)
			limit = 10.0 # J/mol/K
			self.assertTrue(err<limit,"GA (%.1f) and NASA (%.1f) differ by more than %s J/mol/K at %dK"%(ga,nasa,limit,T))
			
	def testEntropy(self):
		"""Check the NASA S matches the GA S for propane.
		
		Uses Propane as a test-case. atoms=11, rotors=2, linear=False
		"""
		
		propane = structure.Structure(SMILES='CCC')
		propane.updateAtomTypes()
		GAthermoData = species.getThermoData(propane,required_class=thermo.ThermoGAData)
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=11, rotors=2, linear=False)
		NASAthermoData = thermo.convertWilhoitToNASA(WilhoitData)
		
		Tlist = thermo.ThermoGAData.CpTlist # just check at defined data points
		for T in Tlist:
			ga = GAthermoData.getEntropy(T)
			nasa = NASAthermoData.getEntropy(T)
			err = abs(ga-nasa)
			limit = 4.0 # J/mol/K
			self.assertTrue(err<limit,"GA (%.1f) and NASA (%.1f) differ by more than %s J/mol/K at %dK"%(ga,nasa,limit,T))

	def testEnthalpy(self):
		"""Check the NASA H matches the GA H for propane.
		
		Uses Propane as a test-case. atoms=11, rotors=2, linear=False
		"""
		
		propane = structure.Structure(SMILES='CCC')
		propane.updateAtomTypes()
		GAthermoData = species.getThermoData(propane,required_class=thermo.ThermoGAData)
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=11, rotors=2, linear=False)
		NASAthermoData = thermo.convertWilhoitToNASA(WilhoitData)
		
		Tlist = thermo.ThermoGAData.CpTlist # just check at defined data points
		for T in Tlist:
			ga = GAthermoData.getEnthalpy(T)
			nasa = NASAthermoData.getEnthalpy(T)
			err = abs(ga-nasa)
			limit = 2000.0 # J/mol  # the wilhoit should be more accurate then trapezoid integration of GA, so wouldn't want them to be exactly the same
			self.assertTrue(err<limit,"GA (%.1f) and NASA (%.1f) differ by more than %s J/mol at %dK"%(ga,nasa,limit,T))

	def testObjectiveFunctionForOxygen(self):
		"""Check weighted objective function encountered during Wilhoit-to-NASA conversion for molecular oxygen)
		"""

		wilhoit = thermo.ThermoWilhoitData(3.5,4.5,-2.343,32.54,-79.26,47.75,8951,-18.19, B=0.5) #this is the scaled version
		q=thermo.TintOpt_objFun(1.0, wilhoit, .298, 6.0, 1, 3)#these are also scaled values
		expectedVal = 0.00018295170781357228 #taken from running in pure-python mode
		relErr = abs(q-expectedVal)/expectedVal
		limit = 0.01 #relative error limit (0.01=1%)
		self.assertTrue(relErr<limit,"Actual (%.8f) and expected (%.8f) differ by more than %s"%(q,expectedVal,limit*expectedVal))

	def testOxygenFromGA(self):
		"""Check conversion of GA values for molecular oxygen to NASA form
		"""

		oxygen = structure.Structure(SMILES='O=O')
		oxygen.updateAtomTypes()
		GAthermoData = species.getThermoData(oxygen,required_class=thermo.ThermoGAData)
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=2, rotors=0, linear=True)
		print WilhoitData
		NASAthermoData = thermo.convertWilhoitToNASA(WilhoitData)
		#if not working properly, this may produce an error in logging

	def testOxygenFromGA2(self):
		"""Check conversion of current GA values for molecular oxygen to Wilhoit form

		This is an attempt to test a discrepancy between results on Greg's Windows computer (same in both C- and pure-python modes) and Josh's Linux computer (C-mode)
		Discrepancy could be during the numerically intensive step (which this tests) or the computation of the group values (which a unit test below tests)
		"""
		#the following values come from Greg's Windows computer; Josh's scaled result for thermoWilhoitData is thermo.ThermoWilhoitData(3.5,4.5,-2.343,32.54,-79.26,47.75,8951,-18.19, B=0.5)
		GAthermoData = thermo.ThermoGAData(H298=0.0, S298=205.026767175, Cp=[29.288, 30.208480000000002, 31.128960000000003, 32.007600000000004, 33.764880000000005, 34.936399999999999, 36.484480000000005], index="+1096+++0+1096+++0")
		WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=2, rotors=0, linear=True)
		limit = 0.01 #relative error limit (0.01=1%)
		a0e = -0.9324
		a0re = abs(WilhoitData.a0-a0e)/abs(a0e)
		self.assertTrue(a0re<limit,"Actual (%.8f) and expected (%.8f) differ by more than %s"%(WilhoitData.a0,a0e,limit*a0e))
		a1e = 26.18
		a1re = abs(WilhoitData.a1-a1e)/abs(a1e)
		self.assertTrue(a1re<limit,"Actual (%.8f) and expected (%.8f) differ by more than %s"%(WilhoitData.a1,a1e,limit*a1e))
		a2e = -70.47
		a2re = abs(WilhoitData.a2-a2e)/abs(a2e)
		self.assertTrue(a2re<limit,"Actual (%.8f) and expected (%.8f) differ by more than %s"%(WilhoitData.a2,a2e,limit*a2e))
		a3e = 44.12
		a3re = abs(WilhoitData.a3-a3e)/abs(a3e)
		self.assertTrue(a3re<limit,"Actual (%.8f) and expected (%.8f) differ by more than %s"%(WilhoitData.a3,a3e,limit*a3e))

	def testOxygenGA(self):
		"""Check molecular oxygen GA heat capacity values

		"""

		oxygen = structure.Structure(SMILES='O=O')
		oxygen.updateAtomTypes()
		GAthermoData = species.getThermoData(oxygen,required_class=thermo.ThermoGAData)
		#the following values come from Greg's Windows computer
		Cpe = [29.288, 30.208480000000002, 31.128960000000003, 32.007600000000004, 33.764880000000005, 34.936399999999999, 36.484480000000005]
		self.assertTrue(GAthermoData.Cp[0]==Cpe[0],"Actual (%.8f) and expected (%.8f) are different"%(GAthermoData.Cp[0],Cpe[0]))
		self.assertTrue(GAthermoData.Cp[1]==Cpe[1],"Actual (%.8f) and expected (%.8f) are different"%(GAthermoData.Cp[1],Cpe[1]))
		self.assertTrue(GAthermoData.Cp[2]==Cpe[2],"Actual (%.8f) and expected (%.8f) are different"%(GAthermoData.Cp[2],Cpe[2]))
		self.assertTrue(GAthermoData.Cp[3]==Cpe[3],"Actual (%.8f) and expected (%.8f) are different"%(GAthermoData.Cp[3],Cpe[3]))
		self.assertTrue(GAthermoData.Cp[4]==Cpe[4],"Actual (%.8f) and expected (%.8f) are different"%(GAthermoData.Cp[4],Cpe[4]))
		self.assertTrue(GAthermoData.Cp[5]==Cpe[5],"Actual (%.8f) and expected (%.8f) are different"%(GAthermoData.Cp[5],Cpe[5]))
		self.assertTrue(GAthermoData.Cp[6]==Cpe[6],"Actual (%.8f) and expected (%.8f) are different"%(GAthermoData.Cp[6],Cpe[6]))

	#in the below integral tests, the expected values come from running in pure-python mode
	def testWilhoitIntegralTM1(self):
		"""Check Wilhoit.IntegralTM1

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral_TM1(6.0) - w.integral_TM1(.298)
		self.assertAlmostEqual(ans, 14.05184933947314, 14)

	def testWilhoitIntegralT0(self):
		"""Check Wilhoit.IntegralT0

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral_T0(6.0) - w.integral_T0(.298)
		self.assertAlmostEqual(ans, 30.156466842356952, 14)

	def testWilhoitIntegralT1(self):
		"""Check Wilhoit.IntegralT1

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral_T1(6.0) - w.integral_T1(.298)
		self.assertAlmostEqual(ans, 100.16557858116997, 16)

	def testWilhoitIntegralT2(self):
		"""Check Wilhoit.IntegralT2

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral_T2(6.0) - w.integral_T2(.298)
		self.assertAlmostEqual(ans, 409.65039004796273, 16)

	def testWilhoitIntegralT3(self):
		"""Check Wilhoit.IntegralT3

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral_T3(6.0) - w.integral_T3(.298)
		self.assertAlmostEqual(ans, 1859.5032978322745, 16)

	def testWilhoitIntegralT4(self):
		"""Check Wilhoit.IntegralT4

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral_T4(6.0) - w.integral_T4(.298)
		self.assertAlmostEqual(ans, 8965.5578894745959, 16)

	def testWilhoitIntegral2TM1(self):
		"""Check Wilhoit.Integral2_TM1

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral2_TM1(6.0) - w.integral2_TM1(.298)
		self.assertAlmostEqual(ans, 67.864846008944539, 13)

	def testWilhoitIntegral2T0(self):
		"""Check Wilhoit.Integral2_T0

		"""
		w = thermo.ThermoWilhoitData(3.0,6.0,1.0,-2.0,3.0,-4.0,1234.5,789.0,B=0.5)
		ans = w.integral2_T0(6.0) - w.integral2_T0(.298)
		self.assertAlmostEqual(ans, 161.77012800697503, 13)

	def testNASAIntegral2TM1(self):
		"""Check NASA.Integral2_TM1

		"""
		n = thermo.ThermoNASAPolynomial(T_range=[.298,1.000], coeffs = [1, -.2, .3, -.4, .5, -6, 7])
		ans = n.integral2_TM1(1.0) - n.integral2_TM1(.298)
		self.assertAlmostEqual(ans, 1.1958658319514555, 15)

	def testNASAIntegral2T0(self):
		"""Check NASA.Integral2_T0

		"""
		n = thermo.ThermoNASAPolynomial(T_range=[.298,1.000], coeffs = [1, -.2, .3, -.4, .5, -6, 7])
		ans = n.integral2_T0(1.0) - n.integral2_T0(.298)
		self.assertAlmostEqual(ans, 0.71887383097545454, 15)

class ThermoCpToNASACheck(unittest.TestCase):
	"""Test conversion from Cp/R to NASA polynomials (using numerical integrals)"""
	def testNASAfromCp(self):
		"""Can we make NASA polynomial data from an arbitrary Cp function"""

		CpObject = thermo.ThermoNASAPolynomial(T_range=[0,8000], coeffs = [2.5, 3.0/1000, 7.0/1000000, 0.0, 0.0, 0, 0])
		NASAthermoData = thermo.convertCpToNASA(CpObject, 1.0, 2.0)
		#print NASAthermoData
		self.assertAlmostEqual(NASAthermoData.getEnthalpy(298.15), 1.0, 4)
		self.assertAlmostEqual(NASAthermoData.getEntropy(298.15), 2.0, 4)
		self.assertAlmostEqual(NASAthermoData.polynomials[0].c0, 2.5, 4)
		self.assertAlmostEqual(NASAthermoData.polynomials[0].c1, 3.0/1000, 7)
		self.assertAlmostEqual(NASAthermoData.polynomials[0].c2, 7.0/1000000, 10)
		self.assertAlmostEqual(NASAthermoData.polynomials[0].c3, 0.0, 4)
		self.assertAlmostEqual(NASAthermoData.polynomials[0].c4, 0.0, 4)
		self.assertAlmostEqual(NASAthermoData.polynomials[1].c0, 2.5, 4)
		self.assertAlmostEqual(NASAthermoData.polynomials[1].c1, 3.0/1000, 7)
		self.assertAlmostEqual(NASAthermoData.polynomials[1].c2, 7.0/1000000, 10)
		self.assertAlmostEqual(NASAthermoData.polynomials[1].c3, 0.0, 4)
		self.assertAlmostEqual(NASAthermoData.polynomials[1].c4, 0.0, 4)



################################################################################

		
		
# Run this only if being run independently
if __name__ == '__main__':	
	from timeit import Timer
	startup = """gc.enable() # enable garbage collection in timeit
import sys
sys.path.append('../source')
import rmg.thermo as thermo
import rmg.species as species
import rmg.structure as structure
from rmg.structure import Structure
propane = structure.Structure(SMILES='CCC')
propane.updateAtomTypes()
GAthermoData = species.getThermoData(propane,required_class=thermo.ThermoGAData)
WilhoitData = thermo.convertGAtoWilhoit(GAthermoData, atoms=11, rotors=2, linear=False)
NASAthermoData = thermo.convertWilhoitToNASA(WilhoitData)
"""
	test1 = "for T in range(600,1100,100):\n\tGAthermoData.getEnthalpy(T)"
	test2 = "for T in range(600,1100,100):\n\tWilhoitData.getEnthalpy(T)"
	test3 = "for T in range(600,1100,100):\n\tNASAthermoData.getEnthalpy(T)"
	print "****"
	print "Timing getEnthalpy:"
	t = Timer(test1,startup)
	times = t.repeat(repeat=5,number=1000)
	print " ThermoGAData    took   %.3f milliseconds (%s)"%(min(times), times)
	t = Timer(test2,startup)
	times = t.repeat(repeat=3,number=1000)
	print " ThermoWilhoitData took %.3f milliseconds (%s)"%(min(times), times)
	t = Timer(test3,startup)
	times = t.repeat(repeat=3,number=1000)
	print " ThermoNASAData   took  %.3f milliseconds (%s)"%(min(times), times)
	print "****\n\nContinuing with tests..."
	
	# Show debug messages (as databases are loading)
	logging.initialize(10,'RMG.log')
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
