#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.append('.')

import os.path
import shutil
import unittest
import pylab
import numpy
import math
import quantities as pq
import logging

import Cantera
import Cantera.Reactor

import rmg.cantera_loader
import rmg.model

from simulationtest import initializeRMGSimulation, runRMGSimulation, postprocessRMGOutput

################################################################################

def initializeCanteraSimulation(filepath, T, P, tf):
	"""
	Initialize a Cantera simulation at temperature `T`, pressure `P`, and solution
	termination time `tf`.
	"""
	
	# we change into the scratch folder for this because it creates xml files
	absfilepath = os.path.abspath(filepath)
	oldwd	= os.getcwd()
	try:
		os.chdir(rmg.constants.scratchDirectory)
		# load the mechanism into gas
		gas = Cantera.importPhase(absfilepath,'chem')
		# create the environment
		gasAir = Cantera.Air()
	finally:
		os.chdir(oldwd)
	
	# set the inital gas and environment conditions
	gas.set(T=T, P=P)
	gasAir.set(T=T, P=P)
	
	# create a reactor for the batch reactor
	# and a reservoir for the environment
	reactor = Cantera.Reactor.Reactor(gas, volume = 1.0)
	environment = Cantera.Reactor.Reservoir(gasAir)
	
	# Define a wall between the reactor and the environment, and
	# make it flexible, so that the pressure in the reactor is held
	# at the environment pressure, and conductive so the temperature likewise
	wall = Cantera.Reactor.Wall(reactor,environment)
	wall.set(K = 1.0e12)   # set expansion parameter. dV/dt = KA(P_1 - P_2)
	wall.set(A = 1.0)
	wall.setHeatTransferCoeff(1.0e15) # W/m2/K
	
	# put the reactor into a reactor network
	sim = Cantera.Reactor.ReactorNet([reactor]) 
	sim.setInitialTime(0.0)
	maxtime = tf
	
	return gas, gasAir, reactor, environment, wall, sim, maxtime

def runCanteraSimulation(times, gas, gasAir, reactor, environment, wall, sim, maxtime):
	
	t = []
	y = []
	
	sim.step(maxtime)
	for time in times:
		if time: # don't integrate if time==0
			sim.advance(time) # now try to hit the rmg_t spot on!
		t.append(sim.time())
		output = [gas.pressure(), reactor.volume(), gas.temperature()]
		molarDensity = float( pq.quantity.Quantity(gas.molarDensity(),'kmol/m^3').simplified )
		# total_molar_density now in RMG units of mol/m^3
		output.extend([gas.moleFraction(i)*molarDensity for i in range(gas.nSpecies())])
		y.append(output)
	
	# Reshape y into a matrix rather than a list of lists
	y0 = numpy.zeros((len(t), len(y[0])), float)
	for i, u in enumerate(y):
		for j, v in enumerate(u):
			y0[i,j] = v
	
	return t, y0

def postprocessCanteraOutput(t, y):
	"""Make concentration plot and show"""
	#pylab.figure()
	pylab.plot(t[1:], y[1:,3:])
	pylab.xlabel('Time (s)')
	pylab.ylabel('Concentration (mol/m^3)')
	pylab.show()

################################################################################

class CanteraLoaderCheck(unittest.TestCase):
	
	testfolder = 'canteraloadertest'
	
	T = 1000.0 # [=] K
	P = 1.0e5 # [=] Pa
	tf = 10.0 # [=] s
	
	
	def failUnlessAlmostEqual(self, first, second, accuracy=5, msg=None):
		"""Fail if the two objects are unequal as determined by the log10 of
		   their difference divided by their average and comparing to the 
		   required accuracy (default 5).  
		   Accuracy is like significant figures, but needn't be an integer.
		"""
		import math
		a=first
		b=second
		if a==b: return
		if a*b<0 or -math.log10(abs((a-b)*2.0/(a+b))) < accuracy:
			raise self.failureException, \
			   (msg or '%r != %r with accuracy %g' % (first, second, accuracy))
	assertAlmostEqual = assertAlmostEquals = failUnlessAlmostEqual

	def setUp(self):
		"""setUp gets called before each test"""
		from rmg import constants, log
		import rmg.cantera_loader
		import ctml_writer as cti
		
		constants.scratchDirectory = os.path.join(self.testfolder,'temp')
		if os.path.isdir(constants.scratchDirectory): shutil.rmtree(constants.scratchDirectory)
		os.makedirs(constants.scratchDirectory)
		
		if cti._species:
			reload(cti)
		if rmg.cantera_loader._species:
			reload(rmg.cantera_loader)
		log.initialize(verbose=20, fstr='RMG.log')
		
		#pylab.figure(1)
	
	def tearDown(self):
		"""tearDown gets called after each test"""
		pass

	def test1CanteraLoad(self):
		filename = 'canteraHXD13.cti'
		filepath = os.path.join(self.testfolder,filename)
		model = rmg.cantera_loader.loadCanteraFile(filepath)
		self.assertTrue(len(model.core.species) == 21)
		self.assertTrue(len(model.core.reactions) == 33)

	def test2ChemkinLoad(self):
		filename = 'chemkinHXD13.inp'
		filepath = os.path.join(self.testfolder,filename)
		model = rmg.cantera_loader.loadChemkinFile(filepath)
		self.assertTrue(len(model.core.species) == 21)
		self.assertTrue(len(model.core.reactions) == 33)

	def test3CanteraVsRMGSimple(self):

		# run it in RMG
		filename = 'canteraA=B_2B.cti'
		filepath = os.path.join(self.testfolder,filename)

		model = rmg.cantera_loader.loadCanteraFile(filepath)
		system = rmg.model.BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, 
									model=model, system=system)
		
		speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[speciesA] = 1.0/molarVolume
		
		rmg_t, rmg_y = runRMGSimulation(model, system)
		postprocessRMGOutput(rmg_t, rmg_y)
		
		#run it in Cantera
		gas, gasAir, reactor, environment, wall, sim, maxtime = \
			initializeCanteraSimulation(filepath=filepath, T=self.T, P=self.P, tf=self.tf)
		
		speciesA_index = gas.speciesIndex('A')
		gas.setMoleFractions("A:1")

		cantera_t, cantera_y = runCanteraSimulation(rmg_t, gas, gasAir, 
								reactor, environment, wall, sim, maxtime)
		postprocessCanteraOutput(cantera_t, cantera_y)

		#check the results are the same shape
		self.assert_( cantera_y.shape == rmg_y.shape )

		#compare pressures
		for i in rmg_y[:,0]/cantera_y[:,0]:
			self.assertAlmostEqual(i,1.0,3)
		#compare volumes
		for i in rmg_y[:,1]/cantera_y[:,1]:
			self.assertAlmostEqual(i,1.0,3)
		#compare temperatures
		for i in rmg_y[:,2]/cantera_y[:,2]:
			self.assertAlmostEqual(i,1.0,3)

		#compare concentration profiles
		ratio = rmg_y[1:,3:] / cantera_y[1:,3:]
		for i in ratio.reshape(ratio.size,1):
			self.assertAlmostEqual(i,1.0,3)

	def test4CanteraVsRMGComplex(self):

		# run it in RMG
		filename = 'canteraHXD13.cti'
		filepath = os.path.join(self.testfolder,filename)

		model = rmg.cantera_loader.loadCanteraFile(filepath)
		system = rmg.model.BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, model=model, system=system)

		HXD13 = rmg.cantera_loader._speciesByName['HXD13(1)'].getRmgSpecies()
		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[HXD13] = 1.0/molarVolume

		logging.info('Running RMG simulation...')

		rmg_t, rmg_y = runRMGSimulation(model, system)
		postprocessRMGOutput(rmg_t, rmg_y, model)

		gas, gasAir, reactor, environment, wall, sim, maxtime = \
			initializeCanteraSimulation(filepath=filepath, T=self.T, P=self.P, tf=self.tf)

		speciesA = gas.speciesIndex('HXD13(1)')
		gas.setMoleFractions("HXD13(1):1")

		logging.info('Running cantera simulation...')

		cantera_t, cantera_y = runCanteraSimulation(rmg_t, gas, gasAir, 
									reactor, environment, wall, sim, maxtime)
		postprocessCanteraOutput(cantera_t, cantera_y)
		
		#check the results are the same shape
		self.assert_( cantera_y.shape == rmg_y.shape )
		
		#compare times
		for i in range(1,len(rmg_t)): #skip time[0]
			# NB. assertAlmostEqual compares *decimal places* not significant 
			# figures so instead of aAE(A,B,3) do aAE(A/B,1.0,3) 
			self.assertAlmostEqual(rmg_t[i] / cantera_t[i], 1.0, 3, 
				"Times diverged at %g sec"%rmg_t[i] )
		#compare pressures
		for i in rmg_y[:,0]/cantera_y[:,0]:
			self.assertAlmostEqual(i,1.0,2)
		#compare volumes
		for i in range(len(rmg_t)):
			self.assertAlmostEqual(rmg_y[i,1]/cantera_y[i,1], 1.0, 2, 
				"Volumes diverged at %g sec"%rmg_t[i] )
		#compare temperatures
		for i in rmg_y[:,2]/cantera_y[:,2]:
			self.assertAlmostEqual(i,1.0,2)
		
		#compare concentration profiles
		ratio = rmg_y[1:,3:] / cantera_y[1:,3:]
		for i in ratio.reshape(ratio.size,1):
			self.assertAlmostEqual(i,1.0,2)

################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
	##suite = unittest.TestLoader().loadTestsFromTestCase(CanteraLoaderCheck)
	##unittest.TextTestRunner(verbosity=2).run(suite)
	# this is handy to manually interrupt a test so you can 
	# play with the PDB debugger
	# CanteraLoaderCheck('test4CanteraVsRMGComplex').debug()