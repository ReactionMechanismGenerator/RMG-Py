#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('.')

import pylab
import numpy
import math
		
from rmg.species import *
from rmg.reaction import *
from rmg.model import *
from rmg.thermo.model import *
from rmg.system.batch import *

################################################################################

def initializeRMGSimulation(T, P, tf, model, system):
	"""
	Initialize an RMG simulation at temperature `T`, pressure `P`, and solution
	termination time `tf`. Returns the initialized reaction model `model` and
	reaction system `system`.
	"""

	model.termination.append(TerminationTime(tf))
	model.absoluteTolerance = 1e-24
	model.relativeTolerance = 1e-12

	system.equationOfState = IdealGas()
	system.temperatureModel = TemperatureModel()
	system.temperatureModel.setIsothermal(T)
	system.pressureModel = PressureModel()
	system.pressureModel.setIsobaric(P)

	return model, system

def runRMGSimulation(model, system):

	t, y, dydt, valid, species = system.simulate(model)
	
	# Reshape y into a matrix rather than a list of lists
	y0 = numpy.zeros((len(t), len(y[0])), float)
	for i, u in enumerate(y):
		for j, v in enumerate(u):
			y0[i,j] = v

	return t, y0

def postprocessRMGOutput(t, y, model=None):

	conc = numpy.zeros((len(t), len(y[0])-3), float)

	for i in range(y.shape[0]):
		conc[i,:] = y[i,3:] / y[i,1]

	# Make concentration plot and show
	pylab.figure()
	pylab.plot(t[1:], conc[1:,:])
	pylab.xlabel('Time (s)')
	pylab.ylabel('Concentration (mol/m^3)')
	if model:
		pylab.legend([s.label for s in model.core.species])
	pylab.show()
	

################################################################################

class SimulationCheck(unittest.TestCase):

	T = 1000.0 # [=] K

	P = 1.0e5 # [=] Pa

	tf = 10.0 # [=] s

	def testIrreversibleAtoB(self):
		"""
		A simple isomerization reaction A --> B, with the
		thermodynamics designed for an equilibrium of all B. This occurs in an
		isothermal, isobaric, homogeneous batch reactor.
		"""

		model = CoreEdgeReactionModel()
		system = BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, model=model, system=system)

		speciesA = Species(1, 'A')
		speciesA.thermoData = ThermoGAModel(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesA)

		speciesB = Species(2, 'B')
		speciesB.thermoData = ThermoGAModel(-500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesB)

		reactionAB = Reaction([speciesA], [speciesB])
		reactionAB.kinetics = ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)
		model.addReactionToCore(reactionAB)

		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[speciesA] = 1.0 / molarVolume

		t, y = runRMGSimulation(model, system)
		postprocessRMGOutput(t, y, model)
		pylab.title('A --> B, irreversible')

		# Check equilibrium
		self.assertTrue(y[-1,3] < 0.0001 * y[0,3])
		self.assertTrue(y[-1,4] > 0.9999 * y[0,3])
		self.assertAlmostEqual(y[-1,1], y[0,1], 3)

		# Check kinetics
		for i in range(len(t)):
			if abs(t[i] - 1.0) < 0.0001:
				self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0), 4)


	def testEquimolarAtoB(self):
		"""
#		A simple isomerization reaction A --> B, with the
#		thermodynamics designed for an equimolar equilibrium. This occurs in an
#		isothermal, isobaric, homogeneous batch reactor.
#		"""

		model = CoreEdgeReactionModel()
		system = BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, model=model, system=system)

		speciesA = Species(1, 'A')
		speciesA.thermoData = ThermoGAModel(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesA)

		speciesB = Species(2, 'B')
		speciesB.thermoData = ThermoGAModel(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesB)

		reactionAB = Reaction([speciesA], [speciesB])
		reactionAB.kinetics = ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)
		model.addReactionToCore(reactionAB)

		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[speciesA] = 1.0 / molarVolume

		t, y = runRMGSimulation(model, system)
		postprocessRMGOutput(t, y, model)
		pylab.title('A --> B, equimolar')

		# Check equilibrium
		self.assertAlmostEqual(y[-1,3], y[-1,4], 4)
		self.assertAlmostEqual(y[-1,1], y[0,1], 3)

	def testReversibleAtoB(self):
		"""
		A simple isomerization reaction A --> B, with the
		thermodynamics designed for an equimolar equilibrium. This occurs in an
		isothermal, isobaric, homogeneous batch reactor.
		"""

		model = CoreEdgeReactionModel()
		system = BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, model=model, system=system)

		speciesA = Species(1, 'A')
		speciesA.thermoData = ThermoGAModel(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesA)

		speciesB = Species(2, 'B')
		speciesB.thermoData = ThermoGAModel(494237.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesB)

		reactionAB = Reaction([speciesA], [speciesB])
		reactionAB.kinetics = ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)
		model.addReactionToCore(reactionAB)

		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[speciesA] = 1.0 / molarVolume

		t, y = runRMGSimulation(model, system)
		postprocessRMGOutput(t, y)
		pylab.title('A --> B, reversible')

		# Check equilibrium
		self.assertAlmostEqual(2.0 * y[-1,3], y[-1,4], 3)
		self.assertAlmostEqual(y[-1,1], y[0,1], 3)

	def testIrreversible2AtoB(self):
		"""
		A simple association reaction 2A --> B, with the
		thermodynamics designed for an equilibrium of all B. This occurs in an
		isothermal, isobaric, homogeneous batch reactor.
		"""

		model = CoreEdgeReactionModel()
		system = BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, model=model, system=system)

		speciesA = Species(1, 'A')
		speciesA.thermoData = ThermoGAModel(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesA)

		speciesB = Species(2, 'B')
		speciesB.thermoData = ThermoGAModel(-500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesB)

		reactionAB = Reaction([speciesA, speciesA], [speciesB])
		reactionAB.kinetics = ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)
		model.addReactionToCore(reactionAB)

		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[speciesA] = 1.0 / molarVolume

		t, y = runRMGSimulation(model, system)
		postprocessRMGOutput(t, y)
		pylab.title('2A --> B, irreversible')

		# Check equilibrium
		self.assertTrue(y[-1,3] < 0.01 * y[0,3])
		self.assertTrue(y[-1,4] > 0.49 * y[0,3] and y[-1,4] < 0.51 * y[0,3])
		self.assertTrue(y[-1,1] > 0.49 * y[0,1] and y[-1,1] < 0.51 * y[0,1])

	def testEquimolar2AtoB(self):
		"""
		A simple association reaction 2A --> B, with the thermodynamics  
		designed for an equimolar mixture of A and B. This occurs in an
		isothermal, isobaric, homogeneous batch reactor.
		"""

		model = CoreEdgeReactionModel()
		system = BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, model=model, system=system)

		speciesA = Species(1, 'A')
		speciesA.thermoData = ThermoGAModel(250000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesA)

		speciesB = Species(2, 'B')
		speciesB.thermoData = ThermoGAModel(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesB)

		reactionAB = Reaction([speciesA, speciesA], [speciesB])
		reactionAB.kinetics = ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)
		model.addReactionToCore(reactionAB)

		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[speciesA] = 1.0 / molarVolume

		t, y = runRMGSimulation(model, system)
		postprocessRMGOutput(t, y)
		pylab.title('2A --> B, equimolar')

	def testIrreversibleAto2B(self):
		"""
		A simple dissociation reaction A --> 2B, with the
		thermodynamics designed for an equilibrium of all B. This occurs in an
		isothermal, isobaric, homogeneous batch reactor.
		"""

		model = CoreEdgeReactionModel()
		system = BatchReactor()
		initializeRMGSimulation(T=self.T, P=self.P, tf=self.tf, model=model, system=system)

		speciesA = Species(1, 'A')
		speciesA.thermoData = ThermoGAModel(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesA)

		speciesB = Species(2, 'B')
		speciesB.thermoData = ThermoGAModel(-500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		model.addSpeciesToCore(speciesB)

		reactionAB = Reaction([speciesA], [speciesB, speciesB])
		reactionAB.kinetics = ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)
		model.addReactionToCore(reactionAB)

		molarVolume = system.equationOfState.getVolume(T=self.T, P=self.P, N=1.0)
		system.initialConcentration[speciesA] = 1.0 / molarVolume

		t, y = runRMGSimulation(model, system)
		postprocessRMGOutput(t, y)
		pylab.title('A --> 2B, irreversible')

		# Check equilibrium
		self.assertTrue(y[-1,3] < 0.01 * y[0,3])
		self.assertTrue(y[-1,4] > 1.99 * y[0,3] and y[-1,4] < 2.01 * y[0,3])
		self.assertTrue(y[-1,1] > 1.99 * y[0,1] and y[-1,1] < 2.01 * y[0,1])

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )