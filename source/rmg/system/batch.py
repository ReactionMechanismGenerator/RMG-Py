#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains classes that represent batch reaction systems.
"""

import os
import os.path
import math
import rmg.log as logging
import numpy

import rmg.settings as settings
import rmg.ctml_writer as ctml_writer
import rmg.model as modelmodule
import rmg.species as species
import rmg.reaction as reaction
from base import ReactionSystem

################################################################################

class BatchReactor(ReactionSystem):
	"""
	A reaction system consisting of a single well-mixed batch reactor. The
	reactor can be configured to be isothermal, adiabatic, or having a known
	overall heat transfer coefficient; isobaric, isochoric, or having a known
	coefficient of expansion.

	.. math:: \\rho \\hat{C}_\\mathrm{p} \\frac{dT}{dt} = U A \\left( T_1 - T_2 \\right)

	.. math:: \\frac{dV}{dt} = K A \\left( P_1 - P_2 \\right)

	The attributes are:

	======================  ====================================================
	Attribute               Description
	======================  ====================================================
	`volume`                The reactor volume in m^3
	`area`                  The total surface area of the reactor in m^2
	`heatTransferCoeff`     The overall heat transfer coefficient of the reactor
	                        in W/m^2*K
	`expansionCoeff`        The coefficient of expansion in m/s*Pa
	`initialTemperature`    The initial temperature in K
	`initialPressure`       The initial pressure in Pa
	`initialMoleFraction`   Dictionary of initial mole fraction for each
	                        species; those not in dictionary have initial
	                        mole fraction of zero
	`reservoirTemperature`  The reservoir temperature in K
	`reservoirPressure`     The reservoir pressure in Pa
	======================  ====================================================

	"""

	def __init__(self):
		ReactionSystem.__init__(self)

	def fromXML(self, document, rootElement, speciesDict):

		# Read volume
		self.volume = float(document.getChildQuantity(rootElement, 'volume', required=True).simplified)

		# Read area
		self.area = float(document.getChildQuantity(rootElement, 'area', required=True).simplified)

		# Read temperature model
		isothermal = document.getChildElement(rootElement, 'isothermal', required=False)
		adiabatic = document.getChildElement(rootElement, 'adiabatic', required=False)
		heatTransferCoefficient = document.getChildElement(rootElement, 'heatTransferCoefficient', required=False)
		if isothermal:
			self.setIsothermal()
		elif adiabatic:
			self.setAdiabatic()
		elif heatTransferCoefficient:
			self.heatTransferCoeff = float(document.getQuantity(heatTransferCoefficient).simplified)
		else:
			raise Exception('No temperature model specified; should use <isothermal>, <adiabatic>, or <heatTransferCoefficient>.')

		# Read pressure/volume model
		isobaric = document.getChildElement(rootElement, 'isobaric', required=False)
		isochoric = document.getChildElement(rootElement, 'isochoric', required=False)
		expansionCoefficient = document.getChildElement(rootElement, 'expansionCoefficient', required=False)
		if isobaric:
			self.setIsobaric()
		elif isochoric:
			self.setIsochoric()
		elif expansionCoefficient:
			self.expansionCoeff = float(document.getQuantity(expansionCoefficient).simplified)
		else:
			raise Exception('No pressure model specified; should use <isobaric>, <isochoric>, or <expansionCoefficient>.')

		# Read initial conditions
		self.initialMoleFraction = {}
		initialConditions = document.getChildElement(rootElement, 'initialConditions', required=True)
		self.initialTemperature = float(document.getChildQuantity(initialConditions, 'temperature', required=True).simplified)
		self.initialPressure = float(document.getChildQuantity(initialConditions, 'pressure', required=True).simplified)
		moleFractions = document.getChildElements(initialConditions, 'moleFraction', required=True)
		for moleFraction in moleFractions:
			specID = str(document.getAttribute(moleFraction, 'speciesID', required=True))
			value = float(document.getElementText(moleFraction))
			self.initialMoleFraction[speciesDict[specID]] = value

		# Read reservoir conditions
		reservoirConditions = document.getChildElement(rootElement, 'reservoirConditions', required=True)
		self.reservoirTemperature = float(document.getChildQuantity(reservoirConditions, 'temperature', required=True).simplified)
		self.reservoirPressure = float(document.getChildQuantity(reservoirConditions, 'pressure', required=True).simplified)

		# Read physical property model
		element = document.getChildElement(rootElement, 'physicalPropertyModel', required=True)
		propModelType = document.getAttribute(element, 'type', required=True)
		if propModelType.lower() == 'idealgas':
			# Set the reaction system's pressure model to isobaric
			self.equationOfState = modelmodule.IdealGas()
		elif propModelType.lower() == 'incompressibleliquid':
			molarVolume = float(document.getChildQuantity(element, 'molarVolume', required=True).simplified)
			self.equationOfState = modelmodule.IncompressibleLiquid(
				P = self.initialPressure,
				T = self.initialTemperature,
				Vmol = molarVolume)
		else:
			raise InvalidInputFileException('Invalid physical property model type "' + propModelType + '".')

	def initializeCantera(self):
		"""Creata a Cantera instance. Call this once"""
		ctml_writer.units(length = "m", time = "s", quantity = "mol", act_energy = "J/mol")
		phase = ctml_writer.ideal_gas(name = "chem",
		      elements = " C H O N Ar He Si",
		      species = "all",
		      reactions = "all",
		      initial_state = ctml_writer.state(temperature = self.initialTemperature ,
		                        pressure = self.initialPressure )    )
		def has_species(sp):
			"""Return 1 is a species with name 's' belongs to the phase,
			or 0 otherwise. Redefined because the ctml_writer one doesn't work for us"""
			if sp in ctml_writer._speciesnames: return 1
			return 0
		phase.has_species = has_species
		#self._cantera = phase  # if we save it, we have to pickle it, and we can't pickle the has_species function
		ctml_writer.validate() # turns on validation

	def runCantera(self, model):
		"""
		Execute a simulation of the reaction system in Cantera. The procedure:
		(1) write a CTML (Cantera) file, (2) read it into Cantera, (3) create
		the reactor in Cantera, and (4) return the simulation results.
		"""

		# Create a folder in the scratch directory for Cantera files if needed
		cantera_folder = os.path.join(settings.scratchDirectory,'cantera')
		os.path.exists(cantera_folder) or os.mkdir(cantera_folder)

		# Write the CTML file to scratch/cantera/ folder
		cti_file = os.path.join(cantera_folder, 'cantera_input_%03d' % len(model.core.species))
		logging.debug("Writing CTML file %s" % cti_file)
		ctml_writer.dataset(cti_file) # change name
		# update the T and P used to convert PdepRate coefficients into rate constants
		ctml_writer._temperature = self.initialTemperature
		ctml_writer._pressure = self.initialTemperature
		ctml_writer.write()

		import Cantera
		import Cantera.Reactor

		# Load the CTML file into Cantera
		logging.info("Preparing Cantera simulation %d" % len(model.core.species))
		Cantera.reset()
		gas = Cantera.importPhase('%s.xml' % cti_file, 'chem', loglevel=1)

		# Set initial concentrations
		moleFractions = numpy.zeros(len(model.core.species))
		for spec, conc in self.initialMoleFraction.iteritems():
			moleFractions[gas.speciesIndex(str(spec))] = conc
		gas.setMoleFractions(moleFractions) # it normalises it to 1

		# Set initial temperature and pressure
		gas.set(T=self.initialTemperature, P=self.initialPressure)

		# create a batch reactor
		if self.heatTransferCoeff == 1.0e100:
			reactor = Cantera.Reactor.Reactor(gas, volume=self.volume, energy='off')
		else:
			reactor = Cantera.Reactor.Reactor(gas, volume=self.volume)

		# set the inital environment conditions
		gasAir = Cantera.Air()
		gasAir.set(T=self.reservoirTemperature, P=self.reservoirPressure)

		# create a reservoir for the environment
		environment = Cantera.Reactor.Reservoir(gasAir)

		# Define a wall between the reactor and the environment, and
		# make it flexible, so that the pressure in the reactor is held
		# at the environment pressure, and conductive so the temperature likewise
		wall = Cantera.Reactor.Wall(reactor, environment)
		wall.set(K=self.expansionCoeff)
		wall.set(A=self.area)
		wall.setHeatTransferCoeff(self.heatTransferCoeff) # W/m2/K

		# put reactor in a reactor network so it can be integrated
		sim = Cantera.Reactor.ReactorNet([reactor])
		sim.setTolerances(atol=model.absoluteTolerance, rtol=model.relativeTolerance)

		#import pdb; pdb.set_trace()
		return sim, gas

	def setIsothermal(self):
		"""
		Sets the reactor to use an isothermal model for temperature. This is
		currently achieved by simply choosing a very large heat transfer
		coefficient.
		"""
		self.heatTransferCoeff = 1.0e100

	def setAdiabatic(self):
		"""
		Sets the reactor to use an adiabatic model for temperature. This is
		achieved by setting the overall heat transfer coefficient to zero.
		"""
		self.heatTransferCoeff = 0.0

	def setIsobaric(self):
		"""
		Sets the reactor to use an isobaric model for pressure. This is
		currently achieved by simply choosing a very large coefficient of
		expansion.
		"""
		self.expansionCoeff = 1.0e12

	def setIsochoric(self):
		"""
		Sets the reactor to use an isochoric model for volume. This is
		achieved by setting the  coefficient of expansion to zero.
		"""
		self.expansionCoeff = 0.0

	def simulate(self, model):
		"""
		Conduct a simulation of the current reaction system using the core-edge
		reaction model `model`.

		Edge species fluxes are tracked, relative to the characteristic core
		flux at that time, throughout the simulation.
		If one exceeds `model.fluxToleranceInterrupt` the simulation
		is interrupted, and that species is returned.
		The highest relative flux reached by each species during the simulation
		is stored for later analysis.
		If one or more of these exceed `model.fluxToleranceMoveToCore` then the
		species with the highest will be returned.

		If the simulation completes without interruption, then any that fall
		below `model.fluxToleranceKeepInEdge` will be removed from the
		edge, along with the reactions that involve them.

		Returns:
		(tlist, ylist, dydtlist, valid?, Edge_species_with_highest_flux)
		"""

		# try writing cantera file
		sim,gas = self.runCantera(model)

		# Assemble stoichiometry matrix for all core and edge species
		# Rows are species by id, columns are reactions by id
		stoichiometry = model.getStoichiometryMatrix()

		tlist = []; ylist = []; dydtlist = []
		maxRelativeSpeciesFluxes = numpy.zeros(species.speciesCounter, float)

		maxRelativeNetworkLeakFluxes = numpy.zeros(len(model.unirxnNetworks), float)

		endtime = 10.0 # default. check for user value:
		for target in model.termination:
			if target.__class__ == modelmodule.TerminationTime:
				endtime = target.time

		# Set up initial conditions
		P = gas.pressure()
		V = sim.reactors()[0].volume()
		T = gas.temperature()
		# Note that, for molar density, Cantera thinks in kmol/m^3, while
		# RMG thinks in mol/m^3
		Ni = gas.molarDensity()*1000.0 * gas.moleFractions() * V
		y = [P, V, T]; y.extend(Ni)
		Ni0 = Ni
		y0 = y

#		# Output information about simulation at current time
		header = 'Time          '
		for target in model.termination:
			if target.__class__ == modelmodule.TerminationConversion: header += 'Conv        '
		header += 'Char. flux     Max. rel. flux to edge'
		logging.debug(header)
#		self.printSimulationStatus(model, 0, y, y0, criticalFlux, maxSpeciesFlux, maxSpecies)
#		tlist.append(0.0); ylist.append(y0)
#		dydtlist.append(self.getResidual(0.0, y0, model, stoichiometry))

		done = False
		time = 0.0
		while not done:

			# Conduct integration
			# Uses a fixed (on a log scale) time step
			nexttime = min(endtime,time*1.2589254117941673)
			# advance cantera one step
			if sim.step(endtime) < endtime:
				# didn't get to endtime, so take another step
				if sim.step(endtime) < nexttime:
					# still didn't get to endtime, so advance to nextime
					sim.advance(nexttime)
#			# Uses the same time steps that the Cantera solver used
#			sim.step(endtime)

			# Get state at current time
			time = sim.time()
			P = gas.pressure()
			V = sim.reactors()[0].volume()
			T = gas.temperature()
			# Note that, for molar density, Cantera thinks in kmol/m^3, while
			# RMG thinks in mol/m^3
			Ni = gas.molarDensity()*1000.0 * gas.moleFractions() * V
			y = [P, V, T]; y.extend(Ni)

			# Calculate species fluxes of all core and edge species at the
			# current time
			rxnRates = self.getReactionRates(P, V, T, Ni, model)
			dNidt = stoichiometry * rxnRates
			
			# Determine characteristic species flux
			charFlux = self.getCharacteristicFlux(model, stoichiometry, rxnRates)

			# Store the highest relative flux for each species
			speciesList = model.core.species[:]; speciesList.extend(model.edge.species)
			for spec in speciesList:
				i = spec.id - 1
				if maxRelativeSpeciesFluxes[i] < abs(dNidt[i])/charFlux:
					maxRelativeSpeciesFluxes[i] = abs(dNidt[i])/charFlux

			# Test for model validity
			criticalFlux = charFlux * model.fluxToleranceInterrupt
			edgeValid, maxSpecies, maxSpeciesFlux = self.isModelValid(model, dNidt, criticalFlux)

			# Test leak fluxes of unimolecular networks
			if settings.unimolecularReactionNetworks:
				maxNetwork = None; maxNetworkFlux = 0.0
				# Get current leak fluxes of all unimolecular reaction networks
				networkLeakFluxes = self.getNetworkLeakFluxes(model, P, V, T, Ni, criticalFlux)
				for i in range(len(networkLeakFluxes)):
					if maxRelativeNetworkLeakFluxes[i] < abs(networkLeakFluxes[i]) / criticalFlux:
						maxRelativeNetworkLeakFluxes[i] = abs(networkLeakFluxes[i]) / criticalFlux
					if networkLeakFluxes[i] > maxNetworkFlux or maxNetwork is None:
						maxNetwork = model.unirxnNetworks[i]
						maxNetworkFlux = networkLeakFluxes[i]
				networksValid = (maxNetworkFlux <= criticalFlux)

			else:
				networksValid = True
				maxNetwork = None
				maxNetworkFlux = 0.0

			# Output information about simulation at current time
			self.printSimulationStatus(model, time, y, y0, charFlux, maxSpeciesFlux/charFlux, maxSpecies)
			tlist.append(time); ylist.append(y)
			#dydtlist.append(self.getResidual(time, solver.y, model, stoichiometry))

			# Exit simulation if model is not valid (exceeds interruption criterion)
			if not edgeValid or not networksValid:
				print gas
				logging.info('')
				# Choose the item with the maximum flux and act on it
				if maxSpeciesFlux >= maxNetworkFlux:
					logging.info('At t = %s, an edge species flux exceeds the critical flux for simulation interruption' % (time))
					logging.info('\tCharacteristic flux: %s' % (charFlux))
					logging.info('\tCritical flux: %s (%s times charFlux)' % (criticalFlux, model.fluxToleranceInterrupt))
					logging.info('\tSpecies flux for %s: %s (%.2g times charFlux)' % (maxSpecies, maxSpeciesFlux, maxSpeciesFlux/charFlux))
					return tlist, ylist, dydtlist, False, maxSpecies
				else:
					logging.info('At t = %s, a network leak flux exceeds the critical flux for simulation interruption' % (time))
					logging.info('\tCharacteristic flux: %s' % (charFlux))
					logging.info('\tCritical flux: %s (%s times charFlux)' % (criticalFlux, model.fluxToleranceInterrupt))
					logging.info('\tNetwork leak flux for %s: %s (%.2g times charFlux)' % (maxNetwork, maxNetworkFlux, maxNetworkFlux/charFlux))
					return tlist, ylist, dydtlist, False, maxNetwork

			# Test for simulation completion
			for target in model.termination:
				if target.__class__ == modelmodule.TerminationConversion:
					index = model.core.species.index(target.species) + 3
					conversion = 1.0 - y[index] / y0[index]
					if conversion > target.conversion: done = True
				elif target.__class__ == modelmodule.TerminationTime:
					if time > target.time: done = True

		logging.info(str(gas))

		# Compare maximum species fluxes
		maxRelativeSpeciesFlux = 0.0; maxSpecies = None
		speciesToRemove = []; maxRelativeFluxes_dict = {}
		for spec in model.edge.species:
			i = spec.id - 1
			# pick out the single highest-flux edge species
			if maxRelativeSpeciesFluxes[i] > maxRelativeSpeciesFlux:
				maxRelativeSpeciesFlux = maxRelativeSpeciesFluxes[i]
				maxSpecies = spec
			# mark for removal those species whose flux is always too low
			if maxRelativeSpeciesFluxes[i] < model.fluxToleranceKeepInEdge:
				speciesToRemove.append(spec)
			# put max relative flux in dictionary
			maxRelativeFluxes_dict[spec] = maxRelativeSpeciesFluxes[i]
		edgeValid = maxRelativeSpeciesFlux <= model.fluxToleranceMoveToCore

		# Compare maximum network leak fluxes
		maxRelativeNetworkLeakFlux = 0.0; maxNetwork = None
		if settings.unimolecularReactionNetworks:
			# Compare maximum species fluxes
			for i in range(len(model.unirxnNetworks)):
				# pick out the single highest-flux edge species
				if maxRelativeNetworkLeakFluxes[i] > maxRelativeNetworkLeakFlux or maxNetwork is None:
					maxRelativeNetworkLeakFlux = maxRelativeNetworkLeakFluxes[i]
					maxNetwork = model.unirxnNetworks[i]
			networksValid = maxRelativeNetworkLeakFlux < model.fluxToleranceMoveToCore
		else:
			networksValid = True

		def removalSortKey(sp):
			return maxRelativeFluxes_dict[sp]
		speciesToRemove.sort(key=removalSortKey)

		# If model is not valid at these criteria, then return
		if not edgeValid or not networksValid:

			# trim the edge according to fluxToleranceKeepInEdge
			logging.info("Removing from edge %d/%d species whose relative flux never exceeded %s"%(
				len(speciesToRemove),len(model.edge.species),model.fluxToleranceKeepInEdge ) )
			logging.info("Max. rel. flux.\tSpecies")
			for sp in speciesToRemove:
				logging.info("%-10.3g    \t%s"%(maxRelativeFluxes_dict[sp], sp))
				model.removeSpeciesFromEdge(sp)

			# trim the edge according to maximumEdgeSpecies
			if len(model.edge.species)> model.maximumEdgeSpecies:
				logging.info("Removing from edge %d/%d species to reach maximum edge size of %s species"%(
					len(model.edge.species)-model.maximumEdgeSpecies,
					len(model.edge.species),
					model.maximumEdgeSpecies ) )
				edgeSpeciesCopy = model.edge.species[:]
				edgeSpeciesCopy.sort(key=removalSortKey)
				logging.info("Max. rel. flux.\tSpecies")
				while len(model.edge.species)>model.maximumEdgeSpecies:
					sp = edgeSpeciesCopy.pop(0)
					logging.info("%-10.3g    \t%s"%(maxRelativeFluxes_dict[sp], sp))
					model.removeSpeciesFromEdge(sp)

			criticalFlux = charFlux * model.fluxToleranceMoveToCore
			print gas
			logging.info('')
			# Choose the item with the maximum flux and act on it
			if maxSpeciesFlux >= maxNetworkFlux:
				logging.info('At some time the species flux for %s exceeded the critical flux\nrelative to the characteristic core flux at that time' % (maxSpecies))
				logging.info('\tCharacteristic flux: %s' % (charFlux))
				logging.info('\tCritical flux: %s (%s times charFlux)' % (criticalFlux, model.fluxToleranceMoveToCore))
				logging.info('\tSpecies flux for %s: %s (%.2g times charFlux)' % (maxSpecies, maxSpeciesFlux, maxSpeciesFlux/charFlux))
				return tlist, ylist, dydtlist, False, maxSpecies
			else:
				logging.info('At some time the network leak flux for %s exceeded the critical flux\nrelative to the characteristic core flux at that time' % (maxNetwork))
				logging.info('\tCharacteristic flux: %s' % (charFlux))
				logging.info('\tCritical flux: %s (%s times charFlux)' % (criticalFlux, model.fluxToleranceMoveToCore))
				logging.info('\tNetwork leak flux for %s: %s (%.2g times charFlux)' % (maxNetwork, maxNetworkFlux, maxNetworkFlux/charFlux))
				return tlist, ylist, dydtlist, False, maxNetwork

		return tlist, ylist, dydtlist, True, None

	def printSimulationStatus(self, model, t, y, y0, charFlux, maxSpeciesFlux, maxSpecies):
		"""
		Log a line of text describing the current status of the simulation. The
		information logged is the current time `t`, the current conversion of
		all species being monitored for conversion targets, the characteristic
		flux `charFlux`, the maximum species flux `maxSpeciesFlux`, and the
		species with that flux `maxSpecies`.
		"""

		# Output information about simulation at current time
		status = '%8.4e' % (t)
		for target in model.termination:
			if target.__class__ == modelmodule.TerminationConversion:
				index = model.core.species.index(target.species) + 3
				conversion = 1.0 - y[index] / y0[index]
				status += '    %8.4g' % (conversion)
		status += '    %8.4e    %8.4g  %s' % (charFlux, maxSpeciesFlux, maxSpecies)
		logging.verbose(status)

	def postprocess(self, model, t, y, dydt, label=''):
		"""
		Postprocess the results of a simulation. This function generates a
		number of plots: temperature, pressure, volume, and concentration
		profiles. The parameters are the reaction `model`, the list of times
		`t`, the list of state vectors `y`, and an optional `label` for the
		reaction system.
		"""
		# Only do if the option for plot generation has been set
		if not settings.generatePlots:
			return
		import pylab
		# Reshape y into a matrix rather than a list of lists
		y0 = numpy.zeros((len(t), len(y[0])), float)
		for i, u in enumerate(y):
			for j, v in enumerate(u):
				y0[i,j] = v
		## Reshape dydt into a matrix rather than a list of lists
		#dydt0 = numpy.zeros((len(t), len(dydt[0])), float)
		#for i, u in enumerate(dydt):
		#	for j, v in enumerate(u):
		#		dydt0[i,j] = v

		# Create the legend for the concentration profile
		legend = []
		for spec in model.core.species:
			legend.append(str(spec))

		# Make pressure plot and save to file
		pylab.semilogx(t[1:], y0[1:,0])
		pylab.xlabel('Time (s)')
		pylab.ylabel('Pressure (Pa)')
		pylab.title('Pressure profile for reaction system ' + label)
		pylab.savefig(settings.outputDirectory + '/plot/pressureProfile' + label + '.svg')
		pylab.clf()

		# Make volume plot and save to file
		pylab.semilogx(t[1:], y0[1:,1])
		pylab.xlabel('Time (s)')
		pylab.ylabel('Volume (m^3)')
		pylab.title('Volume profile for reaction system ' + label)
		pylab.savefig(settings.outputDirectory + '/plot/volumeProfile' + label + '.svg')
		pylab.clf()

		# Make temperature plot and save to file
		pylab.semilogx(t[1:], y0[1:,2])
		pylab.xlabel('Time (s)')
		pylab.ylabel('Temperature (K)')
		pylab.title('Temperature profile for reaction system ' + label)
		pylab.savefig(settings.outputDirectory + '/plot/temperatureProfile' + label + '.svg')
		pylab.clf()

		# Make concentration plot and save to file
		pylab.loglog(t[1:], y0[1:,3:])
		pylab.xlabel('Time (s)')
		pylab.ylabel('Concentration (mol/m^3)')
		pylab.title('Concentration profiles for reaction system ' + label)
		pylab.legend(legend)
		pylab.savefig(settings.outputDirectory + '/plot/concentrationProfile' + label + '.svg')
		pylab.clf()

## Make species flux plot and save to file
		#try:
		#	pylab.loglog(t[1:], abs(dydt0[1:,3:len(model.core.species)+3]))
		#	pylab.xlabel('Time (s)')
		#	pylab.ylabel('Species flux (mol/m^3*s)')
		#	pylab.title('Species flux profiles for reaction system ' + label)
		#	pylab.legend(legend)
		#	pylab.savefig(settings.outputDirectory + '/plot/fluxProfile' + label + '.svg')
		#except OverflowError:
		#	pass

		pylab.clf()

	def getReactionRates(self, P, V, T, Ni, model):
		"""
		Evaluate the reaction rates for all reactions in the model (core and
		edge).
		"""
		Ci = {}
		for i, spec in enumerate(model.core.species):
			Ci[spec] = Ni[i] / V

		return model.getReactionRates(T, P, Ci)

	def getCharacteristicFlux(self, model, stoichiometry, rxnRates0):
		"""
		Determine the characteristic flux for the system given a `model` at
		the specified pressure `P`, volume `V`, temperature `T`, and numbers
		of moles `Ni`. The `stoichiometry` parameter is the stoichiometry
		matrix for the model. The characteristic flux is the root mean square
		of the *net* flux in mol/m^3*s of each core species due to the core
		reactions only.
		"""

		# Generate core species fluxes based on core reactions only
		rxnRates = numpy.zeros(rxnRates0.shape, numpy.float64)
		for rxn in model.core.reactions:
			j = rxn.id - 1
			rxnRates[j] = rxnRates0[j]
		speciesRates = stoichiometry * rxnRates

		# Characteristic flux is root mean square of core species fluxes
		charFlux = 0.0
		for spec in model.core.species:
			flux = speciesRates[spec.id - 1]
			charFlux += flux * flux
		return math.sqrt(charFlux)

	def isModelValid(self, model, dNidt, criticalFlux):
		"""
		Returns :data:`True` if `model` is valid given the set of species fluxes
		`dNidt` and the critical flux `criticalFlux`.
		Also returns the edge species whose flux is greatest, and that flux.
		"""

		if len(model.edge.species) == 0:
			# Nothing in edge - model must be valid!
			return True, None, 0.0

		maxSpecies = None
		maxSpeciesFlux = 0.0
		for spec in model.edge.species:
			j = spec.id - 1
			if dNidt[j] > maxSpeciesFlux:
				maxSpecies = spec
				maxSpeciesFlux = dNidt[j]
		return (maxSpeciesFlux < criticalFlux), maxSpecies, maxSpeciesFlux

		# maxSpeciesFlux, maxSpeciesIndex = max([ (value, i+len(model.core.species)) for i, value in enumerate(dNidt[len(model.core.species):]) ])
		# return (maxSpeciesFlux < criticalFlux), speciesList[maxSpeciesIndex], maxSpeciesFlux

	def getNetworkLeakFluxes(self, model, P, V, T, Ni, criticalFlux):
		"""
		Returns :data:`True` if `model` is valid given the set of species fluxes
		`dNidt` and the critical flux `criticalFlux`.
		Also returns the edge species whose flux is greatest, and that flux.
		"""

		import rmg.reaction as reaction

		conc = {}; totalConc = 0.0
		for i, spec in enumerate(model.core.species):
			conc[spec] = Ni[i] / V
			totalConc += conc[spec]

		leakFluxes = numpy.zeros(len(model.unirxnNetworks), numpy.float64)
		for i, network in enumerate(model.unirxnNetworks):
			leakFluxes[i] = network.getLeakFlux(T, P, conc, totalConc)
		return leakFluxes

