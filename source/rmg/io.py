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

import xml.dom.minidom
import quantities as pq
import logging
import os

import constants
import settings
import model
import structure
import data
import species
import reaction
import thermo

"""
Contains functions for manipulation of RMG input and output files.
"""

################################################################################

class InvalidInputFileException(Exception):
	"""
	An exception used when parsing an RMG input file to indicate that the input
	file is invalid. The `msg` parameter is used to specify what about the file
	caused the exception to be raised.
	"""	

	def __init__(self, msg):
		self.msg = msg
	
	def __str__(self):
		return 'Invalid XML for RMG input file: ' + self.msg

################################################################################

def getFirstChildElement(parent, name):
	"""
	Returns the first child element of the XML `parent` element with a matching
	`name`. The `parent` parameter comes from a node of a DOM tree as produced
	via the :data:`xml.dom.minidom` module.
	"""
	elements = getElements(parent, name)
	if len(elements) > 0:
		return elements[0]
	else:
		return None

def getElements(parent, name):
	"""
	Returns a list of all child elements of the XML `parent` element with a 
	matching `name`. The `parent` parameter comes from a node of a DOM tree as 
	produced via the :data:`xml.dom.minidom` module.
	"""
	return parent.getElementsByTagName(name)

def getElementText(element):
	"""
	Returns the string text that is found between the open and close tags of
	`element`, or an empty string if no text is found. The `element` parameter 
	comes from a node of a DOM tree as produced via the 
	:data:`xml.dom.minidom` module.
	"""
	for child in element.childNodes:
		if child.nodeType == xml.dom.Node.TEXT_NODE:
			return child.data
	return ''

################################################################################

def readInputFile(fstr):
	"""
	Parse an RMG input file at the location `fstr`. If successful, this 
	function returns a :class:`rmg.model.CoreEdgeReactionModel` object and a 
	list of one or more :class:`rmg.model.ReactionSystem` objects.
	"""

	try:
		
		# Parse the RMG input XML file into a DOM tree
		dom = xml.dom.minidom.parse(fstr)
		
		
		# Process root element (must be a rmginput element)
		root = dom.documentElement
		if root.tagName != 'rmginput':
			raise InvalidInputFileException('Incorrect root element. Should be <rmginput>')
		
		# Process units
		element = getFirstChildElement(root, 'units')
		if element is None: units = 'si'
		else: units = getElementText(element)
		pq.set_default_units(units)

		# Read draw molecules option
		element = getFirstChildElement(root, 'drawMolecules')
		drawMolecules = getElementText(element).lower()
		if drawMolecules == 'on' or drawMolecules == 'true' or drawMolecules == 'yes':
			settings.drawMolecules = True
		else:
			settings.drawMolecules = False
		
		
		# Read generate plots option
		element = getFirstChildElement(root, 'generatePlots')
		generatePlots = getElementText(element).lower()
		if generatePlots == 'on' or generatePlots == 'true' or generatePlots == 'yes':
			settings.generatePlots = True
		else:
			settings.generatePlots = False
		
		
		# Process databases
		databases = []
		elements = getElements(root, 'database')
		for element in elements:
			
			# Get database type
			databaseType = element.getAttribute('type').lower()
			if databaseType == '': databaseType = 'general'
			if databaseType != 'general':
				raise InvalidInputFileException('Invalid database type "' + databaseType + '".')
			
			# Get database name and form path
			databaseName = getElementText(element)
			databasePath = os.path.dirname(__file__)
			databasePath = os.path.abspath(databasePath + '/../../data/' + databaseName)
			if not os.path.exists(databasePath):
				raise InvalidInputFileException('Database "' + databaseName + '" not found.')
			
			databases.append([databaseName, databaseType, databasePath])
		
		# Output info about databases
		if len(databases) == 1:
			logging.info('Found ' + str(len(databases)) + ' database')
		else:
			logging.info('Found ' + str(len(databases)) + ' databases')
		
		# Load databases
		generalDatabaseCount = 0
		for database in databases:
			if database[1] == 'general':
				if generalDatabaseCount == 0:
					logging.debug('General database: ' + database[2])
					# Load thermo databases
					species.thermoDatabase = species.ThermoDatabaseSet()
					species.thermoDatabase.load(database[2] + '/')
					# Load forbidden structures
					thermo.forbiddenStructures = data.Dictionary()
					thermo.forbiddenStructures.load(database[2] + '/forbiddenStructure.txt')
					thermo.forbiddenStructures.toStructure()
					# Load kinetic databases (reaction families)
					reaction.kineticsDatabase = reaction.ReactionFamilySet()
					reaction.kineticsDatabase.load(database[2] + '/')
				generalDatabaseCount += 1
				
		logging.debug('')
		
		# Check that exactly one general database was specified
		if generalDatabaseCount == 0:
			raise InvalidInputFileException('No general database specified; one must be present.')
		elif generalDatabaseCount > 1:
			raise InvalidInputFileException('Multiple general databases specified; only one is allowed.')
	
		# Process species
		coreSpecies = []; speciesDict = {}
		elements = getElements(root, 'species')
		logging.info('Found ' + str(len(elements)) + ' species')
		for element in elements:
			
			# Attributes of the species element
			sid = element.getAttribute('id')
			label = element.getAttribute('label')
			reactive = element.getAttribute('reactive').lower()
			if reactive == 'no' or reactive == 'false' or reactive == 'n':
				reactive = False
			else:
				reactive = True		
			
			# Load structure
			struct = structure.Structure()
			
			cml = getFirstChildElement(element, 'cml')
			inchi = getFirstChildElement(element, 'inchi')
			smiles = getFirstChildElement(element, 'smiles')
			if cml is not None:
				cmlstr = str(getFirstChildElement(cml, 'molecule').toxml())
				struct.fromCML(cmlstr)
			elif inchi is not None:
				inchistr = str(getElementText(inchi))
				struct.fromInChI(inchistr)
			elif smiles is not None:
				smilesstr = str(getElementText(smiles))
				struct.fromSMILES(smilesstr)
			else:
				raise InvalidInputFileException('Species '+label+' missing structure information.' )
			
			# Create a new species and append the species to the core
			spec = species.makeNewSpecies(struct, label, reactive)
			coreSpecies.append(spec)
			
			# Add to local species dictionary (for matching with other parts of file)
			speciesDict[sid] = spec

		logging.debug('')
		
		# Create an empty reaction model
		reactionModel = model.CoreEdgeReactionModel()
		
		# Read model flux tolerance
		fluxtolerance = getFirstChildElement(root, 'fluxTolerance')
		element = getFirstChildElement(fluxtolerance, 'keepInEdge')
		reactionModel.fluxTolerance_keepInEdge = float(getElementText(element))
		element = getFirstChildElement(fluxtolerance, 'moveToCore')
		reactionModel.fluxTolerance_moveToCore = float(getElementText(element))
		element = getFirstChildElement(fluxtolerance, 'interruptSimulation')
		reactionModel.fluxTolerance_interruptSimulation = float(getElementText(element))
		
		logging.debug('Model flux tolerances set to:')
		logging.debug('\tKeep in edge:         %s' % (reactionModel.fluxTolerance_keepInEdge) )
		logging.debug('\tMove to core:         %s' % (reactionModel.fluxTolerance_moveToCore) )
		logging.debug('\tInterrupt simulation: %s' % (reactionModel.fluxTolerance_interruptSimulation) )
		logging.debug('')
		
		# Read maximum model size
		maxmodelsize = getFirstChildElement(root, 'maximumModelSize')
		if maxmodelsize is None:
			logging.debug('Maximum model size is not set')
		else:
			element = getFirstChildElement(maxmodelsize, 'edgeSpecies')
			reactionModel.maxModelSize_EdgeSpecies = int(getElementText(element))
			logging.debug('Maximum model size set to:')
			logging.debug('\tEdge species:         %s' % (reactionModel.maxModelSize_EdgeSpecies) )
		logging.debug('')

		# Read dynamic simulator
		element = getFirstChildElement(root, 'simulator')
		reactionModel.absoluteTolerance = float(element.getAttribute('atol'))
		reactionModel.relativeTolerance = float(element.getAttribute('rtol'))
		logging.info('Read dynamic simulator')
		logging.debug('Simulator:')
		logging.debug('\tAbsolute tolerance set to %s' % (reactionModel.absoluteTolerance))
		logging.debug('\tRelative tolerance set to %s' % (reactionModel.relativeTolerance))
		logging.debug('')

		# Read termination targets
		termination = getFirstChildElement(root, 'termination')
		elements = getElements(termination, 'target')
		for element in elements:

			targetType = element.getAttribute('type')
			if targetType == 'conversion':
				sid = element.getAttribute('speciesID')
				spec = speciesDict[sid]
				conv = float(getElementText(element))
				if conv < 0.0 or conv > 1.0:
					raise InvalidInputFileException('Invalid value for termination fractional conversion.')
				reactionModel.termination.append(model.TerminationConversion(spec, conv))
			elif targetType == 'time':
				units = str(element.getAttribute('units'))
				time = float(getElementText(element))
				time = pq.Quantity(time, units); time = float(time.simplified)
				if time < 0.0:
					raise InvalidInputFileException('Invalid value for termination time.')
				reactionModel.termination.append(model.TerminationTime(time))
			else:
				raise InvalidInputFileException('Invalid termination target type "'+targetType+'".')
		if len(reactionModel.termination) == 0:
			raise InvalidInputFileException('No termination targets specified.')

		# Output info about termination targets
		if len(reactionModel.termination) == 1:
			logging.info('Found ' + str(len(reactionModel.termination)) + ' termination target')
		else:
			logging.info('Found ' + str(len(reactionModel.termination)) + ' termination targets')
		for index, target in enumerate(reactionModel.termination):
			string = '\tTermination target #' + str(index+1) + ': '
			if target.__class__ == model.TerminationConversion:
				string += 'conversion ' + str(target.species) + ' ' + str(target.conversion)
			elif target.__class__ == model.TerminationTime:
				string += 'time ' + str(target.time)
			logging.debug(string)	
				
		logging.debug('')

		# Process reaction systems
		reactionSystems = []
		elements = getElements(root, 'reactionSystem')
		for element in elements:
		
			# Create a new reaction system
			rsType = element.getAttribute('type')
			if rsType == 'batch':
				reactionSystem = model.BatchReactor()
			else:
				raise InvalidInputFileException('Invalid reaction system type "' + rsType + '".')
			
			# Temperature model
			temperatureModel = getFirstChildElement(element, 'temperatureModel')
			tempModelType = temperatureModel.getAttribute('type')
			if tempModelType == 'isothermal':
				
				# Read the (constant) temperature from the file
				temperature = getFirstChildElement(temperatureModel, 'temperature')
				value = float(getElementText(temperature))
				units = str(temperature.getAttribute('units'))
				T = pq.Quantity(value, units); T = float(T.simplified)
				
				# Set the reaction system's temperature model to isothermal
				reactionSystem.temperatureModel = model.TemperatureModel()
				reactionSystem.temperatureModel.setIsothermal(T)
				
			else:
				raise InvalidInputFileException('Invalid temperature model type "' + tempModelType + '".')
			
			# Pressure model
			pressureModel = getFirstChildElement(element, 'pressureModel')
			pressModelType = pressureModel.getAttribute('type')
			if pressModelType == 'isobaric':
				
				# Read the (constant) pressure from the file
				pressure = getFirstChildElement(pressureModel, 'pressure')
				value = float(getElementText(pressure))
				units = str(pressure.getAttribute('units'))
				P = pq.Quantity(value, units); P = float(P.simplified)
				
				# Set the reaction system's pressure model to isobaric
				reactionSystem.pressureModel = model.PressureModel()
				reactionSystem.pressureModel.setIsobaric(P)
				
			else:
				raise InvalidInputFileException('Invalid pressure model type "' + pressModelType + '".')

			# Physical property model
			propModel = getFirstChildElement(element, 'physicalPropertyModel')
			propModelType = propModel.getAttribute('type')
			if propModelType.lower() == 'idealgas':
				# Set the reaction system's pressure model to isobaric
				reactionSystem.equationOfState = model.IdealGas()
			elif propModelType.lower() == 'incompressibleliquid':
				molarVolume = getFirstChildElement(element, 'molarVolume')
				value = float(getElementText(molarVolume))
				units = str(molarVolume.getAttribute('units'))
				Vmol = float(pq.Quantity(value, units).simplified); 
				
				reactionSystem.equationOfState = model.IncompressibleLiquid( 
					P = reactionSystem.pressureModel.getPressure(),
					T = reactionSystem.temperatureModel.getTemperature(),
					Vmol = Vmol
					)
			else:
				raise InvalidInputFileException('Invalid physical property model type "' + propModelType + '".')

			# Get total concentration
			T = reactionSystem.temperatureModel.getTemperature(0)
			P = reactionSystem.pressureModel.getPressure(0)
			totalConc = 1.0 / reactionSystem.equationOfState.getVolume(T, P, [1.0])

			# Initialize all initial concentrations to zero
			for spec in coreSpecies:
				reactionSystem.initialConcentration[spec] = 0.0
			
			# List of initial concentrations
			moleFractions = getElements(element, 'moleFraction')
			for moleFraction in moleFractions:
			
				# Read the concentration from the file
				value = float(getElementText(moleFraction))
				sid = moleFraction.getAttribute('speciesID')
				
				reactionSystem.initialConcentration[speciesDict[sid]] = value * totalConc
			
			# Append to list of reaction systems
			reactionSystems.append(reactionSystem)
		
		# Output info about reaction system
		if len(reactionSystems) == 1:
			logging.info('Found ' + str(len(reactionSystems)) + ' reaction system')
		else:
			logging.info('Found ' + str(len(reactionSystems)) + ' reaction systems')
		for index, reactionSystem in enumerate(reactionSystems):
			logging.debug('Reaction system #' + str(index+1) + ':')
			logging.debug('\t' + str(reactionSystem.temperatureModel))
			logging.debug('\t' + str(reactionSystem.pressureModel))
			for spec, conc in reactionSystem.initialConcentration.iteritems():
				if spec.reactive:
					logging.debug('\tInitial concentration of ' + str(spec) + ': ' + str(conc))
				else:
					logging.debug('\tConstant concentration of ' + str(spec) + ': ' + str(conc))
				
				
		logging.debug('')
			
		# Unlink the DOM tree when finished
		dom.unlink()
	except InvalidInputFileException, e:
		logging.exception(str(e))
		raise e
	except IOError, e:
		logging.exception('Input file "' + e.filename + '" not found.')
		raise e
	except xml.parsers.expat.ExpatError, e:
		logging.exception('Invalid XML file: '+e.message+'\n')
		raise InvalidInputFileException('Invalid XML file: '+e.message)
		
	return reactionModel, coreSpecies, reactionSystems

################################################################################

def writeOutputFile(fstr, reactionModel, reactionSystems):
	"""
	Write an RMG output file at the location `fstr` based on the provided
	`reactionModel` and `reactionSystems`.
	"""

	# Create document
	dom = xml.dom.minidom.Document()

	# Process root element
	root = dom.createElement('rmgoutput')
	dom.appendChild(root)

	# Process core species list
	speciesList = dom.createElement('speciesList')
	root.appendChild(speciesList)
	for spec in reactionModel.core.species:

		element = dom.createElement('species')
		element.setAttribute('id', str(spec.id))
		element.setAttribute('label', spec.label)
		element.setAttribute('reactive', 'yes' if spec.reactive else 'no')
		speciesList.appendChild(element)

		# Output the structure using CML
		cml = dom.createElement('cml')
		element.appendChild(cml)
		dom0 = xml.dom.minidom.parseString(spec.toCML())
		cml.appendChild(dom0.documentElement)

		# Output the thermo data
		spec.thermoData.toXML(dom, element)

	# Process core reactions list
	reactionList = dom.createElement('reactionList')
	root.appendChild(reactionList)
	for rxn in reactionModel.core.reactions:

		element = dom.createElement('reaction')
		element.setAttribute('family', rxn.family.label)
		reactionList.appendChild(element)

		for reac in rxn.reactants:
			reactant = dom.createElement('reactant')
			reactant.setAttribute('id', str(reac.id))
			element.appendChild(reactant)
		for prod in rxn.products:
			product = dom.createElement('product')
			product.setAttribute('id', str(prod.id))
			element.appendChild(product)

		kinetics = dom.createElement('kinetics')
		element.appendChild(kinetics)
		for k in rxn.kinetics:
			k.toXML(dom, kinetics)

	# Write output file
	f = open(fstr, 'w')
	f.write('\n'.join([l for l in dom.toprettyxml().split('\n') if l.strip()]))
	f.close()

	# Print to log
	logging.info('')
	logging.info('Output written to ' + fstr)
	
