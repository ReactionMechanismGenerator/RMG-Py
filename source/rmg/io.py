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

import model
import chem
import data
import species
import reaction

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
					species.forbiddenStructures = data.Dictionary()
					species.forbiddenStructures.load(database[2] + '/forbiddenStructure.txt')
					species.forbiddenStructures.toStructure()
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
			structure = chem.Structure()
			
			cml = getFirstChildElement(element, 'cml')
			inchi = getFirstChildElement(element, 'inchi')
			smiles = getFirstChildElement(element, 'smiles')
			if cml is not None:
				cmlstr = str(getFirstChildElement(cml, 'molecule').toxml())
				structure.fromCML(cmlstr)
			elif inchi is not None:
				inchistr = str(getElementText(inchi))
				structure.fromInChI(inchistr)
			elif smiles is not None:
				smilesstr = str(getElementText(smiles))
				structure.fromSMILES(smilesstr)
			else:
				raise InvalidInputFileException('Species '+label+' missing structure information.' )
			
			# Create a new species and append the species to the core
			spec = species.makeNewSpecies(structure, label, reactive)
			coreSpecies.append(spec)
		
			# Add to local species dictionary (for matching with other parts of file)
			speciesDict[sid] = spec

		logging.debug('')
		
		# Create an empty reaction model
		reactionModel = model.CoreEdgeReactionModel()
		
		# Read model flux tolerance
		element = getFirstChildElement(root, 'fluxTolerance')
		reactionModel.fluxTolerance = float(getElementText(element))
		logging.debug('Model flux tolerance set to %s' % (reactionModel.fluxTolerance))
		logging.debug('')

		# Read dynamic simulator
		element = getFirstChildElement(root, 'simulator')
		model.absoluteTolerance = float(element.getAttribute('atol'))
		model.relativeTolerance = float(element.getAttribute('rtol'))
		logging.info('Read dynamic simulator')
		logging.debug('Simulator:')
		logging.debug('\tAbsolute tolerance set to %s' % (reactionModel.absoluteTolerance))
		logging.debug('\tRelative tolerance set to %s' % (reactionModel.relativeTolerance))
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
			
			else:
				raise InvalidInputFileException('Invalid physical property model type "' + propModelType + '".')


			# Initialize all initial concentrations to zero
			for spec in coreSpecies:
				reactionSystem.initialConcentration[spec] = 0.0
			
			# List of initial concentrations
			concentrations = getElements(element, 'concentration')
			for concentration in concentrations:
			
				# Read the concentration from the file
				value = float(getElementText(concentration))
				sid = concentration.getAttribute('speciesID')
				units = str(concentration.getAttribute('units'))
				C = pq.Quantity(value, units); C = float(C.simplified)

				reactionSystem.initialConcentration[speciesDict[sid]] = C
			
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
			
	except InvalidInputFileException, e:
		logging.exception(str(e))
	except IOError, e:
		logging.exception('Input file "' + e.filename + '" not found.')
	finally:
		# Unlink the DOM tree when finished
		dom.unlink()
		
	return reactionModel, coreSpecies, reactionSystems

################################################################################

