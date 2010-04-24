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
import log as logging
import os
import math

import constants
import settings
import model
import structure
import data
import species
import reaction
import thermo.data
import spectral.data

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

class InvalidXMLError(Exception):
	"""
	An exception to be raised when invalid XML is encountered. The
	exception should be passed a string `msg` containing details about the
	error.
	"""

	def __init__(self, msg=''):
		self.msg = msg

	def __str__(self):
		return 'Invalid XML: %s' % self.msg

################################################################################

def loadKineticsDatabase(dstr):
	"""
	Load the RMG kinetics database located at `dstr` into the global variable
	`rmg.reaction.kineticsDatabase`.
	"""
	reaction.kineticsDatabase = reaction.ReactionFamilySet()
	reaction.kineticsDatabase.load(dstr)

def loadFrequencyDatabase(dstr):
	"""
	Load the RMG thermo database located at `dstr` into the global variable
	`rmg.spectral.data.frequencyDatabase`.
	"""
	frequenciesDatabasePath = os.path.join(dstr,'frequencies_groups')
	spectral.data.frequencyDatabase = spectral.data.FrequencyDatabase()
	spectral.data.frequencyDatabase = spectral.data.loadFrequencyDatabase(frequenciesDatabasePath)

################################################################################

class XML:
	"""
	A class for manipulating an XML DOM tree as created by the xml.dom.minidom
	package. The DOM is stored in the `document` attribute.
	"""

	def __init__(self, document=None, path=''):
		self.document = document
		if path != '': self.load(path)

	def load(self, path):
		"""
		Load a DOM from an XML file located at `path`.
		"""
		self.document = xml.dom.minidom.parse(path)

	def save(self, path):
		"""
		Save a DOM to an XML file located at `path`.
		"""
		f = open(path, 'w')
		self.document.writexml(f, indent='', addindent='\t', newl='\n')
		f.close()

	def cleanup(self):
		"""
		Clean up a DOM tree by unlinking. This should only be done when you
		are finished working with the tree.
		"""
		self.document.unlink()

	def getRootElement(self):
		"""
		Return the root element of the DOM tree.
		"""
		return self.document.documentElement

	def initialize(self, elementName):
		"""
		Initializes a new DOM tree with the root element having the name
		`elementName`. Returns the root element.
		"""
		self.document = xml.dom.minidom.getDOMImplementation().createDocument(None, elementName, None)
		return self.getRootElement()

	def getAttribute(self, element, name, required=True, default=None):
		"""
		Return the value of the attribute `name` of element `root`. If no
		attribute is found and `required` is :data:`True`, an 
		:class:`InvalidXMLError` is raised. If no attribute is found and
		`required` is :data:`False`, `default` is returned.
		"""
		attribute = element.getAttribute(name)
		if attribute != '':
			return attribute
		elif not required:
			return default
		else:
			raise InvalidXMLError('No "%s" attribute found for  <%s> element.' % (name, element.tagName))

	def getChildElement(self, parentElement, name, required=True):
		"""
		Return the child element of `root` with the name `name`. If no
		suitable child element is found, :data:`None` is returned if `required`
		is :data:`False` or an :class:`InvalidXMLError` is raised if `required`
		is :data:`True`.
		"""
		elements = parentElement.getElementsByTagName(name)
		if len(elements) == 1:
			return elements[0]
		elif len(elements) > 1:
			raise InvalidXMLError('Multiple <%s> elements found for <%s> element.' % (name, parentElement.tagName))
		elif not required:
			return None
		else:
			raise InvalidXMLError('No <%s> element found in  <%s> element.' % (name, parentElement.tagName))

	def getChildElements(self, parentElement, name, required=True):
		"""
		Return a list of child elements of `root` with the name `name`. If no
		suitable child elements is found, an an :class:`InvalidXMLError` is
		raised if `required` is :data:`True`.
		"""
		elements = parentElement.getElementsByTagName(name)
		if len(elements) > 0 or not required:
			return elements
		else:
			raise InvalidXMLError('No <%s> elements found in  <%s> element.' % (name, parentElement.tagName))

	def getChildElementText(self, parentElement, name, required=True, default=None):
		"""
		Return the text of the child element of `root` with the name
		`name`. If multiple child elements are found, no text is found, or no
		child element is found and `required` is :data:`True`, an 
		:class:`InvalidXMLError` is raised.
		"""
		element = self.getChildElement(parentElement, name, required)
		if element:
			return self.getElementText(element)
		elif not required:
			return default
		else:
			raise InvalidXMLError('No <%s> element found in  <%s> element.' % (name, parentElement.tagName))

	def getElementText(self, element):
		"""
		Return the text of the element `element`. If no text is found, an
		:class:`InvalidXMLError` is raised.
		"""
		for child in element.childNodes:
			if child.nodeType == xml.dom.Node.TEXT_NODE:
				return child.data
		raise InvalidXMLError('No text found in <%s> element.' % (element.tagName))

	def getQuantity(self, element):
		"""
		Process an element of the form 
		
			<elementName units="units">value</elementName>
			
		and return the corresponding :class:`quantities.Quantity` object.
		"""
		units = str(self.getAttribute(element, 'units', required=False, default=''))
		value = str(self.getElementText(element)).split()
		value = [float(v) for v in value]
		if len(value) == 1: value = value[0]
		return pq.Quantity(value, units)

	def getChildQuantity(self, parentElement, name, required=True, default=0.0):
		"""
		Process an element of the form

			<elementName units="units">value</elementName>

		and return the corresponding :class:`quantities.Quantity` object.
		"""
		element = self.getChildElement(parentElement, name, required)
		if element:
			return self.getQuantity(element)
		else:
			return default
		
	def createElement(self, elementName, parentElement):
		"""
		Create the element `elementName` as an immediate child of `parentElement`
		in the DOM tree `document`. The created element is returned.
		"""
		element = self.document.createElement(elementName)
		parentElement.appendChild(element)
		return element

	def createTextElement(self, elementName, parentElement, text):
		"""
		Create the element `elementName` as an immediate child of `parentElement`
		in the DOM tree `document`. Set the text of the created element to `text`,
		which is passed through a string conversion function along the way. The
		created element is returned.
		"""

		element = self.createElement(elementName, parentElement)
		textNode = self.document.createTextNode(str(text))
		element.appendChild(textNode)
		return element

	def createAttribute(self, attributeName, parentElement, attributeValue):
		"""
		Create the attribute `attributeName` as an attribute of `parentElement`
		with the value `attributeValue` in the DOM tree `document`.
		"""
		parentElement.setAttribute(attributeName, attributeValue)

	def createQuantity(self, elementName, parentElement, value, units=''):
		"""
		Create an element representing a quantity of a certain `value` (scalar
		or list) and `units` (string) with the name `elementName` as a child of
		`parentElement`.
		"""

		if isinstance(value, list):
			quantityElement = self.createTextElement(elementName, parentElement, ' '.join([str(v) for v in value]))
		else:
			quantityElement = self.createTextElement(elementName, parentElement, str(value))
		if units != '':
			self.createAttribute('units', quantityElement, units)

		return quantityElement

################################################################################

def readDatabaseList(xml0, rootElement):
	"""
	Parse the portion of the RMG input file that tells where the RMG database(s)
	are that should be used, i.e. the <databaseList> element.
	"""

	# Process databases
	databases = []
	databaseList = xml0.getChildElement(rootElement, 'databaseList')
	databaseElements = xml0.getChildElements(databaseList, 'database')
	for element in databaseElements:

		# Get database type
		databaseType = xml0.getAttribute(element, 'type', required=True)
		databaseType = databaseType.lower()
		if databaseType not in ['general', 'seedmechanism']:
			raise InvalidInputFileException('Invalid database type "' + databaseType + '"; valid types are "general".')

		# Get database name and form path
		databaseName = xml0.getElementText(element).strip()
		databasePath = os.path.dirname(__file__)
		databasePath = os.path.join(databasePath, '..')
		databasePath = os.path.join(databasePath, '..')
		databasePath = os.path.join(databasePath, 'data')
		databasePath = os.path.join(databasePath, databaseName)
		if not os.path.exists(databasePath):
			raise InvalidInputFileException('Database "%s" not found.' % databaseName)

		databases.append([databaseName, databaseType, databasePath])

	# Output info about databases
	logging.info('Found %s database%s' % (len(databases), 's' if len(databases) > 1 else ''))

	# Check that exactly one general database was specified
	generalDatabaseCount = sum([1 for database in databases if database[1] == 'general'])
	if generalDatabaseCount == 0:
		raise InvalidInputFileException('No general database specified; one must be present.')
	elif generalDatabaseCount > 1:
		raise InvalidInputFileException('Multiple general databases specified; only one is allowed.')

	return databases

def readInputFile(fstr):
	"""
	Parse an RMG input file at the location `fstr`. If successful, this 
	function returns a :class:`rmg.model.CoreEdgeReactionModel` object and a 
	list of one or more :class:`rmg.model.ReactionSystem` objects.
	"""

	try:
		
		# Parse the RMG input XML file into a DOM tree
		xml0 = XML(path=fstr)

		# Make sure root element is a <rmginput> element
		rootElement = xml0.getRootElement()
		if rootElement.tagName != 'rmginput':
			raise InvalidInputFileException('Incorrect root element; should be <rmginput>.')
		
		# Process option list
		optionList = xml0.getChildElement(rootElement, 'optionList')
		# Process restart option
		saveRestart = xml0.getChildElement(optionList, 'saveRestart', required=False)
		if saveRestart:
			frequency = str(xml0.getAttribute(saveRestart, 'frequency', required=False, default='always')).lower()
			if frequency not in ['always', 'hourly', 'daily', 'weekly', 'monthly']:
				raise InvalidInputFileException('Invalid value for attribute frequency in <saveRestart> element; should be one of "always", "hourly", "daily", "weekly", "monthly".')
			iterations = int(xml0.getAttribute(saveRestart, 'iterations', required=False, default='0'))
			settings.saveRestart = [frequency, iterations, 0, 0]
		else:
			settings.saveRestart = None
		# Process units option
		units = xml0.getChildElementText(optionList, 'units', required=False, default='si')
		pq.set_default_units(units)
		# Read draw molecules option
		drawMolecules = xml0.getChildElement(optionList, 'drawMolecules', required=False)
		settings.drawMolecules = (drawMolecules is not None)
		# Read generate plots option
		generatePlots = xml0.getChildElement(optionList, 'generatePlots', required=False)
		settings.generatePlots = (generatePlots is not None)
		# Read spectral data estimation option
		spectralDataEstimation = xml0.getChildElement(optionList, 'spectralDataEstimation', required=False)
		settings.spectralDataEstimation = (spectralDataEstimation is not None)

		# Read unimolecular reaction network option
		unirxnNetworks = xml0.getChildElement(optionList, 'unimolecularReactionNetworks', required=False)
		if unirxnNetworks is not None:

			# Read method
			method = str(xml0.getChildElementText(unirxnNetworks, 'method', required=True))
			allowed = ['modifiedstrongcollision', 'reservoirstate']
			if method.lower() not in allowed:
				raise InvalidInputFileException('Invalid unimolecular reaction networks method "%s"; allowed values are %s.' % (method, allowed))

			# Read grain size
			grainSize = xml0.getChildQuantity(unirxnNetworks, 'grainSize', required=False,
				default=pq.Quantity(0.0, 'J/mol'))
			grainSize = float(grainSize.simplified)
			# Read number of grains
			numberOfGrains = int(xml0.getChildElementText(unirxnNetworks, 'numberOfGrains', required=False, default=0))
			if grainSize == 0.0 and numberOfGrains == 0:
				raise InvalidInputFileException('Must specify a grain size or number of grains for unimolecular reaction networks calculations.')

			# Read interpolation model
			interpolationModel = xml0.getChildElement(unirxnNetworks, 'interpolationModel', required=True)
			modelType = str(xml0.getAttribute(interpolationModel, 'type', required=True))
			allowed = ['none', 'chebyshev', 'pdeparrhenius']
			if modelType.lower() not in allowed:
				raise InvalidInputFileException('Invalid unimolecular reaction networks interpolation model "%s"; allowed values are %s.' % (method, allowed))
			if modelType.lower() == 'chebyshev':
				numTPolys = int(xml0.getChildElementText(interpolationModel, 'numberOfTemperaturePolynomials', required=False, default='4'))
				numPPolys = int(xml0.getChildElementText(interpolationModel, 'numberOfPressurePolynomials', required=False, default='4'))
				interpolationModel = (modelType, numTPolys, numPPolys)
			else:
				interpolationModel = (modelType)

			temperatures = None
			pressures = None
			
			# Read temperature range
			temperatureRange = xml0.getChildElement(unirxnNetworks, 'temperatureRange', required=False)
			if temperatureRange:
				Tnumber = int(xml0.getAttribute(temperatureRange, 'number', required=True))
				units = str(xml0.getAttribute(temperatureRange, 'units', required=True))
				Tmin = float(xml0.getAttribute(temperatureRange, 'min', required=True))
				Tmin = pq.Quantity(Tmin, units)
				Tmin = float(Tmin.simplified)
				Tmax = float(xml0.getAttribute(temperatureRange, 'max', required=True))
				Tmax = pq.Quantity(Tmax, units)
				Tmax = float(Tmax.simplified)

			# Read pressure range
			pressureRange = xml0.getChildElement(unirxnNetworks, 'pressureRange', required=False)
			if pressureRange:
				Pnumber = int(xml0.getAttribute(pressureRange, 'number', required=True))
				units = str(xml0.getAttribute(pressureRange, 'units', required=True))
				Pmin = float(xml0.getAttribute(pressureRange, 'min', required=True))
				Pmin = pq.Quantity(Pmin, units)
				Pmin = float(Pmin.simplified)
				Pmax = float(xml0.getAttribute(pressureRange, 'max', required=True))
				Pmax = pq.Quantity(Pmax, units)
				Pmax = float(Pmax.simplified)

			# Determine actual temperatures to use based on interpolation model
			# For Chebyshev polynomials a Gauss-Chebyshev grid should be used
			# Otherwise an even spacing on a T^-1 basis is used
			if interpolationModel[0].lower() == 'chebyshev':
				# The formula for the Gauss-Chebyshev points was taken from
				# the Chemkin theory manual
				temperatures = [-math.cos((2 * x - 1) * math.pi / (2 * Tnumber)) for x in range(1, Tnumber+1)]
				temperatures = [2.0 / ((1.0/Tmax - 1.0/Tmin) * t + 1.0/Tmax + 1.0/Tmin) for t in temperatures]
			else:
				slope = (1.0/Tmax - 1.0/Tmin) / (Tnumber - 1)
				temperatures = [1.0/(slope * x + 1.0/Tmin) for x in range(Tnumber)]

			# Determine actual pressures to use based on interpolation model
			# For Chebyshev polynomials a Gauss-Chebyshev grid should be used
			# Otherwise an even spacing on a log(P) basis is used
			if interpolationModel[0].lower() == 'chebyshev':
				# The formula for the Gauss-Chebyshev points was taken from
				# the Chemkin theory manual
				pressures = [-math.cos((2 * x - 1) * math.pi / (2 * Pnumber)) for x in range(1, Pnumber+1)]
				pressures = [10**(0.5 * ((math.log10(Pmax) - math.log10(Pmin)) * p + math.log10(Pmax) + math.log10(Pmin))) for p in pressures]
			else:
				slope = (math.log10(Pmax) - math.log10(Pmin)) / (Pnumber - 1)
				pressures = [10**(slope * x + math.log10(Pmin)) for x in range(Pnumber)]

			# Read temperatures (overriding <temperatureRange>)
			temperatureList = xml0.getChildQuantity(unirxnNetworks, 'temperatures', required=False, default=None)
			if temperatureList:
				temperatures = [float(T.simplified) for T in temperatureList]
				Tmin = min(temperatures)
				Tmax = max(temperatures)
			# Read pressures (overriding <pressureRange>)
			pressureList = xml0.getChildQuantity(unirxnNetworks, 'pressures', required=False, default=None)
			if pressureList:
				pressures = [float(P.simplified) for P in pressureList]
				Pmin = min(pressures)
				Pmax = max(pressures)

			if not temperatures:
				raise InvalidInputFileException('Must specify a <temperatureRange> or a <temperatures> element as a child of <unimolecularReactionNetworks>.')
			if not pressures:
				raise InvalidInputFileException('Must specify a <pressureRange> or a <pressures> element as a child of <unimolecularReactionNetworks>.')

			settings.unimolecularReactionNetworks = (method, Tmin, Tmax, temperatures, Pmin, Pmax, pressures, grainSize, numberOfGrains, interpolationModel)
		else:
			settings.unimolecularReactionNetworks = None

		# Create an empty reaction model
		reactionModel = model.CoreEdgeReactionModel()

		# Load databases
		databases = readDatabaseList(xml0, rootElement)
		for database in databases:
			if database[1] == 'general':
				logging.verbose('General database: ' + database[2])
				# Load all databases
				thermo.data.loadThermoDatabase(database[2] + os.sep)
				loadKineticsDatabase(database[2] + os.sep)
				loadFrequencyDatabase(database[2])
			elif database[1] == 'seedmechanism':
				logging.verbose('Seed mechanism: ' + database[2])
				reactionModel.loadSeedMechanism(database[2])
			logging.verbose('')
			
		# Process species
		coreSpecies = []; speciesDict = {}
		speciesList = xml0.getChildElement(rootElement, 'speciesList')
		speciesElements = xml0.getChildElements(speciesList, 'species')
		logging.info('Found ' + str(len(speciesElements)) + ' species')
		for element in speciesElements:
			
			# Load species ID
			sid = str(xml0.getAttribute(element, 'id', required=True))

			# Load the species data from the file
			spec = species.Species()
			spec.fromXML(xml0, element)

			# Check that the species isn't already in the core (e.g. from a seed mechanism)
			existingSpecies = None
			for s in reactionModel.core.species:
				if s.isIsomorphic(spec):
					existingSpecies = s
					break

			if existingSpecies is not None:
				# Point to existing species rather than newly created species
				# This means that any information about the species in the
				# input file will be discarded in favor of the existing species
				# data
				spec = existingSpecies
			else:
				# Handle other aspects of RMG species creation
				logging.verbose('Creating new species %s' % str(spec))
				species.processNewSpecies(spec)

			# All species in RMG input file are immediately added to the core
			coreSpecies.append(spec)

			# Add to local species dictionary (for matching with other parts of file)
			speciesDict[sid] = spec

		logging.verbose('')
		
		# Read model flux tolerance
		fluxTolerance = xml0.getChildElement(rootElement, 'fluxTolerance')
		reactionModel.fluxToleranceKeepInEdge = float(xml0.getChildElementText(fluxTolerance, 'keepInEdge'))
		reactionModel.fluxToleranceMoveToCore = float(xml0.getChildElementText(fluxTolerance, 'moveToCore'))
		reactionModel.fluxToleranceInterrupt = float(xml0.getChildElementText(fluxTolerance, 'interruptSimulation'))
		
		logging.debug('Model flux tolerances set to:')
		logging.debug('\tKeep in edge:         %s' % (reactionModel.fluxToleranceKeepInEdge) )
		logging.debug('\tMove to core:         %s' % (reactionModel.fluxToleranceMoveToCore) )
		logging.debug('\tInterrupt simulation: %s' % (reactionModel.fluxToleranceInterrupt) )
		logging.debug('')
		
		# Read maximum model size
		maxModelSize = xml0.getChildElement(rootElement, 'maximumModelSize')
		if maxModelSize is None:
			logging.debug('Maximum model size is not set')
		else:
			reactionModel.maximumEdgeSpecies = int(xml0.getChildElementText(maxModelSize, 'edgeSpecies'))
			logging.debug('Maximum model size set to:')
			logging.debug('\tEdge species:         %s' % (reactionModel.maximumEdgeSpecies) )
		logging.debug('')

		# Read dynamic simulator
		element = xml0.getChildElement(rootElement, 'simulator')
		reactionModel.absoluteTolerance = float(xml0.getAttribute(element, 'atol'))
		reactionModel.relativeTolerance = float(xml0.getAttribute(element, 'rtol'))
		logging.info('Read dynamic simulator')
		logging.debug('Simulator:')
		logging.debug('\tAbsolute tolerance set to %s' % (reactionModel.absoluteTolerance))
		logging.debug('\tRelative tolerance set to %s' % (reactionModel.relativeTolerance))
		logging.debug('')

		# Read termination targets
		termination = xml0.getChildElement(rootElement, 'termination')
		targetElements = xml0.getChildElements(termination, 'target')
		for element in targetElements:

			targetType = xml0.getAttribute(element, 'type')
			if targetType == 'conversion':
				sid = xml0.getAttribute(element, 'speciesID')
				spec = speciesDict[sid]
				conv = float(xml0.getElementText(element))
				if conv < 0.0 or conv > 1.0:
					raise InvalidInputFileException('Invalid value for termination fractional conversion.')
				reactionModel.termination.append(model.TerminationConversion(spec, conv))
			elif targetType == 'time':
				units = str(xml0.getAttribute(element, 'units'))
				time = float(xml0.getElementText(element))
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

		# Get list of available reaction systems
		import system as systemModule
		availableSystems = systemModule.getAvailableReactionSystems()
		
		# Process reaction systems
		reactionSystems = []
		reactionSystemList = xml0.getChildElement(rootElement, 'reactionSystemList')
		systemElements = xml0.getChildElements(reactionSystemList, 'reactionSystem')
		for systemElement in systemElements:
		
			# Determine the class of reaction system
			rsClass = xml0.getAttribute(systemElement, 'class')
			if rsClass not in availableSystems:
				raise InvalidInputFileException('Reaction system class "%s" not available.' % (rsClass))
		
			# Declare the reaction system and populate it with info
			reactionSystem = availableSystems[rsClass]()
			reactionSystem.fromXML(xml0, systemElement, speciesDict)
			reactionSystem.initializeCantera()
			
			# Append to the list of reaction systems
			reactionSystems.append(reactionSystem)

		# Output info about reaction system
		if len(reactionSystems) == 1:
			logging.info('Found ' + str(len(reactionSystems)) + ' reaction system')
		else:
			logging.info('Found ' + str(len(reactionSystems)) + ' reaction systems')
		for index, reactionSystem in enumerate(reactionSystems):
			logging.debug('Reaction system #%i: %s' % (index+1, reactionSystem))
				
		logging.debug('')
			
		# Cleanup the DOM tree when finished
		xml0.cleanup()
		
	except InvalidInputFileException, e:
		logging.exception(str(e))
		raise e
	except InvalidXMLError, e:
		logging.exception(str(e))
		raise InvalidInputFileException(e.msg)
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

	# Create document with root element <rmgoutput>
	document = XML()
	rootElement = document.initialize('rmgoutput')

	# Process core species list
	speciesList = document.createElement('speciesList', rootElement)
	for spec in reactionModel.core.species:

		element = document.createElement('species', speciesList)
		element.setAttribute('id', str(spec.id))
		element.setAttribute('label', spec.label)
		element.setAttribute('reactive', 'yes' if spec.reactive else 'no')
		
		# Output the structure using CML
		cml = document.createElement('cml', element)
		dom0 = xml.dom.minidom.parseString(spec.toCML())
		cml.appendChild(dom0.documentElement)

		# Output the thermo data
		spec.thermoData.toXML(document, element)

	# Process core reactions list
	reactionList = document.createElement('reactionList', rootElement)
	for rxn in reactionModel.core.reactions:

		if isinstance(rxn.family, reaction.ReactionFamily):
			element = document.createElement('reaction', reactionList)
			element.setAttribute('family', rxn.family.label)
		
		for reac in rxn.reactants:
			reactant = document.createElement('reactant', element)
			reactant.setAttribute('id', str(reac.id))
		for prod in rxn.products:
			product = document.createElement('product', element)
			product.setAttribute('id', str(prod.id))
		
		kinetics = document.createElement('kinetics', element)
		if isinstance(rxn, reaction.PDepReaction):
			pass
		else:
			for k in rxn.kinetics:
				k.toXML(document, kinetics, len(rxn.reactants))

	# Write output file
	document.save(fstr)
	
	# Print to log
	logging.info('')
	logging.info('Output written to ' + fstr)
	
