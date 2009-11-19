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
		self.document.writexml(f, indent='', addindent='\t')
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
			raise InvalidXMLError('No "%s" attribute found for  <%s> element.' % (attribute, element.tagName))

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
		value = float(self.getElementText(element))
		return pq.Quantity(value, units)

	def getChildQuantity(self, parentElement, name, required=True):
		"""
		Process an element of the form

			<elementName units="units">value</elementName>

		and return the corresponding :class:`quantities.Quantity` object.
		"""
		element = self.getChildElement(parentElement, name, required)
		return self.getQuantity(element)
		
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

################################################################################

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
		# Process units option
		units = xml0.getChildElementText(optionList, 'units', required=False, default='si')
		pq.set_default_units(units)
		# Read draw molecules option
		drawMolecules = xml0.getChildElementText(optionList, 'drawMolecules', required=False, default='off')
		drawMolecules = drawMolecules.lower()
		settings.drawMolecules = (drawMolecules == 'on' or drawMolecules == 'true' or drawMolecules == 'yes')
		# Read generate plots option
		generatePlots = xml0.getChildElementText(optionList, 'generatePlots', required=False, default='off')
		generatePlots = generatePlots.lower()
		settings.generatePlots = (generatePlots == 'on' or generatePlots == 'true' or generatePlots == 'yes')
		# Read spectral data estimation option
		spectralDataEstimation = xml0.getChildElementText(optionList, 'spectralDataEstimation', required=False, default='off')
		spectralDataEstimation = spectralDataEstimation.lower()
		settings.spectralDataEstimation = (spectralDataEstimation == 'on' or spectralDataEstimation == 'true' or spectralDataEstimation == 'yes')

		# Process databases
		databases = []
		databaseList = xml0.getChildElement(rootElement, 'databaseList')
		databaseElements = xml0.getChildElements(databaseList, 'database')
		for element in databaseElements:
			
			# Get database type
			databaseType = xml0.getAttribute(element, 'type', required=True)
			databaseType = databaseType.lower()
			if databaseType != 'general':
				raise InvalidInputFileException('Invalid database type "' + databaseType + '"; valid types are "general".')
			
			# Get database name and form path
			databaseName = xml0.getElementText(element)
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

		# Load databases
		for database in databases:
			if database[1] == 'general':
				logging.debug('General database: ' + database[2])
				# Load thermo databases
				species.thermoDatabase = species.ThermoDatabaseSet()
				species.thermoDatabase.load(database[2] + os.sep)
				# Load forbidden structures
				thermo.forbiddenStructures = data.Dictionary()
				thermo.forbiddenStructures.load(database[2] + os.sep + 'forbiddenStructure.txt')
				thermo.forbiddenStructures.toStructure()
				# Load kinetic databases (reaction families)
				reaction.kineticsDatabase = reaction.ReactionFamilySet()
				reaction.kineticsDatabase.load(database[2] + os.sep)
				
		logging.debug('')
		
		# Process species
		coreSpecies = []; speciesDict = {}
		speciesList = xml0.getChildElement(rootElement, 'speciesList')
		speciesElements = xml0.getChildElements(speciesList, 'species')
		logging.info('Found ' + str(len(speciesElements)) + ' species')
		for element in speciesElements:
			
			# Attributes of the species element
			sid = xml0.getAttribute(element, 'id')
			label = xml0.getAttribute(element, 'label')
			reactive = xml0.getAttribute(element, 'reactive', required=False, default='yes')
			reactive = reactive.lower()
			reactive = not (reactive == 'no' or reactive == 'false' or reactive == 'n')
			
			# Load structure
			struct = structure.Structure()
			
			cml = xml0.getChildElement(element, 'cml', required=False)
			inchi = xml0.getChildElement(element, 'inchi', required=False)
			smiles = xml0.getChildElement(element, 'smiles', required=False)
			if cml is not None:
				cmlstr = str(xml0.getChildElement(cml, 'molecule', required=True).toxml())
				struct.fromCML(cmlstr)
			elif inchi is not None:
				inchistr = str(xml0.getElementText(inchi))
				struct.fromInChI(inchistr)
			elif smiles is not None:
				smilesstr = str(xml0.getElementText(smiles))
				struct.fromSMILES(smilesstr)
			else:
				raise InvalidInputFileException('Species "%s" missing structure information.' % label)
			
			# Create a new species and append the species to the core
			spec = species.makeNewSpecies(struct, label, reactive)
			coreSpecies.append(spec)
			
			# Add to local species dictionary (for matching with other parts of file)
			speciesDict[sid] = spec

		logging.debug('')
		
		# Create an empty reaction model
		reactionModel = model.CoreEdgeReactionModel()
		
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
	
