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
import pybel

import model
import chem

"""
Contains functions for manipulation of RMG input and output files.
"""

################################################################################

class InvalidInputFileException(Exception):
	"""
	An exception used when parsing an RMG input file to indicate that the input
	file is invalid. The msg parameter is used to specify what about the file
	caused the exception to be raised.
	"""	

	def __init__(self, msg):
		self.msg = msg
	
	def __str__(self):
		return 'Invalid XML for RMG input file: ' + self.msg

################################################################################

def getFirstElement(parent, name):
	elements = getElements(parent, name)
	if len(elements) > 0:
		return elements[0]
	else:
		return None

def getElements(parent, name):
	return parent.getElementsByTagName(name)

def getNodeText(node):
	for child in node.childNodes:
		if child.nodeType == xml.dom.Node.TEXT_NODE:
			return child.data
	return ''

################################################################################

def readInputFile(fstr):
	"""
	Parse an RMG input file.
	"""

	try:
		
		# Parse the RMG input XML file into a DOM tree
		dom = xml.dom.minidom.parse(fstr)

		# Process root element (must be a rmginput element)
		root = dom.documentElement
		if root.tagName != 'rmginput':
			raise InvalidInputFileException('Incorrect root element.')
		
		# Initialize the reaction model
		reactionModel = model.CoreEdgeReactionModel()
		
		# Process species
		speciesDict = {}
		elements = getElements(root, 'species')
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
			mol = None
			cmlMolecule = getFirstElement(element, 'cml:molecule')
			inchi = getFirstElement(element, 'inchi')
			smiles = getFirstElement(element, 'smiles')
			if cmlMolecule is not None:
				molecule = str(cmlMolecule.toxml())
				#molecule = molecule.replace('cml:', '')
				#print molecule
				#mol = pybel.readstring('cml', molecule)
			elif inchi is not None:
				inchi = str(getNodeText(inchi))
				mol = pybel.readstring('inchi', inchi)
			elif smiles is not None:
				smiles = str(getNodeText(smiles))
				mol = pybel.readstring('smiles', smiles)
			else:
				raise InvalidInputFileException('Species missing structure information.')
			
			# Create a new structure from the OBMol
			structure = chem.Structure()
			structure.fromOBMol(mol.OBMol)
			
			# Create a new species and append the species to the core
			species = chem.Species(label, structure, reactive)
			reactionModel.core.species.append(species)
			
			# Add to local species dictionary (for matching with other parts of file)
			speciesDict[sid] = species
		
		# Output info about reaction system
		logging.info('Read ' + str(len(reactionModel.core.species)) + ' species')
		for species in reactionModel.core.species:
			logging.debug('Species ' + str(species) + ':')
			if logging.getLogger().isEnabledFor(logging.DEBUG):
				logging.debug('\t' + species.toInChI())
		logging.debug('')
		
		# Process reaction systems
		reactionSystems = []
		elements = getElements(root, 'reactionSystem')
		for element in elements:
		
			# Create a new reaction system
			reactionSystem = model.ReactionSystem()
			
			# Temperature model
			temperatureModel = getFirstElement(element, 'temperatureModel')
			tempModelType = temperatureModel.getAttribute('type')
			if tempModelType == 'isothermal':
				
				# Read the (constant) temperature from the file
				temperature = getFirstElement(temperatureModel, 'temperature')
				value = float(getNodeText(temperature))
				units = str(temperature.getAttribute('units'))
				T = pq.Quantity(value, units); T.units = 'K'
				
				# Set the reaction system's temperature model to isothermal
				reactionSystem.temperatureModel.setIsothermal(T)
				
			else:
				raise InvalidInputFileException('Invalid temperature model type "' + tempModelType + '".')
			
			# Pressure model
			pressureModel = getFirstElement(element, 'pressureModel')
			pressModelType = pressureModel.getAttribute('type')
			if pressModelType == 'isobaric':
				
				# Read the (constant) pressure from the file
				pressure = getFirstElement(pressureModel, 'pressure')
				value = float(getNodeText(pressure))
				units = str(pressure.getAttribute('units'))
				P = pq.Quantity(value, units); P.units = 'Pa'
				
				# Set the reaction system's pressure model to isobaric
				reactionSystem.pressureModel.setIsobaric(P)
				
			else:
				raise InvalidInputFileException('Invalid pressure model type "' + pressModelType + '".')
			
			# Initialize all initial concentrations to zero
			for species in reactionModel.core.species:
				reactionSystem.initialConcentration[species] = pq.Quantity(0.0, 'mol/m**3')
			
			# List of initial concentrations
			concentrations = getElements(element, 'concentration')
			for concentration in concentrations:
			
				# Read the concentration from the file
				value = float(getNodeText(concentration))
				sid = concentration.getAttribute('speciesID')
				units = str(concentration.getAttribute('units'))
				C = pq.Quantity(value, units); C.units = 'mol/m**3'
			
				reactionSystem.initialConcentration[speciesDict[sid]] = C
			
			# Append to list of reaction systems
			reactionSystems.append(reactionSystem)
		
		# Output info about reaction system
		if len(reactionSystems) == 1:
			logging.info('Read ' + str(len(reactionSystems)) + ' reaction system')
		else:
			logging.info('Read ' + str(len(reactionSystems)) + ' reaction systems')
		for index, reactionSystem in enumerate(reactionSystems):
			logging.debug('Reaction system #' + str(index+1) + ':')
			logging.debug('\t' + str(reactionSystem.temperatureModel))
			logging.debug('\t' + str(reactionSystem.pressureModel))
			for species, conc in reactionSystem.initialConcentration.iteritems():
				if species.reactive:
					logging.debug('\tInitial concentration of ' + str(species) + ': ' + str(conc))
				else:
					logging.debug('\tConstant concentration of ' + str(species) + ': ' + str(conc))
				
				
		logging.debug('')
			
	
	except InvalidInputFileException, e:
		logging.exception(str(e))
	except IOError, e:
		logging.exception('Input file "' + e.filename + '" not found.')
	finally:
		# Unlink the DOM tree when finished
		dom.unlink()
		
	return reactionModel, reactionSystems