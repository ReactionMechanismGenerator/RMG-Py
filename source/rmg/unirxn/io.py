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
This module handles the reading and writing of input and output files involving
standalone called to the rmg.unirxn module.
"""

from rmg.io import *
from rmg.species import Species
from rmg.reaction import Reaction

from network import Network, Isomer

def readInputFile(fstr):
	"""
	Parse an input file found at `fstr`. If successful, this function returns
	a :class:`Network` object containing the unimolecular reaction network and
	a tuple of options.
	"""

	try:

		# Parse the RMG input XML file into a DOM tree
		document = XML(path=fstr)

		# Make sure root element is a <rmginput> element
		rootElement = document.getRootElement()
		if rootElement.tagName != 'rmginput':
			raise InvalidInputFileException('Incorrect root element; should be <rmginput>.')

		# Load databases
		databases = readDatabaseList(document, rootElement)
		for database in databases:
			if database[1] == 'general':
				logging.debug('General database: ' + database[2])
				# Load only frequency database
				loadFrequencyDatabase(database[2] + os.sep)
		logging.debug('')

		# Create Network object
		network = Network()

		# Process species
		speciesDict = {}
		speciesListElement = document.getChildElement(rootElement, 'speciesList')
		speciesElements = document.getChildElements(speciesListElement, 'species')
		for element in speciesElements:
			# Load species ID
			sid = str(document.getAttribute(element, 'id', required=True))
			# Load the species data from the file
			species = Species()
			species.fromXML(document, element)
			# Add to local species dictionary (for matching with other parts of file)
			speciesDict[sid] = species
			logging.debug('\tCreated species "%s"' % species.label)
		logging.info('Found ' + str(len(speciesElements)) + ' species')
		
		# Process reactions
		reactionListElement = document.getChildElement(rootElement, 'reactionList')
		reactionElements = document.getChildElements(reactionListElement, 'reaction')
		for reactionElement in reactionElements:
			reaction = Reaction()
			reaction.fromXML(document, reactionElement)
			network.pathReactions.append(reaction)
			logging.debug('\tCreated reaction "%s"' % reaction.id)
		logging.info('Found %i reactions' % (len(reactionElements)))
		
		# Create isomers
		for reaction in network.pathReactions:

			# Create isomer for the reactant
			isomer = None
			for isom in network.isomers:
				if all([species in isom.species for species in reaction.reactants]):
					isomer = isom
			if not isomer:
				isomer = Isomer(reaction.reactants)
				if isomer.isUnimolecular():
					network.isomers.insert(network.numUniIsomers(), isomer)
				else:
					network.isomers.append(isomer)
				logging.debug('\tCreated isomer "%s"' % (' + '.join(isomer.species)))
			
			# Create isomer for the product
			isomer = None
			for isom in network.isomers:
				if all([species in isom.species for species in reaction.products]):
					isomer = isom
			if not isomer:
				isomer = Isomer(reaction.products)
				if isomer.isUnimolecular():
					network.isomers.insert(network.numUniIsomers(), isomer)
				else:
					network.isomers.append(isomer)
				logging.debug('\tCreated isomer "%s"' % (' + '.join(isomer.species)))
			
		logging.info('Found %i isomers' % (len(network.isomers)))

		# Convert string list to species list
		for isomer in network.isomers:
			isomer.species = [speciesDict[r] for r in isomer.species]
		for reaction in network.pathReactions:
			reaction.reactants = [speciesDict[s] for s in reaction.reactants]
			reaction.products = [speciesDict[s] for s in reaction.products]
		
		# Determine isomer ground-state energies
		for isomer in network.isomers:
			isomer.E0 = sum([species.E0 for species in isomer.species])
		# Determine transition state ground-state energies of the reactions
		for reaction in network.pathReactions:
			E0 = sum([species.E0 for species in reaction.reactants])
			reaction.E0 = E0 + reaction.kinetics[0].Ea
		
		# Read bath gas
		bathGasListElement = document.getChildElement(rootElement, 'bathGasList', required=True)
		bathGasElements = document.getChildElements(bathGasListElement, 'bathGas', required=True)
		for bathGasElement in bathGasElements:
			ref = str(document.getAttribute(bathGasElement, 'ref', required=True))
			network.bathGas = speciesDict[ref]
		if not network.bathGas:
			raise io.InvalidXMLError('No bath gas specified.')
		logging.info('Set species "%s" as bath gas' % (network.bathGas.id))



		# Read option list
		optionListElement = document.getChildElement(rootElement, 'optionList', required=True)

		# Read <temperatures>
		temperatures = document.getChildQuantity(optionListElement, 'temperatures', required=True)
		Tlist = [float(T.simplified) for T in temperatures]

		# Read <pressures>
		pressures = document.getChildQuantity(optionListElement, 'pressures', required=True)
		Plist = [float(P.simplified) for P in pressures]

		logging.info('Read %i temperatures and %i pressures' % (len(Tlist), len(Plist)))

		# Read <grainSize>
		grainSize = document.getChildQuantity(optionListElement, 'grainSize', required=False, default=0.0)
		grainSize = float(grainSize.simplified)

		# Read <numberOfGrains>
		numGrains = int(document.getChildQuantity(optionListElement, 'numberOfGrains', required=False, default=0))

		# Read <method>
		method = str(document.getChildElementText(optionListElement, 'method', required=True))

		# Read <interpolationModel>
		model = str(document.getChildElementText(optionListElement, 'interpolationModel', required=True))

		# Cleanup the DOM tree when finished
		document.cleanup()

		return network, Tlist, Plist, grainSize, numGrains, method, model

	#except InvalidInputFileException, e:
	#	logging.exception(str(e))
	#	raise e
	#except InvalidXMLError, e:
	#	logging.exception(str(e))
	#	raise InvalidInputFileException(e.msg)
	#except IOError, e:
	#	logging.exception('Input file "' + e.filename + '" not found.')
	#	raise e
	except xml.parsers.expat.ExpatError, e:
		logging.exception('Invalid XML file: '+e.message+'\n')
		raise InvalidInputFileException('Invalid XML file: '+e.message)

################################################################################
