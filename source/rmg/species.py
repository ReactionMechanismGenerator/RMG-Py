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
Contains classes describing chemical species.
"""

import log as logging
import math
import pybel
import os
import xml.sax.saxutils
import quantities as pq
import numpy
		
import constants
import settings
import structure
import thermo
import data
import spectral.modes
import ctml_writer

################################################################################

class LennardJones:
	"""
	A set of Lennard-Jones parameters.

	==============  ============================================================
	Attribute       Description
	==============  ============================================================
	`sigma`         The Lennard-Jones sigma parameter in m
	`epsilon`       The Lennard-Jones epsilon parameter in J
	==============  ============================================================

	"""

	def __init__(self, sigma=0.0, epsilon=0.0):
		self.sigma = sigma
		self.epsilon = epsilon

	def fromXML(self, document, rootElement):
		"""
		Convert a <lennardJones> element from a standard RMG-style XML input
		file into a LennardJones object. `document` is an :class:`io.XML`
		class representing the XML DOM tree, and `rootElement` is the
		<lennardJones> element in that tree.
		"""

		# Read <sigma> element
		self.sigma = document.getChildQuantity(rootElement, 'sigma', required=True)
		self.sigma = float(self.sigma.simplified)

		# Read <epsilon> element
		self.epsilon = document.getChildQuantity(rootElement, 'epsilon', required=True)
		self.epsilon = float(self.epsilon.simplified)

	def toXML(self, document, rootElement):
		"""
		Add a <lennardJones> element as a child of `rootElement` using
		RMG-style XML. `document` is an :class:`io.XML` class representing the
		XML DOM tree.
		"""

		# Create <lennardJones> element as child of rootElement
		ljElement = document.createElement('lennardJones', rootElement)
		document.createQuantity('sigma', ljElement, self.sigma, 'm')
		document.createQuantity('epsilon', ljElement, self.epsilon, 'J')

################################################################################

class ThermoSnapshot:
	"""
	A set of thermodynamic state variables for a given state under a single
	set of conditions. During dynamic simulations, this class enables the
	simulator to update the thermodynamics of each species one per iteration.
	The attributes are:
	
	============== ========================================================
	Attribute      Description
	============== ========================================================
	`temperature`  the temperature at which the snapshot was taken in K
	`heatCapacity` the heat capacity at the current conditions in J/mol*K
	`enthalpy`     the enthalpy at the current conditions in J/mol
	`entropy`      the entropy at the current conditions in J/mol*K
	`freeEnergy`   the Gibbs free energy at the current conditions in J/mol
	============== ========================================================
	
	"""

	def __init__(self, temperature=0.0, heatCapacity=0.0, enthalpy=0.0, entropy=0.0):
		self.temperature = temperature
		self.heatCapacity = heatCapacity
		self.enthalpy = enthalpy
		self.entropy = entropy
		self.freeEnergy = enthalpy - temperature * entropy

	def isValid(self, temperature):
		"""
		Return :data:`True` if the current snapshot is still valid at the new
		`temperature` (in K), or :data:`False` if it needs to be recalculated.
		The snapshot is considered valid if the	temperatures are equal to within
		five signficiant figures. The primary benefit is a speedup of simulations
		that are isothermal or steady-state in temperature.
		"""
		if self.temperature == 0: return False
		if temperature==self.temperature: return True
		
		ratio =  temperature / self.temperature
		return (ratio > 0.99999 and ratio < 1.00001)
		

	def update(self, temperature, thermoData):
		"""
		Update the thermodynamics snapshot to a new `temperature` (in K) using
		the information in `thermoData`. The latter can be any object that
		provides the necessary getHeatCapacity(), getEnthalpy(), getEntropy(),
		and getFreeEnergy() methods.
		"""
		self.temperature = temperature
		self.heatCapacity = thermoData.getHeatCapacity(temperature)
		self.enthalpy = thermoData.getEnthalpy(temperature)
		self.entropy = thermoData.getEntropy(temperature)
		self.freeEnergy = self.enthalpy - self.temperature * self.entropy

################################################################################

class Species:
	"""
	Represent a chemical species (including all of its resonance forms). The
	attributes are:

	================  ==========================================================
	Attributes        Description
	================  ==========================================================
	`id`              A unique integer identifier
	`label`           A more descriptive (but not necessarily unique) string
	                  label
	`lennardJones`    The Lennard-Jones parameters for the species
	`reactive`        :data:`True` if the species is reactive, :data:`False` if
	                  inert
	`spectralData`    The spectral data (degrees of freedom) for the species
	`E0`              The ground-state energy of the species in J/mol
	`structure`       A list of :class:`structure.Structure` objects
	                  representing the set of resonance isomers
	`thermoData`      The thermodynamic parameters for the species, always a
	                  derived class of `thermo.ThermoData`
	`thermoSnapshot`  The thermodynamic parameters at the current point in the
	                  simulation
	================  ==========================================================

	"""
	def __repr__(self):
		"""How it looks on the console"""
		return "<Species %d '%s'>"%(self.id, self.label)

	def __init__(self, id=0, label='', structure=None, reactive=True, SMILES=None):
		"""
		Initialize a Species object.
		"""
		self.id = id
		
		self.label = label
		if structure is not None:
			structure.simplifyAtomTypes()
			self.structure = [structure]
		else:
			self.structure = []
		self.reactive = reactive

		self.thermoData = None
		self.thermoSnapshot = ThermoSnapshot()
		self.lennardJones = None
		self.spectralData = None
		self.E0 = None
		
		if SMILES is not None:
			self.fromSMILES(SMILES)

	def __cmp__(self, other):
		"""
		A comparison function that can be used to sort lists of Species objects.
		Currently the sorting method is by increasing ID.
		"""
		return cmp(self.id, other.id)

	def __hash__(self):
		"""
		A hash function that allows for use in dictionaries et al. Currently the
		species ID is used.
		"""
		return self.id
		
	def __str__(self):
		"""
		Return a string representation of the species, in the form 'label(id)'.
		"""
		return self.label + '(' + str(self.id) + ')'
		
	def toCantera(self):
		"""Return a Cantera ctml_writer instance"""
		# contrivedly get a list like ['C', '3', 'H', '9', 'Si', '1']
		atoms = self.structure[0].toOBMol().GetSpacedFormula().split()
		# magically turn that lst into a string like 'C:3 H:9 Si:1'
		atoms = ' '.join([i+':'+j for i,j in zip(*[iter(atoms)]*2)])
		return ctml_writer.species(name = str(self),
		    atoms = " %s "%atoms,
		    thermo = self.thermoData.toCantera(),
			# this escaping should really be done by ctml_writer, but it doesn't do it
		    note = xml.sax.saxutils.escape("%s (%s)"%(self.label,self.thermoData.comment))
		       )
		
	def getFormula(self):
		"""
		Return the chemical formula for the species.
		"""
		return self.structure[0].getFormula()

	def fromXML(self, document, rootElement):
		"""
		Convert a <species> element from a standard RMG-style XML input file
		into a Species object. `document` is an :class:`io.XML` class
		representing the XML DOM tree, and `rootElement` is the <structure>
		element in that tree.
		"""

		# Read id attribute
		self.id = str(document.getAttribute(rootElement, 'id', required=True))

		# Read label attribute
		self.label = str(document.getAttribute(rootElement, 'label', required=False, default=self.id))

		# Read reactive attribute
		self.reactive = str(document.getAttribute(rootElement, 'reactive', required=False, default='yes')).lower()
		self.reactive = not (self.reactive == 'no' or self.reactive == 'false' or self.reactive == 'n')

		# Read <structure> element
		self.structure = []
		structureElement = document.getChildElement(rootElement, 'structure', required=False)
		if structureElement:
			self.structure = [structure.Structure()]
			self.structure[0].fromXML(document, structureElement)

		# Read <thermoData> element
		self.thermoData = None
		thermoElement = document.getChildElement(rootElement, 'thermoData', required=False)
		if thermoElement:
			format = str(document.getAttribute(thermoElement, 'format', required=True)).lower()
			if format == 'group additivity':
				self.thermoData = thermo.ThermoGAData()
				self.thermoData.fromXML(document, thermoElement)
			else:
				raise io.InvalidInputFileException('Invalid format "%s" for thermoData element; allowed values are "group additivity".' % format)

		# Read <groundStateEnergy> element
		E0 = document.getChildQuantity(rootElement, 'groundStateEnergy', required=False, default=pq.Quantity(0.0))
		self.E0 = float(E0.simplified)

		# Read <spectralData> element
		self.spectralData = None
		spectralElement = document.getChildElement(rootElement, 'spectralData', required=False)
		if spectralElement:
			self.spectralData = spectral.modes.SpectralData()
			self.spectralData.fromXML(document, spectralElement)

		# Read <lennardJones> element
		self.lennardJones = None
		ljElement = document.getChildElement(rootElement, 'lennardJones', required=False)
		if ljElement:
			self.lennardJones = LennardJones()
			self.lennardJones.fromXML(document, ljElement)

		# Read <expDownParam> element
		expDownElement = document.getChildElement(rootElement, 'expDownParam', required=False)
		if expDownElement:
			self.expDownParam = float(document.getQuantity(expDownElement).simplified)

	def toXML(self, document, rootElement):
		"""
		Create a <species> element as a child of `rootElement` in the XML DOM
		tree represented by `document`, an :class:`io.XML` class. The format
		matches the format of the :meth:`Species.fromXML()` function.
		"""

		# Create <species> element with id, label, and reactive attributes
		speciesElement = document.createElement('species', rootElement)
		document.createAttribute('id', speciesElement, self.id)
		document.createAttribute('label', speciesElement, self.label)
		document.createAttribute('reactive', speciesElement, 'yes' if self.reactive else 'no')

		# Write structure (only the first resonance form if resonance is present)
		self.structure[0].toXML(document, speciesElement)

		# Write thermo data - the format attribute is written in toXML() of the
		# corresponding ThermoData class
		if self.thermoData: self.thermoData.toXML(document, speciesElement)
		
		# Write ground-state energy
		try:
			if self.E0: document.createQuantity('groundStateEnergy', speciesElement, self.E0 / 1000.0, 'kJ/mol')
		except AttributeError:
			pass
		
		# Write spectral data
		if self.spectralData:
			self.spectralData.toXML(document, speciesElement)

		# Write Lennard-Jones parameters
		if self.lennardJones:
			self.lennardJones.toXML(document, speciesElement)
		
		# Write exponential down parameter
		try:
			if self.expDownParam: document.createQuantity('expDownParam', speciesElement, self.expDownParam / 1000.0, 'kJ/mol')
		except AttributeError:
			pass

	def fromAdjacencyList(self, adjstr):
		"""
		Convert an adjacency list string `adjstr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromAdjacencyList(adjstr)

	def fromCML(self, cmlstr):
		"""
		Convert a string of CML `cmlstr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromCML(cmlstr)
	
	def fromInChI(self, inchistr):
		"""
		Convert an InChI string `inchistr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromInChI(inchistr)

	def fromSMILES(self, smilesstr):
		"""
		Convert a SMILES string `smilesstr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromSMILES(smilesstr)

	def fromOBMol(self, obmol):
		"""
		Convert an OpenBabel OBMol object `obmol` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromOBMol(obmol)

	def toCML(self):
		"""
		Convert a Species object to CML.
		"""
		return self.structure[0].toCML()

	def toInChI(self):
		"""
		Convert a Species object to an InChI string.
		"""
		return self.structure[0].toInChI()

	def toOBMol(self):
		"""
		Convert a Species object to an OBMol object.
		"""
		return self.structure[0].toOBMol()

	def toSMILES(self):
		"""
		Convert a Species object to an InChI string.
		"""
		return self.structure[0].toSMILES()

	def toAdjacencyList(self):
		"""
		Convert a Species object to an adjacency list.
		"""
		return str(self) + '\n' + self.structure[0].toAdjacencyList()

	def getResonanceIsomers(self):
		"""
		Generate all of the resonance isomers of this species. The isomers are
		stored as a list in the `structure` attribute. If the length of
		`structure` is already greater than one, it is assumed that all of the
		resonance isomers have already been generated.
		"""

		if len(self.structure) != 1:
			return

		structure = self.structure[0]

		isomers = [structure]

		# Radicals
		if structure.getRadicalCount() > 0:
			# Iterate over resonance isomers
			index = 0
			while index < len(isomers):
				isomer = isomers[index]

				newIsomers = isomer.getAdjacentResonanceIsomers()
				
				for newIsomer in newIsomers:

					# Append to isomer list if unique
					found = False
					for isom in isomers:
						if isom.isIsomorphic(newIsomer): found = True
					if not found:
						isomers.append(newIsomer)

				# Move to next resonance isomer
				index += 1

		for isomer in isomers:
			isomer.updateAtomTypes()

		self.structure = isomers

	def getThermoData(self):
		"""
		Get thermodynamic data for the species, generating it if necessary.
		
		The specific type of thermo data returned is unknown, but it will be 
		a subclass of :class:`thermo.ThermoData`.  If `self.thermoData` does 
		not yet exist, it is created using :func:`generateThermoData`.
		"""
		if not self.thermoData:
			self.generateThermoData()
		return self.thermoData
		
	def generateThermoData(self):
		"""
		Generate thermodynamic data for the species using the thermo database.
		
		Generates the thermo data for each structure (resonance isomer), 
		picks that with lowest H298 value, and saves it to `self.thermoData`.
		"""
		
		thermoData = []
		for structure in self.structure:
			structure.updateAtomTypes()
			thermoData.append(getThermoData(structure))
		
		# If multiple resonance isomers are present, use the thermo data of
		# the most stable isomer (i.e. one with lowest enthalpy of formation)
		# as the thermo data of the species
		self.thermoData = thermoData[0]
		lowestH298 = self.thermoData.getEnthalpy(298)
		for tdata in thermoData[1:]:
			thisH298 = tdata.getEnthalpy(298)
			if thisH298 < lowestH298:
				self.thermoData = tdata
				lowestH298 = thisH298

		# Put the most stable structure first in the list of structures
		i = thermoData.index(self.thermoData)
		s = self.structure.pop(i)
		self.structure.insert(0, s)
		
		return self.thermoData

	def generateSpectralData(self):
		"""
		Generate spectral data for species using the frequency database. The
		spectral data comprises the molecular degrees of freedom of the
		molecule.
		"""
		# Calculate the thermo data if necessary
		if not self.thermoData:
			self.generateThermoData()
		# Generate the spectral data
		spectral.fit.generateSpectralData(self.structure[0], self.thermoData)

	def calculateDensityOfStates(self, Elist):
		"""
		Calculate and return the density of states in mol/J of the species at
		the specified list of energies `Elist` in J/mol.
		"""

		import unirxn.states as states

		# Return None if the species consists of only one atom
		if len(self.structure[0].atoms()) < 2 or self.spectralData is None:
			return None

		# Make sure we have the necessary information to calculate the
		# density of states
		if not self.spectralData:
			self.generateSpectralData()
			
		# Initialize density of states
		densStates = numpy.zeros(len(Elist), numpy.float64)

		# Create energies in cm^-1 at which to evaluate the density of states
		conv = constants.h * constants.c * 100.0 * constants.Na # [=] J/mol/cm^-1
		Emin = min(Elist) / conv
		Emax = max(Elist) / conv
		dE = (Elist[1] - Elist[0]) / conv
		Elist0 = numpy.arange(Emin, Emax+dE/2, dE)

		# Prepare inputs for density of states function
		vib = numpy.array([mode.frequency for mode in self.spectralData.modes if isinstance(mode, spectral.modes.HarmonicOscillator)])
		rot = numpy.array([mode.frequencies for mode in self.spectralData.modes if isinstance(mode, spectral.modes.RigidRotor)])
		hind = numpy.array([[mode.frequency, mode.barrier] for mode in self.spectralData.modes if isinstance(mode, spectral.modes.HinderedRotor)])
		if len(hind) == 0: hind = numpy.zeros([0,2],numpy.float64)
		linear = 1 if self.structure[0].isLinear() else 0
		symm = self.spectralData.symmetry

		# Calculate the density of states
		densStates, msg = states.densityofstates(Elist0, vib, rot, hind, symm, linear)
		msg = msg.strip()
		if msg != '':
			raise Exception('Error while calculating the density of states for species %s: %s' % (self, msg))

		# Convert density of states from (cm^-1)^-1 to mol/J
		densStates /= conv

		return densStates

	def getHeatCapacity(self, T):
		"""
		Return the heat capacity of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.heatCapacity
		#if self.thermoData is None: self.getThermoData()
		return self.thermoData.getHeatCapacity(T)

	def getEnthalpy(self, T):
		"""
		Return the enthalpy of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.enthalpy
		#if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEnthalpy(T)

	def getEntropy(self, T):
		"""
		Return the entropy of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.entropy
		#if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEntropy(T)

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.freeEnergy
		#if self.thermoData is None: self.getThermoData()
		return self.thermoData.getFreeEnergy(T)

	def getMolecularWeight(self):
		"""
		Return the molecular weight of the species in kg/mol.
		"""
		return self.structure[0].getMolecularWeight()

	def isIsomorphic(self, other):
		"""
		Returns :data:`True` if the two species are isomorphic and data:`False`
		otherwise.
		"""
		if other.__class__ == Species:
			for struct1 in self.structure:
				for struct2 in other.structure:
					if struct1.isIsomorphic(struct2):
						return True
		elif other.__class__ == structure.Structure:
			for struct1 in self.structure:
				if struct1.isIsomorphic(other):
					return True
		return False
	
	def isSubgraphIsomorphic(self, other):
		"""
		Returns :data:`True` if the species is isomorphic with the other
		functional group and data:`False` otherwise.
		"""
		for struct1 in self.structure:
			if struct1.isSubgraphIsomorphic(other):
				return True
		return False

	def findSubgraphIsomorphisms(self, other):
		"""
		Returns a list of the subgraph matches between the species and the
		functional group `other`.
		"""
		maps12 = []; maps21 = []
		for struct1 in self.structure:
			ismatch, map21, map12 = struct1.findSubgraphIsomorphisms(other)
			maps12.extend(map12)
			maps21.extend(map21)
		return (len(maps12) > 0), maps21, maps12

	def calculateLennardJonesParameters(self):
		"""
		Calculate the Lennard-Jones collision parameters for the species. The
		parameters are not returned, but are instead stored in the 
		`lennardJones` attribute.
		"""
		sigma, epsilon = self.structure[0].calculateLennardJonesParameters()
		self.lennardJones = LennardJones(sigma, epsilon)

################################################################################

# The global list of species created at any point during RMG execution
# The list is stored in reverse of the order in which the species are created;
# when searching the list, it is more likely to match a recently created species
# than an older species
#: Global list of species currently in memory.
speciesList = []

# A cache of recently visited species
speciesCache = []
speciesCacheMaxSize = 4

global speciesCounter 
#: Used to label species uniquely. Incremented each time a new species is made.
speciesCounter = 0 

def makeNewSpecies(structure, label='', reactive=True):
	"""
	Attempt to make a new species based on a chemical `structure`, which is a
	:class:`Structure` object.

	The proposed species is checked against the list of existing species; if the
	species already exists, this function returns
	the existing species. If the species does not exist, a :class:`Species`
	object is created and returned after being appended to the global species
	list.
	"""
	global speciesCounter
#	# Recalculate atom types for proposed structure (hopefully not necessary)
#	structure.simplifyAtomTypes()
#	structure.updateAtomTypes()

	# First check cache and return if species is found
	for i, spec in enumerate(speciesCache):
		if spec.isIsomorphic(structure):
			speciesCache.pop(i)
			speciesCache.insert(0, spec)
			return spec

	# Return an existing species if a match is found
	for spec in speciesList:
		if spec.isIsomorphic(structure):
			speciesCache.insert(0, spec)
			if len(speciesCache) > speciesCacheMaxSize: speciesCache.pop()
			return spec

	# Return None if the species has a forbidden structure
	if thermo.forbiddenStructures is not None:
		for lbl, struct in thermo.forbiddenStructures.iteritems():
			if structure.isSubgraphIsomorphic(struct): return None

	# Otherwise make a new species
	if label == '':
#		label = structure.getFormula()
#		for atom in structure.atoms():
#			if atom.hasFreeElectron(): label += 'J'
		label = structure.toSMILES()
	
	# Note in the log
	spec = Species(speciesCounter+1, label, structure, reactive)
	logging.debug('Creating new species %s' % str(spec))
	return processNewSpecies(spec)

def processNewSpecies(spec):
	"""
	Once a species `spec` has been created (e.g. via :meth:`makeNewSpecies`),
	this function handles other aspects	of preparing it for RMG.
	"""
	global speciesCounter

	speciesList.insert(0, spec)
	speciesCounter += 1
	spec.id = speciesCounter
	
	spec.getResonanceIsomers()
	if thermoDatabase is not None:
		spec.getThermoData()

	# Generate spectral data
	if settings.spectralDataEstimation and spec.thermoData and spec.reactive:
		import spectral
		spec.spectralData = spectral.data.generateSpectralData(spec.structure[0], spec.thermoData)
		
	# Generate Lennard-Jones parameters
	spec.calculateLennardJonesParameters()

	# Draw species
	if settings.drawMolecules:
		mol = pybel.Molecule(spec.toOBMol())
		mol.draw(False, os.path.join(settings.outputDirectory, 'species/' + str(spec) + '.png'))

	# Return the newly created species
	speciesCache.insert(0, spec)
	if len(speciesCache) > speciesCacheMaxSize: speciesCache.pop()
	return spec

################################################################################

class ThermoDatabase(data.Database):
	"""
	Represent an RMG thermodynamics database.
	"""

	def __init__(self):
		"""Call the generic `data.Database.__init__()` method.

		This in turn creates
			* self.dictionary = Dictionary()
			* self.library = Library()
			* self.tree = Tree()
		"""
		data.Database.__init__(self)


	def load(self, dictstr, treestr, libstr):
		"""
		Load a thermodynamics group additivity database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. The tree is optional, and should
		be set to '' if not desired.
		"""

		# Load dictionary, library, and (optionally) tree
		data.Database.load(self, dictstr, treestr, libstr)

		# Convert data in library to ThermoData objects or lists of
		# [link, comment] pairs
		for label, item in self.library.iteritems():

			if item is None:
				pass
			elif not item.__class__ is tuple:
				raise data.InvalidDatabaseException('Kinetics library should be tuple at this point. Instead got %r'%data)
			else:
				index,item = item # break apart tuple, recover the 'index' - the beginning of the line in the library file.
				# Is't it dangerous having a local variable with the same name as a module?
				# what if we want to raise another data.InvalidDatabaseException() ?
				if not ( item.__class__ is str or item.__class__ is unicode) :
					raise data.InvalidDatabaseException('Kinetics library data format is unrecognized.')

				items = item.split()
				try:
					thermoData = []; comment = ''
					# First 12 entries are thermo data
					for i in range(0, 12):
						thermoData.append(float(items[i]))
					# Remaining entries are comment
					for i in range(12, len(items)):
						comment += items[i] + ' '

					thermoGAData = thermo.ThermoGAData()
					thermoGAData.fromDatabase(thermoData, comment)
					thermoGAData.index = index
					self.library[label] = thermoGAData
				except (ValueError, IndexError), e:
					# Split data into link string and comment string; store
					# as list of length 2
					link = items[0]
					comment = item[len(link)+1:].strip()
					self.library[label] = [link, comment]

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (dictstr))

		#self.library.removeLinks()

	def toXML(self):
		"""
		Return an XML representation of the thermo database.
		"""
		import xml.dom.minidom
		dom = xml.dom.minidom.parseString('<database type="thermodynamics"></database>')
		root = dom.documentElement

		data.Database.toXML(self, dom, root)

		return dom.toprettyxml()

	def getThermoData(self, structure, atom):
		"""
		Determine the group additivity thermodynamic data for the atom `atom`
		in the structure `structure`.
		"""

		node = self.descendTree(structure, atom, None)

		if node not in self.library:
			# No data present (e.g. bath gas)
			data = thermo.ThermoGAData()
		else:
			data = self.library[node]

		while data.__class__ != thermo.ThermoGAData and data is not None:
			if data[0].__class__ == str or data[0].__class__ == unicode:
				data = self.library[data[0]]

			# This code prints the hierarchy of the found node; useful for debugging
		#result = ''
		#while node is not None:
		#	result = ' -> ' + node + result
		#	node = self.tree.parent[node]
		#print result[4:]
	    #
		return data

	def contains(self, structure):
		"""
		Search the dictionary for the specified `structure`. If found, the label
		corresponding to the structure is returned. If not found, :data:`None`
		is returned.
		"""
		for label, struct in self.dictionary.iteritems():
			if structure.isIsomorphic(struct):
				return label
		return None

################################################################################

class ThermoDatabaseSet:
	"""
	A set of thermodynamics group additivity databases, consisting of a primary
	database of functional groups and a number of secondary databases to provide
	corrections for 1,5-interactions, gauche interactions, radicals, rings,
	and other functionality. There is also a primary library containing data for
	individual species.
	"""


	def __init__(self):
		self.groupDatabase = ThermoDatabase()
		self.int15Database = ThermoDatabase()
		self.gaucheDatabase = ThermoDatabase()
		self.otherDatabase = ThermoDatabase()
		self.radicalDatabase = ThermoDatabase()
		self.ringDatabase = ThermoDatabase()
		self.primaryDatabase = ThermoDatabase()

	def load(self, datapath):
		"""
		Load a set of thermodynamics group additivity databases from the general
		database specified at `datapath`.
		"""

		logging.debug('\tThermodynamics databases:')

		self.groupDatabase.load(datapath + 'thermo/Group_Dictionary.txt', \
			datapath + 'thermo/Group_Tree.txt', \
			datapath + 'thermo/Group_Library.txt')
		logging.debug('\t\tFunctional groups')

		self.int15Database.load(datapath + 'thermo/15_Dictionary.txt', \
			datapath + 'thermo/15_Tree.txt', \
			datapath + 'thermo/15_Library.txt')
		logging.debug('\t\t1,5 interactions')

		self.gaucheDatabase.load(datapath + 'thermo/Gauche_Dictionary.txt', \
			datapath + 'thermo/Gauche_Tree.txt', \
			datapath + 'thermo/Gauche_Library.txt')
		logging.debug('\t\tGauche interactions')

		self.otherDatabase.load(datapath + 'thermo/Other_Dictionary.txt', \
			datapath + 'thermo/Other_Tree.txt', \
			datapath + 'thermo/Other_Library.txt')
		logging.debug('\t\tOther corrections')

		self.radicalDatabase.load(datapath + 'thermo/Radical_Dictionary.txt', \
			datapath + 'thermo/Radical_Tree.txt', \
			datapath + 'thermo/Radical_Library.txt')
		logging.debug('\t\tRadical corrections')

		self.ringDatabase.load(datapath + 'thermo/Ring_Dictionary.txt', \
			datapath + 'thermo/Ring_Tree.txt', \
			datapath + 'thermo/Ring_Library.txt')
		logging.debug('\t\tRing corrections')

		self.primaryDatabase.load(datapath + 'thermo/Primary_Dictionary.txt', \
			'', \
			datapath + 'thermo/Primary_Library.txt')
		logging.debug('\t\tPrimary thermo database')

	def saveXML(self, datapath):
		"""
		Save the loaded databases in the set to XML files.
		"""

		f = open(datapath + 'thermo/group.xml', 'w')
		f.write(self.groupDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/1,5-interactions.xml', 'w')
		f.write(self.int15Database.toXML())
		f.close()

		f = open(datapath + 'thermo/gauche-interactions.xml', 'w')
		f.write(self.gaucheDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/other.xml', 'w')
		f.write(self.otherDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/radicals.xml', 'w')
		f.write(self.radicalDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/ring.xml', 'w')
		f.write(self.ringDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/primary.xml', 'w')
		f.write(self.primaryDatabase.toXML())
		f.close()

	def getThermoData(self, struct):
		"""
		Determine the group additivity thermodynamic data for the structure
		`structure`.
		"""

		import chem
		import structure

		# First check to see if structure is in primary thermo library
		label = self.primaryDatabase.contains(struct)
		if label is not None:
			return self.primaryDatabase.library[label]

		thermoData = thermo.ThermoGAData()

		if struct.getRadicalCount() > 0:
			# Make a copy of the structure so we don't change the original
			saturatedStruct = struct.copy()
			# Saturate structure by replacing all radicals with bonds to
			# hydrogen atoms
			added = {}
			for atom in saturatedStruct.atoms():
				for i in range(0, atom.getFreeElectronCount()):
					H = chem.Atom('H', '0')
					bond = chem.Bond([atom, H], 'S')
					saturatedStruct.addAtom(H)
					saturatedStruct.addBond(bond)
					atom.decreaseFreeElectron()
					if atom not in added:
						added[atom] = []
					added[atom].append(bond)
			# Update the atom types of the saturated structure (not sure why
			# this is necessary, because saturating with H shouldn't be
			# changing atom types, but it doesn't hurt anything and is not
			# very expensive, so will do it anyway)
			saturatedStruct.simplifyAtomTypes()
			saturatedStruct.updateAtomTypes()
			# Get thermo estimate for saturated form of structure
			thermoData = self.getThermoData(saturatedStruct)
			# For each radical site, get radical correction
			# Only one radical site should be considered at a time; all others
			# should be saturated with hydrogen atoms
			for atom in added:

				# Remove the added hydrogen atoms and bond and restore the radical
				for bond in added[atom]:
					H = bond.atoms[1]
					saturatedStruct.removeBond(bond)
					saturatedStruct.removeAtom(H)
					atom.increaseFreeElectron()

				thermoData += self.radicalDatabase.getThermoData(saturatedStruct, {'*':atom})

				# Re-saturate
				for bond in added[atom]:
					H = bond.atoms[1]
					saturatedStruct.addAtom(H)
					saturatedStruct.addBond(bond)
					atom.decreaseFreeElectron()

			# Subtract the enthalpy of the added hydrogens
			thermoData_H = self.primaryDatabase.library['H']
			for bond in added[atom]:
				thermoData.H298 -= thermoData_H.H298
				#thermoData.S298 -= thermoData_H.S298

			# Correct the entropy for the symmetry number

		else:
			# Generate estimate of thermodynamics
			for atom in struct.atoms():
				# Iterate over heavy (non-hydrogen) atoms
				if atom.isNonHydrogen():
					# Get initial thermo estimate from main group database
					thermoData += self.groupDatabase.getThermoData(struct, {'*':atom})
					# Correct for gauche and 1,5- interactions
					thermoData += self.gaucheDatabase.getThermoData(struct, {'*':atom})
					thermoData += self.int15Database.getThermoData(struct, {'*':atom})
					thermoData += self.otherDatabase.getThermoData(struct, {'*':atom})

			# Do ring corrections separately because we only want to match
			# each ring one time; this doesn't work yet
			rings = struct.getSmallestSetOfSmallestRings()
			for ring in rings:

				# Make a temporary structure containing only the atoms in the ring
				ringStructure = structure.Structure()
				for atom in ring: ringStructure.addAtom(atom)
				for atom1 in ring:
					for atom2 in ring:
						if struct.hasBond(atom1, atom2):
							ringStructure.addBond(struct.getBond(atom1, atom2))

				# Get thermo correction for this ring
				thermoData += self.ringDatabase.getThermoData(ringStructure, {})

		return thermoData

thermoDatabase = None
forbiddenStructures = None

################################################################################

def getThermoData(struct, required_class=thermo.ThermoNASAData): # ThermoGAData
	"""
	Get the thermodynamic data associated with `structure` by looking in the
	loaded thermodynamic database.

	`required_class` is the class of thermo object you want returning; default
	is :class:`ThermoNASAData`
	"""
	import constants
	import math

	GAthermoData = thermoDatabase.getThermoData(struct)

	# Correct entropy for symmetry number
	struct.calculateSymmetryNumber()
	GAthermoData.S298 -= constants.R * math.log(struct.symmetryNumber)

	if required_class==thermo.ThermoGAData:
		return GAthermoData  # return here because Wilhoit conversion not wanted

	# Convert to Wilhoit
	rotors = struct.calculateNumberOfRotors()
	atoms = len(struct.atoms())
	linear = struct.isLinear()
	WilhoitData = thermo.convertGAtoWilhoit(GAthermoData,atoms,rotors,linear)

	if required_class==thermo.ThermoWilhoitData:
		return WilhoitData

	# Convert to NASA
	NASAthermoData = thermo.convertWilhoitToNASA(WilhoitData)

	# compute the error for the entire conversion, printing it as info or warning (if it is sufficiently high)
	rmsErr = NASAthermoData.rmsErr(GAthermoData)
	if(rmsErr > 0.35):
	    logging.warning("Poor overall GA-to-NASA fit: Overall RMS error in heat capacity fit = %.3f*R;"%(rmsErr) + " Struct: "+str(struct)+" GAthermoData: "+ str(GAthermoData))
	else:
	    logging.info("Overall RMS error in heat capacity fit = %.3f*R;"%(rmsErr)+ " Struct: "+str(struct)+" GAthermoData: "+ str(GAthermoData))

	if required_class==thermo.ThermoNASAData:
		return NASAthermoData

	# Still not returned?
	raise Exception("Cannot convert themo data into class %r"%(required_class))

################################################################################

if __name__ == '__main__':

	import os.path
	import main
	main.initializeLog(logging.DEBUG)

	datapath = '../data/RMG_database/'

	logging.debug('General database: ' + os.path.abspath(datapath))
	thermoDatabase = ThermoDatabaseSet()
	thermoDatabase.load(datapath)

	for label, struct in thermoDatabase.groupDatabase.dictionary.iteritems():
		label = label.replace('/', '_').replace('\\', '_')
		f = open('data/%s.txt' % (label), 'w')
		f.write(struct.toDOT('graphname'))
		f.close()

#	adjlist = \
#"""
#1 C 0 {2,S} {4,S} {5,S} {6,S}
#2 C 0 {1,S} {3,S} {7,S} {8,S}
#3 C 0 {2,S} {4,S} {9,S} {10,S}
#4 C 0 {3,S} {1,S} {11,S} {12,S}
#5 H 0 {1,S}
#6 H 0 {1,S}
#7 H 0 {2,S}
#8 H 0 {2,S}
#9 H 0 {3,S}
#10 H 0 {3,S}
#11 H 0 {4,S}
#12 H 0 {4,S}
#"""
#
#	structure = structure.Structure()
#	structure.fromAdjacencyList(adjlist)
#	structure.updateAtomTypes()
#
#	print structure.toInChI()
#
#	thermoData = getThermoData(structure)
#
#	T = 1000.0
#	print 'Heat capacity at ' + str(T) + ': ' + str(thermoData.getHeatCapacity(T))
#	print 'Enthalpy at ' + str(T) + ': ' + str(thermoData.getEnthalpy(T))
#	print 'Entropy at ' + str(T) + ': ' + str(thermoData.getEntropy(T))
#	print 'Free energy at ' + str(T) + ': ' + str(thermoData.getFreeEnergy(T))
#
