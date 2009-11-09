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

import logging
import math
import pybel

import constants
import settings
import structure
import data
import thermo
import os

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

	def getFormula(self):
		"""
		Return the chemical formula for the species.
		"""
		return self.structure[0].getFormula()

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
		return self.thermoData

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

	def __str__(self):
		"""
		Return a string representation of the species, in the form 'label(id)'.
		"""
		return self.label + '(' + str(self.id) + ')'

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
	
	speciesCounter += 1
	spec = Species(speciesCounter, label, structure, reactive)
	speciesList.insert(0, spec)
	
	spec.getResonanceIsomers()
	if thermoDatabase is not None:
		spec.getThermoData()
	
	# Draw species
	if settings.drawMolecules:
		mol = pybel.Molecule(spec.toOBMol())
		mol.draw(False, os.path.join(settings.outputDirectory, 'species/' + str(spec) + '.png'))

	# Note in the log
	logging.debug('Created new species ' + str(spec) )# + ': ' + spec.toInChI())
	
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
		import logging

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
