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
Contains classes describing chemical species and their thermodynamics.
"""

import quantities as pq
import logging
import math
import pybel

import constants
import data
import chem
import structure
import thermo

################################################################################

class ThermoSnapshot:
	"""
	A set of thermodynamic state variables for a given state under a single
	set of conditions. During dynamic simulations, this class enables the
	simulator to update the thermodynamics of each species one per iteration.
	The attributes are:
	 ============== ========================================================
	 Attribute      Meaning
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
		
		## This test ensures the temperature is not more than 
		## 100,000 times too low or 100,000 times too high!!!
		## Probably not what was intended  :-)
		#return abs(math.log10(temperature / self.temperature)) < 5
		## I think this is better:
		#return math.log10(abs(1-temperature / self.temperature)) < -5
		## however, this is faster (and this function is called a LOT!)
		ratio =  temperature / self.temperature
		return (ratio > 0.99999 and ratio < 1.00001)
		
# >>> import timeit
# >>> timeit.Timer('math.log10(abs(1-500.0/500.05))<-5', 'import math').timeit()
# 0.78184604644775391
# >>> timeit.Timer('(abs(1-500.0/500.05))<0.00001', 'import math').timeit()
# 0.34678912162780762
# >>> timeit.Timer('ratio=500.0/500.05; (ratio<1.00001 and ratio>0.99999)', 'import math').timeit()
# 0.27086305618286133



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
	Represent a chemical species (including all of its resonance forms). Each
	species has a unique integer `id` assigned automatically by RMG and a
	not-necessarily unique string `label`. The `structure` variable contains a
	list of :class:`Structure` objects representing each resonance form. The
	`reactive` flag is :data:`True` if the species can react and :data:`False`
	if it is inert.
	"""
	def __repr__(self):
		"""How it looks on the console"""
		return "<Species %d '%s'>"%(self.id, self.label)

	def __init__(self, id=0, label='', structure=None, reactive=True):
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
		Generate thermodynamic data for the species by use of the thermo
		database.
		"""
		
		thermoData = []
		for structure in self.structure:
			structure.updateAtomTypes()
			thermoData.append(thermo.getThermoData(structure))
		
		# If multiple resonance isomers are present, use the thermo data of
		# the most stable isomer (i.e. one with lowest enthalpy of formation)
		# as the thermo data of the species
		self.thermoData = thermoData[0]
		for tdata in thermoData[1:]:
			if tdata.H298 < self.thermoData.H298:
				self.thermoData = tdata

	def getHeatCapacity(self, T):
		"""
		Return the heat capacity of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.heatCapacity
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getHeatCapacity(T)

	def getEnthalpy(self, T):
		"""
		Return the enthalpy of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.enthalpy
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEnthalpy(T)

	def getEntropy(self, T):
		"""
		Return the entropy of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.entropy
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEntropy(T)

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy of the species at temperature `T`.
		"""
		if self.thermoSnapshot.isValid(temperature=T):
			return self.thermoSnapshot.freeEnergy
		if self.thermoData is None: self.getThermoData()
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
speciesList = []

# A cache of recently visited species
speciesCache = []
speciesCacheMaxSize = 4

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
	spec = Species(len(speciesList)+1, label, structure, reactive)
	speciesList.insert(0, spec)
	
	spec.getResonanceIsomers()
	if thermo.thermoDatabase is not None:
		spec.getThermoData()
	
	# Draw species in core
	if constants.drawMolecules:
		mol = pybel.Molecule(spec.toOBMol())
		mol.draw(False, constants.outputDir + '/species/' + str(spec) + '.png')

	# Note in the log
	logging.debug('Created new species ' + str(spec) + ': ' + spec.toInChI())
	
	# Return the newly created species
	speciesCache.insert(0, spec)
	if len(speciesCache) > speciesCacheMaxSize: speciesCache.pop()
	return spec

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
