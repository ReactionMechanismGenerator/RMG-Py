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

import chem
import thermo

################################################################################

class Species:
	"""
	Represent a chemical species (including all of its resonance forms). Each
	species has a unique integer `id` assigned automatically by RMG and a
	not-necessarily unique string `label`. The *structure* variable contains a
	list of :class:`Structure` objects representing each resonance form. The
	`reactive` flag is :data:`True` if the species can react and :data:`False`
	if it is inert.
	"""

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

		self.thermo = None
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

	def toInChI(self):
		"""
		Convert a Species object to an InChI string.
		"""
		return self.structure[0].toInChI()

	def toSMILES(self):
		"""
		Convert a Species object to an InChI string.
		"""
		return self.structure[0].toSMILES()

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
						ismatch, map12, map21 = isom.isIsomorphic(newIsomer)
						if ismatch: found = True
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
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getHeatCapacity(T)

	def getEnthalpy(self, T):
		"""
		Return the enthalpy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEnthalpy(T)

	def getEntropy(self, T):
		"""
		Return the entropy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEntropy(T)

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy of the species at temperature `T`.
		"""
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
					ismatch, map12, map21 = struct1.isIsomorphic(struct2)
					if ismatch:
						return True
		elif other.__class__ == chem.Structure:
			for struct1 in self.structure:
				ismatch, map12, map21 = struct1.isIsomorphic(other)
				if ismatch:
					return True
		return False
	
	def isSubgraphIsomorphic(self, other):
		"""
		Returns :data:`True` if the species is isomorphic with the other
		functional group and data:`False` otherwise.
		"""
		for struct1 in self.structure:
			ismatch, map12, map21 = struct1.isSubgraphIsomorphic(other)
			if ismatch:
				return True
		return False

	def findSubgraphIsomorphisms(self, other):
		"""
		Returns a list of the subgraph matches between the species and the
		functional group `other`.
		"""
		maps12 = []; maps21 = []
		for struct1 in self.structure:
			ismatch, map12, map21 = struct1.findSubgraphIsomorphisms(other)
			maps12.extend(map12)
			maps21.extend(map21)
		return (len(maps12) > 0), maps12, maps21

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

	# Return an existing species if a match is found
	for spec in speciesList:
		if spec.isIsomorphic(structure):
			return spec

	# Otherwise make a new species
	if label == '':
		label = structure.getFormula()
		for atom in structure.atoms():
			if atom.hasFreeElectron(): label += 'J'
	spec = Species(len(speciesList)+1, label, structure, reactive)
	speciesList.insert(0, spec)

	spec.getResonanceIsomers()
	spec.getThermoData()

	# Note in the log
	logging.debug('Created new species ' + str(spec) + ': ' + spec.toInChI())

	# Return the newly created species
	return spec

################################################################################

if __name__ == '__main__':

	structure1 = chem.Structure()
	structure1.fromSMILES('C=CC=CC=CC=C[CH]C')
	#structure1.fromSMILES('C=CC=CC=C[CH]C=CC')
	#structure1.fromSMILES('C=CC=C[CH]C=CC=CC')
	#structure1.fromSMILES('C=C[CH]C=CC=CC=CC')
	#structure1.fromSMILES('[CH2]C=CC=CC=CC=CC')
#	isomers = structure1.getAdjacentResonanceIsomers()
#	for isomer in isomers:
#		print isomer.toSMILES()
	
	species1 = Species('', structure1, True)
	
	print len(species1.structure)

