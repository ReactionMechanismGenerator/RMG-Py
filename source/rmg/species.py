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

import chem
import thermo

################################################################################

speciesDictionary = {}

class Species:
	"""
	Represent a chemical species (including all of its resonance forms). Each
	species has a unique integer `id` assigned automatically by RMG and a
	not-necessarily unique string `label`. The *structure* variable contains a
	list of :class:`Structure` objects representing each resonance form. The
	`reactive` flag is :data:`True` if the species can react and :data:`False`
	if it is inert.
	"""

	# A static counter for the number of species created since the RMG job began.
	numSpecies = 0

	def __init__(self, label='', structure=None, reactive=True):
		"""
		Initialize a Species object.
		"""
		Species.numSpecies += 1
		self.id = Species.numSpecies
		speciesDictionary[self.id] = self

		self.label = label
		if structure:
			structure.simplifyAtomTypes()
			self.structure = [structure]
		else:
			self.structure = []
		self.reactive = reactive

		self.thermo = None
		self.lennardJones = None
		self.spectralData = None

		if structure is not None:
			self.__generateData()

	def __generateData(self):
		"""
		Generate supplemental parameters and information about the species:
		resonance isomers, thermodynamic data, etc.
		"""
		self.getResonanceIsomers()
		#self.getThermoData()

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

def makeNewSpecies(structure):

	# Return an existing species if a match is found
	for id, species in speciesDictionary.iteritems():
		if species.isIsomorphic(structure):
			return species

	# Otherwise make a new species
	return Species(structure.getFormula(), structure)

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

