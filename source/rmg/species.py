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
		self.label = label
		self.structure = [structure]
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
		self.generateResonanceIsomers()
		#self.getThermoData()

	def generateResonanceIsomers(self):
		"""
		Generate a list of all of the resonance isomers of this species.
		"""

		if len(self.structure) != 1:
			return

		structure = self.structure[0]

		isomers = [structure]

		# Radicals
		if structure.radicalCount() > 0:
			# Iterate over resonance isomers
			index = 0
			while index < len(isomers):
				isomer = isomers[index]
				# Iterate over radicals in structure; use indices because the
				# pointer to isomer (but not its contents) may be changing
				for i in range(0, len(isomer.atoms())):
					atom = isomer.atoms()[i]
					paths = isomer.findAllDelocalizationPaths(atom)
					for path in paths:

						atom1, atom2, atom3, bond12, bond23 = path

						# Skip if invalid
						if atom1.electronState.order < 1 or \
								not bond12.canIncreaseOrder() or \
								not bond23.canDecreaseOrder():
							continue


						# Make a copy of isomer
						oldIsomer = isomer.copy()
						isomers[index] = oldIsomer
						newIsomer = isomer
						isomer = oldIsomer

						# Adjust to (potentially) new resonance isomer
						atom1.decreaseFreeElectron()
						atom3.increaseFreeElectron()
						bond12.increaseOrder()
						bond23.decreaseOrder()

						# Append to isomer list if unique
						found = False
						for isom in isomers:
							ismatch, map12, map21 = isom.isIsomorphic(newIsomer)
							if ismatch: found = True
						if not found:
							isomers.append(newIsomer)

				# Move to next resonance isomer
				index += 1

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

	def heatCapacity(self, T):
		"""
		Return the heat capacity of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.heatCapacity(T)

	def enthalpy(self, T):
		"""
		Return the enthalpy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.enthalpy(T)

	def entropy(self, T):
		"""
		Return the entropy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.entropy(T)

	def freeEnergy(self, T):
		"""
		Return the Gibbs free energy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.freeEnergy(T)

	def isIsomorphic(self, other):
		"""
		Returns :data:`True` if the two species are isomorphic and data:`False`
		otherwise.
		"""
		for struct1 in self.structure:
			for struct2 in other.structure:
				ismatch, map12, map21 = struct1.isIsomorphic(struct2)
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
			print ismatch, map12, map21
		return (len(maps12) > 0), maps12, maps21

	def __str__(self):
		"""
		Return a string representation of the species, in the form 'label(id)'.
		"""
		return self.label + '(' + str(self.id) + ')'

################################################################################

