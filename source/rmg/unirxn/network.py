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
Contains classes that define an internal representation of a unimolecular
reaction network.
"""

################################################################################

class Isomer:
	"""
	An isomer well or product channel of a unimolecular reaction network. This
	represents a single local minimum on the corresponding potential energy
	surface, and is composed of one or two chemical species. The attributes are:

	==============  ============================================================
	Attribute       Description
	==============  ============================================================
	`species`       A list of the species that make up the isomer; items in the
	                list are :class:`Species` objects
	`E0`            The ground-state energy of the isomer
	--------------  ------------------------------------------------------------
	`collFreq`      For unimolecular isomers, the collision frequency in s^-1
	`densStates`    The density of states for the isomer in mol/J
	`eqDist`        The equilibrium distribution for the isomer
	`Q`             The partition function for the isomer
	==============  ============================================================

	The entries in the lower half of the above table are used for temporary
	storage while conducting a pressure-dependent rate coefficient calculation.
	"""

	def __init__(self, species=None):
		self.species = species or None
		self.E0 = None
		self.clean()

	def clean(self):
		"""
		Free some memory by deleting temporary variables used only during a
		single phenomenological rate coefficient estimation. This should be
		called at the end of that calculation.
		"""
		self.collFreq = None
		self.densStates = None
		self.eqDist = None
		self.Q = 0.0

	def __str__(self):
		return ' + '.join(self.species)

	def isUnimolecular(self):
		"""
		Return :data:`True` if the isomer represents a unimolecular isomer well
		or :data:`False` otherwise.
		"""
		return len(self.species) == 1

	def isMultimolecular(self):
		"""
		Return :data:`True` if the isomer represents a multimolecular reactant
		or product channel or :data:`False` otherwise.
		"""
		return len(self.species) > 1

	def calculateDensityOfStates(self, Elist):
		"""
		Calculate the density of states in mol/J of the isomer at the specified
		list of energies `Elist` in J/mol.
		"""

		# Initialize density of states
		self.densStates = numpy.zeros(len(Elist), numpy.float64)

		# Calculate the density of states for each species, convolving when
		# multiple species are present
		for species in self.species:
			densStates = species.calculateDensityOfStates(Elist)
			states.convolve(self.densStates, densStates, Elist)

		# Shift to appropriate energy grain using E0
		index = -1
		for i, E in enumerate(Elist):
			if self.E0 < E and index < 0:
				index = i
		if index >= 0:
			self.densStates[index:len(densStates)] = densStates[0:len(densStates)-index]

	def calculateEqDist(self, Elist, T):
		"""
		Calculate the equilibrium (Boltzmann) distribution of the isomer using
		the stored density of states which corresponds to the list of energies
		`Elist` in J/mol and at the temperature `T` in K. The partition function
		at that temperature is also determined. The equilibrium distribution is
		not returned, but is instead stored in the `eqDist` attribute.
		"""
		dE = Elist[1] - Elist[0]
		self.eqDist = self.densStates * numpy.exp(-Elist / constants.R / T) * dE
		self.Q = sum(self.eqDist)
		self.eqDist /= self.Q

	def calculateCollisionFrequency(self, T, P, bathGas):
		"""
		Calculate the collision frequency for the isomer (unimolecular only)
		at the given temperature `T` in K, pressure `P` in Pa, and using
		species `bathGas` as the bath gas.
		"""

		if not self.isUnimolecular():
			return

		collisionIntegral = 1.16145 * T**(-0.14874) + 0.52487 * math.exp(-0.77320 * T) + 2.16178 * math.exp(-2.43787 * T) \
			-6.435/10000 * T**0.14874 * math.sin(18.0323 * T**(-0.76830) - 7.27371)

		gasConc = P / constants.kB / T
		mu = 1 / (1/self.species[0].molWt + 1/bathGas.molWt) / 6.022e23
		sigma = 0.5 * (self.species[0].lennardJones.sigma + bathGas.lennardJones.sigma)
		#epsilon = 0.5 * (self.species[0].lennardJones.epsilon + bathGas.lennardJones.epsilon)

		# Evaluate collision frequency
		self.collFreq = collisionIntegral * \
			math.sqrt(8 * constants.kB * T / math.pi / mu) * math.pi * sigma**2 * gasConc

################################################################################
