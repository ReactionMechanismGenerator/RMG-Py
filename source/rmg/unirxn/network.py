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

import math
import logging
import numpy

import rmg.constants as constants
import states

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
		densStates = numpy.zeros(len(Elist), numpy.float64)

		# Calculate the density of states for each species, convolving when
		# multiple species are present
		for species in self.species:
			densStates0 = species.calculateDensityOfStates(Elist)
			if densStates0 is not None:
				states.convolve(densStates, densStates0, Elist)

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
		mu = 1 / (1/self.species[0].getMolecularWeight() + 1/bathGas.getMolecularWeight()) / 6.022e23
		sigma = 0.5 * (self.species[0].lennardJones.sigma + bathGas.lennardJones.sigma)
		#epsilon = 0.5 * (self.species[0].lennardJones.epsilon + bathGas.lennardJones.epsilon)

		# Evaluate collision frequency
		self.collFreq = collisionIntegral * \
			math.sqrt(8 * constants.kB * T / math.pi / mu) * math.pi * sigma**2 * gasConc

	def __efficiency(self, T, E0, alpha, densStates, E):

		dE = E[1] - E[0]

		FeNum = 0; FeDen = 0
		Delta1 = 0; Delta2 = 0; DeltaN = 0; Delta = 1

		for r in range(0, len(E)):
			value = densStates[r] * math.exp(-E[r] / constants.R / T)
			if E[r] > E0:
				FeNum += value * dE
				if FeDen == 0:
					FeDen = value * constants.R * T
		if FeDen == 0: return 1.0
		Fe = FeNum / FeDen

		for r in range(0, len(E)):
			value = densStates[r] * math.exp(-E[r] / constants.R / T)
			# Delta
			if E[r] < E0:
				Delta1 += value * dE
				Delta2 += value * dE * math.exp(-(E0 - E[r]) / (Fe * constants.R * T))
			DeltaN += value * dE

		Delta1 /= DeltaN
		Delta2 /= DeltaN

		Delta = Delta1 - (Fe * constants.R * T) / (alpha + Fe * constants.R * T) * Delta2

		beta = (alpha / (alpha + Fe * constants.R * T))**2 / Delta

		if beta < 0 or beta > 1:
			logging.warning('Invalid collision efficiency %s calculated at %s K.' % (beta, T))
			if beta < 0: beta = 0
			elif beta > 1: beta = 1
		
		return beta

	def calculateCollisionEfficiency(self, T, reactions, dEdown, Elist):

		if not self.isUnimolecular():
			return

		# Determine the "barrier height" as the minimum transition state energy connected to that well
		E0 = 1.0e100
		for reaction in reactions:
			if reaction.reactant == self or reaction.product == self:
				if reaction.E0 < E0:
					E0 = reaction.E0

		# Ensure that the barrier height is sufficiently above the ground state
		# Otherwise invalid efficiencies are observed
		if E0 - self.E0 < 100000:
			E0 = self.E0 + 100000

		# Calculate efficiency
		return self.__efficiency(T, E0, dEdown, self.densStates, Elist)

	def getActiveSpaceEnergy(self, reactions):

		Eres = 1.0e100

		for reaction in reactions:
			if reaction.reactant == self or reaction.product == self:
				if reaction.E0 < Eres:
					Eres = reaction.E0

		return Eres

################################################################################

class Network:
	"""
	A representation of a unimolecular reaction network. The attributes are:

	=============== ============================================================
	Attribute       Description
	=============== ============================================================
	`isomers`       A list of :class:`Isomer` objects that make up the network
	`pathReactions` A list of :class:`Reaction` objects that connect adjacent
	                isomers (the "path" reactions)
	`netReactions`  A list of :class:`Reaction` objects that connect any pair of
	                isomers (the "net" reactions)
	`valid`         :data:`True` if the net reaction kinetics are up to date;
	                :data:`False` if the kinetics are out of date and need to
	                be recomputed
	=============== ============================================================

	"""

	def __init__(self):
		self.isomers = []
		self.pathReactions = []
		self.netReactions = []
		self.valid = True

	def numUniIsomers(self):
		"""
		Return the number of unimolecular isomers in the network.
		"""
		return sum([1 for isomer in self.isomers if isomer.isUnimolecular()])

	def numMultiIsomers(self):
		"""
		Return the number of multimolecular isomers in the network.
		"""
		return sum([1 for isomer in self.isomers if isomer.isMultimolecular()])

	def indexOf(self, object):
		"""
		Return the integer index associated with a given `object`, which can be
		either an instance of :class:`Isomer` or :class:`Reaction`.
		"""
		if isinstance(object, Isomer):
			for (index, isomer) in enumerate(self.isomers):
				if isomer is object:
					return index
		elif isinstance(object, Reaction):
			for (index, reaction) in enumerate(self.pathReactions):
				if reaction is object:
					return index
		return -1

	def containsSpecies(self, species):
		"""
		Return :data:`True` if the `species` is a *unimolecular* isomer in the
		network, and :data:`False` if not.
		"""
		for rxn in self.pathReactions:
			if rxn.isIsomerization() or rxn.isDissociation():
				if rxn.reactants[0] is species: return True
			if rxn.isIsomerization() or rxn.isAssociation():
				if rxn.products[0] is species: return True
		return False

	def addPathReaction(self, rxn):
		"""
		Add the :class:`Reaction` object `rxn` to the network as a path 
		reaction.
		"""
		# Add the reaction to the list of path reactions
		self.pathReactions.append(rxn)
		self.invalidate()
	
	def merge(self, other):
		"""
		Merge another :class:`Network` object `other` into this object.
		"""
		self.isomers.extend(other.isomers)
		self.pathReactions.extend(other.pathReactions)
		self.invalidate()
	
	def invalidate(self):
		"""
		Mark a network as in need of a new pressure-dependence calculation.
		"""
		self.valid = False

	def calculateDensitiesOfStates(self, Elist):
		"""
		Calculate the density of states for all isomers in the network whose
		species have spectral data. This is required for all unimolecular
		isomers and for all multimolecular isomers that are reactants of a
		path reaction with the forward reaction defined as the association.
		"""

		# Calculate the density of states for each isomer
		Elist0 = self.__getEnergyGrains(min(Elist), max(Elist), Elist[1] - Elist[0], 0)
		for isomer in self.isomers:
			isomer.calculateDensityOfStates(Elist0)

	def shiftToZeroEnergy(self):
		"""
		Adjust the energies of the system (i.e. the `E0` attribute of isomers)
		such that the lowest energy in the system is set to zero.
		"""
		Emin = 10000000.0
		for isomer in self.isomers:
			isomer.E0 = sum([species.E0 for species in isomer.species])
			if isomer.E0 < Emin:
				Emin = isomer.E0
		for isomer in self.isomers:
			isomer.E0 -= Emin
		for reaction in self.pathReactions:
			reaction.E0 -= Emin

	def __getEnergyGrains(self, Emin, Emax, dE0, nGrains0):
		"""
		Return an array of energy grains that have a minimum of `Emin`, a
		maximum of `Emax`, and either a spacing of `dE0` or have number of
		grains `nGrains0`. The first three parameters are in J/mol.
		"""
		if nGrains0 == 0:
			nGrains = int((Emax - Emin) / dE0) + 1
			dE = dE0
		else:
			nGrains = nGrains0
			dE = (Emax - Emin) / (nGrains0 - 1)

		return numpy.arange(Emin, Emax + dE, dE, numpy.float64)

	def determineEnergyGrains(self, grainSize, numGrains, Tmax):
		"""
		Select a suitable list of energies to use for subsequent calculations.
		The procedure is to (1) calculate the equilibrium distribution of the
		highest-energy isomer at the largest temperature of interest (to get the
		broadest distribution), (2) calculate the energy at which the tail of
		the distribution is some fraction of the maximum, and (3) add the
		difference between the ground-state energy of the isomer and the
		highest ground-state energy in the system (either isomer or transition
		state). Parameters are the desired grain size `grainSize` and number of
		grains `numGrains` and the temperature to use for the equilibrium
		calculation `Tmax`, which should be the highest temperature of interest.
		"""

		# For the purposes of finding the maximum energy we will use 4001 grains
		nE = 4001
		dE = 0.0

		# Determine minimum energy and isomer with minimum ground-state energy
		isomer = None
		for i in self.isomers:
			if isomer is None: isomer = i
			elif i.E0 < isomer.E0:
				isomer = i
		Emin = math.floor(isomer.E0)

		# Determine maximum energy and isomer with maximum ground-state energy
		isomer = None
		for i in self.isomers:
			if i.isUnimolecular():
				if isomer is None: isomer = i
				elif i.E0 > isomer.E0: isomer = i
		Emax0 = isomer.E0

		# (Try to) purposely overestimate Emax using arbitrary multiplier
		# This is to (hopefully) avoid multiple density of states calculations
		mult = 50
		done = False
		while not done:

			Emax = math.ceil(Emax0 + mult * constants.R * Tmax)

			Elist = self.__getEnergyGrains(Emin, Emax, dE, nE)
			isomer.calculateDensityOfStates(Elist)
			isomer.calculateEqDist(Elist, Tmax)

			# Find maximum of distribution
			maxIndex = 0
			value = 0.0
			for r, E in enumerate(Elist):
				if isomer.eqDist[r] > value:
					value = isomer.eqDist[r]
					maxIndex = r

			# If tail of distribution is much lower than the maximum, then we've found bounds for Emax
			tol = 1e-8
			if isomer.eqDist[-1] / value < tol:
				r = nE - 1
				while r > 0 and not done:
					if isomer.eqDist[r] / value > tol:
						done = True
					else:
						r -= 1
				Emax = Elist[r] + max([rxn.getBestKinetics(Tmax).Ea for rxn in self.pathReactions]) - Emax0
			else:
				mult += 50

		# Add difference between isomer ground-state energy and highest
		# transition state or isomer energy
		Emax0_iso = max([isomer.E0 for isomer in self.isomers])
		Emax0_rxn = max([reaction.E0 for reaction in self.pathReactions])
		Emax += max([Emax0_iso, Emax0_rxn]) - isomer.E0

		# Round Emax up to nearest integer
		Emax = math.ceil(Emax)

		# Return the chosen energy grains
		return self.__getEnergyGrains(Emin, Emax, grainSize, numGrains)

	def calculateRateCoefficients(self, Tlist, Plist, Elist, method):
		"""
		Calculate the phenomenological rate coefficients for the network.
		"""

		K = numpy.zeros([len(Tlist), len(Plist),\
			len(self.isomers), len(self.isomers)], numpy.float64)

		for t, T in enumerate(Tlist):

			# Calculate equilibrium distributions
			for isomer in self.isomers:
				if isomer.densStates is not None:
					isomer.calculateEqDist(Elist, T)

#			# DEBUG: Plot equilibrium distributions
#			import pylab
#			for isomer in self.isomers:
#				if isomer.densStates is not None:
#					pylab.plot(Elist / 1000.0, isomer.eqDist, '-')
#			pylab.xlabel('Energy (kJ/mol)')
#			pylab.ylabel('Equilibrium distribution')
#			pylab.show()

			# Calculate microcanonical rates k(E)
			# It might seem odd that this is dependent on temperature, and it
			# isn't -- unless the Arrhenius expression has a negative n
			for reaction in self.pathReactions:
				reaction.kf, reaction.kb = reaction.calculateMicrocanonicalRate(Elist,
					T, reaction.reactant.densStates, reaction.product.densStates)
				
#			# DEBUG: Plot microcanonical rates
#			import pylab
#			for reaction in self.pathReactions:
#				if reaction.isIsomerization():
#					pylab.semilogy(Elist / 1000.0, reaction.kf, '-')
#					pylab.semilogy(Elist / 1000.0, reaction.kb, '--')
#			pylab.xlabel('Energy (kJ/mol)')
#			pylab.ylabel('Microcanonical rate')
#			pylab.show()

			for p, P in enumerate(Plist):

				# Calculate collision frequencies
				for isomer in self.isomers:
					if isomer.isUnimolecular():
						isomer.calculateCollisionFrequency(T, P, self.bathGas)

				# Determine phenomenological rate coefficients using approximate
				# method
				K[t,p,:,:] = self.applyApproximateMethod(T, P, Elist, method)

		return K

	def applyApproximateMethod(self, T, P, Elist, method):
		"""
		Apply the approximate method specified in `method` to estimate the
		phenomenological rate coefficients for the network. This function
		expects that all preparations have already been made, as in the
		:meth:`calculateRateCoefficients` method.
		"""

		# Matrix and vector size indicators
		nIsom = self.numUniIsomers()
		nProd = self.numMultiIsomers()
		nGrains = len(Elist)

		# Equilibrium distribution of eash isomer
		eqDist = numpy.zeros([nIsom,nGrains], numpy.float64)
		for i in range(nIsom): eqDist[i,:] = self.isomers[i].eqDist

		# Active-state energy of each isomer
		Eres = numpy.zeros([nIsom+nProd], numpy.float64)
		for i in range(nIsom+nProd):
			Eres[i] = self.isomers[i].getActiveSpaceEnergy(self.pathReactions)

		# Isomerization, dissociation, and association microcanonical rate
		# coefficients, respectively
		Kij = numpy.zeros([nIsom,nIsom,nGrains], numpy.float64)
		Gnj = numpy.zeros([nProd,nIsom,nGrains], numpy.float64)
		Fim = numpy.zeros([nIsom,nProd,nGrains], numpy.float64)
		for reaction in self.pathReactions:
			i = self.indexOf(reaction.reactant)
			j = self.indexOf(reaction.product)
			if reaction.isIsomerization():
				Kij[j,i,:] = reaction.kf
				Kij[i,j,:] = reaction.kb
			elif reaction.isDissociation():
				Gnj[j-nIsom,i,:] = reaction.kf
				Fim[i,j-nIsom,:] = reaction.kb
			elif reaction.isAssociation():
				Fim[j,i-nIsom,:] = reaction.kf
				Gnj[i-nIsom,j,:] = reaction.kb

		if method.lower() == 'modifiedstrongcollision':

			# Modified collision frequency of each isomer
			collFreq = numpy.zeros([nIsom], numpy.float64)
			for i in range(nIsom): collFreq[i] = self.isomers[i].collFreq * \
				self.isomers[i].calculateCollisionEfficiency(T, self.pathReactions, self.bathGas.expDownParam, Elist)

			# Apply modified strong collision method
			import msc
			K, msg = msc.estimateratecoefficients(T, P, Elist, collFreq, eqDist, Eres,
				Kij, Fim, Gnj)
			msg = msg.strip()
			if msg != '':
				raise Exception('Unable to apply modified strong collision method: %s' % msg)

		elif method.lower() == 'reservoirstate':

			# Average energy transferred in a deactivating collision
			dEdown = self.bathGas.expDownParam

			# Ground-state energy for each isomer
			E0 = numpy.zeros([nIsom], numpy.float64)
			for i in range(nIsom): E0[i] = self.isomers[i].E0

			# The full collision matrix for each isomer
			import mastereqn
			Mcoll = numpy.zeros([nIsom,nGrains,nGrains], numpy.float64)
			for i in range(nIsom):
				collFreq = self.isomers[i].collFreq
				densStates = self.isomers[i].densStates
				Mcoll[i,:,:], msg = mastereqn.collisionmatrix(T, P, Elist, collFreq, densStates, E0[i], dEdown)
				msg = msg.strip()
				if msg != '':
					raise Exception('Unable to determine collision matrix for isomer %i: %s' % (i, msg))

			# Apply reservoir state method
			import rs
			logging.debug('\tCalculating phenomenological rate coefficients %g K, %g bar...' % (T, P / 1e5))
			K, msg = rs.estimateratecoefficients(T, P, Elist, Mcoll, eqDist, E0, Eres,
				Kij, Fim, Gnj, dEdown)
			msg = msg.strip()
			if msg != '':
				raise Exception('Unable to apply reservoir state method: %s' % msg)

		return K

################################################################################
