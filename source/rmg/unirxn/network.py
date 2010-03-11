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
import rmg.log as logging
import numpy

import rmg.constants as constants
import states

################################################################################

class UnirxnNetworkException(Exception):
	"""
	An general exception used when an error is encountered while dealing with a
	unimolecular reaction network. The `msg` attribute can be used to store
	a string description of the cause of the exception.
	"""

	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return self.msg

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
		return ' + '.join([str(spec) for spec in self.species])

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
		else:
			raise UnirxnNetworkException('Invalid energy range %s-%s kJ/mol for density of states calculation for %s with ground-state energy %s kJ/mol.'
				% (min(Elist) / 1000.0, max(Elist) / 1000.0, self, self.E0 / 1000.0))

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

		Eres = None

		for reaction in reactions:
			if reaction.reactant == self or reaction.product == self:
				if Eres is None:
					Eres = reaction.E0
				elif reaction.E0 < Eres:
					Eres = reaction.E0
		if Eres is None:
			for rxn in reactions: print rxn
			raise UnirxnNetworkException('Unable to determine active space energy cutoff for isomer %s.' % (self))

		return Eres

################################################################################

class Network:
	"""
	A representation of a unimolecular reaction network. The attributes are:

	=============== ============================================================
	Attribute       Description
	=============== ============================================================
	`id`            A unique identifier for the network
	`isomers`       A list of :class:`Isomer` objects that make up the network
	`pathReactions` A list of :class:`Reaction` objects that connect adjacent
	                isomers (the "path" reactions)
	`netReactions`  A list of :class:`Reaction` objects that connect any pair of
	                isomers (the "net" reactions)
	`valid`         :data:`True` if the net reaction kinetics are up to date;
	                :data:`False` if the kinetics are out of date and need to
	                be recomputed
	`explored`      A list of all of the fully-explored unimolecular isomers
	=============== ============================================================

	"""

	def __init__(self, id=0):
		self.id = id
		self.isomers = []
		self.pathReactions = []
		self.netReactions = []
		self.valid = True
		self.explored = []

	def __repr__(self):
		return '<Network "%s">' % ([str(isomer.species[0]) for isomer in self.isomers if isomer.isUnimolecular()])

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
		either an instance of :class:`Isomer` or :class:`Reaction`. Raises a
		KeyError if the object is not found in the network; thus you can assume
		that the returned index is valid.
		"""
		if isinstance(object, Isomer):
			for (index, isomer) in enumerate(self.isomers):
				if isomer is object:
					return index
		elif isinstance(object, Reaction):
			for (index, reaction) in enumerate(self.pathReactions):
				if reaction is object:
					return index
		raise KeyError('%s not found in network %s.' % (object, self))

	def containsSpecies(self, species):
		"""
		Return :data:`True` if the `species` is a *unimolecular* isomer in the
		network, and :data:`False` if not.
		"""
		for rxn in self.pathReactions:
			if rxn.isIsomerization() or rxn.isDissociation():
				if rxn.reactants[0] == species: return True
			if rxn.isIsomerization() or rxn.isAssociation():
				if rxn.products[0] == species: return True
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
		self.explored.extend(other.explored)
		self.invalidate()
	
	def invalidate(self):
		"""
		Mark a network as in need of a new pressure-dependence calculation.
		"""
		self.valid = False

	def getSpeciesList(self):
		"""
		Return a list of all of the species in the network (excluding the bath
		gas), in no particular order.
		"""
		speciesList = []
		for reaction in self.pathReactions:
			for reactant in reaction.reactants:
				if reactant not in speciesList: speciesList.append(reactant)
			for product in reaction.products:
				if product not in speciesList: speciesList.append(product)
		return speciesList

	def getLeakFlux(self, T, P, conc, totalConc=None):
		"""
		Return the leak flux of the network: the forward flux to all unexplored
		unimolecular isomers in the network.
		"""	
		# Get leak fluxes of all unexplored [unimolecular] isomers
		self.leakFluxes = {}
		for rxn in self.netReactions:
			rate = rxn.getRate(T, P, conc, totalConc)
			if rxn.isIsomerization() or rxn.isDissociation():
				spec = rxn.reactants[0]
				if spec not in self.explored:
					if spec in self.leakFluxes:
						self.leakFluxes[spec] -= rate
					else:
						self.leakFluxes[spec] = -rate
			if rxn.isIsomerization() or rxn.isAssociation():
				spec = rxn.reactants[0]
				if spec not in self.explored:
					if spec in self.leakFluxes:
						self.leakFluxes[spec] += rate
					else:
						self.leakFluxes[spec] = rate
		
		return sum(self.leakFluxes.values())

	def getMaximumLeakSpecies(self):
		"""
		Get the unexplored (unimolecular) isomer with the maximum leak flux.
		"""
		if len(self.leakFluxes) == 0:
			raise UnirxnNetworkException('No unimolecular isomers left to explore!')

		# Choose species with maximum leak flux
		maxSpeciesFlux = 0.0; maxSpecies = None
		for spec, flux in self.leakFluxes.iteritems():
			if maxSpecies is None or flux > maxSpeciesFlux:
				maxSpecies = spec; maxSpeciesFlux = flux

		# Return the species and flux
		return maxSpecies, maxSpeciesFlux

	def calculateDensitiesOfStates(self, Elist):
		"""
		Calculate the density of states for all isomers in the network whose
		species have spectral data. This is required for all unimolecular
		isomers and for all multimolecular isomers that are reactants of a
		path reaction with the forward reaction defined as the association.
		"""
		# Calculate the density of states for each isomer
		for isomer in self.isomers:
			isomer.calculateDensityOfStates(Elist)

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

	def getEnergyGrains(self, Emin, Emax, dE0, nGrains0):
		"""
		Return an array of energy grains that have a minimum of `Emin`, a
		maximum of `Emax`, and either a spacing of `dE0` or have number of
		grains `nGrains0`. The first three parameters are in J/mol.
		"""
		useGrainSize = False; useNumGrains = False
		
		if nGrains0 <= 0 and dE0 != 0.0:
			useGrainSize = False; useNumGrains = True
		elif nGrains0 != 0 and dE0 <= 0.0:
			useGrainSize = True; useNumGrains = False
		else:
			# Choose the tighter constraint
			dE1 = (Emax - Emin) / (nGrains0 - 1)
			useNumGrains = (dE1 < dE0); useGrainSize = not useNumGrains

		if useNumGrains:
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

		# For the purposes of finding the maximum energy we will use 401 grains
		nE = 401
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

			Elist = self.getEnergyGrains(Emin, Emax, dE, nE)
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
		Emax0_iso = max([i.E0 for i in self.isomers])
		Emax0_rxn = max([r.E0 for r in self.pathReactions])
		Emax += max([Emax0_iso, Emax0_rxn])
		
		# Round Emax up to nearest integer
		Emax = math.ceil(Emax)

		# Return the chosen energy grains
		return self.getEnergyGrains(Emin, Emax, grainSize, numGrains)

	def calculateRateCoefficients(self, Tlist, Plist, Elist, method, errorCheck=True, nProd=0):
		"""
		Calculate the phenomenological rate coefficients for the network.
		"""

		K = numpy.zeros([len(Tlist), len(Plist),\
			len(self.isomers), len(self.isomers)], numpy.float64)

		try:

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
					
#				# DEBUG: Plot microcanonical rates
#				import pylab
#				for i, reaction in enumerate(self.pathReactions):
#					if reaction.isIsomerization() or reaction.isDissociation():
#						pylab.semilogy(Elist / 1000.0, reaction.kf, '-')
#					if reaction.isIsomerization() or reaction.isAssociation():
#						pylab.semilogy(Elist / 1000.0, reaction.kb, '--')
#				pylab.xlabel('Energy (kJ/mol)')
#				pylab.ylabel('Microcanonical rate')
#				pylab.show()

				for p, P in enumerate(Plist):

					# Calculate collision frequencies
					for isomer in self.isomers:
						if isomer.isUnimolecular():
							isomer.calculateCollisionFrequency(T, P, self.bathGas)
							
					# Determine phenomenological rate coefficients using approximate
					# method
					K[t,p,:,:] = self.applyApproximateMethod(T, P, Elist, method, errorCheck, nProd)

		except UnirxnNetworkException, e:

			if method.lower() == 'reservoirstate':
				logging.warning(e.msg)
			else:
				logging.error(e.msg)

			# Save network to file for later testing
			fstr = 'unirxn_input.xml'
			logging.info('Troublesome network saved to %s.' % fstr)
			import io
			io.writeInputFile(fstr, self, Tlist, Plist, Elist, method, 'none')

			if method.lower() == 'reservoirstate':
				logging.info('Falling back to modified strong collision for this network.')
				return self.calculateRateCoefficients(Tlist, Plist, Elist, 'modifiedstrongcollision')
			else:
				raise e

		return K

	def applyApproximateMethod(self, T, P, Elist, method, errorCheck=True, nProd=0):
		"""
		Apply the approximate method specified in `method` to estimate the
		phenomenological rate coefficients for the network. This function
		expects that all preparations have already been made, as in the
		:meth:`calculateRateCoefficients` method.
		"""

		logging.debug('Applying %s method at %g K, %g bar...' % (method, T, P*1e-5))

		# Matrix and vector size indicators
		nIsom = self.numUniIsomers()
		nReac = self.numMultiIsomers() - nProd
		nGrains = len(Elist)

		dE = Elist[1] - Elist[0]

		# Density of states per partition function (i.e. normalized density of
		# states with respect to Boltzmann weighting factor) for each isomer
		densStates = numpy.zeros([nIsom,nGrains], numpy.float64)
		for i in range(nIsom): densStates[i,:] = self.isomers[i].densStates * dE / self.isomers[i].Q
		
		# If there are no product channels, we must temporarily create a fake
		# one; this is because f2py can't handle matrices with a dimension of zero
		usingFakeProduct = False
		if nReac == 0 and nProd == 0:
			nReac = 1
			usingFakeProduct = True

		# Active-state energy of each isomer
		Eres = numpy.zeros([nIsom+nReac+nProd], numpy.float64)
		for i, isomer in enumerate(self.isomers):
			Eres[i] = isomer.getActiveSpaceEnergy(self.pathReactions)
		
		# Isomerization, dissociation, and association microcanonical rate
		# coefficients, respectively
		Kij = numpy.zeros([nIsom,nIsom,nGrains], numpy.float64)
		Gnj = numpy.zeros([nReac+nProd,nIsom,nGrains], numpy.float64)
		Fim = numpy.zeros([nIsom,nReac,nGrains], numpy.float64)
		for reaction in self.pathReactions:
			i = self.indexOf(reaction.reactant)
			j = self.indexOf(reaction.product)
			if reaction.isIsomerization():
				Kij[j,i,:] = reaction.kf
				Kij[i,j,:] = reaction.kb
			elif reaction.isDissociation():
				Gnj[j-nIsom,i,:] = reaction.kf
				if j - nIsom < nReac:
					Fim[i,j-nIsom,:] = reaction.kb
			elif reaction.isAssociation():
				if i - nIsom < nReac:
					Fim[j,i-nIsom,:] = reaction.kf
				Gnj[i-nIsom,j,:] = reaction.kb

		if method.lower() == 'modifiedstrongcollision':

			# Modified collision frequency of each isomer
			collFreq = numpy.zeros([nIsom], numpy.float64)
			for i in range(nIsom): collFreq[i] = self.isomers[i].collFreq * \
				self.isomers[i].calculateCollisionEfficiency(T, self.pathReactions, self.bathGas.expDownParam, Elist)

			# Apply modified strong collision method
			import msc
			K, msg = msc.estimateratecoefficients_msc(T, P, Elist, collFreq, densStates, Eres,
				Kij, Fim, Gnj, nisom=nIsom, nreac=nReac, nprod=nProd, ngrains=nGrains)
			msg = msg.strip()
			if msg != '':
				raise UnirxnNetworkException('Unable to apply modified strong collision method: %s' % msg)

		elif method.lower() == 'reservoirstate' or method.lower() == 'chemicaleigenvalues':

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
				densStates0 = self.isomers[i].densStates
				Mcoll[i,:,:], msg = mastereqn.collisionmatrix(T, P, Elist, collFreq, densStates0, E0[i], dEdown)
				msg = msg.strip()
				if msg != '':
					raise UnirxnNetworkException('Unable to determine collision matrix for isomer %i: %s' % (i, msg))

			if method.lower() == 'reservoirstate':

				# Apply reservoir state method
				import rs
				K, msg = rs.estimateratecoefficients_rs(T, P, Elist, Mcoll, densStates, E0, Eres,
					Kij, Fim, Gnj, dEdown, nisom=nIsom, nreac=nReac, nprod=nProd, ngrains=nGrains)
				msg = msg.strip()
			
			elif method.lower() == 'chemicaleigenvalues':

				# Ground-state energy for each isomer
				E0 = numpy.zeros([nIsom+nReac+nProd], numpy.float64)
				for i in range(nIsom+nReac+nProd): E0[i] = self.isomers[i].E0

				# Use free energy to determine equilibrium ratios of each isomer and product channel
				eqRatios = numpy.zeros(nIsom+nReac+nProd, numpy.float64)
				for i, isom in enumerate(self.isomers):
					G = sum([spec.getFreeEnergy(T) for spec in isom.species])
					eqRatios[i] = math.exp(-G / constants.R / T)
				eqRatios /= numpy.sum(eqRatios)
				
				# Apply chemically-significant eigenvalue method
				import cse
				K, msg = cse.estimateratecoefficients_cse(T, P, Elist, Mcoll, E0,
					densStates, eqRatios, Kij, Fim, Gnj, nisom=nIsom, nreac=nReac, nprod=nProd, ngrains=nGrains)
				msg = msg.strip()

				#print K
				#quit()

		if errorCheck:
			if msg == '':
				if not numpy.isfinite(K).all():
					print K
					msg = 'Non-finite rate constant returned at %s K, %s Pa.' % (T, P)

			if msg != '':
				raise UnirxnNetworkException('Unable to apply method %s: %s' % (method, msg))

		# If we had to create a temporary (fake) product channel, then don't
		# return the last row and column of the rate coefficient matrix
		if usingFakeProduct:
			return K[:-1,:-1]
		else:
			return K

################################################################################
