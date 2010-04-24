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
Contains classes describing chemical reactions:

* :class:`Reaction` - A general chemical reaction

* :class:`ReactionRecipe` - A set of actions to take when applying a reaction

* :class:`ReactionFamily` - A database of a general family of related reactions

* :class:`ReactionFamilySet` - A set of reaction families

"""

import math
import log as logging
import os.path

import constants

from kinetics.model import *

import ctml_writer

################################################################################

class Reaction:
	"""
	Represent a generic chemical reaction. The attributes are:

	================= ==========================================================
	Attribute         Description
	================= ==========================================================
	`id`              A unique integer identifier
	`atomLabels`      A dictionary with the keys representing the labels and the
	                  values the atoms of the reactant or product species that
	                  were labeled at the time the reaction was generated
	`bestKinetics`    The best kinetics for the reaction, always a derived class
	                  of :class:`kinetics.Kinetics`
	`family`          The reaction family that this reaction represents, as a
	                  pointer to a :class:`ReactionFamily` object
	`kinetics`        A list of all of the valid sets of kinetics for the reaction
	`multiplier`      A multiplier to use to increase the reaction rate for cases
	                  when the reaction is generated multiple times due to
	                  different parts of the reactants yielding the same behavior
	`products`        A list of the species that are produced by this reaction
	`reactants`       A list of the species that are consumed by this reaction
	`reverse`         A pointer to the reverse reaction, also a :class:`Reaction`
	                  object
	`canteraReaction` A pointer to the corresponding reaction instance in
	                  Cantera
	`thirdBody`       :data:`True` if the reaction kinetics imply a third body,
	                  :data:`False` if not
	================= ==========================================================

	By convention, the forward reaction is taken to be that for which the
	provided kinetics apply; the reverse kinetics are taken from thermodynamic
	reversibility. Lists of reactions in a model or mechanism should therefore
	only store the forward reaction. Note that the reverse reaction exists as a
	:class:`Reaction` object in the `reverse` attribute because this is a
	convenient way to represent the reverse reaction.
	"""

	def __init__(self, id=0, reactants=None, products=None, family=None, kinetics=None, thirdBody=False):
		"""Initialize a reaction object."""
		self.id = id
		self.reactants = reactants or []
		self.products = products or []
		self.family = family
		self.reverse = None
		self.kinetics = kinetics or []
		self.multiplier = 1.0
		self.thirdBody = thirdBody

		self.E0 = None

		# A cache for the best kinetics for this reaction
		self.bestKinetics = None

		# A dictionary of the labeled atoms for the reactants
		self.atomLabels = {}

		# A pointer to the corresponding reaction in Cantera, if one exists
		self.canteraReaction = None

	def __getstate__(self):
		"""
		Used to specify what should be pickled. In this case, we pickle
		everything except the instance of the Cantera reaction, because that
		isn't pickled properly (nor should it be).
		"""
		pickleMe = self.__dict__.copy()
		if 'canteraReaction' in pickleMe: del pickleMe['canteraReaction']
		return pickleMe

	def __str__(self):
		"""
		Return a string representation of the reaction.
		"""
		#string = ' + '.join([str(s) for s in self.reactants])
		string = ''
		for reactant in self.reactants:
			string += str(reactant) + ' + '
		string = string[0:-3] + ' <=> '
		for product in self.products:
			string += str(product) + ' + '
		string = string[0:-3]
		return string

	def toCantera(self, T=1000, P=1.0e5):
		"""Return a Cantera ctml_writer instance"""
		#  Made up. Unimolecular rate constant 1/s
		#reaction(  "A <=> B",  [1.00000E+00, 0, 0])
		rxnstring = ' + '.join([str(sp) for sp in self.reactants])
		rxnstring += ' <=> '
		rxnstring += ' + '.join([str(sp) for sp in self.products])
		k = self.getBestKinetics(T, P)
		A = k.A
		Ea= k.Ea
		n = k.n
		#import pdb; pdb.set_trace()

		options = []

		makeCanteraReaction = True
		try:
			if self.canteraReaction is not None:
				makeCanteraReaction = False
		except AttributeError:
			pass

		if not makeCanteraReaction:
			# If we're updating this reaction, then remove the original version
			ctml_writer._reactions.remove(self.canteraReaction)
			# If the old reaction was a duplicate, then the new one is too
			if 'duplicate' in self.canteraReaction._options:
				options.append('duplicate')
		else:
			# If we're making this reaction for the first time then we need to
			# check for duplicate reactions
			# Get ID of each reactant and product of this reaction
			reactants = [str(r) for r in self.reactants]; reactants.sort()
			products = [str(p) for p in self.products]; products.sort()
			# Remove any IDs that appear in both the reactant and product lists
			# This is because Cantera treats A --> B + C and A + D --> B + C + D
			# as requiring the duplicate tag
			speciesToRemove = []
			for spec in reactants:
				if spec in products: speciesToRemove.append(spec)
			speciesToRemove = list(set(speciesToRemove))
			for spec in speciesToRemove:
				reactants.remove(spec)
				products.remove(spec)
			# Iterate over all existing Cantera reactions
			for rxn in ctml_writer._reactions:
				# Get ID of each reactant and product
				reac = []; prod = []
				for r, v in rxn._r.iteritems():
					for i in range(int(v)): reac.append(r)
				for p, v in rxn._p.iteritems():
					for i in range(int(v)): prod.append(p)
				reac.sort(); prod.sort()
				# Remove any IDs that appear in both the reactant and product lists
				speciesToRemove = []
				for spec in reac:
					if spec in prod: speciesToRemove.append(spec)
				speciesToRemove = list(set(speciesToRemove))
				for spec in speciesToRemove:
					reac.remove(spec)
					prod.remove(spec)
				# Compare with reactants and products of this reaction
				if (reactants == reac and products == prod) or (reactants == prod and products == reac):
					if 'duplicate' not in options or 'duplicate' not in rxn._options:
						logging.debug('Marking reaction %s as duplicate' % (self))
					if 'duplicate' not in options:
						options.append('duplicate')
					if 'duplicate' not in rxn._options:
						rxn._options.append('duplicate')

		self.canteraReaction = ctml_writer.reaction(rxnstring, ctml_writer.Arrhenius(A, n, Ea), options=options)
		return self.canteraReaction

	def fromXML(self, document, rootElement):
		"""
		Convert a <reaction> element from a standard RMG-style XML input file
		into a Reaction object. `document` is an :class:`io.XML` class
		representing the XML DOM tree, and `rootElement` is the <reaction>
		element in that tree.
		"""

		# Read label attribute
		self.id = str(document.getAttribute(rootElement, 'id', required=True))

		# Read <reactant> elements
		self.reactants = []
		reactantElements = document.getChildElements(rootElement, 'reactant', required=True)
		for reactantElement in reactantElements:
			spec = str(document.getAttribute(reactantElement, 'ref'))
			self.reactants.append(spec)

		# Read <product> elements
		self.products = []
		productElements = document.getChildElements(rootElement, 'product', required=True)
		for productElement in productElements:
			spec = str(document.getAttribute(productElement, 'ref'))
			self.products.append(spec)

		# Read <kinetics> element
		self.kinetics = []
		kineticsElement = document.getChildElement(rootElement, 'kinetics', required=False)
		if kineticsElement:
			format = str(document.getAttribute(kineticsElement, 'type', required=True)).lower()
			if format == 'arrhenius':
				self.kinetics = [ArrheniusModel()]
				self.kinetics[0].fromXML(document, kineticsElement)
			else:
				raise io.InvalidInputFileException('Invalid type "%s" for kinetics element; allowed values are "Arrhenius".' % format)

		# Read <groundStateEnergy> element
		groundStateEnergyElement = document.getChildElement(rootElement, 'groundStateEnergy', required=False)
		if groundStateEnergyElement:
			E0 = document.getChildQuantity(rootElement, 'groundStateEnergy', required=False, default=pq.Quantity(0.0))
			self.E0 = float(E0.simplified)
		else:
			self.E0 = None

	def toXML(self, document, rootElement):
		"""
		Create a <reaction> element as a child of `rootElement` in the XML DOM
		tree represented by `document`, an :class:`io.XML` class. The format
		matches the format of the :meth:`Reaction.fromXML()` function.
		"""

		# Create <species> element with id attribute
		reactionElement = document.createElement('reaction', rootElement)
		document.createAttribute('id', reactionElement, self.id)

		# Write reactants
		for reactant in self.reactants:
			reactantElement = document.createElement('reactant', reactionElement)
			document.createAttribute('ref', reactantElement, reactant.id)

		# Write products
		for product in self.products:
			productElement = document.createElement('product', reactionElement)
			document.createAttribute('ref', productElement, product.id)

		# Write kinetics; format is to be added in kinetics class
		for kinetics in self.kinetics:
			kinetics.toXML(document, reactionElement, len(self.reactants))

	def hasTemplate(self, reactants, products):
		"""
		Return :data:`True` if the reaction matches the template of `reactants`
		and `products`, which are both lists of :class:`species.Species`
		objects.
		"""
		return ((all([spec in self.reactants for spec in reactants]) and
			all([spec in self.products for spec in products])) or
			(all([spec in self.products for spec in reactants]) and
			all([spec in self.reactants for spec in products])))

	def isUnimolecular(self):
		"""
		Return :data:`True` if the forward reaction has one reactant and
		:data:`False` otherwise.
		"""
		return len(self.reactants) == 1

	def isBimolecular(self):
		"""
		Return :data:`True` if the forward reaction has two reactants and
		:data:`False` otherwise.
		"""
		return len(self.reactants) == 2

	def equivalent(self, other):
		"""
		Return :data:`True` if the two reactions are equivalent (i.e. they have
		the same reactants and products and are of the same reaction family) and
		:data:`False` otherwise.
		"""

		if len(self.reactants) != len(other.reactants) or \
		  len(self.products) != len(other.products):
			return False
		elif self.family is not other.family:
			return False

		reactantsMatch = False
		if len(self.reactants) == 1:
			indices = [[0]]
		elif len(self.reactants) == 2:
			indices = [[0, 1], [1, 0]]
		elif len(self.reactants) == 3:
			indices = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
		for index in indices:
			if reactantsMatch: break
			match = True
			for i in range(len(index)):
				if not self.reactants[i].isIsomorphic(other.reactants[index[i]]):
					match = False
			if match:
				reactantsMatch = True

		productsMatch = False
		if len(self.products) == 1:
			indices = [[0]]
		elif len(self.products) == 2:
			indices = [[0, 1], [1, 0]]
		elif len(self.products) == 3:
			indices = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
		for index in indices:
			if productsMatch: break
			match = True
			for i in range(len(index)):
				if not self.products[i].isIsomorphic(other.products[index[i]]):
					match = False
			if match:
				productsMatch = True

		return reactantsMatch and productsMatch

	def getEnthalpyOfReaction(self, T):
		"""
		Return the enthalpy of reaction evaluated at temperature `T`.
		"""
		dHrxn = -self.reactants[0].getEnthalpy(T)
		for reactant in self.reactants[1:]:
			dHrxn -= reactant.getEnthalpy(T)
		for product in self.products:
			dHrxn += product.getEnthalpy(T)
		return dHrxn

	def getEntropyOfReaction(self, T):
		"""
		Return the entropy of reaction evaluated at temperature `T`.
		"""
		dSrxn = -self.reactants[0].getEntropy(T)
		for reactant in self.reactants[1:]:
			dSrxn -= reactant.getEntropy(T)
		for product in self.products:
			dSrxn += product.getEntropy(T)
		return dSrxn

	def getFreeEnergyOfReaction(self, T):
		"""
		Return the Gibbs free energy of reaction evaluated at temperature `T`.
		"""
		dGrxn = -self.reactants[0].getFreeEnergy(T)
		for reactant in self.reactants[1:]:
			dGrxn -= reactant.getFreeEnergy(T)
		for product in self.products:
			dGrxn += product.getFreeEnergy(T)
		return dGrxn

	def getEquilibriumConstant(self, T):
		"""
		Return the equilibrium constant K(T) with respect to concentration
		:math:`K_\\mathrm{c}` evaluated at temperature `T` in K. The returned
		equilibrium constant has units related to concentration (mol, m^3).
		"""
		dGrxn = self.getFreeEnergyOfReaction(T)
		K = math.exp(-dGrxn / constants.R / T)
		# Convert from Ka to Kc; conc is the reference concentration
		conc = 1e5 / constants.R / T
		K *= conc ** (len(self.products) - len(self.reactants))
		return K


	def getBestKinetics(self, T, P=1.0e5):
		"""
		Return the best set of ArrheniusKinetics parameters for the forward reaction
		evaluated at the temperature `T`. This function follows the convention
		that the forward reaction is the one for which we are using the kinetic
		expression, and that the reverse rate constant is evaluated using
		thermochemical equilibrium.
		Evans-Polyani ArrheniusEPKinetics are converted to ArrheniusKinetics
		using dHrxn(298K)
		"""

		# Check cache first
		if self.bestKinetics is not None:
			if self.bestKinetics.isTemperatureInRange(T):
				return self.bestKinetics

		# Check that self.kinetics is storing a list and not a single object
		# If the latter, use that as the best kinetics without any other
		# checking
		if self.kinetics.__class__ != list:
			dHrxn = self.getEnthalpyOfReaction(T)
			self.bestKinetics = self.kinetics.getArrhenius(dHrxn)
			return self.bestKinetics

		kinetics = self.kinetics[:]

		# Prune out all kinetic data not valid at desired temperature
		i = 0
		while i < len(kinetics):
			k = kinetics[i]
			if not k.isTemperatureInRange(T): kinetics.remove(k)
			else: i += 1

		# If no kinetic parameters are left to choose from, print a warning
		# The reaction rate for the reactions is set to zero
		# This may not be the best course of action
#		if len(kinetics) == 0:
#			#logging.warning('Warning: No kinetics available for reaction ' + str(self) + ' at ' + str(T) + ' K.')
#			kinetics = ArrheniusModel(A=0.0, Ea=0.0, n=0.0)
#			kinetics.Trange = [0.0, 100000.0]
#			return kinetics

		# If no kinetic parameters are left to choose from, ignore the
		# temperature ranges
		if len(kinetics) == 0:
			kinetics = self.kinetics[:]

		# Choose kinetics based on rank (i.e. lowest non-zero rank)
		bestRank = kinetics[0].rank
		bestKinetics = kinetics[0]
		for k in kinetics[1:]:
			if k.rank < bestRank and k.rank != 0:
				bestRank = k.rank
				bestKinetics = k
		if isinstance(bestKinetics, ArrheniusEPModel):
			# Convert to ArrheniusKinetics
			# Use T = 298 K to calculate enthalpy and free energy of reaction
			T = 298.0
			dHrxn = self.getEnthalpyOfReaction(T)
			bestKinetics = bestKinetics.getArrhenius(dHrxn)

		self.bestKinetics = bestKinetics
		return self.bestKinetics

	def getRateConstant(self, T, P=1.0e5):
		"""
		Return the value of the rate constant k(T) at the temperature `T`. The
		pressure `P` in Pa is not required.
		"""
		kinetics = self.getBestKinetics(T)
		if kinetics is None:
			raise Exception('Unable to determine the rate constant of reaction ' + str(self) + '.')
		return kinetics.getRateConstant(T)

	def getStoichiometricCoefficient(self, spec):
		"""
		Return the stoichiometric coefficient of species `spec` in the reaction.
		The stoichiometric coefficient is increased by one for each time `spec`
		appears as a product and decreased by one for each time `spec` appears
		as a reactant.
		"""
		stoich = 0
		for reactant in self.reactants:
			if reactant is spec: stoich -= 1
		for product in self.products:
			if product is spec: stoich += 1
		return stoich

	def getRate(self, T, P, conc, totalConc=None):
		"""
		Return the net rate of reaction at temperature `T` and pressure `P`. The
		parameter `conc` is a map with species as keys and concentrations as
		values. A reactant not found in the `conc` map is treated as having zero
		concentration.

		If passed a `totalConc`, it won't bother recalculating it.
		"""

		# Calculate total concentration
		if totalConc is None:
			totalConc=sum( conc.values() )

		# Evaluate rate constant
		rateConstant = self.getRateConstant(T, P)
		if self.thirdBody: rateConstant *= totalConc

		# Evaluate equilibrium constant
		equilibriumConstant = self.getEquilibriumConstant(T)

		# Evaluate forward concentration product
		forward = 1.0
		for reactant in self.reactants:
			if reactant in conc:
				forward = forward * conc[reactant]
			else:
				forward = 0.0
				break

		# Evaluate reverse concentration product
		reverse = 1.0
		for product in self.products:
			if product in conc:
				reverse = reverse * conc[product]
			else:
				reverse = 0.0
				break

		# Return rate
		return rateConstant * (forward - reverse / equilibriumConstant)

	def isIsomerization(self):
		"""
		Return :data:`True` if the reaction is an isomerization, i.e. has the
		form :math:`\\mathrm{A} \\rightleftharpoons \\mathrm{B}`.
		Returns :data:`False` otherwise.
		"""
		return len(self.reactants) == 1 and len(self.products) == 1

	def isDissociation(self):
		"""
		Return :data:`True` if the reaction is a dissocition, i.e. has the
		form :math:`\\mathrm{A} \\rightleftharpoons \\mathrm{B} + \\mathrm{C}`.
		Returns :data:`False` otherwise.
		"""
		return len(self.reactants) == 1 and len(self.products) > 1

	def isAssociation(self):
		"""
		Return :data:`True` if the reaction is an association, i.e. has the
		form :math:`\\mathrm{A} + \\mathrm{B} \\rightleftharpoons \\mathrm{C}`.
		Returns :data:`False` otherwise.
		"""
		return len(self.reactants) > 1 and len(self.products) == 1

	def calculateMicrocanonicalRate(self, Elist, T, reacDensStates, prodDensStates=None):
		"""
		Calculate and return the microcanonical rate coefficients k(E) for the
		forward and reverse reactions from the high-pressure limit canonical
		rate coefficient k(T) using the inverse Laplace transform method. The
		reverse microcanonical rate is calculated using the canonical
		equilibrium constant to ensure that the thermodynamics is correct, even
		when there are small errors in the density of states. For dissociation
		reactions the reverse rate coefficient is actually the product of the
		reverse rate and the product equilibrium distribution.
		"""

		dE = Elist[1] - Elist[0]

		import numpy
		kf = numpy.zeros(len(Elist), numpy.float64)
		kb = numpy.zeros(len(Elist), numpy.float64)

		kinetics = self.getBestKinetics(T)

		reacQ = numpy.sum(reacDensStates * numpy.exp(-Elist / constants.R / T))
		if prodDensStates is not None:
			prodQ = numpy.sum(prodDensStates * numpy.exp(-Elist / constants.R / T))

		Keq = self.getEquilibriumConstant(T)

		if self.isIsomerization():

			kf = kineticsInverseLaplaceTransform(kinetics, self.E0, reacDensStates, Elist, T)
			for r in range(len(Elist)):
				if prodDensStates[r] != 0:
					kb[r] = kf[r] / Keq * (reacDensStates[r] / reacQ) / (prodDensStates[r] / prodQ)

		elif self.isDissociation():

			kf = kineticsInverseLaplaceTransform(kinetics, self.E0, reacDensStates, Elist, T)

			for r in range(len(Elist)):
				kb[r] = kf[r] / Keq * reacDensStates[r] * math.exp(-Elist[r] / constants.R / T) / reacQ

		elif self.isAssociation():

			if self.reactant.densStates is None:
				raise Exception('Unable to process association reaction; no density of states available for the reactant isomers.')

			kf = kineticsInverseLaplaceTransform(kinetics, self.E0, reacDensStates, Elist, T)
			for r in range(len(Elist)):
				kf[r] *= reacDensStates[r] * math.exp(-Elist[r] / constants.R / T) / reacQ

			for r in range(len(Elist)):
				bn = prodDensStates[r] * math.exp(-Elist[r] / constants.R / T) / prodQ
				if bn != 0:
					kb[r] = kf[r] / Keq / bn

		return kf, kb

################################################################################

class PDepReaction(Reaction):
	"""
	A reaction with kinetics that depend on both temperature and pressure. Much
	of the functionality is inherited from :class:`Reaction`, with one
	exception: as the kinetics are explicitly functions of pressure, the
	pressure must be specified when calling :meth:`getRateConstant` and
	:meth:`getBestKinetics`.
	"""

	def __init__(self, id=0, reactants=None, products=None, network=None, kinetics=None):
		Reaction.__init__(self, id=id, reactants=reactants, products=products, family=None)
		self.kinetics = kinetics
		self.network = network

	def getRateConstant(self, T, P):
		"""
		Return the value of the rate constant k(T) at the temperature `T` in K
		and pressure `P` in Pa.
		"""
		return self.kinetics.getRateConstant(T, P)

	def getBestKinetics(self, T, P):
		"""
		Return the best set of ArrheniusKinetics parameters for the forward
		reaction evaluated at the temperature `T` and pressure `P`. Currently
		this simply sets the prefactor to the value of :math:`k(T,P)` and
		sets the other Arrhenius parameters to zero.
		"""
		if isinstance(self.kinetics, PDepArrheniusModel):
			return self.kinetics.getArrhenius(P)
		else:
			k = float(self.getRateConstant(T, P))
			return ArrheniusModel(A=k, n=0.0, Ea=0.0)

	def toCantera(self,T=1000,P=1.0e5):
		"""Add this to Cantera ctml_writer"""
		# create a cantera Reaction
		self.canteraReaction = Reaction.toCantera(self,T,P)
		# replace the forward rate coefficient
		rate_function_of_T_P = self.getRateConstant #(T, P)
		self.canteraReaction._kf = ctml_writer.PdepRate(rate_function_of_T_P)
		return self.canteraReaction

################################################################################

def kineticsInverseLaplaceTransform(kinetics, E0, densStates, Elist, T):
	"""
	Apply the inverse Laplace transform to the modified Arrhenius expression
	`kinetics` for a reaction with reactant density of states `densStates` to
	determine the microcanonical rate coefficient over the range of energies
	`Elist`. For :math:`n = 0` the rate is given exactly by the equation

	.. math:: k(E) = A \\frac{\\rho(E - E_\\mathrm{a})}{\\rho(E)}

	For :math:`n > 0` an exact expression also exists, although it is the more
	complicated equation

	.. math:: k(E) = \\frac{A}{R \\Gamma(n)} \\frac{1}{\\rho(E)} \\int_0^E (E - x)^{n-1} \\rho(x) dx

	For :math:`n < 0` we opt to use the very approximate equation

	.. math:: k(E) = (A T^n) \\frac{\\rho(E - E_\\mathrm{a})}{\\rho(E)}

	For this last equation we also need the temperature `T` of interest.
	"""

	import numpy
	import spectral.states as states

	# Extract the Arrhenius parameters (so we only look them up once; also
	# useful in case we modify them (e.g. for negative Ea))
	A = kinetics.A
	n = kinetics.n
	Ea = kinetics.Ea

	# The ILT method cannot handle negative activation energies, so we move
	# this contribution from the exponential to the preexponential
	if Ea < 0.0:
		A *= math.exp(-Ea / constants.R / T)
		Ea = 0.0

	dE = Elist[1] - Elist[0]
	k = numpy.zeros(len(Elist), numpy.float64)

	if n == 0.0:
		# Determine the microcanonical rate directly
		s = int(math.floor(Ea / dE))
		for r in range(len(Elist)):
			if Elist[r] > E0 and densStates[r] != 0:
				k[r] = A * densStates[r - s] / densStates[r]

	elif n > 0.0:
		import scipy.special
		# Evaluate the inverse Laplace transform of the T**n piece, which only
		# exists for n >= 0
		phi = numpy.zeros(len(Elist), numpy.float64)
		for i, E in enumerate(Elist):
			if E == 0.0:
				phi[i] = 0.0
			else:
				phi[i] = E**(n-1) / (constants.R**n * scipy.special.gamma(n))
		# Evaluate the convolution
		states.convolve(phi, densStates, Elist)
		# Apply to determine the microcanonical rate
		s = int(math.floor(Ea / dE))
		for r in range(len(Elist)):
			if Elist[r] > E0 and densStates[r] != 0:
				k[r] = A * phi[r - s] / densStates[r]

	else:
		# Use the cheating method for n < 0
		s = int(math.floor(Ea / dE))
		for r in range(len(Elist)):
			if Elist[r] > E0 and densStates[r] != 0:
				k[r] = A * T**n * densStates[r - s] / densStates[r]

	return k

################################################################################

# The global list of reactions created at any point during RMG execution
# The list is stored in reverse of the order in which the reactions are created;
# when searching the list, it is more likely to match a recently created
# reaction than an older reaction
reactionDict = {'seed':{}}

global reactionCounter
#: Used to label reactions uniquely. Incremented each time a new reaction is made.
reactionCounter = 0

def checkForExistingReaction(rxn):
	"""
	Check to see if an existing reaction has the same reactants, products, and
	family as `rxn`. Returns :data:`True` or :data:`False` and the matched
	reaction (if found).
	"""
	
	# Get the short-list of reactions with the same family, reactant1 and reactant2
	r1 = rxn.reactants[0]
	if len(rxn.reactants)==1: r2 = None
	else: r2 = rxn.reactants[1]
	try:
		my_reactionList = reactionDict[rxn.family][r1][r2]
	except KeyError: # no such short-list: must be new, unless in seed.
		my_reactionList = []
	
	# Now use short-list to check for matches. All should be in same forward direction.
	for rxn0 in my_reactionList:
		if (rxn0.reactants == rxn.reactants and rxn0.products == rxn.products):
			return True, rxn0
		
	# Now check seed reactions.
	# First check seed short-list in forward direction
	try:
		my_reactionList = reactionDict['seed'][r1][r2]
	except KeyError:
		my_reactionList = []
	for rxn0 in my_reactionList:
		if (rxn0.reactants == rxn.reactants and rxn0.products == rxn.products) or \
			(rxn0.reactants == rxn.products and rxn0.products == rxn.reactants):
			return True, rxn0
	# Now get the seed short-list of the reverse reaction
	r1 = rxn.products[0]
	if len(rxn.products)==1: r2 = None
	else: r2 = rxn.products[1]
	try:
		my_reactionList = reactionDict['seed'][r1][r2]
	except KeyError:
		my_reactionList = []
	for rxn0 in my_reactionList:
		if (rxn0.reactants == rxn.reactants and rxn0.products == rxn.products) or \
			(rxn0.reactants == rxn.products and rxn0.products == rxn.reactants):
			return True, rxn0		
	
	return False, None

def makeNewReaction(forward, checkExisting=True):
	"""
	Make a new reaction given a :class:`Reaction` object `forward`. The kinetics
	of the reaction are estimated and the reaction is added to the global list
	of reactions. Returns the reaction in the direction that corresponds to the
	estimated kinetics, along with whether or not the reaction is new to the
	global reaction list.
	
	The forward direction is determined using the "is_reverse" attribute of the 
	reaction's family.  If the reaction family is its own reverse, then it is 
	made such that the forward reaction is exothermic at 298K.
	"""
	
	# switch it around if it's in the reverse direction
	if forward.family.is_reverse:
		reverse = forward
		forward = reverse.reverse
	else:
		reverse = forward.reverse
	# now forward family is not marked as a reverse family. 
	# if backwards isn't either then it's its own reverse, 
	# and we should ensure the forward is exothermic.
	if not reverse.family.is_reverse:
		if forward.getEnthalpyOfReaction(298) > 0:
			reverse = forward
			forward = reverse.reverse

	if checkExisting:
		found, rxn = checkForExistingReaction(forward)
		if found: return rxn, False

	# Note in the log
	logging.verbose('Creating new %s reaction %s' % (forward.family, forward))
	
	def prepareStructures(forward, reverse, speciesList, atomLabels):
	
		speciesList = speciesList[:]

		for spec in speciesList:
			for struct in spec.structure:
				struct.clearLabeledAtoms()

		structures = []
		for labels in atomLabels:
			found = False
			for spec in speciesList:
				for struct in spec.structure:
					atom = labels.values()[0]
					if isinstance(atom, list):
						for a in atom:
							if not found and a in struct.atoms():
								structures.append(struct)
								found = True
					else:
						if not found and atom in struct.atoms():
							structures.append(struct)
							found = True
				if found:
					speciesList.remove(spec)
					break
		if len(speciesList) != 0:
			raise UndeterminableKineticsException(forward)
	
		if len(structures) == 2 and structures[0] == structures[1]:
			structures[1], map = structures[1].copy(returnMap=True)
			for label, atom in atomLabels[1].iteritems():
				atomLabels[1][label] = map[atom]

		# Apply atom labels to structures
		for labels in atomLabels:
			for label, atom in labels.iteritems():
				atom.label = label
		
		return structures
	
	# By convention, we only work with the reaction in the direction for which
	# we have assigned kinetics from the kinetics database; the kinetics of the
	# reverse of that reaction come from thermodynamics
	
	forwardAtomLabels = [labels.copy() for labels in forward.atomLabels]
	reactantStructures = prepareStructures(forward, reverse, forward.reactants, forwardAtomLabels)
	forwardKinetics = forward.family.getKinetics(forward, reactantStructures)

	if forwardKinetics is None:
		raise UndeterminableKineticsException(forward)
		
	forward.kinetics = forwardKinetics
	
	return processNewReaction(forward)

def processNewReaction(rxn):
	"""
	Once a reaction `rxn` has been created (e.g. via :meth:`makeNewReaction`),
	this function handles other aspects	of preparing it for RMG.
	
	It increments the global reactionCounter, assigns this to the id of the 
	forward and reverse reactions, and stores the forward reaction in the 
	global list of reactions, which is now a dictionary of short-lists 
	(to speed up searching through it).
	"""
	
	# Update counter
	global reactionCounter
	reactionCounter += 1
	rxn.id = rxn.reverse.id = reactionCounter
	
	# Add to the global dict/list of existing reactions (a list broken down by family, r1, r2)
	# identify r1 and r2
	r1 = rxn.reactants[0]
	if len(rxn.reactants)==1: r2 = None
	else: r2 = rxn.reactants[1]
	# make dictionary entries if necessary
	if not reactionDict[rxn.family].has_key(r1):
		reactionDict[rxn.family][r1] = dict()
	if not reactionDict[rxn.family][r1].has_key(r2):
		reactionDict[rxn.family][r1][r2] = list()
	# store this reaction at the top of the relevent short-list
	reactionDict[rxn.family][r1][r2].insert(0, rxn)

	# Return newly created reaction
	return rxn, True

def removeFromGlobalList(rxn):
	"""
	Remove a reaction from the global list of reactions (eg. because we 
	just discarded it from the edge)
	"""
	# Get the short-list of reactions with the same family, reactant1 and reactant2
	r1 = rxn.reactants[0]
	if len(rxn.reactants)==1: r2 = None
	else: r2 = rxn.reactants[1]
	try:
		my_reactionList = reactionDict[rxn.family][r1][r2]
		reactionDict[rxn.family][r1][r2].remove(rxn)
	except KeyError, ValueError: 
		raise Exception("Reaction %s wasn't in the global reaction list to be removed"%rxn)

def makeNewPDepReaction(reactants, products, network, kinetics):
	"""
	Make a new pressure-dependent reaction based on a list of `reactants` and a
	list of `products`. The reaction belongs to the specified `network` and
	has pressure-dependent kinetics given by `kinetics`.

	No checking for existing reactions is made here. The returned PDepReaction
	object is not added to the global list of reactions, as that is intended
	to represent only the high-pressure-limit set. The reactionCounter is
	incremented, however, since the returned reaction can and will exist in
	the model edge and/or core.
	"""

	global reactionCounter

	# Sort reactants and products (to make comparisons easier/faster)
	reactants.sort()
	products.sort()

	forward = PDepReaction(id=reactionCounter+1, reactants=reactants, products=products, network=network, kinetics=kinetics)
	reverse = PDepReaction(id=reactionCounter+1, reactants=products, products=reactants, network=network, kinetics=None)
	forward.reverse = reverse
	reverse.reverse = forward

	reactionCounter += 1
	return forward

