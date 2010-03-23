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
import os
import os.path

import constants
import settings
import chem
import data
import structure
import species
from kinetics import *

import ctml_writer

################################################################################

class Reaction:
	"""
	Represent a generic chemical reaction. The attributes are:
	
	================= ==========================================================
	Attribute         Description
	================= ==========================================================
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
	
	def __init__(self, reactants=None, products=None, family=None, kinetics=None, thirdBody=False):
		"""Initialize a reaction object."""
		self.reactants = reactants or []
		self.products = products or []
		self.family = family
		self.reverse = None
		self.kinetics = kinetics or []
		self.multiplier = 1.0
		self.thirdBody = thirdBody

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
				self.kinetics = [ArrheniusKinetics()]
				self.kinetics[0].fromXML(document, kineticsElement)
			else:
				raise io.InvalidInputFileException('Invalid type "%s" for kinetics element; allowed values are "Arrhenius".' % format)

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
			if reactantsMatch: continue
			match = True
			for i in range(len(self.reactants)):
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
			if productsMatch: continue
			match = True
			for i in range(len(self.reactants)):
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

	def getEquilibriumConstant(self, T, conc):
		"""
		Return the equilibrium constant K(T) evaluated at temperature `T` in a
		system with total concentration `conc`.
		"""
		dGrxn = self.getFreeEnergyOfReaction(T)
		K = math.exp(-dGrxn / constants.R / T)
		# Convert from Ka to Kc
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
#			kinetics = ArrheniusKinetics(A=0.0, Ea=0.0, n=0.0)
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
		if isinstance(bestKinetics, ArrheniusEPKinetics):
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
		#totalConc = 0.0
		#for spec in conc:
		#	totalConc += conc[spec]
		if totalConc is None:
			totalConc=sum( conc.values() )

		# Evaluate rate constant
		rateConstant = self.getRateConstant(T, P)
		if self.thirdBody: rateConstant *= totalConc

		# Evaluate equilibrium constant
		equilibriumConstant = self.getEquilibriumConstant(T, totalConc)

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

		Keq = self.getEquilibriumConstant(T, conc=1.0)

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

	def __init__(self, reactants, products, network, kinetics):
		Reaction.__init__(self, reactants, products, None)
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
		if isinstance(self.kinetics, PDepArrheniusKinetics):
			return self.kinetics.getArrhenius(P)
		else:
			k = float(self.getRateConstant(T, P))
			return ArrheniusKinetics(A=k, n=0.0, Ea=0.0)

################################################################################

class InvalidActionException(Exception):
	"""
	An exception to be raised when an invalid action is encountered in a
	reaction recipe.
	"""

	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return self.msg

################################################################################

class ReactionRecipe:
	"""
	Represent a list of actions that, when executed, result in the conversion
	of a set of reactants to a set of products. There are currently five such
	actions:

	============  =======================  =====================================
	Action Name   Arguments                Action
	============  =======================  =====================================
	CHANGE_BOND   center1, order, center2  change the bond order of the bond between center1 and center2 by order; do not break or form bonds
	FORM_BOND     center1, order, center2  form a new bond between center1 and center2 of type order
	BREAK_BOND    center1, order, center2  break the bond between center1 and center2, which should be of type order
	GAIN_RADICAL  center, radical          increase the number of free electrons on center by radical
	LOSE_RADICAL  center, radical          decrease the number of free electrons on center by radical
	============  =======================  =====================================

	The actions are stored as a list in the `actions` attribute. Each action is
	a list of items; the first is the action name, while the rest are the
	action parameters as indicated above.
	"""

	def __init__(self, actions=None):
		self.actions = actions or []

	def addAction(self, action):
		"""
		Add an action to the reaction recipe.
		"""
		self.actions.append(action)

	def getReverse(self):
		"""
		Generate a reaction recipe that, when applied, does the opposite of
		what the current recipe does, i.e., it is the recipe for the reverse
		of the reaction that this is the recipe for.
		"""
		other = ReactionRecipe()
		for action in self.actions:
			if action[0] == 'CHANGE_BOND':
				other.addAction(['CHANGE_BOND', action[1], str(-int(action[2])), action[3]])
			elif action[0] == 'FORM_BOND':
				other.addAction(['BREAK_BOND', action[1], action[2], action[3]])
			elif action[0] == 'BREAK_BOND':
				other.addAction(['FORM_BOND', action[1], action[2], action[3]])
			elif action[0] == 'LOSE_RADICAL':
				other.addAction(['GAIN_RADICAL', action[1], action[2]])
			elif action[0] == 'GAIN_RADICAL':
				other.addAction(['LOSE_RADICAL', action[1], action[2]])
		return other

	def __apply(self, struct, doForward, unique):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. The `doForward` parameter is used to indicate
		whether the forward or reverse recipe should be applied. The atoms in
		the structure should be labeled with the appropriate atom centers.
		"""

		try:

			for action in self.actions:
				if action[0] == 'CHANGE_BOND' or action[0] == 'FORM_BOND' or action[0] == 'BREAK_BOND':

					label1, info, label2 = action[1:]

					# Find associated atoms
					atom1 = struct.getLabeledAtom(label1)
					atom2 = struct.getLabeledAtom(label2)
					if atom1 is None or atom2 is None or atom1 is atom2:
						raise InvalidActionException('Invalid atom labels encountered.')

					# Find associated bond, if present
					bond = struct.getBond(atom1, atom2)

					# If found, apply action
					if action[0] == 'CHANGE_BOND':
						if bond is None:
							raise InvalidActionException('Attempted to change the bond order of a nonexistent bond.')

						info = int(info)
						if abs(info) != 1:
							raise InvalidActionException('Attempted to change the bond order of a bond by more than one at a time.')

						if (info > 0 and doForward) or (info < 0 and not doForward):
							# Increment bond order
							bond.increaseOrder()
							atom1.incrementBond(struct.getBonds(atom1), unique)
							atom2.incrementBond(struct.getBonds(atom2), unique)
						else:
							# Decrement bond order
							bond.decreaseOrder()
							atom1.decrementBond(struct.getBonds(atom1), unique)
							atom2.decrementBond(struct.getBonds(atom2), unique)

					elif (action[0] == 'FORM_BOND' and doForward) or \
						(action[0] == 'BREAK_BOND' and not doForward):
						if bond is not None:
							raise InvalidActionException('Attempted to form a bond that already exists.')
						# Form single bond
						bond = chem.Bond([atom1, atom2], 'S')
						struct.addBond(bond)
						atom1.formBond(struct.getBonds(atom1), unique)
						atom2.formBond(struct.getBonds(atom2), unique)

						if doForward:
							bond = chem.Bond([atom1, atom2], info)
							struct.addBond(bond)
						else:
							bond = struct.getBond(atom1, atom2)
							if bond is None: raise InvalidActionException('Attempted to remove a nonexistent bond.')
							struct.removeBond(bond)

					elif (action[0] == 'BREAK_BOND' and doForward) or \
						(action[0] == 'FORM_BOND' and not doForward):
						if bond is None: raise InvalidActionException('Attempted to remove a nonexistent bond.')
						# Break single bond
						struct.removeBond(bond)
						atom1.breakBond(struct.getBonds(atom1), unique)
						atom2.breakBond(struct.getBonds(atom2), unique)

				elif action[0] == 'LOSE_RADICAL' or action[0] == 'GAIN_RADICAL':

					label, change = action[1:]

					# Find associated atoms
					atom = struct.getLabeledAtom(label)
					if atom is None:
						raise Exception('Invalid atom labels found while attempting to execute reaction recipe.')
					
					change = int(change)
					
					for i in range(change):
						if (action[0] == 'GAIN_RADICAL' and doForward) or \
							(action[0] == 'LOSE_RADICAL' and not doForward):
							# Increment radical count
							atom.increaseFreeElectron()
						elif (action[0] == 'LOSE_RADICAL' and doForward) or \
							(action[0] == 'GAIN_RADICAL' and not doForward):
							# Decrement radical count
							atom.decreaseFreeElectron()

				else:
					raise InvalidActionException('Unknown action "' + action[0] + '" encountered.')

		except InvalidActionException, e:
			logging.warning('Warning: ' + e.msg)
			return False

		return True

	def applyForward(self, struct, unique=True):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures.
		"""
		return self.__apply(struct, True, unique)

	def applyReverse(self, struct, unique=True):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures.
		"""
		return self.__apply(struct, False, unique)

################################################################################

class ReactionFamily(data.Database):
	"""
	Represent a reaction family: a set of reactions with similar chemistry, and
	therefore similar reaction rates. Besides the dictionary, tree, and library
	inherited from :class:`data.Database`, the attributes are:
	
	===========  ===============================================================
	Attribute    Description
	===========  ===============================================================
	`label`      The name of the reaction family
	`template`   A :class:`Reaction` object representing the forward reaction
	             template
	`recipe`     A :class:`ReactionRecipe` object representing the steps to 
	             take when applying the reaction to a set of reactants
	`forbidden`  (Optional) A dictionary of forbidden product structures
	`reverse`    A pointer to the reverse reaction family (or :data:`None` if
	             the family is its own reverse
	===========  ===============================================================

	"""

	def __init__(self, label='', template='', recipe=None):
		data.Database.__init__(self)
		# calling  data.Database.__init__(self) sets:
		# 	self.dictionary = Dictionary()
		# 	self.library = Library()
		# 	self.tree = Tree()
		self.label = label
		self.template = template
		self.recipe = recipe
		self.forbidden = None
		self.reverse = None
		
	def __str__(self):
		return '<ReactionFamily(%s) from %s>'%(self.label,os.path.basename(self._path))
	
	def __str__(self):
		# append the path folder name to the reaction family label, and to the reverse
		return self.label + ' [%s]'%(os.path.basename(self._path))

	def getTemplateLists(self):
		"""
		Return lists containing the top-level nodes of each tree representing
		the reactants and the products. These lists are flattened versions of
		the lists available in self.template.reactants and
		self.template.products.
		"""
		forward = []
		for reactant in self.template.reactants:
			if type(reactant) is list:	forward.extend(reactant)
			else:						forward.append(reactant)
		reverse = []
		for product in self.template.products:
			if type(product) is list:	reverse.extend(product)
			else:						reverse.append(product)
		return forward, reverse

	def load(self, path):
		"""
		Load a reaction family located in the directory `path`. The family
		consists of the files::
		
			dictionary.txt
			tree.txt
			library.txt
			template.txt
			forbiddenGroups.txt
		
		"""
		import re
		# Generate paths to files in the database
		dictstr = os.path.join(path, 'dictionary.txt')
		treestr = os.path.join(path, 'tree.txt')
		libstr  = os.path.join(path, 'library.txt')
		tempstr = os.path.join(path, 'template.txt')
		forbstr = os.path.join(path, 'forbiddenGroups.txt')

		#: The path of the database that was loaded.
		self._path = path

		# Load the dictionary and tree using the generic methods
		# We can't use the generic method to load the library because it has
		# the type ('Arrhenius_EP') as the first meaningful line
		data.Database.load(self, dictstr, treestr, '')
		
		# Load the forbidden groups if the file 'forbiddenGroups.txt' is present
		# This file has the form of a standard dictionary so we can use the
		# standard dictionary loading function
		if os.path.exists(forbstr):
			self.forbidden = data.Dictionary()
			self.forbidden.load(forbstr)
			self.forbidden.toStructure()

		# Load the reaction template information and generate the reverse family
		# This requires that the dictionary and tree be loaded
		self.loadTemplate(tempstr)

		# Process the data in the library
		lines = self.library.load(libstr)
		
		# pop off the first line and check that it's 'Arrhenius_EP'
		line = lines.pop(0)
		if line != 'Arrhenius_EP':
			raise data.InvalidDatabaseException("Was expecting 'Arrhenius_EP' as first line, but got %s, in %s"%(line,libstr))
		
		#figure out how many labels there are
		test_line = lines[0].split()
		for token_no, token in enumerate(test_line):
			# skip token_no=0 because it's the label (but may match the regular expression)
			if token_no and re.match('^[0-9\-.]*$',token):
				# found the Temperature range at token_no
				number_of_groups = token_no-1
				logging.debug("Deduced there are %d groups %s in %s"%(number_of_groups,test_line[1:token_no],libstr))
				break
		else: # didn't break
			raise data.InvalidDatabaseException("Unable to figure out how many groups in %s using line %s"%(libstr,' '.join(test_line)))
		
		self.library.parse(lines, number_of_groups )
		self.processLibraryData()

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (path))

		# Fill in missing nodes in library via an averaging scheme of the
		# existing data; note that this disregards all temperature range
		# information
		forwardTemplate, reverseTemplate = self.getTemplateLists()

		#self.generateMissingEntriesFromBelow(forwardTemplate)
		#self.generateMissingEntriesFromAbove(forwardTemplate)
		
		# give the path to the reverse family too
		if self.reverse:
			self.reverse._path = self._path


	def prune(self, template=None):
		"""
		Remove nodes from the tree and dictionary that are not referred to in
		the library. Nodes that are not referred to in the library but that have
		one or more descendants that are referred to in the library are
		retained. The `template` parameter is a list of the nodes at which to
		begin, e.g. the template lists returned from :func:`getTemplateLists()`.
		"""

		if template is None:
			template, reverseTemplate = self.getTemplateLists()

		# Get lists of all nodes, sorted by top-level node
		nodeLists = []
		for parent in template:
			nodeList = []
			for node in self.tree.children:
				temp = node
				while temp is not None and temp not in template:
					temp = self.tree.parent[temp]
				if temp == parent:
					nodeList.append(node)
			nodeLists.append(nodeList)

		pruneList = []

		for i in range(len(template)):

			nodes = nodeLists[:]
			for node in nodes[i]:

				nodeList = [node]
				children = self.tree.children[node]
				nodeList.extend(children)

				while len(children) > 0:
					temp = children[:]
					children = []
					for child in temp:
						children.extend(self.tree.children[child])
					nodeList.extend(children)

				nodes[i] = nodeList

				nodesList = data.getAllCombinations(nodes)

				# If none of the possible combinations have data, the node
				# is a candidate for pruning
				candidate = True
				for nodeList in nodesList:
					if self.library.getData(nodeList) is not None:
						candidate = False

				# However, we will keep all top-level nodes because they are
				# needed for the template
				parent = self.tree.parent[node]
				if parent is None:
					candidate = False
				# We will also keep a node that is the child of a union in case
				# unions need the children explicitly defined in the tree
				elif self.dictionary[parent] == 'union':
						candidate = False

				if candidate:
					pruneList.append(node)

		# Complete pruning
		for node in pruneList:
			del self.dictionary[node]
			del self.tree.parent[node]
			del self.tree.children[node]

		# Also prune unneeded structures from dictionary
		pruneList = []
		for label, struct in self.dictionary.iteritems():
			if label not in self.tree.children:
				pruneList.append(label)
		for node in pruneList:
			del self.dictionary[node]

	def generateMissingEntriesFromBelow(self, nodes):
		"""
		Generate a nonexisting entry in the library based on an averaging
		scheme.
		"""

		# Get all possible combinations of child nodes
		children = []
		for node in nodes:
			nodeList = [node]; nodeList.extend(self.tree.children[node])
			children.append(nodeList)
		childNodesList = data.getAllCombinations(children)
		# Remove current nodes
		for i in range(childNodesList.count(nodes)):
			childNodesList.remove(nodes)

		# Recursively descend child nodes
		kinetics = []
		for childNodes in childNodesList:
			k = self.generateMissingEntriesFromBelow(childNodes)
			if k is not None: kinetics.append(k)

		# Only generate new entry if data does not exist or rank is zero;
		# otherwise return existing value
		if self.library.getData(nodes) is not None:
			if self.library.getData(nodes).rank > 0:
				return self.library.getData(nodes)

		# Otherwise average all of the kinetics parameters found above
		if len(kinetics) > 0:
			kin = self.averageKinetics(kinetics)
			kin.rank = 5
			self.library.add(nodes, kin)
			print nodes, kin


	def averageKinetics(self, kinetics):
		"""
		Return the average kinetic parameters for the list of kinetic data
		`kinetics`.
		"""
		if len(kinetics) == 0:
			return None

		# Use geometric average of parameters
		lnA = 0.0; E0 = 0.0; n = 0.0; alpha = 0.0
		for k in kinetics:
			lnA += math.log(k.A)
			E0 += k.E0
			n += k.n
			alpha += k.alpha

		lnA /= len(kinetics)
		E0 /= len(kinetics)
		n /= len(kinetics)
		alpha /= len(kinetics)

		kin = ArrheniusEPKinetics(math.exp(lnA), E0, n, alpha)
		kin.Trange = [0.0, 0.0]
		return kin

	def drawFullGraphOfTree(self):
		"""
		Create a PyDOT representation of the current tree.
		"""

		import pydot

		graph = pydot.Dot(size='10,8', page='10,8' ,  rankdir='LR',
				graph_type='digraph', simplify=True, fontsize=10,
				overlap='true', dpi='85',center="True")

		forwardTemplate, reverseTemplate = self.getTemplateLists()

		nodeLists = [[] for top in forwardTemplate]
		for node in self.tree.parent:
			ancestors = self.tree.ancestors(node)
			index = -1
			if len(ancestors) > 0:
				try:
					index = forwardTemplate.index(ancestors[-1])
				except ValueError:
					pass
			elif node in forwardTemplate:
				index = forwardTemplate.index(node)
			if index >= 0 and index < len(forwardTemplate):
				nodeLists[index].append(node)

		nodesList = data.getAllCombinations(nodeLists)

		# Create vertices of graph
		for nodes in nodesList:
			label = ';'.join(nodes)
			node = pydot.Node(label)
			if self.library.getData(nodes) is not None:
				node.set_style('filled')
				node.set_fillcolor('#000000FF')
				node.set_fontcolor('#FFFFFFFF')
			graph.add_node(node)

		# Create edges of graph
		for nodes in nodesList:
			label = ';'.join(nodes)
			for i, node in enumerate(nodes):
				parent = nodes[:]
				parent[i] = self.tree.parent[node]
				if None not in parent:
					parentLabel = ';'.join(parent)
					graph.add_edge(pydot.Edge(parentLabel,label))

		return graph


	def drawGraphOfTree(self, nodes):
		"""draw a graph of the tree"""
		import pydot, re
		g=pydot.Dot(size='10,8', page='10,8' ,  rankdir='LR',
				graph_type='digraph', simplify=True, fontsize=10,
				overlap='true', dpi='85',center="True")

		rates =  self.library.keys() # known reaction rates in library
		# add reactionrate nodes to graph
		for rate in rates:
			g.add_node(pydot.Node(rate))
		# add edges from each node to each of its ancestors
		for rate in rates:
			nodes=rate.split(';')
			for trialRate in rates: # the one we are testing for ancestry
				if trialRate == rate: continue # reject if itself
				for i,node in enumerate(nodes):
					trialNode = trialRate.split(';')[i]
					if trialNode == node: continue # ok if equal
					if trialNode not in self.tree.ancestors(node): break
					# if stopped through break, then it will not run the else clause
				else:
					# loop fell through all nodes without breaking:
					# trialRate must be ancestor of rate
					g.add_edge(pydot.Edge(trialRate,rate))
		g.set('fontsize','10')
		format='svg'
		prog='dot'
		f=open(self.label+'.dot','w')
		f.write(g.to_string())
		f.close()
		filename=self.label+'.'+format
		if format=='svg':  # annoyingly, dot creates svg's without units on the font size attribute.
			st=g.create_svg(prog=prog)
			st=re.sub(r"(font\-size\:[0-9]+\.*[0-9]*)([^p])",r"\1pt\2",st)
			f=open(filename,'w')
			f.write(st)
			f.close()
		else:
			g.write(filename,format=format,prog=prog)



	def generateMissingEntriesFromAbove(self, nodes):
		"""
		Generate a nonexisting entry in the library based on an averaging
		scheme.
		"""

		# Generate list of all sets of nodes that should have entries in the
		# library
		nodeLists = []
		for parent in nodes:
			nodeList = []
			for node in self.tree.children:
				temp = node
				while temp is not None and temp not in nodes:
					temp = self.tree.parent[temp]
				if temp == parent:
					nodeList.append(node)
			nodeLists.append(nodeList)
		nodesList = data.getAllCombinations(nodeLists)

		data = []
		for nodeList in nodesList:
			k = self.generateMissingEntryFromAbove(nodeList)
			if k is not None:
				data.append((nodeList, k))

		for nodeList, kinetics in data:
			self.library.add(nodeList, kinetics)
			#print nodeList, kinetics

	def generateMissingEntryFromAbove(self, nodes):

		# If an entry is already present, return it
		if self.library.getData(nodes) is not None:
			return self.library.getData(nodes)

		# Generate list of parents
		parentNodeLists = []
		for node in nodes:
			node0 = node
			parentNodeList = []
			while node0 is not None:
				parentNodeList.append(node0)
				node0 = self.tree.parent[node0]
			parentNodeLists.append(parentNodeList)
		parentNodesList = data.getAllCombinations(parentNodeLists)

		# Generate list of existing kinetic parameters along with "distance"
		# values (smaller is better)
		# This assumes that the amount of specificity indicated by moving from
		# parent to child in each tree is approximately equal
		minDistance = 1000000
		kinetics = []
		for parentNodes in parentNodesList:

			# Calculate distance
			distance = 0
			for i, node in enumerate(parentNodes):
				distance += parentNodeLists[i].index(node)

			# Don't bother if we've already found kinetics for a smaller
			# distance than the current distance, as we're just going to ignore
			# it anyway
			if distance < minDistance:
				# Get kinetics and append if available
				k = self.library.getData(parentNodes)
				if k is not None:
					kinetics.append((distance, k))
					if distance < minDistance: minDistance = distance

		# Fail if no kinetics found
		if len(kinetics) == 0:
			logging.warning('Unable to estimate missing kinetics for nodes %s in reaction family %s.' % (nodes, self.label))
			return None

		# Prune all entries with distances higher than the minimum
		kinetics = [k for d, k in kinetics if d == minDistance]
		if len(kinetics) == 0:
			logging.warning('Unable to estimate missing kinetics for nodes %s in reaction family %s.' % (nodes, self.label))
			return None

		# Average remaining kinetics
		kin = self.averageKinetics(kinetics)
		return kin

	def processLibraryData(self):
		"""
		Convert the data in the library from a string/unicode object to either
		an :class:`ArrheniusEPKinetics` object or a list of [link, comment]
		string pairs. This function is generally called in the course of
		loading a database from files.
		"""
		
		for label, item in self.library.iteritems():
			
			if item is None:
				pass
			elif not item.__class__ is tuple:
				raise data.InvalidDatabaseException('Kinetics library should be tuple at this point. Instead got %r'%data) 
			else: 
				index,data = item # break apart tuple, recover the 'index' - the beginning of the line in the library file.
				# Is't it dangerous having a local variable with the same name as a module?
				# what if we want to raise another data.InvalidDatabaseException() ?
				if not ( data.__class__ is str or data.__class__ is unicode) :
					raise data.InvalidDatabaseException('Kinetics library data format is unrecognized.')
				
				items = data.split()
				try:
					kineticData = [];
					# First item is temperature range
					kineticData.extend(items[0].split('-'))
					if len(kineticData) == 2:
						kineticData[0] = float(kineticData[0])
						kineticData[1] = float(kineticData[1])
					elif len(kineticData) == 1:
						kineticData = [float(items[0]), float(items[0])]
					# Middle items are Arrhenius + Evans-Polanyi data
					for i in range(1, 5):
						kineticData.append(float(items[i]))
					for i in range(5, 9):
						# Convert multiplicative uncertainties to additive
						# uncertainties if needed
						if items[i][0] == '*':
							kineticData.append((float(items[i][1:]) - 1.0) * float(items[i-4]))
						else:
							kineticData.append(float(items[i]))
					# Final item before comment is quality
					kineticData.append(int(items[9]))
					# Everything else is a comment
					comment = ' '.join(items[10:])

					kinetics = ArrheniusEPKinetics()
					kinetics.fromDatabase(kineticData, comment, len(self.template.reactants))
					kinetics.family = self
					kinetics.label = label
					kinetics.index = index
					#kinetics.comment = self.label + ' ' + label + ' ' + kinetics.comment
					self.library[label] = kinetics
					
				except (ValueError, IndexError), e:
					# Split data into link string and comment string; store
					# as list of length 2
					link = items[0]
					comment = data[len(link)+1:].strip()
					self.library[label] = [link, comment]
					
	def loadTemplate(self, path):
		"""
		Load and process a reaction template file located at `path`. This file
		is part of every reaction family.
		"""
		
		# Process the template file, removing comments and empty lines
		info = ''
		try:
			frec = open(path, 'r')
			for line in frec:
				line = data.removeCommentFromLine(line).strip()
				if len(line) > 0:
					info += line + '\n'
		except data.InvalidDatabaseException, e:
			logging.exception(str(e))
			return
		except IOError, e:
			logging.exception('Database file "' + e.filename + '" not found.')
			return
		finally:
			frec.close()
			
		lines = info.splitlines()
		
		# First line is 'Forward: <name of forward reaction>
		# Second line is 'Reverse: <name of reverse reaction>
		forward = ''; reverse = ''
		if lines[0].find('Forward:') > -1:
			self.label = lines[0][9:].strip()
		if lines[1].find('Reverse:') > -1:
			reverse = lines[1][9:].strip()
		
		# Third line is reaction template, of the form
		# A1 A2 ... + B1 B2 ... + ... <---> C1 C2 ... + D1 D2 ... + ...
		# A, B, ... are reactants; C, D, ... are products
		# A1, A2, ... represent different trees for the same species
		# The first tree of each species is always used to identify the 
		# reactants, so it should have all of the labeled atoms that are in that
		# reactant
		# The other trees can be used to provide functional group trees for
		# different parts of the molecule
		reactants = []; products = []; species = []
		atArrow = False
		items = lines[2].split(); items.extend('+')
		for item in items:
			if item == '+' or (item[0] == '<' and item[-1] == '>' and item.find('-') > -1):
				if len(species) == 1: species = species[0]
				if atArrow:		products.append(species)
				else:			reactants.append(species)
				species = []
				if item[0] == '<' and item[-1] == '>' and item.find('-') > -1:
					atArrow = True
			else:
				# Check that all reactant structures are in dictionary
				# The product structures are generated automatically and need
				# not be included
				if item not in self.dictionary and not atArrow:
					raise data.InvalidDatabaseException('Reaction family template contains an unknown structure.')
				species.append(item)

		# Set template reaction
		self.template = Reaction(reactants, products)

		# Remaining lines are reaction recipe for forward reaction
		self.recipe = ReactionRecipe()
		for line in lines[3:]:
			line = line.strip()
			
			# First item is the name of the action
			items = line.split()
			action = [ items[0].upper() ]
			if items[0] != 'CHANGE_BOND' and items[0] != 'FORM_BOND' and \
			  items[0] != 'BREAK_BOND' and items[0] != 'GAIN_RADICAL' and \
			  items[0] != 'LOSE_RADICAL':
				print items[0]

			# Remaining items are comma-delimited list of parameters enclosed by
			# {}, which we will split into individual parameters
			action.extend(line[len(items[0]):].strip()[1:-1].split(','))

			self.recipe.addAction(action)

		# Generate the reverse template
		if reverse != self.label:
			template = Reaction(self.template.products, self.template.reactants)
			self.reverse = ReactionFamily(reverse, template, self.recipe.getReverse())
			self.reverse.dictionary = self.dictionary
			self.reverse.tree = self.tree
			self.reverse.library = data.Library()
			self.reverse.forbidden = self.forbidden
			self.reverse.reverse = self

		# If necessary, generate the product template structure(s)
		# Don't need to do this if family is its own reverse
		if reverse != self.label:
			self.generateProductTemplate()

	def generateProductTemplate(self):
		"""
		Generate the product structures by applying the reaction template to
		the top-level nodes. For reactants defined by multiple structures, only
		the first is used here; it is assumed to be the most generic.
		"""

		# First, generate a list of reactant structures that are actual
		# structures, rather than unions
		reactantStructures = []
		for reactant in self.template.reactants:
			if isinstance(reactant, list):	reactants = [reactant[0]]
			else:							reactants = [reactant]
			unionFound = True
			while unionFound:
				unionFound = False; reactantsToRemove = []; reactantsToAdd = []
				for s in reactants:
					if self.dictionary[s] == 'union':
						reactantsToRemove.append(s)
						reactantsToAdd.extend(self.tree.children[s])
						unionFound = True
				for s in reactantsToRemove: reactants.remove(s)
				reactants.extend(reactantsToAdd)
			reactantStructures.append([self.dictionary[s] for s in reactants])
		
		# Second, get all possible combinations of reactant structures
		reactantStructures = data.getAllCombinations(reactantStructures)

		# Third, generate all possible product structures by applying the
		# recipe to each combination of reactant structures
		# Note that bimolecular products are split by labeled atoms
		productStructures = []
		for reactantStructure in reactantStructures:
			productStructure = self.applyRecipe(reactantStructure, unique=False)
			productStructures.append(productStructure)
		
		# Fourth, remove duplicates from the lists
		productStructureList = [[] for i in range(len(productStructures[0]))]
		for productStructure in productStructures:
			for i, struct in enumerate(productStructure):
				found = False
				for s in productStructureList[i]:
					if s.isIsomorphic(struct): found = True
				if not found:
					productStructureList[i].append(struct)

		# Fifth, associate structures with product template
		for i in range(len(self.template.products)):
			if len(productStructureList[i]) == 1:
				self.dictionary[self.template.products[i]] = productStructureList[i][0]
				self.tree.parent[self.template.products[i]] = None
				self.tree.children[self.template.products[i]] = []
			else:
				self.dictionary[self.template.products[i]] = 'union'
				children = []
				for j in range(len(productStructureList[i])):
					label = '%s_%i' % (self.template.products[i], j+1)
					self.dictionary[label] = productStructureList[i][j]
					children.append(label)
					self.tree.parent[label] = self.template.products[i]
					self.tree.children[label] = []

				self.tree.parent[self.template.products[i]] = None
				self.tree.children[self.template.products[i]] = children

	def reactantMatch(self, reactant, templateReactant):
		"""
		Return :data:`True` if the provided reactant matches the provided
		template reactant and :data:`False` if not.
		"""
		maps12 = []; maps21 = []
		if templateReactant.__class__ == list: templateReactant = templateReactant[0]
		struct = self.dictionary[templateReactant]
		if struct.__class__ == str or struct.__class__ == unicode:
			if struct.lower() == 'union':
				for child in self.tree.children[templateReactant]:
					ismatch, map21, map12 = self.reactantMatch(reactant, child)
					if ismatch:
						maps12.extend(map12); maps21.extend(map21)
		elif struct.__class__ == structure.Structure:
			return reactant.findSubgraphIsomorphisms(struct)

		return len(maps12) > 0, maps21, maps12

	def applyRecipe(self, reactantStructures, unique=True):
		"""
		Apply the recipe for this reaction family to the list of
		:class:`structure.Structure` objects `reactantStructures`. The atoms
		of the reactant structures must already be tagged with the appropriate
		labels. Returns a list of structures corresponding to the products
		after checking that the correct number of products was produced.
		"""

		# There is some hardcoding of reaction families in this function, so
		# we need the label of the reaction family for this
		label = self.label.lower()

		# Merge reactant structures into single structure
		# Also copy structures so we don't modify the originals
		# Since the tagging has already occurred, both the reactants and the
		# products will have tags
		reactantStructure = structure.Structure()
		for s in reactantStructures:
			reactantStructure = reactantStructure.merge(s.copy())

		# Hardcoding of reaction family for radical recombination (colligation)
		# because the two reactants are identical, they have the same tags
		# In this case, we must change the labels from '*' and '*' to '*1' and
		# '*2'
		if label == 'colligation':
			identicalCenterCounter = 0
			for atom in reactantStructure.atoms():
				if atom.label == '*':
					identicalCenterCounter += 1
					atom.label = '*' + str(identicalCenterCounter)
			if identicalCenterCounter != 2:
				raise Exception('Unable to change labels from "*" to "*1" and "*2" for reaction family %s.' % (label))

		# Generate the product structure by applying the recipe
		if not self.recipe.applyForward(reactantStructure, unique):
			return None
		productStructure = reactantStructure.copy()

		# Hardcoding of reaction family for reverse of radical recombination
		# (Unimolecular homolysis)
		# Because the two products are identical, they should the same tags
		# In this case, we must change the labels from '*1' and '*2' to '*' and
		# '*'
		if label == 'unimolecular homolysis':
			for atom in productStructure.atoms():
				if atom.label == '*1' or atom.label == '*2': atom.label = '*'
			
		
		# If reaction family is its own reverse, relabel atoms
		if not self.reverse:
			# Get atom labels for products
			atomLabels = {}
			for atom in productStructure.atoms():
				if atom.label != '':
					atomLabels[atom.label] = atom
			
			# This is hardcoding of reaction families (bad!)
			label = self.label.lower()
			if label == 'h abstraction':
				# '*2' is the H that migrates
				# it moves from '*1' to '*3'
				atomLabels['*1'].label = '*3'
				atomLabels['*3'].label = '*1'
			
			elif label == 'intra h migration':
				# '*3' is the H that migrates
				# swap the two ends between which the H moves
				atomLabels['*1'].label = '*2'
				atomLabels['*2'].label = '*1'
				# reverse all the atoms in the chain between *1 and *2
				# i.e. swap *4 with the highest, *5 with the second-highest
				highest = len(atomLabels)
				if highest>4:
					for i in range(4,highest+1):
						atomLabels['*%d'%i].label = '*%d'%(4+highest-i)
		
		
		# Split product structure into multiple species if necessary
		if len(self.template.products) > 1:
			productStructures = productStructure.split()
		else:
			productStructures = [productStructure]
		
		# Make sure we've made the expected number of products
		if len(self.template.products) != len(productStructures):
			# We have a different number of products than expected by the template.
			# It might be because we found a ring-opening using a homolysis template
			if (label=='unimolecular homolysis'
			 and len(productStructures) == 1
			 and len(reactantStructures) == 1):
				# just be absolutely sure (maybe slow, but safe)
				rs = reactantStructures[0]
				if ( rs.graph.isVertexInCycle(rs.getLabeledAtom('*1'))
				 and rs.graph.isVertexInCycle(rs.getLabeledAtom('*2'))):
					# both *1 and *2 are in cycles (probably the same one)
					# so it's pretty safe to just fail quietly,
					# and try the next reaction
					return None
			
			# no other excuses, raise an exception
			message = 'Application of reaction recipe failed; expected %s product(s), but %s found.\n' % (len(self.template.products), len(productStructures))
			message += "Reaction family: %s \n"%str(self)
			message += "Reactant structures: %s \n"%reactantStructures
			message += "Product structures: %s \n"%productStructures
			message += "Template: %s"%self.template
			logging.error(message)
			return None # don't fail!!! muhahaha
			raise Exception(message)

		# If there are two product structures, place the one containing '*1' first
		if len(productStructures) == 2:
			if not productStructures[0].containsLabeledAtom('*1') and \
				productStructures[1].containsLabeledAtom('*1'):
				productStructures.reverse()

		# reset any cached structure information because it is now invalid
		for struct in productStructures:
			struct.resetCachedStructureInfo()

		# Return the product structures
		return productStructures
	
	def makeReaction(self, reactants, reactantStructures, maps):
		"""
		Create a reaction involving a list of `reactants`. The `reactantStructures`
		parameter is a list of structures in the order the reactants are stored
		in the reaction family template, and the `maps` parameter is a list of
		mappings of the top-level tree node of each template reactant to the
		corresponding structure.
		"""
		
		# Clear any previous atom labeling from all reactant structures
		for struct in reactantStructures: struct.clearLabeledAtoms()
		
		# Tag atoms with labels
		for map in maps:
			for templateAtom, reactantAtom in map.iteritems():
				reactantAtom.label = templateAtom.label
		
		# Generate the product structures by applying the forward reaction recipe
		try:
			productStructures = self.applyRecipe(reactantStructures)
			if not productStructures: return None
		except chem.InvalidChemicalActionException, e:
			print 'Unable to apply reaction recipe!'
			print 'Reaction family is %s' % self
			print 'Reactant structures are:'
			for struct in reactantStructures:
				print struct.toAdjacencyList()
			raise e

		# Check that reactant and product structures are allowed in this family
		# If not, then stop
		if self.forbidden is not None:
			for label, struct2 in self.forbidden.iteritems():
				for struct in reactantStructures:
					if struct.isSubgraphIsomorphic(struct2): return None
				for struct in productStructures:
					if struct.isSubgraphIsomorphic(struct2): return None
		
		# Convert structure(s) to products
		products = []
		for product in productStructures:
			spec = species.makeNewSpecies(product)
			# Don't make a new reaction if no species was returned from
			# makeNewSpecies() (e.g. due to forbidden structure)
			if spec is None: return None
			products.append(spec)
		
		# Create reaction and add if unique
		rxn, isNew = makeNewReaction(reactants[:], products, reactantStructures, productStructures, self)
		if isNew:	return rxn
		else:		return None
	
	def getReactionList(self, reactants):
		"""
		Generate a list of all of the possible reactions of this family between
		the list of `reactants`.
		"""
		rxnList = []
		# If the number of reactants provided does not match the number of
		# reactants in the template, return False
		if len(reactants) == 1 and self.template.isUnimolecular():

			# Iterate over all resonance isomers of the reactant
			for structure in reactants[0].structure:

				ismatch, map21, map12 = self.reactantMatch(structure, self.template.reactants[0])
				if ismatch:
					for map in map12:
						rxn = self.makeReaction(reactants, [structure], [map])
						if rxn is not None:
							rxnList.append(rxn)
			
		# Bimolecular reactants: A + B --> products
		elif len(reactants) == 2 and self.template.isBimolecular():
			
			# Make copies of the structure lists of the two reactants
			# This is a workaround for an issue in which the two reactant
			# structure lists were getting swapped around, resulting in
			# unbalanced reactions
			# The copy is needed for cases where A and B are the same
			structuresA = []; structuresB = []
			for structureA in reactants[0].structure:
				structuresA.append(structureA.copy())
			for structureB in reactants[1].structure:
				structuresB.append(structureB.copy())
			
			# Iterate over all resonance isomers of the reactant
			for structureA in structuresA:
				for structureB in structuresB:
				
					# Reactants stored as A + B
					ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[0])
					ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[1])
					
					# Iterate over each pair of matches (A, B)
					if ismatch_A and ismatch_B:
						for mapA in map12_A:
							for mapB in map12_B:
								rxn = self.makeReaction(reactants, [structureA, structureB], [mapA, mapB])
								if rxn is not None:
									rxnList.append(rxn)
									
					# Only check for swapped reactants if they are different
					if reactants[0].id != reactants[1].id:
						
						# Reactants stored as B + A
						ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[1])
						ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[0])
						
						# Iterate over each pair of matches (A, B)
						if ismatch_A and ismatch_B:
							for mapA in map12_A:
								for mapB in map12_B:
									rxn = self.makeReaction(reactants, [structureB, structureA], [mapB, mapA])
									if rxn is not None:
										rxnList.append(rxn)
		
		return rxnList
	
	def getKinetics(self, reaction, structures):
		"""
		Determine the appropriate kinetics for `reaction` which involves the
		labeled atoms in `atoms`.
		"""
		
		# Get forward reaction template and remove any duplicates
		forwardTemplate, reverseTemplate = self.getTemplateLists()
		#forwardTemplate = list(set(forwardTemplate)) # this can shuffle the order!
		temporary=[]
		symmetric_tree=False
		for node in forwardTemplate:
			if node not in temporary:
				temporary.append(node)
			else: 
				# duplicate node found at top of tree
				# eg. R_recombination: ['Y_rad', 'Y_rad']
				assert len(forwardTemplate)==2 , 'Can currently only do symmetric trees with nothing else in them'
				symmetric_tree = True
		forwardTemplate = temporary
		
		# Descend reactant trees as far as possible
		template = []
		for forward in forwardTemplate:
			# 'forward' is a head node that should be matched.
			# Get labeled atoms of forward
			node = forward
			group = self.dictionary[node]
			# to sort out "union" groups:
			# descends to the first child that's not a union
			while isinstance(group,str) or isinstance(group,unicode):
				node = self.tree.children[node][0]
				group = self.dictionary[node]
			# ...but this child may not match the structure.
			# eg. an R3 ring node will not match an R4 ring structure.
			# (but at least the first such child will contain fewest labels - we hope)
			
			atomList = group.getLabeledAtoms() # list of atom labels in highest non-union node
			
			for struct in structures:
				# Match labeled atoms
				# Check this structure has each of the atom labels in this group
				has_all_atom_labels = True
				for label in atomList:
					if not struct.containsLabeledAtom(label):
						has_all_atom_labels = False
						#print "%s does not contain %s atom labels"%(struct.toSMILES(),node)
						break
						# this structure (reactant/product) does not contain this group's atom labels
				if not has_all_atom_labels: continue # don't try to match this structure - the atoms aren't there!
				
				# Match structures
				atoms = struct.getLabeledAtoms()
				matched_node = self.descendTree(struct, atoms, forward)
				if matched_node is not None:
					template.append(matched_node)
				else:
					logging.warning("Couldn't find match for %s in %s"%(forward,atomList))
					logging.warning( struct.toAdjacencyList() )
					
		# Get fresh templates (with duplicate nodes back in)
		forwardTemplate, reverseTemplate = self.getTemplateLists()
		
		# Check that we were able to match the template.
		# template is a list of the actual matched nodes
		# forwardTemplate is a list of the top level nodes that should be matched
		if len(template) != len(forwardTemplate):
			logging.warning('Warning: Unable to find matching template for reaction %s in reaction family %s' % (str(reaction), str(self)) )
			logging.warning(" Trying to match",forwardTemplate)
			logging.warning(" Matched",template)
			raise UndeterminableKineticsException(reaction)
			print str(self), template, forwardTemplate, reverseTemplate
			for reactant in reaction.reactants:
				print reactant.toAdjacencyList() + '\n'
			for product in reaction.products:
				print product.toAdjacencyList() + '\n'
				
			## If unable to match template, use the most general template
			#template = forwardTemplate
		
		
#		k = self.library.getData(template)
#		print template, k
#		if k is not None: return [k]
#		else: return None
		
		
		# climb the tree finding ancestors
		nodeLists = []
		for temp in template:
			nodeList = []
			while temp is not None:
				nodeList.append(temp)
				temp = self.tree.parent[temp]
			nodeLists.append(nodeList)
		
		# Generate all possible combinations of nodes
		items = data.getAllCombinations(nodeLists)
		
		# Generate list of kinetics at every node
		#logging.debug("   Template contains %s"%forwardTemplate)
		kinetics = []
		for item in items:
			itemData = self.library.getData(item)
			#logging.debug("   Looking for %s found %r"%(item, itemData))
			if itemData is not None:
				kinetics.append(itemData)
				
			if symmetric_tree: # we might only store kinetics the other way around
				item.reverse()
				itemData = self.library.getData(item)
				#logging.debug("   Also looking for %s found %r"%(item, itemData))
				if itemData is not None:
					kinetics.append(itemData)
					
		if len(kinetics) == 0: return None
		
		return kinetics
	
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
	import unirxn.states as states

	if kinetics.Ea < 0.0:
		logging.warning('Negative activation energy of %s kJ/mol encountered during unirxn calculation; setting to zero.' % (kinetics.Ea / 1000.0))
		Ea = 0.0
	else:
		Ea = kinetics.Ea

	dE = Elist[1] - Elist[0]
	k = numpy.zeros(len(Elist), numpy.float64)

	if kinetics.n == 0.0:
		# Determine the microcanonical rate directly
		s = int(math.floor(Ea / dE))
		for r in range(len(Elist)):
			if Elist[r] > E0 and densStates[r] != 0:
				k[r] = kinetics.A * densStates[r - s] / densStates[r]

	elif kinetics.n > 0.0:
		import scipy.special
		# Evaluate the inverse Laplace transform of the T**n piece, which only
		# exists for n >= 0
		phi = numpy.zeros(len(Elist), numpy.float64)
		for i, E in enumerate(Elist):
			if E == 0.0:
				phi[i] = 0.0
			else:
				phi[i] = E**(kinetics.n-1) / (constants.R**kinetics.n * scipy.special.gamma(kinetics.n))
		# Evaluate the convolution
		states.convolve(phi, densStates, Elist)
		# Apply to determine the microcanonical rate
		s = int(math.floor(Ea / dE))
		for r in range(len(Elist)):
			if Elist[r] > E0 and densStates[r] != 0:
				k[r] = kinetics.A * phi[r - s] / densStates[r]

	else:
		# Use the cheating method for n < 0
		s = int(math.floor(Ea / dE))
		for r in range(len(Elist)):
			if Elist[r] > E0 and densStates[r] != 0:
				k[r] = kinetics.A * T**kinetics.n * densStates[r - s] / densStates[r]

	return k

################################################################################

class ReactionFamilySet:
	"""
	Represent a set of reaction families. The `families` attribute stores a
	dictionary of :class:`ReactionFamily` objects representing the families in
	the set.
	"""
	
	def __init__(self):
		self.families = {}
	
	def load(self, datapath, only_families=False):
		"""
		Load a set of reaction families from the general database
		specified at `datapath`. If only_families is present, families not in
		this list will not be loaded (e.g. only_families=['H_Abstraction'] )
		"""
		
		datapath = os.path.abspath(datapath)

		logging.info('Loading reaction family databases from %s...' % datapath)
		
		# Load the families from kinetics/families.txt
		familyList = []
		try:
			ffam = open(os.path.join(datapath,'kinetics','families.txt'), 'r')
			for line in ffam:
				line = data.removeCommentFromLine(line).strip()
				if len(line) > 0:
					items = line.split()
					items[0] = int(items[0])
					familyList.append(items)
		except data.InvalidDatabaseException, e:
			logging.exception(str(e))
			return
		except IOError, e:
			logging.exception('Database file "' + e.filename + '" not found.')
			return
		finally:
			ffam.close()
		
		# Load the reaction families (if they exist and status is 'on')
		self.families = {}
		for index, status, label in familyList:
			path = os.path.join(datapath, 'kinetics', label)
			if os.path.isdir(path) and status.lower() == 'on':

				# skip families not in only_families, if it's set
				if only_families and label not in only_families: continue

				logging.info('Loading reaction family %s from %s...' % (label, datapath))
				family = ReactionFamily(label)
				family.load(path)
				self.families[family.label] = family
				if family.reverse is not None:
					self.families[family.reverse.label] = family.reverse

	def getReactions(self, species):
		"""
		Generate a list of reactions that involve a list of one or two `species`
		as a reactant or product.
		"""
		
		rxnList = []

		# Don't bother if any or all of the species are marked as nonreactive
		if not all([spec.reactive for spec in species]):
			return rxnList

		log_text = ' + '.join([str(spec) for spec in species])
		
		logging.info('Looking for reactions of %s'%(log_text))
		
		for key, family in self.families.iteritems():
			rxnList.extend(family.getReactionList(species))

		if len(rxnList) == 1:
			logging.info('Found %s reaction for %s'%(len(rxnList), log_text))
		else:
			logging.info('Found %s reactions for %s'%(len(rxnList), log_text))
		
		return rxnList


kineticsDatabase = None

################################################################################

class ReactionException(Exception):
	"""
	An base exception for reactions.
	Takes a reaction object, and optional message
	"""
	def __init__(self, reaction, message=''):
		self.reaction = reaction
		self.message = message
		
	def __str__(self):
		string = "Reaction: "+str(self.reaction) + '\n'
		string += "Reaction Family: "+str(self.reaction.family) + '\n'
		for reactant in self.reaction.reactants:
			string += reactant.toAdjacencyList() + '\n'
		for product in self.reaction.products:
			string += product.toAdjacencyList() + '\n'
		if self.message: string += "Message: "+self.message
		return string


class UndeterminableKineticsException(ReactionException):
	"""
	An exception raised when attempts to select appropriate kinetic parameters
	for a chemical reaction are unsuccessful.
	"""
	def __init__(self, reaction, message=''):
		new_message = 'Kinetics could not be determined. '+message
		ReactionException.__init__(self,reaction,new_message)


################################################################################

# The global list of reactions created at any point during RMG execution
# The list is stored in reverse of the order in which the reactions are created;
# when searching the list, it is more likely to match a recently created
# reaction than an older reaction
reactionList = []

def makeNewReaction(reactants, products, reactantStructures, productStructures, family):
	"""
	Attempt to make a new reaction based on a list of `reactants` and a list of
	`products`. The combination of these and a reaction `family` string uniquely
	identifies a reaction. The reactant and product lists must contain 
	:class:`Species` objects, not :class:`Structure` objects.

	The proposed reaction is checked against the list of
	existing reactions; if the reaction already exists, this function returns
	the existing reaction. If the reaction does not exist, a :class:`Reaction`
	object is created and returned after being appended to the global reaction
	list.
	"""
	
	# Sort reactants and products (to make comparisons easier/faster)
	reactants.sort()
	products.sort()
	
	# Check that the reaction actually results in a different set of species
	if reactants == products:
		return None, False
	
	# If a species appears in both the reactants and products, then remove it
	speciesToRemove = []
	for spec in reactants:
		if spec in products:
			speciesToRemove.append(spec)
	for spec in speciesToRemove:
		reactants.remove(spec)
		products.remove(spec)
	if len(reactants) == 0 or len(products) == 0:
		return None, False

	# Check that the reaction is unique
	matchReaction = None
	for rxn in reactionList:
		if isinstance(rxn.family, ReactionFamily) and (
			rxn.family.label != family.label and rxn.family.reverse.label != family.label):
			# rxn is not from seed, and families are different
			continue # not a match, try next rxn
		if (rxn.reactants == reactants and rxn.products == products) or \
			(rxn.reactants == products and rxn.products == reactants):
			matchReaction = rxn
			break # found a match so stop checking other rxn
	
	# If a match was found, take an
	if matchReaction is not None:
		#matchReaction.multiplier += 1.0
		return matchReaction, False
	
	# If this point is reached, the proposed reaction is new, so make new
	# Reaction objects for forward and reverse reaction
	forward = Reaction(reactants, products, family)
	reverseFamily = None
	if family is not None: reverseFamily = family.reverse or family
	reverse = Reaction(products, reactants, reverseFamily)
	forward.reverse = reverse
	reverse.reverse = forward
	
	# Get atom labels of reactants
	reactantLabels = {}; productLabels = {}
	for structure in reactantStructures:
		for atom in structure.atoms():
			if atom.label == '*': 
				if atom.label in reactantLabels: 
					reactantLabels[atom.label].append(atom)
				else:
					reactantLabels[atom.label] = [atom]
			elif atom.label != '': 
				reactantLabels[atom.label] = atom
	# Get atom labels of products
	for structure in productStructures:
		for atom in structure.atoms():
			if atom.label == '*':
				if atom.label in productLabels:
					productLabels[atom.label].append(atom)
				else:
					productLabels[atom.label] = [atom]
			elif atom.label != '':
				productLabels[atom.label] = atom
			
	# Dictionaries containing the labeled atoms for the reactants and products
	forward.atomLabels = reactantLabels
	reverse.atomLabels = productLabels
	
	if forward.family is None or reverse.family is None:
		reactionList.insert(0, forward)
		return forward, True
	
	# Attempt to get the kinetics of the forward and reverse reactions
	forwardKinetics = forward.family.getKinetics(forward, reactantStructures)
	reverseKinetics = reverse.family.getKinetics(reverse, productStructures)
	
	# By convention, we only work with the reaction in the direction for which
	# we have assigned kinetics from the kinetics database; the kinetics of the
	# reverse of that reaction come from thermodynamics
	# If we have assigned kinetics in both directions, then (for now) choose the
	# kinetics for the forward reaction
	rxn = forward
	if forwardKinetics is not None and reverseKinetics is not None:
		rxn = forward
		reverseKinetics = []
	elif forwardKinetics is not None:
		rxn = forward
	elif reverseKinetics is not None:
		rxn = reverse
	else:
		raise UndeterminableKineticsException(forward)
		return None, False
	
	forward.kinetics = forwardKinetics
	reverse.kinetics = reverseKinetics

	# Note in the log
	logging.verbose('Creating new %s reaction %s' % (rxn.family, rxn))

	return processNewReaction(rxn)

def processNewReaction(rxn):
	"""
	Once a reaction `rxn` has been created (e.g. via :meth:`makeNewReaction`),
	this function handles other aspects	of preparing it for RMG.
	"""

	reactionList.insert(0, rxn)

	# Return newly created reaction
	return rxn, True

################################################################################

if __name__ == '__main__':
	
	import os.path
	import main
	main.initializeLog(logging.DEBUG)

	import thermo

	datapath = '../data/RMG_database/'

	logging.debug('General database: ' + os.path.abspath(datapath))
	species.thermoDatabase = species.ThermoDatabaseSet()
	species.thermoDatabase.load(datapath)
	kineticsDatabase = ReactionFamilySet()
	kineticsDatabase.load(datapath)

	structure1 = structure.Structure()
	structure1.fromAdjacencyList("""HXD13
1 C 0 {2,D} {7,S} {8,S}
2 C 0 {1,D} {3,S} {9,S}
3 C 0 {2,S} {4,D} {10,S}
4 C 0 {3,D} {5,S} {11,S}
5 C 1 {4,S} {6,S} {12,S}
6 C 0 {5,S} {14,S} {15,S} {16,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {3,S}
11 H 0 {4,S}
12 H 0 {5,S}
14 H 0 {6,S}
15 H 0 {6,S}
16 H 0 {6,S}
""")

	species1 = species.makeNewSpecies(structure1, 'C6H9J', True)

	#print len(species1.structure)

	structure2 = structure.Structure()
	structure2.fromSMILES('[H][H]')
	species2 = species.makeNewSpecies(structure2, 'H2', True)

	rxnList = kineticsDatabase.getReactions([species1])
	#rxnList = kineticsDatabase.getReactions([species1, species2])
	for rxn in rxnList:
		print rxn
