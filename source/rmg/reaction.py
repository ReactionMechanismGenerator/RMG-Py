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
Contains classes describing chemical reactions.
"""

import quantities as pq
import logging
import os

################################################################################

class Reaction:
	"""
	Represent a generic chemical reaction. A reaction is defined by the set of
	reactants, set of products, and the transition state; thus both the forward
	and reverse reaction are represented by a single object since they share 
	the same transition state. By convention, the forward reaction is taken to
	be that for which the provided kinetics apply; the reverse kinetics are 
	taken from thermodynamic reversibility.
	"""
	
	def __init__(self, reactants=None, products=None, family=None, kinetics=None):
		"""Initialize a reaction object."""
		self.reactants = reactants or []
		self.products = products or []
		self.family = family
		self.kinetics = kinetics
	
	def __str__(self):
		"""
		Return a string representation of the reaction.
		"""
		string = ''
		for reactant in self.reactants:
			string += str(reactant) + ' + '
		string = string[0:-3] + ' <---> '
		for product in self.products:
			string += str(product) + ' + '
		string = string[0:-3]
		return string

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

################################################################################

# The global list of reactions created at any point during RMG execution
# The list is stored in reverse of the order in which the reactions are created;
# when searching the list, it is more likely to match a recently created
# reaction than an older reaction
reactionList = []

def makeNewReaction(reactants, products, family=''):
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
	if len(reactants) == len(products):
		match = True
		for i in range(len(reactants)):
			if reactants[i] != products[i]: match = False
		if match: return None, False

	# Check that the reaction is unique
	for rxn in reactionList:
		# Check forward reaction for match
		if rxn.family is family or rxn.family is None or family is None:
			if len(rxn.reactants) == len(reactants) and len(rxn.products) == len(products):
				match = True
				for i in range(len(reactants)):
					if rxn.reactants[i] != reactants[i]: match = False
				for i in range(len(products)):
					if rxn.products[i] != products[i]: match = False
				if match: return rxn, False
		# Check reverse reaction for match
		if rxn.family.reverse is family or rxn.family.reverse is None or family is None:
			if len(rxn.reactants) == len(products) and len(rxn.products) == len(reactants):
				match = True
				for i in range(len(reactants)):
					if rxn.products[i] != reactants[i]: match = False
				for i in range(len(products)):
					if rxn.reactants[i] != products[i]: match = False
				if match: return rxn, False

	# If this point is reached, the proposed reaction is new, so make new
	# Reaction object and append to global reaction list
	rxn = Reaction(reactants, products, family)
	reactionList.insert(0, rxn)
	
	# Note in the log
	logging.debug('Created new reaction ' + str(rxn))

	# Return newly created reaction
	return rxn, True

################################################################################

if __name__ == '__main__':
	
	reactant1 = chem.Species('CH4')
	
	product1 = chem.Species('CH3')
	product2 = chem.Species('H')
	
	reaction = Reaction([reactant1], [product1, product2], None)
	print reaction
	