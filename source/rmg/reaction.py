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
	
	def __init__(self, reactants=[], products=[], kinetics=None):
		"""Initialize a reaction object."""
		self.reactants = reactants
		self.products = products
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

################################################################################

if __name__ == '__main__':
	
	reactant1 = chem.Species('CH4')
	
	product1 = chem.Species('CH3')
	product2 = chem.Species('H')
	
	reaction = Reaction([reactant1], [product1, product2], None)
	print reaction
	