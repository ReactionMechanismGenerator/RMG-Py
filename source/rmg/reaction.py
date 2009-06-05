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

################################################################################

class ArrheniusEPKinetics:
	"""
	Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
	kinetic expression has the form
	
	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_\mathrm{a}}{RT} \\right)
		
	The parameter :math:`\\alpha` is used to correct the activation energy 
	:math:`E_\\mathrm{a}` via the Evans-Polanyi formula
	
	.. math:: E_\\mathrm{a} = E_\\mathrm{a}^0 - (\\alpha - 1) \\Delta H_\\mathrm{rxn} 
	
	"""

	def __init__(self, A=0.0, Ea=0.0, n=0.0, alpha=0.0):
		self.A = A
		self.Ea = Ea
		self.n = n
		self.alpha = alpha
	
	def fromDatabase(self, data, comment):
		"""
		Process a list of numbers `data` and associated description `comment`
		generated while reading from a kinetics database.
		"""
	
		if len(data) != 11:
			raise Exception('Invalid list of kinetic data; should be a list of numbers of length 10.')
		
		Tmin, Tmax, A, n, alpha, Ea, dA, dn, dalpha, dEa, rank = data
		
		self.Trange = pq.Quantity([Tmin, Tmax], 'K')
		
		self.A = pq.UncertainQuantity(A, 's^-1', dA)
		self.Ea = pq.UncertainQuantity(Ea, 'kcal/mol', dEa)
		#self.Ea.units = 'J/mol'
		self.n = pq.UncertainQuantity(n, '', dn)
		self.alpha = pq.UncertainQuantity(alpha, 's^-1', dalpha)
		
		self.rank = rank
		self.comment = comment

################################################################################

if __name__ == '__main__':
	
	reactant1 = chem.Species('CH4')
	
	product1 = chem.Species('CH3')
	product2 = chem.Species('H')
	
	reaction = Reaction([reactant1], [product1, product2], None)
	print reaction
	