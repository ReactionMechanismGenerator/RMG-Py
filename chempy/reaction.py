#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains classes and functions for working with chemical reactions.

From the `IUPAC Compendium of Chemical Terminology 
<http://dx.doi.org/10.1351/goldbook>`_, a chemical reaction is "a process that 
results in the interconversion of chemical species".

In ChemPy, a chemical reaction is called a Reaction object and is represented in
memory as an instance of the :class:`Reaction` class.
"""

################################################################################

class Reaction:
	"""
	A chemical reaction.
	
	=============== =========================== ================================
	Attribute       Type                        Description
	=============== =========================== ================================
	`index`         :class:`int`                A unique nonnegative integer index
	`reactants`     :class:`list`               The reactant species (as :class:`Species` objects)
	`products`      :class:`list`               The product species (as :class:`Species` objects)
	=============== =========================== ================================
	
	"""
	
	def __init__(self, index=-1, reactants=None, products=None):
		self.index = index
		self.reactants = reactants
		self.products = products
	
	def __repr__(self):
		"""
		Return a string representation of the reaction, suitable for console output.
		"""
		return "<Reaction %i '%s'>" % (self.index, str(self.label))
	
	def __str__(self):
		"""
		Return a string representation of the reaction, in the form 'A + B <=> C + D'.
		"""
		return ' <=> '.join([' + '.join([str(s) for s in self.reactants]), ' + '.join([str(s) for s in self.products])])
	
