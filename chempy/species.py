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
This module contains classes and functions for working with chemical species.

From the `IUPAC Compendium of Chemical Terminology 
<http://dx.doi.org/10.1351/goldbook>`_, a chemical species is "an 
ensemble of chemically identical molecular entities that can explore the same 
set of molecular energy levels on the time scale of the experiment". This
definition is purposefully vague to allow the user flexibility in application.

In ChemPy, a chemical species is called a Species object and is represented in
memory as an instance of the :class:`Species` class.
"""

################################################################################

class Species:
	"""
	A chemical species.
	
	=============== =========================== ================================
	Attribute       Type                        Description
	=============== =========================== ================================
	`index`         :class:`int`                A unique nonnegative integer index
	`label`         :class:`str`                A descriptive string label
	=============== =========================== ================================
	
	"""
	
	def __init__(self, index=-1, label=''):
		self.index = index
		self.label = label
	
	def __repr__(self):
		"""
		Return a string representation of the species, suitable for console output.
		"""
		return "<Species %i '%s'>" % (self.index, self.label)
	
	def __str__(self):
		"""
		Return a string representation of the species, in the form 'label(id)'.
		"""
		return '%s(%i)' % (self.label, self.index)
	