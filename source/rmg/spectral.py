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

import data

################################################################################

class CharacteristicFrequency:
	"""
	Represent a characteristic frequency in the frequency database. The
	characteristic frequency has a real lower bound `lower`, a real upper bound
	`upper`, and an integer `degeneracy`.
	"""

	def __init__(self, lower=0.0, upper=0.0, degeneracy=1):
		self.lower = lower
		self.upper = upper
		self.degeneracy = degeneracy

	def generateFrequencies(self, count=1):
		"""
		Generate a set of frequencies. The number of frequencies returned is
		self.degeneracy * count, and these are distributed linearly between
		`self.lower` and `self.upper`.
		"""
		numFreqs = self.degeneracy * count

################################################################################

class FrequencyDatabase(data.Database):
	"""
	Represent an RMG frequency database.
	"""

	def __init__(self):
		"""
		Call the generic `data.Database.__init__()` method. This in turn creates
			* self.dictionary = Dictionary()
			* self.library = Library()
			* self.tree = Tree()
		"""
		data.Database.__init__(self)


	def load(self, dictstr, treestr, libstr):
		"""
		Load a group frequency database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. The tree is optional, and should
		be set to '' if not desired.
		"""

		# Load dictionary, library, and (optionally) tree
		data.Database.load(self, dictstr, treestr, libstr)

		# Convert data in library to ThermoData objects or lists of
		# [link, comment] pairs
		for label, item in self.library.iteritems():

			if item is None:
				pass
			elif not item.__class__ is tuple:
				raise data.InvalidDatabaseException('Frequencies library should be tuple at this point. Instead got %r'%data)
			else:
				index, item = item
				if not (item.__class__ is str or item.__class__ is unicode):
					raise data.InvalidDatabaseException('Frequencies library data format is unrecognized.')

				items = item.split()

				# Items should be a multiple of three (no comments allowed at the moment)
				if len(items) % 3 != 0:
					raise data.InvalidDatabaseException('Unexpected number of items encountered in frequencies library.')

				# Convert list of data into a list of characteristic frequencies
				frequencies = []
				count = len(items) / 3
				for i in range(count):
					frequency = CharacteristicFrequency(
						lower=float(items[3*i+1]),
						upper=float(items[3*i+2]),
						degeneracy=int(items[3*i]))
					frequencies.append(frequency)

				self.library[label] = frequencies

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (dictstr))

frequencyDatabase = None

################################################################################
