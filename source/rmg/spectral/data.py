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

import logging

import rmg.data as data
import rmg.constants as constants
import rmg.thermo as thermo

from modes import *

################################################################################

class CharacteristicFrequency:
	"""
	Represent a characteristic frequency in the frequency database. The
	characteristic frequency has a real lower bound `lower`, a real upper bound
	`upper`, an integer `degeneracy`.
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
		from numpy import linspace
		number = self.degeneracy * count
		if number == 1:
			return [(self.lower + self.upper) / 2.0]
		else:
			return linspace(self.lower, self.upper, number, endpoint=True)

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

				# First item is the symmetry correction 
				frequencies = [int(items.pop(0))]

				# Items should be a multiple of three (no comments allowed at the moment)
				if len(items) % 3 != 0:
					raise data.InvalidDatabaseException('Unexpected number of items encountered in frequencies library.')

				# Convert list of data into a list of characteristic frequencies
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

def loadFrequencyDatabase(databasePath):
	"""
	Create and load the frequencies databases.
	"""
	import os.path
	import logging

	databasePath = os.path.join(databasePath, 'frequencies')

	# Create and load thermo databases
	database = FrequencyDatabase()
	logging.debug('\tFrequencies database')
	database.load(
		dictstr=os.path.join(databasePath, 'Dictionary.txt'),
		treestr=os.path.join(databasePath, 'Tree.txt'),
		libstr=os.path.join(databasePath, 'Library.txt'))

	return database

################################################################################

def generateSpectralData(struct, thermoData):
	"""
	Generate the spectral data for a :class:`structure.Structure` object
	`struct` with corresponding :class:`thermo.ThermoGAData` object `thermo`.
	The group frequency method is used to do so; this method has two steps:

	1.	Search the structure for certain functional groups for which
		characteristic frequencies are known, and use those frequencies.
	2.	For any remaining degrees of freedom, fit the parameters such that they
		replicate the heat capacity.

	This method only fits the internal modes (i.e. vibrations and hindered
	rotors).
	"""

	# No spectral data for single atoms
	if len(struct.atoms()) < 2:
		return None

	# Determine linearity of structure
	linear = struct.isLinear()

	# Determine the number of internal modes for this structure
	numModes = 3 * len(struct.atoms()) - (5 if linear else 6)

	# Determine the number of vibrational and hindered rotor modes for this
	# structure
	numRotors = struct.countInternalRotors()
	numVibrations = numModes - numRotors
	#print 'For %s, I found %i internal rotors and %i vibrations for a total of %i modes' % (struct, numRotors, numVibrations, numModes)

	# For each group in library, find all subgraph isomorphisms
	groupCount = {}
	for node, data in frequencyDatabase.library.iteritems():
		ismatch, map12List, map21List = struct.findSubgraphIsomorphisms(frequencyDatabase.dictionary[node])
		if ismatch:
			count = len(map12List)
		else:
			count = 0
		if count % data[0] != 0:
			raise Exception('Incorrect number of matches of node "%s" while estimating frequencies of %s.' % (node, struct))
		groupCount[node] = count / data[0]

	# For debugging, print a list of the groups found
	#print 'Groups found:'
	#for node, count in groupCount.iteritems():
	#	if count != 0: print '\t', node, count

	# Get characteristic frequencies
	frequencies = []
	for node, count in groupCount.iteritems():
		for charFreq in frequencyDatabase.library[node][1:]:
			frequencies.extend(charFreq.generateFrequencies(count))

	# Check that we have the right number of degrees of freedom specified
	if len(frequencies) > numVibrations:
		# We have too many vibrational modes
		difference = len(frequencies) - numVibrations
		# First try to remove hindered rotor modes until the proper number of modes remains
		if numRotors >= difference:
			numRotors -= difference
			logging.warning('For structure %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i internal rotors to compensate.' % (struct, difference))
		# If that doesn't work, turn off functional groups until the problem is underspecified again
		else:
			groupsRemoved = 0
			freqsRemoved = 0
			freqCount = len(frequencies)
			while freqCount > numVibrations:
				minDegeneracy, minNode = min([(sum([charFreq.degeneracy for charFreq in frequencyDatabase.library[node][1:]]), node) for node in groupCount])
				if groupCount[minNode] > 1:
					groupCount[minNode] -= 1
				else:
					del groupCount[minNode]
				groupsRemoved += 1
				freqsRemoved += minDegeneracy
				freqCount -= minDegeneracy
			# Log warning
			logging.warning('For structure %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i groups (%i frequencies) to compensate.' % (struct, groupsRemoved, freqsRemoved))
			# Regenerate characteristic frequencies
			frequencies = []
			for node, count in groupCount.iteritems():
				for charFreq in frequencyDatabase.library[node][1:]:
					frequencies.extend(charFreq.generateFrequencies(count))

	# Create spectral data object with characteristic frequencies
	spectralData = SpectralData()
	for freq in frequencies:
		spectralData.modes.append(HarmonicOscillator(frequency=freq))

	# Subtract out contributions to heat capacity from the characteristic modes
	import numpy
	Cp = [thermoData.getHeatCapacity(T) for T in thermo.ThermoGAData.CpTlist]
	Cv = numpy.array(Cp) / constants.R
	Tlist = numpy.array(thermo.ThermoGAData.CpTlist)
	for mode in spectralData.modes:
		Cv -= mode.heatCapacity(Tlist)
	# Subtract out translational modes
	Cv -= 1.5
	# Subtract out external rotational modes
	Cv -= (1.5 if not linear else 1.0)
	# Subtract out PV term (Cp -> Cv)
	Cv -= 1.0
	# Check that all Cv values are still positive (should we do this?)
	#for C in Cv:
	#	if C <= 0.0: raise Exception('Remaining heat capacity is negative.')
	
	# Fit remaining frequencies and hindered rotors to the heat capacity data
	import fit
	freeVibrations = numVibrations - len(frequencies)
	if freeVibrations > 0 and numRotors > 0:
		vib, hind = fit.fitspectraldata(Cv, Tlist, freeVibrations, numRotors)
		for v in vib:
			spectralData.modes.append(HarmonicOscillator(frequency=v))
		for v, b in hind:
			spectralData.modes.append(HinderedRotor(frequency=v, barrier=b))
	elif freeVibrations > 0 and numRotors == 0:
		vib = fit.fitspectraldatanorotors(Cv, Tlist, freeVibrations)
		for v in vib:
			spectralData.modes.append(HarmonicOscillator(frequency=v))
	elif freeVibrations == 0 and numRotors > 0:
		hind = fit.fitspectraldatanooscillators(Cv, Tlist, numRotors)
		for v, b in hind:
			spectralData.modes.append(HinderedRotor(frequency=v, barrier=b))

	return spectralData

################################################################################
