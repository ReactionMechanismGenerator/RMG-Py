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

import os.path

import rmg.log as logging
import rmg.data as data
import rmg.constants as constants

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

		# Convert data in library to SpectralData objects or lists of
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

				# First item is the number of vibrations
				numVibrations = int(items.pop(0))

				# Next are the vibrational frequencies in cm^-1
				vibrations = [float(items.pop(0)) for i in range(numVibrations)]
					
				# Next item is the number of hindered rotors
				numRotors = int(items.pop(0))
				
				# Next are the hindered rotor frequency-degeneracy pairs, both in cm^-1
				rotors = [[float(items.pop(0)), float(items.pop(0))] for i in range(numRotors)]

				# Create the SpectralData object
				spectralData = SpectralData()
				for freq in vibrations:
					spectralData.modes.append(HarmonicOscillator(frequency=freq))
				for freq, barr in rotors:
					spectralData.modes.append(HinderedRotor(frequency=freq, barrier=barr))

				self.library[label] = spectralData

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (dictstr))

	def contains(self, structure):
		"""
		Search the dictionary for the specified `structure`. If found, the label
		corresponding to the structure is returned. If not found, :data:`None`
		is returned.
		"""
		for label, struct in self.dictionary.iteritems():
			if structure.isIsomorphic(struct):
				return label
		return None

frequencyDatabase = None

################################################################################

class GroupFrequencyDatabase(data.Database):
	"""
	Represent an RMG group frequency database.
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

		# Convert data in library to CharacteristicFrequency objects or lists of
		# [link, comment] pairs
		for label, item in self.library.iteritems():

			if item is None:
				pass
			elif not item.__class__ is tuple:
				raise data.InvalidDatabaseException('Group frequencies library should be tuple at this point. Instead got %r'%data)
			else:
				index, item = item
				if not (item.__class__ is str or item.__class__ is unicode):
					raise data.InvalidDatabaseException('Group frequencies library data format is unrecognized.')

				items = item.split()

				# First item is the symmetry correction
				frequencies = [int(items.pop(0))]

				# Items should be a multiple of three (no comments allowed at the moment)
				if len(items) % 3 != 0:
					raise data.InvalidDatabaseException('Unexpected number of items encountered in group frequencies library.')

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

class FrequencyDatabaseSet:
	"""
	A set of frequency databases, consisting of a primary database of species
	whose frequencies are known and a functional group database for estimating
	frequencies using characteristic values.
	"""

	def __init__(self):
		self.groupDatabase = GroupFrequencyDatabase()
		self.primaryDatabase = FrequencyDatabase()

	def load(self, datapath):
		"""
		Load a set of thermodynamics group additivity databases from the general
		database specified at `datapath`.
		"""
		datapath = os.path.abspath(datapath)
		def DTLpaths(prefix):
			"""Get a tuple of Dictionary, Tree, and Library paths for a given prefix"""
			dict_path = os.path.join(datapath, prefix+'_Dictionary.txt')
			tree_path = os.path.join(datapath, prefix+'_Tree.txt')
			libr_path = os.path.join(datapath, prefix+'_Library.txt')
			return dict_path, tree_path, libr_path

		logging.info('Loading frequency databases from %s...' % datapath)

		logging.verbose('Loading functional group frequency database from %s...' % datapath)
		self.groupDatabase.load(*DTLpaths('Group')) # the '*' unpacks the tuple into three separate arguments

		logging.verbose('Loading primary frequency database from %s...' % datapath)
		self.primaryDatabase.load(os.path.join(datapath, 'Primary_Dictionary.txt'),
			'',
			os.path.join(datapath, 'Primary_Library.txt'))

	def getSpectralData(self, struct, thermoData):
		"""
		Generate the spectral data for a :class:`structure.Structure` object
		`struct` with corresponding :class:`ThermoModel` object `thermoData`.
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

		# First try primary database
		label = self.primaryDatabase.contains(struct)
		if label is not None:
			return self.primaryDatabase.library[label]

		# If not in primary database, estimate using group database

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
		for node, data in self.groupDatabase.library.iteritems():
			ismatch, map12List, map21List = struct.findSubgraphIsomorphisms(self.groupDatabase.dictionary[node])
			if ismatch:
				count = len(map12List)
			else:
				count = 0
			if count % data[0] != 0:
				raise Exception('Incorrect number of matches of node "%s" while estimating frequencies of %s; expected a multiple of %s, got %s.' % (node, struct, data[0], count))
			groupCount[node] = count / data[0]

		# For debugging, print a list of the groups found
	#	print 'Groups found:'
	#	for node, count in groupCount.iteritems():
	#		if count != 0: print '\t', node, count

		# Get characteristic frequencies
		frequencies = []
		for node, count in groupCount.iteritems():
			for charFreq in self.groupDatabase.library[node][1:]:
				frequencies.extend(charFreq.generateFrequencies(count))

		# Check that we have the right number of degrees of freedom specified
		if len(frequencies) > numVibrations:
			# We have too many vibrational modes
			difference = len(frequencies) - numVibrations
			# First try to remove hindered rotor modes until the proper number of modes remains
			if numRotors >= difference:
				numRotors -= difference
				numVibrations = len(frequencies)
				logging.warning('For %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i internal rotors to compensate.' % (struct, difference))
			# If that won't work, turn off functional groups until the problem is underspecified again
			else:
				groupsRemoved = 0
				freqsRemoved = 0
				freqCount = len(frequencies)
				while freqCount > numVibrations:
					minDegeneracy, minNode = min([(sum([charFreq.degeneracy for charFreq in self.groupDatabase.library[node][1:]]), node) for node in groupCount if groupCount[node] > 0])
					if groupCount[minNode] > 1:
						groupCount[minNode] -= 1
					else:
						del groupCount[minNode]
					groupsRemoved += 1
					freqsRemoved += minDegeneracy
					freqCount -= minDegeneracy
				# Log warning
				logging.warning('For %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i groups (%i frequencies) to compensate.' % (struct, groupsRemoved, freqsRemoved))
				# Regenerate characteristic frequencies
				frequencies = []
				for node, count in groupCount.iteritems():
					for charFreq in self.groupDatabase.library[node][1:]:
						frequencies.extend(charFreq.generateFrequencies(count))

		# Create spectral data object with characteristic frequencies
		spectralData = SpectralData()
		for freq in frequencies:
			spectralData.modes.append(HarmonicOscillator(frequency=freq))

		# Subtract out contributions to heat capacity from the characteristic modes
		import numpy
		Tlist = numpy.arange(300.0, 1501.0, 100.0, numpy.float64)
		Cp = [thermoData.getHeatCapacity(T) for T in Tlist]
		Cv = numpy.array(Cp) / constants.R
		for mode in spectralData.modes:
			Cv -= mode.getHeatCapacity(Tlist)
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
		try:
			vib, hind = fit.fitSpectralDataToHeatCapacity(struct, Tlist, Cv, numVibrations - len(frequencies), numRotors)
			for freq, degen in vib:
				for d in range(degen):
					spectralData.modes.append(HarmonicOscillator(frequency=freq, degeneracy=1))
			for freq, barr, degen in hind:
				for d in range(degen):
					spectralData.modes.append(HinderedRotor(frequency=freq, barrier=barr, degeneracy=1))

			# Check: Does the fitted data reproduce the Cv data?
			# We use root mean square error per data point as our basis for judging
			Cp_fit = spectralData.getHeatCapacity(Tlist) + (3.0 if not linear else 2.5)
			Cp_data = [thermoData.getHeatCapacity(T) / 8.314472 for T in Tlist]
			rms = 0.0
			for i in range(len(Tlist)):
				rms += (Cp_fit[i] - Cp_data[i]) * (Cp_fit[i] - Cp_data[i])
			rms /= len(Tlist)
			if rms > 1.0:
				logging.warning('Fitted spectral data may not reproduce heat capacity data to within tolerance. RMS/point = %s' % rms)
				#if rms > 3.0:
				#	raise fit.SpectralFitException('Fitted spectral data does not reproduce heat capacity data to within tolerance. RMS = %s\nModel heat capacity is %s\nData heat capacity is %s' % (rms, Cp_fit, Cp_data))

		except fit.SpectralFitException, e:
			e.msg += '\nThe species I was fitting spectral data to was %s' % str(struct)
			raise

		return spectralData

################################################################################

def loadFrequencyDatabase(dstr):
	"""
	Load the RMG thermo database located at `dstr` into the global variable
	`rmg.spectral.data.frequencyDatabase`.
	"""
	global frequencyDatabase

	frequenciesDatabasePath = os.path.join(dstr,'frequencies_groups')
	
	# Create and load frequency databases
	frequencyDatabase = FrequencyDatabaseSet()
	logging.debug('\tFrequencies database from '+frequenciesDatabasePath)
	frequencyDatabase.load(frequenciesDatabasePath)

################################################################################

def generateSpectralData(struct, thermoData):
	return frequencyDatabase.getSpectralData(struct, thermoData)

################################################################################
