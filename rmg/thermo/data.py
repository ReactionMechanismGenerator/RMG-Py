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

import os
import math

import rmg.constants as constants
import rmg.data as data
import rmg.log as logging
import rmg.chem as chem
import rmg.structure as structure

from model import *
from converter import *

################################################################################

class ThermoDatabase(data.Database):
	"""
	Represent an RMG thermodynamics database.
	"""

	def __init__(self):
		"""Call the generic `data.Database.__init__()` method.

		This in turn creates
			* self.dictionary = Dictionary()
			* self.library = Library()
			* self.tree = Tree()
		"""
		data.Database.__init__(self)


	def load(self, dictstr, treestr, libstr):
		"""
		Load a thermodynamics group additivity database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. The tree is optional, and should
		be set to '' if not desired.
		"""

		# Load dictionary, library, and (optionally) tree
		data.Database.load(self, dictstr, treestr, libstr)

		# Convert data in library to ThermoGAModel objects or lists of
		# [link, comment] pairs
		for label, item in self.library.iteritems():

			if item is None:
				pass
			elif not item.__class__ is tuple:
				raise data.InvalidDatabaseException('Thermo library should be tuple at this point. Instead got %r'%data)
			else:
				index,item = item # break apart tuple, recover the 'index' - the beginning of the line in the library file.
				# Is't it dangerous having a local variable with the same name as a module?
				# what if we want to raise another data.InvalidDatabaseException() ?
				if not ( item.__class__ is str or item.__class__ is unicode) :
					raise data.InvalidDatabaseException('Thermo library data format is unrecognized.')

				items = item.split()
				try:
					thermoData = []; comment = ''
					# First 12 entries are thermo data
					for i in range(0, 12):
						thermoData.append(float(items[i]))
					# Remaining entries are comment
					for i in range(12, len(items)):
						comment += items[i] + ' '

					thermoGAData = ThermoGAModel()
					thermoGAData.fromDatabase(thermoData, comment)
					thermoGAData.index = index
					self.library[label] = thermoGAData
				except (ValueError, IndexError), e:
					# Split data into link string and comment string; store
					# as list of length 2
					link = items[0]
					comment = item[len(link)+1:].strip()
					self.library[label] = [link, comment]

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (dictstr))

		#self.library.removeLinks()

	def toXML(self):
		"""
		Return an XML representation of the thermo database.
		"""
		import xml.dom.minidom
		dom = xml.dom.minidom.parseString('<database type="thermodynamics"></database>')
		root = dom.documentElement

		data.Database.toXML(self, dom, root)

		return dom.toprettyxml()

	def getThermoData(self, structure, atom):
		"""
		Determine the group additivity thermodynamic data for the atom `atom`
		in the structure `structure`.
		"""

		node = self.descendTree(structure, atom, None)

		if node not in self.library:
			# No data present (e.g. bath gas)
			data = ThermoGAModel()
		else:
			data = self.library[node]

		while data.__class__ != ThermoGAModel and data is not None:
			if data[0].__class__ == str or data[0].__class__ == unicode:
				data = self.library[data[0]]

			# This code prints the hierarchy of the found node; useful for debugging
		#result = ''
		#while node is not None:
		#	result = ' -> ' + node + result
		#	node = self.tree.parent[node]
		#print result[4:]
	    #
		return data

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

################################################################################

class ThermoDatabaseSet:
	"""
	A set of thermodynamics group additivity databases, consisting of a primary
	database of functional groups and a number of secondary databases to provide
	corrections for 1,5-interactions, gauche interactions, radicals, rings,
	and other functionality. There is also a primary library containing data for
	individual species.
	"""


	def __init__(self):
		self.groupDatabase = ThermoDatabase()
		self.int15Database = ThermoDatabase()
		self.gaucheDatabase = ThermoDatabase()
		self.otherDatabase = ThermoDatabase()
		self.radicalDatabase = ThermoDatabase()
		self.ringDatabase = ThermoDatabase()
		self.primaryDatabase = ThermoDatabase()

	def load(self, datapath):
		"""
		Load a set of thermodynamics group additivity databases from the general
		database specified at `datapath`.
		"""
		datapath = os.path.join(datapath,'thermo_groups')
		datapath = os.path.abspath(datapath)
		def DTLpaths(prefix):
			"""Get a tuple of Dictionary, Tree, and Library paths for a given prefix"""
			dict_path = os.path.join(datapath, prefix+'_Dictionary.txt')
			tree_path = os.path.join(datapath, prefix+'_Tree.txt')
			libr_path = os.path.join(datapath, prefix+'_Library.txt')
			return dict_path, tree_path, libr_path

		logging.info('Loading thermodynamics databases from %s...' % datapath)

		logging.verbose('Loading functional group thermo database from %s...' % datapath)
		self.groupDatabase.load(*DTLpaths('Group')) # the '*' unpacks the tuple into three separate arguments

		logging.verbose('Loading 1,5 interactions thermo database from %s...' % datapath)
		self.int15Database.load(*DTLpaths('15'))

		logging.verbose('Loading gauche interactions thermo database from %s...' % datapath)
		self.gaucheDatabase.load(*DTLpaths('Gauche'))

		logging.verbose('Loading radical corrections thermo database from %s...' % datapath)
		self.radicalDatabase.load(*DTLpaths('Radical'))

		logging.verbose('Loading ring corrections thermo database from %s...' % datapath)
		self.ringDatabase.load(*DTLpaths('Ring'))

		logging.verbose('Loading other corrections thermo database from %s...' % datapath)
		self.otherDatabase.load(*DTLpaths('Other'))

		logging.verbose('Loading primary thermo database from %s...' % datapath)
		self.primaryDatabase.load(os.path.join(datapath, 'Primary_Dictionary.txt'), \
			'', \
			os.path.join(datapath, 'Primary_Library.txt'))

	def saveXML(self, datapath):
		"""
		Save the loaded databases in the set to XML files.
		"""

		f = open(datapath + 'thermo/group.xml', 'w')
		f.write(self.groupDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/1,5-interactions.xml', 'w')
		f.write(self.int15Database.toXML())
		f.close()

		f = open(datapath + 'thermo/gauche-interactions.xml', 'w')
		f.write(self.gaucheDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/other.xml', 'w')
		f.write(self.otherDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/radicals.xml', 'w')
		f.write(self.radicalDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/ring.xml', 'w')
		f.write(self.ringDatabase.toXML())
		f.close()

		f = open(datapath + 'thermo/primary.xml', 'w')
		f.write(self.primaryDatabase.toXML())
		f.close()

	def getThermoData(self, struct):
		"""
		Determine the group additivity thermodynamic data for the structure
		`structure`.
		"""

		# First check to see if structure is in primary thermo library
		label = self.primaryDatabase.contains(struct)
		if label is not None:
			return self.primaryDatabase.library[label]

		thermoData = ThermoGAModel()

		if struct.getRadicalCount() > 0:
			# Make a copy of the structure so we don't change the original
			saturatedStruct = struct.copy()
			# Saturate structure by replacing all radicals with bonds to
			# hydrogen atoms
			added = {}
			for atom in saturatedStruct.atoms():
				for i in range(0, atom.getFreeElectronCount()):
					H = chem.Atom('H', '0')
					bond = chem.Bond([atom, H], 'S')
					saturatedStruct.addAtom(H)
					saturatedStruct.addBond(bond)
					atom.decreaseFreeElectron()
					if atom not in added:
						added[atom] = []
					added[atom].append(bond)
			# Update the atom types of the saturated structure (not sure why
			# this is necessary, because saturating with H shouldn't be
			# changing atom types, but it doesn't hurt anything and is not
			# very expensive, so will do it anyway)
			saturatedStruct.simplifyAtomTypes()
			saturatedStruct.updateAtomTypes()
			# Get thermo estimate for saturated form of structure
			thermoData = self.getThermoData(saturatedStruct)
			# For each radical site, get radical correction
			# Only one radical site should be considered at a time; all others
			# should be saturated with hydrogen atoms
			for atom in added:

				# Remove the added hydrogen atoms and bond and restore the radical
				for bond in added[atom]:
					H = bond.atoms[1]
					saturatedStruct.removeBond(bond)
					saturatedStruct.removeAtom(H)
					atom.increaseFreeElectron()

				thermoData += self.radicalDatabase.getThermoData(saturatedStruct, {'*':atom})

				# Re-saturate
				for bond in added[atom]:
					H = bond.atoms[1]
					saturatedStruct.addAtom(H)
					saturatedStruct.addBond(bond)
					atom.decreaseFreeElectron()

			# Subtract the enthalpy of the added hydrogens
			thermoData_H = self.primaryDatabase.library['H']
			for bond in added[atom]:
				thermoData.H298 -= thermoData_H.H298
				#thermoData.S298 -= thermoData_H.S298

			# Correct the entropy for the symmetry number

		else:
			# Generate estimate of thermodynamics
			for atom in struct.atoms():
				# Iterate over heavy (non-hydrogen) atoms
				if atom.isNonHydrogen():
					# Get initial thermo estimate from main group database
					thermoData += self.groupDatabase.getThermoData(struct, {'*':atom})
					# Correct for gauche and 1,5- interactions
					thermoData += self.gaucheDatabase.getThermoData(struct, {'*':atom})
					thermoData += self.int15Database.getThermoData(struct, {'*':atom})
					thermoData += self.otherDatabase.getThermoData(struct, {'*':atom})

			# Do ring corrections separately because we only want to match
			# each ring one time; this doesn't work yet
			rings = struct.getSmallestSetOfSmallestRings()
			for ring in rings:

				# Make a temporary structure containing only the atoms in the ring
				ringStructure = structure.Structure()
				for atom in ring: ringStructure.addAtom(atom)
				for atom1 in ring:
					for atom2 in ring:
						if struct.hasBond(atom1, atom2):
							ringStructure.addBond(struct.getBond(atom1, atom2))

				# Get thermo correction for this ring
				thermoData += self.ringDatabase.getThermoData(ringStructure, {})

		return thermoData

################################################################################

thermoDatabase = None

forbiddenStructures = None

def loadThermoDatabase(dstr):
	"""
	Load the RMG thermo database located at `dstr` into the global variable
	`rmg.species.thermoDatabase`. Also loads the forbidden structures into
	`rmg.thermo.forbiddenStructures`.
	"""
	global thermoDatabase
	global forbiddenStructures
	
	thermoDatabase = ThermoDatabaseSet()
	thermoDatabase.load(dstr)

	forbiddenStructures = data.Dictionary()
	forbiddenStructures.load(os.path.join(dstr, 'ForbiddenStructures.txt'))
	forbiddenStructures.toStructure()

################################################################################

def generateThermoData(struct, thermoClass=NASAModel):
	"""
	Get the thermodynamic data associated with `structure` by looking in the
	loaded thermodynamic database. The parameter `thermoClass` is the class of
	thermo object you want returning; default is :class:`NASAModel`.
	"""

	GAthermoData = thermoDatabase.getThermoData(struct)

	# Correct entropy for symmetry number
	struct.calculateSymmetryNumber()
	GAthermoData.S298 -= constants.R * math.log(struct.symmetryNumber)

	logging.debug('Group-additivity thermo data: %s' % GAthermoData)

	if thermoClass == ThermoGAModel:
		return GAthermoData  # return here because Wilhoit conversion not wanted

	# Convert to Wilhoit
	rotors = struct.calculateNumberOfRotors()
	atoms = len(struct.atoms())
	linear = struct.isLinear()
	WilhoitData = convertGAtoWilhoit(GAthermoData,atoms,rotors,linear)

	logging.debug('Wilhoit thermo data: %s' % WilhoitData)

	if thermoClass == WilhoitModel:
		return WilhoitData

	# Convert to NASA
	NASAthermoData = convertWilhoitToNASA(WilhoitData)

	logging.debug('NASA thermo data: %s' % NASAthermoData)

	# compute the error for the entire conversion, printing it as info or warning (if it is sufficiently high)
	rmsErr = NASAthermoData.rmsErr(GAthermoData)
	if(rmsErr > 0.35):
	    logging.warning("Poor overall GA-to-NASA fit: Overall RMS error in heat capacity fit = %.3f*R." % (rmsErr))
	else:
	    logging.debug("Overall RMS error in heat capacity fit = %.3f*R" % (rmsErr))

	if thermoClass == NASAModel:
		return NASAthermoData

	# Still not returned?
	raise Exception("Cannot convert themo data into class %r"%(required_class))

################################################################################

def isStructureForbidden(struct):
	"""
	Return :data:`True` if the structure `struct` contains any of the subgraphs
	listed in the forbidden structures database, or :data:`False` if not.
	"""

	if forbiddenStructures:
		for lbl, s in forbiddenStructures.iteritems():
			if struct.isSubgraphIsomorphic(s):
				return True

	return False
