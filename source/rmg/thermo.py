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
A module for working with the thermodynamics of chemical species. This module
seeks to provide functionality for answering the question, "Given a species,
what are its thermodynamics?"
"""

import quantities as pq
import logging
import math

import data
import chem

################################################################################

class ThermoGAData:
	"""
	A set of thermodynamic parameters as determined from Benson's group
	additivity data. The attributes are:

	- `H298` = the standard enthalpy of formation at 298 K in J/mol

	- `S298` = the standard entropy of formation at 298 K in J/mol*K

	- `Cp300` = the standard heat capacity at 300 K in J/mol*K

	- `Cp400` = the standard heat capacity at 400 K in J/mol*K

	- `Cp500` = the standard heat capacity at 500 K in J/mol*K

	- `Cp600` = the standard heat capacity at 600 K in J/mol*K

	- `Cp800` = the standard heat capacity at 800 K in J/mol*K

	- `Cp1000` = the standard heat capacity at 1000 K in J/mol*K

	- `Cp1500` = the standard heat capacity at 1500 K in J/mol*K

	- `comment` = a string describing the source of the data
	"""

	CpTlist = pq.Quantity([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K')

	def __init__(self, H298=pq.UncertainQuantity(0.0, 'J/mol', 0.0), \
			S298=pq.UncertainQuantity(0.0, 'J/(mol*K)', 0.0),
			Cp=pq.UncertainQuantity([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'J/(mol*K)', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), \
			comment=''):
		"""Initialize a set of group additivity thermodynamic data."""

		self.H298 = H298
		self.S298 = S298
		self.Cp = Cp
		self.comment = comment

	def fromDatabase(self, data, comment):
		"""
		Process a list of numbers `data` and associated description `comment`
		generated while reading from a thermodynamic database.
		"""

		if len(data) != 12:
			raise Exception('Invalid list of thermo data; should be a list of numbers of length 12.')

		H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, \
		  dH, dS, dCp = data
		self.H298 = pq.UncertainQuantity(H298, 'kcal/mol', dH)
		self.H298.units = 'J/mol'
		self.S298 = pq.UncertainQuantity(S298, 'cal/(mol*K)', dS)
		self.S298.units = 'J/(mol*K)'
		self.Cp = pq.UncertainQuantity([Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500], 'cal/(mol*K)', [dCp, dCp, dCp, dCp, dCp, dCp, dCp])
		self.Cp.units = 'J/(mol*K)'

		self.comment = comment

	def __add__(self, other):
		"""
		Add two sets of thermodynamic data together. All parameters are
		additive.
		"""
		if other is None: other = ThermoGAData()
		
		new = ThermoGAData()
		new.H298 = self.H298 + other.H298
		new.S298 = self.S298 + other.S298
		new.Cp = []
		for i in range(len(self.Cp)):
			new.Cp.append(self.Cp[i] + other.Cp[i])
		if self.comment == '': new.comment = other.comment
		elif other.comment == '': new.comment = self.comment
		else: new.comment = self.comment + '+ ' + other.comment
		return new

	def toXML(self, dom, root):

		thermo = dom.createElement('thermo')
		root.appendChild(thermo)

		self.valueToXML(dom, thermo, 'enthalpyOfFormation', self.H298, '298 K')
		self.valueToXML(dom, thermo, 'entropyOfFormation', self.S298, '298 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[0], '300 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[1], '400 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[2], '500 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[3], '600 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[4], '800 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[5], '1000 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[6], '1500 K')

		element = dom.createElement('comment')
		thermo.appendChild(element)
		comment = dom.createTextNode(self.comment)
		element.appendChild(comment)

	def valueToXML(self, dom, root, name, value, temp):
		element = dom.createElement(name)
		root.appendChild(element)

		units = str(value.units).split()[1]

		element.setAttribute('temperature', temp)
		element.setAttribute('units', units)
		element.setAttribute('uncertainty', str(value.uncertainty))

		valueNode = dom.createTextNode(str(float(value)))
		element.appendChild(valueNode)

	def __str__(self):
		"""
		Return a string summarizing the thermodynamic data.
		"""
		string = ''
		string += 'Enthalpy of formation: ' + str(self.H298) + '\n'
		string += 'Entropy of formation: ' + str(self.S298) + '\n'
		for i in range (0, len(self.Cp)):
			string += 'Heat capacity at ' + str(ThermoGAData.CpTlist[i]) + ': ' + str(self.Cp[i]) + '\n'
		string += 'Comment: ' + str(self.comment)
		return string

	def heatCapacity(self, T):
		"""
		Return the heat capacity at temperature `T`.
		"""
		if T < ThermoGAData.CpTlist[0]:
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from group additivity.')
		elif T > ThermoGAData.CpTlist[-1]:
			return self.Cp[-1]
		else:
			for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
					ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
				if Tmin <= T and T < Tmax:
					return (Cpmax - Cpmin) * ((T - Tmin) / (Tmax - Tmin)) + Cpmin

	def enthalpy(self, T):
		"""
		Return the enthalpy at temperature `T`.
		"""
		H = self.H298
		if T < ThermoGAData.CpTlist[0]:
			raise data.TemperatureOutOfRangeException('Invalid temperature for enthalpy estimation from group additivity.')
		for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
				ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
			if T > Tmin:
				slope = (Cpmax - Cpmin) / (Tmax - Tmin)
				intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
				if T < Tmax:	H += 0.5 * slope * (T**2 - Tmin**2) + intercept * (T - Tmin)
				else:			H += 0.5 * slope * (Tmax**2 - Tmin**2) + intercept * (Tmax - Tmin)
		if T > ThermoGAData.CpTlist[-1]:
			H += self.Cp[-1] * (T - ThermoGAData.CpTlist[-1])
		return H

	def entropy(self, T):
		"""
		Return the entropy at temperature `T`.
		"""
		S = self.S298
		if T < ThermoGAData.CpTlist[0]:
			raise data.TemperatureOutOfRangeException('Invalid temperature for entropy estimation from group additivity.')
		for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
				ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
			if T > Tmin:
				slope = (Cpmax - Cpmin) / (Tmax - Tmin)
				intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
				if T < Tmax:	S += slope * (T - Tmin) + intercept * math.log(T/Tmin)
				else:			S += slope * (Tmax - Tmin) + intercept * math.log(Tmax/Tmin)
		if T > ThermoGAData.CpTlist[-1]:
			S += self.Cp[-1] * math.log(T / ThermoGAData.CpTlist[-1])
		return S

	def freeEnergy(self, T):
		"""
		Return the Gibbs free energy at temperature `T`.
		"""
		if T < ThermoGAData.CpTlist[0]:
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from group additivity.')
		return self.enthalpy(T) - T * self.entropy(T)


################################################################################

class ThermoDatabase(data.Database):
	"""
	Represent an RMG thermodynamics database.
	"""

	def __init__(self):
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

		# Convert data in library to ThermoData objects or lists of
		# [link, comment] pairs
		for label, item in self.library.iteritems():

			if item is None:
				pass
			elif item.__class__ == str or item.__class__ == unicode:
				items = item.split()
				try:
					thermoData = []; comment = ''
					# First 12 entries are thermo data
					for i in range(0, 12):
						thermoData.append(float(items[i]))
					# Remaining entries are comment
					for i in range(12, len(items)):
						comment += items[i] + ' '

					thermoGAData = ThermoGAData()
					thermoGAData.fromDatabase(thermoData, comment)
					self.library[label] = thermoGAData
				except (ValueError, IndexError), e:
					# Split data into link string and comment string; store
					# as list of length 2
					link = items[0]
					comment = item[len(link)+1:].strip()
					self.library[label] = [link, comment]

			else:
				raise data.InvalidDatabaseException('Thermo library data format is unrecognized.')

		#self.library.removeLinks()

	def toXML(self):
		"""
		Return an XML representation of the thermo database.
		"""

		dom = xml.dom.minidom.parseString('<database type="thermodynamics"></database>')
		root = dom.documentElement

		Database.toXML(self, dom, root)

		return dom.toprettyxml()

	def getThermoData(self, structure, atom):
		"""
		Determine the group additivity thermodynamic data for the atom `atom`
		in the structure `structure`.
		"""

		node = self.descendTree(structure, atom, None)
		#print node
		
		if node not in self.library:
			# No data present (e.g. bath gas)
			data = ThermoGAData()
		else:
			data = self.library[node]

		while data.__class__ != ThermoGAData and data is not None:
			if data[0].__class__ == str or data[0].__class__ == unicode:
				data = self.library[data[0]]

#		result = ''
#		while node is not None:
#			result = ' -> ' + node + result
#			node = self.tree.parent[node]
#		print result[4:]

		return data

	def matchNodeToStructure(self, node, structure, atom):
		"""
		Return :data:`True` if the `structure` centered at `atom` matches the
		structure at `node` in the dictionary.
		"""
		group = self.dictionary[node]
		if group.__class__ == str or group.__class__ == unicode:
			if group.lower() == 'union':
				match = False
				for child in self.tree.children[node]:
					if self.matchNodeToStructure(child, structure, atom):
						match = True
				return match
			else:
				return False
		else:
			center = group.getCenter()
			if atom.equivalent(center):
				match, map12, map21 = structure.isSubgraphIsomorphic(group, {center:atom}, {atom:center})
				return match
			else:
				return False

	def descendTree(self, structure, atom, root=None):
		"""
		Descend the tree in search of the functional group node that best
		matches the local structure around `atom` in `structure`.
		"""

		if root is None:
			root = self.tree.top[0]
			if not self.matchNodeToStructure(root, structure, atom):
				return None

		next = None
		for child in self.tree.children[root]:
			if self.matchNodeToStructure(child, structure, atom):
				next = child

		if next is not None:
			return self.descendTree(structure, atom, next)
		else:
			return root

	def contains(self, structure):
		"""
		Search the dictionary for the specified `structure`. If found, the label
		corresponding to the structure is returned. If not found, :data:`None`
		is returned.
		"""
		for label, struct in self.dictionary.iteritems():
			match, map12, map21 = structure.isIsomorphic(struct)
			if match:
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

		logging.debug('\tThermodynamics databases:')

		self.groupDatabase.load(datapath + 'thermo/Group_Dictionary.txt', \
			datapath + 'thermo/Group_Tree.txt', \
			datapath + 'thermo/Group_Library.txt')
		logging.debug('\t\tFunctional groups')

		self.int15Database.load(datapath + 'thermo/15_Dictionary.txt', \
			datapath + 'thermo/15_Tree.txt', \
			datapath + 'thermo/15_Library.txt')
		logging.debug('\t\t1,5 interactions')

		self.gaucheDatabase.load(datapath + 'thermo/Gauche_Dictionary.txt', \
			datapath + 'thermo/Gauche_Tree.txt', \
			datapath + 'thermo/Gauche_Library.txt')
		logging.debug('\t\tGauche interactions')

		self.otherDatabase.load(datapath + 'thermo/Other_Dictionary.txt', \
			datapath + 'thermo/Other_Tree.txt', \
			datapath + 'thermo/Other_Library.txt')
		logging.debug('\t\tOther corrections')

		self.radicalDatabase.load(datapath + 'thermo/Radical_Dictionary.txt', \
			datapath + 'thermo/Radical_Tree.txt', \
			datapath + 'thermo/Radical_Library.txt')
		logging.debug('\t\tRadical corrections')

		self.ringDatabase.load(datapath + 'thermo/Ring_Dictionary.txt', \
			datapath + 'thermo/Ring_Tree.txt', \
			datapath + 'thermo/Ring_Library.txt')
		logging.debug('\t\tRing corrections')

		self.primaryDatabase.load(datapath + 'thermo/Primary_Dictionary.txt', \
			'', \
			datapath + 'thermo/Primary_Library.txt')
		logging.debug('\t\tPrimary thermo database')

	def saveXML(self, datapath):

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

	def getThermoData(self, structure):
		"""
		Determine the group additivity thermodynamic data for the structure
		`structure`.
		"""

		# First check to see if structure is in primary thermo library
		label = self.primaryDatabase.contains(structure)
		if label is not None:
			return self.primaryDatabase.library[label]

		thermoData = ThermoGAData()

		if structure.radicalCount() > 0:
			# Make a copy of the structure so we don't change the original
			struct = structure.copy()
			# Saturate structure by replacing all radicals with bonds to 
			# hydrogen atoms
			added = {}
			for atom in struct.atoms():
				for i in range(0, atom.freeElectronCount()):
					H = chem.Atom('H', '0')
					bond = chem.Bond([atom, H], 'S')
					struct.addAtom(H)
					struct.addBond(bond)
					atom.decreaseFreeElectron()
					if atom not in added:
						added[atom] = []
					added[atom].append(bond)
			# Get thermo estimate for saturated form of structure
			thermoData = self.getThermoData(struct)
			# For each radical site, get radical correction
			# Only one radical site should be considered at a time; all others
			# should be saturated with hydrogen atoms
			for atom in added:

				# Remove the added hydrogen atoms and bond and restore the radical
				for bond in added[atom]:
					H = bond.atoms[1]
					struct.removeBond(bond)
					struct.removeAtom(H)
					atom.increaseFreeElectron()

				thermoData += self.radicalDatabase.getThermoData(structure, atom)
				
				# Re-saturate
				for bond in added[atom]:
					H = bond.atoms[1]
					struct.addAtom(H)
					struct.addBond(bond)
					atom.decreaseFreeElectron()

			# Subtract the enthalpy of the added hydrogens
			
			# Correct the entropy for the symmetry number

		else:
			# Generate estimate of thermodynamics
			for atom in structure.atoms():
				# Iterate over heavy (non-hydrogen) atoms
				if atom.isNonHydrogen():
					# Get initial thermo estimate from main group database
					thermoData += self.groupDatabase.getThermoData(structure, atom)
					# Correct for gauche and 1,5- interactions
					thermoData += self.gaucheDatabase.getThermoData(structure, atom)
					thermoData += self.int15Database.getThermoData(structure, atom)
					thermoData += self.otherDatabase.getThermoData(structure, atom)

			# Do ring corrections separately because we only want to match
			# each ring one time; this doesn't work yet
#			atoms = structure.atoms()[:]
#			for atom in atoms:
#				# Iterate over heavy (non-hydrogen) atoms
#				if atom.isNonHydrogen():
#					newData = self.ringDatabase.getThermoData(structure, atom)
#					if newData is not None:
#						thermoData += nestructure.atoms()wData
#						atoms.remove(atom)

		return thermoData

database = None

################################################################################

def getThermoData(structure):
	"""
	Get the thermodynamic data associated with `structure` by looking in the
	loaded thermodynamic database.
	"""
	return database.getThermoData(structure)

################################################################################

if __name__ == '__main__':

	import os.path
	import main
	main.initializeLog(logging.DEBUG)

	datapath = '../data/RMG_database/'

	logging.debug('General database: ' + os.path.abspath(datapath))
	database = ThermoDatabaseSet()
	database.load(datapath)

#	adjlist = \
#"""
#1 C 0 {2,D} {7,S} {8,S}
#2 C 0 {1,D} {3,S} {9,S}
#3 C 0 {2,S} {4,D} {10,S}
#4 C 0 {3,D} {5,S} {11,S}
#5 C 0 {4,S} {6,S} {12,S} {13,S}
#6 C 0 {5,S} {14,S} {15,S} {16,S}
#7 H 0 {1,S}
#8 H 0 {1,S}
#9 H 0 {2,S}
#10 H 0 {3,S}
#11 H 0 {4,S}
#12 H 0 {5,S}
#13 H 0 {5,S}
#14 H 0 {6,S}
#15 H 0 {6,S}
#16 H 0 {6,S}
#"""

#	adjlist = \
#"""
#1 H 0 {2,S}
#2 H 0 {1,S}
#"""


	adjlist = \
"""
1 C 0 {2,S} {4,S} {5,S} {6,S}
2 C 0 {1,S} {3,S} {7,S} {8,S}
3 C 0 {2,S} {4,S} {9,S} {10,S}
4 C 0 {3,S} {1,S} {11,S} {12,S}
5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {2,S}
8 H 0 {2,S}
9 H 0 {3,S}
10 H 0 {3,S}
11 H 0 {4,S}
12 H 0 {4,S}
"""

	structure = chem.Structure()
	structure.fromAdjacencyList(adjlist)
	structure.updateAtomTypes()

	print structure.toInChI()

	thermoData = getThermoData(structure)

	T = pq.Quantity(1000.0, 'K')
	print 'Heat capacity at ' + str(T) + ': ' + str(thermoData.heatCapacity(T))
	print 'Enthalpy at ' + str(T) + ': ' + str(thermoData.enthalpy(T))
	print 'Entropy at ' + str(T) + ': ' + str(thermoData.entropy(T))
	print 'Free energy at ' + str(T) + ': ' + str(thermoData.freeEnergy(T))


