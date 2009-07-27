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
Contains classes describing chemical species and their thermodynamics.
"""

import quantities as pq
import logging
import math
import pybel

import constants
import data
import chem
import structure

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

	CpTlist = list(pq.Quantity([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K').simplified)
	for i in range(len(CpTlist)): CpTlist[i] = float(CpTlist[i])

	def __init__(self, H298=0.0, S298=0.0, Cp=None, comment=''):
		"""Initialize a set of group additivity thermodynamic data."""

		self.H298 = H298
		self.S298 = S298
		self.Cp = Cp or [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
		
		self.H298 = float(pq.Quantity(H298, 'kcal/mol').simplified)
		self.S298 = float(pq.Quantity(S298, 'cal/(mol*K)').simplified)
		self.Cp = list(pq.Quantity([Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500], 'cal/(mol*K)').simplified)
		for i in range(len(self.Cp)): self.Cp[i] = float(self.Cp[i])
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
		"""
		Generate an XML representation of the thermodynamic data using the
		:data:`xml.dom.minidom` package. The `dom` and `root` parameters
		represent the DOM and the element in the DOM used as the parent of
		the generated XML.
		"""

		thermo = dom.createElement('thermodynamics')
		thermo.setAttribute('comment', self.comment)
		root.appendChild(thermo)

		enthalpy = dom.createElement('enthalpyOfFormation')
		enthalpy.setAttribute('temperature', '298.0 K')
		thermo.appendChild(enthalpy)
		data.createXMLQuantity(dom, enthalpy, self.H298, 'J/mol')

		entropy = dom.createElement('entropyOfFormation')
		entropy.setAttribute('temperature', '298.0 K')
		thermo.appendChild(entropy)
		data.createXMLQuantity(dom, entropy, self.S298, 'J/(mol*K)')

		for i, Cp in enumerate(self.Cp):
			heatCapacity = dom.createElement('heatCapacity')
			heatCapacity.setAttribute('temperature', '%s K' % (self.CpTlist[i]) )
			thermo.appendChild(heatCapacity)
			data.createXMLQuantity(dom, heatCapacity, Cp, 'J/(mol*K)')

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

	def getHeatCapacity(self, T):
		"""
		Return the heat capacity at temperature `T`.
		"""
		if T < ThermoGAData.CpTlist[0] and T != 298.0:
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from group additivity.')
		elif T > ThermoGAData.CpTlist[-1]:
			return self.Cp[-1]
		else:
			for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
					ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
				if Tmin <= T and T <= Tmax:
					return (Cpmax - Cpmin) * ((T - Tmin) / (Tmax - Tmin)) + Cpmin

	def getEnthalpy(self, T):
		"""
		Return the enthalpy at temperature `T`.
		"""
		H = self.H298
		if T < ThermoGAData.CpTlist[0] and T != 298.0:
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

	def getEntropy(self, T):
		"""
		Return the entropy at temperature `T`.
		"""
		S = self.S298
		if T < ThermoGAData.CpTlist[0] and T != 298.0:
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

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy at temperature `T`.
		"""
		if T < ThermoGAData.CpTlist[0] and T != T != 298.0:
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from group additivity.')
		return self.getEnthalpy(T) - T * self.getEntropy(T)
###############
class ThermoNASAPolynomial:
	"""A single NASA polynomial for thermochemistry. 
	See __init__ docstring for more details"""

	def __init__(self, T_range=None, coeffs=None, comment=''):
		"""Initialize a NASA polynomial data.
		
		T_range is [Tmin, Tmax]
		coeffs is [a1, a2, a3, a4, a5, a6, a7]
		
 		The thermodynamic properties are given by 
		 Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
		 H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
		 S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
		which was copied from http://www.me.berkeley.edu/gri-mech/data/nasa_plnm.html
		"""
		self.T_range = T_range or [0.0]*2
		self.coeffs = coeffs or [0.0]*7
		self.comment = comment
	def isValidTemperature(self,temperature):
		"""Is the given temperature within the range of this polynomial?"""
		if temperature<self.T_range[0]: return False
		if temperature>self.T_range[1]: return False
		return True
		
	def getHeatCapacity(self, T):
		"""
		Return the heat capacity at temperature `T`.
		assumes T in [K] and returns Cp in [J/mol/K]
		"""
		c = self.coeffs
		T2 = T*T
		# Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
#		HeatCapacityOverR = c[0] + c[1]*T + c[2]*T*T + c[3]*T*T*T + c[4]*T*T*T*T
		HeatCapacityOverR = c[0] + T*(c[1] + T*(c[2] + T*(c[3] + c[4]*T)))
		Cp = HeatCapacityOverR * constants.R
		return Cp
		
	def getEnthalpy(self, T):
		"""
		Return the enthalpy at temperature `T`.
		"""
		c = self.coeffs
		T2 = T*T
		T4 = T2*T2
		# H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
		EnthalpyOverR = ( c[0]*T + c[1]*T2/2 + c[2]*T2*T/3 + c[3]*T4/4 +
						  c[4]*T4*T/5 + c[5] )
		H = EnthalpyOverR * constants.R
		return H
		
	def getEntropy(self, T):
		"""
		Return the entropy at temperature `T`.
		"""
		# S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
		c = self.coeffs
		T2 = T*T
		T4 = T2*T2
		EntropyOverR = ( c[0]*math.log(T) + c[1]*T + c[2]*T2/2 +
					c[3]*T2*T/3 + c[4]*T4/4 + c[6] )
		S = EntropyOverR * constants.R
		return S
		
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy at temperature `T`.
		"""
		return self.getEnthalpy(T) - T * self.getEntropy(T)

	def toXML(self, dom, root):
#   <polynomial>
#      <validRange>
#         <bound kind="lower" property="temperature" units="K">300.000</bound>
#         <bound kind="upper" property="temperature" units="K">1000</bound>
#      </validRange>
#      <coefficient id="1" label="a1">4.0733</coefficient>
#      <coefficient id="2" label="a2">0.011308</coefficient>
#      <coefficient id="3" label="a3">-1.6565e-005</coefficient>
#      <coefficient id="4" label="a4">1.1784e-008</coefficient>
#      <coefficient id="5" label="a5">-3.3006e-012</coefficient>
#      <coefficient id="6" label="a6">-19054.23</coefficient>
#      <coefficient id="7" label="a7">4.4083</coefficient>
#   </polynomial>
		pass
		
class ThermoNASAData:
	"""
	A set of thermodynamic parameters given by NASA polynomials.
	- `comment` = a string describing the source of the data
	- `NASA_polynomials` = a list of ThermoNASAPolynomial() objects
	"""
	def __init__(self, NASA_polynomials = None, comment=''):
		"""Initialize a set of NASA thermodynamic data."""
		self.NASA_polynomials = NASA_polynomials or [] # a list of ThermoNASAPolynomial objects
		self.comment = comment
		# do some tests?
	
	def addPolynomial(self, polynomial):
		if not isinstance(polynomial,ThermoNASAPolynomial):
			raise TypeError("polynomial should be instance of ThermoNASAPolynomial")
		self.NASA_polynomials.append(polynomial)

	def selectPolynomialForT(self, temperature):
		for poly in self.NASA_polynomials:
			if poly.isValidTemperature(temperature): break
		else: 
			raise data.TemperatureOutOfRangeException("No polynomial found for T=%s"%temperature)
		return poly

	def toXML(self, dom, root):
# <thermodynamicPolynomials type="nasa7">
#   <referenceState>
#      <Tref units="K">298.15</Tref>
#      <Pref units="Pa">100000</Pref>
#   </referenceState>
#   <dfH units="J/mol">-145177.5731</dfH>
#   <polynomial> FROM NASA_polynomials </polynomial>
# </thermodynamicPolynomials>
		pass

	def __str__(self):
		"""
		Return a string summarizing the thermodynamic data.
		"""
		string = 'NASA polynomials'
		return string
		
	def getHeatCapacity(self, T):
		"""
		Return the heat capacity at temperature `T`.
		"""
		poly = self.selectPolynomialForT(T)
		return poly.getHeatCapacity(T)
		
	def getEnthalpy(self, T):
		"""
		Return the enthalpy at temperature `T`.
		"""
		poly = self.selectPolynomialForT(T)
		return poly.getEnthalpy(T)
		
	def getEntropy(self, T):
		"""
		Return the entropy at temperature `T`.
		"""
		poly = self.selectPolynomialForT(T)
		return poly.getEntropy(T)

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy at temperature `T`.
		"""
		poly = self.selectPolynomialForT(T)
		return poly.getFreeEnergy(T)


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

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (path))

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

	def contains(self, structure):
		"""
		Search the dictionary for the specified `structure`. If found, the label
		corresponding to the structure is returned. If not found, :data:`None`
		is returned.
		"""
		for label, struct in self.dictionary.iteritems():
			match, map21, map12 = structure.isIsomorphic(struct)
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

		if structure.getRadicalCount() > 0:
			# Make a copy of the structure so we don't change the original
			struct = structure.copy()
			# Saturate structure by replacing all radicals with bonds to
			# hydrogen atoms
			added = {}
			for atom in struct.atoms():
				for i in range(0, atom.getFreeElectronCount()):
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

				thermoData += self.radicalDatabase.getThermoData(struct, {'*':atom})

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
					thermoData += self.groupDatabase.getThermoData(structure, {'*':atom})
					# Correct for gauche and 1,5- interactions
					thermoData += self.gaucheDatabase.getThermoData(structure, {'*':atom})
					thermoData += self.int15Database.getThermoData(structure, {'*':atom})
					thermoData += self.otherDatabase.getThermoData(structure, {'*':atom})

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

thermoDatabase = None
forbiddenStructures = None

################################################################################

def getThermoData(structure):
	"""
	Get the thermodynamic data associated with `structure` by looking in the
	loaded thermodynamic database.
	"""
	return thermoDatabase.getThermoData(structure)

################################################################################

class Species:
	"""
	Represent a chemical species (including all of its resonance forms). Each
	species has a unique integer `id` assigned automatically by RMG and a
	not-necessarily unique string `label`. The *structure* variable contains a
	list of :class:`Structure` objects representing each resonance form. The
	`reactive` flag is :data:`True` if the species can react and :data:`False`
	if it is inert.
	"""

	def __init__(self, id=0, label='', structure=None, reactive=True):
		"""
		Initialize a Species object.
		"""
		self.id = id
		
		self.label = label
		if structure is not None:
			structure.simplifyAtomTypes()
			self.structure = [structure]
		else:
			self.structure = []
		self.reactive = reactive

		self.thermoData = None
		self.lennardJones = None
		self.spectralData = None

	def __cmp__(self, other):
		"""
		A comparison function that can be used to sort lists of Species objects.
		Currently the sorting method is by increasing ID.
		"""
		return cmp(self.id, other.id)

	def __hash__(self):
		"""
		A hash function that allows for use in dictionaries et al. Currently the
		species ID is used.
		"""
		return self.id

	def getFormula(self):
		"""
		Return the chemical formula for the species.
		"""
		return self.structure[0].getFormula()

	def fromAdjacencyList(self, adjstr):
		"""
		Convert an adjacency list string `adjstr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromAdjacencyList(adjstr)

	def fromCML(self, cmlstr):
		"""
		Convert a string of CML `cmlstr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromCML(cmlstr)
	
	def fromInChI(self, inchistr):
		"""
		Convert an InChI string `inchistr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromInChI(inchistr)

	def fromSMILES(self, smilesstr):
		"""
		Convert a SMILES string `smilesstr` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromSMILES(smilesstr)

	def fromOBMol(self, obmol):
		"""
		Convert an OpenBabel OBMol object `obmol` to a Species object.
		"""
		self.structure = [structure.Structure()]
		self.structure[0].fromOBMol(obmol)

	def toCML(self):
		"""
		Convert a Species object to CML.
		"""
		return self.structure[0].toCML()

	def toInChI(self):
		"""
		Convert a Species object to an InChI string.
		"""
		return self.structure[0].toInChI()

	def toOBMol(self):
		"""
		Convert a Species object to an OBMol object.
		"""
		return self.structure[0].toOBMol()

	def toSMILES(self):
		"""
		Convert a Species object to an InChI string.
		"""
		return self.structure[0].toSMILES()

	def toAdjacencyList(self):
		"""
		Convert a Species object to an adjacency list.
		"""
		return str(self) + '\n' + self.structure[0].toAdjacencyList()

	def getResonanceIsomers(self):
		"""
		Generate all of the resonance isomers of this species. The isomers are
		stored as a list in the `structure` attribute. If the length of
		`structure` is already greater than one, it is assumed that all of the
		resonance isomers have already been generated.
		"""

		if len(self.structure) != 1:
			return

		structure = self.structure[0]

		isomers = [structure]

		# Radicals
		if structure.getRadicalCount() > 0:
			# Iterate over resonance isomers
			index = 0
			while index < len(isomers):
				isomer = isomers[index]

				newIsomers = isomer.getAdjacentResonanceIsomers()
				
				for newIsomer in newIsomers:

					# Append to isomer list if unique
					found = False
					for isom in isomers:
						ismatch, map21, map12 = isom.isIsomorphic(newIsomer)
						if ismatch: found = True
					if not found:
						isomers.append(newIsomer)

				# Move to next resonance isomer
				index += 1

		for isomer in isomers:
			isomer.updateAtomTypes()

		self.structure = isomers

	def getThermoData(self):
		"""
		Generate thermodynamic data for the species by use of the thermo
		database.
		"""

		thermoData = []
		for structure in self.structure:
			structure.updateAtomTypes()
			thermoData.append(getThermoData(structure))

		# If multiple resonance isomers are present, use the thermo data of
		# the most stable isomer (i.e. one with lowest enthalpy of formation)
		# as the thermo data of the species
		self.thermoData = thermoData[0]
		for tdata in thermoData[1:]:
			if tdata.H298 < self.thermoData.H298:
				self.thermoData = tdata

	def getHeatCapacity(self, T):
		"""
		Return the heat capacity of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getHeatCapacity(T)

	def getEnthalpy(self, T):
		"""
		Return the enthalpy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEnthalpy(T)

	def getEntropy(self, T):
		"""
		Return the entropy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getEntropy(T)

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy of the species at temperature `T`.
		"""
		if self.thermoData is None: self.getThermoData()
		return self.thermoData.getFreeEnergy(T)

	def isIsomorphic(self, other):
		"""
		Returns :data:`True` if the two species are isomorphic and data:`False`
		otherwise.
		"""
		if other.__class__ == Species:
			for struct1 in self.structure:
				for struct2 in other.structure:
					ismatch, map21, map12 = struct1.isIsomorphic(struct2)
					if ismatch:
						return True
		elif other.__class__ == structure.Structure:
			for struct1 in self.structure:
				ismatch, map21, map12 = struct1.isIsomorphic(other)
				if ismatch:
					return True
		return False
	
	def isSubgraphIsomorphic(self, other):
		"""
		Returns :data:`True` if the species is isomorphic with the other
		functional group and data:`False` otherwise.
		"""
		for struct1 in self.structure:
			ismatch, map21, map12 = struct1.isSubgraphIsomorphic(other)
			if ismatch:
				return True
		return False

	def findSubgraphIsomorphisms(self, other):
		"""
		Returns a list of the subgraph matches between the species and the
		functional group `other`.
		"""
		maps12 = []; maps21 = []
		for struct1 in self.structure:
			ismatch, map21, map12 = struct1.findSubgraphIsomorphisms(other)
			maps12.extend(map12)
			maps21.extend(map21)
		return (len(maps12) > 0), maps21, maps12

	def __str__(self):
		"""
		Return a string representation of the species, in the form 'label(id)'.
		"""
		return self.label + '(' + str(self.id) + ')'

################################################################################

# The global list of species created at any point during RMG execution
# The list is stored in reverse of the order in which the species are created;
# when searching the list, it is more likely to match a recently created species
# than an older species
speciesList = []

# A cache of recently visited species
speciesCache = []
speciesCacheMaxSize = 4

def makeNewSpecies(structure, label='', reactive=True):
	"""
	Attempt to make a new species based on a chemical `structure`, which is a
	:class:`Structure` object.

	The proposed species is checked against the list of existing species; if the
	species already exists, this function returns
	the existing species. If the species does not exist, a :class:`Species`
	object is created and returned after being appended to the global species
	list.
	"""

#	# Recalculate atom types for proposed structure (hopefully not necessary)
#	structure.simplifyAtomTypes()
#	structure.updateAtomTypes()

	# First check cache and return if species is found
	for i, spec in enumerate(speciesCache):
		if spec.isIsomorphic(structure):
			speciesCache.pop(i)
			speciesCache.insert(0, spec)
			return spec

	# Return an existing species if a match is found
	for spec in speciesList:
		if spec.isIsomorphic(structure):
			speciesCache.insert(0, spec)
			if len(speciesCache) > speciesCacheMaxSize: speciesCache.pop()
			return spec

	# Return None if the species has a forbidden structure
	if forbiddenStructures is not None:
		for lbl, struct in forbiddenStructures.iteritems():
			match, map21, map12 = structure.isSubgraphIsomorphic(struct)
			if match: return None

	# Otherwise make a new species
	if label == '':
		label = structure.getFormula()
		for atom in structure.atoms():
			if atom.hasFreeElectron(): label += 'J'
	spec = Species(len(speciesList)+1, label, structure, reactive)
	speciesList.insert(0, spec)
	
	spec.getResonanceIsomers()
	if thermoDatabase is not None:
		spec.getThermoData()

	# Draw species in core
	if constants.drawMolecules:
		mol = pybel.Molecule(spec.toOBMol())
		mol.draw(False, constants.outputDir + '/species/' + str(spec) + '.png')

	# Note in the log
	logging.debug('Created new species ' + str(spec) + ': ' + spec.toInChI())
	
	# Return the newly created species
	speciesCache.insert(0, spec)
	if len(speciesCache) > speciesCacheMaxSize: speciesCache.pop()
	return spec

################################################################################

if __name__ == '__main__':

	import os.path
	import main
	main.initializeLog(logging.DEBUG)

	datapath = '../data/RMG_database/'

	logging.debug('General database: ' + os.path.abspath(datapath))
	thermoDatabase = ThermoDatabaseSet()
	thermoDatabase.load(datapath)

	for label, struct in thermoDatabase.groupDatabase.dictionary.iteritems():
		label = label.replace('/', '_').replace('\\', '_')
		f = open('data/%s.txt' % (label), 'w')
		f.write(struct.toDOT('graphname'))
		f.close()

#	adjlist = \
#"""
#1 C 0 {2,S} {4,S} {5,S} {6,S}
#2 C 0 {1,S} {3,S} {7,S} {8,S}
#3 C 0 {2,S} {4,S} {9,S} {10,S}
#4 C 0 {3,S} {1,S} {11,S} {12,S}
#5 H 0 {1,S}
#6 H 0 {1,S}
#7 H 0 {2,S}
#8 H 0 {2,S}
#9 H 0 {3,S}
#10 H 0 {3,S}
#11 H 0 {4,S}
#12 H 0 {4,S}
#"""
#
#	structure = structure.Structure()
#	structure.fromAdjacencyList(adjlist)
#	structure.updateAtomTypes()
#
#	print structure.toInChI()
#
#	thermoData = getThermoData(structure)
#
#	T = 1000.0
#	print 'Heat capacity at ' + str(T) + ': ' + str(thermoData.getHeatCapacity(T))
#	print 'Enthalpy at ' + str(T) + ': ' + str(thermoData.getEnthalpy(T))
#	print 'Entropy at ' + str(T) + ': ' + str(thermoData.getEntropy(T))
#	print 'Free energy at ' + str(T) + ': ' + str(thermoData.getFreeEnergy(T))
#
