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

import data

################################################################################

class ThermoData:
	"""
	A base class for all forms of thermodynamic data used by RMG. The common
	attributes are:

	========= ============================================================
	Attribute Meaning
	========= ============================================================
	`Trange`  a 2-tuple containin the minimum and maximum temperature in K
	`comment` a string describing the source of the data
	========= ============================================================
	"""

	def __init__(self, Trange=None, comment=''):
		self.Trange = Trange or (0.0, 0.0)
		self.comment = comment

	def isTemperatureValid(self, T):
		"""
		Return :data:`True` if the temperature `T` in K is within the valid
		temperature range of the thermodynamic data, or :data:`False` if not.
		"""
		if self.Trange is None:
			return True
		else:
			return self.Trange[0] <= T and T <= self.Trange[1]

################################################################################

ThermoGAData_CpTlist = list(pq.Quantity([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K').simplified)
ThermoGAData_CpTlist = [float(T) for T in ThermoGAData_CpTlist]

class ThermoGAData(ThermoData):
	"""
	A set of thermodynamic parameters as determined from Benson's group
	additivity data. The attributes are:

	========= ========================================================
	Attribute Meaning
	========= ========================================================
	`H298`    the standard enthalpy of formation at 298 K in J/mol
	`S298`    the standard entropy of formation at 298 K in J/mol*K
	`Cp`      the standard heat capacity in J/mol*K at 300, 400, 500,
		      600, 800, 1000, and 1500 K
	========= ========================================================
	"""

	def __init__(self, H298=0.0, S298=0.0, Cp=None, comment=''):
		"""Initialize a set of group additivity thermodynamic data."""
		ThermoData.__init__(self, Trange=(298.0, 2500.0), comment=comment)
		self.H298 = H298
		self.S298 = S298
		self.Cp = Cp or [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	def __add__(self, other):
		"""
		Add two sets of thermodynamic data together. All parameters are
		considered additive. Returns a new :class:`ThermoGAData` object that is
		the sum of the two sets of thermodynamic data.
		"""
		new = ThermoGAData()
		new.H298 = self.H298 + other.H298
		new.S298 = self.S298 + other.S298
		new.Cp = [self.Cp[i] + other.Cp[i] for i in range(len(self.Cp))]
		if self.comment == '': new.comment = other.comment
		elif other.comment == '': new.comment = self.comment
		else: new.comment = self.comment + '+ ' + other.comment
		return new

	def __str__(self):
		"""
		Return a string summarizing the thermodynamic data.
		"""
		string = ''
		string += 'Enthalpy of formation: %s J/mol\n' % (self.H298)
		string += 'Entropy of formation: %s J/mol*K\n' % (self.S298)
		for T, Cp in zip(ThermoGAData.Tlist, self.Cp):
			string += 'Heat capacity at %s K: %s J/mol*K' % (T, Cp)
		string += 'Comment: %s' % (self.comment)
		return string

	def getHeatCapacity(self, T):
		"""
		Return the heat capacity in J/mol*K at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from group additivity.')
		if T < 300.0:
			return self.Cp[0]
		elif T > ThermoGAData_CpTlist[-1]:
			return self.Cp[-1]
		else:
			for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData_CpTlist[:-1], \
					ThermoGAData_CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
				if Tmin <= T and T <= Tmax:
					return (Cpmax - Cpmin) * ((T - Tmin) / (Tmax - Tmin)) + Cpmin

	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		H = self.H298
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for enthalpy estimation from group additivity.')
		for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData_CpTlist[:-1], \
				ThermoGAData_CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
			if T > Tmin:
				slope = (Cpmax - Cpmin) / (Tmax - Tmin)
				intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
				if T < Tmax:	H += 0.5 * slope * (T**2 - Tmin**2) + intercept * (T - Tmin)
				else:			H += 0.5 * slope * (Tmax**2 - Tmin**2) + intercept * (Tmax - Tmin)
		if T > ThermoGAData_CpTlist[-1]:
			H += self.Cp[-1] * (T - ThermoGAData_CpTlist[-1])
		return H

	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		import math
		S = self.S298
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for entropy estimation from group additivity.')
		for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData_CpTlist[:-1], \
				ThermoGAData_CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
			if T > Tmin:
				slope = (Cpmax - Cpmin) / (Tmax - Tmin)
				intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
				if T < Tmax:	S += slope * (T - Tmin) + intercept * math.log(T/Tmin)
				else:			S += slope * (Tmax - Tmin) + intercept * math.log(Tmax/Tmin)
		if T > ThermoGAData_CpTlist[-1]:
			S += self.Cp[-1] * math.log(T / ThermoGAData.Tlist[-1])
		return S

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from group additivity.')
		return self.getEnthalpy(T) - T * self.getEntropy(T)

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
			heatCapacity.setAttribute('temperature', '%s K' % (ThermoGAData_CpTlist[i]) )
			thermo.appendChild(heatCapacity)
			data.createXMLQuantity(dom, heatCapacity, Cp, 'J/(mol*K)')

################################################################################

class ThermoNASAPolynomial(ThermoData):
	"""
	A single NASA polynomial for thermodynamic data. The `coeffs` attribute
	stores the seven polynomial coefficients
	:math:`\\mathbf{a} = \\left[a_1\\ a_2\\ a_3\\ a_4\\ a_5\\ a_6\\ a_7 \\right]`
	from which the relevant thermodynamic parameters are evaluated via the
	expressions

	.. math:: \\frac{C_\\mathrm{p}(T)}{R} = a_1 + a_2 T + a_3 T^2 + a_4 T^3 + a_5 T^4

	.. math:: \\frac{H(T)}{RT} = a_1 + \\frac{1}{2} a_2 T + \\frac{1}{3} a_3 T^2 + \\frac{1}{4} a_4 T^3 + \\frac{1}{5} a_5 T^4 + \\frac{a_6}{T}

	.. math:: \\frac{S(T)}{R} = a_1 \\ln T + a_2 T + \\frac{1}{2} a_3 T^2 + \\frac{1}{3} a_4 T^3 + \\frac{1}{4} a_5 T^4 + a_7

	The above was adapted from `this page <http://www.me.berkeley.edu/gri-mech/data/nasa_plnm.html>`_.
	"""

	def __init__(self, T_range=None, coeffs=None, comment=''):
		ThermoData.__init__(self, Trange=T_range, comment=comment)
		self.coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

	def getHeatCapacity(self, T):
		"""
		Return the heat capacity in J/mol*K at temperature `T` in K.
		"""
		import constants
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from NASA polynomial.')
		c = self.coeffs
		T2 = T*T
		# Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
#		HeatCapacityOverR = c[0] + c[1]*T + c[2]*T*T + c[3]*T*T*T + c[4]*T*T*T*T
		HeatCapacityOverR = c[0] + T*(c[1] + T*(c[2] + T*(c[3] + c[4]*T)))
		Cp = HeatCapacityOverR * constants.R
		return Cp

	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		import constants
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for enthalpy estimation from NASA polynomial.')
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
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		import math
		import constants
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for entropy estimation from NASA polynomial.')
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
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from NASA polynomial.')
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

################################################################################

class ThermoNASAData(ThermoData):
	"""
	A set of thermodynamic parameters given by NASA polynomials. This class
	stores a list of :class:`ThermoNASAPolynomial` objects in the `polynomials`
	attribute. When evaluating a thermodynamic quantity, a polynomial that
	contains the desired temperature within its valid range will be used.
	"""

	def __init__(self, polynomials=None, comment=''):
		ThermoData.__init__(self, Trange=None, comment=comment)
		self.polynomials = polynomials or []

	def addPolynomial(self, polynomial):
		if not isinstance(polynomial,ThermoNASAPolynomial):
			raise TypeError("Polynomial attribute should be instance of ThermoNASAPolynomial")
		self.polynomials.append(polynomial)

	def selectPolynomialForTemperature(self, T):
		for poly in self.polynomials:
			if poly.isTemperatureValid(T): break
		else:
			raise data.TemperatureOutOfRangeException("No polynomial found for T=%s" % T)
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
		Return the heat capacity in J/mol*K at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
		return poly.getHeatCapacity(T)

	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
		return poly.getEnthalpy(T)

	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
		return poly.getEntropy(T)

	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
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
			data = ThermoGAData()
		else:
			data = self.library[node]

		while data.__class__ != ThermoGAData and data is not None:
			if data[0].__class__ == str or data[0].__class__ == unicode:
				data = self.library[data[0]]

		# This code prints the hierarchy of the found node; useful for debugging
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
		import logging

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

	def getThermoData(self, struct):
		"""
		Determine the group additivity thermodynamic data for the structure
		`structure`.
		"""

		import chem
		import structure

		# First check to see if structure is in primary thermo library
		label = self.primaryDatabase.contains(struct)
		if label is not None:
			return self.primaryDatabase.library[label]

		thermoData = ThermoGAData()

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

thermoDatabase = None
forbiddenStructures = None

################################################################################

def getThermoData(struct):
	"""
	Get the thermodynamic data associated with `structure` by looking in the
	loaded thermodynamic database.
	"""
	import constants
	import math
		
	thermoData = thermoDatabase.getThermoData(struct)

	# Correct entropy for symmetry number
	struct.calculateSymmetryNumber()
	thermoData.S298 -= constants.R * math.log(struct.symmetryNumber)

	return thermoData

################################################################################

if __name__ == '__main__':

	pass