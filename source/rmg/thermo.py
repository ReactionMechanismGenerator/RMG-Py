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
Contains classes relating to thermochemistry.
"""

import quantities as pq
import logging

import data

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
		
	def __init__(self, H298=0.0, S298=0.0, Cp300=0.0, Cp400=0.0, Cp500=0.0, \
	             Cp600=0.0, Cp800=0.0, Cp1000=0.0, Cp1500=0.0, comment=''):
		"""Initialize a set of group additivity thermodynamic data."""
		
		self.H298 = H298
		self.S298 = S298
		self.Cp = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
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
		#self.H298.units = 'J/mol'
		self.S298 = pq.UncertainQuantity(S298, 'cal/(mol*K)', dS)
		#self.S298.units = 'J/(mol*K)'
		self.Cp = pq.UncertainQuantity([Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500], 'kcal/(mol*K)', [dCp, dCp, dCp, dCp, dCp, dCp, dCp])
		#self.Cp.units = 'J/(mol*K)'
		
		#self.Cp[0] = pq.UncertainQuantity(Cp300, 'kcal/(mol*K)', dCp)
		#self.Cp[0].units = 'J/(mol*K)'
		#self.Cp[1] = pq.UncertainQuantity(Cp400, 'kcal/(mol*K)', dCp)
		#self.Cp[1].units = 'J/(mol*K)'
		#self.Cp[2] = pq.UncertainQuantity(Cp500, 'kcal/(mol*K)', dCp)
		#self.Cp[2].units = 'J/(mol*K)'
		#self.Cp[3] = pq.UncertainQuantity(Cp600, 'kcal/(mol*K)', dCp)
		#self.Cp[3].units = 'J/(mol*K)'
		#self.Cp[4] = pq.UncertainQuantity(Cp800, 'kcal/(mol*K)', dCp)
		#self.Cp[4].units = 'J/(mol*K)'
		#self.Cp[5] = pq.UncertainQuantity(Cp1000, 'kcal/(mol*K)', dCp)
		#self.Cp[5].units = 'J/(mol*K)'
		#self.Cp[6] = pq.UncertainQuantity(Cp1500, 'kcal/(mol*K)', dCp)
		#self.Cp[6].units = 'J/(mol*K)'
		self.comment = comment
	
	def __add__(self, other):
		"""
		Add two sets of thermodynamic data together. All parameters are 
		additive.
		"""
		new = ThermoGAData()
		new.H298 = self.H298 + other.H298
		new.S298 = self.S298 + other.S298
		new.Cp = []
		for i in range(len(self.Cp)):
			new.Cp.append(self.Cp + other.Cp)
		new.comment = self.comment + '; ' + other.comment
		return new
	
	def heatCapacity(self, T):
		"""
		Return the heat capacity at the specified temperature `T`. This is done
		via linear interpolation between the provided values.
		"""
		T.units = 'K'; T = float(T)
		if T < 300.0:
			raise TemperatureOutOfRangeException('No thermodynamic data available for T < 300 K.')
		# Use Cp(1500 K) if T > 1500 K
		elif T > 1500.0: T = 1500.0
		
		Cpfun = scipy.interpolate.interp1d(ThermoGAData.CpTlist, self.Cp)
		return pq.Quantity(Cpfun(T), 'J/(mol*K)')
		
	def enthalpy(self, T):
		"""
		Return the enthalpy of formation at the specified temperature `T`.
		"""
		T.units = 'K'; T = float(T)
		if T < 300.0:
			raise TemperatureOutOfRangeException('No thermodynamic data available for T < 300 K.')
	
	def getCpLinearization(self):
		slope = []; intercept = []
		for i in range(0, len(self.Cp)-1):
			slope.append((self.Cp[i+1] - self.Cp[i]) / (ThermoGAData.CpTlist[i+1] - ThermoGAData.CpTlist[i]))
			print slope, ThermoGAData.CpTlist[i], slope * ThermoGAData.CpTlist[i]
			quit()
			intercept.append(self.Cp[i] - slope[i] * ThermoGAData.CpTlist[i])
		return slope, intercept
	
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
		
################################################################################

class ThermoDatabase(data.Database):
	"""
	Represent an RMG thermodynamics database. 
	"""

	def __init__(self):
		data.Database.__init__(self)
		
	def load(self, dictstr, treestr, libstr, isfgroup=True):
		"""
		Load a thermodynamics group additivity database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. The tree is optional, and should
		be set to '' if not desired.
		"""
		
		# Load dictionary, library, and (optionally) tree
		data.Database.load(self, dictstr, treestr, libstr, isfgroup)
		
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
				raise InvalidDatabaseException('Thermo library data format is unrecognized.')

		#self.library.removeLinks()
	
	def toXML(self):
		"""
		Return an XML representation of the thermo database.
		"""
		
		dom = xml.dom.minidom.parseString('<database type="thermodynamics"></database>')
		root = dom.documentElement
		
		Database.toXML(self, dom, root)
	
		return dom.toprettyxml()
	
################################################################################

def loadOldThermoDatabases(datapath):
	"""
	Load a set of thermodynamics group additivity databases from the general
	database specified at `datapath`. The loaded databases are: 
	
	- 1,5-interactions
	
	- Gauche interactions
	
	- Acyclic functional groups
	
	- Other functional groups
	
	- Radical functional groups
	
	- Cyclic functional groups
	
	- Primary thermo library
	"""

	logging.debug('\tThermodynamics databases:')
	
	int15Database = ThermoDatabase()
	int15Database.load(datapath + 'thermo/15_Dictionary.txt', \
		datapath + 'thermo/15_Tree.txt', \
		datapath + 'thermo/15_Library.txt', True)
	
	logging.debug('\t\t1,5 interactions')
	
	gaucheDatabase = ThermoDatabase()
	gaucheDatabase.load(datapath + 'thermo/Gauche_Dictionary.txt', \
		datapath + 'thermo/Gauche_Tree.txt', \
		datapath + 'thermo/Gauche_Library.txt', True)
	
	logging.debug('\t\tGauche interactions')
	
	groupDatabase = ThermoDatabase()
	groupDatabase.load(datapath + 'thermo/Group_Dictionary.txt', \
		datapath + 'thermo/Group_Tree.txt', \
		datapath + 'thermo/Group_Library.txt', True)
	
	logging.debug('\t\tAcyclic functional groups')
	
	otherDatabase = ThermoDatabase()
	otherDatabase.load(datapath + 'thermo/Other_Dictionary.txt', \
		datapath + 'thermo/Other_Tree.txt', \
		datapath + 'thermo/Other_Library.txt', True)
	
	logging.debug('\t\tOther functional groups')
	
	radicalDatabase = ThermoDatabase()
	radicalDatabase.load(datapath + 'thermo/Radical_Dictionary.txt', \
		datapath + 'thermo/Radical_Tree.txt', \
		datapath + 'thermo/Radical_Library.txt', True)
	
	logging.debug('\t\tRadical functional groups')
	
	ringDatabase = ThermoDatabase()
	ringDatabase.load(datapath + 'thermo/Ring_Dictionary.txt', \
		datapath + 'thermo/Ring_Tree.txt', \
		datapath + 'thermo/Ring_Library.txt', True)
	
	logging.debug('\t\tCyclic functional groups')
	
	primaryDatabase = ThermoDatabase()
	primaryDatabase.load(datapath + 'thermo/Primary_Dictionary.txt', \
		'', \
		datapath + 'thermo/Primary_Library.txt', False)

	logging.debug('\t\tPrimary thermo database')
	
	return int15Database, gaucheDatabase, groupDatabase, otherDatabase, \
	        radicalDatabase, ringDatabase, primaryDatabase

################################################################################

def saveNewThermoDatabases(datapath, int15Database, gaucheDatabase, \
		groupDatabase, otherDatabase, radicalDatabase, ringDatabase, \
		primaryDatabase):
	
	f = open(datapath + 'thermo/1,5-interactions.xml', 'w')
	f.write(int15Database.toXML())
	f.close()
	
	f = open(datapath + 'thermo/gauche-interactions.xml', 'w')
	f.write(gaucheDatabase.toXML())
	f.close()
	
	f = open(datapath + 'thermo/acyclics.xml', 'w')
	f.write(groupDatabase.toXML())
	f.close()
	
	f = open(datapath + 'thermo/other.xml', 'w')
	f.write(otherDatabase.toXML())
	f.close()
	
	f = open(datapath + 'thermo/radicals.xml', 'w')
	f.write(radicalDatabase.toXML())
	f.close()
	
	f = open(datapath + 'thermo/cyclics.xml', 'w')
	f.write(ringDatabase.toXML())
	f.close()
	
	f = open(datapath + 'thermo/primary.xml', 'w')
	f.write(primaryDatabase.toXML())
	f.close()
	
	
################################################################################

if __name__ == '__main__':
	
	datapath = '../data/RMG_database/'

	int15Database, gaucheDatabase, groupDatabase, otherDatabase, \
	        radicalDatabase, ringDatabase, primaryDatabase = \
			loadOldThermoDatabases(datapath)
			
	#saveNewThermoDatabases(datapath, int15Database, gaucheDatabase, \
	        #groupDatabase, otherDatabase, radicalDatabase, ringDatabase, \
			#primaryDatabase)
	
	
	#abrahamDatabase = loadDatabase(datapath + 'thermo/Abraham_Dictionary.txt', \
		#datapath + 'thermo/Abraham_Tree.txt', \
		#datapath + 'thermo/Abraham_Library.txt', True)
	
	#unifacDatabase = loadDatabase(datapath + 'thermo/Unifac_Dictionary.txt', \
		#datapath + 'thermo/Unifac_Tree.txt', \
		#datapath + 'thermo/Unifac_Library.txt', True)
	
	# Thermo data for H2
	#data = [0.0, 31.233, 6.895, 6.975, 6.994, 7.009, 7.081, 7.219, 7.720, 0, 0.0007, 0]
	#thermo = ThermoGAData()
	#thermo.fromDatabase(data, '')
	#Cp = thermo.heatCapacity(pq.Quantity(700, 'K')); Cp.units = 'kcal/(mol*K)'
	#print Cp

	