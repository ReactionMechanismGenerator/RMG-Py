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
Contains classes and functions for working with the various RMG databases.
"""

import logging
import quantities as pq
import scipy.interpolate

import chem
import fgroup

pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol') 

################################################################################

class InvalidDatabaseException(Exception):
	"""
	An exception used when parsing an RMG database to indicate that the
	database is invalid. The `msg` parameter is used to specify what about the 
	database caused the exception to be raised.
	"""	

	def __init__(self, msg):
		self.msg = msg
	
	def __str__(self):
		return 'Invalid format for RMG database: ' + self.msg

class TemperatureOutOfRangeException(Exception):
	"""
	An exception used when parsing an RMG database to indicate that a
	temperature is outside the allowed range. At the time of this writing the
	allowed range was T >= 300 K.
	"""	

	def __init__(self, msg):
		self.msg = msg
	
	def __str__(self):
		return self.msg

################################################################################

class Node:
	"""
	Represent a node in an RMG database. As the database is represented by an
	n-ary tree, each object has a `parent` and a list of `children`. Each node
	in a database should also have a unique `label`, an associated chemical
	or functional group `structure`, a list of `data` values, and a `comment`
	describing the source of the data values (or a place to go to find that
	information).
	"""

	def __init__(self, parent=None, label='', structure=None, data=None, \
	             comment='', children=[]):
		"""
		Initialize an RMG database node.
		"""
		self.parent = parent
		self.label = label
		self.structure = structure
		self.data = data
		self.comment = comment
		self.children = children
        
################################################################################

class Database:
	"""
	Represent an RMG database. An RMG database is structured as an n-ary tree,
	were each node has a chemical or functional group `structure` (that is a
	superset of its parent node) and a library of `data` values. This class is 
	intended to be generic and can be subclassed for specific databases if
	desired.

	Each Database object maintains its own internal dictionary of the nodes
	in the tree. Thus it is strongly recommended that you utilize the add()
	and remove() functions for manipulating the database.
	"""

	def __init__(self):
		"""
		Initialize an RMG database.
		"""
		self.top = None
		self.nodes = {}

	def addNodeToTree(self, node, parent):
		"""
		Add `node` to the database as a child of `parent`, another node
		already in the database.
		"""
		
		node.children = []
		
		if parent is None:
			if self.top is None:
				self.top = node
			else:
				raise InvalidDatabaseException('Database can only have one parent node.')
		elif parent.label not in self.nodes:
			raise InvalidDatabaseException('Parent node is not in database.')
		else:
			parent.children.append(node)
			
	def removeNode(self, node):
		"""
		Remove `node` and all of its children from the database.
		"""
		
		if node not in self.nodes:
			raise InvalidDatabaseException('Node to be removed is not in the database.')
		else:
			while len(node.children) > 0:
				self.remove(node.children[0])
			del self.nodes[node.label]
			node.parent.children.remove(node)


	def readFromFiles(self, dictstr, treestr, libstr):
		"""
		Parse an RMG database. On disk each database is made up of three files: 
		
		1. The *dictionary*, a list of key-value pairs of a string label and a 
		string record; each record is separated by at least one empty line
		
		2. The *tree*, an n-ary tree representing the hierarchy of items in the
		dictionary
		
		3. The *library*, associating each dictionary item with a set of parameter
		values.
		
		The file at `dictstr` is treated as the dictionary, `treestr` the 
		corresponding tree, and `libstr` the corresponding library.
		"""

		# The current record
		record = ''
		# The dictionary as read from the file
		records = {}
		# An array of parents used when forming the tree
		parents = []
		
		# Process the dictionary
		try:
			
			fdict = open(dictstr, 'r')
			for line in fdict:
				
				line = removeCommentFromLine(line).strip()
				
				# If at blank line, end of record has been found
				if len(line) == 0 and len(record) > 0:
					
					# Label is first line of record
					lines = record.splitlines()
					label = lines[0]
			
					# Add record to temporary dictionary
					node = Node(None, label, record, None, '')
					self.nodes[label] = node
					
					# Clear record in preparation for next iteration
					record = ''
					
				# Otherwise append line to record (if not empty)
				elif len(line) > 0:
					record += line + '\n'
					
		except InvalidDatabaseException, e:
			logging.exception(str(e))
			return
		except IOError, e:
			logging.exception('Database dictionary file "' + e.filename + '" not found.')
			return
		finally:	
			fdict.close()

		# Process the library
		try:
		
			flib = open(libstr, 'r')
			for line in flib:
				line = removeCommentFromLine(line).strip()
				if len(line) > 0:
					
					# Extract level
					info = line.split()
					label = info[1]
					comment = ''
					
					if len(info) == 3:
						# Data is label to other node that contains the data to use
						data = info[2]
					else:
						# Parse the data the node contains; numbers are added 
						# from info[2:] to data until an entry is found that 
						# cannot be converted to a number; all subsequent 
						# entries are treated as supplemental comments
						data = []
						index = 2
						try:
							while index < len(info):
								data.append(float(info[index]))
								index += 1
						except ValueError, e:
							pass
						for i in range(index, len(info)):
							comment += info[i] + ' '
						if len(comment) > 0:
							comment = comment[:-1]
					
					# Check that tree has associated record from dictionary
					if label not in self.nodes:
						raise InvalidDatabaseException('Node "' + label + \
							'" in library, but not in dictionary.')
						
					# Set data of node (no further processing here, since at
					# this point we don't know anything about what the data
					# represents)
					self.nodes[label].data = data				
					self.nodes[label].comment = comment				
					
		except InvalidDatabaseException, e:
			logging.exception(str(e))
			print dictstr, treestr, libstr
			quit()
		except IOError, e:
			logging.exception('Database library file "' + e.filename + '" not found.')
		finally:	
			flib.close()
			
		# Process the tree (optional)
		if treestr != '':
		
			try:
			
				ftree = open(treestr, 'r')
				for line in ftree:
					line = removeCommentFromLine(line).strip()
					if len(line) > 0:
						
						# Extract level
						data = line.split()
						level = int(data[0].replace('L', '').replace(':', ''))
						label = data[1]
						
						# Check that tree has associated record from dictionary
						if label not in self.nodes:
							raise InvalidDatabaseException('Node "' + label + \
								'" in tree, but not in dictionary.')
								
						# Create new node
						node = self.nodes[label]
						
						# Find immediate parent of the new node
						if len(parents) < level:
							raise InvalidDatabaseException('Invalid level specified.')
						else:
							while len(parents) > level:
								parents.remove(parents[-1])
							if len(parents) > 0:
								node.parent = parents[level-1]
						
						# Formally add node to tree
						self.addNodeToTree(node, node.parent)
						
						parents.append(node)
						
									
			except InvalidDatabaseException, e:
				logging.exception(str(e))
			except IOError, e:
				logging.exception('Database tree file "' + e.filename + '" not found.')
			finally:	
				ftree.close()

		else:
			
			for label, record in records.iteritems():
				node = Node(None, label, record, None, '')
				self.nodes[label] = node

		


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
		self.H298.units = 'J/mol'
		self.S298 = pq.UncertainQuantity(S298, 'cal/(mol*K)', dS)
		self.S298.units = 'J/(mol*K)'
		self.Cp[0] = pq.UncertainQuantity(Cp300, 'kcal/(mol*K)', dCp)
		self.Cp[0].units = 'J/(mol*K)'
		self.Cp[1] = pq.UncertainQuantity(Cp400, 'kcal/(mol*K)', dCp)
		self.Cp[1].units = 'J/(mol*K)'
		self.Cp[2] = pq.UncertainQuantity(Cp500, 'kcal/(mol*K)', dCp)
		self.Cp[2].units = 'J/(mol*K)'
		self.Cp[3] = pq.UncertainQuantity(Cp600, 'kcal/(mol*K)', dCp)
		self.Cp[3].units = 'J/(mol*K)'
		self.Cp[4] = pq.UncertainQuantity(Cp800, 'kcal/(mol*K)', dCp)
		self.Cp[4].units = 'J/(mol*K)'
		self.Cp[5] = pq.UncertainQuantity(Cp1000, 'kcal/(mol*K)', dCp)
		self.Cp[5].units = 'J/(mol*K)'
		self.Cp[6] = pq.UncertainQuantity(Cp1500, 'kcal/(mol*K)', dCp)
		self.Cp[6].units = 'J/(mol*K)'
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
		
		
################################################################################

class ThermoDatabase(Database):
	"""
	Represent an RMG thermodynamics database. 
	"""

	def __init__(self):
		Database.__init__(self)
		

	def load(self, dictstr, treestr, libstr, isfgroup=True):
		"""
		Load a thermodynamics group additivity database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. The tree is optional, and should
		be set to '' if not desired.
		"""
		
		self.readFromFiles(dictstr, treestr, libstr)
		for label, node in self.nodes.iteritems():
			
			# Convert adjacency list to structure
			if isfgroup:
				structure = fgroup.Structure()
			else:
				structure = chem.Structure()
			structure.fromAdjacencyList(node.structure)
			node.structure = structure
			
			# Convert list of data to object
			thermoData = ThermoGAData()
			# Some nodes have string labels that act as pointers to other
			# nodes' data lists; these are navigated now so that they don't 
			# have to be in the future
			if node.data.__class__ == str or node.data.__class__ == unicode:
				dataNode = node
				while dataNode.data.__class__ == str or dataNode.data.__class__ == unicode:
					if dataNode.data not in self.nodes:
						raise InvalidDatabaseException('Node references parameters from a node that is not in the library.')
					if dataNode.data == dataNode.label:
						raise InvalidDatabaseException('Node references parameters from itself.')
					dataNode = self.nodes[dataNode.data]
				if dataNode.data.__class__ != ThermoGAData:
					self.nodeToThermoGAData(dataNode)
				node.data = dataNode.data
			# Some nodes have no data specified; this is not an exceptional
			# circumstance and should simply be ignored
			elif node.data.__class__ != ThermoGAData:
				self.nodeToThermoGAData(node)
	
	def nodeToThermoGAData(self, node):
		"""
		Convert the data of `node` from a list to a :class:`ThermoGAData` 
		object. The node must not contain a string label pointing to another
		node; it can contain either a list of 12 numbers or :data:`None`.
		"""
		if node.data is None:
			pass
		elif len(node.data) != 12:
			pass
		else:
			thermoData = ThermoGAData()
			thermoData.fromDatabase(node.data, node.comment)
			node.data = thermoData
	
################################################################################

def removeCommentFromLine(line):
	"""
	Remove a C++/Java style comment from a line of text.
	"""

	index = line.find('//')
	if index >= 0:
		line = line[0:index]
	return line

################################################################################

def loadThermoDatabases(datapath):
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

	int15Database = ThermoDatabase()
	int15Database.load(datapath + 'thermo/15_Dictionary.txt', \
		datapath + 'thermo/15_Tree.txt', \
		datapath + 'thermo/15_Library.txt', True)
	
	logging.debug('\t1,5 interactions: ' + str(len(int15Database.nodes)) + ' nodes')
	
	gaucheDatabase = ThermoDatabase()
	gaucheDatabase.load(datapath + 'thermo/Gauche_Dictionary.txt', \
		datapath + 'thermo/Gauche_Tree.txt', \
		datapath + 'thermo/Gauche_Library.txt', True)
	
	logging.debug('\tGauche interactions: ' + str(len(gaucheDatabase.nodes)) + ' nodes')
	
	groupDatabase = ThermoDatabase()
	groupDatabase.load(datapath + 'thermo/Group_Dictionary.txt', \
		datapath + 'thermo/Group_Tree.txt', \
		datapath + 'thermo/Group_Library.txt', True)
	
	logging.debug('\tAcyclic functional groups: ' + str(len(groupDatabase.nodes)) + ' nodes')
	
	otherDatabase = ThermoDatabase()
	otherDatabase.load(datapath + 'thermo/Other_Dictionary.txt', \
		datapath + 'thermo/Other_Tree.txt', \
		datapath + 'thermo/Other_Library.txt', True)
	
	logging.debug('\tOther functional groups: ' + str(len(otherDatabase.nodes)) + ' nodes')
	
	radicalDatabase = ThermoDatabase()
	radicalDatabase.load(datapath + 'thermo/Radical_Dictionary.txt', \
		datapath + 'thermo/Radical_Tree.txt', \
		datapath + 'thermo/Radical_Library.txt', True)
	
	logging.debug('\tRadical functional groups: ' + str(len(radicalDatabase.nodes)) + ' nodes')
	
	ringDatabase = ThermoDatabase()
	ringDatabase.load(datapath + 'thermo/Ring_Dictionary.txt', \
		datapath + 'thermo/Ring_Tree.txt', \
		datapath + 'thermo/Ring_Library.txt', True)
	
	logging.debug('\tCyclic functional groups: ' + str(len(ringDatabase.nodes)) + ' nodes')
	
	primaryDatabase = ThermoDatabase()
	primaryDatabase.load(datapath + 'thermo/Primary_Dictionary.txt', \
		'', \
		datapath + 'thermo/Primary_Library.txt', False)

	logging.debug('\tPrimary thermo database: ' + str(len(primaryDatabase.nodes)) + ' nodes')
	
	return int15Database, gaucheDatabase, groupDatabase, otherDatabase, \
	        radicalDatabase, ringDatabase, primaryDatabase

################################################################################

if __name__ == '__main__':

	#int15Database, gaucheDatabase, groupDatabase, otherDatabase, \
	        #radicalDatabase, ringDatabase, primaryDatabase = \
			#loadThermoDatabases('../data/RMG_database/')
			
	#print int15Database.top.children
	#for key, value in int15Database.nodes.iteritems():
		#print key, value
		#for child in value.children:
			#print child.label
	#quit()
	
	#abrahamDatabase = loadDatabase(datapath + 'thermo/Abraham_Dictionary.txt', \
		#datapath + 'thermo/Abraham_Tree.txt', \
		#datapath + 'thermo/Abraham_Library.txt', True)
	
	#unifacDatabase = loadDatabase(datapath + 'thermo/Unifac_Dictionary.txt', \
		#datapath + 'thermo/Unifac_Tree.txt', \
		#datapath + 'thermo/Unifac_Library.txt', True)
	
	# Thermo data for H2
	data = [0.0, 31.233, 6.895, 6.975, 6.994, 7.009, 7.081, 7.219, 7.720, 0, 0.0007, 0]
	thermo = ThermoGAData()
	thermo.fromDatabase(data, '')
	Cp = thermo.heatCapacity(pq.Quantity(700, 'K')); Cp.units = 'kcal/(mol*K)'
	print Cp
	