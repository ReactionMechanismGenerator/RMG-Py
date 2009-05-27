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
import chem
import fgroup

################################################################################

class InvalidDatabaseException(Exception):
	"""
	An exception used when parsing an RMG database to indicate that the
	databas is invalid. The `msg` parameter is used to specify what about the 
	database caused the exception to be raised.
	"""	

	def __init__(self, msg):
		self.msg = msg
	
	def __str__(self):
		return 'Invalid format for RMG database: ' + self.msg

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

	def addNode(self, node, parent):
		"""
		Add `node` to the database as a child of `parent`, another node
		already in the database.
		"""
		
		if parent is None:
			if self.top is None:
				self.top = node
			else:
				raise InvalidDatabaseException('Database can only have one parent node.')
		elif parent.label not in self.nodes:
			raise InvalidDatabaseException('Parent node is not in database.')
		
		self.nodes[node.label] = node
		if parent is not None:
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
					
					# Label is first line
					lines = record.splitlines()
					label = lines[0]
			
					# Add record to temporary dictionary
					records[label] = record
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

		# Process the tree
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
					if label not in records:
						raise InvalidDatabaseException('Node "' + label + \
							'" in tree, but not in database.')
							
					# Create new node
					node = Node(None, label, records[label], None, '')
					
					# Find immediate parent of the new node
					if len(parents) < level:
						raise InvalidDatabaseException('Invalid level specified.')
					else:
						while len(parents) > level:
							parents.remove(parents[-1])
						if len(parents) > 0:
							node.parent = parents[level-1]
						
					# Formally add node to tree
					self.addNode(node, node.parent)
					parents.append(node)
					
								
		except InvalidDatabaseException, e:
			logging.exception(str(e))
		except IOError, e:
			logging.exception('Database tree file "' + e.filename + '" not found.')
		finally:	
			ftree.close()

		# Process the library
		try:
		
			flib = open(libstr, 'r')
			for line in flib:
				line = removeCommentFromLine(line).strip()
				if len(line) > 0:
					
					# Extract level
					info = line.split()
					label = info[1]
					
					if len(info) == 3:
						# Data is label to other node that contains the data to use
						data = info[2]
						comment = ''
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
						comment = comment[:-1]
					
					# Check that tree has associated record from dictionary
					if label not in self.nodes:
						raise InvalidDatabaseException('Node "' + label + \
							'" in library, but not in tree or database.')
					
					# Set data of node (no further processing here, since at
					# this point we don't know anything about the data)
					self.nodes[label].data = data				
					self.nodes[label].comment = comment				
					
		except InvalidDatabaseException, e:
			logging.exception(str(e))
		except IOError, e:
			logging.exception('Database library file "' + e.filename + '" not found.')
		finally:	
			flib.close()


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

def parseAdjacencyList(adjlist, isfgroup=True):
	"""
	Convert a string adjacency list `adjlist` into a functional group object.
	The optional argument `isfgroup` specifies whether the adjacency tree
	represents a functional group or a full chemical species.
	"""
	
	atoms = []; bonds = []; atomdict = {}
			
	lines = adjlist.splitlines()
	label = lines[0]
	for line in lines[1:]:
		
		data = line.split()
		
		# First item is index for atom
		aid = int(data[0])
		
		# If second item is '*', the atom is the center atom
		center = False
		index = 1
		if data[1] == '*':
			center = True
			index = 2
		
		# Next is the element or functional group element
		element = data[index]
		if element == '{Cd,CO}':
			element = 'Cd'
			
		# Next is the electron state
		elecState = data[index+1]
		
		# Create a new atom based on the above information
		if isfgroup:
			atom = fgroup.Atom(element, elecState, 0, center)
		else:
			atom = chem.Atom(element, elecState, 0)
			
		# Add the atom to the list
		atoms.append(atom)
		atomdict[aid] = atom
		
		# Process list of bonds
		for datum in data[index+2:]:
			aid2, btype = datum.replace('{', '').replace('}', '').replace(',',' ').split()
			aid2 = int(aid2)
			if aid2 in atomdict:
				if isfgroup:
					bond = fgroup.Bond([atomdict[aid], atomdict[aid2]], btype)
				else:
					bond = chem.Bond([atomdict[aid], atomdict[aid2]], btype)
				bonds.append(bond)
		
	# Create and return functional group or species
	if isfgroup:
		structure = fgroup.Structure(atoms, bonds)
		return fgroup.FunctionalGroup(label, structure)
	else:
		structure = chem.Structure(atoms, bonds)
		return chem.Species(label, structure)

################################################################################

if __name__ == '__main__':

	database = Database()
	
	database.readFromFiles('../data/RMG_database/thermo/Group_Dictionary.txt', \
		'../data/RMG_database/thermo/Group_Tree.txt', \
		'../data/RMG_database/thermo/Group_Library.txt')
	
	for label, node in database.nodes.iteritems():
		node.structure = parseAdjacencyList(node.structure)
		