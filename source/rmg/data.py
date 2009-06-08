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

import os
import logging
import quantities as pq
import xml.dom.minidom

import chem
import reaction

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

class Dictionary(dict):
	"""
	An RMG dictionary class, extended from the Python dictionary class to
	include	functions for loading the dictionary from a file.
	"""

	def load(self, path):
		"""
		Parse an RMG database dictionary located at path. An RMG
		dictionary is a list of key-value pairs of a string	label and a string
		record. Each record is separated by at least one empty line.
		"""
		
		# The current record
		record = ''
		# The dictionary as read from the file
		records = {}
		
		# Process the dictionary
		try:
			
			fdict = open(path, 'r')
			for line in fdict:
				
				line = removeCommentFromLine(line).strip()
				
				# If at blank line, end of record has been found
				if len(line) == 0 and len(record) > 0:
					
					# Label is first line of record
					lines = record.splitlines()
					label = lines[0]
			
					# Add record to dictionary
					self[label] = record
					
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
		
	def toStructure(self):
		"""
		Convert the values stored in the dictionary from adjacency list strings
		to :class:`chem.Structure` objects. If a record is a union, it is stored
		as the string 'union', and automatically uses all immediate children of
		the node as the union.
		"""
	
		for label, record in self.iteritems():
		
			# If record is a union, then store as string 'union'
			# By definition a union includes all of its immediate children
			if record.lower().find('union') > -1:
				self[label] = 'union'
			# Otherwise convert adjacency list to structure
			else:
				structure = chem.Structure()
				structure.fromAdjacencyList(record)
				self[label] = structure

	def toXML(self, dom, root):
		"""
		Return an XML representation of the dictionary.
		"""
			
		dictionary = dom.createElement('dictionary')
		root.appendChild(dictionary)
	
		for label, structure in self.iteritems():
			entry = dom.createElement('entry')
			entry.setAttribute('label', label)
			dictionary.appendChild(entry)
	
			structure.toXML(dom, entry)		

################################################################################

class Tree:
	"""
	An implementation of an n-ary tree used for representing a hierarchy of
	data. Each instance contains dictionaries `parent` and `children` which use
	nodes as keys; the values are the node's parent and a list of the node's
	children, respectively.
	"""
	
	def __init__(self):
		self.top = []
		self.parent = {}
		self.children = {}
	
	def add(self, node, parent):
		"""
		Add `node` to the tree as a child of `parent`, another node	already in 
		the tree.
		"""
		
		# Set parent of node
		self.parent[node] = parent
		
		# Initialize list of children of node
		self.children[node] = []
		
		# Set node as child of parent
		if parent is not None:
			self.children[parent].append(node)
		
		# Set top if needed
		if parent is None:
			self.top.append(node)
			
	def remove(self, node):
		"""
		Remove `node` and all of its children from the tree.
		"""
		
		while len(self.children[node]) > 0:
			self.remove(self.children[node][0])
		if self.parent[node] is not None:
			self.children[self.parent[node]].remove(node)
		del self.parent[node]
		del self.children[node]
		
	def load(self, path):
		"""
		Parse an RMG database tree located at `path`. An RMG tree is an 
		n-ary tree representing the hierarchy of items in the dictionary.
		"""
		
		# An array of parents used when forming the tree
		parents = [None]
		
		# Process the tree (optional)
		try:
		
			ftree = open(path, 'r')
			for line in ftree:
				line = removeCommentFromLine(line).strip()
				if len(line) > 0:
					
					# Extract level
					data = line.split()
					level = int(data[0].replace('L', '').replace(':', ''))
					label = data[1]
					
					# Find immediate parent of the new node
					parent = None
					if len(parents) < level:
						raise InvalidDatabaseException('Invalid level specified.')
					else:
						while len(parents) > level:
							parents.remove(parents[-1])
						if len(parents) > 0:
							parent = parents[level-1]
					
					# Add node to tree
					self.add(label, parent)
					
					# Add node to list of parents
					parents.append(label)
					
								
		except InvalidDatabaseException, e:
			logging.exception(str(e))
		except IOError, e:
			logging.exception('Database tree file "' + e.filename + '" not found.')
		finally:	
			ftree.close()

	def toXML(self, dom, root):
		"""
		Return an XML representation of the tree.
		"""
		
		if self.top is None:
			return
		
		tree = dom.createElement('tree')
		root.appendChild(tree)
	
		nodes = {}
		for label in self.parent.keys():
			element = dom.createElement('node')
			element.setAttribute('label', label)
			nodes[label] = element
		
		for key, values in self.children.iteritems():
			for value in values:
				nodes[key].appendChild(nodes[value])
				
		tree.appendChild(nodes[self.top])
		
################################################################################

class Library(dict):
	"""
	An RMG database library class, extended from the base RMG dictionary class
	to be able to handle an array of string labels.
	"""
	
	def add(self, labels, data):
		"""
		Add an item of `data` to the library based on the value of the list
		of `labels`.
		"""
		names = self.hashLabels(labels)
		for name in names:
			self[name] = data
		
	def remove(self, labels):
		"""
		Remove an item of data from the library based on the value of the list
		of `labels`.
		"""
		names = self.hashLabels(labels)
		for name in names:
			del self[name]
	
	def hashLabels(self, labels):
		"""
		Convert a list of string `labels` to a list of single strings that
		represent permutations of the individual strings in the `labels` list.
		"""
		names = []
		if len(labels) == 1:
			names.append(labels[0])
		elif len(labels) == 2:
			names.append(labels[0] + ';' + labels[1])
			names.append(labels[1] + ';' + labels[0])
		elif len(labels) == 3:
			names.append(labels[0] + ';' + labels[1] + ';' + labels[2])
			names.append(labels[0] + ';' + labels[2] + ';' + labels[1])
			names.append(labels[1] + ';' + labels[0] + ';' + labels[2])
			names.append(labels[1] + ';' + labels[2] + ';' + labels[0])
			names.append(labels[2] + ';' + labels[0] + ';' + labels[1])
			names.append(labels[2] + ';' + labels[1] + ';' + labels[0])
		return names
	
	#def __getitem__(self, key):
		#if key.__class__ == str or key.__class__ == unicode:
			#return self[key]
		#else:
			#names = self.hashLabels(labels)
			#return self[names[0]]
	
	def load(self, path, numLabels=1):
		"""
		Parse an RMG database library located at `path`.
		"""
		
		# Process the library
		try:
		
			flib = open(path, 'r')
			for line in flib:
				line = removeCommentFromLine(line).strip()
				if len(line) > 0:
					
					info = line.split()
					
					# Skip if the number of items on the line is invalid
					if len(info) < 2:
						continue
					
					# Extract label(s)
					labels = []
					for i in range(0, numLabels):
						labels.append(info[i+1])
					comment = ''
					
					data = ''
					for i in range(numLabels+1, len(info)):
						data += info[i] + ' '
					
					self.add(labels, data)
					
		except InvalidDatabaseException, e:
			logging.exception(str(e))
			print dictstr, treestr, libstr
			quit()
		except IOError, e:
			logging.exception('Database library file "' + e.filename + '" not found.')
		finally:	
			flib.close()
			
	def removeLinks(self):
		"""
		Ensure all values in the library are either a 
		:class:`chem.ThermoGAValue` object or None by following and replacing
		all string links.
		"""
		
		for label, data in self.iteritems():
		
			dataLabel = ''
			while data.__class__ == str or data.__class__ == unicode:
				if data not in self:
					raise InvalidDatabaseException('Node references parameters from a node that is not in the library.')
				if dataLabel == data:
					raise InvalidDatabaseException('Node references parameters from itself.')
				dataLabel = data; data = self[dataLabel]
			
			self[label] = data

	def toXML(self, dom, root):
		"""
		Return an XML representation of the library. 
		"""
		
		library = dom.createElement('library')
		root.appendChild(library)
	
		for label, data in self.iteritems():
			element = dom.createElement('data')
			element.setAttribute('label', label)
			library.appendChild(element)

			if data.__class__ == str or data.__class__ == unicode:
				link = dom.createElement('link')
				link.setAttribute('target', data)
				element.appendChild(link)
			else:
				data.toXML(dom, element)
			
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
		self.dictionary = Dictionary()
		self.library = Library()
		self.tree = Tree()
		
	def load(self, dictstr, treestr, libstr):
		"""
		Load a dictionary-tree-library based database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. The tree is optional, and should
		be set to '' if not desired.
		"""
		
		# Load dictionary, library, and (optionally) tree
		self.dictionary.load(dictstr)
		self.dictionary.toStructure()
		self.library.load(libstr, 1)
		if treestr != '': self.tree.load(treestr)

	def toXML(self, dom, root):
		"""
		Return an XML representation of the database.
		"""
		self.dictionary.toXML(dom, root)
		self.tree.toXML(dom, root)
		self.library.toXML(dom, root)

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

if __name__ == '__main__':
	pass
	