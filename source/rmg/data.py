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
Contains classes and functions for working with the various RMG databases. In
RMG a database is composed of three parts:

* A *dictionary*, associating a string identifier with a chemical structure.

* A *tree*, providing a hierarchical and extensible organization of the chemical structures.

* A *library*, associating a chemical structure with some sort of physical information.

These are implemented via the :class:`Dictionary`, :class:`Tree`, and
:class:`Library` classes, respectively. A base class for databases is
implemented in the :class:`Database` class, which is often overloaded to
provide specific functionality for individual databases.
"""

import os
import log as logging
import quantities as pq
import re
import structure

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
	include functions for loading the dictionary from a file. The keys of the
	dictionary are strings that represent unique identifiers, while the
	corresponding values are either :class:`chem.Structure` objects representing
	a chemical structure or the string 'union' to indicate that the structure
	is a union of all of its child nodes in the corresponding tree.
	"""

	def load(self, path):
		"""
		Parse an RMG database dictionary located at path. An RMG
		dictionary is a list of key-value pairs of a string	label and a string
		record. Each record is separated by at least one empty line.
		"""
		
		# The current record
		record = ''
		
		fdict=None
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
			
			# process the last record! (after end of for loop)
			# Label is first line of record
			if record:
				lines = record.splitlines()
				label = lines[0]
				# Add record to dictionary
				self[label] = record
			
		except InvalidDatabaseException, e:
			logging.exception(str(e))
			raise
		except IOError, e:
			logging.exception('Database dictionary file "' + e.filename + '" not found.')
			raise
		finally:	
			if fdict: fdict.close()
		
	def toStructure(self, addH=False):
		"""
		Convert the values stored in the dictionary from adjacency list strings
		to :class:`structure.Structure` objects. If a record is a union, it is 
		stored as the string 'union', and automatically uses all immediate
		children of the node as the union.
		If a record is a logical node, it is converted into the appropriate class.
		"""
	
		for label, record in self.iteritems():
		
			# If record is a union, 
			# If record is a logical node, make it into one.
			if re.match('(?i)\s*OR|AND|NOT|UNION',record.splitlines()[1] ):
				self[label] = makeLogicNode(' '.join(record.splitlines()[1:]) )
			# Otherwise convert adjacency list to structure
			else:
				try:
					struct = structure.Structure()
					struct.fromAdjacencyList(record, addH)
					self[label] = struct
				except structure.InvalidAdjacencyListException, e:
					logging.error('\t\t\t' + str(e))
					raise

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

def getAllCombinations(nodeLists):
	"""
	Generate a list of all possible combinations of items in the list of
	lists `nodeLists`. Each combination takes one item from each list
	contained within `nodeLists`. The order of items in the returned lists 
	reflects the order of lists in `nodeLists`. For example, if `nodeLists` was
	[[A, B, C],	[N], [X, Y]], the returned combinations would be
	[[A, N, X], [A, N, Y], [B, N, X], [B, N, Y], [C, N, X], [C, N, Y]].
	"""

	items = [[]]
	for nodeList in nodeLists:
		items = [ item + [node] for node in nodeList for item in items ]

	return items

################################################################################

class Tree:
	"""
	An implementation of an n-ary tree used for representing a hierarchy of
	data. The tree is represented as a pair of dictionaries:

	* The `parent` dictionary. For a given identifier, returns the identifier
	  corresponding to the parent.

	* The `children` dictionary. For a given identifier, returns a list of
	  identifiers corresponding to the children.

	The top-level node(s) of the tree are stored in a list in the `top`
	attribute.

	In constructing the tree, it is important to develop a hierarchy such that
	siblings are mutually exclusive, to ensure that there is a unique path of
	descent down a tree for each structure. If non-mutually exclusive siblings
	are encountered, a warning is raised and the parent of the siblings is
	returned.

	There is no requirement that the children of a node span the range of
	more specific permutations of the parent. As the database gets more complex,
	attempting to maintain complete sets of children for each parent in each
	database rapidly becomes untenable, and is against the spirit of
	extensibility behind the database development.
	"""
	
	def __init__(self):
		self.top = []
		self.parent = {}
		self.children = {}
		
	def ancestors(self, node):
		"""
		Returns all the ancestors of a node, climbing up the tree to the top.
		"""
		parent=self.parent[node]
		if parent is None:
			return list()
		else:
			ancestors = [parent]
			ancestors.extend(self.ancestors(parent))
			return ancestors
	
	def descendants(self, node):
		"""
		Returns all the descendants of a node, climbing down the tree to the bottom.
		"""
		children = self.children[node][:]
		temp0 = children[:]

		while len(temp0) > 0:
			temp = []
			for child in temp0:
				children.extend(self.children[child])
				temp.extend(self.children[child])
			temp0 = temp

		return children

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
		# Recursively remove the children of node
		while len(self.children[node]) > 0:
			self.remove(self.children[node][0])
			
		# Remove the current node from the list of its parent's children
		if self.parent[node] is not None:
			self.children[self.parent[node]].remove(node)

		# Delete the node from the parent and children dictionaries
		del self.parent[node]
		del self.children[node]
		
	def load(self, path):
		"""
		Parse an RMG database tree located at `path`. An RMG tree is an 
		n-ary tree representing the hierarchy of items in the dictionary.
		"""

		# An array of parents used when forming the tree
		parents = [None]
		
		import re
		parser = re.compile('^\s*L(?P<level>\d+)\s*:\s*(?P<label>\S+)')
		# should match '  L3 : foo_bar '  and 'L3:foo_bar'
		
		# Process the tree (optional)
		try:
			ftree = open(path, 'r')
			for line in ftree:
				line = removeCommentFromLine(line).strip()
				if len(line) > 0:
					# Extract level
					match = parser.match(line)
					if not match: 
						raise InvalidDatabaseException("Couldn't parse line '%s'"%line.strip() ) 
					level = int(match.group('level'))
					label = match.group('label')
					
					# Find immediate parent of the new node
					parent = None
					if len(parents) < level:
						raise InvalidDatabaseException("Invalid level specified in line '%s'"%line.strip() )
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

	def toDOT(self):
		"""
		Create a DOT representation of the tree.
		"""
		
		import pydot

		graph = pydot.Dot(size='10,8', page='10,8' ,  rankdir='LR',
				graph_type='digraph', simplify=True, fontsize=10,
				overlap='true', dpi='85',center="True")

		# Create vertices of digraph
		for label in self.parent:
			node = pydot.Node(label)
			graph.add_node(node)
		
		# Create edges of digraph
		for key, value in self.parent.iteritems():
			if value is not None:
				graph.add_edge(pydot.Edge(value,key))

		return graph

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

	def write(self, children):
		"""
		Write the tree to a string in the syntax used by the RMG database. The
		`children` parameter is a list of the current children, used to enable
		recursive writing.
		"""

		string = ''

		for child in children:

			# Determine level
			level = 1
			temp = child
			while self.parent[temp] is not None:
				level += 1
				temp = self.parent[temp]

			# Write current node
			for i in range(level):
				string += '\t'
			string += 'L%s: %s\n' % (str(level), child)

			# Recursively descend children
			string += self.write(self.children[child])

		return string



################################################################################

class Library(dict):
	"""
	An RMG database library class, extended from the base RMG dictionary class
	to be able to handle an array of string labels.
	"""
	
	def add(self, index, labels, data):
		"""
		Add an item of `data` to the library based on the value of the list
		of `labels`. Only add and return True if there is not preexisting data 
		with those labels, else return False.
		"""
		if self.getData(labels) is not None:
			logging.warning("There was already something labelled %s in the database. Ignoring '%s' (%s)"%(labels,index, data))
			return False
		names = self.hashLabels(labels)
		
		for name in names:
			self[name] = (index, data )

		return True
		
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
		represent permutations of the individual strings in the `labels` list::
		
			>>> hashLabels(['a','b'])
			['a;b', 'b;a']
		"""
		names = []
		if len(labels) == 1:
			names.append(labels[0])
		elif len(labels) == 2:
			names.append(labels[0] + ';' + labels[1])
		#	names.append(labels[1] + ';' + labels[0])
		elif len(labels) == 3:
			names.append(labels[0] + ';' + labels[1] + ';' + labels[2])
		#	names.append(labels[0] + ';' + labels[2] + ';' + labels[1])
		#	names.append(labels[1] + ';' + labels[0] + ';' + labels[2])
		#	names.append(labels[1] + ';' + labels[2] + ';' + labels[0])
		#	names.append(labels[2] + ';' + labels[0] + ';' + labels[1])
		#	names.append(labels[2] + ';' + labels[1] + ';' + labels[0])
		return names
	
	def getData(self, key):
		"""
		Return the data in the library associated with the label or list of
		labels denoted by `key`.
		"""
		if key.__class__ == str or key.__class__ == unicode:
			if key in self: return self[key]
		else:
			names = self.hashLabels(key)
			for name in names:
				if name in self: return self[name]
		return None
	
	def load(self, path):
		"""
		Parse an RMG database library located at `path`.
		"""
		lines = []
		
		# Process the library
		flib = None
		try:
			flib = open(path, 'r')
			for line in flib:
				line = removeCommentFromLine(line).strip()
				if len(line) > 0:
					lines.append(line)
					
		except IOError, e:
			logging.exception('Database library file "' + e.filename + '" not found.')
			raise
		finally:	
			if flib: flib.close()
		
		return lines

	def parse(self, lines, numLabels=1):
		"""
		Parse an RMG database library located at `path`.
		
		It splits lines on whitespace then treats tokens as::
		
			<index> <label1> ... <labelN> <data1>  <data2> ...
		
		`numLabels` determines how  many labels are assumed.
		All the data are concatenated to a single string with single spaces between items.
		"""
		
		# Process the library
		try:
			for line in lines:
				info = line.split()
				
				# Skip if the number of items on the line is invalid
				if len(info) < 2:
					continue
				
				# Extract label(s)
				labels = []
				for i in range(0, numLabels):
					labels.append(info[i+1])
				
				data = ''
				for i in range(numLabels+1, len(info)):
					data += info[i] + ' '
				
				index = info[0]
				
				self.add(index, labels, data)
					
		except InvalidDatabaseException, e:
			logging.exception(str(e))
			for s in (dictstr, treestr, libstr):
				logging.exception(s)
			quit()
		
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
					raise InvalidDatabaseException('Node %s references parameters from a node %s that is not in the library.'%(label,data))
				if dataLabel == data:
					raise InvalidDatabaseException('Node %s references parameters from itself.'%label )
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

def makeLogicNode(string):
	"""
	Creates and returns a node in the tree which is a logic node.
	
	String should be of the form:
	
	* OR{}
	* AND{}
	* NOT OR{}
	* NOT AND{} 
	
	And the returned object will be of class LogicOr or LogicAnd
	"""
	
	match = re.match("(?i)\s*(NOT)?\s*(OR|AND|UNION)\s*(.*)",string)  # the (?i) makes it case-insensitive
	if not match:
		raise Exception("Unexpected string for Logic Node: %s"%string)
		
	if match.group(1): invert = True
	else: invert = False
	
	logic = match.group(2)  # OR or AND
	
	contents = match.group(3).strip()
	while contents.startswith('{'):
		if not contents.endswith('}'):
			raise Exception("Unbalanced braces in Logic Node: %s"%string)
		contents = contents[1:-1]
	
	items=[]
	chars=[]
	brace_depth = 0
	for character in contents:
		if character == '{':
			brace_depth += 1
		if character == '}':
			brace_depth -= 1
		if character == ',' and brace_depth == 0:
			items.append(''.join(chars).lstrip().rstrip() )
			chars = []
		else: 
			chars.append(character)
	if chars: # add last item
		items.append(''.join(chars).lstrip().rstrip() ) 
	if brace_depth != 0: raise Exception("Unbalanced braces in Logic Node: %s"%string)
		
	logging.debug("Creating %s Logic Node with items: %s"%(logic,items))
	
	if logic.upper() in ['OR', 'UNION']:
		return LogicOr(items, invert)
	if logic == 'AND':
		return LogicAnd(items, invert)
		
	raise Exception("Could not create Logic Node from %s"%string)
		

class LogicNode():
	"""
	A base class for AND and OR logic nodes. 
	"""
	symbol="<TBD>" # To be redefined by subclass
	def __init__(self,items,invert):
		self.components = []
		for item in items:
			if re.match('(?i)\s*OR|AND|NOT|UNION',item):
				component = makeLogicNode(item)
			else:
				component = item
			self.components.append(component)
		self.invert = bool(invert)
		logging.debug("Created Logic Node: "+str(self) )
	def __str__(self):
		result = ''
		if self.invert: self += 'NOT '
		result += self.symbol
		result += "{%s}"%(', '.join([str(c) for c in self.components]))
		return result

class LogicOr(LogicNode):
	"""
	A logical OR node. Structure can match any component.
	
	Initialize with a list of component items and a boolean instruction to invert the answer.
	"""
	symbol = "OR"
	def matchToStructure(self,database,structure,atoms):
		"""
		Does this node in the given database match the given structure with the labeled atoms?
		"""
		for node in self.components:
			if isinstance(node,LogicNode):
				match = node.matchToStructure(database,structure,atoms)
			else:
				match = database.matchNodeToStructure(node, structure, atoms)
			if match:
				return True != self.invert
		return False != self.invert
	
class LogicAnd(LogicNode):
	"""A logical AND node. Structure must match all components."""
	symbol = "AND"
	def matchToStructure(self,database,structure,atoms):
		"""
		Does this node in the given database match the given structure with the labeled atoms?
		"""
		for node in self.components:
			if isinstance(node,LogicNode):
				match = node.matchToStructure(database,structure,atoms)
			else:
				match = database.matchNodeToStructure(node, structure, atoms)
			if not match:
				return False != self.invert
		return True != self.invert 
	
		
			
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
		assert len(self.dictionary)>0
		
		if treestr != '':
			self.tree.load(treestr)
			# Check that all nodes in tree are also in dictionary
			for node in self.tree.children:
				if node not in self.dictionary:
					if node.startswith('Others-'):
						logging.warning("Node %s found in tree but not dictionary %s. Letting it through for now because it's an 'Others-' node!"%(node,os.path.abspath(treestr)))
					else: 
						raise InvalidDatabaseException('Node "' + node + '" found in tree "' + os.path.abspath(treestr) + '", but not in corresponding dictionary "' + os.path.abspath(dictstr) + '".')
			# Sort children by decreasing size; the algorithm returns the first
			# match of each children, so this makes it less likely to miss a
			# more detailed functional group
			# First determine if we can do the sort (that is, all children have
			# one structure.Structure)
			canSort = True
			for node, children in self.tree.children.iteritems():
				for child in children:
					if self.dictionary[child].__class__ != structure.Structure: canSort = False
			if canSort:
				for node, children in self.tree.children.iteritems():
					children.sort(lambda x, y: cmp(len(self.dictionary[x].atoms()), len(self.dictionary[y].atoms())))
		if libstr != '':
			lines = self.library.load(libstr)
			self.library.parse(lines, 1)


	def save(self, dictstr, treestr, libstr):
		"""
		Save a dictionary-tree-library based database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. If a file is not desired to be
		saved, its path should be set to ''.
		"""

		# Save dictionary
		if dictstr != '':
			f = open(dictstr, 'w')
			f.write('////////////////////////////////////////////////////////////////////////////////\n')
			f.write('//\n')
			f.write('//\tDictionary\n')
			f.write('//\n')
			f.write('////////////////////////////////////////////////////////////////////////////////\n')
			f.write('\n')
			for label, struct in self.dictionary.iteritems():
				f.write(label + '\n')
				if isinstance(struct, structure.Structure):
					f.write(struct.toAdjacencyList() + '\n')
				elif struct == 'union':
					union = 'Union {'
					children = [child for child in self.tree.children[label]]
					union += ','.join(children)
					union += '}'
					f.write(union + '\n\n')
				else:
					raise InvalidDatabaseException('Unexpected item with label %s encountered in dictionary while attempting to save.' % label)
			f.close()

		# Save tree
		if treestr != '':
			logging.warning('Tree saving not yet implemented.')

		# Save library
		if libstr != '':
			logging.warning('Library saving not yet implemented.')


	def toXML(self, dom, root):
		"""
		Return an XML representation of the database.
		"""
		self.dictionary.toXML(dom, root)
		self.tree.toXML(dom, root)
		self.library.toXML(dom, root)

	def isWellFormed(self):
		"""
		Return :data:`True` if the database is well-formed. A well-formed
		database has an entry in the dictionary for every entry in the tree, and
		an entry in the tree for every entry in the library. If no tree is
		present (e.g. the primary libraries), then every entry in the library
		must have an entry in the dictionary. Finally, each entry in the
		library must have the same number of nodes as the number of top-level
		nodes in the tree, if the tree is present; this is for databases with
		multiple trees, e.g. the kinetics databases.
		"""

		wellFormed = True

		# Make list of all nodes in library
		libraryNodes = []
		for nodes in self.library:
			libraryNodes.extend(nodes.split(';'))
		libraryNodes = list(set(libraryNodes))
		
		
		for node in libraryNodes:
			
			# All nodes in library must be in dictionary
			try:
				if node not in self.dictionary:
					raise InvalidDatabaseException('Node "%s" in library is not present in dictionary.' % (node))
			except InvalidDatabaseException, e:
				wellFormed = False
				logging.error(e.msg)

			# If a tree is present, all nodes in library must be in tree
			if len(self.tree.parent) > 0:
				try:
					if node not in self.tree.parent:
						raise InvalidDatabaseException('Node "%s" in library is not present in tree.' % (node))
				except InvalidDatabaseException, e:
					wellFormed = False
					logging.error(e.msg)

		# If a tree is present, all nodes in tree must be in dictionary
		if self.tree is not None:
			for node in self.tree.parent:
				try:
					if node not in self.dictionary:
						raise InvalidDatabaseException('Node "%s" in tree is not present in dictionary.' % (node))
				except InvalidDatabaseException, e:
					wellFormed = False
					logging.error(e.msg)

		return wellFormed


	def matchNodeToStructure(self, node, structure, atoms):
		"""
		Return :data:`True` if the `structure` centered at `atom` matches the
		structure at `node` in the dictionary. The structure at `node` should
		have atoms with the appropriate labels because they are set on loading
		and never change. However, the atoms in `structure` may not have the
		correct labels, hence the `atoms` parameter. The `atoms` parameter may
		include extra labels, and so we only require that every labeled atom in
		the functional group represented by `node` has an equivalent labeled
		atom in `structure`.
		"""
	#	if node in ['CH2_triplet', 'Y_1centerbirad']:
	#		logging.debug("Trying to match structure with node %s"%node)
		group = self.dictionary[node]
		if isinstance(group, LogicNode):
			return group.matchToStructure(self, structure, atoms)
			
		#if group.__class__ == str or group.__class__ == unicode:
		#	if group.lower() == 'union':
		#		match = False
		#		for child in self.tree.children[node]:
		#			if self.matchNodeToStructure(child, structure, atoms):
		#				match = True
		#		return match
		#	else:
		#		return False
		else:
			# try to pair up labeled atoms
			centers = group.getLabeledAtoms()
			map12_0 = {}; map21_0 = {}
			for label in centers.keys():
				# Make sure the labels are in both group and structure.
				if not (label in centers and label in atoms):
					return False 
				center = centers[label]
				atom = atoms[label]
				# Make sure labels actually point to atoms.
				if center is None or atom is None:
					return False
				# Semantic check #1: atoms with same label are equivalent
				elif not atom.isSpecificCaseOf(center):
					return False
				# Semantic check #2: labeled atoms that share bond in the group (node)
				# also share equivalent (or more specific) bond in the structure
				for atom1, atom2 in map21_0.iteritems():
					if group.hasBond(center, atom1) and structure.hasBond(atom, atom2):
						bond1 = group.getBond(center, atom1)   # bond1 is group
						bond2 = structure.getBond(atom, atom2) # bond2 is structure
						if not bond2.isSpecificCaseOf(bond1):
							return False
					elif group.hasBond(center, atom1): # but structure doesn't 
						return False
					# but we don't mind if...
					elif structure.hasBond(atom, atom2): # but group doesn't
						logging.debug("We don't mind that structure %s has bond"+
							" but group %s doesn't"%(str(structure),node))
				# Passed semantic checks, so add to maps of already-matched atoms
				map21_0[center] = atom; map12_0[atom] = center
			# use mapped (labeled) atoms to try to match subgraph
			return structure.isSubgraphIsomorphic(group, map12_0, map21_0)

	def descendTree(self, structure, atoms, root=None):
		"""
		Descend the tree in search of the functional group node that best
		matches the local structure around `atoms` in `structure`.
		
		If root=None then uses the first matching top node.
		
		Returns None if there is no matching root.
		"""
		
		if root is None:
			for root in self.tree.top:
				if self.matchNodeToStructure(root, structure, atoms):
					break # We've found a matching root
			else: # didn't break - matched no top nodes
				return None
		elif not self.matchNodeToStructure(root, structure, atoms):
			return None
			
		next = []
		for child in self.tree.children[root]:
			if self.matchNodeToStructure(child, structure, atoms):
				next.append(child)

		if len(next) == 1:
			return self.descendTree(structure, atoms, next[0])
		elif len(next) == 0:
			return root
		else:
			#print structure.toAdjacencyList()
			#raise InvalidDatabaseException('For structure %s, a node %s with non-mutually-exclusive children %s was encountered in tree with top level nodes %s.' % (structure.getFormula(), root, next, self.tree.top))
			logging.warning('For %s, a node %s with overlapping children %s was encountered in tree with top level nodes %s.' % (structure, root, next, self.tree.top))
			return root
		
################################################################################

def removeCommentFromLine(line):
	"""
	Remove a C++/Java style comment from a line of text. This refers
	particularly to comments that begin with a double-slash '//' and continue
	to the end of the line.
	"""
	
	index = line.find('//')
	if index >= 0:
		line = line[0:index]
	return line

################################################################################

if __name__ == '__main__':
	pass
	