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
Contains classes describing chemical reactions:

* :class:`Reaction` - A general chemical reaction

* :class:`ReactionRecipe` - A set of actions to take when applying a reaction

* :class:`ReactionFamily` - A database of a general family of related reactions

* :class:`ReactionFamilySet` - A set of reaction families

"""

import math
import os
import os.path

import rmg.chem as chem
import rmg.data as data
import rmg.log as logging
import rmg.structure as structure
import rmg.species as species
import rmg.reaction as reaction

from model import *

################################################################################

class ReactionException(Exception):
	"""
	An base exception for reactions.
	Takes a reaction object, and optional message
	"""
	def __init__(self, reaction, message=''):
		self.reaction = reaction
		self.message = message

	def __str__(self):
		string = "Reaction: "+str(self.reaction) + '\n'
		string += "Reaction Family: "+str(self.reaction.family) + '\n'
		for reactant in self.reaction.reactants:
			string += reactant.toAdjacencyList() + '\n'
		for product in self.reaction.products:
			string += product.toAdjacencyList() + '\n'
		if self.message: string += "Message: "+self.message
		return string

class UndeterminableKineticsException(ReactionException):
	"""
	An exception raised when attempts to select appropriate kinetic parameters
	for a chemical reaction are unsuccessful.
	"""
	def __init__(self, reaction, message=''):
		new_message = 'Kinetics could not be determined. '+message
		ReactionException.__init__(self,reaction,new_message)

class InvalidActionException(Exception):
	"""
	An exception to be raised when an invalid action is encountered in a
	reaction recipe.
	"""

	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return self.msg

################################################################################

class ReactionRecipe:
	"""
	Represent a list of actions that, when executed, result in the conversion
	of a set of reactants to a set of products. There are currently five such
	actions:

	============  =======================  =====================================
	Action Name   Arguments                Action
	============  =======================  =====================================
	CHANGE_BOND   center1, order, center2  change the bond order of the bond between center1 and center2 by order; do not break or form bonds
	FORM_BOND     center1, order, center2  form a new bond between center1 and center2 of type order
	BREAK_BOND    center1, order, center2  break the bond between center1 and center2, which should be of type order
	GAIN_RADICAL  center, radical          increase the number of free electrons on center by radical
	LOSE_RADICAL  center, radical          decrease the number of free electrons on center by radical
	============  =======================  =====================================

	The actions are stored as a list in the `actions` attribute. Each action is
	a list of items; the first is the action name, while the rest are the
	action parameters as indicated above.
	"""

	def __init__(self, actions=None):
		self.actions = actions or []

	def addAction(self, action):
		"""
		Add an action to the reaction recipe.
		"""
		self.actions.append(action)

	def getReverse(self):
		"""
		Generate a reaction recipe that, when applied, does the opposite of
		what the current recipe does, i.e., it is the recipe for the reverse
		of the reaction that this is the recipe for.
		"""
		other = ReactionRecipe()
		for action in self.actions:
			if action[0] == 'CHANGE_BOND':
				other.addAction(['CHANGE_BOND', action[1], str(-int(action[2])), action[3]])
			elif action[0] == 'FORM_BOND':
				other.addAction(['BREAK_BOND', action[1], action[2], action[3]])
			elif action[0] == 'BREAK_BOND':
				other.addAction(['FORM_BOND', action[1], action[2], action[3]])
			elif action[0] == 'LOSE_RADICAL':
				other.addAction(['GAIN_RADICAL', action[1], action[2]])
			elif action[0] == 'GAIN_RADICAL':
				other.addAction(['LOSE_RADICAL', action[1], action[2]])
		return other

	def __apply(self, struct, doForward, unique):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. The `doForward` parameter is used to indicate
		whether the forward or reverse recipe should be applied. The atoms in
		the structure should be labeled with the appropriate atom centers.
		"""

		try:

			for action in self.actions:
				if action[0] == 'CHANGE_BOND' or action[0] == 'FORM_BOND' or action[0] == 'BREAK_BOND':

					label1, info, label2 = action[1:]

					# Find associated atoms
					atom1 = struct.getLabeledAtom(label1)
					atom2 = struct.getLabeledAtom(label2)
					if atom1 is None or atom2 is None or atom1 is atom2:
						raise InvalidActionException('Invalid atom labels encountered.')

					# Find associated bond, if present
					bond = struct.getBond(atom1, atom2)

					# If found, apply action
					if action[0] == 'CHANGE_BOND':
						if bond is None:
							raise InvalidActionException('Attempted to change the bond order of a nonexistent bond.')

						info = int(info)
						if abs(info) != 1:
							raise InvalidActionException('Attempted to change the bond order of a bond by more than one at a time.')

						if (info > 0 and doForward) or (info < 0 and not doForward):
							# Increment bond order
							bond.increaseOrder()
							atom1.incrementBond(struct.getBonds(atom1), unique)
							atom2.incrementBond(struct.getBonds(atom2), unique)
						else:
							# Decrement bond order
							bond.decreaseOrder()
							atom1.decrementBond(struct.getBonds(atom1), unique)
							atom2.decrementBond(struct.getBonds(atom2), unique)

					elif (action[0] == 'FORM_BOND' and doForward) or \
						(action[0] == 'BREAK_BOND' and not doForward):
						if bond is not None:
							raise InvalidActionException('Attempted to form a bond that already exists.')
						# Form single bond
						bond = chem.Bond([atom1, atom2], 'S')
						struct.addBond(bond)
						atom1.formBond(struct.getBonds(atom1), unique)
						atom2.formBond(struct.getBonds(atom2), unique)

						if doForward:
							bond = chem.Bond([atom1, atom2], info)
							struct.addBond(bond)
						else:
							bond = struct.getBond(atom1, atom2)
							if bond is None: raise InvalidActionException('Attempted to remove a nonexistent bond.')
							struct.removeBond(bond)

					elif (action[0] == 'BREAK_BOND' and doForward) or \
						(action[0] == 'FORM_BOND' and not doForward):
						if bond is None: raise InvalidActionException('Attempted to remove a nonexistent bond.')
						# Break single bond
						struct.removeBond(bond)
						atom1.breakBond(struct.getBonds(atom1), unique)
						atom2.breakBond(struct.getBonds(atom2), unique)

				elif action[0] == 'LOSE_RADICAL' or action[0] == 'GAIN_RADICAL':

					label, change = action[1:]

					# Find associated atoms
					atom = struct.getLabeledAtom(label)
					if atom is None:
						raise Exception('Invalid atom labels found while attempting to execute reaction recipe.')

					change = int(change)

					for i in range(change):
						if (action[0] == 'GAIN_RADICAL' and doForward) or \
							(action[0] == 'LOSE_RADICAL' and not doForward):
							# Increment radical count
							atom.increaseFreeElectron()
						elif (action[0] == 'LOSE_RADICAL' and doForward) or \
							(action[0] == 'GAIN_RADICAL' and not doForward):
							# Decrement radical count
							atom.decreaseFreeElectron()

				else:
					raise InvalidActionException('Unknown action "' + action[0] + '" encountered.')

		except InvalidActionException, e:
			logging.warning(e.msg)
			return False

		return True

	def applyForward(self, struct, unique=True):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures.
		"""
		return self.__apply(struct, True, unique)

	def applyReverse(self, struct, unique=True):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures.
		"""
		return self.__apply(struct, False, unique)

################################################################################

class ReactionFamily(data.Database):
	"""
	Represent a reaction family: a set of reactions with similar chemistry, and
	therefore similar reaction rates. Besides the dictionary, tree, and library
	inherited from :class:`data.Database`, the attributes are:

	===========  ===============================================================
	Attribute    Description
	===========  ===============================================================
	`label`      The name of the reaction family
	`template`   A :class:`Reaction` object representing the forward reaction
	             template
	`recipe`     A :class:`ReactionRecipe` object representing the steps to
	             take when applying the reaction to a set of reactants
	`forbidden`  (Optional) A dictionary of forbidden product structures
	`reverse`    A pointer to the reverse reaction family (or :data:`None` if
	             the family is its own reverse
	===========  ===============================================================

	"""

	def __init__(self, label='', template='', recipe=None):
		data.Database.__init__(self)
		# calling  data.Database.__init__(self) sets:
		# 	self.dictionary = Dictionary()
		# 	self.library = Library()
		# 	self.tree = Tree()
		self.label = label
		self.template = template
		self.recipe = recipe
		self.forbidden = None
		self.reverse = None
		self.is_reverse = False

	def __str__(self):
		return '<ReactionFamily(%s) from %s>'%(self.label,os.path.basename(self._path))

	def __str__(self):
		# append the path folder name to the reaction family label, and to the reverse
		return self.label + ' [%s]'%(os.path.basename(self._path))

	def getTemplateLists(self):
		"""
		Return lists containing the top-level nodes of each tree representing
		the reactants and the products. These lists are flattened versions of
		the lists available in self.template.reactants and
		self.template.products.
		"""
		forward = []
		for reactant in self.template.reactants:
			if type(reactant) is list:	forward.extend(reactant)
			else:						forward.append(reactant)
		reverse = []
		for product in self.template.products:
			if type(product) is list:	reverse.extend(product)
			else:						reverse.append(product)
		return forward, reverse

	def load(self, path):
		"""
		Load a reaction family located in the directory `path`. The family
		consists of the files::

			dictionary.txt
			tree.txt
			library.txt
			template.txt
			forbiddenGroups.txt

		"""
		import re
		# Generate paths to files in the database
		dictstr = os.path.join(path, 'dictionary.txt')
		treestr = os.path.join(path, 'tree.txt')
		libstr  = os.path.join(path, 'library.txt')
		tempstr = os.path.join(path, 'template.txt')
		forbstr = os.path.join(path, 'forbiddenGroups.txt')

		#: The path of the database that was loaded.
		self._path = path

		# Load the dictionary and tree using the generic methods
		# We can't use the generic method to load the library because it has
		# the type ('Arrhenius_EP') as the first meaningful line
		data.Database.load(self, dictstr, treestr, '')

		# Load the forbidden groups if the file 'forbiddenGroups.txt' is present
		# This file has the form of a standard dictionary so we can use the
		# standard dictionary loading function
		if os.path.exists(forbstr):
			self.forbidden = data.Dictionary()
			self.forbidden.load(forbstr)
			self.forbidden.toStructure()

		# Load the reaction template information and generate the reverse family
		# This requires that the dictionary and tree be loaded
		self.loadTemplate(tempstr)

		# Process the data in the library
		lines = self.library.load(libstr)

		# pop off the first line and check that it's 'Arrhenius_EP'
		line = lines.pop(0)
		if line != 'Arrhenius_EP':
			raise data.InvalidDatabaseException("Was expecting 'Arrhenius_EP' as first line, but got %s, in %s"%(line,libstr))

		#figure out how many labels there are
		test_line = lines[0].split()
		for token_no, token in enumerate(test_line):
			# skip token_no=0 because it's the label (but may match the regular expression)
			if token_no and re.match('^[0-9\-.]*$',token):
				# found the Temperature range at token_no
				number_of_groups = token_no-1
				logging.debug("Deduced there are %d groups %s in %s"%(number_of_groups,test_line[1:token_no],libstr))
				break
		else: # didn't break
			raise data.InvalidDatabaseException("Unable to figure out how many groups in %s using line %s"%(libstr,' '.join(test_line)))

		self.library.parse(lines, number_of_groups )
		self.processLibraryData()

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (path))

		# Fill in missing nodes in library via an averaging scheme of the
		# existing data; note that this disregards all temperature range
		# information
		forwardTemplate, reverseTemplate = self.getTemplateLists()

		#self.generateMissingEntriesFromBelow(forwardTemplate)
		#self.generateMissingEntriesFromAbove(forwardTemplate)

		# give the path to the reverse family too
		if self.reverse:
			self.reverse._path = self._path


	def prune(self, template=None):
		"""
		Remove nodes from the tree and dictionary that are not referred to in
		the library. Nodes that are not referred to in the library but that have
		one or more descendants that are referred to in the library are
		retained. The `template` parameter is a list of the nodes at which to
		begin, e.g. the template lists returned from :func:`getTemplateLists()`.
		"""

		if template is None:
			template, reverseTemplate = self.getTemplateLists()

		# Get lists of all nodes, sorted by top-level node
		nodeLists = []
		for parent in template:
			nodeList = []
			for node in self.tree.children:
				temp = node
				while temp is not None and temp not in template:
					temp = self.tree.parent[temp]
				if temp == parent:
					nodeList.append(node)
			nodeLists.append(nodeList)

		pruneList = []

		for i in range(len(template)):

			nodes = nodeLists[:]
			for node in nodes[i]:

				nodeList = [node]
				children = self.tree.children[node]
				nodeList.extend(children)

				while len(children) > 0:
					temp = children[:]
					children = []
					for child in temp:
						children.extend(self.tree.children[child])
					nodeList.extend(children)

				nodes[i] = nodeList

				nodesList = data.getAllCombinations(nodes)

				# If none of the possible combinations have data, the node
				# is a candidate for pruning
				candidate = True
				for nodeList in nodesList:
					if self.library.getData(nodeList) is not None:
						candidate = False

				# However, we will keep all top-level nodes because they are
				# needed for the template
				parent = self.tree.parent[node]
				if parent is None:
					candidate = False
				# We will also keep a node that is the child of a union in case
				# unions need the children explicitly defined in the tree
				elif isinstance(self.dictionary[parent],data.LogicNode):
						candidate = False
				
				if candidate:
					pruneList.append(node)

		# Complete pruning
		for node in pruneList:
			del self.dictionary[node]
			del self.tree.parent[node]
			del self.tree.children[node]

		# Also prune unneeded structures from dictionary
		pruneList = []
		for label, struct in self.dictionary.iteritems():
			if label not in self.tree.children:
				pruneList.append(label)
		for node in pruneList:
			del self.dictionary[node]

	def generateMissingEntriesFromBelow(self, nodes):
		"""
		Generate a nonexisting entry in the library based on an averaging
		scheme.
		"""

		# Get all possible combinations of child nodes
		children = []
		for node in nodes:
			nodeList = [node]; nodeList.extend(self.tree.children[node])
			children.append(nodeList)
		childNodesList = data.getAllCombinations(children)
		# Remove current nodes
		for i in range(childNodesList.count(nodes)):
			childNodesList.remove(nodes)

		# Recursively descend child nodes
		kinetics = []
		for childNodes in childNodesList:
			k = self.generateMissingEntriesFromBelow(childNodes)
			if k is not None: kinetics.append(k)

		# Only generate new entry if data does not exist or rank is zero;
		# otherwise return existing value
		if self.library.getData(nodes) is not None:
			if self.library.getData(nodes).rank > 0:
				return self.library.getData(nodes)

		# Otherwise average all of the kinetics parameters found above
		if len(kinetics) > 0:
			kin = self.averageKinetics(kinetics)
			kin.rank = 5
			self.library.add(nodes, kin)
			print nodes, kin


	def averageKinetics(self, kinetics):
		"""
		Return the average kinetic parameters for the list of kinetic data
		`kinetics`.
		"""
		if len(kinetics) == 0:
			return None

		# Use geometric average of parameters
		lnA = 0.0; E0 = 0.0; n = 0.0; alpha = 0.0
		for k in kinetics:
			lnA += math.log(k.A)
			E0 += k.E0
			n += k.n
			alpha += k.alpha

		lnA /= len(kinetics)
		E0 /= len(kinetics)
		n /= len(kinetics)
		alpha /= len(kinetics)

		kin = ArrheniusEPModel(math.exp(lnA), E0, n, alpha)
		kin.Trange = [0.0, 0.0]
		return kin

	def drawFullGraphOfTree(self):
		"""
		Create a PyDOT representation of the current tree.
		"""

		import pydot

		graph = pydot.Dot(size='10,8', page='10,8' ,  rankdir='LR',
				graph_type='digraph', simplify=True, fontsize=10,
				overlap='true', dpi='85',center="True")

		forwardTemplate, reverseTemplate = self.getTemplateLists()

		nodeLists = [[] for top in forwardTemplate]
		for node in self.tree.parent:
			ancestors = self.tree.ancestors(node)
			index = -1
			if len(ancestors) > 0:
				try:
					index = forwardTemplate.index(ancestors[-1])
				except ValueError:
					pass
			elif node in forwardTemplate:
				index = forwardTemplate.index(node)
			if index >= 0 and index < len(forwardTemplate):
				nodeLists[index].append(node)

		nodesList = data.getAllCombinations(nodeLists)

		# Create vertices of graph
		for nodes in nodesList:
			label = ';'.join(nodes)
			node = pydot.Node(label)
			if self.library.getData(nodes) is not None:
				node.set_style('filled')
				node.set_fillcolor('#000000FF')
				node.set_fontcolor('#FFFFFFFF')
			graph.add_node(node)

		# Create edges of graph
		for nodes in nodesList:
			label = ';'.join(nodes)
			for i, node in enumerate(nodes):
				parent = nodes[:]
				parent[i] = self.tree.parent[node]
				if None not in parent:
					parentLabel = ';'.join(parent)
					graph.add_edge(pydot.Edge(parentLabel,label))

		return graph


	def drawGraphOfTree(self, nodes):
		"""draw a graph of the tree"""
		import pydot, re
		g=pydot.Dot(size='10,8', page='10,8' ,  rankdir='LR',
				graph_type='digraph', simplify=True, fontsize=10,
				overlap='true', dpi='85',center="True")

		rates =  self.library.keys() # known reaction rates in library
		# add reactionrate nodes to graph
		for rate in rates:
			g.add_node(pydot.Node(rate))
		# add edges from each node to each of its ancestors
		for rate in rates:
			nodes=rate.split(';')
			for trialRate in rates: # the one we are testing for ancestry
				if trialRate == rate: continue # reject if itself
				for i,node in enumerate(nodes):
					trialNode = trialRate.split(';')[i]
					if trialNode == node: continue # ok if equal
					if trialNode not in self.tree.ancestors(node): break
					# if stopped through break, then it will not run the else clause
				else:
					# loop fell through all nodes without breaking:
					# trialRate must be ancestor of rate
					g.add_edge(pydot.Edge(trialRate,rate))
		g.set('fontsize','10')
		format='svg'
		prog='dot'
		f=open(self.label+'.dot','w')
		f.write(g.to_string())
		f.close()
		filename=self.label+'.'+format
		if format=='svg':  # annoyingly, dot creates svg's without units on the font size attribute.
			st=g.create_svg(prog=prog)
			st=re.sub(r"(font\-size\:[0-9]+\.*[0-9]*)([^p])",r"\1pt\2",st)
			f=open(filename,'w')
			f.write(st)
			f.close()
		else:
			g.write(filename,format=format,prog=prog)



	def generateMissingEntriesFromAbove(self, nodes):
		"""
		Generate a nonexisting entry in the library based on an averaging
		scheme.
		"""

		# Generate list of all sets of nodes that should have entries in the
		# library
		nodeLists = []
		for parent in nodes:
			nodeList = []
			for node in self.tree.children:
				temp = node
				while temp is not None and temp not in nodes:
					temp = self.tree.parent[temp]
				if temp == parent:
					nodeList.append(node)
			nodeLists.append(nodeList)
		nodesList = data.getAllCombinations(nodeLists)

		data = []
		for nodeList in nodesList:
			k = self.generateMissingEntryFromAbove(nodeList)
			if k is not None:
				data.append((nodeList, k))

		for nodeList, kinetics in data:
			self.library.add(nodeList, kinetics)
			#print nodeList, kinetics

	def generateMissingEntryFromAbove(self, nodes):

		# If an entry is already present, return it
		if self.library.getData(nodes) is not None:
			return self.library.getData(nodes)

		# Generate list of parents
		parentNodeLists = []
		for node in nodes:
			node0 = node
			parentNodeList = []
			while node0 is not None:
				parentNodeList.append(node0)
				node0 = self.tree.parent[node0]
			parentNodeLists.append(parentNodeList)
		parentNodesList = data.getAllCombinations(parentNodeLists)

		# Generate list of existing kinetic parameters along with "distance"
		# values (smaller is better)
		# This assumes that the amount of specificity indicated by moving from
		# parent to child in each tree is approximately equal
		minDistance = 1000000
		kinetics = []
		for parentNodes in parentNodesList:

			# Calculate distance
			distance = 0
			for i, node in enumerate(parentNodes):
				distance += parentNodeLists[i].index(node)

			# Don't bother if we've already found kinetics for a smaller
			# distance than the current distance, as we're just going to ignore
			# it anyway
			if distance < minDistance:
				# Get kinetics and append if available
				k = self.library.getData(parentNodes)
				if k is not None:
					kinetics.append((distance, k))
					if distance < minDistance: minDistance = distance

		# Fail if no kinetics found
		if len(kinetics) == 0:
			logging.warning('Unable to estimate missing kinetics for nodes %s in reaction family %s.' % (nodes, self.label))
			return None

		# Prune all entries with distances higher than the minimum
		kinetics = [k for d, k in kinetics if d == minDistance]
		if len(kinetics) == 0:
			logging.warning('Unable to estimate missing kinetics for nodes %s in reaction family %s.' % (nodes, self.label))
			return None

		# Average remaining kinetics
		kin = self.averageKinetics(kinetics)
		return kin

	def processLibraryData(self):
		"""
		Convert the data in the library from a string/unicode object to either
		an :class:`ArrheniusEPModel` object or a list of [link, comment]
		string pairs. This function is generally called in the course of
		loading a database from files.
		"""

		for label, item in self.library.iteritems():

			if item is None:
				pass
			elif not item.__class__ is tuple:
				raise data.InvalidDatabaseException('Kinetics library should be tuple at this point. Instead got %r'%data)
			else:
				index,data = item # break apart tuple, recover the 'index' - the beginning of the line in the library file.
				# Is't it dangerous having a local variable with the same name as a module?
				# what if we want to raise another data.InvalidDatabaseException() ?
				if not ( data.__class__ is str or data.__class__ is unicode) :
					raise data.InvalidDatabaseException('Kinetics library data format is unrecognized.')

				items = data.split()
				try:
					kineticData = [];
					# First item is temperature range
					kineticData.extend(items[0].split('-'))
					if len(kineticData) == 2:
						kineticData[0] = float(kineticData[0])
						kineticData[1] = float(kineticData[1])
					elif len(kineticData) == 1:
						kineticData = [float(items[0]), float(items[0])]
					# Middle items are Arrhenius + Evans-Polanyi data
					for i in range(1, 5):
						kineticData.append(float(items[i]))
					for i in range(5, 9):
						# Convert multiplicative uncertainties to additive
						# uncertainties if needed
						if items[i][0] == '*':
							kineticData.append((float(items[i][1:]) - 1.0) * float(items[i-4]))
						else:
							kineticData.append(float(items[i]))
					# Final item before comment is quality
					kineticData.append(int(items[9]))
					# Everything else is a comment
					comment = ' '.join(items[10:])

					kinetics = ArrheniusEPModel()
					kinetics.fromDatabase(kineticData, comment, len(self.template.reactants))
					kinetics.family = self
					kinetics.label = label
					kinetics.index = index
					#kinetics.comment = self.label + ' ' + label + ' ' + kinetics.comment
					self.library[label] = kinetics

				except (ValueError, IndexError), e:
					# Split data into link string and comment string; store
					# as list of length 2
					link = items[0]
					comment = data[len(link)+1:].strip()
					self.library[label] = [link, comment]

	def loadTemplate(self, path):
		"""
		Load and process a reaction template file located at `path`. This file
		is part of every reaction family.
		"""

		# Process the template file, removing comments and empty lines
		info = ''
		frec = None
		try:
			frec = open(path, 'r')
			for line in frec:
				line = data.removeCommentFromLine(line).strip()
				if len(line) > 0:
					info += line + '\n'
		except data.InvalidDatabaseException, e:
			logging.exception(str(e))
			return
		except IOError, e:
			logging.exception('Database file "' + e.filename + '" not found.')
			raise
		finally:
			if frec: frec.close()

		lines = info.splitlines()

		# First line is 'Forward: <name of forward reaction>
		# Second line is 'Reverse: <name of reverse reaction>
		forward = ''; reverse = ''
		if lines[0].find('Forward:') > -1:
			self.label = lines[0][9:].strip()
		if lines[1].find('Reverse:') > -1:
			reverse = lines[1][9:].strip()

		# Third line is reaction template, of the form
		# A1 A2 ... + B1 B2 ... + ... <---> C1 C2 ... + D1 D2 ... + ...
		# A, B, ... are reactants; C, D, ... are products
		# A1, A2, ... represent different trees for the same species
		# The first tree of each species is always used to identify the
		# reactants, so it should have all of the labeled atoms that are in that
		# reactant
		# The other trees can be used to provide functional group trees for
		# different parts of the molecule
		reactants = []; products = []; species = []
		atArrow = False
		items = lines[2].split(); items.extend('+')
		for item in items:
			if item == '+' or (item[0] == '<' and item[-1] == '>' and item.find('-') > -1):
				if len(species) == 1: species = species[0]
				if atArrow:		products.append(species)
				else:			reactants.append(species)
				species = []
				if item[0] == '<' and item[-1] == '>' and item.find('-') > -1:
					atArrow = True
			else:
				# Check that all reactant structures are in dictionary
				# The product structures are generated automatically and need
				# not be included
				if item not in self.dictionary and not atArrow:
					raise data.InvalidDatabaseException('Reaction family template contains an unknown structure.')
				species.append(item)

		# Set template reaction
		self.template = reaction.Reaction(reactants=reactants, products=products)

		# Remaining lines are reaction recipe for forward reaction
		self.recipe = ReactionRecipe()
		for line in lines[3:]:
			line = line.strip()

			# First item is the name of the action
			items = line.split()
			action = [ items[0].upper() ]
			if items[0] != 'CHANGE_BOND' and items[0] != 'FORM_BOND' and \
			  items[0] != 'BREAK_BOND' and items[0] != 'GAIN_RADICAL' and \
			  items[0] != 'LOSE_RADICAL':
				print items[0]
			assert action[0] in ['CHANGE_BOND','FORM_BOND','BREAK_BOND','GAIN_RADICAL','LOSE_RADICAL']

			# Remaining items are comma-delimited list of parameters enclosed by
			# {}, which we will split into individual parameters
			action.extend(line[len(items[0]):].strip()[1:-1].split(','))

			self.recipe.addAction(action)

		# Generate the reverse template
		if reverse != self.label:
			template = reaction.Reaction(reactants=self.template.products, products=self.template.reactants)
			self.reverse = ReactionFamily(reverse, template, self.recipe.getReverse())
			self.reverse.dictionary = self.dictionary
			self.reverse.tree = self.tree
			self.reverse.library = data.Library()
			self.reverse.forbidden = self.forbidden
			self.reverse.reverse = self
			self.reverse.is_reverse = True

		# If necessary, generate the product template structure(s)
		# Don't need to do this if family is its own reverse
		if reverse != self.label:
			self.generateProductTemplate()

	def generateProductTemplate(self):
		"""
		Generate the product structures by applying the reaction template to
		the top-level nodes. For reactants defined by multiple structures, only
		the first is used here; it is assumed to be the most generic.
		"""

		# First, generate a list of reactant structures that are actual
		# structures, rather than unions
		reactantStructures = []
		
		logging.debug("Generating template for products.")
		for reactant in self.template.reactants:
			if isinstance(reactant, list):	reactants = [reactant[0]]
			else:							reactants = [reactant]
			
			logging.debug("Reactants:%s"%reactants)
			for s in reactants: #
				struct = self.dictionary[s]
				#logging.debug("Reactant %s is class %s"%(str(s),struct.__class__))
				if isinstance(struct, data.LogicNode):
					all_structures = struct.getPossibleStructures(self.dictionary) 
					logging.debug('Expanding node %s to %s'%(s, all_structures))
					reactantStructures.append(all_structures)
				else:
					reactantStructures.append([struct])
		
		# Second, get all possible combinations of reactant structures
		reactantStructures = data.getAllCombinations(reactantStructures)

		# Third, generate all possible product structures by applying the
		# recipe to each combination of reactant structures
		# Note that bimolecular products are split by labeled atoms
		productStructures = []
		for reactantStructure in reactantStructures:
			productStructure = self.applyRecipe(reactantStructure, unique=False)
			productStructures.append(productStructure)

		# Fourth, remove duplicates from the lists
		productStructureList = [[] for i in range(len(productStructures[0]))]
		for productStructure in productStructures:
			for i, struct in enumerate(productStructure):
				found = False
				for s in productStructureList[i]:
					if s.isIsomorphic(struct): found = True
				if not found:
					productStructureList[i].append(struct)

		# Fifth, associate structures with product template
		for i in range(len(self.template.products)):
			if len(productStructureList[i]) == 1:
				self.dictionary[self.template.products[i]] = productStructureList[i][0]
				self.tree.parent[self.template.products[i]] = None
				self.tree.children[self.template.products[i]] = []
			else:
				children = []
				for j in range(len(productStructureList[i])):
					label = '%s_%i' % (self.template.products[i], j+1)
					self.dictionary[label] = productStructureList[i][j]
					children.append(label)
					self.tree.parent[label] = self.template.products[i]
					self.tree.children[label] = []

				self.tree.parent[self.template.products[i]] = None
				self.tree.children[self.template.products[i]] = children
				
				self.dictionary[self.template.products[i]] = data.LogicOr(children,invert=False)

	def reactantMatch(self, reactant, templateReactant):
		"""
		Return :data:`True` if the provided reactant matches the provided
		template reactant and :data:`False` if not.
		
		Also returns complete lists of mappings
		"""
		maps12 = []; maps21 = []
		if templateReactant.__class__ == list: templateReactant = templateReactant[0]
		struct = self.dictionary[templateReactant]
		
		if isinstance(struct, data.LogicNode):
			for child_structure in struct.getPossibleStructures(self.dictionary):
				ismatch, map21, map12 = reactant.findSubgraphIsomorphisms(child_structure)
				if ismatch:
					maps12.extend(map12); maps21.extend(map21)
			return len(maps12) > 0, maps21, maps12
			
		elif struct.__class__ == structure.Structure:
			return reactant.findSubgraphIsomorphisms(struct)

	def applyRecipe(self, reactantStructures, unique=True):
		"""
		Apply the recipe for this reaction family to the list of
		:class:`structure.Structure` objects `reactantStructures`. The atoms
		of the reactant structures must already be tagged with the appropriate
		labels. Returns a list of structures corresponding to the products
		after checking that the correct number of products was produced.
		"""

		# There is some hardcoding of reaction families in this function, so
		# we need the label of the reaction family for this
		label = self.label.lower()

		# Merge reactant structures into single structure
		# Also copy structures so we don't modify the originals
		# Since the tagging has already occurred, both the reactants and the
		# products will have tags
		reactantStructure = structure.Structure()
		for s in reactantStructures:
			reactantStructure = reactantStructure.merge(s.copy())

		# Hardcoding of reaction family for radical recombination (colligation)
		# because the two reactants are identical, they have the same tags
		# In this case, we must change the labels from '*' and '*' to '*1' and
		# '*2'
		if label == 'colligation':
			identicalCenterCounter = 0
			for atom in reactantStructure.atoms():
				if atom.label == '*':
					identicalCenterCounter += 1
					atom.label = '*' + str(identicalCenterCounter)
			if identicalCenterCounter != 2:
				raise Exception('Unable to change labels from "*" to "*1" and "*2" for reaction family %s.' % (label))

		# Generate the product structure by applying the recipe
		if not self.recipe.applyForward(reactantStructure, unique):
			return None
		productStructure = reactantStructure

		# Hardcoding of reaction family for reverse of radical recombination
		# (Unimolecular homolysis)
		# Because the two products are identical, they should the same tags
		# In this case, we must change the labels from '*1' and '*2' to '*' and
		# '*'
		if label == 'unimolecular homolysis':
			for atom in productStructure.atoms():
				if atom.label == '*1' or atom.label == '*2': atom.label = '*'


		# If reaction family is its own reverse, relabel atoms
		if not self.reverse:
			# Get atom labels for products
			atomLabels = {}
			for atom in productStructure.atoms():
				if atom.label != '':
					atomLabels[atom.label] = atom

			# This is hardcoding of reaction families (bad!)
			label = self.label.lower()
			if label == 'h abstraction':
				# '*2' is the H that migrates
				# it moves from '*1' to '*3'
				atomLabels['*1'].label = '*3'
				atomLabels['*3'].label = '*1'

			elif label == 'intra h migration':
				# '*3' is the H that migrates
				# swap the two ends between which the H moves
				atomLabels['*1'].label = '*2'
				atomLabels['*2'].label = '*1'
				# reverse all the atoms in the chain between *1 and *2
				# i.e. swap *4 with the highest, *5 with the second-highest
				highest = len(atomLabels)
				if highest>4:
					for i in range(4,highest+1):
						atomLabels['*%d'%i].label = '*%d'%(4+highest-i)


		# Split product structure into multiple species if necessary
		if len(self.template.products) > 1:
			productStructures = productStructure.split()
		else:
			productStructures = [productStructure]

		# Make sure we've made the expected number of products
		if len(self.template.products) != len(productStructures):
			# We have a different number of products than expected by the template.
			# It might be because we found a ring-opening using a homolysis template
			if (label=='unimolecular homolysis'
			 and len(productStructures) == 1
			 and len(reactantStructures) == 1):
				# just be absolutely sure (maybe slow, but safe)
				rs = reactantStructures[0]
				if ( rs.graph.isVertexInCycle(rs.getLabeledAtom('*1'))
				 and rs.graph.isVertexInCycle(rs.getLabeledAtom('*2'))):
					# both *1 and *2 are in cycles (probably the same one)
					# so it's pretty safe to just fail quietly,
					# and try the next reaction
					return None

			# no other excuses, raise an exception
			message = 'Application of reaction recipe failed; expected %s product(s), but %s found.\n' % (len(self.template.products), len(productStructures))
			message += "Reaction family: %s \n"%str(self)
			message += "Reactant structures: %s \n"%reactantStructures
			message += "Product structures: %s \n"%productStructures
			message += "Template: %s"%self.template
			logging.error(message)
			return None # don't fail!!! muhahaha
			raise Exception(message)

		# If there are two product structures, place the one containing '*1' first
		if len(productStructures) == 2:
			if not productStructures[0].containsLabeledAtom('*1') and \
				productStructures[1].containsLabeledAtom('*1'):
				productStructures.reverse()

		# reset any cached structure information because it is now invalid
		for struct in productStructures:
			struct.resetCachedStructureInfo()

		# Return the product structures
		return productStructures

	def generateProductStructures(self, reactantStructures, maps):
		"""
		For a given set of `reactantStructures` and a given set of `maps`,
		generate and return the corresponding product structures. The
		`reactantStructures` parameter should be given in the order the
		reactants are stored in the reaction family template. The `maps`
		parameter is a list of mappings of the top-level tree node of each
		*template* reactant to the corresponding *structure*. This function
		returns the product structures, species, and a boolean that tells
		whether any species are new.
		"""

		# Clear any previous atom labeling from all reactant structures
		for struct in reactantStructures: struct.clearLabeledAtoms()

		# If there are two structures and they are the same, then make a copy
		# of the second one and adjust the second map to point to its atoms
		# This is for the case where A + A --> products
		if len(reactantStructures) == 2 and reactantStructures[0] == reactantStructures[1]:
			reactantStructures[1], newMap = reactantStructures[1].copy(returnMap=True)
			maps[1] = dict([(templateAtom,newMap[reactantAtom]) for templateAtom, reactantAtom in maps[1].iteritems()])

		# Tag atoms with labels
		for map in maps:
			for templateAtom, reactantAtom in map.iteritems():
				reactantAtom.label = templateAtom.label

		# Generate the product structures by applying the forward reaction recipe
		try:
			productStructures = self.applyRecipe(reactantStructures)
			if not productStructures: return None
		except chem.InvalidChemicalActionException, e:
			print 'Unable to apply reaction recipe!'
			print 'Reaction family is %s' % self
			print 'Reactant structures are:'
			for struct in reactantStructures:
				print struct.toAdjacencyList()
			raise

		# Check that reactant and product structures are allowed in this family
		# If not, then stop
		if self.forbidden is not None:
			for label, struct2 in self.forbidden.iteritems():
				for struct in reactantStructures:
					if struct.isSubgraphIsomorphic(struct2): return None
				for struct in productStructures:
					if struct.isSubgraphIsomorphic(struct2): return None

		return productStructures

	def createReaction(self, reactants, reactantStructures, productStructures, reactantAtomLabels):
		"""
		Create a :class:`Reaction` object representing the reaction that
		converts the structures in `reactantStructures` corresponding to the
		species in `reactants` to the structures in `productStructures`. The
		atom labels for the reactants should already be known, and they are
		passed in the `reactantAtomLabels` parameter.
		"""

		productAtomLabels = []
		for struct in productStructures:
			productAtomLabels.append(struct.getLabeledAtoms())
		
		# Convert product structures to product species
		products = []
		for i, struct0 in enumerate(productStructures):
			found, spec, struct, map = species.checkForExistingSpecies(struct0)
			if found:
				# Adjust atom labels mapping accordingly
				for label, atom in productAtomLabels[i].iteritems():
					productAtomLabels[i][label] = map[atom]
				# Save struct rather than struct0
				productStructures[i] = struct
				# Append product species to list of products
				products.append(spec)
			else:
				product, isNew = species.makeNewSpecies(struct0, checkExisting=False)
				# Don't make a new reaction if no species was returned from
				# makeNewSpecies() (e.g. due to forbidden structure)
				if product is None: return None
				products.append(product)

		# Sort reactants and products (to make comparisons easier/faster)
		reactants.sort()
		products.sort()

		# Check that the reaction actually results in a different set of species
		if reactants == products:
			return None

		# Create reaction object
		forward = reaction.Reaction(reactants=reactants, products=products, family=self)
		reverse = reaction.Reaction(reactants=products, products=reactants, family=self.reverse or self)
		forward.reverse = reverse
		reverse.reverse = forward

		# Dictionaries containing the labeled atoms for the reactants and products
		forward.atomLabels = reactantAtomLabels
		reverse.atomLabels = productAtomLabels

		# Return the created reaction (forward direction only)
		return forward

	def getReactionList(self, reactants):
		"""
		Generate a list of all of the possible reactions of this family between
		the list of `reactants`. The number of reactants provided must match
		the number of reactants expected by the template, or this function
		will return an empty list.
		"""

		rxnList = []

		# Unimolecular reactants: A --> products
		if len(reactants) == 1 and self.template.isUnimolecular():

			# Iterate over all resonance isomers of the reactant
			for structure in reactants[0].structure:

				ismatch, map21, map12 = self.reactantMatch(structure, self.template.reactants[0])
				if ismatch:
					for map in map12:

						reactantAtomLabels = [{}]
						for atom1, atom2 in map.iteritems():
							reactantAtomLabels[0][atom1.label] = atom2

						reactantStructures = [structure]
						productStructures = self.generateProductStructures(reactantStructures, [map])
						if productStructures:
							rxn = self.createReaction(reactants, reactantStructures, productStructures, reactantAtomLabels)
							if rxn: rxnList.append(rxn)

		# Bimolecular reactants: A + B --> products
		elif len(reactants) == 2 and self.template.isBimolecular():

			structuresA = reactants[0].structure
			structuresB = reactants[1].structure

			# Iterate over all resonance isomers of the reactant
			for structureA in structuresA:
				for structureB in structuresB:

					# Reactants stored as A + B
					ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[0])
					ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[1])

					# Iterate over each pair of matches (A, B)
					if ismatch_A and ismatch_B:
						for mapA in map12_A:
							for mapB in map12_B:

								reactantAtomLabels = [{},{}]
								for atom1, atom2 in mapA.iteritems():
									reactantAtomLabels[0][atom1.label] = atom2
								for atom1, atom2 in mapB.iteritems():
									reactantAtomLabels[1][atom1.label] = atom2

								reactantStructures = [structureA, structureB]
								productStructures = self.generateProductStructures(reactantStructures, [mapA, mapB])
								if productStructures:
									rxn = self.createReaction(reactants, reactantStructures, productStructures, reactantAtomLabels)
									if rxn: rxnList.append(rxn)

					# Only check for swapped reactants if they are different
					if reactants[0].id != reactants[1].id:

						# Reactants stored as B + A
						ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[1])
						ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[0])

						# Iterate over each pair of matches (A, B)
						if ismatch_A and ismatch_B:
							for mapA in map12_A:
								for mapB in map12_B:

									reactantAtomLabels = [{},{}]
									for atom1, atom2 in mapA.iteritems():
										reactantAtomLabels[0][atom1.label] = atom2
									for atom1, atom2 in mapB.iteritems():
										reactantAtomLabels[1][atom1.label] = atom2

									reactantStructures = [structureA, structureB]
									productStructures = self.generateProductStructures(reactantStructures, [mapA, mapB])
									if productStructures:
										rxn = self.createReaction(reactants, reactantStructures, productStructures, reactantAtomLabels)
										if rxn: rxnList.append(rxn)

		# Merge duplicate reactions and increment multiplier
		# In this context we already know that the family and the reactants
		# match, so we only need to check the products
		reactionsToRemove = []
		for i, rxn1 in enumerate(rxnList):
			for j, rxn2 in enumerate(rxnList[i+1:]):
				if rxn2 not in reactionsToRemove:
					if rxn1.products == rxn2.products:
						reactionsToRemove.append(rxn2)
						rxn1.multiplier += 1.0
		for rxn in reactionsToRemove:
			rxnList.remove(rxn)

		# For R_Recombination reactions, the multiplier is twice what it should
		# be, so divide those by two
		# This is hardcoding of reaction families!
		if self.label.lower() == 'unimolecular homolysis':
			for rxn in rxnList:
				rxn.multiplier /= 2

		# Formally make the new reactions
		reactionsToRemove = []
		for i in range(len(rxnList)):
			rxn, isNew = reaction.makeNewReaction(rxnList[i])
			if isNew:
				rxnList[i] = rxn
			else:
				reactionsToRemove.append(rxnList[i])
		for rxn in reactionsToRemove:
			rxnList.remove(rxn)

		return rxnList

	def getKinetics(self, reaction, structures):
		"""
		Determine the appropriate kinetics for `reaction` which involves the
		labeled atoms in `atoms`.
		"""

		# Get forward reaction template and remove any duplicates
		forwardTemplate, reverseTemplate = self.getTemplateLists()
		#forwardTemplate = list(set(forwardTemplate)) # this can shuffle the order!
		temporary=[]
		symmetric_tree=False
		for node in forwardTemplate:
			if node not in temporary:
				temporary.append(node)
			else:
				# duplicate node found at top of tree
				# eg. R_recombination: ['Y_rad', 'Y_rad']
				assert len(forwardTemplate)==2 , 'Can currently only do symmetric trees with nothing else in them'
				symmetric_tree = True
		forwardTemplate = temporary

		# Descend reactant trees as far as possible
		template = []
		for forward in forwardTemplate:
			# 'forward' is a head node that should be matched.
			# Get labeled atoms of forward
			
			#if forward=='Rn':
			#	import pdb
			#	pdb.set_trace()
			node = forward
			group = self.dictionary[node]
			# to sort out "union" groups:
			# descends to the first child that's not a logical node
			if isinstance(group,data.LogicNode):
				all_structures = group.getPossibleStructures(self.dictionary)
				group = all_structures[0]
				node = 'First sub-structure of '+node
			# ...but this child may not match the structure.
			# eg. an R3 ring node will not match an R4 ring structure.
			# (but at least the first such child will contain fewest labels - we hope)

			atomList = group.getLabeledAtoms() # list of atom labels in highest non-union node

			for struct in structures:
				# Match labeled atoms
				# Check this structure has each of the atom labels in this group
				has_all_atom_labels = True
				for label in atomList:
					if not struct.containsLabeledAtom(label):
						has_all_atom_labels = False
						#print "%s does not contain %s atom labels"%(struct.toSMILES(),node)
						break
						# this structure (reactant/product) does not contain this group's atom labels
				if not has_all_atom_labels: continue # don't try to match this structure - the atoms aren't there!

				# Match structures
				atoms = struct.getLabeledAtoms()
				matched_node = self.descendTree(struct, atoms, forward)
				if matched_node is not None:
					template.append(matched_node)
				else:
					logging.warning("Couldn't find match for %s in %s"%(forward,atomList))
					logging.warning( struct.toAdjacencyList() )

		# Get fresh templates (with duplicate nodes back in)
		forwardTemplate, reverseTemplate = self.getTemplateLists()

		# Check that we were able to match the template.
		# template is a list of the actual matched nodes
		# forwardTemplate is a list of the top level nodes that should be matched
		if len(template) != len(forwardTemplate):
			logging.warning('Warning: Unable to find matching template for reaction %s in reaction family %s' % (str(reaction), str(self)) )
			logging.warning(" Trying to match " + str(forwardTemplate))
			logging.warning(" Matched "+str(template))
			raise UndeterminableKineticsException(reaction)
			print str(self), template, forwardTemplate, reverseTemplate
			for reactant in reaction.reactants:
				print reactant.toAdjacencyList() + '\n'
			for product in reaction.products:
				print product.toAdjacencyList() + '\n'

			## If unable to match template, use the most general template
			#template = forwardTemplate


#		k = self.library.getData(template)
#		print template, k
#		if k is not None: return [k]
#		else: return None


		# climb the tree finding ancestors
		nodeLists = []
		for temp in template:
			nodeList = []
			while temp is not None:
				nodeList.append(temp)
				temp = self.tree.parent[temp]
			nodeLists.append(nodeList)

		# Generate all possible combinations of nodes
		items = data.getAllCombinations(nodeLists)

		# Generate list of kinetics at every node
		#logging.debug("   Template contains %s"%forwardTemplate)
		kinetics = []
		for item in items:
			itemData = self.library.getData(item)
			#logging.debug("   Looking for %s found %r"%(item, itemData))
			if itemData is not None:
				kinetics.append(itemData)

			if symmetric_tree: # we might only store kinetics the other way around
				item.reverse()
				itemData = self.library.getData(item)
				#logging.debug("   Also looking for %s found %r"%(item, itemData))
				if itemData is not None:
					kinetics.append(itemData)

		if len(kinetics) == 0: return None

		return kinetics

################################################################################

class ReactionFamilySet:
	"""
	Represent a set of reaction families. The `families` attribute stores a
	dictionary of :class:`ReactionFamily` objects representing the families in
	the set.
	"""

	def __init__(self):
		self.families = {}

	def load(self, datapath, only_families=False):
		"""
		Load a set of reaction families from the general database
		specified at `datapath`. If only_families is present, families not in
		this list will not be loaded (e.g. only_families=['H_Abstraction'] )
		"""

		datapath = os.path.abspath(datapath)

		logging.info('Loading reaction family databases from %s...' % datapath)

		# Load the families from kinetics/families.txt
		familyList = []
		ffam = None
		try:
			ffam = open(os.path.join(datapath,'kinetics_groups','families.txt'), 'r')
			for line in ffam:
				line = data.removeCommentFromLine(line).strip()
				if len(line) > 0:
					items = line.split()
					items[0] = int(items[0])
					familyList.append(items)
		except data.InvalidDatabaseException, e:
			logging.exception(str(e))
			return
		except IOError, e:
			logging.exception('Database file "' + e.filename + '" not found.')
			raise
		finally:
			if ffam: ffam.close()

		# Load the reaction families (if they exist and status is 'on')
		self.families = {}
		for index, status, label in familyList:
			path = os.path.join(datapath, 'kinetics_groups', label)
			if os.path.isdir(path) and status.lower() == 'on':
				# skip families not in only_families, if it's set
				if only_families and label not in only_families: continue

				logging.info('Loading reaction family %s from %s...' % (label, datapath))
				family = ReactionFamily(label)
				family.load(path)
				self.families[family.label] = family
				if family.reverse is not None:
					self.families[family.reverse.label] = family.reverse
			else: logging.info("NOT loading family %s."%label)
			
		# initialize the global reactionDict with family names
		for key in self.families.values():
			reaction.reactionDict[key] = dict()

	def getReactions(self, species):
		"""
		Generate a list of reactions that involve a list of one or two `species`
		as a reactant or product.
		"""

		rxnList = []

		# Don't bother if any or all of the species are marked as nonreactive
		if not all([spec.reactive for spec in species]):
			return rxnList

		log_text = ' + '.join([str(spec) for spec in species])

		logging.info('Looking for reactions of %s'%(log_text))

		for key, family in self.families.iteritems():
			rxnList.extend(family.getReactionList(species))

		if len(rxnList) == 1:
			logging.info('Found %s new reaction for %s'%(len(rxnList), log_text))
		else:
			logging.info('Found %s new reactions for %s'%(len(rxnList), log_text))

		return rxnList

################################################################################

kineticsDatabase = None

def loadKineticsDatabase(dstr):
	"""
	Load the RMG kinetics database located at `dstr` into the global variable
	`rmg.reaction.kineticsDatabase`.
	"""
	global kineticsDatabase

	kineticsDatabase = ReactionFamilySet()
	kineticsDatabase.load(dstr)

################################################################################

