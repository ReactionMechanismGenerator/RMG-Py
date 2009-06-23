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
Contains classes describing chemical reactions.
"""

import quantities as pq
import logging
import os
import os.path

import chem
import data
import species

################################################################################

class ArrheniusEPKinetics:
	"""
	Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
	kinetic expression has the form

	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_\mathrm{a}}{RT} \\right)

	The parameter :math:`\\alpha` is used to correct the activation energy
	:math:`E_\\mathrm{a}` via the Evans-Polanyi formula

	.. math:: E_\\mathrm{a} = E_\\mathrm{a}^0 - (\\alpha - 1) \\Delta H_\\mathrm{rxn}

	"""

	def __init__(self, A=0.0, Ea=0.0, n=0.0, alpha=0.0):
		self.A = A
		self.Ea = Ea
		self.n = n
		self.alpha = alpha

	def fromDatabase(self, data, comment):
		"""
		Process a list of numbers `data` and associated description `comment`
		generated while reading from a kinetics database.
		"""

		if len(data) != 11:
			raise Exception('Invalid list of kinetic data; should be a list of numbers of length 10.')

		Tmin, Tmax, A, n, alpha, Ea, dA, dn, dalpha, dEa, rank = data

		self.Trange = pq.Quantity([Tmin, Tmax], 'K')

		self.A = pq.UncertainQuantity(A, 's^-1', dA)
		self.Ea = pq.UncertainQuantity(Ea, 'kcal/mol', dEa)
		self.Ea.units = 'J/mol'
		self.n = pq.UncertainQuantity(n, '', dn)
		self.alpha = pq.UncertainQuantity(alpha, 's^-1', dalpha)

		self.rank = rank
		self.comment = comment

################################################################################

class ReactionRecipe:
	"""
	Represent a list of actions that, when executed, result in the conversion
	of a set of reactants to a set of products. There are currently five such
	actions:

	- CHANGE_BOND {center1,order,center2} - change the bond order of the bond
	between center1 and center2 by order; do not break or form bonds

	- FORM_BOND {center1,order,center2} - form a new bond between center1 and
	center2 of type order

	- BREAK_BOND {center1,order,center2} - break the bond between center1 and
	center2, which should be of type order

	- GAIN_RADICAL {center,radical} - Increase the number of free electrons on
	center by radical

	- LOSE_RADICAL {center,radical} - Decrease the number of free electrons on
	center by radical

	Each action is a list of items; the first is the action name, while the
	rest are the action parameters as indicated above.
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

	def apply(self, structure, doForward):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. The `doForward` parameter is used to indicate
		whether the forward or reverse recipe should be applied. The atoms in
		the structure should be labeled with the appropriate atom centers.
		"""

		for action in self.actions:
			if action[0] == 'CHANGE_BOND' or action[0] == 'FORM_BOND' or action[0] == 'BREAK_BOND':

				label1, info, label2 = action[1:]

				# Find associated atoms
				atom1 = structure.getLabeledAtom(label1)
				atom2 = structure.getLabeledAtom(label2)
				if atom1 is None or atom2 is None or atom1 is atom2:
					raise Exception('Invalid atom labels found while attempting to execute reaction recipe.')

				# If found, change bond
				if action[0] == 'CHANGE_BOND':
					bond = structure.getBond(atom1, atom2)
					info = int(info)
					if bond is not None:
						for i in range(0, abs(info)):
							if doForward:
								if info > 0:	bond.increaseOrder()
								elif info < 0:	bond.decreaseOrder()
							else:
								if info > 0:	bond.decreaseOrder()
								elif info < 0:	bond.increaseOrder()
				elif action[0] == 'FORM_BOND':
					if doForward:
						bond = chem.Bond([atom1, atom2], info)
						structure.addBond(bond)
					else:
						bond = structure.getBond(atom1, atom2)
						structure.removeBond(bond)
				elif action[0] == 'BREAK_BOND':
					if doForward:
						bond = structure.getBond(atom1, atom2)
						structure.removeBond(bond)
					else:
						bond = chem.Bond([atom1, atom2], info)
						structure.addBond(bond)

			elif action[0] == 'LOSE_RADICAL' or action[0] == 'GAIN_RADICAL':

				label, change = action[1:]
				change = int(change)
				# Find associated atoms
				atom = structure.getLabeledAtom(label)
				if atom is None:
					raise Exception('Invalid atom labels found while attempting to execute reaction recipe.')

				# If found, adjust radical
				for i in range(0, change):
					if doForward:
						if action[0] == 'LOSE_RADICAL':		atom.decreaseFreeElectron()
						elif action[0] == 'GAIN_RADICAL':	atom.increaseFreeElectron()
					else:
						if action[0] == 'LOSE_RADICAL':		atom.increaseFreeElectron()
						elif action[0] == 'GAIN_RADICAL':	atom.decreaseFreeElectron()

			else:
				raise Exception('Unknown action "' + action[0] + '" encountered while attempting to execute reaction recipe.')


	def applyForward(self, structure):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures.
		"""
		return self.apply(structure, True)

	def applyReverse(self, structure):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. 
		"""
		return self.apply(structure, False)

################################################################################

class ReactionFamily(data.Database):
	"""
	Represent a reaction family: a set of reactions with similar chemistry, and
	therefore similar reaction rates. Attributes include the family name `name`,
	the reaction template `template`, a list `actions` of actions to take when
	reacting, and a dictionary-library-tree `database` of rate rules.
	"""

	def __init__(self, label='', template='', recipe=None):
		data.Database.__init__(self)
		self.label = label
		self.template = template
		self.recipe = recipe
		self.forbidden = None
		self.reverse = None

	def __str__(self):
		return self.label

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
		Load a reaction family located in the directory `path`.
		"""

		# Generate paths to files in the database
		dictstr = path + '/dictionary.txt'
		treestr = path + '/tree.txt'
		libstr = path + '/library.txt'
		tempstr = path + '/template.txt'
		forbstr = path + '/forbiddenGroups.txt'

		# Load the dictionary, tree, and library using the generic methods
		data.Database.load(self, dictstr, treestr, '')

		
		# Load the forbidden groups if necessary
		if os.path.exists(forbstr):
			self.forbidden = data.Dictionary()
			self.forbidden.load(forbstr)
			self.forbidden.toStructure()

		# Load the adjlist (must be last, as the reverse family is also generated)
		self.loadTemplate(tempstr)

		# Process the data in the library
		lines = self.library.load(libstr)
		self.library.parse(lines[2:], int(lines[1]))
		self.processLibraryData()

	def processLibraryData(self):
		"""
		Convert the data in the library from a string/unicode object to either
		an :class:`ArrheniusEPKinetics` object or a list of [link, comment]
		string pairs. This function is generally called in the course of
		loading a database from files.
		"""

		for label, data in self.library.iteritems():

			if data is None:
				pass
			elif data.__class__ == str or data.__class__ == unicode:
				items = data.split()
				try:
					kineticData = []; comment = ''
					# First item is temperature range
					kineticData.extend(items[0].split('-'))
					kineticData[0] = float(kineticData[0])
					kineticData[1] = float(kineticData[1])
					# Middle items are Arrhenius + Evans-Polanyi data
					for i in range(1, 9):
						kineticData.append(float(items[i]))
					# Final item before comment is quality
					kineticData.append(int(items[9]))
					# Everything else is a comment
					for i in range(10, len(items)):
						comment += items[i] + ' '

					kinetics = ArrheniusEPKinetics()
					kinetics.fromDatabase(kineticData, comment)
					self.library[label] = kinetics

				except (ValueError, IndexError), e:
					# Split data into link string and comment string; store
					# as list of length 2
					link = items[0]
					comment = data[len(link)+1:].strip()
					self.library[label] = [link, comment]

			else:
				raise data.InvalidDatabaseException('Kinetics library data format is unrecognized.')

	def loadTemplate(self, path):
		"""
		Load and process a reaction template file located at `path`. This file
		is part of every reaction family.
		"""

		# Process the template file, removing comments and empty lines
		info = ''
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
			return
		finally:
			frec.close()

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
		# The first tree of each species is always used to identify reactions,
		# while the others are used to select the appropriate kinetics
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
				# Check that all structures are in dictionary
				if item not in self.dictionary:
					raise data.InvalidDatabaseException('Reaction family template contains an unknown structure.')
				species.append(item)

		# Set template reaction
		self.template = Reaction(reactants, products)

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

			# Remaining items are comma-delimited list of parameters enclosed by
			# {}, which we will split into individual parameters
			action.extend(line[len(items[0]):].strip()[1:-1].split(','))

			self.recipe.addAction(action)

		# Generate the reverse template
		if reverse != self.label:
			template = Reaction(self.template.products, self.template.reactants)
			self.reverse = ReactionFamily(reverse, template, self.recipe.getReverse())
			self.reverse.dictionary = self.dictionary
			self.reverse.tree = self.tree
			self.reverse.library = self.library
			self.reverse.forbidden = self.forbidden
			self.reverse.reverse = self

	def reactantMatch(self, reactant, templateReactant):
		"""
		Return :data:`True` if the provided reactant matches the provided
		template reactant and :data:`False` if not.
		"""
		maps12 = []; maps21 = []
		if templateReactant.__class__ == list: templateReactant = templateReactant[0]
		structure = self.dictionary[templateReactant]
		if structure.__class__ == str or structure.__class__ == unicode:
			if structure.lower() == 'union':
				for child in self.tree.children[templateReactant]:
					ismatch, map12, map21 = self.reactantMatch(reactant, child)
					maps12.extend(map12); maps21.extend(map21)
		elif structure.__class__ == chem.Structure:
			return reactant.findSubgraphIsomorphisms(structure)

		return len(maps12) > 0, maps12, maps21

	def makeReaction(self, reactants, reactantStructures, maps):
		"""
		Create a reaction involving a list of `reactants`. The `reactantStructures`
		parameter is a list of structures in the order the reactants are stored
		in the reaction family template, and the `maps` parameter is a list of
		mappings of the top-level tree node of each template reactant to the
		corresponding structure.
		"""

#		# Find atoms associated with each label
#		atoms = {}
#		for map in maps:
#			atom = {}
#			for key, value in map.iteritems():
#				if key.label != '':
#					atoms[key.label] = value

		# Tag atoms with labels
		for struct in reactantStructures: struct.clearLabeledAtoms()
		for map in maps:
			for key, value in map.iteritems():
				if key.label != '':
					value.label = key.label
					
		# Merge reactant structures into single structure
		structure = chem.Structure()
		for struct in reactantStructures:
			structure = structure.merge(struct)

		# Generate the product structure
		self.recipe.applyForward(structure)
		productStructure = structure.copy()
		self.recipe.applyReverse(structure)
		
		# Split product structure into multiple species if necessary
		if len(self.template.products) > 1:
			productStructures = productStructure.split()
		else:
			productStructures = [productStructure]

		# Convert structure(s) to products
		products = []
		for product in productStructures:
			spec = species.makeNewSpecies(product)
			products.append(spec)

#		# Determine top level nodes
#		reactantTemplate = self.template.reactants[:]
#		productTemplate = self.template.products[:]
#		print reactantTemplate, productTemplate
#
#		# Descend all reactant trees as deeply as possible
#		for i in range(len(reactantStructures)):
#			if reactantTemplate[i].__class__ == list:
#				for j in range(len(reactantTemplate[i])):
#					reactantTemplate[i][j] = self.descendTree(reactantStructures[i], atoms, reactantTemplate[i][j])
#			else:
#				reactantTemplate[i] = self.descendTree(reactantStructures[i], atoms, reactantTemplate[i])
#
#
#		# Descend all product trees as deeply as possible
#		for i in range(len(productStructures)):
#			if productTemplate[i].__class__ == list:
#				for j in range(len(productTemplate[i])):
#					productTemplate[i][j] = self.descendTree(productStructures[i], atoms2, productTemplate[i][j])
#			else:
#				productTemplate[i] = self.descendTree(productStructures[i], atoms2, productTemplate[i])
#
#		print reactantTemplate, productTemplate

		# Create reaction and add if unique
		rxn, isNew = makeNewReaction(reactants, products, self)
		if isNew:	return rxn
		else:		return None

	def getReactionList(self, reactants):
		"""
		Generate a list of all of the possible reactions of this family between
		the list of `reactants`.
		"""
		rxnList = []
		# If the number of reactants provided does not match the number of
		# reactants in the template, return False
		if len(reactants) == 1 and self.template.isUnimolecular():

			# Iterate over all resonance isomers of the reactant
			for structure in reactants[0].structure:

				# Make a copy of structure so we don't modify the original
				structureCopy = structure.copy()

				ismatch, map12, map21 = self.reactantMatch(structureCopy, self.template.reactants[0])
				for map in map12:
					rxn = self.makeReaction(reactants, [structureCopy], [map])
					if rxn is not None:
						rxnList.append(rxn)

		# Bimolecular reactants: A + B --> products
		elif len(reactants) == 2 and self.template.isBimolecular():

			# Iterate over all resonance isomers of the reactant
			for structureA in reactants[0].structure:
				for structureB in reactants[1].structure:

					# Make a copy of structure so we don't modify the original
					structureACopy = structureA.copy()
					structureBCopy = structureB.copy()

					# Reactants stored as A + B
					ismatch_A, map12_A, map21_A = self.reactantMatch(structureACopy, self.template.reactants[0])
					ismatch_B, map12_B, map21_B = self.reactantMatch(structureBCopy, self.template.reactants[1])

					# Iterate over each pair of matches (A, B)
					for mapA in map12_A:
						for mapB in map12_B:
							rxn = self.makeReaction(reactants, [structureACopy, structureBCopy], [mapA, mapB])
							if rxn is not None:
								rxnList.append(rxn)

					# Reactants stored as B + A
					ismatch_0, map12_A, map21_A = self.reactantMatch(structureACopy, self.template.reactants[1])
					ismatch_1, map12_B, map21_B = self.reactantMatch(structureBCopy, self.template.reactants[0])

					# Iterate over each pair of matches (A, B)
					for mapA in map12_A:
						for mapB in map12_B:
							rxn = self.makeReaction(reactants, [structureBCopy, structureACopy], [mapB, mapA])
							if rxn is not None:
								rxnList.append(rxn)

		return rxnList

	def getKinetics(self, reaction):
		"""
		Determine the appropriate kinetics for `reaction` which involves the
		labeled atoms in `atoms`.
		"""
		return None
		#data = self.library.getData(template)
		#return data
	
#		# Fail if reaction family is not self
#		if reaction.family is not self: return None
#		print self.template
#
#		# Descend each tree as far as possible
#		nodes = []
#		for top in self.tree.top:
#
#			for templateReactant in self.template.reactants:
#
#				if top == templateReactant or top in templateReactant:
#					for reactant in reaction.reactants:
#						for atoms in atomsList:
#							node = self.descendTree(reactant.structure[0], atoms, top)
#							if node is not None:
#								nodes.append(node)
#
#		print nodes

#						result = ''
#						while node is not None:
#							result = ' -> ' + node + result
#							node = self.tree.parent[node]
#						print result[4:]

################################################################################

class ReactionFamilySet:
	"""
	Represent a set of reaction families.
	"""

	def __init__(self):
		self.families = {}

	def load(self, datapath):
		"""
		Load a set of reaction families from the general database
		specified at `datapath`.
		"""

		logging.debug('\tReaction families:')

		# Load the families from kinetics/families.txt
		familyList = []
		try:
			ffam = open(datapath + 'kinetics/families.txt', 'r')
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
			return
		finally:
			ffam.close()

		# Load the reaction families (if they exist and status is 'on')
		self.families = {}
		for index, status, label in familyList:
			path = datapath + 'kinetics/' + label
			if os.path.isdir(path) and status.lower() == 'on':
				family = ReactionFamily(label)
				family.load(path)
				self.families[family.label] = family
				logging.debug('\t\t' + family.label)
				if family.reverse is not None:
					self.families[family.reverse.label] = family.reverse
					logging.debug('\t\t' + family.reverse.label)

	def getReactions(self, species):
		"""
		Generate a list of reactions that involve a list of one or two `species`
		as a reactant or product.
		"""

		log = str(species[0])
		for spec in species[1:]: log += ' + ' + str(spec)

		rxnList = []
		for key, family in self.families.iteritems():
			rxnList.extend(family.getReactionList(species))

		if len(rxnList) == 1:
			logging.info('Found %s reaction for %s' % (len(rxnList), log))
		else:
			logging.info('Found %s reactions for %s' % (len(rxnList), log))


		return rxnList

kineticsDatabase = None

################################################################################

class Reaction:
	"""
	Represent a generic chemical reaction. A reaction is defined by the set of
	reactants, set of products, and the transition state; thus both the forward
	and reverse reaction are represented by a single object since they share 
	the same transition state. By convention, the forward reaction is taken to
	be that for which the provided kinetics apply; the reverse kinetics are 
	taken from thermodynamic reversibility.
	"""
	
	def __init__(self, reactants=None, products=None, family=None, kinetics=None):
		"""Initialize a reaction object."""
		self.reactants = reactants or []
		self.products = products or []
		self.family = family
		self.reverse = None
		self.kinetics = kinetics
	
	def __str__(self):
		"""
		Return a string representation of the reaction.
		"""
		string = ''
		for reactant in self.reactants:
			string += str(reactant) + ' + '
		string = string[0:-3] + ' <---> '
		for product in self.products:
			string += str(product) + ' + '
		string = string[0:-3]
		return string

	def isUnimolecular(self):
		"""
		Return :data:`True` if the forward reaction has one reactant and
		:data:`False` otherwise.
		"""
		return len(self.reactants) == 1

	def isBimolecular(self):
		"""
		Return :data:`True` if the forward reaction has two reactants and
		:data:`False` otherwise.
		"""
		return len(self.reactants) == 2

	def equivalent(self, other):
		"""
		Return :data:`True` if the two reactions are equivalent (i.e. they have
		the same reactants and products and are of the same reaction family) and
		:data:`False` otherwise.
		"""

		if len(self.reactants) != len(other.reactants) or \
		  len(self.products) != len(other.products):
			return False
		elif self.family is not other.family:
			return False

		reactantsMatch = False
		if len(self.reactants) == 1:
			indices = [[0]]
		elif len(self.reactants) == 2:
			indices = [[0, 1], [1, 0]]
		elif len(self.reactants) == 3:
			indices = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
		for index in indices:
			if reactantsMatch: continue
			match = True
			for i in range(len(self.reactants)):
				if not self.reactants[i].isIsomorphic(other.reactants[index[i]]):
					match = False
			if match:
				reactantsMatch = True

		productsMatch = False
		if len(self.products) == 1:
			indices = [[0]]
		elif len(self.products) == 2:
			indices = [[0, 1], [1, 0]]
		elif len(self.products) == 3:
			indices = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
		for index in indices:
			if productsMatch: continue
			match = True
			for i in range(len(self.reactants)):
				if not self.products[i].isIsomorphic(other.products[index[i]]):
					match = False
			if match:
				productsMatch = True

		return reactantsMatch and productsMatch

################################################################################

# The global list of reactions created at any point during RMG execution
# The list is stored in reverse of the order in which the reactions are created;
# when searching the list, it is more likely to match a recently created
# reaction than an older reaction
reactionList = []

def makeNewReaction(reactants, products, family):
	"""
	Attempt to make a new reaction based on a list of `reactants` and a list of
	`products`. The combination of these and a reaction `family` string uniquely
	identifies a reaction. The reactant and product lists must contain 
	:class:`Species` objects, not :class:`Structure` objects.

	The proposed reaction is checked against the list of
	existing reactions; if the reaction already exists, this function returns
	the existing reaction. If the reaction does not exist, a :class:`Reaction`
	object is created and returned after being appended to the global reaction
	list.
	"""
	
	# Sort reactants and products (to make comparisons easier/faster)
	reactants.sort()
	products.sort()

	# Check that the reaction actually results in a different set of species
	if len(reactants) == len(products):
		match = True
		for i in range(len(reactants)):
			if reactants[i] != products[i]: match = False
		if match: return None, False

	# Check that the reaction is unique
	for rxn in reactionList:
		# Check forward reaction for match
		if rxn.family is family or rxn.family is None or family is None:
			if len(rxn.reactants) == len(reactants) and len(rxn.products) == len(products):
				match = True
				for i in range(len(reactants)):
					if rxn.reactants[i] != reactants[i]: match = False
				for i in range(len(products)):
					if rxn.products[i] != products[i]: match = False
				if match: return rxn, False
		# Check reverse reaction for match
		if rxn.family.reverse is family or rxn.family.reverse is None or family is None:
			if len(rxn.reactants) == len(products) and len(rxn.products) == len(reactants):
				match = True
				for i in range(len(reactants)):
					if rxn.products[i] != reactants[i]: match = False
				for i in range(len(products)):
					if rxn.reactants[i] != products[i]: match = False
				if match: return rxn, False

	# If this point is reached, the proposed reaction is new, so make new
	# Reaction objects for forward and reverse reaction
	forward = Reaction(reactants, products, family)
	reverse = Reaction(products, reactants, family.reverse or family)
	forward.reverse = reverse
	reverse.reverse = forward

	# Attempt to get the kinetics of the forward and reverse reactions
	forwardKinetics = forward.family.getKinetics(forward)
	reverseKinetics = reverse.family.getKinetics(reverse)
	
	# By convention, we only work with the reaction in the direction for which
	# we have assigned kinetics from the kinetics database; the kinetics of the
	# reverse of that reaction come from thermodynamics
	rxn = forward
	if forwardKinetics is None and reverseKinetics is None:
		pass
		#raise Exception('Unable to determine kinetics of reaction ' + str(forward) + '.')
	elif forwardKinetics is None: rxn = reverse
	
	reactionList.insert(0, rxn)
	
	# Note in the log
	logging.debug('Created new ' + str(rxn.family) + ' reaction ' + str(rxn))

	# Return newly created reaction
	return rxn, True

################################################################################

if __name__ == '__main__':
	
	import os.path
	import main
	main.initializeLog(logging.DEBUG)

	datapath = '../data/RMG_database/'

	logging.debug('General database: ' + os.path.abspath(datapath))
	species.thermoDatabase = species.ThermoDatabaseSet()
	species.thermoDatabase.load(datapath)
	kineticsDatabase = ReactionFamilySet()
	kineticsDatabase.load(datapath)

	structure1 = chem.Structure()
	structure1.fromAdjacencyList("""HXD13
1 C 0 {2,D} {7,S} {8,S}
2 C 0 {1,D} {3,S} {9,S}
3 C 0 {2,S} {4,D} {10,S}
4 C 0 {3,D} {5,S} {11,S}
5 C 1 {4,S} {6,S} {12,S}
6 C 0 {5,S} {14,S} {15,S} {16,S}
7 H 0 {1,S}
8 H 0 {1,S}
9 H 0 {2,S}
10 H 0 {3,S}
11 H 0 {4,S}
12 H 0 {5,S}
14 H 0 {6,S}
15 H 0 {6,S}
16 H 0 {6,S}
""")

	species1 = species.makeNewSpecies(structure1, 'C6H9J', True)

	#print len(species1.structure)

	structure2 = chem.Structure()
	structure2.fromSMILES('[H][H]')
	species2 = species.makeNewSpecies(structure2, 'H2', True)

	rxnList = kineticsDatabase.getReactions([species1])
	#rxnList = kineticsDatabase.getReactions([species1, species2])
	for rxn in rxnList:
		print rxn
