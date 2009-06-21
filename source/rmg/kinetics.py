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

import logging
import quantities as pq
import os.path

import chem
import data
import species
import reaction

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
		#self.Ea.units = 'J/mol'
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
	
	def apply(self, structure, atoms, doForward):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. The `atoms` parameter is a dictionary of the labeled atom
		centers in the structure. The `doForward` parameter is used to indicate
		whether the forward or reverse recipe should be applied.
		"""

		for action in self.actions:
			if action[0] == 'CHANGE_BOND' or action[0] == 'FORM_BOND' or action[0] == 'BREAK_BOND':

				label1, info, label2 = action[1:]
				if label1 not in atoms or label2 not in atoms:
					raise Exception('Invalid atom labels found while attempting to execute reaction recipe.')

				# Find associated atoms
				atom1 = atoms[label1]
				atom2 = atoms[label2]
				if atom1 is None or atom2 is None:
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
				if label not in atoms:
					raise Exception('Invalid atom labels found while attempting to execute reaction recipe.')

				# Find associated atoms
				atom = atoms[label]
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


	def applyForward(self, structure, atoms):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. The `atoms` parameter is a dictionary of the labeled atom
		centers in the structure. The recipe is applied in the forward
		direction.
		"""
		return self.apply(structure, atoms, True)

	def applyReverse(self, structure, atoms):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. The `atoms` parameter is a dictionary of the labeled atom
		centers in the structure. The recipe is applied in the reverse
		direction.
		"""
		return self.apply(structure, atoms, False)

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
		data.Database.load(self, dictstr, treestr, libstr)

		# Process the data in the library
		self.processLibraryData()

		# Load the forbidden groups if necessary
		if os.path.exists(forbstr):
			self.forbidden = data.Dictionary()
			self.forbidden.load(forbstr)
			self.forbidden.toStructure()
		
		# Load the adjlist (must be last, as the reverse family is also generated)
		self.loadTemplate(tempstr)

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
		self.template = reaction.Reaction(reactants, products)
		
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
			template = reaction.Reaction(self.template.products, self.template.reactants)
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

	def makeReaction(self, reactants, structures, maps):
		"""
		Create a reaction.
		"""
		
		# Find atoms associated with each label
		atoms = {}
		for map in maps:
			for key, value in map.iteritems():
				if key.label != '':
					atoms[key.label] = value
		
		# Merge reactant structures into single structure
		structure = chem.Structure()
		for struct in structures:
			structure = structure.merge(struct)
			
		# Generate the product structure
		self.recipe.applyForward(structure, atoms)
		productStructure = structure.copy()
		self.recipe.applyReverse(structure, atoms)

		# Split product structure into multiple species if necessary
		if len(self.template.products) > 1:
			productList = productStructure.split()
		else:
			productList = [productStructure]

		# Convert structure(s) to products
		products = []
		for product in productList:
			spec = species.makeNewSpecies(product)
			products.append(spec)

		# Create reaction and add if unique
		rxn, isNew = reaction.makeNewReaction(reactants, products, self)
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
							rxn = self.makeReaction(reactants, [structureACopy, structureBCopy], [mapA, mapB])
							if rxn is not None:
								rxnList.append(rxn)

		return rxnList

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
		Generate a list of reactions that involve a single `species` as a reactant or product.
		"""
		rxnList = []
		for key, family in self.families.iteritems():
			rxnList.extend(family.getReactionList(species))
		return rxnList

database = None

################################################################################

if __name__ == '__main__':

	import os.path
	import main
	main.initializeLog(logging.DEBUG)

	datapath = '../data/RMG_database/'

	logging.debug('General database: ' + os.path.abspath(datapath))
	database = ReactionFamilySet()
	database.load(datapath)

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

	rxnList = database.getReactions([species1])
	#rxnList = database.getReactions([species1, species2])
	for rxn in rxnList:
		print rxn
	
	