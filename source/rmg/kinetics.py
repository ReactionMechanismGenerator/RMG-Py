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

class ReactionFamily(data.Database):
	"""
	Represent a reaction family: a set of reactions with similar chemistry, and
	therefore similar reaction rates. Attributes include the family name `name`,
	the reaction template `template`, a list `actions` of actions to take when
	reacting, and a dictionary-library-tree `database` of rate rules.
	"""

	def __init__(self, label='', template='', actions='', database=None):
		data.Database.__init__(self)
		self.label = label
		self.template = template
		self.actions = actions
		self.forbidden = None

	def load(self, path):
		"""
		Load a reaction family located in the directory `path`.
		"""

		# Generate paths to files in the database
		dictstr = path + '/dictionary.txt'
		treestr = path + '/tree.txt'
		libstr = path + '/rateLibrary.txt'
		adjstr = path + '/reactionAdjList.txt'
		forbstr = path + '/forbiddenGroups.txt'

		# Load the dictionary, tree, and library using the generic methods
		data.Database.load(self, dictstr, treestr, libstr)

		# Process the data in the library
		self.__processLibraryData()

		# Load the adjlist
		self.__loadAdjList(adjstr)

		# Load the forbidden groups if necessary
		if os.path.exists(forbstr):
			self.forbidden = data.Dictionary()
			self.forbidden.load(forbstr)
			self.forbidden.toStructure()

	def __processLibraryData(self):
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

	def __loadAdjList(self, path):
		"""
		Load and process an adjList file located at `path`. This file is part
		of every reaction family.
		"""

		# Process the adjList file, removing comments and empty lines
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

		# Process adjlist
		lines = info.splitlines()
		# First line is template
		reactants = []; products = []
		atArrow = False
		for item in lines[0].split():
			item = item.strip()
			if item == '+': pass
			elif item.find('->') > -1: atArrow = True
			else:
				if atArrow: products.append(item)
				else: reactants.append(item)
		# Check that all reactants are in dictionary; if so, convert to structures
		for reactant in reactants:
			if reactant not in self.dictionary:
				raise data.InvalidDatabaseException('Reaction family template contains an unknown reactant.')
		# Set template reaction
		self.template = reaction.Reaction(reactants, products)
		# Skip forward/reverse and thermo_consistence information
		index = 1
		while not lines[index].lower().find('actions') and index < len(lines):
			index += 1
		# Read actions
		self.actions = ''
		for line in lines[index:]:
			self.actions += line + '\n'

	def reactantMatch(self, reactant, templateReactant):
		"""
		Return :data:`True` if the provided reactant matches the provided
		template reactant and :data:`False` if not.
		"""
		maps12 = []; maps21 = []
		structure = self.dictionary[templateReactant]
		if structure.__class__ == str or structure.__class__ == unicode:
			if structure.lower() == 'union':
				for child in self.tree.children[templateReactant]:
					ismatch, map12, map21 = self.reactantMatch(reactant, child)
					maps12.extend(map12); maps21.extend(map21)
		elif structure.__class__ == chem.Structure:
			return reactant.findSubgraphIsomorphisms(structure)
		
		return len(maps12) > 0, maps12, maps21
	
	def getReactionList(self, reactants):
		"""
		Generate a list of all of the possible reactions of this family between
		the list of `reactants`.
		"""
		rxnList = []
		# If the number of reactants provided does not match the number of
		# reactants in the template, return False
		if len(reactants) == 1 and self.template.isUnimolecular():
			ismatch, map12, map21 = self.reactantMatch(reactants[0], self.template.reactants[0])
			if len(map12) > 0:
				print self.label, len(map12)
		elif len(reactants) == 2 and self.template.isBimolecular():

			# Reactants stored as A + B
			ismatch_0, map12_0, map21_0 = self.reactantMatch(reactants[0], self.template.reactants[0])
			ismatch_1, map12_1, map21_1 = self.reactantMatch(reactants[1], self.template.reactants[1])
			if len(map12_0) > 0 and len(map12_1) > 0:
				print self.label, len(map12_0), len(map12_1)

			# Reactants stored as B + A
			ismatch_0, map12_0, map21_0 = self.reactantMatch(reactants[0], self.template.reactants[1])
			ismatch_1, map12_1, map21_1 = self.reactantMatch(reactants[1], self.template.reactants[0])
			if len(map12_0) > 0 and len(map12_1) > 0:
				print self.label, len(map12_0), len(map12_1)

		return []

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
				self.families[label] = family
				logging.debug('\t\t' + label)

	def getUnimolecularReactions(self, species):
		"""
		Generate a list of reactions that involve a single `species` as a reactant or product.
		"""
		for key, family in self.families.iteritems():
			rxnList = family.getReactionList([species])
			#print family.label, len(rxnList)

	def getBimolecularReactions(self, species1, species2):
		"""
		Generate a list of reactions that involve a pair of species `species1`
		and `species2` as a set of reactants or products.
		"""
		for key, family in self.families.iteritems():
			rxnList = family.getReactionList([species1, species2])
			#print family.label, len(rxnList)

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

	structure = chem.Structure()
	#structure.fromInChI('InChI=1/C6H10/c1-3-5-6-4-2/h3,5-6H,1,4H2,2H3')
	structure.fromAdjacencyList("""HXD13
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
	species1 = species.Species('HXD13', structure, True)
	structure = chem.Structure()
	structure.fromSMILES('[H][H]')
	species2 = species.Species('H2', structure, True)

	database.getUnimolecularReactions(species1)
	database.getBimolecularReactions(species1, species2)

	