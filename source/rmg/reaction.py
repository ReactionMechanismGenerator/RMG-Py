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
import math
import logging
import os
import os.path

import constants
import chem
import data
import structure
import species

################################################################################

class Kinetics:
	"""
	Represent a set of kinetic data. No matter how the kinetic data is modeled,
	it should be accompanied by a temperature range over which it is valid,
	an integer rank denoting the degree of confidence in the data (1 = high,
	5 = low, 0 = none), and a comment describing the source of the data.
	"""

	def __init__(self, Trange=None, rank=0, comment=''):
		self.Trange = Trange
		self.rank = 0
		self.comment = ''
		
	def isTemperatureInRange(self, T):
		"""
		Return :data:`True` if temperature `T` is within the valid temperature
		range specified by self.Trange and :data:`False` if not. If no
		temperature range is specified in self.Trange, the kinetic data is
		assumed to be valid at all temperatures.
		"""
		if self.Trange is not None:
			if T < self.Trange[0] or T > self.Trange[1]:
				return False
		return True

################################################################################

class ArrheniusKinetics(Kinetics):
	"""
	Represent a set of modified Arrhenius kinetics. The kinetic expression has
	the form

	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_\mathrm{a}}{RT} \\right)

	"""

	def __init__(self, A=0.0, Ea=0.0, n=0.0):
		Kinetics.__init__(self)
		self.A = A
		self.Ea = Ea
		self.n = n

	def __str__(self):
		return 'k(T) = %s * T ** %s * math.exp(-%s / constants.R / T)\t%s < T < %s' % (self.A, self.n, self.Ea, self.Trange[0], self.Trange[1])

	def getRateConstant(self, T):
		"""
		Return the rate constant k(T) at temperature `T` by evaluating the
		Arrhenius expression.
		"""

		# Raise exception if T is outside of valid temperature range
		if not self.isTemperatureInRange(T):
			raise Exception('Attempted to evaluate a rate constant at an invalid temperature.')

		return self.A * (T ** self.n) * math.exp(-self.Ea / constants.R / T)

	def getReverse(self, dHrxn, Keq, T):
		"""
		Generate the reverse of the current kinetics for a reaction with
		standard enthalpy of reaction `dHrxn` and equilibrium constant `Keq` at
		298 K, respectively, defined in the same direction that these kinetics
		are.
		"""
		
		kinetics = ArrheniusKinetics(self.A / Keq * math.exp(-dHrxn / constants.R / T), self.Ea - dHrxn, self.n)
		kinetics.Trange = self.Trange
		kinetics.rank = self.rank
		kinetics.comment = self.comment

		return kinetics

	def toXML(self, dom, root, numReactants):
		"""
		Generate the XML for these kinetics using the :data:`xml.dom.minidom`
		package. The `dom` and `root` parameters refer to the DOM and the
		point within the DOM to place this item.
		"""
		kinetics = dom.createElement('arrheniusKinetics')
		root.appendChild(kinetics)
		kinetics.setAttribute('Trange', '%s-%s K' % (self.Trange[0], self.Trange[1]))
		kinetics.setAttribute('rank', str(self.rank))
		kinetics.setAttribute('comment', self.comment)

		preexponential = dom.createElement('preexponential')
		kinetics.appendChild(preexponential)
		preexponentialUnits = None
		if len(self.reactants) == 1:
			preexponentialUnits = 's^-1'
		else:
			preexponentialUnits = 'm^%s/(mol^%s*s)' % ((numReactants-1)*3, numReactants-1)
		data.createXMLQuantity(dom, preexponential, self.A, preexponentialUnits)

		exponent = dom.createElement('exponent')
		kinetics.appendChild(exponent)
		data.createXMLQuantity(dom, exponent, self.n, '')

		activationEnergy = dom.createElement('activationEnergy')
		kinetics.appendChild(activationEnergy)
		data.createXMLQuantity(dom, activationEnergy, self.Ea, 'J/mol')

################################################################################

class ArrheniusEPKinetics(Kinetics):
	"""
	Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
	kinetic expression has the form

	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_\mathrm{a}}{RT} \\right)

	The parameter :math:`\\alpha` is used to correct the activation energy
	:math:`E_\\mathrm{a}` via the Evans-Polanyi formula

	.. math:: E_\\mathrm{a} = E_0 + \\alpha \\Delta H_\\mathrm{rxn}

	"""

	def __init__(self, A=0.0, E0=0.0, n=0.0, alpha=0.0):
		Kinetics.__init__(self)
		self.A = A
		self.E0 = E0
		self.n = n
		self.alpha = alpha

	def __str__(self):
		return 'k(T) = %s * T ** %s * math.exp(-(%s + %s * DHrxn) / constants.R / T)\t%s < T < %s' % (self.A, self.n, self.E0, self.alpha, self.Trange[0], self.Trange[1])

	def getActivationEnergy(self, dHrxn):
		"""
		Return the activation energy using the enthalpy of reaction `dHrxn`.
		"""
		return self.E0 + self.alpha * dHrxn
		
	def getArrhenius(self, dHrxn):
		"""
		Return the Arrhenius form of k(T) at temperature `T` by correcting E0
		to Ea using the enthalpy of reaction `dHrxn`.
		"""

		Ea = self.getActivationEnergy(dHrxn)

		kinetics = ArrheniusKinetics(self.A, Ea, self.n)
		kinetics.Trange = self.Trange
		kinetics.rank = self.rank
		kinetics.comment = self.comment

		return kinetics

	def getRateConstant(self, T, dHrxn):
		"""
		Return the rate constant k(T) at temperature `T` by evaluating the
		Arrhenius expression. The reaction has an enthalpy of reaction `dHrxn`.
		"""

		# Raise exception if T is outside of valid temperature range
		#if not self.isTemperatureInRange(T):
		#	raise Exception('Attempted to evaluate a rate constant at an invalid temperature.')

		Ea = self.getActivationEnergy(dHrxn)

		return self.A * T ** self.n * math.exp(-Ea / constants.R / T)

	def fromDatabase(self, data, comment, numReactants):
		"""
		Process a list of numbers `data` and associated description `comment`
		generated while reading from a kinetics database. The `numReactants`
		parameter is used to interpret the units of the preexponential.
		"""

		if len(data) != 11:
			raise Exception('Invalid list of kinetic data; should be a list of numbers of length 11.')

		Tmin, Tmax, A, n, alpha, E0, dA, dn, dalpha, dE0, rank = data

		originalUnits = 's^-1'; desiredUnits = 's^-1'
		if numReactants == 2:
			originalUnits = 'cm^3/(mol*s)'
		elif numReactants > 2:
			originalUnits = 'cm^%s/(mol^%s*s)' % ((numReactants-1)*3, numReactants-1)
		
		self.Trange = pq.Quantity([Tmin, Tmax], 'K').simplified
		self.Trange = [float(self.Trange[0]), float(self.Trange[1])]

		self.A = float(pq.Quantity(A, originalUnits).simplified)
		self.E0 = float(pq.Quantity(E0, 'kcal/mol').simplified)
		self.n = float(pq.Quantity(n, '').simplified)
		self.alpha = float(pq.Quantity(alpha, '').simplified)

		self.rank = rank
		self.comment = comment

	def toXML(self, dom, root, numReactants):
		"""
		Generate the XML for these kinetics using the :data:`xml.dom.minidom`
		package. The `dom` and `root` parameters refer to the DOM and the
		point within the DOM to place this item.
		"""
		kinetics = dom.createElement('arrheniusEPKinetics')
		root.appendChild(kinetics)
		kinetics.setAttribute('Trange', '%s-%s K' % (self.Trange[0], self.Trange[1]))
		kinetics.setAttribute('rank', str(self.rank))
		kinetics.setAttribute('comment', self.comment)


		preexponential = dom.createElement('preexponential')
		kinetics.appendChild(preexponential)
		preexponentialUnits = None
		if numReactants == 1:
			preexponentialUnits = 's^-1'
		else:
			preexponentialUnits = 'm^%s/(mol^%s*s)' % ((numReactants-1)*3, numReactants-1)
		data.createXMLQuantity(dom, preexponential, self.A, preexponentialUnits)

		exponent = dom.createElement('exponent')
		kinetics.appendChild(exponent)
		data.createXMLQuantity(dom, exponent, self.n, '')

		alpha = dom.createElement('evansPolanyiSlope')
		kinetics.appendChild(alpha)
		data.createXMLQuantity(dom, alpha, self.alpha, '')

		activationEnergy = dom.createElement('evansPolanyiIntercept')
		kinetics.appendChild(activationEnergy)
		data.createXMLQuantity(dom, activationEnergy, self.E0, 'J/mol')

################################################################################

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

	- CHANGE_BOND {center1,order,center2} - change the bond order of the bond between center1 and center2 by order; do not break or form bonds

	- FORM_BOND {center1,order,center2} - form a new bond between center1 and center2 of type order

	- BREAK_BOND {center1,order,center2} - break the bond between center1 and center2, which should be of type order

	- GAIN_RADICAL {center,radical} - Increase the number of free electrons on center by radical

	- LOSE_RADICAL {center,radical} - Decrease the number of free electrons on center by radical

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

	def apply(self, struct, doForward):
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

					# If found, change bond
					if action[0] == 'CHANGE_BOND':
						bond = struct.getBond(atom1, atom2)
						info = int(info)
						if bond is None: raise InvalidActionException('Attempted to change the bond order of a nonexistent bond.')
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
							struct.addBond(bond)
						else:
							bond = struct.getBond(atom1, atom2)
							if bond is None: raise InvalidActionException('Attempted to remove a nonexistent bond.')
							struct.removeBond(bond)
					elif action[0] == 'BREAK_BOND':
						if doForward:
							bond = struct.getBond(atom1, atom2)
							if bond is None: raise InvalidActionException('Attempted to remove a nonexistent bond.')
							struct.removeBond(bond)
						else:
							bond = chem.Bond([atom1, atom2], info)
							struct.addBond(bond)

				elif action[0] == 'LOSE_RADICAL' or action[0] == 'GAIN_RADICAL':

					label, change = action[1:]
					change = int(change)
					# Find associated atoms
					atom = struct.getLabeledAtom(label)
					if atom is None:
						raise Exception('Invalid atom labels found while attempting to execute reaction recipe.')

					# If found, adjust radical
					for i in range(0, change):
						if doForward:
							if action[0] == 'LOSE_RADICAL':
								if not atom.canDecreaseFreeElectron(): raise InvalidActionException('Attempted to decrease the number of free electrons below the minimum.')
								atom.decreaseFreeElectron()
							elif action[0] == 'GAIN_RADICAL':
								if not atom.canIncreaseFreeElectron(): raise InvalidActionException('Attempted to increase the number of free electrons above the maximum.')
								atom.increaseFreeElectron()
						else:
							if action[0] == 'LOSE_RADICAL':
								if not atom.canIncreaseFreeElectron(): raise InvalidActionException('Attempted to increase the number of free electrons above the maximum.')
								atom.increaseFreeElectron()
							elif action[0] == 'GAIN_RADICAL':
								if not atom.canDecreaseFreeElectron(): raise InvalidActionException('Attempted to decrease the number of free electrons below the minimum.')
								atom.decreaseFreeElectron()

				else:
					raise InvalidActionException('Unknown action "' + action[0] + '" encountered.')
		
		except InvalidActionException, e:
			logging.warning('Warning: ' + e.msg)
			return False

		return True

	def applyForward(self, struct):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures.
		"""
		return self.apply(struct, True)

	def applyReverse(self, struct):
		"""
		Apply the reaction recipe to the set of molecules contained in
		`structure`, a single Structure object that contains one or more
		structures. 
		"""
		return self.apply(struct, False)

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

	def getAllNodeCombinations(self, nodeLists):
		"""
		Generate a list of all possible combinations of reactant nodes.
		"""

		items = [[]]
		for nodeList in nodeLists:
			items = [ item + [node] for node in nodeList for item in items ]

		return items

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

		# Fill in missing nodes in library via an averaging scheme of the
		# existing data; note that this disregards all temperature range
		# information
		forwardTemplate, reverseTemplate = self.getTemplateLists()
		#self.generateMissingEntriesFromBelow(forwardTemplate)
		#self.generateMissingEntriesFromAbove(forwardTemplate)

	def generateMissingEntriesFromBelow(self, nodes):
		"""
		Generate a nonexisting entry in the library based on an averaging
		scheme.
		"""

		# Get all possible combinations of child nodes
		children = []
		for node in nodes:
			children.append(self.tree.children[node])
		childNodesList = self.getAllNodeCombinations(children)

		# Only generate new entry if data does not exist
		if self.library.getData(nodes) is None:

			# Check all combinations of child nodes for existing kinetics
			kinetics = []
			for childNodes in childNodesList:
				k = self.library.getData(childNodes)
				if k is not None:
					kinetics.append(k)

			# Average all of the kinetics parameters found above
			if len(kinetics) > 0:
				kin = self.averageKinetics(kinetics)
				self.library.add(nodes, kin)
				#print nodes, kin

		# Recursively descend child nodes
		for childNodes in childNodesList:
			self.generateMissingEntriesFromBelow(childNodes)

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

		kin = ArrheniusEPKinetics(math.exp(lnA), E0, n, alpha)
		kin.Trange = [0.0, 0.0]
		return kin

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
		nodesList = self.getAllNodeCombinations(nodeLists)

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
		parentNodesList = self.getAllNodeCombinations(parentNodeLists)

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
					for i in range(10, len(items)):
						comment += items[i] + ' '
					
					kinetics = ArrheniusEPKinetics()
					kinetics.fromDatabase(kineticData, comment, len(self.template.reactants))
					kinetics.comment = self.label + ' ' + label + ' ' + kinetics.comment
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
		struct = self.dictionary[templateReactant]
		if struct.__class__ == str or struct.__class__ == unicode:
			if struct.lower() == 'union':
				for child in self.tree.children[templateReactant]:
					ismatch, map21, map12 = self.reactantMatch(reactant, child)
					maps12.extend(map12); maps21.extend(map21)
		elif struct.__class__ == structure.Structure:
			return reactant.findSubgraphIsomorphisms(struct)

		return len(maps12) > 0, maps21, maps12

	def makeReaction(self, reactants, structures, maps):
		"""
		Create a reaction involving a list of `reactants`. The `reactantStructures`
		parameter is a list of structures in the order the reactants are stored
		in the reaction family template, and the `maps` parameter is a list of
		mappings of the top-level tree node of each template reactant to the
		corresponding structure.
		"""

		# Clear any previous atom labeling from all structures
		for struct in structures: struct.clearLabeledAtoms()

		# Tag atoms with labels
		counter = 0
		for map in maps:
			for key, value in map.iteritems():
				# For reactions involving two identical centers, both labeled
				# '*', label one as '*1' and the other as '*2'
				# An example reaction family is colligation (reverse of
				# unimolecular homolysis)
				if key.label == '*':
					counter += 1; value.label = '*' + str(counter)
				# Otherwise pass label as given
				elif key.label != '':
					value.label = key.label

		# Copy structures so we don't modify the originals
		# Do this after tagging the originals so both reactants and products
		# have tags
		reactantStructures = []
		for struct in structures:
			reactantStructures.append(struct.copy())

		# Merge reactant structures into single structure
		struct = structure.Structure()
		for s in reactantStructures:
			struct = struct.merge(s)

		# Generate the product structure
		if not self.recipe.applyForward(struct):
			return None

		productStructure = struct.copy()
		self.recipe.applyReverse(struct)
		
		# Restore original atom labels of the reactants if they were changed
		# before
		if counter > 0:
			#for struct in structures:
			for s in reactantStructures:
				for atom in s.atoms():
					if atom.label != '':
						atom.label = '*'

		# Remove numbers from labeled atoms of products if they are the same
		# top-level node
		if len(self.template.products) > 1:
			if self.template.products[0] == self.template.products[1]:
				for atom in productStructure.atoms():
					if atom.label != '': atom.label = '*'

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
				atomLabels['*1'].label = '*3'
				atomLabels['*3'].label = '*1'
			elif label == 'intra h migration' or label == 'alkyl hydroperoxyl intra OH migration':
				atomLabels['*1'].label = '*2'
				atomLabels['*2'].label = '*1'


		# Split product structure into multiple species if necessary
		if len(self.template.products) > 1:
			productStructures = productStructure.split()
		else:
			productStructures = [productStructure]

		# Recalculate atom types of product structures, since they may have
		# changed as a result of the reaction
		for struct in productStructures:
			struct.simplifyAtomTypes()
			struct.updateAtomTypes()

		# Check that reactant and product structures are allowed in this family
		# If not, then stop
		if self.forbidden is not None:
			for label, struct2 in self.forbidden.iteritems():
				for struct in reactantStructures:
					match, map21, map12 = struct.isSubgraphIsomorphic(struct2)
					if match: return None
				for struct in productStructures:
					match, map21, map12 = struct.isSubgraphIsomorphic(struct2)
					if match: return None

		# Convert structure(s) to products
		products = []
		for product in productStructures:
			spec = species.makeNewSpecies(product)
			# Don't make a new reaction if no species was returned from
			# makeNewSpecies() (e.g. due to forbidden structure)
			if spec is None: return None
			products.append(spec)

		# Create reaction and add if unique
		rxn, isNew = makeNewReaction(reactants, products, reactantStructures, productStructures, self)
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

				ismatch, map21, map12 = self.reactantMatch(structure, self.template.reactants[0])
				for map in map12:
					rxn = self.makeReaction(reactants, [structure], [map])
					if rxn is not None:
						rxnList.append(rxn)

		# Bimolecular reactants: A + B --> products
		elif len(reactants) == 2 and self.template.isBimolecular():

			# Make copies of the structure lists of the two reactants
			# This is a workaround for an issue in which the two reactant 
			# structure lists were getting swapped around, resulting in 
			# unbalanced reactions
			# The copy is needed for cases where A and B are the same
			structuresA = []; structuresB = []
			for structureA in reactants[0].structure:
				structuresA.append(structureA.copy())
			for structureB in reactants[1].structure:
				structuresB.append(structureB.copy())

			# Iterate over all resonance isomers of the reactant
			for structureA in structuresA:
				for structureB in structuresB:

					# Reactants stored as A + B
					ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[0])
					ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[1])

					# Iterate over each pair of matches (A, B)
					for mapA in map12_A:
						for mapB in map12_B:
							rxn = self.makeReaction(reactants, [structureA, structureB], [mapA, mapB])
							if rxn is not None:
								rxnList.append(rxn)

					# Only check for swapped reactants if they are different
					if reactants[0].id != reactants[1].id:

						# Reactants stored as B + A
						ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[1])
						ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[0])

						# Iterate over each pair of matches (A, B)
						for mapA in map12_A:
							for mapB in map12_B:
								rxn = self.makeReaction(reactants, [structureB, structureA], [mapB, mapA])
								if rxn is not None:
									rxnList.append(rxn)

		return rxnList

	def getKinetics(self, reaction, structures):
		"""
		Determine the appropriate kinetics for `reaction` which involves the
		labeled atoms in `atoms`.
		"""

		# Get forward reaction template and remove any duplicates
		forwardTemplate, reverseTemplate = self.getTemplateLists()
		forwardTemplate = list(set(forwardTemplate))
		
		# Descend reactant trees as far as possible
		template = []
		for forward in forwardTemplate:
			# Get labeled atoms of forward
			node = forward;	group = self.dictionary[node]
			while group.__class__ == str or group.__class__ == unicode:
				node = self.tree.children[node][0]
				group = self.dictionary[node]

			atomList = group.getLabeledAtoms()
			for structure in structures:
				# Match labeled atoms
				match = True
				for label in atomList:
					if not structure.containsLabeledAtom(label):
						match = False
				# Match structures
				atoms = {}
				node = self.descendTree(structure, atoms, forward)
				if match and node is not None:
					template.append(node)

		# Check that we were able to match the template
		if len(template) == 0:
			logging.warning('Warning: Unable to find matching template for reaction %s in reaction family %s; using the most general reaction template.' % (str(reaction), str(self)))
			# If unable to match template, use the most general template
			forwardTemplate, reverseTemplate = self.getTemplateLists()

			print str(self), template, forwardTemplate, reverseTemplate
			for reactant in reaction.reactants:
				print reactant.toAdjacencyList() + '\n'
			for product in reaction.products:
				print product.toAdjacencyList() + '\n'
			
			template = forwardTemplate

		k = self.library.getData(template)
		print template, k
		if k is not None: return [k]
		else: return None
		
		nodeLists = []
		for temp in template:
			nodeList = []
			while temp is not None:
				nodeList.append(temp)
				temp = self.tree.parent[temp]
			nodeLists.append(nodeList)

		# Generate all possible combinations of nodes
		items = self.getAllNodeCombinations(nodeLists)
		
		# Generate list of kinetics at every node
		kinetics = []
		for item in items:
			data = self.library.getData(item)
			if data is not None:
				kinetics.append(data)

		if len(kinetics) == 0: return None
		
		return kinetics

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
	taken from thermodynamic reversibility. A multiplier attribute is used when
	reactions involving equivalent reactive sites are used.
	"""
	
	def __init__(self, reactants=None, products=None, family=None, kinetics=None):
		"""Initialize a reaction object."""
		self.reactants = reactants or []
		self.products = products or []
		self.family = family
		self.reverse = None
		self.kinetics = kinetics or []

		# A cache for the best kinetics for this reaction
		self.bestKinetics = None

		# A dictionary of the labeled atoms for the reactants
		self.atomLabels = {}
		
		self.multiplier = 1.0

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

	def getEnthalpyOfReaction(self, T):
		"""
		Return the enthalpy of reaction evaluated at temperature `T`.
		"""
		dHrxn = -self.reactants[0].getEnthalpy(T)
		for reactant in self.reactants[1:]:
			dHrxn -= reactant.getEnthalpy(T)
		for product in self.products:
			dHrxn += product.getEnthalpy(T)
		return dHrxn

	def getEntropyOfReaction(self, T):
		"""
		Return the entropy of reaction evaluated at temperature `T`.
		"""
		dSrxn = -self.reactants[0].getEntropy(T)
		for reactant in self.reactants[1:]:
			dSrxn -= reactant.getEntropy(T)
		for product in self.products:
			dSrxn += product.getEntropy(T)
		return dSrxn

	def getFreeEnergyOfReaction(self, T):
		"""
		Return the Gibbs free energy of reaction evaluated at temperature `T`.
		"""
		dGrxn = -self.reactants[0].getFreeEnergy(T)
		for reactant in self.reactants[1:]:
			dGrxn -= reactant.getFreeEnergy(T)
		for product in self.products:
			dGrxn += product.getFreeEnergy(T)
		return dGrxn

	def getEquilibriumConstant(self, T, conc):
		"""
		Return the equilibrium constant K(T) evaluated at temperature `T` in a
		system with total concentration `conc`.
		"""
		dGrxn = self.getFreeEnergyOfReaction(T)
		K = math.exp(-dGrxn / constants.R / T)
		# Convert from Ka to Kc
		K *= conc ** (len(self.products) - len(self.reactants))
		return K

	def getBestKinetics(self, T):
		"""
		Return the best set of kinetic parameters for the forward reaction
		evaluated at the temperature `T`. This function follows the convention
		that the forward reaction is the one for which we are using the kinetic
		expression, and that the reverse rate constant is evaluated using
		thermochemical equilibrium.
		"""

		# Check cache first
		if self.bestKinetics is not None:
			if self.bestKinetics.isTemperatureInRange(T):
				return self.bestKinetics

		kinetics = self.kinetics[:]
		
		# Prune out all kinetic data not valid at desired temperature
		i = 0
		while i < len(kinetics):
			k = kinetics[i]
			if not k.isTemperatureInRange(T): kinetics.remove(k)
			else: i += 1
		
		# If no kinetic parameters are left to choose from, print a warning
		# The reaction rate for the reactions is set to zero
		# This may not be the best course of action
		if len(kinetics) == 0:
			#logging.warning('Warning: No kinetics available for reaction ' + str(self) + ' at ' + str(T) + ' K.')
			kinetics = ArrheniusKinetics(0.0, 0.0, 0.0)
			kinetics.Trange = [0.0, 100000.0]
			return kinetics

		# Choose kinetics based on rank (i.e. lowest non-zero rank)
		bestRank = kinetics[0].rank; bestKinetics = kinetics[0]
		for k in kinetics[1:]:
			if k.rank < bestRank and k.rank != 0:
				bestRank = k.rank
				bestKinetics = k

		# Convert to ArrheniusKinetics and return
		# Use T = 298 K to calculate enthalpy and free energy of reaction
		T = 298.0
		dHrxn = self.getEnthalpyOfReaction(T)
		self.bestKinetics = bestKinetics.getArrhenius(dHrxn)
		
		return self.bestKinetics

	def getRateConstant(self, T):
		"""
		Return the value of the rate constant k(T) at the temperature `T`.
		"""
		kinetics = self.getBestKinetics(T)
		if kinetics is None:
			raise Exception('Unable to determine the rate constant of reaction ' + str(self) + '.')
		return kinetics.getRateConstant(T)

	def getStoichiometricCoefficient(self, spec):
		"""
		Return the stoichiometric coefficient of species `spec` in the reaction.
		The stoichiometric coefficient is increased by one for each time `spec`
		appears as a product and decreased by one for each time `spec` appears
		as a reactant.
		"""
		stoich = 0
		for reactant in self.reactants:
			if reactant is spec: stoich -= 1
		for product in self.products:
			if product is spec: stoich += 1
		return stoich

	def getRate(self, T, P, conc):
		"""
		Return the net rate of reaction at temperature `T` and pressure `P`. The
		parameter `conc` is a map with species as keys and concentrations as
		values. A reactant not found in the `conc` map is treated as having zero
		concentration.
		"""

		# Calculate total concentration
		totalConc = None
		for spec in conc:
			if totalConc is None: totalConc = conc[spec]
			else: totalConc += conc[spec]

		# Evaluate rate constant
		rateConstant = self.getRateConstant(T)

		# Evaluate equilibrium constant
		equilibriumConstant = self.getEquilibriumConstant(T, totalConc)

		# Evaluate forward concentration product
		forward = 1.0
		for reactant in self.reactants:
			if reactant in conc:
				forward = forward * conc[reactant]
			else:
				forward = forward * 0.0

		# Evaluate reverse concentration product
		reverse = 1.0
		for product in self.products:
			if product in conc:
				reverse = reverse * conc[product]
			else:
				reverse = reverse * 0.0

		# Return rate
		return rateConstant * (forward - reverse / equilibriumConstant)

################################################################################

class UndeterminableKineticsException(Exception):
	"""
	An exception raised when attempts to select appropriate kinetic parameters
	for a chemical reaction are unsuccessful.
	"""

	def __init__(self, reaction):
		self.reaction = reaction
		
	def __str__(self):
		string = str(self.reaction) + '\n'
		string += str(self.reaction.family) + '\n'
		for reactant in self.reaction.reactants:
			string += reactant.toAdjacencyList() + '\n'
		for product in self.reaction.products:
			string += product.toAdjacencyList() + '\n'

		return string


################################################################################

# The global list of reactions created at any point during RMG execution
# The list is stored in reverse of the order in which the reactions are created;
# when searching the list, it is more likely to match a recently created
# reaction than an older reaction
reactionList = []

def makeNewReaction(reactants, products, reactantStructures, productStructures, family):
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

	# Get atom labels of reactants
	reactantLabels = {}; productLabels = {}
	for structure in reactantStructures:
		for atom in structure.atoms():
			if atom.label == '*': 
				if atom.label in reactantLabels: 
					reactantLabels[atom.label].append(atom)
				else:
					reactantLabels[atom.label] = [atom]
			elif atom.label != '': reactantLabels[atom.label] = atom
	# Get atom labels of products
	for structure in productStructures:
		for atom in structure.atoms():
			if atom.label == '*':
				if atom.label in productLabels:
					productLabels[atom.label].append(atom)
				else:
					productLabels[atom.label] = [atom]
			elif atom.label != '': productLabels[atom.label] = atom

	# Check that the reaction actually results in a different set of species
	if len(reactants) == len(products):
		match = True
		for i in range(len(reactants)):
			if reactants[i] != products[i]: match = False
		if match: return None, False

	# Check that the reaction is unique
	matchReaction = None
	for rxn in reactionList:
		# Check forward reaction for match
		if rxn.family is family or rxn.family is None or family is None:
			if len(rxn.reactants) == len(reactants) and len(rxn.products) == len(products):
				match = True
				for i in range(len(reactants)):
					if rxn.reactants[i] != reactants[i]: match = False
				for i in range(len(products)):
					if rxn.products[i] != products[i]: match = False
				if match: matchReaction = rxn
		# Check reverse reaction for match
		if rxn.reverse.family is family or rxn.reverse.family is None or family is None:
			if len(rxn.reactants) == len(products) and len(rxn.products) == len(reactants):
				match = True
				for i in range(len(reactants)):
					if rxn.products[i] != reactants[i]: match = False
				for i in range(len(products)):
					if rxn.reactants[i] != products[i]: match = False
				if match: matchReaction = rxn

	# If a match was found, take an
	if matchReaction is not None:
		#matchReaction.multiplier += 1.0
		return matchReaction, False

	# If this point is reached, the proposed reaction is new, so make new
	# Reaction objects for forward and reverse reaction
	forward = Reaction(reactants, products, family)
	reverseFamily = None
	if family is not None: reverseFamily = family.reverse or family
	reverse = Reaction(products, reactants, reverseFamily)
	forward.reverse = reverse
	reverse.reverse = forward

	# Dictionaries containing the labeled atoms for the reactants and products
	forward.atomLabels = reactantLabels
	reverse.atomLabels = productLabels

	if forward.family is None or reverse.family is None:
		reactionList.insert(0, forward)
		return forward, True

	# Attempt to get the kinetics of the forward and reverse reactions
	forwardKinetics = forward.family.getKinetics(forward, reactantStructures)
	reverseKinetics = reverse.family.getKinetics(reverse, productStructures)
	
	# By convention, we only work with the reaction in the direction for which
	# we have assigned kinetics from the kinetics database; the kinetics of the
	# reverse of that reaction come from thermodynamics
	# If we have assigned kinetics in both directions, then (for now) choose the
	# kinetics for the forward reaction
	rxn = forward
	if forwardKinetics is not None and reverseKinetics is not None:
		rxn = forward
		reverseKinetics = []
	elif forwardKinetics is not None:
		rxn = forward
	elif reverseKinetics is not None:
		rxn = reverse
	else:
		raise UndeterminableKineticsException(forward)
		return None, False
	
	forward.kinetics = forwardKinetics
	reverse.kinetics = reverseKinetics

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

	structure1 = structure.Structure()
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

	structure2 = structure.Structure()
	structure2.fromSMILES('[H][H]')
	species2 = species.makeNewSpecies(structure2, 'H2', True)

	rxnList = kineticsDatabase.getReactions([species1])
	#rxnList = kineticsDatabase.getReactions([species1, species2])
	for rxn in rxnList:
		print rxn
