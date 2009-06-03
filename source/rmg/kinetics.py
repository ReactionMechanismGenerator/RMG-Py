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

import data

################################################################################

class ArrheniusEPKinetics:
	"""
	Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
	kinetic expression has the form
	
	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_\mathrm{a}}{RT} \\right)
		
	The parameter :math:`\\alpha` is used to correct the activation energy 
	:math:`E_\\mathrm{a}` via the Evans-Polanyi formula
	
	.. math:: E_\\mathrm{a} = E_\\mathrm{a}^0 + (\\alpha - 1) \\Delta H_\\mathrm{rxn} 
	
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

	def __init__(self, name='', template='', actions='', database=None):
		data.Database.__init__(self)
		self.name = name
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
		data.Database.load(self, dictstr, treestr, libstr, True)
		
		# Process the data in the library
		self.__processLibraryData()
		
		# Load the adjlist
		self.__loadAdjList(adjstr)
		
		# Load the forbidden groups if necessary
		if os.path.exists(forbstr):
			self.forbidden = data.Dictionary()
			self.forbidden.load(forbstr)
			self.forbidden.toStructure(False)
		
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
		self.template = lines[0]
		# Skip forward/reverse and thermo_consistence information
		index = 1
		while not lines[index].lower().find('actions') and index < len(lines):
			index += 1
		# Read actions
		self.actions = ''
		for line in lines[index:]:
			self.actions += line + '\n'

################################################################################

def loadReactionFamilies(datapath):
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
	families = {}
	for index, status, label in familyList: 
		path = datapath + 'kinetics/' + label
		if os.path.isdir(path) and status.lower() == 'on':
			family = ReactionFamily()
			family.load(path)
			families[label] = family
			logging.debug('\t\t' + label)
	
	return families
	
################################################################################

if __name__ == '__main__':
	
	#reactant1 = chem.Species('CH4')
	
	#product1 = chem.Species('CH3')
	#product2 = chem.Species('H')
	
	#reaction = Reaction([reactant1], [product1, product2], None)
	#print reaction
	
	datapath = '../data/RMG_database/'
	families = loadReactionFamilies(datapath)
	