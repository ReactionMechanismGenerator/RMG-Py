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

import pydot
import re

import rmg.main as main
import rmg.data as data
import rmg.species as species
import rmg.reaction as reaction

################################################################################

def loadThermoDatabases(databasePath):
	"""
	Create and load the thermodynamics databases.
	"""
	# Create and load thermo databases
	species.thermoDatabase = species.ThermoDatabaseSet()
	species.thermoDatabase.load(databasePath + '/')

	# Create and load forbidden structures
	species.forbiddenStructures = data.Dictionary()
	species.forbiddenStructures.load(databasePath + '/forbiddenStructure.txt')
	species.forbiddenStructures.toStructure()

def loadKineticsDatabases(databasePath):
	"""
	Create and load the kinetics databases (reaction families).
	"""
	reaction.kineticsDatabase = reaction.ReactionFamilySet()
	reaction.kineticsDatabase.load(databasePath + '/')

################################################################################

def drawKineticsTrees():
	"""
	For each reaction family, output a DOT file containing a combined form of
	the various trees that make up the reaction family.
	"""

	# Iterate through reaction families: key is a string containing the family
	# name, family is the corresponding ReactionFamily object
	for key, family in reaction.kineticsDatabase.families.iteritems():

		print '\tCreating DOT object...'

		graph = family.drawFullGraphOfTree()

		graph.set('fontsize','10')
		format='svg'
		prog='dot'

		print '\tCreating DOT file...'
		f=open(key+'.dot','w')
		f.write(graph.to_string())
		f.close()

		print '\tCreating SVG file...'
		filename=key+'.'+format
		if format=='svg':  # annoyingly, dot creates svg's without units on the font size attribute.
			st=graph.create_svg(prog=prog)
			st=re.sub(r"(font\-size\:[0-9]+\.*[0-9]*)([^p])",r"\1pt\2",st)
			f=open(filename,'w')
			f.write(st)
			f.close()
		else:
			graph.write(filename,format=format,prog=prog)

		print 'Created DOT for reaction family %s' % (key)

################################################################################

if __name__ == '__main__':

	# Show debug messages (as databases are loading)
	main.initializeLog(10)

	# Load databases
	databasePath = '../data/RMG_database'
	#loadThermoDatabases(databasePath)
	loadKineticsDatabases(databasePath)

#	# Prune kinetics dictionaries and trees
#	for key, family in reaction.kineticsDatabase.families.iteritems():
#		forwardTemplate, reverseTemplate = family.getTemplateLists()
#		family.prune(forwardTemplate)
#
#		# Print dictionary
#		print '/////////////////////////////////////////////////////////////////////////////////'
#		print '//'
#		print '// Structure dictionary for %s' % family.label
#		print '//'
#		print '/////////////////////////////////////////////////////////////////////////////////'
#		print ''
#		for label, struct in family.dictionary.iteritems():
#			if struct.__class__ == str:
#				union = label + '\nUnion {' + ','.join(family.tree.children[label]) + '}\n'
#				print union
#			else:
#				print struct.toAdjacencyList(label)
#
#		# Print tree
#		print '/////////////////////////////////////////////////////////////////////////////////'
#		print '//'
#		print '// Structure tree for %s' % family.label
#		print '//'
#		print '/////////////////////////////////////////////////////////////////////////////////'
#		print ''
#		print family.tree.write(family.tree.top)

	# Draw kinetics trees
	#drawKineticsTrees()

	# DECOUPLE LIBRARY DATA
	family = reaction.kineticsDatabase.families['H abstraction']

	# Get available data
	nodesList = []; dataList = []
	for key, value in family.library.iteritems():
		nodes = key.split(';')
		if value not in dataList:
			nodesList.append(nodes)
			dataList.append(value)

	# Get set of all nodes; remove those with no data
	nodeSet = family.tree.parent.keys()
#	nodesToRemove = []
#	for node in nodeSet:
#		hasData = False
#		children = family.tree.descendants(node)
#		for nodes in nodesList:
#			if node in nodes: hasData = True
#			for child in children:
#				if child in nodes: hasData = True
#		if not hasData:
#			nodesToRemove.append(node)
#	for node in nodesToRemove:
#		nodeSet.remove(node)
	
	# Determine size of matrix and vector
	ncols = len(nodeSet) + 1
	nrows = 0
	for nodes in nodesList:
		nodeLists = []
		for node in nodes:
			nodeList = []
			temp = node
			while temp is not None:
				nodeList.append(temp)
				temp = family.tree.parent[temp]
			nodeLists.append(nodeList)
		nodeLists = data.getAllCombinations(nodeLists)
		nrows += len(nodeLists)

	import numpy
	import numpy.linalg
	import math

	# Initialize matrix and vector
	A = numpy.zeros((nrows, ncols), float)
	b = numpy.zeros((nrows, 4), float)

	# Fill in matrix and vector
	rowcount = 0
	for nodes, kinetics in zip(nodesList, dataList):
		nodeLists = []
		for node in nodes:
			nodeList = []
			temp = node
			while temp is not None:
				nodeList.append(temp)
				temp = family.tree.parent[temp]
			nodeLists.append(nodeList)
		nodeLists = data.getAllCombinations(nodeLists)
		for nodeList in nodeLists:
			for i, node in enumerate(nodeList):
				temp = node
				while temp is not None:
					A[rowcount,nodeSet.index(temp)] += 1
					temp = family.tree.parent[temp]

			A[rowcount,len(nodeSet)] += 1
			b[rowcount,0] = math.log(kinetics.A)
			b[rowcount,1] = kinetics.n
			b[rowcount,2] = kinetics.alpha
			b[rowcount,3] = kinetics.E0
			rowcount += 1

	x, residues, rank, s = numpy.linalg.lstsq(A,b)

	for i in range(len(nodeSet)):
		print nodeSet[i], x[i,:]
	print x[len(nodeSet),:]
#	print x
#	print residues
#	print rank, nrows, ncols, 4
#	print s
