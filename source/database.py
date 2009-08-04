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

import logging

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
	import math
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
	
	# Get set of all nodes
	node_set = family.tree.parent.keys()
	non_top_nodes = [ node for node in node_set if family.tree.parent[node] ]
	top_nodes = [ node for node in node_set if not family.tree.parent[node] ]

	A_list = []
	b_list = []
	# Get available data
	nodesList = []; dataList = []
	for key, kinetics in family.library.iteritems():
		nodes = key.split(';')
		# example:
		#  nodes = ['A11', 'B11']
		#  kinetics = <rmg.reaction.ArrheniusEPKinetics instance>
		
		b_row = [ math.log(kinetics.A),
				  kinetics.n,
				  kinetics.alpha,
				  kinetics.E0 ]
			
		all_ancestors=[]
		for node in nodes:
			# start with the most specific - the node itself
			# then add the ancestors
			ancestors = [node]
			ancestors.extend( family.tree.ancestors(node) )
			# append to the list of lists
			all_ancestors.append(ancestors)
		
		# example
		#  all_ancestors = [['A11','A1','A'], ['B11','B1','B']]
		
		all_combinations = data.getAllCombinations(all_ancestors)
		
		# example:
		#  all_combinations = 
		#  [['A11', 'B11'],
		#   ['A1', 'B11'],
		#   ['A', 'B11'],
		#   ['A11', 'B1'],
		#   ['A1', 'B1'],
		#   ['A', 'B1'],
		#   ['A11', 'B'],
		#   ['A1', 'B'],
		#   ['A', 'B']]
		
		for combination in all_combinations:
			# Create a row of the A matrix. Each column is for a non_top_node
			# It contains 1 if that node exists in combination, else 0
			A_row = [int(node in combination) for node in non_top_nodes]
			# Add on a column at the end for constant C which is always there
			A_row.append(1)
			
			A_list.append(A_row)
			b_list.append(b_row)
			
	import numpy
	import numpy.linalg

	A = numpy.array(A_list)
	b = numpy.array(b_list)
	
	x, residues, rank, s = numpy.linalg.lstsq(A,b)
	
	group_values=dict()

	for i in range(len(non_top_nodes)):
		group_values[non_top_nodes[i]] = tuple(x[i,:])
	group_values['Constant'] = tuple(x[len(non_top_nodes),:])
	family.tree.children['Constant']=[]
	for node in top_nodes:
		group_values[node] = (0,0,0,0)
	
	def print_node_tree(node,indent=0):
		print (' '*indent +
				node.ljust(17-indent) + 
				"\t%.2g\t%.2g\t%.2g   \t%.2g"%group_values[node] )
		children = family.tree.children[node]
		if children:
			children.sort()
			for child in children:
				print_node_tree(child,indent+1)
	
	for node in top_nodes:
		print_node_tree(node)
	print_node_tree('Constant')


if False: # The following is Josh's old code:
	# Get set of all nodes
	nodeSet = family.tree.parent.keys()
	# remove those with no data
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
