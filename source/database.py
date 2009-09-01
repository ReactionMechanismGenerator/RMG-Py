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
	species.thermoDatabase.load(databasePath )

	# Create and load forbidden structures
	species.forbiddenStructures = data.Dictionary()
	species.forbiddenStructures.load(os.path.join(databasePath, 'forbiddenStructure.txt'))
	species.forbiddenStructures.toStructure()

def loadKineticsDatabases(databasePath, only_families=False):
	"""
	Create and load the kinetics databases (reaction families).
	If only_families is a list like ['H_Abstraction'] then only families in this
	list will be loaded.
	"""
	reaction.kineticsDatabase = reaction.ReactionFamilySet()
	reaction.kineticsDatabase.load(databasePath, only_families=only_families)

################################################################################

def drawKineticsTrees():
	"""
	For each reaction family, output a DOT file containing a combined form of
	the various trees that make up the reaction family.
	"""
	import pydot
	
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
		
		
###################
def fit_groups(family_names = None):
	"""Decouples a nested tree and fits values to groups for each seperate tree.
	   If given a list of family names, only does those families.
	"""
	import os
	import math
	import numpy
	import numpy.linalg
	import pylab
	
	if not family_names: 
		family_names = reaction.kineticsDatabase.families.keys()
		
	for family_name in family_names:
		family = reaction.kineticsDatabase.families[family_name]
		print 
		if not family.library:
			logging.debug("Family '%s' has no data in the library."%family_name)
			if family.reverse.library:
				logging.debug("(but its reverse '%s' does)"%family.reverse.label)
			continue
		
		logging.info("Fitting groups for reaction family: %s (%s)"%(family_name,
			os.path.basename(os.path.abspath(family._path))) )
		
		# Get set of all nodes
		node_set = family.tree.parent.keys()
		non_top_nodes = [ node for node in node_set if family.tree.parent[node] ]
		top_nodes = [ node for node in node_set if not family.tree.parent[node] ]
		group_names = [node for node in non_top_nodes] # poor man's copy
		group_names.append("Constant")
		family.tree.children['Constant']=[]
		
		A_list = []
		b_list = []
		# Get available data
		for key, kinetics in family.library.iteritems():
			nodes = key.split(';')
			# example:
			#  nodes = ['A11', 'B11']
			#  kinetics = <rmg.reaction.ArrheniusEPKinetics instance>
			#b_row = [ math.log(kinetics.A),
			#		  kinetics.n,
			#		  kinetics.alpha,
			#		  kinetics.E0 ]
			if kinetics.alpha:
				logging.info("Warning: %s has EP alpha = %g"%(nodes,kinetics.alpha))
			Ts = [300, 500, 1000, 1500]
			Hrxn=0
			b_row = [ math.log10(kinetics.getRateConstant(T,Hrxn)) for T in Ts ]
				
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
			#  [['A11', 'B11'], ['A1', 'B11'], ['A', 'B11'],  ['A11', 'B1'],
			#   ['A1', 'B1'], ['A', 'B1'],  ['A11', 'B'], ['A1', 'B'], ['A', 'B']]
			
			for combination in all_combinations:
				# Create a row of the A matrix. Each column is for a non_top_node
				# It contains 1 if that node exists in combination, else 0
				A_row = [int(node in combination) for node in non_top_nodes]
				# Add on a column at the end for constant C which is always there
				A_row.append(1)
				
				A_list.append(A_row)
				b_list.append(b_row)
				
		A = numpy.array(A_list)
		b = numpy.array(b_list)
		
		logging.info("Library contained %d rates"%len(family.library))
		logging.info("Matrix for inversion is %d x %d"%A.shape)
		
		x, residues, rank, s = numpy.linalg.lstsq(A,b)
		
		fitted_b = numpy.dot(A,x)
		errors = fitted_b - b
		errors_sum_squared = numpy.sum(errors*errors, axis=1)
		
		
		group_values=dict()
		group_error=dict()
		group_count=dict()
		group_error_MAD_by_T=dict()
		
		for node in top_nodes:
			group_values[node] = (0,0,0,0)
			group_error[node] = 0
			group_count[node] = 0
			group_error_MAD_by_T[node] = (0,0,0,0)
			
		for i in range(len(x)):
			group_values[group_names[i]] = tuple(x[i,:])
			
		
		for i in range(len(x)): # for each group
			rates_in_group = A[:,i] 
			group_error[group_names[i]] = numpy.sqrt(
				sum(rates_in_group * errors_sum_squared)  /
					 sum(rates_in_group) / len(Ts)   )
			group_count[group_names[i]] = sum(rates_in_group)
			group_error_MAD_by_T[group_names[i]] = tuple( 
				numpy.dot(rates_in_group, abs(errors)) /
				 sum(rates_in_group)  )
		 #   if group_names[i]=='Cs_rad':
		 #   	print "foobar"
			
		
		def print_node_tree(node,indent=0):
			print (' '*indent +
					node.ljust(17-indent) + 
					"\t%7.2g\t%7.2g\t%7.2g\t%7.2g"%group_values[node]  +
					"\t%6.2g\t%d"%(group_error[node],group_count[node]) + 
					"\t%7.3g\t%7.3g\t%7.3g\t%7.3g"%group_error_MAD_by_T[node]
				)
			children = family.tree.children[node]
			if children:
				children.sort()
				for child in children:
					print_node_tree(child,indent+1)
					
		print ("Log10(k) at T=   \t%7g\t%7g\t%7g\t%7g"%tuple(Ts) + 
				'\t RMS\tcount' + 
				"\tMAD @ %d\tMAD @ %d\tMAD @ %d\tMAD @ %d"%tuple(Ts) 
			)
		print_node_tree('Constant')
		for node in top_nodes:
			print_node_tree(node)
		print
		
		pylab.figure(1)
		pylab.semilogx(rates_in_group)
		#http://matplotlib.sourceforge.net/users/event_handling.html
			
	#	graph = family.drawFullGraphOfTree()
	#	return graph

def write_xml(family_names = None):
	"""Writes the library to xml files
	"""
	import os
	import xml.dom.minidom
	
	# Create document
	dom = xml.dom.minidom.Document()

	# Process root element
	root = dom.createElement('rmgoutput')
	dom.appendChild(root)
	
	
	if not family_names: 
		family_names = reaction.kineticsDatabase.families.keys()
		
	for family_name in family_names:
		family = reaction.kineticsDatabase.families[family_name]
		print 
		if not family.library:
			logging.debug("Family '%s' has no data in the library."%family_name)
			if family.reverse.library:
				logging.debug("(but its reverse '%s' does)"%family.reverse.label)
			continue
		
		logging.info("Writing xml for reaction family: %s (%s)"%(family_name,
			os.path.basename(os.path.abspath(family._path))) )
			
		
		
		family.library.toXML(dom,root)
		print dom.toprettyxml()

################################################################################

if __name__ == '__main__':
	import math
	# Show debug messages (as databases are loading)
	main.initializeLog(10)

	# Load databases
	databasePath = '../data/RMG_database'
	#loadThermoDatabases(databasePath)
	loadKineticsDatabases(databasePath,only_families=['H_Abstraction'])

#	fit_groups(['H abstraction'])	
	#graph = fit_groups()
	
	write_xml()
	
#	for node in graph.get_node_list():	
#		node.set_style('filled')
#		node.set_fontcolor('#FFFFFFFF')
#		node.set_fillcolor('#000000FF')
	
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
