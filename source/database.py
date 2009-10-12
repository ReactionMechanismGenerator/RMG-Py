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
import rmg.structure as structure
import rmg.reaction as reaction
import rmg.thermo as thermo

import logging

################################################################################

def loadThermoDatabases(databasePath):
	"""
	Create and load the thermodynamics databases.
	"""
	import os.path
	databasePath += '/'

	# Create and load thermo databases
	thermo.thermoDatabase = thermo.ThermoDatabaseSet()
	thermo.thermoDatabase.load(databasePath )

	# Create and load forbidden structures
	thermo.forbiddenStructures = data.Dictionary()
	thermo.forbiddenStructures.load(os.path.join(databasePath, 'forbiddenStructure.txt'))
	thermo.forbiddenStructures.toStructure()

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


################################################################################

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

def pruneKineticsTrees():

	# Prune kinetics dictionaries and trees
	for key, family in reaction.kineticsDatabase.families.iteritems():
		forwardTemplate, reverseTemplate = family.getTemplateLists()
		family.prune(forwardTemplate)

		# Print dictionary
		print '/////////////////////////////////////////////////////////////////////////////////'
		print '//'
		print '// Structure dictionary for %s' % family.label
		print '//'
		print '/////////////////////////////////////////////////////////////////////////////////'
		print ''
		for label, struct in family.dictionary.iteritems():
			if struct.__class__ == str:
				union = label + '\nUnion {' + ','.join(family.tree.children[label]) + '}\n'
				print union
			else:
				print struct.toAdjacencyList(label)

		# Print tree
		print '/////////////////////////////////////////////////////////////////////////////////'
		print '//'
		print '// Structure tree for %s' % family.label
		print '//'
		print '/////////////////////////////////////////////////////////////////////////////////'
		print ''
		print family.tree.write(family.tree.top)

################################################################################

def findCatchallNodes():

	import rmg.structure as structure

	families = { 'group': thermo.thermoDatabase.groupDatabase, \
				 'radical': thermo.thermoDatabase.radicalDatabase, \
				 '1,5-interations': thermo.thermoDatabase.int15Database, \
				 'gauche': thermo.thermoDatabase.gaucheDatabase, \
				 'ring': thermo.thermoDatabase.ringDatabase, \
				 'other': thermo.thermoDatabase.otherDatabase }

	#families = reaction.kineticsDatabase.families


	for key, family in families.iteritems():
		for parent, children in family.tree.children.iteritems():
			for child in children:
				try:
					print family.library[parent], family.library[child]
					if family.library[parent] == family.library[child]:
						print key, parent, child
				except KeyError:
					pass
#				struct1 = family.dictionary[parent]
#				struct2 = family.dictionary[child]
#				if isinstance(struct1, structure.Structure) and isinstance(struct2, structure.Structure):
#
#					struct1.updateAtomTypes()
#					struct2.updateAtomTypes()
#
#					match = True
#					if len(struct1.atoms()) != len(struct2.atoms()):
#						match = False
#					else:
#						for i in range(len(struct1.atoms())):
#							if struct1.atoms()[i].atomType != struct2.atoms()[i].atomType or \
#								struct1.atoms()[i].electronState != struct2.atoms()[i].electronState:
#								match = False

#					labeled1 = struct1.getLabeledAtoms().values()[0]
#					labeled2 = struct2.getLabeledAtoms().values()[0]
#					map21 = {labeled2: labeled1}
#					map12 = {labeled1: labeled2}
#					match, map21List, map12List = struct1.findSubgraphIsomorphisms(struct2, map12, map21)
#					match = False
#					for map in map21List:
#						found = True
#						for k, v in map.iteritems():
#							if k.atomType != v.atomType: found = False
#						if found: match = True
#
#					if match:
#						print key, parent, child
	print ''
		
################################################################################

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
		
		#: a dictionary of lists. key = node, value = list of kinetics items which contributed to that node
		kinetics_used_in={'Constant':[]}
		for node in node_set: # must initialise in loop so each has a separate list instance!
			kinetics_used_in[node] = list() 
		Ts = [300, 500, 1000, 1500]
		
		def rates_string(k):
			"""Return a string representing the rates of :class:`kinetics` object k
			
			log10 of the k at a bunch of T values"""
			string = "%5.2f "*len(Ts)
			return string%tuple([ math.log10(k.getRateConstant(T,Hrxn)) for T in Ts ])

		
		A_list = []
		b_list = []
		# Get available data
		
		to_delete=[]
		for key, kinetics in family.library.iteritems():
			if kinetics.alpha:
				logging.warning("Warning: %s %s has EP alpha = %g"%(kinetics.index, kinetics.label, kinetics.alpha))
				to_delete.append(key)
				
			#if re.search('O2b',kinetics.label): 
			#	logging.warning("Removing %s %s because I don't like O2b"%(kinetics.index, kinetics.label))
			#	to_delete.append(key)
			#	
		for key in to_delete:
			del family.library[key]
			logging.warning("Deleting %s from kinetics library!"%key)
				
				
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
				logging.warning("Warning: %s has EP alpha = %g"%(nodes,kinetics.alpha))
			
			Hrxn=0
			
			b_row = [ math.log10(kinetics.getRateConstant(T,Hrxn)) for T in Ts ]
				
			all_ancestors=list()
			kinetics.used_in_groups = list()
			kinetics.used_in_combinations = list()
			for node in nodes:
				# start with the most specific - the node itself
				# then add the ancestors
				ancestors = [node]
				ancestors.extend( family.tree.ancestors(node) )
				# append to the list of lists
				all_ancestors.append(ancestors)
				# add to the list 
				kinetics.used_in_groups.extend(ancestors)
				
				for ancestor in ancestors:
					kinetics_used_in[ancestor].append(kinetics)
			kinetics_used_in['Constant'].append(kinetics)
			
			# example
			#  all_ancestors = [['A11','A1','A'], ['B11','B1','B']]
			#  kinetics.used_in_groups = [ 'A11','A1','A','B11','B1','B' ]
			#  kinetics_used_in['A11'] = kinetics_used_in['A1'] ... = [... <kinetics>]
			
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
				
				kinetics.used_in_combinations.append(len(A_list))
				A_list.append(A_row)
				b_list.append(b_row)
				
		A = numpy.array(A_list)
		b = numpy.array(b_list)
		
		logging.info("Library contained %d rates"%len(family.library))
		logging.info("Matrix for inversion is %d x %d"%A.shape)
		
		x, residues, rank, s = numpy.linalg.lstsq(A,b)
		
		fitted_b = numpy.dot(A,x)
		errors = fitted_b - b
		#: squared and SUMMED over temperatures, not averaged
		errors_sum_squared = numpy.sum(errors*errors, axis=1)
		
		group_values=dict()
		group_error=dict()
		group_count=dict()
		group_error_MAD_by_T=dict()
		
		for node in top_nodes:
			group_values[node] = tuple([0 for i in Ts]) # eg. (0 0 0 0 0)
			group_error[node] = 0
			group_count[node] = 0
			group_error_MAD_by_T[node] = tuple([0 for i in Ts]) # eg. (0 0 0 0 0)
			
		for i in range(len(x)):
			group_values[group_names[i]] = tuple(x[i,:])
			
		for i in range(len(x)): # for each group
			#: vector of 1s and 0s, one for each rate-group
			rates_in_group = A[:,i]  
			#: number of data points training this group (each measured rate may be counted many times)
			group_count[group_names[i]] = sum(rates_in_group)
			#: RMS error for this group (where M = mean over temperatures and rates training the group) 
			group_error[group_names[i]] = numpy.sqrt(
				sum(rates_in_group * errors_sum_squared)  /
					 sum(rates_in_group) / len(Ts)   )
			#: Mean Absolute Deviation, reported by Temperature (as tuple)
			group_error_MAD_by_T[group_names[i]] = tuple( 
				numpy.dot(rates_in_group, abs(errors)) /
				 sum(rates_in_group)  )
				
		for key, kinetics in family.library.iteritems():
			rows = kinetics.used_in_combinations
			#: RMS error for this rate (where M = mean over temperatures and group combinations it's estimated by)
			kinetics.RMS_error = numpy.sqrt( 
				sum([errors_sum_squared[i] for i in rows])
				 / len(rows) / len(Ts) 
				)
			kinetics.key = key
		rates = family.library.values()
		rates.sort(cmp=lambda x,y: cmp(x.RMS_error, y.RMS_error))
		print "Rate expressions sorted by how well they are predicted by their group combinations"
		
		rates_1000 = []
		rates_err = []
		for k in rates:
			print "%-5s %-30s\tRMS error: %.2f  Rates: %s  %.30s"%(k.index, k.key, k.RMS_error, rates_string(k), k.comment )
			rates_1000.append( math.log10(k.getRateConstant(1000,Hrxn)) )
			rates_err.append( k.RMS_error )  # [Ts.index(T)]
		rates_1000 = numpy.array(rates_1000)
		rates_err = numpy.array(rates_err)
		
		fig_number = family_names.index(family_name)
		fig1 = pylab.figure( fig_number )
		pylab.plot(rates_1000, rates_err, 'o')
		pylab.xlabel('log10(k) at 1000K')
		pylab.ylabel('RMSE')
		pylab.show()
		
		def print_node_tree(node,indent=0):
			print (' '*indent +
					node.ljust(17-indent) + 
					("\t%7.2g"*len(group_values[node])) % group_values[node]  +
					"\t%6.2g\t%d"%(group_error[node],group_count[node]) + 
					("\t%7.3g"*len(group_error_MAD_by_T[node])) % group_error_MAD_by_T[node] 
				)
			children = family.tree.children[node]
			if children:
				children.sort()
				for child in children:
					# recurse!
					print_node_tree(child,indent+1)
					
		print ("Log10(k) at T=   " + ("\t%7g"*len(Ts)) % tuple(Ts) + 
				'\t RMS\tcount' + 
				("\tMAD @ %d"*len(Ts)) % tuple(Ts) 
			)
			
		print_node_tree('Constant')
		for node in top_nodes:
			print_node_tree(node)
		print
		
		
		fig = pylab.figure( 100 + fig_number )
		
		xvals = numpy.array([ group_count[group] for group in group_names ])
		yvals = numpy.array([ group_error[group] for group in group_names ])
		pylab.semilogx(xvals,yvals,'o',picker=5) # 5 points tolerance
		pylab.title(family_name)
		
		def onpick(event):
			thisline = event.artist
			xdata = thisline.get_xdata()
			ydata = thisline.get_ydata()
			for ind in event.ind:
				group_name = group_names[ind]
				print "#%d Name: %s \tRates:%d \tNode-Rates:%d \tRMS error: %g"%(ind, group_name, len(kinetics_used_in[group_name]) , xvals[ind], yvals[ind])
				print "MAD errors:"+("  %.2f"*len(Ts))%group_error_MAD_by_T[group_name]
				print "Kinetics taken from:"
				rates = kinetics_used_in[group_name]
				rates.sort(cmp=lambda x,y: cmp(x.RMS_error, y.RMS_error))
				for k in rates:
					print "%s\tIndex:%s \t%s "%(k.key,k.index,repr(k))
					print "RMS error: %.2f"%(k.RMS_error), 
					print "Rates: ",rates_string(k)
					for combo in k.used_in_combinations:
						#print "A[%d,%d] ="%(combo,ind),A[combo,ind]
						if not A[combo,ind]:
							#print "Rate didn't use the node in question (presumably used an ancestor)"
							continue
						print "Using",
						used_nodes = [ group_names[i] for i in A[combo,:].nonzero()[0] ]
						used_nodes.remove(group_name)
						print group_name + ' with ' + ' + '.join(used_nodes) + '\t',
						rms = numpy.sqrt( errors_sum_squared[combo] / len(Ts) )
						print "RMSE: %.2f  Err(T):"%(rms), errors[combo]
					print 
				#print 'check %g:'%ind, zip(xdata[ind], ydata[ind])
				
		connection_id = fig.canvas.mpl_connect('pick_event', onpick)
		# disconnect with: fig.canvas.mpl_disconnect(connection_id) 
		pylab.show()
		#http://matplotlib.sourceforge.net/users/event_handling.html
		
		import pdb; pdb.set_trace()
			
	#	graph = family.drawFullGraphOfTree()
	#	return graph

################################################################################

def saveDOTAndImage(graph, root, format='svg', prog='dot'):
	"""
	Save a DOT file and an image file from the Dot class `graph` to the location
	specified by `root`, a directory and filename without extension.
	"""

	#print '\tCreating DOT file...'
	f=open(root+'.dot','w')
	f.write(graph.to_string())
	f.close()

	#print '\tCreating SVG file...'
	filename=root+'.'+format
	if format=='svg':  # annoyingly, dot creates svg's without units on the font size attribute.
		st=graph.create_svg(prog=prog)
		st=re.sub(r"(font\-size\:[0-9]+\.*[0-9]*)([^p])",r"\1pt\2",st)
		f=open(filename,'w')
		f.write(st)
		f.close()
	else:
		graph.write(filename,format=format,prog=prog)


def drawAllTrees(root):
	"""
	Draws all of the trees in the thermo and kinetics database. The trees are
	placed in the folder specified by `root`.
	"""

	import os

	thermoDatabases = {'group': thermo.thermoDatabase.groupDatabase,
		'1,5-interactions': thermo.thermoDatabase.int15Database,
		'gauche': thermo.thermoDatabase.gaucheDatabase,
		'other': thermo.thermoDatabase.otherDatabase,
		'radical': thermo.thermoDatabase.radicalDatabase,
		'ring': thermo.thermoDatabase.ringDatabase}
	
	
	# Create directories
	try:
		os.makedirs(root)
		os.makedirs(root+os.sep+'thermo')
		os.makedirs(root+os.sep+'kinetics')
		for key in thermoDatabases:
			os.makedirs(root+os.sep+'thermo'+os.sep+key)
		for key in reaction.kineticsDatabase.families:
			os.makedirs(root+os.sep+'kinetics'+os.sep+key)

	except OSError:
		raise

	# Process thermo databases
	for key, thermoDatabase in thermoDatabases.iteritems():
		# Process all structures in dictionary
		for label, node in thermoDatabase.dictionary.iteritems():
			if isinstance(node, structure.Structure):
				print '\t'+ label
				graph = node.toDOT()
				graph.set('fontsize','10')
				saveDOTAndImage(graph, root+os.sep+'thermo'+os.sep+key+os.sep+label, 'svg', 'neato')
		# Process tree itself
		print '\t'+key
		graph = thermoDatabase.tree.toDOT()
		graph.set('fontsize','10')
		saveDOTAndImage(graph, root+os.sep+'thermo'+os.sep+key, 'svg', 'dot')
		print 'Created DOT for thermo database %s' % (key)

	# Process kinetics databases
	for key, family in reaction.kineticsDatabase.families.iteritems():
		# Process all structures in dictionary
		for label, node in family.dictionary.iteritems():
			if isinstance(node, structure.Structure):
				print '\t'+ label
				graph = node.toDOT()
				graph.set('fontsize','10')
				saveDOTAndImage(graph, root+os.sep+'thermo'+os.sep+key+os.sep+label, 'svg', 'neato')
		# Process tree itself
		print '\t'+key
		graph = family.tree.toDOT()
		graph.set('fontsize','10')
		saveDOTAndImage(graph, root+os.sep+'kinetics'+os.sep+key, 'svg', 'neato')
		print 'Created DOT for kinetics database %s' % (key)

#
#
#	import rmg.structure
#	s = rmg.structure.Structure()
#	s.fromSMILES('C=CC=CCC')
#	s.updateAtomTypes()
#	graph = s.toDOT()
#
#	key = 'hxd13'
#	print '\tCreating DOT file...'
#	f=open(key+'.dot','w')
#	f.write(graph.to_string())
#	f.close()
#
#	print '\tCreating SVG file...'
#	format='svg'
#	prog='neato'
#	filename=key+'.'+format
#	if format=='svg':  # annoyingly, dot creates svg's without units on the font size attribute.
#		st=graph.create_svg(prog=prog)
#		st=re.sub(r"(font\-size\:[0-9]+\.*[0-9]*)([^p])",r"\1pt\2",st)
#		f=open(filename,'w')
#		f.write(st)
#		f.close()
#	else:
#		graph.write(filename,format=format,prog=prog)
#
#	quit()


################################################################################

if __name__ == '__main__':

	import math
	# Show debug messages (as databases are loading)
	main.initializeLog(10)

	# Load databases
	databasePath = '../data/RMG_database'
#	loadThermoDatabases(databasePath)
	loadKineticsDatabases(databasePath,only_families=['H_Abstraction'])

	# findCatchallNodes()
	fit_groups()	
#	fit_groups(['H abstraction'])
#	graph = fit_groups()
#	write_xml()
	
#	for node in graph.get_node_list():	
#		node.set_style('filled')
#		node.set_fontcolor('#FFFFFFFF')
#		node.set_fillcolor('#000000FF')
	
	# Prune kinetics trees
	#pruneKineticsTrees()

	# Draw kinetics trees
	#drawKineticsTrees()

	# Draw trees in database
	drawAllTrees('database')

#	for key, family in reaction.kineticsDatabase.families.iteritems():
#
#		print '\tCreating DOT object...'
#
#		graph = family.tree.toDOT()
#
#		graph.set('fontsize','10')
#		format='svg'
#		prog='dot'
#
#		print '\tCreating DOT file...'
#		f=open(key+'.dot','w')
#		f.write(graph.to_string())
#		f.close()
#
#		print '\tCreating SVG file...'
#		filename=key+'.'+format
#		if format=='svg':  # annoyingly, dot creates svg's without units on the font size attribute.
#			st=graph.create_svg(prog=prog)
#			st=re.sub(r"(font\-size\:[0-9]+\.*[0-9]*)([^p])",r"\1pt\2",st)
#			f=open(filename,'w')
#			f.write(st)
#			f.close()
#		else:
#			graph.write(filename,format=format,prog=prog)
#
#		print 'Created DOT for reaction family %s' % (key)

