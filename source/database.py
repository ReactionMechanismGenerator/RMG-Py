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

	# Prune kinetics dictionaries and trees
	for key, family in reaction.kineticsDatabase.families.iteritems():
		forwardTemplate, reverseTemplate = family.getTemplateLists()
		family.prune(forwardTemplate)
		
	# Draw kinetics trees
	#drawKineticsTrees()

