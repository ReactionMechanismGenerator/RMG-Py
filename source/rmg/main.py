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
Contains the main function for RMG execution and several helper functions that
support it.
"""

import os
import time
import logging
import io
import sys
import numpy
import pylab
import pybel

import model

################################################################################

def execute(inputFile, outputDir, scratchDir, libraryDir, verbose):
	"""
	Generate a reaction model for the set of reaction systems specified in an 
	input file at location `inputFile`. Output and temporary files will be 
	placed in the directories `outputDir` and `scratchDir`, respectively; if
	either of these are empty strings, the corresponding files will be placed
	in the same directory as the input file. The `libraryDir` parameter can be
	used to specify additional directories to search for RMG databases. The
	`verbose` parameter is an integer specifying the amount of log text seen
	at the console; the levels correspond to those of the :data:`logging` module.
	"""

	# Set up log (uses stdout)
	initializeLog(verbose)
	
	# Log start timestamp
	logging.info('RMG execution initiated at ' + time.asctime() + '\n')
	
	# Print out RMG header
	printRMGHeader()
	
	# Make output subdirectories
	plotDir = outputDir + os.sep + 'plot'
	if os.path.exists(plotDir):
		for f in os.listdir(plotDir):
			os.remove(plotDir + '/' + f)
		os.rmdir(plotDir)
	os.mkdir(plotDir)

	specDir = outputDir + os.sep + 'species'
	if os.path.exists(specDir):
		for f in os.listdir(specDir):
			os.remove(specDir + '/' + f)
		os.rmdir(specDir)
	os.mkdir(specDir)

	# Read input file
	reactionModel, coreSpecies, reactionSystems = io.readInputFile(inputFile)
	
	# Initialize reaction model
	reactionModel.initialize(coreSpecies)

	# Main RMG loop
	done = False
	while not done:

		done = True
		for index, reactionSystem in enumerate(reactionSystems):
			
			# Conduct simulation
			t, y, valid, species = reactionSystem.simulate(reactionModel)

			y0 = numpy.zeros((len(t), len(y[0])), float)
			for i, u in enumerate(y):
				for j, v in enumerate(u):
					y0[i,j] = v

			legend = []
			for spec in reactionModel.core.species:
				legend.append(str(spec))

			# Make plot and save to file
			pylab.semilogx(t[1:], y0[1:,0])
			pylab.xlabel('Time (s)')
			pylab.ylabel('Pressure (Pa)')
			pylab.title('Pressure profile for reaction system #' + str(index+1))
			pylab.savefig(outputDir + '/plot/pressureProfile' + str(index+1) + '.svg')
			pylab.clf()

			# Make plot and save to file
			pylab.semilogx(t[1:], y0[1:,1])
			pylab.xlabel('Time (s)')
			pylab.ylabel('Volume (m^3)')
			pylab.title('Volume profile for reaction system #' + str(index+1))
			pylab.savefig(outputDir + '/plot/volumeProfile' + str(index+1) + '.svg')
			pylab.clf()

			# Make plot and save to file
			pylab.semilogx(t[1:], y0[1:,2])
			pylab.xlabel('Time (s)')
			pylab.ylabel('Temperature (K)')
			pylab.title('Temperature profile for reaction system #' + str(index+1))
			pylab.savefig(outputDir + '/plot/temperatureProfile' + str(index+1) + '.svg')
			pylab.clf()

			# Make plot and save to file
			pylab.loglog(t[1:], y0[1:,3:])
			pylab.xlabel('Time (s)')
			pylab.ylabel('Concentration (mol/m^3)')
			pylab.title('Concentration profiles for reaction system #' + str(index+1))
			pylab.legend(legend)
			pylab.savefig(outputDir + '/plot/concentrationProfile' + str(index+1) + '.svg')
			pylab.clf()

			# Draw species in core
			#for spec in reactionModel.core.species:
			#	mol = pybel.Molecule(spec.toOBMol())
			#	mol.draw(False, outputDir + '/species/' + str(spec) + '.svg')

			# Enlarge reaction model if simulation is invalid
			if not valid:
				reactionModel.enlarge(species)
				done = False


	# Write output file
	io.writeOutputFile(outputDir + '/output.xml', reactionModel, reactionSystems)

	# Log end timestamp
	logging.info('\nRMG execution terminated at ' + time.asctime())
	
################################################################################

def initializeLog(verbose):
	"""
	Set up a logger for RMG to use to print output to stdout. The
	`verbose` parameter is an integer specifying the amount of log text seen
	at the console; the levels correspond to those of the :data:`logging` module.
	"""

	# Create logger
	logger = logging.getLogger()
	logger.setLevel(verbose)
	
	# Create console handler and set level to debug; send everything to stdout
	# rather than stderr
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(verbose)
	
	# Create formatter and add to console handler
	#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
	formatter = logging.Formatter('%(message)s', '%Y-%m-%d %H:%M:%S')
	ch.setFormatter(formatter)
	
	# Add ch to logger
	logger.addHandler(ch)
	
################################################################################

def printRMGHeader():
	"""
	Output a header containing identifying information about RMG to the log.
	"""

	logging.info('################################################################')
	logging.info('#                                                              #')
	logging.info('#              RMG - Reaction Mechanism Generator              #')
	logging.info('#                    Python Version 0.0.1                      #')
	logging.info('#                         14 May 2009                          #')
	logging.info('#                                                              #')
	logging.info('#                 http://rmg.sourceforge.net/                  #')
	logging.info('#                                                              #')
	logging.info('################################################################\n')
	