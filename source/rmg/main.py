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
import os.path
import time
import logging
import io
import sys

import constants

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

	# Set directories
	constants.outputDir = outputDir
	constants.scratchDir = scratchDir
	constants.libraryDir = libraryDir

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

	# RMG execution statistics
	coreSpeciesCount = []
	coreReactionCount = []
	edgeSpeciesCount = []
	edgeReactionCount = []
	execTime = []
	restartSize = []
	memoryUse = []

	# Main RMG loop
	done = False
	while not done:

		done = True
		speciesToAdd = []
		for index, reactionSystem in enumerate(reactionSystems):
			
			# Conduct simulation
			logging.info('Conducting simulation of reaction system %s...' % (index+1))
			t, y, dydt, valid, species = reactionSystem.simulate(reactionModel)

			# Postprocess results
			logging.info('Saving simulation results for reaction system %s...' % (index+1))
			reactionSystem.postprocess(reactionModel, t, y, dydt, str(index+1))

			# If simulation is invalid, note which species should be added to
			# the core
			if not valid:
				speciesToAdd.append(species)
				done = False

		# Add the notes species to the core
		speciesToAdd = list(set(speciesToAdd))
		for species in speciesToAdd:
			reactionModel.enlarge(species)

		# Save the restart file
		logging.info('Saving restart file...')
		import pickle
		f = open(outputDir + '/restart.pkl', 'wb')
		pickle.dump(reactionModel, f)
		pickle.dump(reactionSystems, f)
		f.close()

		if not done:
			logging.info('Updating RMG execution statistics...')
			# Update RMG execution statistics
			coreSpeciesCount.append(len(reactionModel.core.species))
			coreReactionCount.append(len(reactionModel.core.reactions))
			edgeSpeciesCount.append(len(reactionModel.edge.species))
			edgeReactionCount.append(len(reactionModel.edge.reactions))
			execTime.append(time.clock())
			from guppy import hpy
			hp = hpy()
			memoryUse.append(hp.heap().size / 1.0e6)
			restartSize.append(os.path.getsize(outputDir + '/restart.pkl') / 1.0e6)
			generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)

		logging.info('')


	# Write output file
	logging.info('MODEL GENERATION COMPLETED')
	logging.info('')
	logging.info('The final model core has %s species and %s reactions' % (len(reactionModel.core.species), len(reactionModel.core.reactions)))
	logging.info('The final model edge has %s species and %s reactions' % (len(reactionModel.edge.species), len(reactionModel.edge.reactions)))
	io.writeOutputFile(outputDir + '/output.xml', reactionModel, reactionSystems)

	# Log end timestamp
	logging.info('RMG execution terminated at ' + time.asctime())
	
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
	
	# remove old handlers!
	while logger.handlers:
		logger.removeHandler(logger.handlers[0])
	
	# Add ch to logger
	logger.addHandler(ch)
	
################################################################################

def generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize):

	if not constants.generatePlots:
		return

	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.semilogx(execTime, coreSpeciesCount, 'o-b')
	ax1.set_xlabel('Execution time (s)')
	ax1.set_ylabel('Number of core species')
	ax2 = ax1.twinx()
	ax2.semilogx(execTime, coreReactionCount, 'o-r')
	ax2.set_ylabel('Number of core reactions')
	plt.savefig(constants.outputDir + '/plot/coreSize.svg')
	plt.clf()

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.loglog(execTime, edgeSpeciesCount, 'o-b')
	ax1.set_xlabel('Execution time (s)')
	ax1.set_ylabel('Number of edge species')
	ax2 = ax1.twinx()
	ax2.loglog(execTime, edgeReactionCount, 'o-r')
	ax2.set_ylabel('Number of edge reactions')
	plt.savefig(constants.outputDir + '/plot/edgeSize.svg')
	plt.clf()

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.semilogx(execTime, memoryUse, 'o-k')
	ax1.semilogx(execTime, restartSize, 'o-g')
	ax1.set_xlabel('Execution time (s)')
	ax1.set_ylabel('Memory (MB)')
	ax1.legend(['RAM', 'Restart file'], loc=2)
	plt.savefig(constants.outputDir + '/plot/memoryUse.svg')
	plt.clf()

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
	