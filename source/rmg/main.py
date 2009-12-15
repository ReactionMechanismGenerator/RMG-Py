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
import settings
import species
import reaction
import unirxn.network

################################################################################

def execute(inputFile, options):
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
	settings.outputDirectory = options.outputDirectory
	settings.scratchDirectory = options.scratchDirectory
	settings.libraryDirectory = options.libraryDirectory

	# Save initialization time
	settings.initializationTime = time.time()

	# Set up log (uses stdout and a file)
	initializeLog(options.verbose)
	
	# Log start timestamp
	logging.info('RMG execution initiated at ' + time.asctime() + '\n')
	
	# Print out RMG header
	printRMGHeader()
	
	# Make output subdirectories
	plotDir = os.path.join(settings.outputDirectory, 'plot')
	if os.path.exists(plotDir):
		for f in os.listdir(plotDir):
			os.remove(plotDir + '/' + f)
		os.rmdir(plotDir)
	os.mkdir(plotDir)

	specDir = os.path.join(settings.outputDirectory, 'species')
	if os.path.exists(specDir):
		for f in os.listdir(specDir):
			os.remove(specDir + '/' + f)
		os.rmdir(specDir)
	os.mkdir(specDir)

	# Read input file
	reactionModel, coreSpecies, reactionSystems = io.readInputFile(inputFile)
	
	#reactionModel.loadSeedMechanism('/files/research/rmgpy/devel/data/GRI-Mech3.0')
	
	# Initialize reaction model
	if options.restart:
		import cPickle
		import ctml_writer
		logging.info('Loading previous restart file...')
		f = open(os.path.join(settings.outputDirectory,'restart.pkl'), 'rb')
		species.speciesList, species.speciesCounter, reaction.reactionList, \
			reactionModel, reactionSystems = cPickle.load(f)
		f.close()
		# Cantera stuff
		reload(ctml_writer) # ensure new empty ctml_writer._species and ._reactions lists
		for reactor in reactionSystems:
			# initialise the ctml_writer thing
			reactor.initializeCantera()
		for spec in reactionModel.core.species:
			# add species to ctml_writer._species list
			spec.toCantera() 
		for rxn in reactionModel.core.reactions:
			# add reaction to ctml_writer._reactions list
			rxn.toCantera()
		#print "enter 'c' to continue"; import pdb; pdb.set_trace()
		options.restart = False # have already restarted
	else:
		reactionModel.initialize(coreSpecies)

	# RMG execution statistics
	coreSpeciesCount = []
	coreReactionCount = []
	edgeSpeciesCount = []
	edgeReactionCount = []
	execTime = []
	restartSize = []
	memoryUse = []

	# Handle unimolecular (pressure dependent) reaction networks
	if settings.unimolecularReactionNetworks:
		reactionModel.updateUnimolecularReactionNetworks()
		logging.info('')

	# Main RMG loop
	done = False
	while not done:

		done = True
		objectsToEnlarge = []
		for index, reactionSystem in enumerate(reactionSystems):
			
			# Conduct simulation
			logging.info('Conducting simulation of reaction system %s...' % (index+1))
			t, y, dydt, valid, obj = reactionSystem.simulate(reactionModel)

			# Postprocess results
			logging.info('')
			logging.info('Saving simulation results for reaction system %s...' % (index+1))
			reactionSystem.postprocess(reactionModel, t, y, dydt, str(index+1))

			# If simulation is invalid, note which species should be added to
			# the core
			if not valid:
				objectsToEnlarge.append(obj)
				done = False

		# Enlarge objects identified by the simulation for enlarging
		# These should be Species or Network objects
		logging.info('')
		objectsToEnlarge = list(set(objectsToEnlarge))
		for object in objectsToEnlarge:
			reactionModel.enlarge(object)
		
		# Handle unimolecular (pressure dependent) reaction networks
		if settings.unimolecularReactionNetworks:
			reactionModel.updateUnimolecularReactionNetworks()
			logging.info('')

		# Save the restart file
		# In order to get all the references preserved, you must pickle all of
		# the objects in one concerted dump; this also has the added benefits
		# of using less space and running faster
		import cPickle
		logging.info('Saving restart file...')
		f = open(os.path.join(settings.outputDirectory,'restart.pkl'), 'wb')
		cPickle.dump((
			species.speciesList,
			species.speciesCounter,
			reaction.reactionList,
			reactionModel,
			reactionSystems),
			f)
		f.close()

		if not done:
			logging.info('Updating RMG execution statistics...')
			# Update RMG execution statistics
			coreSpeciesCount.append(len(reactionModel.core.species))
			coreReactionCount.append(len(reactionModel.core.reactions))
			edgeSpeciesCount.append(len(reactionModel.edge.species))
			edgeReactionCount.append(len(reactionModel.edge.reactions))
			execTime.append(time.time() - settings.initializationTime)
			from guppy import hpy
			hp = hpy()
			memoryUse.append(hp.heap().size / 1.0e6)
			logging.debug('Execution time: %s s' % (execTime[-1]))
			logging.debug('Memory used: %s MB' % (memoryUse[-1]))
			restartSize.append(os.path.getsize(os.path.join(settings.outputDirectory,'restart.pkl')) / 1.0e6)
			generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)

		logging.info('')


	# Write output file
	logging.info('MODEL GENERATION COMPLETED')
	logging.info('')
	logging.info('The final model core has %s species and %s reactions' % (len(reactionModel.core.species), len(reactionModel.core.reactions)))
	logging.info('The final model edge has %s species and %s reactions' % (len(reactionModel.edge.species), len(reactionModel.edge.reactions)))
	io.writeOutputFile(os.path.join(settings.outputDirectory,'output.xml'), reactionModel, reactionSystems)

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
	
	# create file handler
	log_file_name = os.path.join(settings.outputDirectory,'RMG.log')
	if os.path.exists(log_file_name):
		backup_name = '_backup'.join(os.path.splitext(log_file_name))
		if os.path.exists(backup_name):
			print "Removing old file %s"%backup_name
			os.remove(backup_name)
		print "Renaming %s to %s"%(log_file_name, backup_name)
		os.rename(log_file_name, backup_name)
	fh = logging.FileHandler(filename=log_file_name) #, backupCount=3)
	fh.setLevel(logging.DEBUG) # always verbose in the file
	fh.setFormatter(formatter)
	# notice that STDERR does not get saved to the log file
	# so errors from underlying libraries (eg. openbabel) etc. that report
	# on stderr will not be logged to disk.
	
	# remove old handlers!
	while logger.handlers:
		logger.removeHandler(logger.handlers[0])
	
	# Add ch to logger
	logger.addHandler(ch)
	logger.addHandler(fh)
	
################################################################################

def generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount, \
	edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize):
	"""
	Generate a number of plots describing the statistics of the RMG job, 
	including the reaction model core and edge size and memory use versus
	execution time. These will be placed in the output directory in the plot/
	folder.
	"""

	# Only generate plots if that flag is turned on (in input file)
	if not settings.generatePlots:
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
	plt.savefig(os.path.join(settings.outputDirectory, 'plot/coreSize.svg'))
	plt.clf()

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.loglog(execTime, edgeSpeciesCount, 'o-b')
	ax1.set_xlabel('Execution time (s)')
	ax1.set_ylabel('Number of edge species')
	ax2 = ax1.twinx()
	ax2.loglog(execTime, edgeReactionCount, 'o-r')
	ax2.set_ylabel('Number of edge reactions')
	plt.savefig(os.path.join(settings.outputDirectory, 'plot/edgeSize.svg'))
	plt.clf()

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.semilogx(execTime, memoryUse, 'o-k')
	ax1.semilogx(execTime, restartSize, 'o-g')
	ax1.set_xlabel('Execution time (s)')
	ax1.set_ylabel('Memory (MB)')
	ax1.legend(['RAM', 'Restart file'], loc=2)
	plt.savefig(os.path.join(settings.outputDirectory, 'plot/memoryUse.svg'))
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
	