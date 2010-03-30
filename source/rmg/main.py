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
import rmg.log as logging
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

	# Set wall time
	if options.wallTime == '0': settings.wallTime = 0
	else:
		try:
			data = options.wallTime.split(':')
			if len(data) == 1:
				settings.wallTime = int(data[-1])
			elif len(data) == 2:
				settings.wallTime = int(data[-1]) + 60 * int(data[-2])
			elif len(data) == 3:
				settings.wallTime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3])
			elif len(data) == 4:
				settings.wallTime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3]) + 86400 * int(data[-4])
			else:
				raise ValueError('Invalid format for wall time; should be HH:MM:SS.')
		except ValueError, e:
			raise ValueError('Invalid format for wall time; should be HH:MM:SS.')

	# Save initialization time
	settings.initializationTime = time.time()

	# Set up log (uses stdout and a file)
	logging.initialize(options.verbose, os.path.join(settings.outputDirectory,'RMG.log'))

	# Log start timestamp
	logging.info('RMG execution initiated at ' + time.asctime() + '\n')
	
	# Print out RMG header
	logging.logHeader()
	
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
	
	# Initialize reaction model
	if options.restart:
		import gzip
		import cPickle
		import ctml_writer
		logging.info('Loading previous restart file...')
		f = gzip.GzipFile(os.path.join(settings.outputDirectory,'restart.pkl'), 'rb')
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

		if not done:
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
			# We also compress the restart file to save space (and lower the
			# disk read/write time)
			import gzip
			import cPickle
			logging.info('Saving restart file...')
			f = gzip.GzipFile(os.path.join(settings.outputDirectory,'restart.pkl'), 'wb')
			cPickle.dump((
				species.speciesList,
				species.speciesCounter,
				reaction.reactionList,
				reactionModel,
				reactionSystems),
				f)
			f.close()

			# Update RMG execution statistics
			logging.info('Updating RMG execution statistics...')
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
			saveExecutionStatistics(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)
			generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)

		logging.info('')
		
		# Consider stopping gracefully if the next iteration might take us
		# past the wall time
		if settings.wallTime > 0 and len(execTime) > 1:
			t = execTime[-1]
			dt = execTime[-1] - execTime[-2]
			if t + 2 * dt > settings.wallTime:
				logging.info('MODEL GENERATION TERMINATED')
				logging.info('')
				logging.info('There is not enough time to complete the next iteration before the wall time is reached.')
				logging.info('The output model may be incomplete.')
				logging.info('')
				logging.info('The current model core has %s species and %s reactions' % (len(reactionModel.core.species), len(reactionModel.core.reactions)))
				logging.info('The current model edge has %s species and %s reactions' % (len(reactionModel.edge.species), len(reactionModel.edge.reactions)))
				io.writeOutputFile(os.path.join(settings.outputDirectory,'output.xml'), reactionModel, reactionSystems)
				return

	# Write output file
	logging.info('MODEL GENERATION COMPLETED')
	logging.info('')
	logging.info('The final model core has %s species and %s reactions' % (len(reactionModel.core.species), len(reactionModel.core.reactions)))
	logging.info('The final model edge has %s species and %s reactions' % (len(reactionModel.edge.species), len(reactionModel.edge.reactions)))
	io.writeOutputFile(os.path.join(settings.outputDirectory,'output.xml'), reactionModel, reactionSystems)

	# Log end timestamp
	logging.info('')
	logging.info('RMG execution terminated at ' + time.asctime())
	
################################################################################

def saveExecutionStatistics(execTime, coreSpeciesCount, coreReactionCount, \
	edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize):
	"""
	Save the statistics of the RMG job to an Excel spreadsheet for easy viewing
	after the run is complete. The statistics are saved to the file
	`statistics.xls` in the output directory. The ``xlwt`` package is used to
	create the spreadsheet file; if this package is not installed, no file is
	saved.
	"""

	# Attempt to import the xlwt package; return if not installed
	try:
		import xlwt
	except ImportError:
		logging.warning('Package xlwt not found. Unable to save execution statistics.')
		return

	# Create workbook and sheet for statistics to be places
	workbook = xlwt.Workbook()
	sheet = workbook.add_sheet('Statistics')

	# First column is execution time
	sheet.write(0,0,'Execution time (s)')
	for i, etime in enumerate(execTime):
		sheet.write(i+1,0,etime)

	# Second column is number of core species
	sheet.write(0,1,'Core species')
	for i, count in enumerate(coreSpeciesCount):
		sheet.write(i+1,1,count)

	# Third column is number of core reactions
	sheet.write(0,2,'Core reactions')
	for i, count in enumerate(coreReactionCount):
		sheet.write(i+1,2,count)

	# Fourth column is number of edge species
	sheet.write(0,3,'Edge species')
	for i, count in enumerate(edgeSpeciesCount):
		sheet.write(i+1,3,count)

	# Fifth column is number of edge reactions
	sheet.write(0,4,'Edge reactions')
	for i, count in enumerate(edgeReactionCount):
		sheet.write(i+1,4,count)

	# Sixth column is memory used
	sheet.write(0,5,'Memory used (MB)')
	for i, memory in enumerate(memoryUse):
		sheet.write(i+1,5,memory)

	# Seventh column is restart file size
	sheet.write(0,6,'Restart file size (MB)')
	for i, memory in enumerate(restartSize):
		sheet.write(i+1,6,memory)

	# Save workbook to file
	fstr = os.path.join(settings.outputDirectory, 'statistics.xls')
	workbook.save(fstr)

################################################################################

def generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount,
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

