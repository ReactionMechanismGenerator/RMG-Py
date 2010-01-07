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
RMG is an automatic chemical mechanism generator. It is awesomely awesome.
"""
import os.path

import time
import logging

from rmg.main import printRMGHeader, initializeLog
import rmg.settings as settings
import rmg.kinetics as kinetics
import rmg.unirxn.io as io
import rmg.reaction as reaction

################################################################################

def execute(inputFile, options):

	# Set directories
	settings.libraryDirectory = options.libraryDirectory

	# Set up log (uses stdout and a file)
	initializeLog(options.verbose)

	# Log start timestamp
	logging.info('RMG execution initiated at ' + time.asctime() + '\n')

	# Print out RMG header
	printRMGHeader()

	# Read input file
	network, Tlist, Plist, grainSize, numGrains, method, model = io.readInputFile(inputFile)

	# Shift network such that lowest-energy isomer has a ground state of 0.0
	logging.info('Zeroing lowest energy isomer...')
	network.shiftToZeroEnergy()

	# Determine energy grains
	logging.info('Determining energy grains...')
	Elist = network.determineEnergyGrains(grainSize, numGrains, max(Tlist))

	# Calculate density of states for all isomers in network
	logging.info('Calculating densities of states...')
	network.calculateDensitiesOfStates(Elist)

#	# DEBUG: Plot densities of states
#	import pylab
#	legend = []
#	for isomer in network.isomers:
#		if isomer.densStates is not None:
#			pylab.semilogy(Elist / 1000.0, isomer.densStates, '-')
#			#legend.append(str(isomer))
#	pylab.xlabel('Energy (kJ/mol)')
#	pylab.ylabel('Density of states (mol/J)')
#	#pylab.legend(legend, loc=4)
#	pylab.show()

	# Determine phenomenological rate coefficients
	logging.info('Calculating phenomenological rate coefficients...')
	K = network.calculateRateCoefficients(Tlist, Plist, Elist, method)

	# Create net reaction objects
	network.netReactions = []
	index = 0
	for i, reactantIsomer in enumerate(network.isomers):
		for j, productIsomer in enumerate(network.isomers[0:i]):
			netReaction = reaction.PDepReaction(reactantIsomer.species, productIsomer.species, network, None)
			netReaction.id = 'netReaction%i' % (index+1)
			index += 1
			network.netReactions.append(netReaction)
			# Fit interpolation model
			if model[0].lower() == 'chebyshev':
				modelType, degreeT, degreeP = model
				chebyshev = kinetics.ChebyshevKinetics()
				chebyshev.fitToData(Tlist, Plist, K[:,:,j,i], degreeT, degreeP)
				netReaction.kinetics = [chebyshev]
			elif model[0].lower() == 'pdeparrhenius':
				pass
			else:
				pass

	# Save results to file
	logging.info('Saving results...')
	outputFile = os.path.join(os.path.dirname(inputFile), 'output.xml')
	io.writeOutputFile(outputFile, network, Tlist, Plist, Elist, method, model)

	logging.info('')

	# Log end timestamp
	logging.info('RMG execution terminated at ' + time.asctime())

################################################################################

if __name__ == '__main__':

	# Command-line options
	description = 'RMG is an automatic chemical reaction mechanism ' + \
				  'generator that constructs kinetic models composed of ' + \
				  'elementary chemical reaction steps using a general ' + \
				  'understanding of how molecules react.'

	# Initialize command-line option parser
	import optparse
	parser = optparse.OptionParser(usage='usage: %prog [options] FILE',
								   version="RMG v0.0.1",
								   description=description)

	# Add options
	parser.add_option('-q', '--quiet',
					  action='store_const', const=30, default=20, dest='verbose',
					  help='quiet mode; only log errors and warnings')
	parser.add_option('-v', '--verbose',
					  action='store_const', const=10, default=20, dest='verbose',
					  help='verbose mode; log debug info')
	parser.add_option('-l', '--library-directory', default='',
					  action="store", type="string", dest="libraryDirectory",
					  help='use DIR as library directory', metavar='DIR')
	
	# Parse the command-line arguments
	options, args = parser.parse_args()

	# There should be exactly one positional argument: the input file
	# If this is not the case, print the help information and stop
	if len(args) != 1:
		parser.parse_args(['-h'])
		quit()

	execute(args[0], options)