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

import optparse
import rmg

################################################################################

if __name__ == '__main__':

	# Command-line options
	description = 'RMG is an automatic chemical reaction mechanism ' + \
	              'generator that constructs kinetic models composed of ' + \
				  'elementary chemical reaction steps using a general ' + \
				  'understanding of how molecules react.'
	
	# Initialize command-line option parser
	parser = optparse.OptionParser(usage='usage: %prog [options] FILE', \
	                               version="RMG v0.0.1", \
								   description=description)
	
	# Add options
	parser.add_option('-o', '--output-directory', default='', \
		              action="store", type="string", dest="outputDirectory", \
					  help='use DIR as output directory', metavar='DIR')
	parser.add_option('-s', '--scratch-directory', default='', \
	                  action="store", type="string", dest="scratchDirectory", \
					  help='use DIR as scratch directory', metavar='DIR')
	parser.add_option('-l', '--library-directory', default='', \
	                  action="store", type="string", dest="libraryDirectory", \
					  help='use DIR as library directory', metavar='DIR')
	parser.add_option('-v', '--verbose', default=0, \
	                  action="store", type="int", dest="verbose", \
					  help='set verbosity level to VALUE', metavar='VALUE')
	
	# Parse the command-line arguments
	options, args = parser.parse_args()
	
	# There should be exactly one positional argument: the input file
	# If this is not the case, print the help information and stop
	if len(args) != 1:
		parser.parse_args(['-h'])
		quit()
	
	# Execute RMG
	#rmg.core.execute(args[0], options.outputDirectory, \
	                 #options.scratchDirectory, options.libraryDirectory, \
					 #options.verbose)

