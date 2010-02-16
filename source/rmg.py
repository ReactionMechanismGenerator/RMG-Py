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

################################################################################

def breakpoint():
	"""To insert a breakpoint in the code, just call breakpoint()"""
	import pdb; pdb.set_trace()
	
class tee:
	"""A simple tee to create a stream which prints to many streams"""
	def __init__(self, *fileobjects):
		self.fileobjects=fileobjects
	def write(self, string):
		for fileobject in self.fileobjects:
			fileobject.write(string)


def process_stats(stats_file, log_file):
		out_stream = tee(sys.stdout,open(log_file,'a')) # print to screen AND append to RMG.log
		stats = pstats.Stats(stats_file,stream=out_stream)
		stats.strip_dirs()
		print >>out_stream, "Sorted by internal time"
		stats.sort_stats('time')
		stats.print_stats(25)
		stats.print_callers(25)
		print >>out_stream, "Sorted by cumulative time"
		stats.sort_stats('cumulative')
		stats.print_stats(25)
		stats.print_callers(25)
		stats.print_callees(25)

if __name__ == '__main__':

	from guppy import hpy
	hp = hpy()
	hp.heap()

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
					  action='store_const', const=15, default=20, dest='verbose',
					  help='verbose mode; add extra details to log')
	parser.add_option('-d', '--debug',
					  action='store_const', const=10, default=20, dest='verbose',
					  help='debug mode; log debug info (very verbose)')
	parser.add_option('-o', '--output-directory', default='',
					  action="store", type="string", dest="outputDirectory",
					  help='use DIR as output directory', metavar='DIR')
	parser.add_option('-s', '--scratch-directory', default='',
					  action="store", type="string", dest="scratchDirectory",
					  help='use DIR as scratch directory', metavar='DIR')
	parser.add_option('-l', '--library-directory', default='',
					  action="store", type="string", dest="libraryDirectory",
					  help='use DIR as library directory', metavar='DIR')
	parser.add_option('-r', '--restart',
						action="store_true", dest="restart", default=False,
						help="restart an incomplete job." )
	parser.add_option('-p', '--profile',
						action="store_true", dest="profile", default=False,
						help="run under cProfile to gather profiling statistics, and postprocess them if job completes." )
	parser.add_option('-P', '--postprocess',
					action="store_true", dest="postprocess_only", default=False,
					help="postprocess profiling statistics from previous [failed] run. Does not run the simulation.")
	parser.add_option('-t', '--walltime',
					action="store", type="string", dest="wallTime", default="0",
					help="set the maximum execution time (HH:MM:SS)")

	# Parse the command-line arguments
	options, args = parser.parse_args()

	# There should be exactly one positional argument: the input file
	# If this is not the case, print the help information and stop
	if len(args) != 1:
		parser.parse_args(['-h'])
		quit()

	# For output and scratch directories, if they are empty strings, set them
	# to match the input file location
	import os.path
	inputDirectory = os.path.abspath(os.path.dirname(args[0]))
	if options.outputDirectory == '':
		options.outputDirectory = inputDirectory
	if options.scratchDirectory == '':
		options.scratchDirectory = inputDirectory

	# Execute RMG
	import rmg.main

	if options.postprocess_only:
		print "Postprocessing the profiler statistics (will be appended to RMG.log)"
		options.profile = True


	if options.profile:
		import cProfile, sys, pstats, os
		global_vars = {}
		local_vars = {'args': args, 'options':options, 'rmg':rmg}
		command = """rmg.main.execute( args[0], options )"""
		stats_file = os.path.join(options.outputDirectory,'RMG.profile')
		print("Running under cProfile")
		if not options.postprocess_only:
			# actually run the program!
			cProfile.runctx(command, global_vars, local_vars, stats_file)
		# postprocess the stats
		log_file = os.path.join(options.outputDirectory,'RMG.log')
		process_stats(stats_file, log_file)

	else:
		rmg.main.execute( args[0], options )

