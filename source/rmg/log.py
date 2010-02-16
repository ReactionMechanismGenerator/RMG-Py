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
Contains functions useful for manipulating the logging information that RMG
saves during execution. In RMG we define several levels of information:

* CRITICAL - For log messages pertaining to critical (fatal) errors; rarely used
  in RMG, as we'd usually rather raise an exception

* ERROR - For log messages pertaining to normal (non-fatal) errors; rarely used
  in RMG, as we'd usually rather raise an exception

* WARNING - For log messages pertaining to warnings, which refer to unexpected
  but perhaps not incorrect conditions

* INFO - For log messages that summarize the steps taken by the code

* VERBOSE - For log messages that give more details on the steps taken by the
  code

* DEBUG - For log messages that give lots and lots of details on the steps taken
  by the code

Each level above has a corresponding function to add a message to the log at
that level.
"""

import os
import os.path
import sys

from logging import critical, error, warning, info, debug, log, \
	CRITICAL, ERROR, WARNING, INFO, DEBUG, \
	addLevelName, getLogger, StreamHandler, Formatter, FileHandler

################################################################################

VERBOSE = 15

addLevelName(CRITICAL, 'Critical: ')
addLevelName(ERROR, 'Error: ')
addLevelName(WARNING, 'Warning: ')
addLevelName(INFO, '')
addLevelName(VERBOSE, '')
addLevelName(DEBUG, '')

#def verbose(msg, *args, **kwargs):
#	log(15, msg, args, kwargs)

def verbose(msg):
	log(15, msg)

################################################################################

def initialize(verbose, fstr):
	"""
	Set up a logger for RMG to use to print output to stdout. The
	`verbose` parameter is an integer specifying the amount of log text seen
	at the console; the levels correspond to those of the :data:`logging` module.
	"""
	# Create logger
	logger = getLogger()
	logger.setLevel(verbose)

	# Create console handler and set level to debug; send everything to stdout
	# rather than stderr
	ch = StreamHandler(sys.stdout)
	ch.setLevel(verbose)

	# Create formatter and add to console handler
	#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
	#formatter = Formatter('%(message)s', '%Y-%m-%d %H:%M:%S')
	formatter = Formatter('%(levelname)s%(message)s')
	ch.setFormatter(formatter)

	# create file handler
	log_file_name = fstr
	if os.path.exists(log_file_name):
		backup_name = '_backup'.join(os.path.splitext(log_file_name))
		if os.path.exists(backup_name):
			print "Removing old file %s"%backup_name
			os.remove(backup_name)
		print "Renaming %s to %s"%(log_file_name, backup_name)
		os.rename(log_file_name, backup_name)
	fh = FileHandler(filename=log_file_name) #, backupCount=3)
	fh.setLevel(VERBOSE) # always verbose in the file
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

def logHeader(level=INFO):
	"""
	Output a header containing identifying information about RMG to the log.
	"""

	log(level, '################################################################')
	log(level, '#                                                              #')
	log(level, '#              RMG - Reaction Mechanism Generator              #')
	log(level, '#                    Python Version 0.0.1                      #')
	log(level, '#                         14 May 2009                          #')
	log(level, '#                                                              #')
	log(level, '#                 http://rmg.sourceforge.net/                  #')
	log(level, '#                                                              #')
	log(level, '################################################################\n')
