#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains the :class:`CanTherm` class, the main class used to run
CanTherm.
"""

import os
import os.path
import sys
import logging
import argparse
import time
try:
    import matplotlib
    matplotlib.rc('mathtext', default='regular')
except ImportError:
    pass

from rmgpy.cantherm.input import loadInputFile

################################################################################

class CanTherm:
    """
    The :class:`CanTherm` class represents an instance of CanTherm, a tool for
    computing properties of chemical species and reactions. The attributes are:
    
    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `jobList`           A list of the jobs to execute
    `inputFile`         The path of the input file defining the jobs to execute
    `outputDirectory`   The directory in which to write the output files
    `verbose`           The level of detail in the generated logging messages
    =================== ========================================================
    
    The output directory defaults to the same directory as the input file if
    not explicitly specified.
    
    To use this class programmatically, create an instance and set its
    attributes using either the :meth:`__init__()` method or by directly
    accessing the attributes, and then invoke the :meth:`execute()` method.
    You can also populate the attributes from the command line using the
    :meth:`parseCommandLineArguments()` method before running :meth:`execute()`.
    """
    
    def __init__(self, inputFile=None, outputDirectory=None, verbose=logging.INFO):
        self.jobList = []
        self.inputFile = inputFile
        self.outputDirectory = outputDirectory
        self.verbose = verbose
    
    def parseCommandLineArguments(self):
        """
        Parse the command-line arguments being passed to CanTherm. This uses the
        :mod:`argparse` module, which ensures that the command-line arguments are
        sensible, parses them, and returns them.
        """
    
        parser = argparse.ArgumentParser(description=
        """
        CanTherm is a Python toolkit for computing chemical reaction rates and other
        properties used in detailed kinetics models using various methodologies
        and theories.
        """)
        parser.add_argument('file', metavar='FILE', type=str, nargs=1,
            help='a file describing the job to execute')
    
        # Options for controlling the amount of information printed to the console
        # By default a moderate level of information is printed; you can either
        # ask for less (quiet), more (verbose), or much more (debug)
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-q', '--quiet', action='store_const', const=logging.WARNING, default=logging.INFO, dest='verbose', help='only print warnings and errors')
        group.add_argument('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO, dest='verbose', help='print more verbose output')
        group.add_argument('-d', '--debug', action='store_const', const=0, default=logging.INFO, dest='verbose', help='print debug information')
    
        # Add options for controlling what directories files are written to
        parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
            metavar='DIR', help='use DIR as output directory')

        # Add options for controlling generation of plots
        parser.add_argument('-p', '--plot', action='store_true', default=False, help='generate plots of results')

        args = parser.parse_args()
        
        # Extract the input file
        self.inputFile = args.file[0] 
        
        # Extract the log verbosity
        self.verbose = args.verbose
        
        # Extract the plot settings
        self.plot = args.plot
        
        # Determine the output directory
        # By default the directory containing the input file is used, unless an
        # alternate directory is specified using the -o flag
        if args.output_directory and os.path.isdir(args.output_directory[0]):
            self.outputDirectory = os.path.abspath(args.output_directory[0])
        else:
            self.outputDirectory = os.path.dirname(os.path.abspath(args.file[0]))
    
    def initializeLog(self, verbose=logging.INFO, logFile=None):
        """
        Set up a logger for CanTherm to use to print output to stdout. The
        `verbose` parameter is an integer specifying the amount of log text seen
        at the console; the levels correspond to those of the :data:`logging` module.
        """
        # Create logger
        logger = logging.getLogger()
        logger.setLevel(verbose)
    
        # Use custom level names for cleaner log output
        logging.addLevelName(logging.CRITICAL, 'Critical: ')
        logging.addLevelName(logging.ERROR, 'Error: ')
        logging.addLevelName(logging.WARNING, 'Warning: ')
        logging.addLevelName(logging.INFO, '')
        logging.addLevelName(logging.DEBUG, '')
        logging.addLevelName(0, '')
    
        # Create formatter and add to handlers
        formatter = logging.Formatter('%(levelname)s%(message)s')
        
        # Remove old handlers before adding ours
        while logger.handlers:
            logger.removeHandler(logger.handlers[0])
       
        # Create console handler; send everything to stdout rather than stderr
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(verbose)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    
        # Create file handler; always be at least verbose in the file
        if logFile:
            fh = logging.FileHandler(filename=logFile)
            fh.setLevel(min(logging.DEBUG,verbose))
            fh.setFormatter(formatter)
            logger.addHandler(fh)
    
    def logHeader(self, level=logging.INFO):
        """
        Output a header containing identifying information about CanTherm to the log.
        """
        logging.log(level, 'CanTherm execution initiated at {0}'.format(time.asctime()))
        logging.log(level, '')
    
        logging.log(level, '###############################################################')
        logging.log(level, '#                                                             #')
        logging.log(level, '#                          CanTherm                           #')
        logging.log(level, '#                                                             #')
        logging.log(level, '#   Version: 0.1.0 (14 May 2009)                              #')
        logging.log(level, '#   Authors: RMG Developers (rmg_dev@mit.edu)                 #')
        logging.log(level, '#   P.I.:    William H. Green (whgreen@mit.edu)               #')
        logging.log(level, '#   Website: http://rmg.sourceforge.net/                      #')
        logging.log(level, '#                                                             #')
        logging.log(level, '###############################################################')
        logging.log(level, '')
    
    def logFooter(self, level=logging.INFO):
        """
        Output a footer to the log.
        """
        logging.log(level, '')
        logging.log(level, 'CanTherm execution terminated at {0}'.format(time.asctime()))

    def loadInputFile(self, inputFile):
        """
        Load a set of jobs from the given `inputFile` on disk. Returns the
        loaded set of jobs as a list.
        """
        self.inputFile = inputFile
        self.jobList = loadInputFile(self.inputFile)
        logging.info('')
        return self.jobList
        
    def execute(self):
        """
        Execute, in order, the jobs found in input file specified by the
        `inputFile` attribute.
        """
        
        # Initialize the logging system (both to the console and to a file in the
        # output directory)
        self.initializeLog(self.verbose, os.path.join(self.outputDirectory, 'cantherm.log'))
        
        # Print some information to the beginning of the log
        self.logHeader()
        
        # Load the input file for the job
        self.jobList = self.loadInputFile(self.inputFile)
        logging.info('')
        
        # Initialize (and clear!) the output files for the job
        if self.outputDirectory is None:
            self.outputDirectory = os.path.dirname(os.path.abspath(self.inputFile))
        outputFile = os.path.join(self.outputDirectory, 'output.py')
        with open(outputFile, 'w') as f:
            pass
        chemkinFile = os.path.join(self.outputDirectory, 'chem.inp')
        with open(chemkinFile, 'w') as f:
            pass
        
        # Run the jobs
        for job in self.jobList:
            job.execute(outputFile=outputFile, plot=self.plot)
        
        # Print some information to the end of the log
        self.logFooter()
