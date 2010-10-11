#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   CanTherm - 
#    
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import argparse
import logging
import time
import os.path
import numpy

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to CanTherm.
    """

    parser = argparse.ArgumentParser(description="""
        CanTherm is a tool for computing the thermodynamics properties of
        chemical species and the high-pressure-limit rate coefficients for
        chemical reactions using the results of quantum chemistry calculations.
    """)
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
        help='a job to run')
    
    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group2.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')

    return parser.parse_args()

################################################################################

def initializeLogging(args):
    """
    Initialize the logging system. The level of information printed is 
    determined by looking at the ``args.quiet`` and ``args.verbose`` attributes
    to see if either of the corresponding flags were set. The parameter `args`
    is an object returned by the ``argparse`` module.
    """

    # Set up logging system
    level = logging.INFO
    if args.quiet: 
        level = logging.WARNING
    elif args.verbose: 
        level = logging.DEBUG

    # Reassign the level names so that they look better on printing
    logging.addLevelName(logging.CRITICAL, 'CRITICAL: ')
    logging.addLevelName(logging.ERROR, 'ERROR: ')
    logging.addLevelName(logging.WARNING, 'Warning: ')
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, '')

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Create console handler and set level to debug
    # Also send everything to stdout rather than stderr
    import sys
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)

    # Create formatter and add to console handler
    formatter = logging.Formatter('%(levelname)s%(message)s')
    ch.setFormatter(formatter)

    # Remove any old handlers that might exist
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Add ch to logger
    logger.addHandler(ch)

################################################################################

def execute(path):
    """
    Execute the CanTherm job located at `path` on disk.
    """
    logging.info('Executing job "%s"...' % os.path.abspath(path))

################################################################################

if __name__ == '__main__':
    
    # Parse the command-line arguments
    args = parseCommandLineArguments()
    
    # Initialize the logging system
    initializeLogging(args)
    
    # Log start timestamp
    logging.info('CanTherm execution initiated at ' + time.asctime() + '\n')
    
    # Execute job
    execute(args.file[0])

    # Log end timestamp
    logging.info('')
    logging.info('CanTherm execution terminated at ' + time.asctime())
    