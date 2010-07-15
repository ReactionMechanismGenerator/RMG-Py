#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   MEASURE - Master Equation Automatic Solver for Unimolecular REactions
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

"""
This is the primary MEASURE module. To run MEASURE, invoke this script
via ::

$ python measure.py FILE

where ``FILE`` is the path to a valid MEASURE input file describing the job
to be run and providing the necessary information about the unimolecular
reaction network. Other command-line arguments control the level of 
verbosity of information printed to the console.
"""

import argparse
import logging
import time

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to MEASURE. These are
    described in the module docstring.
    """

    parser = argparse.ArgumentParser(description='Master Equation Automatic Solver for Unimolecular REactions.')
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
        help='a file containing information about the network')
    
    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')

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

def logHeader(level=logging.INFO):
    """
    Output a header containing identifying information about RMG to the log.
    """

    logging.log(level, '###############################################################')
    logging.log(level, '# Master Equation Automatic Solver for Unimolecular REactions #')
    logging.log(level, '# (MEASURE)                                                   #')
    logging.log(level, '# Release: 0.1.0 (7 July 2010)                                #')
    logging.log(level, '# Author: Joshua W. Allen (jwallen@mit.edu)                   #')
    logging.log(level, '# Website: http://jwallen.github.com/MEASURE                  #')
    logging.log(level, '###############################################################\n')

################################################################################

if __name__ == '__main__':
    
    # Parse the command-line arguments
    args = parseCommandLineArguments()
    
    # Initialize the logging system
    initializeLogging(args)
    
    # Log start timestamp
    logging.info('MEASURE execution initiated at ' + time.asctime() + '\n')
    
    # Log header
    logHeader()
    
    # Load input file
    from measure.input import readInput
    network, Tlist, Plist, Elist, method = readInput(args.file[0])
    
    # Only proceed if the input network is valid
    if network is not None:
        
        # Automatically choose a suitable set of energy grains if they were not
        # explicitly specified in the input file
        if len(Elist) == 2:
            logging.info('Automatically determining energy grains...')
            Tmax = max(Tlist)
            grainSize, Ngrains = Elist
            Elist = network.autoGenerateEnergyGrains(Tmax=Tmax, grainSize=grainSize, Ngrains=Ngrains)
            logging.debug('Using %i energy grains from %g to %g kJ/mol in steps of %g kJ/mol' % (len(Elist), Elist[0] / 1000, Elist[-1] / 1000, (Elist[1] - Elist[0]) / 1000))
            logging.debug('')
        
        # Calculate the rate coefficients
        K = network.calculateRateCoefficients(Tlist, Plist, Elist, method)
        
    # Log end timestamp
    logging.info('')
    logging.info('MEASURE execution terminated at ' + time.asctime())
    