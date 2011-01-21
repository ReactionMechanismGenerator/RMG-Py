#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   CanTherm
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
    
    # Option for controlling the output file location
    parser.add_argument('-o', '--output', metavar='OUTFILE', type=str, nargs=1,
        help='specify location of output file')

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

def execute(path, output):
    """
    Execute the CanTherm job located at `path` on disk, saving the output to
    `output` on disk.
    """
    
    try:
        f = open(path)
    except IOError, e:
        logging.error('The input file "%s" could not be opened.' % (path))
        logging.error('Check that the file exists and that you have read access.')
        return

    logging.info('Executing job "%s"...' % os.path.abspath(path))

    # Change directory to that of the job file
    # This way all relative paths given in the job file will be valid
    root, ext = os.path.split(path)
    os.chdir(root)

    from cantherm.input import setOutputFile, setModelChemistry, hinderedRotor, loadSpecies, loadTransitionState, loadReaction, generateStates, generateThermo, generateKinetics
    setOutputFile(output)    
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'range': range,
        'modelChemistry': setModelChemistry,
        'HinderedRotor': hinderedRotor,
        'species': loadSpecies,
        'transitionState': loadTransitionState,
        'reaction': loadReaction,
        'states': generateStates,
        'thermo': generateThermo,
        'kinetics': generateKinetics,
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "%s" was invalid:' % (path))
        raise
    finally:
        f.close()

################################################################################

if __name__ == '__main__':
    
    # Parse the command-line arguments
    args = parseCommandLineArguments()
    outputDirectory = os.path.dirname(os.path.abspath(args.file[0]))
    out = os.path.join(outputDirectory, 'output.py')
    if args.output:
        out = os.path.abspath(args.output[0])
            
    # Initialize the logging system
    initializeLogging(args)
    
    # Log start timestamp
    logging.info('CanTherm execution initiated at ' + time.asctime() + '\n')
    
    # Execute job
    execute(args.file[0], out)

    # Log end timestamp
    logging.info('')
    logging.info('CanTherm execution terminated at ' + time.asctime())
    