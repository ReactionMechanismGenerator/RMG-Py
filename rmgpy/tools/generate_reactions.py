#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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
This script is used to generate all the possible reactions involving a given
set of reactants, including pressure-dependent effects if desired. This is
effectively the first step in the RMG rate-based mechanism generation algorithm.

The input file is a subset of that used with regular RMG jobs. 
"""

import os.path
import argparse
import logging

from rmgpy.rmg.main import initializeLog, RMG
from rmgpy.chemkin import ChemkinWriter
from rmgpy.rmg.output import OutputHTMLWriter

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to RMG Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser(description=
    """
    Reaction Mechanism Generator (RMG) is an automatic chemical reaction
    mechanism generator that constructs kinetic models composed of
    elementary chemical reaction steps using a general understanding of
    how molecules react.
    """)
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
        help='a file describing the job to execute')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')

    # Add options for controlling what directories files are written to
    parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as output directory')
    parser.add_argument('-s', '--scratch-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as scratch directory')
    parser.add_argument('-l', '--library-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as library directory')
    
    args = parser.parse_args()
    args.walltime = '0'
    args.restart = False

    return args


def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.
    """
    # Parse the command-line arguments (requires the argparse module)
    args = parseCommandLineArguments()

    # For output and scratch directories, if they are empty strings, set them
    # to match the input file location
    inputFile = args.file[0]

    inputDirectory = os.path.abspath(os.path.dirname(inputFile))

    if args.output_directory == '':
        args.output_directory = inputDirectory
    if args.scratch_directory == '':
        args.scratch_directory = inputDirectory
    
    # Initialize the logging system (resets the RMG.log file)
    level = logging.INFO
    if args.debug: level = 0
    elif args.verbose: level = logging.DEBUG
    elif args.quiet: level = logging.WARNING

    kwargs = {
            'scratch_directory': args.scratch_directory,
            'restart': args.restart,
            'walltime': args.walltime,
            'log': level,
            }

    initializeLog(level, os.path.join(args.output_directory,'RMG.log'))

    rmg = RMG()

    # Add output listeners:
    rmg.attach(ChemkinWriter())
    rmg.attach(OutputHTMLWriter())

    execute(rmg, inputFile, args.output_directory, **kwargs)

def execute(rmg, inputFile, output_directory, **kwargs):
    """

    Generates all the possible reactions involving a given
    set of reactants, including pressure-dependent effects if desired. This is
    effectively the first step in the RMG rate-based mechanism generation algorithm.

    The input file is a subset of that used with regular RMG jobs. 

    Returns an RMG object.
    """   

    rmg.initialize(inputFile, output_directory, **kwargs)
    
    # Show all core and edge species and reactions in the output
    rmg.reactionModel.outputSpeciesList.extend(rmg.reactionModel.edge.species)
    rmg.reactionModel.outputReactionList.extend(rmg.reactionModel.edge.reactions)
            
    rmg.saveEverything()

    rmg.finish()
    
    return rmg