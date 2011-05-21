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

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to MEASURE. These are
    described in the module docstring.
    """

    parser = argparse.ArgumentParser(description="""
        Master Equation Automatic Solver for Unimolecular REactions (MEASURE):
        A tool for estimating pressure-dependent phenomenological rate
        coefficients k(T,P) for unimolecular reaction networks of arbitrary
        size and complexity using the master equation. Multiple methods of
        varying accuracy, speed, and robustness are available for determining
        the k(T,P) values. The output is a set of k(T,P) functions suitable for
        use in chemical kinetics mechanisms.
    """)
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
        help='a file containing information about the network')
    
    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument('-d', '--draw', metavar='IMGFILE', type=str, nargs=1,
        help='draw potential energy surface and exit')
    group1.add_argument('-o', '--output', metavar='OUTFILE', type=str, nargs=1,
        help='specify location of output file')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group2.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')

    return parser.parse_args()

################################################################################

if __name__ == '__main__':
    
    # Parse the command-line arguments
    args = parseCommandLineArguments()
    
    quiet = args.quiet
    verbose = args.verbose
    
    inputFile = args.file[0]
    outputFile = args.output[0] if args.output is not None else None
    drawFile = args.draw[0] if args.draw is not None else None
    
    # Run MEASURE
    from rmgpy.measure.main import execute
    execute(inputFile, outputFile=outputFile, drawFile=drawFile, logFile=None, quiet=quiet, verbose=verbose)
