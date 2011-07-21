#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert a FAME input file to a MEASURE input file.
"""

import argparse
import numpy
import os.path

from rmgpy.molecule import Molecule

from rmgpy.measure.main import MEASURE

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to MEASURE. These are
    described in the module docstring.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('file', metavar='FILE', type=str, nargs='+',
        help='a file to convert')
    parser.add_argument('-d', '--dictionary', metavar='DICTFILE', type=str, nargs=1,
        help='the RMG dictionary corresponding to these files')

    return parser.parse_args()

################################################################################

if __name__ == '__main__':
    
    # Parse the command-line arguments
    args = parseCommandLineArguments()
    
    # Load RMG dictionary if specified
    moleculeDict = {}
    if args.dictionary is not None:
        f = open(args.dictionary[0])
        adjlist = ''; label = ''
        for line in f:
            if len(line.strip()) == 0:
                if len(adjlist.strip()) > 0:
                    molecule = Molecule()
                    molecule.fromAdjacencyList(adjlist)
                    moleculeDict[label] = molecule
                adjlist = ''; label = ''
            else:
                if len(adjlist.strip()) == 0:
                    label = line.strip()
                adjlist += line
                    
        f.close()
    
    method = None

    for fstr in args.file:

        # Construct MEASURE job from FAME input
        measure = MEASURE()
        measure.loadFAMEInput(fstr, moleculeDict)

        # Save MEASURE input file based on the above
        dirname, basename = os.path.split(os.path.abspath(fstr))
        basename, ext = os.path.splitext(basename)
        path = os.path.join(dirname, basename + '.py')
        measure.saveInput(path)
