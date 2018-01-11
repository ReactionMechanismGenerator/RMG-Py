#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script generates an html page of species and reactions from a chemkin
input file and RMG species dictionary.

To use, pass the paths of the Chemkin file and species dictionary on the
command-line, e.g.

    $ python generateChemkinHTML.py /path/to/chem.inp /path/to/species_dictionary.txt [/path/to/output/directory/]

The resulting HTML file and species image folder are placed in the execution
directory, unless an output directory is specified.
"""

import os
import argparse

from rmgpy.chemkin import loadChemkinFile
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.output import saveOutputHTML

################################################################################

def main(chemkin, dictionary, output, foreign):
    model = CoreEdgeReactionModel()
    model.core.species, model.core.reactions = loadChemkinFile(chemkin, dictionary, readComments=not foreign)
    outputPath = os.path.join(output, 'output.html')
    speciesPath = os.path.join(output, 'species')
    if not os.path.isdir(speciesPath):
        os.makedirs(speciesPath)
    saveOutputHTML(outputPath, model)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('chemkin', metavar='CHEMKIN', type=str, nargs=1,
                        help='the Chemkin input file to visualize')
    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, nargs=1,
                        help='the RMG species dictionary')
    parser.add_argument('output', metavar='OUTPUT', type=str, nargs='?', default=None,
                        help='output directory, defaults to current directory')
    parser.add_argument('-f', '--foreign', action='store_true', help='not an RMG generated Chemkin file')

    args = parser.parse_args()

    chemkin = os.path.abspath(args.chemkin[0])
    dictionary = os.path.abspath(args.dictionary[0])
    output = os.path.abspath(args.output[0]) if args.output is not None else os.getcwd()
    foreign = args.foreign

    main(chemkin, dictionary, output, foreign)

