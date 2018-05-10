#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

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

