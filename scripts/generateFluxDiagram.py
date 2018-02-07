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
This script generates a video showing the flux diagram for a given reaction
model as it evolves in time. It takes as its arguments the path to an RMG-Py
input file corresponding to a job that has already been run and the
corresponding Chemkin mechanism and RMG species dictionary files. If a folder
of species images is available, it can be passed as an optional argument. A
Chemkin output file can also be passed as an optional positional argument.
"""

import os
import argparse

from rmgpy.tools.fluxdiagram import createFluxDiagram

################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('input', metavar='INPUT', type=str, help='RMG input file')
    parser.add_argument('chemkin', metavar='CHEMKIN', type=str, help='Chemkin file')
    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, help='RMG dictionary file')
    parser.add_argument('chemkinOutput', metavar='CHEMKIN_OUTPUT', type=str, nargs='?', default=None,
                        help='Chemkin output file')
    parser.add_argument('--java', action='store_true', help='process RMG-Java model')
    parser.add_argument('--no-dlim', dest='dlim', action='store_false', help='Turn off diffusion-limited rates')
    parser.add_argument('-s', '--species', metavar='DIR', type=str, help='Path to folder containing species images')
    parser.add_argument('-f', '--foreign', dest='checkDuplicates', action='store_true',
                        help='Not an RMG generated Chemkin file (will be checked for duplicates)')
    parser.add_argument('-n', '--maxnode', metavar='N', type=int, help='Maximum number of nodes to show in diagram')
    parser.add_argument('-e', '--maxedge', metavar='N', type=int, help='Maximum number of edges to show in diagram')
    parser.add_argument('-c', '--conctol', metavar='TOL', type=float, help='Lowest fractional concentration to show')
    parser.add_argument('-r', '--ratetol', metavar='TOL', type=float, help='Lowest fractional species rate to show')
    parser.add_argument('-t', '--tstep', metavar='S', type=float,
                        help='Multiplicative factor to use between consecutive time points')
    parser.add_argument('--centralSpecies', metavar='s1,s2,...', type=lambda s: [int(idx) for idx in s.split(',')],
                        help='List of indices of central species')
    parser.add_argument('--rad', metavar='R', type=int, help='Graph radius around a central species')
    parser.add_argument('--super', action='store_true', help='Superimpose central species onto normal flux diagram to'
                                                             ' ensure that they appear in diagram')

    args = parser.parse_args()

    inputFile = os.path.abspath(args.input)
    chemkinFile = os.path.abspath(args.chemkin)
    dictFile = os.path.abspath(args.dictionary)
    speciesPath = os.path.abspath(args.species) if args.species is not None else None
    chemkinOutput = os.path.abspath(args.chemkinOutput) if args.chemkinOutput is not None else ''
    useJava = args.java
    dflag = args.dlim
    checkDuplicates = args.checkDuplicates
    centralSpeciesList = args.centralSpecies
    superimpose = args.super

    keys = ('maximumNodeCount', 'maximumEdgeCount', 'concentrationTolerance', 'speciesRateTolerance', 'radius', 'timeStep')
    vals = (args.maxnode, args.maxedge, args.conctol, args.ratetol, args.rad, args.tstep)
    settings = {k: v for k, v in zip(keys, vals) if v is not None}
    
    return inputFile, chemkinFile, dictFile, speciesPath, chemkinOutput, useJava, dflag, checkDuplicates, settings, centralSpeciesList, superimpose

def main():
    inputFile, chemkinFile, dictFile, speciesPath, chemkinOutput, useJava, dflag, checkDuplicates, settings, centralSpeciesList, superimpose = parse_arguments()

    createFluxDiagram(inputFile, chemkinFile, dictFile, speciesPath=speciesPath, java=useJava, settings=settings,
                      chemkinOutput=chemkinOutput, diffusionLimited=dflag, centralSpeciesList=centralSpeciesList,
                      superimpose=superimpose, checkDuplicates=checkDuplicates)

if __name__ == '__main__':
    main()