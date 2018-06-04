#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    parser.add_argument('species', metavar='SPECIES', type=str, nargs='?', default=None, help='Path to species images')
    parser.add_argument('chemkinOutput', metavar='CHEMKIN_OUTPUT', type=str, nargs='?', default=None,
                        help='Chemkin output file')
    parser.add_argument('--java', action='store_true', help='process RMG-Java model')
    parser.add_argument('--no-dlim', dest='dlim', action='store_false', help='Turn off diffusion-limited rates')
    parser.add_argument('-n', '--maxnode', metavar='N', type=int, help='Maximum number of nodes to show in diagram')
    parser.add_argument('-e', '--maxedge', metavar='N', type=int, help='Maximum number of edges to show in diagram')
    parser.add_argument('-c', '--conctol', metavar='TOL', type=float, help='Lowest fractional concentration to show')
    parser.add_argument('-r', '--ratetol', metavar='TOL', type=float, help='Lowest fractional species rate to show')
    parser.add_argument('-t', '--tstep', metavar='S', type=float,
                        help='Multiplicative factor to use between consecutive time points')

    args = parser.parse_args()

    inputFile = os.path.abspath(args.input)
    chemkinFile = os.path.abspath(args.chemkin)
    dictFile = os.path.abspath(args.dictionary)
    speciesPath = os.path.abspath(args.species) if args.species is not None else None
    chemkinOutput = os.path.abspath(args.chemkinOutput) if args.chemkinOutput is not None else ''
    useJava = args.java
    dflag = args.dlim

    keys = ('maximumNodeCount', 'maximumEdgeCount', 'concentrationTolerance', 'speciesRateTolerance', 'timeStep')
    vals = (args.maxnode, args.maxedge, args.conctol, args.ratetol, args.tstep)
    settings = {k: v for k, v in zip(keys, vals) if v is not None}
    
    return inputFile, chemkinFile, dictFile, speciesPath, chemkinOutput, useJava, dflag, settings

def main():
    inputFile, chemkinFile, dictFile, speciesPath, chemkinOutput, useJava, dflag, settings = parse_arguments()

    createFluxDiagram(inputFile, chemkinFile, dictFile, speciesPath=speciesPath, java=useJava, settings=settings,
                      chemkinOutput=chemkinOutput, diffusionLimited=dflag)

if __name__ == '__main__':
    main()