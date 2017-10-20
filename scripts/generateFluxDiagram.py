#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script generates a video showing the flux diagram for a given reaction
model as it evolves in time. It takes as its arguments the path to an RMG-Py
input file corresponding to a job that has already been run and the
corresponding Chemkin mechanism and RMG species dictionary files. If a folder
of species images is available, it can be passed as an optional argument.
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
    parser.add_argument('--java', action='store_true', help='process RMG-Java model')
    parser.add_argument('--no-dlim', dest='dlim', action='store_false', help='Turn off diffusion-limited rates')

    args = parser.parse_args()

    inputFile = os.path.abspath(args.input)
    chemkinFile = os.path.abspath(args.chemkin)
    dictFile = os.path.abspath(args.dictionary)
    speciesPath = os.path.abspath(args.species) if args.species is not None else None
    useJava = args.java
    dflag = args.dlim
    
    return inputFile, chemkinFile, dictFile, speciesPath, useJava, dflag

def main():
    inputFile, chemkinFile, dictFile, speciesPath, useJava, dflag = parse_arguments()

    createFluxDiagram(inputFile, chemkinFile, dictFile, speciesPath=speciesPath, java=useJava, diffusionLimited=dflag)

if __name__ == '__main__':
    main()