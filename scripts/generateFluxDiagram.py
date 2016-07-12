#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script generates a video showing the flux diagram for a given reaction
model as it evolves in time. It takes as its lone required argument the path
to an RMG-Py input file corresponding to a job that has already been run.
This script will automatically read from the necessary output files to extract
the information needed to generate the flux diagram.
"""

import os.path
import argparse

import rmgpy.tools.fluxdiagram as fluxdiagram

################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='the RMG input file to use')
    parser.add_argument('species', metavar='SPECIES', type=str, nargs='?', default=None,
        help='path to species images')
    parser.add_argument('--java', action='store_true', help='process RMG-Java model')

    args = parser.parse_args()

    inputFile = os.path.abspath(args.input[0])
    speciesPath = os.path.abspath(args.species[0]) if args.species is not None else None
    useJava = args.java
    
    return inputFile, speciesPath, useJava

def main():
    # This might not work anymore because functions were modified for use with webserver

    inputFile, speciesPath, useJava = parse_arguments()

    fluxdiagram.run(inputFile, speciesPath, useJava)

if __name__ == '__main__':
    main()