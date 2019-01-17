#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script runs stand-alone sensitivity analysis on an RMG job.
"""

import os.path
import argparse

from rmgpy.tools.sensitivity import runSensitivity

################################################################################

def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='RMG input file')
    parser.add_argument('chemkin', metavar='CHEMKIN', type=str, nargs=1,
        help='Chemkin file')
    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, nargs=1,
        help='RMG dictionary file')
    parser.add_argument('--no-dlim', dest='dlim', action='store_false',
        help='Turn off diffusion-limited rates for LiquidReactor')
    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    chemkinFile = os.path.abspath(args.chemkin[0])
    dictFile = os.path.abspath(args.dictionary[0])
    dflag = args.dlim

    return inputFile, chemkinFile, dictFile, dflag

def main():
    # This might not work anymore because functions were modified for use with webserver

    inputFile, chemkinFile, dictFile, dflag = parse_arguments()

    runSensitivity(inputFile, chemkinFile, dictFile, dflag)

################################################################################

if __name__ == '__main__':
    main()

    
    