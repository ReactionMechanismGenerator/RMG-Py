#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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
This module contains classes and functions for comparing observables between
two RMG generated models.
"""
import logging
import os.path
import argparse

from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.species import Species
from rmgpy.tools.observablesRegression import ObservablesTestCase
from rmgpy.tools.canteraModel import CanteraCondition

observables = []
setups = None
casetitle = ''
initialSpecies = {}
tol = 0.05

def readInputFile(path):
    """
    Read an regression input file at `path` on disk.
    """
    global observables, setups, initialSpecies, casetitle, tol

    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise e

    logging.info('Reading input file "{0}"...'.format(full_path))
    logging.info(f.read())
    f.seek(0)# return to beginning of file

    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'observable': observable,
        'SMILES': SMILES,
        'species': species,
        'adjacencyList': adjacencyList,
        'reactorSetups': reactorSetups,
        'options': options,
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    # convert keys of the initial mole fraction dictionaries from labels into Species objects.
    imfs = setups[3]
    newImfs = []
    for imf in imfs:
        newDict = convert(imf, initialSpecies)
        newImfs.append(newDict)
    setups[3] = newImfs

    return casetitle, observables, setups, tol

def observable(label, structure):
    spc = species(label, structure)
    observables.append(spc)

def species(label, structure):
    spc = Species(label=label, molecule=[structure])
    initialSpecies[label] = spc
    return spc

def reactorSetups(reactorTypes, temperatures, pressures, initialMoleFractionsList, terminationTimes):

    global setups

    terminationTimes = Quantity(terminationTimes)
    temperatures = Quantity(temperatures)
    pressures = Quantity(pressures)

    setups = [reactorTypes, temperatures, pressures, initialMoleFractionsList, terminationTimes]

def SMILES(string):
    return Molecule().fromSMILES(string)

def adjacencyList(string):
    return Molecule().fromAdjacencyList(string)

def options(title = '', tolerance=0.05):

    global casetitle, tol
    casetitle = title
    tol = tolerance

def convert(origDict, initialSpecies):
    """
    Convert the original dictionary with species labels as keys
    into a new dictionary with species objects as keys,
    using the given dictionary of species.
    """
    newDict = {}
    
    for label, value in origDict.iteritems():
        newDict[initialSpecies[label]] = value

    return newDict

def run(benchmarkDir, testDir, title, observables, setups, tol):
    
    case = ObservablesTestCase(
        title = title,
        oldDir = benchmarkDir,
        newDir = testDir,
        observables = {'species': observables}
        )

    reactorTypes, temperatures, pressures, initialMoleFractionsList, terminationTimes = setups
    case.generateConditions(
        reactorTypeList = reactorTypes,
        reactionTimeList = terminationTimes,
        molFracList = initialMoleFractionsList,
        Tlist = temperatures,
        Plist = pressures
        )

    case.compare(tol)

def parseArguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='regression input file')
    parser.add_argument('benchmark', metavar='BENCHMARK', type=str, nargs=1,
        help='folder of the benchmark model')
    parser.add_argument('tested', metavar='TESTED', type=str, nargs=1,
        help='folder of the tested model')

    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    benchmark = os.path.abspath(args.benchmark[0])
    tested = os.path.abspath(args.tested[0])

    return inputFile, benchmark, tested


def main():
    inputFile, benchmark, tested = parseArguments()

    args = readInputFile(inputFile)

    run(benchmark, tested, *args)

if __name__ == '__main__':
    main()
    
