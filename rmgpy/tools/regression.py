#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains classes and functions for comparing observables between
two RMG generated models.
"""
import argparse
import logging
import os.path
import sys

import rmgpy
rmgpy.DISABLE_JULIA = True  # disable Julia to save time and warnings

from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.species import Species
from rmgpy.rmg.input import fragment_adj, fragment_smiles, smiles, adjacency_list
from rmgpy.tools.canteramodel import CanteraCondition
from rmgpy.tools.observablesregression import ObservablesTestCase

observables = []
setups = None
casetitle = ''
initialSpecies = {}
tol = 0.05


def read_input_file(path):
    """
    Read an regression input file at `path` on disk.
    """
    global observables, setups, initialSpecies, casetitle, tol

    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise

    logging.info('Reading input file "{0}"...'.format(full_path))
    logging.info(f.read())
    f.seek(0)  # return to beginning of file

    global_context = {'__builtins__': None}
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'observable': observable,
        'SMILES': smiles,
        'fragment_adj': fragment_adj,
        'fragment_SMILES': fragment_smiles,
        'species': species,
        'adjacencyList': adjacency_list,
        'reactorSetups': reactorSetups,
        'options': options,
    }

    try:
        exec(f.read(), global_context, local_context)
    except (NameError, TypeError, SyntaxError) as e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    # convert keys of the initial mole fraction dictionaries from labels into Species objects.
    imfs = setups[3]
    new_imfs = []
    for imf in imfs:
        new_dict = convert(imf, initialSpecies)
        new_imfs.append(new_dict)
    setups[3] = new_imfs

    # convert keys of the initial surface mole fraction dictionaries (ismfs) from labels into Species objects.
    if setups[5][0]:
        ismfs = setups[5]
        new_ismfs = []
        for ismf in ismfs:
            new_dict = convert(ismf, initialSpecies)
            new_ismfs.append(new_dict)
        setups[5] = new_ismfs

    return casetitle, observables, setups, tol


def observable(label, structure):
    spc = species(label, structure)
    observables.append(spc)


def species(label, structure):
    spc = Species(label=label, molecule=[structure])
    initialSpecies[label] = spc
    return spc


def reactorSetups(reactorTypes, temperatures, pressures, initialMoleFractionsList, terminationTimes, initialSurfaceMoleFractionsList=None):
    global setups

    if initialSurfaceMoleFractionsList is None:
        initialSurfaceMoleFractionsList = [None]

    terminationTimes = Quantity(terminationTimes)
    temperatures = Quantity(temperatures)
    pressures = Quantity(pressures)

    setups = [reactorTypes, temperatures, pressures, initialMoleFractionsList, terminationTimes, initialSurfaceMoleFractionsList]


def options(title='', tolerance=0.05):
    global casetitle, tol
    casetitle = title
    tol = tolerance


def convert(origDict, initialSpecies):
    """
    Convert the original dictionary with species labels as keys
    into a new dictionary with species objects as keys,
    using the given dictionary of species.
    """
    new_dict = {}

    for label, value in origDict.items():
        new_dict[initialSpecies[label]] = value

    return new_dict


def run(benchmarkDir, testDir, title, observables, setups, tol):
    case = ObservablesTestCase(
        title=title,
        old_dir=benchmarkDir,
        new_dir=testDir,
        observables={'species': observables},
        ck2cti=False,
    )

    reactor_types, temperatures, pressures, initial_mole_fractions_list, termination_times, initial_surface_mole_fractions_list = setups
    case.generate_conditions(
        reactor_type_list=reactor_types,
        reaction_time_list=termination_times,
        mol_frac_list=initial_mole_fractions_list,
        surface_mol_frac_list=initial_surface_mole_fractions_list,
        Tlist=temperatures,
        Plist=pressures
    )

    variables_failed = case.compare(tol)
    return variables_failed  # will be None if no failures


def parse_command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
                        help='regression input file')
    parser.add_argument('benchmark', metavar='BENCHMARK', type=str, nargs=1,
                        help='folder of the benchmark model')
    parser.add_argument('tested', metavar='TESTED', type=str, nargs=1,
                        help='folder of the tested model')

    args = parser.parse_args()

    input_file = os.path.abspath(args.input[0])
    benchmark = os.path.abspath(args.benchmark[0])
    tested = os.path.abspath(args.tested[0])

    return input_file, benchmark, tested


def main():
    "Returns the list of variables that failed the regression."
    input_file, benchmark, tested = parse_command_line_arguments()

    args = read_input_file(input_file)  # casetitle, observables, setups, tol

    return run(benchmark, tested, *args)


if __name__ == '__main__':
    variables_failed = main()
    sys.exit(1 if variables_failed else 0)
