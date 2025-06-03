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

import argparse
import logging
import math
import sys

# Disable Julia imports for checkModels.py since they're not needed
import rmgpy
rmgpy.DISABLE_JULIA = True

from rmgpy.tools.diffmodels import execute

logger = logging.getLogger('checkModels')


def parse_command_line_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('name', metavar='NAME', type=str, nargs=1,
                        help='Name of test target model')

    parser.add_argument('benchChemkin', metavar='BENCHCHEMKIN', type=str, nargs=1,
                        help='The path to the the Chemkin file of the benchmark model')
    parser.add_argument('benchSpeciesDict', metavar='BENCHSPECIESDICT', type=str, nargs=1,
                        help='The path to the the species dictionary file of the benchmark model')

    parser.add_argument('testChemkin', metavar='TESTEDCHEMKIN', type=str, nargs=1,
                        help='The path to the the Chemkin file of the tested model')
    parser.add_argument('testSpeciesDict', metavar='TESTEDSPECIESDICT', type=str, nargs=1,
                        help='The path to the the species dictionary file of the tested model')

    args = parser.parse_args()

    return args


def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.

    Returns `True` if there is any error (discrepancy between the two models), `False` otherwise.
    """
    args = parse_command_line_arguments()

    name = args.name[0]
    initialize_log(logging.WARNING, name + '.log')

    bench_chemkin = args.benchChemkin[0]
    bench_species_dict = args.benchSpeciesDict[0]

    test_chemkin = args.testChemkin[0]
    test_species_dict = args.testSpeciesDict[0]

    error = check(name, bench_chemkin, bench_species_dict, test_chemkin, test_species_dict)
    return error


def check(name, benchChemkin, benchSpeciesDict, testChemkin, testSpeciesDict):
    """
    Compare the tested model to the benchmark model.
    """
    kwargs = {
        'web': True,
    }

    test_thermo, bench_thermo = None, None
    common_species, unique_species_orig, unique_species_test, common_reactions, unique_reactions_orig, unique_reactions_test = \
        execute(benchChemkin, benchSpeciesDict, bench_thermo, testChemkin, testSpeciesDict, test_thermo, **kwargs)

    error_model = checkModel(common_species, unique_species_test, unique_species_orig, common_reactions, unique_reactions_test,
                            unique_reactions_orig)

    error_species = checkSpecies(common_species, unique_species_test, unique_species_orig)

    error_reactions = checkReactions(common_reactions, unique_reactions_test, unique_reactions_orig)

    return error_model or error_species or error_reactions


def checkModel(commonSpecies, uniqueSpeciesTest, uniqueSpeciesOrig, commonReactions, uniqueReactionsTest,
               uniqueReactionsOrig):
    """
    Compare the species and reaction count of both models.
    """

    test_model_species = len(commonSpecies) + len(uniqueSpeciesTest)
    orig_model_species = len(commonSpecies) + len(uniqueSpeciesOrig)

    logger.error(f"Original model has {orig_model_species} species.")
    logger.error(f"Test model has {test_model_species} species. " +
                 ('✅' if test_model_species == orig_model_species else '❌'))


    test_model_rxns = len(commonReactions) + len(uniqueReactionsTest)
    orig_model_rxns = len(commonReactions) + len(uniqueReactionsOrig)

    logger.error(f"Original model has {orig_model_rxns} reactions.")
    logger.error(f"Test model has {test_model_rxns} reactions. " +
                 ('✅' if test_model_rxns == orig_model_rxns else '❌'))


    return (test_model_species != orig_model_species) or (test_model_rxns != orig_model_rxns)


def checkSpecies(commonSpecies, uniqueSpeciesTest, uniqueSpeciesOrig):
    error = False

    # check for unique species in one of the models:
    if uniqueSpeciesOrig:
        error = True
        logger.error(f'The original model has {len(uniqueSpeciesOrig)} species '
                     'that the tested model does not have. ❌')
        [printSpecies(spc) for spc in uniqueSpeciesOrig]

    if uniqueSpeciesTest:
        error = True
        logger.error(f'The tested model has {len(uniqueSpeciesTest)} species '
                     'that the original model does not have. ❌')
        [printSpecies(spc) for spc in uniqueSpeciesTest]

    # check for different thermo among common species::
    if commonSpecies:
        for spec1, spec2 in commonSpecies:
            logger.info('    %s', spec1)
            if spec1.thermo and spec2.thermo:
                if not spec1.thermo.is_similar_to(spec2.thermo):
                    error = True
                    logger.error('')
                    logger.error('Non-identical thermo! ❌')
                    logger.error(f'original: `{spec1.label}`')
                    logger.error(f'tested:   `{spec2.label}`')
                    logger.error("|{0:10}|{1:10}|{2:10}|{3:10}|{4:10}|{5:10}|{6:10}|{7:10}|{8:10}|"
                                 .format('Hf(300K)', 'S(300K)', 'Cp(300K)', 'Cp(400K)', 'Cp(500K)', 'Cp(600K)',
                                         'Cp(800K)', 'Cp(1000K)', 'Cp(1500K)')
                                 )
                    logger.error("|----------|----------|----------|----------|----------|----------|----------|----------|----------|")
                    [printThermo(spc) for spc in [spec1, spec2]]
                    logger.error('')
                    if spec1.thermo.comment != spec2.thermo.comment:
                        [printSpeciesComments(spc) for spc in [spec1, spec2]]
                    else:
                        logger.error('Identical thermo comments:')
                        printSpeciesComments(spec1)

    return error


def checkReactions(commonReactions, uniqueReactionsTest, uniqueReactionsOrig):
    error = False

    # check for unique reactions in one of the models:
    if uniqueReactionsOrig:
        error = True

        logger.error(f'The original model has {len(uniqueReactionsOrig)} reactions '
                        'that the tested model does not have. ❌')
        [printReaction(rxn) for rxn in uniqueReactionsOrig]

    if uniqueReactionsTest:
        error = True
        logger.error(f'The tested model has {len(uniqueReactionsTest)} reactions '
                        'that the original model does not have. ❌')
        [printReaction(rxn) for rxn in uniqueReactionsTest]

    if commonReactions:
        for rxn1, rxn2 in commonReactions:
            logger.info('    %s', rxn1)
            if rxn1.kinetics and rxn2.kinetics:
                if not rxn1.kinetics.is_similar_to(rxn2.kinetics):
                    error = True
                    logger.error('')
                    logger.error('Non-identical kinetics! ❌')
                    logger.error('original:')
                    printReaction(rxn1)
                    logger.error('tested:')
                    printReaction(rxn2)

                    logger.error("|{0:7}|{1:7}|{2:7}|{3:7}|{4:7}|{5:7}|{6:7}|{7:7}|{8:7}|"
                                 .format('k(1bar)', '300K', '400K', '500K', '600K', '800K', '1000K', '1500K', '2000K')
                                 )
                    logger.error("|-------|-------|-------|-------|-------|-------|-------|-------|-------|")
                    [printRates(rxn) for rxn in [rxn1, rxn2]]

                    logger.error('')
                    [printKinetics(rxn) for rxn in [rxn1, rxn2]]

                    if rxn1.kinetics.comment != rxn2.kinetics.comment:
                        [printReactionComments(rxn) for rxn in [rxn1, rxn2]]
                    else:
                        logger.error('Identical kinetics comments:')
                        printReactionComments(rxn1)

    return error


def printSpecies(spc):
    """

    """
    logger.error(
        'spc: {}'.format(spc)
    )


def printRates(rxn):
    """
    Print the rate coefficients of a reaction at various temperatures, in a markdown table.
    """
    logger.error("|{0:7}|{1:7.2f}|{2:7.2f}|{3:7.2f}|{4:7.2f}|{5:7.2f}|{6:7.2f}|{7:7.2f}|{8:7.2f}|"
        .format(
        'k(T): ',
        math.log10(rxn.kinetics.get_rate_coefficient(300, 1e5)),
        math.log10(rxn.kinetics.get_rate_coefficient(400, 1e5)),
        math.log10(rxn.kinetics.get_rate_coefficient(500, 1e5)),
        math.log10(rxn.kinetics.get_rate_coefficient(600, 1e5)),
        math.log10(rxn.kinetics.get_rate_coefficient(800, 1e5)),
        math.log10(rxn.kinetics.get_rate_coefficient(1000, 1e5)),
        math.log10(rxn.kinetics.get_rate_coefficient(1500, 1e5)),
        math.log10(rxn.kinetics.get_rate_coefficient(2000, 1e5)),
    ))


def printThermo(spec):
    """
    Print the thermo of a species at various temperatures, in a markdown table.
    """
    logger.error("|{0:10.2f}|{1:10.2f}|{2:10.2f}|{3:10.2f}|{4:10.2f}|{5:10.2f}|{6:10.2f}|{7:10.2f}|{8:10.2f}|"
        .format(
        spec.thermo.get_enthalpy(300) / 4184.,
        spec.thermo.get_entropy(300) / 4.184,
        spec.thermo.get_heat_capacity(300) / 4.184,
        spec.thermo.get_heat_capacity(400) / 4.184,
        spec.thermo.get_heat_capacity(500) / 4.184,
        spec.thermo.get_heat_capacity(600) / 4.184,
        spec.thermo.get_heat_capacity(800) / 4.184,
        spec.thermo.get_heat_capacity(1000) / 4.184,
        spec.thermo.get_heat_capacity(1500) / 4.184,
    ))


def printReaction(rxn):
    logger.error('rxn: `{}`\t\torigin: {}'.format(rxn, rxn.get_source()))


def printReactionComments(rxn):
    logger.error('kinetics: {}'.format(rxn.kinetics.comment))


def printSpeciesComments(spc):
    logger.error('thermo: {}'.format(spc.thermo.comment.replace('\n',' ')))


def printKinetics(rxn):
    logger.error('kinetics: `{}`'.format(rxn.kinetics))


def initialize_log(verbose, log_file_name='checkModels.log'):
    """
    Set up a logger for RMG to use to print output to stdout. The
    `verbose` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """
    # since RDKit 2022.03.1, logging is done using the Python logger instead of the
    # Cout streams. This does not affect running RMG normally, but this testing file
    # only works properly if it is the only logger
    # see https://github.com/rdkit/rdkit/pull/4846 for the changes in RDKit

    logging.basicConfig(
        filename=log_file_name,
        filemode='w',
        format='%(message)s',
        level=verbose,
        force=True  # clear all previous log handlers, including those from RDKit
    )

if __name__ == '__main__':
    error = main()
    sys.exit(1 if error else 0)
