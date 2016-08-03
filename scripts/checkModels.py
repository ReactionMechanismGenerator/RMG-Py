#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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

import sys
import os
import os.path
import math

import logging
import argparse

from rmgpy.tools.diff_models import execute

logger = logging.getLogger('checkModels')


def parseCommandLineArguments():
        
    
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
    """
    args = parseCommandLineArguments()

    name = args.name[0]
    initializeLog(logging.WARNING, name + '.log')

    benchChemkin = args.benchChemkin[0]
    benchSpeciesDict = args.benchSpeciesDict[0]

    testChemkin = args.testChemkin[0]
    testSpeciesDict = args.testSpeciesDict[0]

    check(name, benchChemkin, benchSpeciesDict, testChemkin, testSpeciesDict)

def check(name, benchChemkin, benchSpeciesDict, testChemkin, testSpeciesDict):
    """
    Compare the tested model to the benchmark model.
    """
    kwargs = {
        'web': True,
        }

    testThermo, benchThermo = None, None
    commonSpecies, uniqueSpeciesTest, uniqueSpeciesOrig, commonReactions, uniqueReactionsTest, uniqueReactionsOrig = \
    execute(benchChemkin, benchSpeciesDict, benchThermo, testChemkin, testSpeciesDict, testThermo, **kwargs)

    errorModel = checkModel(commonSpecies, uniqueSpeciesTest, uniqueSpeciesOrig, commonReactions, uniqueReactionsTest, uniqueReactionsOrig)

    errorSpecies = checkSpecies(commonSpecies, uniqueSpeciesTest, uniqueSpeciesOrig)

    errorReactions = checkReactions(commonReactions, uniqueReactionsTest, uniqueReactionsOrig)

def checkModel(commonSpecies, uniqueSpeciesTest, uniqueSpeciesOrig, commonReactions, uniqueReactionsTest, uniqueReactionsOrig):
    """
    Compare the species and reaction count of both models.
    """

    testModelSpecies = len(commonSpecies) + len(uniqueSpeciesTest)
    origModelSpecies = len(commonSpecies) + len(uniqueSpeciesOrig)

    logger.error('Test model has {} species.'.format(testModelSpecies))
    logger.error('Original model has {} species.'.format(origModelSpecies))

    testModelRxns = len(commonReactions) + len(uniqueReactionsTest)
    origModelRxns = len(commonReactions) + len(uniqueReactionsOrig)
    logger.error('Test model has {} reactions.'.format(testModelRxns))
    logger.error('Original model has {} reactions.'.format(origModelRxns))

    return (testModelSpecies != origModelSpecies) or (testModelRxns != origModelRxns)

def checkSpecies(commonSpecies, uniqueSpeciesTest, uniqueSpeciesOrig):

    error = False

    # check for unique species in one of the models:
    if uniqueSpeciesOrig:
        error = True
        logger.error(
            'The original model has {} species that the tested model does not have.'
            .format(len(uniqueSpeciesOrig))
            )
        [printSpecies(spc) for spc in uniqueSpeciesOrig]

    if uniqueSpeciesTest:
        error = True
        logger.error(
            'The tested model has {} species that the original model does not have.'
            .format(len(uniqueSpeciesTest))
            )
        [printSpecies(spc) for spc in uniqueSpeciesTest]

    # check for different thermo among common species::
    if commonSpecies:
        for spec1, spec2 in commonSpecies:
            logger.info('    {0!s}'.format(spec1))
            if spec1.thermo and spec2.thermo:
                if not spec1.thermo.isSimilarTo(spec2.thermo):
                    error = True
                    logger.error('')
                    logger.error('Non-identical thermo!')
                    logger.error('tested:\t{}'.format(spec1.label))
                    logger.error('original:\t{}'.format(spec2.label))
                    logger.error("{0:10}|{1:10}|{2:10}|{3:10}|{4:10}|{5:10}|{6:10}|{7:10}|{8:10}"
                        .format('Hf(300K)','S(300K)','Cp(300K)','Cp(400K)','Cp(500K)','Cp(600K)','Cp(800K)','Cp(1000K)','Cp(1500K)')
                        )

                    [printThermo(spc) for spc in [spec1, spec2]]

                    if spec1.thermo.comment != spec2.thermo.comment:
                        [printSpeciesComments(spc) for spc in [spec1, spec2]]
                    else:
                        logger.error('Identical thermo comments')

    return error

def checkReactions(commonReactions, uniqueReactionsTest, uniqueReactionsOrig):

    error = False

    # check for unique reactions in one of the models:
    if uniqueReactionsOrig:
        error = True
        
        logger.error(
            'The original model has {} reactions that the tested model does not have.'
            .format(len(uniqueReactionsOrig))
            )
        
        [printReaction(rxn) for rxn in uniqueReactionsOrig]

    if uniqueReactionsTest:
        error = True
        
        logger.error(
            'The tested model has {} reactions that the original model does not have.'
            .format(len(uniqueReactionsTest))
            )

        [printReaction(rxn) for rxn in uniqueReactionsTest]

    if commonReactions:
        for rxn1, rxn2 in commonReactions:
            logger.info('    {0!s}'.format(rxn1))
            if rxn1.kinetics and rxn2.kinetics:
                if not rxn1.kinetics.isSimilarTo(rxn2.kinetics):
                    error = True
                    logger.error('')
                    logger.error('Non-identical kinetics!')
                    logger.error('tested:')
                    printReaction(rxn1)
                    logger.error('benchm:')
                    printReaction(rxn2)
                    
                    logger.error("{0:7}|{1:7}|{2:7}|{3:7}|{4:7}|{5:7}|{6:7}|{7:7}|{8:7}"
                        .format('k(1bar)','300K','400K','500K','600K','800K','1000K','1500K','2000K')
                        )

                    logger.error('')
                    [printRates(rxn) for rxn in [rxn1, rxn2]]

                    logger.error('')
                    [printKinetics(rxn) for rxn in [rxn1, rxn2]]

                    if rxn1.kinetics.comment != rxn2.kinetics.comment:
                        [printReactionComments(rxn) for rxn in [rxn1, rxn2]]
                    else:
                        logger.error('Identical kinetics comments')
                            

    return error

def printSpecies(spc):
    """

    """

    logger.error(
        'spc: {}'.format(spc)
        )        

def printRates(rxn):
    """

    """
    logger.error("{0:7}|{1:7.2f}|{2:7.2f}|{3:7.2f}|{4:7.2f}|{5:7.2f}|{6:7.2f}|{7:7.2f}|{8:7.2f}"
        .format(
        'k(T): ',
        math.log10(rxn.kinetics.getRateCoefficient(300, 1e5)),
        math.log10(rxn.kinetics.getRateCoefficient(400, 1e5)),
        math.log10(rxn.kinetics.getRateCoefficient(500, 1e5)),
        math.log10(rxn.kinetics.getRateCoefficient(600, 1e5)),
        math.log10(rxn.kinetics.getRateCoefficient(800, 1e5)),
        math.log10(rxn.kinetics.getRateCoefficient(1000, 1e5)),
        math.log10(rxn.kinetics.getRateCoefficient(1500, 1e5)),
        math.log10(rxn.kinetics.getRateCoefficient(2000, 1e5)),
    ))

def printThermo(spec):
    """

    """
    logger.error("{0:10.2f}|{1:10.2f}|{2:10.2f}|{3:10.2f}|{4:10.2f}|{5:10.2f}|{6:10.2f}|{7:10.2f}|{8:10.2f}"
        .format(
        spec.thermo.getEnthalpy(300) / 4184.,
        spec.thermo.getEntropy(300) / 4.184,
        spec.thermo.getHeatCapacity(300) / 4.184,
        spec.thermo.getHeatCapacity(400) / 4.184,
        spec.thermo.getHeatCapacity(500) / 4.184,
        spec.thermo.getHeatCapacity(600) / 4.184,
        spec.thermo.getHeatCapacity(800) / 4.184,
        spec.thermo.getHeatCapacity(1000) / 4.184,
        spec.thermo.getHeatCapacity(1500) / 4.184,
    ))

def printReaction(rxn):
    logger.error('rxn: {}\t\torigin: {}'.format(rxn, rxn.getSource()))

def printReactionComments(rxn):
    logger.error('kinetics: {}'.format(rxn.kinetics.comment))

def printSpeciesComments(spc):
    logger.error('thermo: {}'.format(spc.thermo.comment))

def printKinetics(rxn):
    logger.error('Kinetics: {}'.format(rxn.kinetics))

def initializeLog(verbose, log_file_name='checkModels.log'):
    """
    Set up a logger for RMG to use to print output to stdout. The
    `verbose` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """
    logging.basicConfig(
        filename=log_file_name,
        filemode='w',
        format='%(name)s=%(message)s',
        level=verbose
        )

if __name__ == '__main__':
    main()