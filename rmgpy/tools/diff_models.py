#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
This script can be used to compare two RMG-generated kinetics models. To use,
pass the chem.inp and species_dictionary.txt files to the script. The syntax
is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 CHEMKIN2 SPECIESDICT2

Optionally, you may use the --thermo1 and/or --thermo2 flags to add separate
thermo chemkin files.

The optional --web flag is used for running this script through the RMG-website

With all the above options the syntax is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 --thermo1 THERMO1 CHEMKIN2 SPECIESDICT2 --thermo2 THERMO2 --web

Further option flags:
======================= ====================================================================================
Flag                    Description
======================= ====================================================================================
--diffOnly              Only show species and reactions which are unique or have different values
--commonDiffOnly        Only show species and reactions present in BOTH models which have different values
======================= ====================================================================================
"""
import os
import math
import numpy
import os.path

import logging
import argparse

from rmgpy.chemkin import loadChemkinFile
from rmgpy.rmg.model import ReactionModel
from rmgpy.rmg.output import saveDiffHTML

################################################################################

def compareModelKinetics(model1, model2):
    """
    Compare the kinetics of :class:`ReactionModel` objects `model1` and 
    `model2`, printing the results to stdout.
    """
    from matplotlib import pylab
    # Determine reactions that both models have in common
    commonReactions = {}
    for rxn1 in model1.reactions:
        for rxn2 in model2.reactions:
            if rxn1.isIsomorphic(rxn2):
                commonReactions[rxn1] = rxn2
                model2.reactions.remove(rxn2)
                break
    uniqueReactions1 = [rxn for rxn in model1.reactions if rxn not in commonReactions.keys()]
    uniqueReactions2 = model2.reactions
    
    logging.info('{0:d} reactions were found in both models:'.format(len(commonReactions)))
    for rxn in commonReactions:
        logging.info('    {0!s}'.format(rxn))
    logging.info('{0:d} reactions were only found in the first model:'.format(len(uniqueReactions1)))
    for rxn in uniqueReactions1:
        logging.info('    {0!s}'.format(rxn))
    logging.info('{0:d} reactions were only found in the second model:'.format(len(uniqueReactions2)))
    for rxn in uniqueReactions2:
        logging.info('    {0!s}'.format(rxn))
    
    from rmgpy.kinetics import Chebyshev
    
    T = 1000; P = 1e5
    kinetics1 = []; kinetics2 = []
    for rxn1, rxn2 in commonReactions.iteritems():
        kinetics1.append(rxn1.getRateCoefficient(T,P))
        if rxn1.isIsomorphic(rxn2, eitherDirection=False):
            kinetics2.append(rxn2.getRateCoefficient(T,P))
        else:
            kinetics2.append(rxn2.getRateCoefficient(T,P) / rxn2.getEquilibriumConstant(T))
    fig = pylab.figure(figsize=(8,6))
    ax = pylab.subplot(1,1,1) 
    pylab.loglog(kinetics1, kinetics2, 'o', picker=5)
    xlim = pylab.xlim()
    ylim = pylab.ylim()
    lim = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
    ax.loglog(lim, lim, '-k')
    pylab.xlabel('Model 1 rate coefficient (SI units)')
    pylab.ylabel('Model 2 rate coefficient (SI units)')
    pylab.title('T = {0:g} K, P = {1:g} bar'.format(T, P/1e5))
    pylab.xlim(lim)
    pylab.ylim(lim)
    
    def onpick(event):
        xdata = event.artist.get_xdata()
        ydata = event.artist.get_ydata()
        for ind in event.ind:
            logging.info(commonReactions.keys()[ind])
            logging.info('k(T,P) = {0:9.2e} from model 1'.format(xdata[ind]))
            logging.info('k(T,P) = {0:9.2e} from model 2'.format(ydata[ind]))
            logging.info('ratio = 10**{0:.2f}'.format(math.log10(xdata[ind] / ydata[ind])))
        
    connection_id = fig.canvas.mpl_connect('pick_event', onpick)
        
        
    pylab.show()

def compareModelSpecies(model1, model2):
    """
    This function compares two RMG models and returns a list of common species (with a nested list containing
    both species objects as elements), as well as a list of unique species for each model.
    """

    commonSpecies = []
    uniqueSpecies1 = model1.species[:]
    uniqueSpecies2 = []
    
    for spec2 in model2.species:
        for spec1 in uniqueSpecies1[:]: # make a copy so you don't remove from the list you are iterating over
            if spec1.isIsomorphic(spec2):
                commonSpecies.append([spec1, spec2])
                uniqueSpecies1.remove(spec1)
                break
        else:
            uniqueSpecies2.append(spec2)
    # Remove species in the mechanism that aren't identified (includes those called out as species
    # but not used)        
    for spec in uniqueSpecies1[:]: # make a copy so you don't remove from the list you are iterating over
        if not len(spec.molecule):
            uniqueSpecies1.remove(spec)
            logging.warning("Removing species {!r} from model 1 because it has no molecule info".format(spec))
    for spec in uniqueSpecies2[:]: # make a copy so you don't remove from the list you are iterating over
        if not spec.molecule:
            uniqueSpecies2.remove(spec)
            logging.warning("Removing species {!r} from model 2 because it has no molecule info".format(spec))
    return commonSpecies, uniqueSpecies1, uniqueSpecies2

def compareModelReactions(model1, model2):
    """
    This function compares two RMG models and returns a list of common reactions (with a nested list containing
    both reaction objects as elements), as well as a list of unique reactions for each model.
    """
    reactionList1 = model1.reactions[:]
    reactionList2 = model2.reactions[:]
    
    # remove reactions that have an unidentified species
    to_remove = []
    for reactionList in (reactionList1, reactionList2):
        for reaction in reactionList:
            for side in (reaction.products, reaction.reactants):
                for species in side:
                    if not species.molecule:
                        to_remove.append((reactionList,reaction))
                        logging.warning("Removing reaction {!r} that had unidentified species {!r}".format(reaction, species))
                        break
    for reactionList, reaction in to_remove:
        reactionList.remove(reaction)
    
    commonReactions = []; uniqueReactions1 = []; uniqueReactions2 = []
    for rxn1 in reactionList1:
        for rxn2 in reactionList2[:]: # make a copy so you don't remove from the list you are iterating over
            if rxn1.isIsomorphic(rxn2):
                commonReactions.append([rxn1, rxn2])
                # Remove reaction 2 from being chosen a second time.
                # Let each reaction only appear only once in the diff comparison.
                # Otherwise this miscounts number of reactions in model 2.
                reactionList2.remove(rxn2)
                break
    for rxn1 in reactionList1:
        for r1, r2 in commonReactions:
            if rxn1 is r1:
                break
        else:
            uniqueReactions1.append(rxn1)
    for rxn2 in reactionList2:
        for r1, r2 in commonReactions:
            if rxn2 is r2:
                break
        else:
            uniqueReactions2.append(rxn2)

    return commonReactions, uniqueReactions1, uniqueReactions2

def saveCompareHTML(outputDir,chemkinPath1,speciesDictPath1,chemkinPath2,speciesDictPath2,readComments1=True,readComments2=True):
    """
    Saves a model comparison HTML file based on two sets of chemkin and species dictionary
    files.
    """
    model1 = ReactionModel()
    model1.species, model1.reactions = loadChemkinFile(chemkinPath1, speciesDictPath1, readComments = readComments1)
    model2 = ReactionModel()
    model2.species, model2.reactions = loadChemkinFile(chemkinPath2, speciesDictPath2, readComments = readComments2)
    commonReactions, uniqueReactions1, uniqueReactions2 = compareModelReactions(model1, model2)
    commonSpecies, uniqueSpecies1, uniqueSpecies2 = compareModelSpecies(model1, model2)
    
    outputPath = outputDir + 'diff.html'            
    saveDiffHTML(outputPath, commonSpecies, uniqueSpecies1, uniqueSpecies2, commonReactions, uniqueReactions1, uniqueReactions2)

def enthalpyDiff(species):
    """
    Returns the enthalpy discrepancy between the same species in the two models
    """
    thermo0 = species[0].thermo
    thermo1 = species[1].thermo
    if thermo0 and thermo1:    
        diff =  species[0].thermo.discrepancy(species[1].thermo)
    else:
        diff = 99999999
    return -1*diff

def kineticsDiff(reaction):
    """
    Returns some measure of the discrepancy between two reactions in a model
    """    
    kinetics0 = reaction[0].kinetics
    kinetics1 = reaction[1].kinetics
    if kinetics0 and kinetics1:
        diff = reaction[0].kinetics.discrepancy(reaction[1].kinetics)        
    else:
        diff = 9999999
    return -1*diff

def identicalThermo(species_pair):
    return species_pair[0].thermo.isIdenticalTo(species_pair[1].thermo)

def identicalKinetics(reaction_pair):
    return reaction_pair[0].kinetics.isIdenticalTo(reaction_pair[1].kinetics)

################################################################################

def parseCommandLineArguments():
        
    
    parser = argparse.ArgumentParser()
    parser.add_argument('chemkin1', metavar='CHEMKIN1', type=str, nargs=1,
        help='the Chemkin file of the first model')
    parser.add_argument('speciesDict1', metavar='SPECIESDICT1', type=str, nargs=1,
        help='the species dictionary file of the first model')
    parser.add_argument('--thermo1', metavar = 'THERMO1', type=str, nargs = 1,
        help = 'the thermo file of the first model')
    parser.add_argument('chemkin2', metavar='CHEMKIN2', type=str, nargs=1,
        help='the Chemkin file of the second model')
    parser.add_argument('speciesDict2', metavar='SPECIESDICT2', type=str, nargs=1,
        help='the species dictionary file of the second model')
    parser.add_argument('--thermo2', metavar = 'THERMO2', type=str, nargs = 1,
        help = 'the thermo file of the second model')
    parser.add_argument('--web', action='store_true', help='Running diff models through the RMG-website')
    parser.add_argument('--diffOnly', action='store_true', help='Do not show identical species thermo or reactions')
    parser.add_argument('--commonDiffOnly', action='store_true',
        help='Only show species and reactions present in BOTH models which have different values')
    
    args = parser.parse_args()

    return args

def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.
    """
    args = parseCommandLineArguments()


    chemkin1 = args.chemkin1[0]
    speciesDict1 = args.speciesDict1[0]
    if args.thermo1: 
        thermo1 = args.thermo1[0]
    else:
        thermo1 = None
    chemkin2 = args.chemkin2[0]
    speciesDict2 = args.speciesDict2[0]
    if args.thermo2: 
        thermo2 = args.thermo2[0]
    else:
        thermo2 = None

    kwargs = {
            'web': args.web,
            'wd': os.getcwd(),
            'diffOnly': args.diffOnly,
            'commonDiffOnly': args.commonDiffOnly,
            }

    execute(chemkin1, speciesDict1, thermo1, chemkin2, speciesDict2, thermo2, **kwargs)

def execute(chemkin1, speciesDict1, thermo1, chemkin2, speciesDict2, thermo2, **kwargs):
    
    model1 = ReactionModel()
    model1.species, model1.reactions = loadChemkinFile(chemkin1, speciesDict1, thermoPath = thermo1)
    model2 = ReactionModel()
    model2.species, model2.reactions = loadChemkinFile(chemkin2, speciesDict2, thermoPath = thermo2)
    
    commonSpecies, uniqueSpecies1, uniqueSpecies2 = compareModelSpecies(model1, model2)
    commonReactions, uniqueReactions1, uniqueReactions2 = compareModelReactions(model1, model2)

    try:
        diffOnly = kwargs['diffOnly']
    except KeyError:
        diffOnly = False

    try:
        commonDiffOnly = kwargs['commonDiffOnly']
    except KeyError:
        commonDiffOnly = False

    if diffOnly or commonDiffOnly:
        commonSpecies = filter(lambda x: not identicalThermo(x), commonSpecies)
        commonReactions = filter(lambda x: not identicalKinetics(x), commonReactions)

    if commonDiffOnly:
        uniqueSpecies1 = []
        uniqueSpecies2 = []
        uniqueReactions1 = []
        uniqueReactions2 = []
    
    try:
        web = kwargs['web']
    except KeyError:
        web = False

    if not web:
        logging.info('{0:d} species were found in both models:'.format(len(commonSpecies)))
        for spec1, spec2 in commonSpecies:
            logging.info('    {0!s}'.format(spec1))
            if spec1.thermo and spec2.thermo:
                spec1.molecule[0].getSymmetryNumber()
                logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f} {8:7.2f}'.format(
                    spec1.thermo.getEnthalpy(300) / 4184.,
                    spec1.thermo.getEntropy(300) / 4.184,
                    spec1.thermo.getHeatCapacity(300) / 4.184,
                    spec1.thermo.getHeatCapacity(400) / 4.184,
                    spec1.thermo.getHeatCapacity(500) / 4.184,
                    spec1.thermo.getHeatCapacity(600) / 4.184,
                    spec1.thermo.getHeatCapacity(800) / 4.184,
                    spec1.thermo.getHeatCapacity(1000) / 4.184,
                    spec1.thermo.getHeatCapacity(1500) / 4.184,
                ))
                logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f} {8:7.2f}'.format(
                    spec2.thermo.getEnthalpy(300) / 4184.,
                    spec2.thermo.getEntropy(300) / 4.184,
                    spec2.thermo.getHeatCapacity(300) / 4.184,
                    spec2.thermo.getHeatCapacity(400) / 4.184,
                    spec2.thermo.getHeatCapacity(500) / 4.184,
                    spec2.thermo.getHeatCapacity(600) / 4.184,
                    spec2.thermo.getHeatCapacity(800) / 4.184,
                    spec2.thermo.getHeatCapacity(1000) / 4.184,
                    spec2.thermo.getHeatCapacity(1500) / 4.184,
                ))
        logging.info('{0:d} species were only found in the first model:'.format(len(uniqueSpecies1)))
        for spec in uniqueSpecies1:
            logging.info('    {0!s}'.format(spec))
        logging.info('{0:d} species were only found in the second model:'.format(len(uniqueSpecies2)))
        for spec in uniqueSpecies2:
            logging.info('    {0!s}'.format(spec))
    
        logging.info('{0:d} reactions were found in both models:'.format(len(commonReactions)))
        for rxn1, rxn2 in commonReactions:
            logging.info('    {0!s}'.format(rxn1))
            if rxn1.kinetics and rxn2.kinetics:
                logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                    math.log10(rxn1.kinetics.getRateCoefficient(300, 1e5)),
                    math.log10(rxn1.kinetics.getRateCoefficient(400, 1e5)),
                    math.log10(rxn1.kinetics.getRateCoefficient(500, 1e5)),
                    math.log10(rxn1.kinetics.getRateCoefficient(600, 1e5)),
                    math.log10(rxn1.kinetics.getRateCoefficient(800, 1e5)),
                    math.log10(rxn1.kinetics.getRateCoefficient(1000, 1e5)),
                    math.log10(rxn1.kinetics.getRateCoefficient(1500, 1e5)),
                    math.log10(rxn1.kinetics.getRateCoefficient(2000, 1e5)),
                ))
                logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                    math.log10(rxn2.kinetics.getRateCoefficient(300, 1e5)),
                    math.log10(rxn2.kinetics.getRateCoefficient(400, 1e5)),
                    math.log10(rxn2.kinetics.getRateCoefficient(500, 1e5)),
                    math.log10(rxn2.kinetics.getRateCoefficient(600, 1e5)),
                    math.log10(rxn2.kinetics.getRateCoefficient(800, 1e5)),
                    math.log10(rxn2.kinetics.getRateCoefficient(1000, 1e5)),
                    math.log10(rxn2.kinetics.getRateCoefficient(1500, 1e5)),
                    math.log10(rxn2.kinetics.getRateCoefficient(2000, 1e5)),
                ))
        logging.info('{0:d} reactions were only found in the first model:'.format(len(uniqueReactions1)))
        for rxn in uniqueReactions1:
            logging.info('    {0!s}'.format(rxn))
        logging.info('{0:d} reactions were only found in the second model:'.format(len(uniqueReactions2)))
        for rxn in uniqueReactions2:
            logging.info('    {0!s}'.format(rxn))

    logging.info("Saving output in diff.html")

    try:
        wd = kwargs['wd']
    except KeyError:
        wd = os.getcwd()

    outputPath = os.path.join(wd, 'diff.html')
    saveDiffHTML(outputPath, commonSpecies, uniqueSpecies1, uniqueSpecies2, commonReactions, uniqueReactions1, uniqueReactions2)
    logging.info("Finished!")

    return commonSpecies, uniqueSpecies1, uniqueSpecies2, commonReactions, uniqueReactions1, uniqueReactions2
