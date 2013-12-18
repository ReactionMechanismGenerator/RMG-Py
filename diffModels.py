#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script can be used to compare two RMG-generated kinetics models. To use,
pass the 
"""

import math
import numpy
from matplotlib import pylab
import os.path
#import matplotlib.pyplot
import logging

from rmgpy.chemkin import loadChemkinFile
from rmgpy.reaction import ReactionModel
from rmgpy.rmg.output import saveDiffHTML

################################################################################

def compareModelKinetics(model1, model2):
    """
    Compare the kinetics of :class:`ReactionModel` objects `model1` and 
    `model2`, printing the results to stdout.
    """
    
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
    
    print '{0:d} reactions were found in both models:'.format(len(commonReactions))
    for rxn in commonReactions:
        print '    {0!s}'.format(rxn)
    print '{0:d} reactions were only found in the first model:'.format(len(uniqueReactions1))
    for rxn in uniqueReactions1:
        print '    {0!s}'.format(rxn)
    print '{0:d} reactions were only found in the second model:'.format(len(uniqueReactions2))
    for rxn in uniqueReactions2:
        print '    {0!s}'.format(rxn)
    
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
            print commonReactions.keys()[ind]
            print 'k(T,P) = {0:9.2e} from model 1'.format(xdata[ind])
            print 'k(T,P) = {0:9.2e} from model 2'.format(ydata[ind])
            print 'ratio = 10**{0:.2f}'.format(math.log10(xdata[ind] / ydata[ind]))
        
    connection_id = fig.canvas.mpl_connect('pick_event', onpick)
        
        
    pylab.show()

def compareModelSpecies(model1, model2):
    """
    This function compares two RMG models and returns a list of common reactions
    as a dictionary, as well as a list of unique reactions for each model.
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
    This function compares two RMG models and returns a list of common reactions
    as a dictionary, as well as a list of unique reactions for each model.
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
################################################################################

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('chemkin1', metavar='CHEMKIN1', type=str, nargs=1,
        help='the Chemkin file of the first model')
    parser.add_argument('speciesDict1', metavar='SPECIESDICT1', type=str, nargs=1,
        help='the species dictionary file of the first model')
    parser.add_argument('thermo1', metavar = 'THERMO1', type=str, nargs = 1,
        help = 'the thermo file of the first model')
    parser.add_argument('chemkin2', metavar='CHEMKIN2', type=str, nargs=1,
        help='the Chemkin file of the second model')
    parser.add_argument('speciesDict2', metavar='SPECIESDICT2', type=str, nargs=1,
        help='the species dictionary file of the second model')
    parser.add_argument('thermo2', metavar = 'THERMO2', type=str, nargs = 1,
        help = 'the thermo file of the second model')
    
    args = parser.parse_args()
    chemkin1 = args.chemkin1[0]
    speciesDict1 = args.speciesDict1[0]
    thermo1 = args.thermo1[0]
    chemkin2 = args.chemkin2[0]
    speciesDict2 = args.speciesDict2[0]
    thermo2 = args.thermo2[0]
    
    model1 = ReactionModel()
    model1.species, model1.reactions = loadChemkinFile(chemkin1, speciesDict1, thermoPath = thermo1)
    model2 = ReactionModel()
    model2.species, model2.reactions = loadChemkinFile(chemkin2, speciesDict2, thermoPath = thermo2)
    
    commonSpecies, uniqueSpecies1, uniqueSpecies2 = compareModelSpecies(model1, model2)
    commonReactions, uniqueReactions1, uniqueReactions2 = compareModelReactions(model1, model2)

    print '{0:d} species were found in both models:'.format(len(commonSpecies))
    for spec1, spec2 in commonSpecies:
        print '    {0!s}'.format(spec1)
        if spec1.thermo and spec2.thermo:
            spec1.molecule[0].calculateSymmetryNumber()
            print '        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f} {8:7.2f}'.format( 
                spec1.thermo.getEnthalpy(300) / 4184.,
                spec1.thermo.getEntropy(300) / 4.184,
                spec1.thermo.getHeatCapacity(300) / 4.184,
                spec1.thermo.getHeatCapacity(400) / 4.184,
                spec1.thermo.getHeatCapacity(500) / 4.184,
                spec1.thermo.getHeatCapacity(600) / 4.184,
                spec1.thermo.getHeatCapacity(800) / 4.184,
                spec1.thermo.getHeatCapacity(1000) / 4.184,
                spec1.thermo.getHeatCapacity(1500) / 4.184,
            )
            print '        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f} {8:7.2f}'.format( 
                spec2.thermo.getEnthalpy(300) / 4184.,
                spec2.thermo.getEntropy(300) / 4.184,
                spec2.thermo.getHeatCapacity(300) / 4.184,
                spec2.thermo.getHeatCapacity(400) / 4.184,
                spec2.thermo.getHeatCapacity(500) / 4.184,
                spec2.thermo.getHeatCapacity(600) / 4.184,
                spec2.thermo.getHeatCapacity(800) / 4.184,
                spec2.thermo.getHeatCapacity(1000) / 4.184,
                spec2.thermo.getHeatCapacity(1500) / 4.184,
            )
    print '{0:d} species were only found in the first model:'.format(len(uniqueSpecies1))
    for spec in uniqueSpecies1:
        print '    {0!s}'.format(spec)
    print '{0:d} species were only found in the second model:'.format(len(uniqueSpecies2))
    for spec in uniqueSpecies2:
        print '    {0!s}'.format(spec)

    print '{0:d} reactions were found in both models:'.format(len(commonReactions))
    for rxn1, rxn2 in commonReactions:
        print '    {0!s}'.format(rxn1)
        if rxn1.kinetics and rxn2.kinetics:
            print '        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                math.log10(rxn1.kinetics.getRateCoefficient(300, 1e5)),
                math.log10(rxn1.kinetics.getRateCoefficient(400, 1e5)),
                math.log10(rxn1.kinetics.getRateCoefficient(500, 1e5)),
                math.log10(rxn1.kinetics.getRateCoefficient(600, 1e5)),
                math.log10(rxn1.kinetics.getRateCoefficient(800, 1e5)),
                math.log10(rxn1.kinetics.getRateCoefficient(1000, 1e5)),
                math.log10(rxn1.kinetics.getRateCoefficient(1500, 1e5)),
                math.log10(rxn1.kinetics.getRateCoefficient(2000, 1e5)),
            )
            print '        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                math.log10(rxn2.kinetics.getRateCoefficient(300, 1e5)),
                math.log10(rxn2.kinetics.getRateCoefficient(400, 1e5)),
                math.log10(rxn2.kinetics.getRateCoefficient(500, 1e5)),
                math.log10(rxn2.kinetics.getRateCoefficient(600, 1e5)),
                math.log10(rxn2.kinetics.getRateCoefficient(800, 1e5)),
                math.log10(rxn2.kinetics.getRateCoefficient(1000, 1e5)),
                math.log10(rxn2.kinetics.getRateCoefficient(1500, 1e5)),
                math.log10(rxn2.kinetics.getRateCoefficient(2000, 1e5)),
            )
    print '{0:d} reactions were only found in the first model:'.format(len(uniqueReactions1))
    for rxn in uniqueReactions1:
        print '    {0!s}'.format(rxn)
    print '{0:d} reactions were only found in the second model:'.format(len(uniqueReactions2))
    for rxn in uniqueReactions2:
        print '    {0!s}'.format(rxn)
    
    #commonSpecies.sort(key = enthalpyDiff)
    #commonReactions.sort(key = kineticsDiff)

    print "Saving output in diff.html"
    outputPath = 'diff.html'
    saveDiffHTML(outputPath, commonSpecies, uniqueSpecies1, uniqueSpecies2, commonReactions, uniqueReactions1, uniqueReactions2)
    print "Finished!"
