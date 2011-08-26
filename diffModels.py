#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script can be used to compare two RMG-generated kinetics models. To use,
pass the 
"""

import math
import numpy
import pylab
#import matplotlib.pyplot

from rmgpy.chemkin import loadChemkinFile
from rmgpy.reaction import ReactionModel

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
                break
    uniqueReactions1 = [rxn for rxn in model1.reactions if rxn not in commonReactions.keys()]
    uniqueReactions2 = [rxn for rxn in model2.reactions if rxn not in commonReactions.values()]
    
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
    
################################################################################

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('chemkin1', metavar='CHEMKIN1', type=str, nargs=1,
        help='the Chemkin file of the first model')
    parser.add_argument('speciesDict1', metavar='SPECIESDICT1', type=str, nargs=1,
        help='the species dictionary file of the first model')
    parser.add_argument('chemkin2', metavar='CHEMKIN2', type=str, nargs=1,
        help='the Chemkin file of the second model')
    parser.add_argument('speciesDict2', metavar='SPECIESDICT2', type=str, nargs=1,
        help='the species dictionary file of the second model')
    
    args = parser.parse_args()
    chemkin1 = args.chemkin1[0]
    speciesDict1 = args.speciesDict1[0]
    chemkin2 = args.chemkin2[0]
    speciesDict2 = args.speciesDict2[0]
    
    model1 = ReactionModel()
    model1.species, model1.reactions = loadChemkinFile(chemkin1, speciesDict1)
    model2 = ReactionModel()
    model2.species, model2.reactions = loadChemkinFile(chemkin2, speciesDict2)
    
    compareModelKinetics(model1, model2)
    