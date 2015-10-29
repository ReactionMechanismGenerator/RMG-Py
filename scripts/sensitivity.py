#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script runs stand-alone sensitivity analysis on an RMG job.
"""

import os.path
import logging
import csv
from time import time

from rmgpy.rmg.main import RMG
from rmgpy.tools.loader import loadRMGJob

################################################################################

def simulate(rmg):
    """
    Simulate the RMG job and run the sensitivity analysis if it is on, generating
    output csv files
    """
        
    for index, reactionSystem in enumerate(rmg.reactionSystems):
            
        if reactionSystem.sensitiveSpecies:
            logging.info('Conducting sensitivity analysis of reaction system %s...' % (index+1))
            
            if rmg.saveSimulationProfiles:
                csvfile = file(os.path.join(rmg.outputDirectory, 'simulation_{0}.csv'.format(index+1)),'w')
                worksheet = csv.writer(csvfile)
            else:
                worksheet = None
                
            sensWorksheet = []
            for spec in reactionSystem.sensitiveSpecies:
                csvfile = file(os.path.join(rmg.outputDirectory, 'sensitivity_{0}_SPC_{1}.csv'.format(index+1, spec.index)),'w')
                sensWorksheet.append(csv.writer(csvfile))
    
            pdepNetworks = []
            for source, networks in rmg.reactionModel.networkDict.items():
                pdepNetworks.extend(networks)
            terminated, obj = reactionSystem.simulate(
                coreSpecies = rmg.reactionModel.core.species,
                coreReactions = rmg.reactionModel.core.reactions,
                edgeSpecies = rmg.reactionModel.edge.species,
                edgeReactions = rmg.reactionModel.edge.reactions,
                toleranceKeepInEdge = 0,
                toleranceMoveToCore = 1,
                toleranceInterruptSimulation = 1,
                pdepNetworks = pdepNetworks,
                worksheet = worksheet,
                absoluteTolerance = rmg.absoluteTolerance,
                relativeTolerance = rmg.relativeTolerance,
                sensitivity = True,
                sensitivityAbsoluteTolerance = rmg.sensitivityAbsoluteTolerance,
                sensitivityRelativeTolerance = rmg.sensitivityRelativeTolerance,
                sensWorksheet = sensWorksheet,
            )                      


################################################################################

def runSensitivity(inputFile, chemkinFile, dictFile):
    # Load the RMG job
    rmg = loadRMGPyJob(inputFile, chemkinFile, dictFile, generateImages=False)    
    
    start_time = time()
    # conduct sensitivity simulation
    simulate(rmg)
    end_time = time()
    time_taken = end_time - start_time
    print "Sensitivity analysis took {0} seconds".format(time_taken)

################################################################################

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='RMG input file')
    parser.add_argument('chemkin', metavar='CHEMKIN', type=str, nargs=1,
        help='Chemkin file')
    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, nargs=1,
        help='RMG dictionary file')
    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    chemkinFile = os.path.abspath(args.chemkin[0])
    dictFile = os.path.abspath(args.dictionary[0])
    
    runSensitivity(inputFile, chemkinFile, dictFile)