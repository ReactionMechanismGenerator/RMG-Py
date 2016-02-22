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

import os.path
import logging
from time import time

from rmgpy.chemkin import getSpeciesIdentifier
from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
import csv
from .loader import loadRMGJob
import rmgpy.util as util 

def simulate(rmg):
    """
    Simulate the RMG job and run the sensitivity analysis if it is on, generating
    output csv files
    """
    util.makeOutputSubdirectory(rmg.outputDirectory, 'solver')

    for index, reactionSystem in enumerate(rmg.reactionSystems):
            
        if reactionSystem.sensitiveSpecies:
            logging.info('Conducting simulation and sensitivity analysis of reaction system %s...' % (index+1))
        
        else:
            logging.info('Conducting simulation of reaction system %s...' % (index+1))
            
        if rmg.saveSimulationProfiles:
            reactionSystem.attach(SimulationProfileWriter(
                rmg.outputDirectory, index, rmg.reactionModel.core.species))   
            reactionSystem.attach(SimulationProfilePlotter(
                    rmg.outputDirectory, index, rmg.reactionModel.core.species))  
        else:
            worksheet = None
            
        sensWorksheet = []
        for spec in reactionSystem.sensitiveSpecies:
            csvfile = file(os.path.join(rmg.outputDirectory, 'solver', 'sensitivity_{0}_SPC_{1}.csv'.format(index+1, spec.index)),'w')
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
            absoluteTolerance = rmg.absoluteTolerance,
            relativeTolerance = rmg.relativeTolerance,
            sensitivity = True if reactionSystem.sensitiveSpecies else False,
            sensitivityAbsoluteTolerance = rmg.sensitivityAbsoluteTolerance,
            sensitivityRelativeTolerance = rmg.sensitivityRelativeTolerance,
            sensWorksheet = sensWorksheet,
        )                      

def runSensitivity(inputFile, chemkinFile, dictFile):
    """
    Runs a standalone simulation of RMG.  Runs sensitivity analysis if sensitive species are given.
    """
    
    rmg = loadRMGJob(inputFile, chemkinFile, dictFile, generateImages=False)    
    
    start_time = time()
    # conduct sensitivity simulation
    simulate(rmg)
    end_time = time()
    time_taken = end_time - start_time
    print "Simulation took {0} seconds".format(time_taken)
