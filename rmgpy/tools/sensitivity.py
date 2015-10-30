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
import csv
from time import time

from rmgpy.rmg.main import RMG
from .loader import loadRMGJob

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
    
    rmg = loadRMGJob(inputFile, chemkinFile, dictFile, generateImages=False)    
    
    start_time = time()
    # conduct sensitivity simulation
    simulate(rmg)
    end_time = time()
    time_taken = end_time - start_time
    print "Sensitivity analysis took {0} seconds".format(time_taken)