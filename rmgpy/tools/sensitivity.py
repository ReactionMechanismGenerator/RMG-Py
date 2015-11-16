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
import numpy as np

from rmgpy.chemkin import getSpeciesIdentifier
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

class SimulationProfileWriter(object):
    """
        SimulationProfileWriter listens to a ReactionSystem subject
        and writes the species mole fractions as a function of the reaction time
        to a csv file.
    """
    def __init__(self, outputDirectory, reaction_sys_index, coreSpecies):
        super(SimulationProfileWriter, self).__init__()
        
        self.outputDirectory = outputDirectory
        self.reaction_sys_index = reaction_sys_index
        self.coreSpecies = coreSpecies

    def update(self, reactionSystem):
        """
        Opens a file with filename referring to:
            - reaction system
            - number of core species

        Writes to a csv file:
            - header row with species names
            - each row with mole fractions of the core species in the given reaction system.
        """

        filename = os.path.join(
            self.outputDirectory,
            'solver',
            'simulation_{0}_{1:d}.csv'.format(
                self.reaction_sys_index + 1, len(self.coreSpecies)
                )
            )

        header = ['Time (s)', 'Volume (m^3)']
        for spc in self.coreSpecies:
            header.append(getSpeciesIdentifier(spc))

        with open(filename, 'w') as csvfile:
            worksheet = csv.writer(csvfile)

            # add header row:
            worksheet.writerow(header) 

            # add mole fractions:
            worksheet.writerows(reactionSystem.snapshots)
            
                
        