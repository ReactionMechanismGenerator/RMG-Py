#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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
from rmgpy.tools.plot import ReactionSensitivityPlot, ThermoSensitivityPlot
from rmgpy.rmg.RMGSettings import ModelSettings

def plotSensitivity(outputDirectory, reactionSystemIndex, sensitiveSpeciesList, number=10, fileformat='.png'):
    """
    A function for plotting the top reaction thermo sensitivities (the number is 
    inputted as the variable `number`) in bar plot format.
    To be called after running a simulation on a particular reactionSystem.
    """
    
    for species in sensitiveSpeciesList:
        csvFile = os.path.join(
            outputDirectory,
            'solver',
            'sensitivity_{0}_SPC_{1}.csv'.format(
                reactionSystemIndex + 1, species.index
                )
            )
        
        reactionPlotFile = os.path.join(
            outputDirectory,
            'solver',
            'sensitivity_{0}_SPC_{1}_reactions'.format(
                reactionSystemIndex + 1, species.index
                ) + fileformat
            )
        
        thermoPlotFile = os.path.join(
            outputDirectory,
            'solver',
            'sensitivity_{0}_SPC_{1}_thermo'.format(
                reactionSystemIndex + 1, species.index
                ) + fileformat
            )

        ReactionSensitivityPlot(csvFile=csvFile, numReactions=number).barplot(reactionPlotFile)
        ThermoSensitivityPlot(csvFile=csvFile, numSpecies=number).barplot(thermoPlotFile)



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
            csvfilePath = os.path.join(rmg.outputDirectory, 'solver', 'sensitivity_{0}_SPC_{1}.csv'.format(index+1, spec.index))
            sensWorksheet.append(csvfilePath)

        pdepNetworks = []
        for source, networks in rmg.reactionModel.networkDict.items():
            pdepNetworks.extend(networks)
        
        modelSettings = ModelSettings(toleranceKeepInEdge = 0, toleranceMoveToCore = 1, toleranceInterruptSimulation = 1)
        simulatorSettings = rmg.simulatorSettingsList[-1]
        
        terminated, obj,sspcs,srxns = reactionSystem.simulate(
            coreSpecies = rmg.reactionModel.core.species,
            coreReactions = rmg.reactionModel.core.reactions,
            edgeSpecies = rmg.reactionModel.edge.species,
            edgeReactions = rmg.reactionModel.edge.reactions,
            surfaceSpecies = [],
            surfaceReactions = [],
            pdepNetworks = pdepNetworks,
            sensitivity = True if reactionSystem.sensitiveSpecies else False,
            sensWorksheet = sensWorksheet,
            modelSettings = modelSettings,
            simulatorSettings = simulatorSettings,
        )
        
        if reactionSystem.sensitiveSpecies:
            plotSensitivity(rmg.outputDirectory, index, reactionSystem.sensitiveSpecies)

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
