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

import logging
import os.path
from time import time

import rmgpy.util as util
from rmgpy.kinetics.diffusionLimited import diffusionLimiter
from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
from rmgpy.rmg.main import initializeLog
from rmgpy.rmg.settings import ModelSettings
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.tools.loader import loadRMGJob
from rmgpy.tools.plot import plot_sensitivity


def simulate(rmg, diffusionLimited=True):
    """
    Simulate the RMG job and run the sensitivity analysis if it is on, generating
    output csv files
    diffusionLimited=True implies that if it is a liquid reactor diffusion limitations will be enforced
    otherwise they will not be in a liquid reactor
    """
    util.make_output_subdirectory(rmg.outputDirectory, 'solver')

    for index, reactionSystem in enumerate(rmg.reactionSystems):

        if reactionSystem.sensitiveSpecies:
            logging.info('Conducting simulation and sensitivity analysis of reaction system %s...' % (index + 1))
            if reactionSystem.sensitiveSpecies == ['all']:
                reactionSystem.sensitiveSpecies = rmg.reactionModel.core.species

        else:
            logging.info('Conducting simulation of reaction system %s...' % (index + 1))

        reactionSystem.attach(SimulationProfileWriter(
            rmg.outputDirectory, index, rmg.reactionModel.core.species))
        reactionSystem.attach(SimulationProfilePlotter(
            rmg.outputDirectory, index, rmg.reactionModel.core.species))

        sens_worksheet = []
        for spec in reactionSystem.sensitiveSpecies:
            csvfile_path = os.path.join(rmg.outputDirectory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(index + 1, spec.index))
            sens_worksheet.append(csvfile_path)

        pdep_networks = []
        for source, networks in rmg.reactionModel.networkDict.items():
            pdep_networks.extend(networks)

        model_settings = ModelSettings(toleranceKeepInEdge=0, toleranceMoveToCore=1, toleranceInterruptSimulation=1)
        simulator_settings = rmg.simulatorSettingsList[-1]

        if isinstance(reactionSystem, LiquidReactor):
            if diffusionLimited:
                rmg.loadDatabase()
                solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
                diffusionLimiter.enable(solvent_data, rmg.database.solvation)

            # Store constant species indices
            if reactionSystem.constSPCNames is not None:
                reactionSystem.get_constSPCIndices(rmg.reactionModel.core.species)
        elif rmg.uncertainty is not None:
            rmg.verboseComments = True
            rmg.loadDatabase()

        reactionSystem.simulate(
            coreSpecies=rmg.reactionModel.core.species,
            coreReactions=rmg.reactionModel.core.reactions,
            edgeSpecies=rmg.reactionModel.edge.species,
            edgeReactions=rmg.reactionModel.edge.reactions,
            surfaceSpecies=[],
            surfaceReactions=[],
            pdepNetworks=pdep_networks,
            sensitivity=True if reactionSystem.sensitiveSpecies else False,
            sensWorksheet=sens_worksheet,
            modelSettings=model_settings,
            simulatorSettings=simulator_settings,
        )

        if reactionSystem.sensitiveSpecies:
            plot_sensitivity(rmg.outputDirectory, index, reactionSystem.sensitiveSpecies)
            rmg.run_uncertainty_analysis()


def run_simulation(inputFile, chemkinFile, dictFile, diffusionLimited=True, checkDuplicates=True):
    """
    Runs a standalone simulation of RMG.  Runs sensitivity analysis if sensitive species are given.
    Also runs uncertainty analysis if uncertainty options block is present in input file.

    diffusionLimited=True implies that if it is a liquid reactor diffusion limitations will be enforced
    otherwise they will not be in a liquid reactor
    """
    output_dir = os.path.abspath(os.path.dirname(inputFile))
    initializeLog(logging.INFO, os.path.join(output_dir, 'simulate.log'))

    rmg = loadRMGJob(inputFile, chemkinFile, dictFile, generateImages=False, checkDuplicates=checkDuplicates)

    start_time = time()
    # conduct simulation
    simulate(rmg, diffusionLimited)
    end_time = time()
    time_taken = end_time - start_time
    logging.info("Simulation took {0} seconds".format(time_taken))
