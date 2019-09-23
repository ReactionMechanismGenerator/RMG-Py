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
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
from rmgpy.rmg.main import initialize_log
from rmgpy.rmg.settings import ModelSettings
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.tools.loader import load_rmg_job
from rmgpy.tools.plot import plot_sensitivity


def simulate(rmg, diffusion_limited=True):
    """
    Simulate the RMG job and run the sensitivity analysis if it is on, generating
    output csv files
    diffusion_limited=True implies that if it is a liquid reactor diffusion limitations will be enforced
    otherwise they will not be in a liquid reactor
    """
    util.make_output_subdirectory(rmg.output_directory, 'solver')

    for index, reaction_system in enumerate(rmg.reaction_systems):

        if reaction_system.sensitive_species:
            logging.info('Conducting simulation and sensitivity analysis of reaction system %s...' % (index + 1))
            if reaction_system.sensitive_species == ['all']:
                reaction_system.sensitive_species = rmg.reaction_model.core.species

        else:
            logging.info('Conducting simulation of reaction system %s...' % (index + 1))

        reaction_system.attach(SimulationProfileWriter(
            rmg.output_directory, index, rmg.reaction_model.core.species))
        reaction_system.attach(SimulationProfilePlotter(
            rmg.output_directory, index, rmg.reaction_model.core.species))

        sens_worksheet = []
        for spec in reaction_system.sensitive_species:
            csvfile_path = os.path.join(rmg.output_directory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(index + 1, spec.index))
            sens_worksheet.append(csvfile_path)

        pdep_networks = []
        for source, networks in rmg.reaction_model.network_dict.items():
            pdep_networks.extend(networks)

        model_settings = ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=1)
        simulator_settings = rmg.simulator_settings_list[-1]

        if isinstance(reaction_system, LiquidReactor):
            if diffusion_limited:
                rmg.load_database()
                solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
                diffusion_limiter.enable(solvent_data, rmg.database.solvation)

            # Store constant species indices
            if reaction_system.const_spc_names is not None:
                reaction_system.get_const_spc_indices(rmg.reaction_model.core.species)
        elif rmg.uncertainty is not None:
            rmg.verbose_comments = True
            rmg.load_database()

        reaction_system.simulate(
            core_species=rmg.reaction_model.core.species,
            core_reactions=rmg.reaction_model.core.reactions,
            edge_species=rmg.reaction_model.edge.species,
            edge_reactions=rmg.reaction_model.edge.reactions,
            surface_species=[],
            surface_reactions=[],
            pdep_networks=pdep_networks,
            sensitivity=True if reaction_system.sensitive_species else False,
            sens_worksheet=sens_worksheet,
            model_settings=model_settings,
            simulator_settings=simulator_settings,
        )

        if reaction_system.sensitive_species:
            plot_sensitivity(rmg.output_directory, index, reaction_system.sensitive_species)
            rmg.run_uncertainty_analysis()


def run_simulation(input_file, chemkin_file, dict_file, diffusion_limited=True, check_duplicates=True):
    """
    Runs a standalone simulation of RMG.  Runs sensitivity analysis if sensitive species are given.
    Also runs uncertainty analysis if uncertainty options block is present in input file.

    diffusion_limited=True implies that if it is a liquid reactor diffusion limitations will be enforced
    otherwise they will not be in a liquid reactor
    """
    output_dir = os.path.abspath(os.path.dirname(input_file))
    initialize_log(logging.INFO, os.path.join(output_dir, 'simulate.log'))

    rmg = load_rmg_job(input_file, chemkin_file, dict_file, generate_images=False, check_duplicates=check_duplicates)

    start_time = time()
    # conduct simulation
    simulate(rmg, diffusion_limited)
    end_time = time()
    time_taken = end_time - start_time
    logging.info("Simulation took {0} seconds".format(time_taken))
