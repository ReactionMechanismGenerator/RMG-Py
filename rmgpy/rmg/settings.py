#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains settings classes for manipulation of RMG run parameters
==================================================================================================================================================
    `atol`                                    The absolute tolerance used in the ODE/DAE solver
    `rtol`                                    The relative tolerance used in the ODE/DAE solver
    `sens_atol`                               The absolute tolerance used in the ODE/DAE solver for the sensitivities
    `sens_rtol`                               The relative tolerance used in the ODE/DAE solver for the sensitivities
    `tol_keep_in_edge`                        The relative species flux below which species are discarded from the edge
    `tol_move_to_core`                        The relative species flux above which species are moved from the edge to the core
    `tol_interrupt_simulation`                The relative species flux above which the simulation will halt
    `tol_move_edge_rxn_to_core`               The dynamics number above which reactions are moved from the edge to the core
    `tol_move_edge_rxn_to_core_interrupt`     The edge dynamics number above which the simulation will halt
    `tol_move_edge_rxn_to_surface`            The dynamics number above which reactions are moved from the edge to the surface
    `tol_move_edge_rxn_to_surface_interrupt`  The edge dynamics number above which simulation will halt
    `tol_move_surface_rxn_to_core`            The surface dynamics number above which reactions are moved from the surface to the bulk core
    `tol_move_surface_spc_to_core`            The relative species flux above which species are moved from the surface to the bulk core
    `maximum_edge_species`                    The maximum number of edge species allowed at any time
    `min_core_size_for_prune`                 Minimum number of core species before pruning is allowed
    `min_species_exist_iterations_for_prune`  Minimum number of iterations a species must exist before it can be pruned
    `filter_reactions`                        Specify whether to filter reactions during model enlarging step
    `filter_threshold`                        Bimolecular reaction filtering threshold rate constant
    `ignore_overall_flux_criterion`           flag indicating that the ordinary flux criterion should be ignored except for pdep purposes
    `max_num_species`                         Number of core species at which a stage/job will terminate
    `maxNumObjPerIter`                        Maximum number of objects that can be sent for enlargement from a single simulation
    `transitory_tol_dict`                     Dictionary mapping species names to transitory sensitivity tolerances
    `transitory_step_period`                  Number of steps between transitory edge analyses
==================================================================================================================================================
"""
import numpy as np

from rmgpy.quantity import Quantity


class ModelSettings(object):
    """
    class for holding the parameters affecting an RMG run
    """

    def __init__(self, tol_move_to_core=None, tol_move_edge_rxn_to_core=np.inf, tol_keep_in_edge=0.0,
                 tol_interrupt_simulation=1.0,
                 tol_move_edge_rxn_to_surface=np.inf, tol_move_surface_spc_to_core=np.inf,
                 tol_move_surface_rxn_to_core=np.inf,
                 tol_move_edge_rxn_to_surface_interrupt=None, tol_move_edge_rxn_to_core_interrupt=None,
                 maximum_edge_species=1000000, min_core_size_for_prune=50,
                 min_species_exist_iterations_for_prune=2, filter_reactions=False, filter_threshold=1e8,
                 ignore_overall_flux_criterion=False, max_num_species=None, max_num_objects_per_iter=1,
                 terminate_at_max_objects=False, thermo_tol_keep_spc_in_edge=np.inf,
                 dynamics_time_scale=Quantity((0.0, 'sec')),
                 tol_branch_rxn_to_core=0.0, branching_index=0.5, branching_ratio_max=1.0,transitory_tol_dict=dict(),
                 transitory_step_period=20):

        self.tol_keep_in_edge = tol_keep_in_edge
        self.tol_move_to_core = tol_move_to_core
        self.tol_move_edge_rxn_to_core = tol_move_edge_rxn_to_core
        self.tol_interrupt_simulation = tol_interrupt_simulation
        self.maximum_edge_species = maximum_edge_species
        self.min_core_size_for_prune = min_core_size_for_prune
        self.min_species_exist_iterations_for_prune = min_species_exist_iterations_for_prune
        self.filter_reactions = filter_reactions
        self.filter_threshold = filter_threshold
        self.ignore_overall_flux_criterion = ignore_overall_flux_criterion
        self.tol_move_edge_rxn_to_surface = tol_move_edge_rxn_to_surface
        self.tol_move_surface_spc_to_core = tol_move_surface_spc_to_core
        self.tol_move_surface_rxn_to_core = tol_move_surface_rxn_to_core
        self.thermo_tol_keep_spc_in_edge = thermo_tol_keep_spc_in_edge
        self.terminate_at_max_objects = terminate_at_max_objects
        self.dynamics_time_scale = dynamics_time_scale.value_si
        self.tol_branch_rxn_to_core = tol_branch_rxn_to_core
        self.branching_index = branching_index
        self.branching_ratio_max = branching_ratio_max
        self.transitory_tol_dict = transitory_tol_dict
        self.transitory_step_period = transitory_step_period

        if tol_interrupt_simulation:
            self.tol_interrupt_simulation = tol_interrupt_simulation
        else:
            self.tol_interrupt_simulation = tol_move_to_core

        if tol_move_edge_rxn_to_surface_interrupt:
            self.tol_move_edge_rxn_to_surface_interrupt = tol_move_edge_rxn_to_surface_interrupt
        else:
            self.tol_move_edge_rxn_to_surface_interrupt = tol_move_edge_rxn_to_surface

        if tol_move_edge_rxn_to_core_interrupt:
            self.tol_move_edge_rxn_to_core_interrupt = tol_move_edge_rxn_to_core_interrupt
        else:
            self.tol_move_edge_rxn_to_core_interrupt = tol_move_edge_rxn_to_core

        if max_num_species:
            self.max_num_species = max_num_species
        else:
            self.max_num_species = np.inf

        if max_num_objects_per_iter <= 0:
            self.max_num_objects_per_iter = np.inf
        else:
            self.max_num_objects_per_iter = max_num_objects_per_iter


class SimulatorSettings(object):
    """
    class for holding the parameters affecting the behavior of the solver
    """

    def __init__(self, atol=1e-16, rtol=1e-8, sens_atol=1e-6, sens_rtol=1e-4):
        self.atol = atol
        self.rtol = rtol
        self.sens_atol = sens_atol
        self.sens_rtol = sens_rtol
