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

cimport numpy as np
from cpython cimport bool
include "settings.pxi"
if DASPK == 1:
    from pydas.daspk cimport DASPK as DASx
else:
    from pydas.dassl cimport DASSL as DASx

################################################################################

cdef class ReactionSystem(DASx):

    # reactor state variables:
    cdef public float t0
    cdef public np.ndarray y0
    cdef public np.ndarray dydt0

    #  variables that determine the dimensions of arrays and matrices:
    cdef public Py_ssize_t num_core_species
    cdef public Py_ssize_t num_core_reactions
    cdef public Py_ssize_t num_edge_species
    cdef public Py_ssize_t num_edge_reactions
    cdef public Py_ssize_t num_pdep_networks
    cdef public Py_ssize_t neq

    # variables that store stoichiometry data
    cdef public dict species_index
    cdef public dict reaction_index
    cdef public np.ndarray reactant_indices
    cdef public np.ndarray product_indices
    cdef public np.ndarray network_indices

    # matrices that cache kinetic and rate data
    cdef public np.ndarray kf  # forward rate coefficients
    cdef public np.ndarray kb  # reverse rate coefficients
    cdef public np.ndarray Keq  # equilibrium constants
    cdef public np.ndarray network_leak_coefficients
    cdef public np.ndarray jacobian_matrix

    cdef public np.ndarray core_species_concentrations
    
    #surface information
    cdef public np.ndarray surface_species_indices
    cdef public np.ndarray surface_reaction_indices
    cdef public np.ndarray valid_layering_indices
    
    # The reaction and species rates at the current time (in mol/m^3*s)
    cdef public np.ndarray core_species_rates
    cdef public np.ndarray core_reaction_rates
    cdef public np.ndarray core_species_production_rates
    cdef public np.ndarray core_species_consumption_rates
    cdef public np.ndarray edge_species_rates
    cdef public np.ndarray edge_reaction_rates

    cdef public np.ndarray network_leak_rates    

    # variables that cache maximum rate (ratio) data
    cdef public np.ndarray max_edge_species_rate_ratios
    cdef public np.ndarray max_network_leak_rate_ratios
    
    #for managing prunable edge species
    cdef public list prunable_species
    cdef public list prunable_networks
    cdef public np.ndarray prunable_species_indices
    cdef public np.ndarray prunable_network_indices
    
    # sensitivity variables
    # cdef public int sensmethod
    cdef public np.ndarray sensitivity_coefficients
    cdef public list sensitive_species
    cdef public double sensitivity_threshold
    # cdef public np.ndarray senpar

    # tolerance settings
    cdef public np.ndarray atol_array
    cdef public np.ndarray rtol_array
    
    cdef public list snapshots

    cdef public list termination
    
    # Trimolecular reactants flag
    cdef public bint trimolecular

    # reaction threshold settings
    cdef public np.ndarray unimolecular_threshold
    cdef public np.ndarray bimolecular_threshold
    cdef public np.ndarray trimolecular_threshold

    cdef public int retry
    cdef public int max_retries

    # methods
    cpdef initialize_model(self, list core_species, list core_reactions, list edge_species, list edge_reactions,
        list surface_species=?, list surface_reactions=?, list pdep_networks=?, atol=?, rtol=?,
        sensitivity=?, sens_atol=?, sens_rtol=?, filter_reactions=?, dict conditions=?)

    cpdef simulate(self, list core_species, list core_reactions, list edge_species, 
        list edge_reactions,list surface_species, list surface_reactions,
        list pdep_networks=?, bool prune=?, bool sensitivity=?, list sens_worksheet=?, object model_settings=?,
        object simulator_settings=?, dict conditions=?)

    cpdef log_rates(self, double char_rate, object species, double species_rate, double max_dif_ln_accum_num, object network, double network_rate)
     
    cpdef log_conversions(self, species_index, y0)
    
    cpdef get_layering_indices(self)
    
    cpdef initialize_surface(self,list core_species,list core_reactions,list surface_species,list surface_reactions)
    
    cpdef add_reactions_to_surface(self,list new_surface_reactions,list new_surface_reaction_inds,list surface_species,list surface_reactions,list edge_species)
