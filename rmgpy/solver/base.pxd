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
    cdef public int numCoreSpecies
    cdef public int numCoreReactions
    cdef public int numEdgeSpecies
    cdef public int numSurfaceSpecies
    cdef public int numSurfaceReactions
    cdef public int numEdgeReactions
    cdef public int numPdepNetworks
    cdef public int neq

    # variables that store stoichiometry data
    cdef public dict speciesIndex
    cdef public dict reactionIndex
    cdef public np.ndarray reactantIndices
    cdef public np.ndarray productIndices
    cdef public np.ndarray networkIndices

    # matrices that cache kinetic and rate data
    cdef public np.ndarray kf  # forward rate coefficients
    cdef public np.ndarray kb  # reverse rate coefficients
    cdef public np.ndarray Keq  # equilibrium constants
    cdef public np.ndarray networkLeakCoefficients
    cdef public np.ndarray jacobianMatrix

    cdef public np.ndarray coreSpeciesConcentrations
    
    #surface information
    cdef public np.ndarray surfaceSpeciesIndices
    cdef public np.ndarray surfaceReactionIndices
    cdef public np.ndarray validLayeringIndices
    
    # The reaction and species rates at the current time (in mol/m^3*s)
    cdef public np.ndarray coreSpeciesRates
    cdef public np.ndarray coreReactionRates
    cdef public np.ndarray coreSpeciesProductionRates
    cdef public np.ndarray coreSpeciesConsumptionRates
    cdef public np.ndarray edgeSpeciesRates
    cdef public np.ndarray edgeReactionRates

    cdef public np.ndarray networkLeakRates    

    # variables that cache maximum rate (ratio) data
    cdef public np.ndarray maxEdgeSpeciesRateRatios
    cdef public np.ndarray maxNetworkLeakRateRatios
    
    #for managing prunable edge species
    cdef public list prunableSpecies
    cdef public list prunableNetworks
    cdef public np.ndarray prunableSpeciesIndices
    cdef public np.ndarray prunableNetworkIndices
    
    # sensitivity variables
    # cdef public int sensmethod
    cdef public np.ndarray sensitivityCoefficients
    cdef public list sensitiveSpecies
    cdef public double sensitivityThreshold
    # cdef public np.ndarray senpar

    # tolerance settings
    cdef public np.ndarray atol_array
    cdef public np.ndarray rtol_array
    
    cdef public list snapshots

    cdef public list termination
    
    # Trimolecular reactants flag
    cdef public bint trimolecular

    # reaction threshold settings
    cdef public np.ndarray unimolecularThreshold
    cdef public np.ndarray bimolecularThreshold
    cdef public np.ndarray trimolecularThreshold

    # methods
    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions,
        list surfaceSpecies=?, list surfaceReactions=?, list pdepNetworks=?, atol=?, rtol=?,
        sensitivity=?, sens_atol=?, sens_rtol=?, filterReactions=?, dict conditions=?)

    cpdef simulate(self, list coreSpecies, list coreReactions, list edgeSpecies, 
        list edgeReactions,list surfaceSpecies, list surfaceReactions,
        list pdepNetworks=?, bool prune=?, bool sensitivity=?, list sensWorksheet=?, object modelSettings=?,
        object simulatorSettings=?, dict conditions=?)

    cpdef logRates(self, double charRate, object species, double speciesRate, double maxDifLnAccumNum, object network, double networkRate)
     
    cpdef logConversions(self, speciesIndex, y0)
    
    cpdef getLayeringIndices(self)
    
    cpdef initialize_surface(self,list coreSpecies,list coreReactions,list surfaceSpecies,list surfaceReactions)
    
    cpdef addReactionsToSurface(self,list newSurfaceReactions,list newSurfaceReactionInds,list surfaceSpecies,list surfaceReactions,list edgeSpecies)