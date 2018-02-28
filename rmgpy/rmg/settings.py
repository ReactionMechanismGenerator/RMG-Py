#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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
    `atol`                                          The absolute tolerance used in the ODE/DAE solver
    `rtol`                                          The relative tolerance used in the ODE/DAE solver
    `sens_atol`                                     The absolute tolerance used in the ODE/DAE solver for the sensitivities
    `sens_rtol`                                     The relative tolerance used in the ODE/DAE solver for the sensitivities
    `toleranceKeepInEdge`                           The relative species flux below which species are discarded from the edge
    `toleranceMoveToCore`                           The relative species flux above which species are moved from the edge to the core
    `toleranceInterruptSimulation`                  The relative species flux above which the simulation will halt
    `toleranceMoveEdgeReactionToCore`               The dynamics number above which reactions are moved from the edge to the core
    `toleranceMoveEdgeReactionToCoreInterrupt`      The edge dynamics number above which the simulation will halt
    `toleranceMoveEdgeReactionToSurface`            The dynamics number above which reactions are moved from the edge to the surface
    `toleranceMoveEdgeReactionToSurfaceInterrupt`   The edge dynamics number above which simulation will halt
    `toleranceMoveSurfaceReactionToCore`            The surface dynamics number above which reactions are moved from the surface to the bulk core
    `toleranceMoveSurfaceSpeciesToCore`             The relative species flux above which species are moved from the surface to the bulk core
    `maximumEdgeSpecies`                            The maximum number of edge species allowed at any time
    `minCoreSizeForPrune`                           Minimum number of core species before pruning is allowed
    `minSpeciesExistIterationsForPrune`             Minimum number of iterations a species must exist before it can be pruned
    `filterReactions`                               Specify whether to filter reactions during model enlarging step
    `ignoreOverallFluxCriterion`                    flag indicating that the ordinary flux criterion should be ignored except for pdep purposes
    `maxNumSpecies`                                 Number of core species at which a stage/job will terminate
    `maxNumObjPerIter`                              Maximum number of objects that can be sent for enlargement from a single simulation
==================================================================================================================================================
"""
import numpy
from rmgpy.quantity import Quantity

class ModelSettings(object):
    """
    class for holding the parameters affecting an RMG run
    """
    def __init__(self,toleranceMoveToCore=None, toleranceMoveEdgeReactionToCore=numpy.inf,toleranceKeepInEdge=0.0, toleranceInterruptSimulation=1.0, 
          toleranceMoveEdgeReactionToSurface=numpy.inf, toleranceMoveSurfaceSpeciesToCore=numpy.inf, toleranceMoveSurfaceReactionToCore=numpy.inf,
          toleranceMoveEdgeReactionToSurfaceInterrupt=None,toleranceMoveEdgeReactionToCoreInterrupt=None, maximumEdgeSpecies=1000000, minCoreSizeForPrune=50, 
          minSpeciesExistIterationsForPrune=2, filterReactions=False, ignoreOverallFluxCriterion=False, maxNumSpecies=None, maxNumObjsPerIter=1,
          terminateAtMaxObjects=False,toleranceThermoKeepSpeciesInEdge=numpy.inf,dynamicsTimeScale = Quantity((0.0,'sec')),
          branching=None):

        
        self.fluxToleranceKeepInEdge = toleranceKeepInEdge
        self.fluxToleranceMoveToCore = toleranceMoveToCore
        self.toleranceMoveEdgeReactionToCore = toleranceMoveEdgeReactionToCore
        self.fluxToleranceInterrupt = toleranceInterruptSimulation
        self.maximumEdgeSpecies = maximumEdgeSpecies
        self.minCoreSizeForPrune = minCoreSizeForPrune
        self.minSpeciesExistIterationsForPrune = minSpeciesExistIterationsForPrune
        self.filterReactions = filterReactions
        self.ignoreOverallFluxCriterion=ignoreOverallFluxCriterion
        self.toleranceMoveEdgeReactionToSurface = toleranceMoveEdgeReactionToSurface
        self.toleranceMoveSurfaceSpeciesToCore = toleranceMoveSurfaceSpeciesToCore
        self.toleranceMoveSurfaceReactionToCore = toleranceMoveSurfaceReactionToCore
        self.toleranceThermoKeepSpeciesInEdge = toleranceThermoKeepSpeciesInEdge
        self.terminateAtMaxObjects = terminateAtMaxObjects
        self.dynamicsTimeScale = dynamicsTimeScale.value_si
        
        self.branching = branching
        
        if toleranceInterruptSimulation:
            self.fluxToleranceInterrupt = toleranceInterruptSimulation
        else:
            self.fluxToleranceInterrupt = toleranceMoveToCore
            
        if toleranceMoveEdgeReactionToSurfaceInterrupt:
            self.toleranceMoveEdgeReactionToSurfaceInterrupt = toleranceMoveEdgeReactionToSurfaceInterrupt
        else:
            self.toleranceMoveEdgeReactionToSurfaceInterrupt = toleranceMoveEdgeReactionToSurface
        
        if toleranceMoveEdgeReactionToCoreInterrupt:
            self.toleranceMoveEdgeReactionToCoreInterrupt = toleranceMoveEdgeReactionToCoreInterrupt
        else:
            self.toleranceMoveEdgeReactionToCoreInterrupt = toleranceMoveEdgeReactionToCore
            
        if maxNumSpecies:
            self.maxNumSpecies = maxNumSpecies
        else:
            self.maxNumSpecies = numpy.inf
        
        if maxNumObjsPerIter <= 0:
            self.maxNumObjsPerIter = numpy.inf
        else:
            self.maxNumObjsPerIter = maxNumObjsPerIter
            
class SimulatorSettings(object):
    """
    class for holding the parameters affecting the behavior of the solver
    """
    def __init__(self,atol=1e-16, rtol=1e-8, sens_atol=1e-6, sens_rtol=1e-4):
        self.atol = atol
        self.rtol = rtol
        self.sens_atol = sens_atol
        self.sens_rtol = sens_rtol