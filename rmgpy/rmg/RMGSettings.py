#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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

"""
This module contains settings classes for manipulation of RMG run parameters
"""
import numpy

class ModelSettings:
    """
    class for holding the parameters affecting an RMG run
    """
    def __init__(self,toleranceMoveToCore=None, toleranceMoveEdgeReactionToCore=numpy.inf,toleranceKeepInEdge=0.0, toleranceInterruptSimulation=1.0, 
          toleranceMoveEdgeReactionToSurface=numpy.inf, toleranceMoveSurfaceSpeciesToCore=numpy.inf, toleranceMoveSurfaceReactionToCore=numpy.inf,
          toleranceMoveEdgeReactionToSurfaceInterrupt=None,toleranceMoveEdgeReactionToCoreInterrupt=None, maximumEdgeSpecies=1000000, minCoreSizeForPrune=50, 
          minSpeciesExistIterationsForPrune=2, filterReactions=False, ignoreOverallFluxCriterion=False, maxNumSpecies=None, maxNumObjsPerIter=1):
        
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
            
class SimulatorSettings:
    """
    class for holding the parameters affecting the behavior of the solver
    """
    def __init__(self,atol=1e-16, rtol=1e-8, sens_atol=1e-6, sens_rtol=1e-4):
        self.absoluteTolerance = atol
        self.relativeTolerance = rtol
        self.sensitivityAbsoluteTolerance = sens_atol
        self.sensitivityRelativeTolerance = sens_rtol