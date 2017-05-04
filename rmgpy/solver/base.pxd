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

cimport numpy

include "settings.pxi"
if DASPK == 1:
    from pydas.daspk cimport DASPK as DASx
else:
    from pydas.dassl cimport DASSL as DASx
    
################################################################################

cdef class ReactionSystem(DASx):

    # reactor state variables:
    cdef public float t0
    cdef public numpy.ndarray y0
    cdef public numpy.ndarray dydt0

    #  variables that determine the dimensions of arrays and matrices:
    cdef public int numCoreSpecies
    cdef public int numCoreReactions
    cdef public int numEdgeSpecies
    cdef public int numEdgeReactions
    cdef public int numPdepNetworks
    cdef public int neq

    # variables that store stoichiometry data
    cdef public dict speciesIndex
    cdef public dict reactionIndex
    cdef public numpy.ndarray reactantIndices
    cdef public numpy.ndarray productIndices
    cdef public numpy.ndarray networkIndices

    # matrices that cache kinetic and rate data
    cdef public numpy.ndarray kf # forward rate coefficients
    cdef public numpy.ndarray kb # reverse rate coefficients
    cdef public numpy.ndarray Keq # equilibrium constants
    cdef public numpy.ndarray networkLeakCoefficients
    cdef public numpy.ndarray jacobianMatrix

    cdef public numpy.ndarray coreSpeciesConcentrations
    
    #surface information
    cdef public numpy.ndarray surfaceSpeciesIndices
    cdef public numpy.ndarray surfaceReactionIndices
    
    # The reaction and species rates at the current time (in mol/m^3*s)
    cdef public numpy.ndarray coreSpeciesRates
    cdef public numpy.ndarray coreReactionRates
    cdef public numpy.ndarray coreSpeciesProductionRates
    cdef public numpy.ndarray coreSpeciesConsumptionRates
    cdef public numpy.ndarray edgeSpeciesRates
    cdef public numpy.ndarray edgeReactionRates

    cdef public numpy.ndarray networkLeakRates    

    # variables that cache maximum rate (ratio) data
    cdef public numpy.ndarray maxCoreSpeciesRates
    cdef public numpy.ndarray maxEdgeSpeciesRates
    cdef public numpy.ndarray maxNetworkLeakRates
    cdef public numpy.ndarray maxEdgeSpeciesRateRatios
    cdef public numpy.ndarray maxNetworkLeakRateRatios

    # sensitivity variables
    # cdef public int sensmethod
    cdef public numpy.ndarray sensitivityCoefficients
    cdef public list sensitiveSpecies
    cdef public double sensitivityThreshold
    # cdef public numpy.ndarray senpar

    # tolerance settings
    cdef public numpy.ndarray atol_array
    cdef public numpy.ndarray rtol_array
    
    cdef public list snapshots

    cdef public list termination
    
    # reaction threshold settings
    cdef public numpy.ndarray unimolecularThreshold
    cdef public numpy.ndarray bimolecularThreshold

    # methods
    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list surfaceSpecies=?,
        list surfaceReactions=?, list pdepNetworks=?, atol=?, rtol=?, sensitivity=?, sens_atol=?, sens_rtol=?, filterReactions=?)

    cpdef simulate(self, list coreSpecies, list coreReactions, list edgeSpecies, list surfaceSpecies, list edgeReactions,list surfaceReactions,
        double toleranceKeepInEdge, double toleranceMoveToCore,  double toleranceInterruptSimulation,double toleranceMoveEdgeReactionToCore=?,double toleranceMoveEdgeReactionToCoreInterrupt=?,
        double toleranceMoveEdgeReactionToSurface=?, double toleranceMoveSurfaceSpeciesToCore=?, double toleranceMoveSurfaceReactionToCore=?,
        double toleranceMoveEdgeReactionToSurfaceInterrupt=?, 
        list pdepNetworks=?,ignoreOverallFluxCriterion=?, absoluteTolerance=?, relativeTolerance=?, sensitivity=?, 
        sensitivityAbsoluteTolerance=?, sensitivityRelativeTolerance=?, sensWorksheet=?,
        filterReactions=?)

    cpdef logRates(self, double charRate, object species, double speciesRate, double maxDifLnAccumNum, object network, double networkRate)
     
    cpdef logConversions(self, speciesIndex, y0)
    
    cpdef maxIndUnderSurfaceLayeringConstraint(self,numpy.ndarray[numpy.float64_t,ndim=1] arr,numpy.ndarray[numpy.int_t,ndim=1] surfSpeciesIndices)
    
    cpdef initialize_surface(self,list coreSpecies,list coreReactions,list surfaceSpecies,list surfaceReactions)
    
    cpdef addReactionsToSurface(self,list newSurfaceReactions,list newSurfaceReactionInds,list surfaceSpecies,list surfaceReactions,list edgeSpecies)