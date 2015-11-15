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

    cdef public numpy.ndarray coreSpeciesConcentrations
    cdef public numpy.ndarray coreSpeciesRates
    cdef public numpy.ndarray coreReactionRates
    cdef public numpy.ndarray edgeSpeciesRates
    cdef public numpy.ndarray edgeReactionRates
    cdef public numpy.ndarray networkLeakRates

    cdef public numpy.ndarray maxCoreSpeciesRates
    cdef public numpy.ndarray maxEdgeSpeciesRates
    cdef public numpy.ndarray maxNetworkLeakRates
    cdef public numpy.ndarray maxEdgeSpeciesRateRatios
    cdef public numpy.ndarray maxNetworkLeakRateRatios
    cdef public numpy.ndarray sensitivityCoefficients
    
    cdef public list snapshots

    cdef public list termination

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list pdepNetworks=?, atol=?, rtol=?, sensitivity=?, sens_atol=?, sens_rtol=?)
    
    cpdef simulate(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions,
        double toleranceKeepInEdge, double toleranceMoveToCore, double toleranceInterruptSimulation,
        list pdepNetworks=?, worksheet=?, absoluteTolerance=?, relativeTolerance=?, sensitivity=?, sensitivityAbsoluteTolerance=?, sensitivityRelativeTolerance=?, sensWorksheet=?)

    cpdef logRates(self, double charRate, object species, double speciesRate, object network, double networkRate)

    cpdef logConversions(self, speciesIndex, y0)
