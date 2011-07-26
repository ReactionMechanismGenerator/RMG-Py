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
Contains the :class:`SimpleReactor` class, providing a reaction system
consisting of a homogeneous, isothermal, isobaric batch reactor.
"""

import numpy
cimport numpy
from pydas cimport DASSL
from base cimport ReactionSystem
cimport cython

from rmgpy.quantity import constants

cdef class SimpleReactor(ReactionSystem):
    """
    A reaction system consisting of a homogeneous, isothermal, isobaric batch
    reactor. These assumptions allow for a number of optimizations that enable
    this solver to complete very rapidly, even for large kinetic models.
    """

    cdef public double T
    cdef public double P
    cdef public dict initialMoleFractions

    cdef numpy.ndarray reactantIndices
    cdef numpy.ndarray productIndices
    cdef numpy.ndarray networkIndices
    cdef numpy.ndarray forwardRateCoefficients
    cdef numpy.ndarray reverseRateCoefficients
    cdef numpy.ndarray networkLeakCoefficients

    def __init__(self, T, P, initialMoleFractions):
        ReactionSystem.__init__(self)
        self.T = T
        self.P = P
        self.initialMoleFractions = initialMoleFractions

        # These are helper variables used within the solver
        self.reactantIndices = None
        self.productIndices = None
        self.networkIndices = None
        self.forwardRateCoefficients = None
        self.reverseRateCoefficients = None

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list pdepNetworks=None, atol=1e-16, rtol=1e-8):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions, pdepNetworks, atol, rtol)

        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks
        cdef int i, j, l, index
        cdef dict speciesIndex, reactionIndex
        cdef numpy.ndarray[numpy.int_t, ndim=2] reactantIndices, productIndices, networkIndices
        cdef numpy.ndarray[numpy.float64_t, ndim=1] forwardRateCoefficients, reverseRateCoefficients, networkLeakCoefficients
        
        pdepNetworks = pdepNetworks or []

        numCoreSpecies = len(coreSpecies)
        numCoreReactions = len(coreReactions)
        numEdgeSpecies = len(edgeSpecies)
        numEdgeReactions = len(edgeReactions)
        numPdepNetworks = len(pdepNetworks)

        # Assign an index to each species (core first, then edge)
        speciesIndex = {}
        for index, spec in enumerate(coreSpecies):
            speciesIndex[spec] = index
        for index, spec in enumerate(edgeSpecies):
            speciesIndex[spec] = index + numCoreSpecies
        # Assign an index to each reaction (core first, then edge)
        reactionIndex = {}
        for index, rxn in enumerate(coreReactions):
            reactionIndex[rxn] = index
        for index, rxn in enumerate(edgeReactions):
            reactionIndex[rxn] = index + numCoreReactions

        # Generate reactant and product indices
        # Generate forward and reverse rate coefficients k(T,P)
        reactantIndices = -numpy.ones((numCoreReactions + numEdgeReactions, 3), numpy.int )
        productIndices = -numpy.ones_like(reactantIndices)
        forwardRateCoefficients = numpy.zeros((numCoreReactions + numEdgeReactions), numpy.float64)
        reverseRateCoefficients = numpy.zeros_like(forwardRateCoefficients)
        for rxnList in [coreReactions, edgeReactions]:
            for rxn in rxnList:
                j = reactionIndex[rxn]
                forwardRateCoefficients[j] = rxn.getRateCoefficient(self.T, self.P)
                reverseRateCoefficients[j] = forwardRateCoefficients[j] / rxn.getEquilibriumConstant(self.T)
                for l, spec in enumerate(rxn.reactants):
                    i = speciesIndex[spec]
                    reactantIndices[j,l] = i
                for l, spec in enumerate(rxn.products):
                    i = speciesIndex[spec]
                    productIndices[j,l] = i

        networkIndices = -numpy.ones((numPdepNetworks, 3), numpy.int )
        networkLeakCoefficients = numpy.zeros((numPdepNetworks), numpy.float64)
        for j, network in enumerate(pdepNetworks):
            networkLeakCoefficients[j] = network.getLeakCoefficient(self.T, self.P)
            for l, spec in enumerate(network.source):
                i = speciesIndex[spec]
                networkIndices[j,l] = i

        self.reactantIndices = reactantIndices
        self.productIndices = productIndices
        self.forwardRateCoefficients = forwardRateCoefficients
        self.reverseRateCoefficients = reverseRateCoefficients
        self.networkIndices = networkIndices
        self.networkLeakCoefficients = networkLeakCoefficients
        
        # Set initial conditions
        t0 = 0.0
        y0 = numpy.zeros((numCoreSpecies), numpy.float64)
        for spec, moleFrac in self.initialMoleFractions.iteritems():
            y0[speciesIndex[spec]] = moleFrac * (self.P / constants.R / self.T)
            self.coreSpeciesConcentrations[speciesIndex[spec]] = y0[speciesIndex[spec]]
        
        # Initialize the model
        dydt0 = - self.residual(t0, y0, numpy.zeros((numCoreSpecies), numpy.float64))[0]
        DASSL.initialize(self, t0, y0, dydt0, atol, rtol)

    cpdef writeWorksheetHeader(self, worksheet):
        """
        Write some descriptive information about the reaction system to the
        first two rows of the given `worksheet`.
        """
        import xlwt
        style0 = xlwt.easyxf('font: bold on')
        worksheet.write(0, 0, 'Simple Reactor', style0)
        worksheet.write(1, 0, 'T = {0:g} K, P = {1:g} bar'.format(self.T, self.P/1e5))

    @cython.boundscheck(False)
    def residual(self, double t, numpy.ndarray[numpy.float64_t, ndim=1] y, numpy.ndarray[numpy.float64_t, ndim=1] dydt):

        """
        Return the residual function for the governing DAE system for the
        simple reaction system.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip, inet
        cdef numpy.ndarray[numpy.float64_t, ndim=1] res, kf, kr, knet
        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks
        cdef int j, first, second, third
        cdef double k, reactionRate
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesRates, coreReactionRates, edgeSpeciesRates, edgeReactionRates, networkLeakRates

        res = numpy.zeros(y.shape[0], numpy.float64)

        ir = self.reactantIndices
        ip = self.productIndices
        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients
        inet = self.networkIndices
        knet = self.networkLeakCoefficients

        numCoreSpecies = len(self.coreSpeciesRates)
        numCoreReactions = len(self.coreReactionRates)
        numEdgeSpecies = len(self.edgeSpeciesRates)
        numEdgeReactions = len(self.edgeReactionRates)
        numPdepNetworks = len(self.networkLeakRates)

        coreSpeciesConcentrations = numpy.zeros_like(self.coreSpeciesConcentrations)
        coreSpeciesRates = numpy.zeros_like(self.coreSpeciesRates)
        coreReactionRates = numpy.zeros_like(self.coreReactionRates)
        edgeSpeciesRates = numpy.zeros_like(self.edgeSpeciesRates)
        edgeReactionRates = numpy.zeros_like(self.edgeReactionRates)
        networkLeakRates = numpy.zeros_like(self.networkLeakRates)

        for j in range(y.shape[0]):
            coreSpeciesConcentrations[j] = y[j]
        
        for j in range(ir.shape[0]):
            k = kf[j]
            if ir[j,0] >= numCoreSpecies or ir[j,1] >= numCoreSpecies or ir[j,2] >= numCoreSpecies:
                reactionRate = 0.0
            elif ir[j,1] == -1: # only one reactant
                reactionRate = k * y[ir[j,0]]
            elif ir[j,2] == -1: # only two reactants
                reactionRate = k * y[ir[j,0]] * y[ir[j,1]]
            else: # three reactants!! (really?)
                reactionRate = k * y[ir[j,0]] * y[ir[j,1]] * y[ir[j,2]]
            k = kr[j]
            if ip[j,0] >= numCoreSpecies or ip[j,1] >= numCoreSpecies or ip[j,2] >= numCoreSpecies:
                pass
            elif ip[j,1] == -1: # only one reactant
                reactionRate -= k * y[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                reactionRate -= k * y[ip[j,0]] * y[ip[j,1]]
            else: # three reactants!! (really?)
                reactionRate -= k * y[ip[j,0]] * y[ip[j,1]] * y[ip[j,2]]

            # Set the reaction and species rates
            if j < numCoreReactions:
                # The reaction is a core reaction
                coreReactionRates[j] = reactionRate

                # Add/substract the total reaction rate from each species rate
                # Since it's a core reaction we know that all of its reactants
                # and products are core species
                first = ir[j,0]
                coreSpeciesRates[first] -= reactionRate
                second = ir[j,1]
                if second != -1:
                    coreSpeciesRates[second] -= reactionRate
                    third = ir[j,2]
                    if third != -1:
                        coreSpeciesRates[third] -= reactionRate
                first = ip[j,0]
                coreSpeciesRates[first] += reactionRate
                second = ip[j,1]
                if second != -1:
                    coreSpeciesRates[second] += reactionRate
                    third = ip[j,2]
                    if third != -1:
                        coreSpeciesRates[third] += reactionRate

            else:
                # The reaction is an edge reaction
                edgeReactionRates[j-numCoreReactions] = reactionRate

                # Add/substract the total reaction rate from each species rate
                # Since it's an edge reaction its reactants and products could
                # be either core or edge species
                # We're only interested in the edge species
                first = ir[j,0]
                if first >= numCoreSpecies: edgeSpeciesRates[first-numCoreSpecies] -= reactionRate
                second = ir[j,1]
                if second != -1:
                    if second >= numCoreSpecies: edgeSpeciesRates[second-numCoreSpecies] -= reactionRate
                    third = ir[j,2]
                    if third != -1:
                        if third >= numCoreSpecies: edgeSpeciesRates[third-numCoreSpecies] -= reactionRate
                first = ip[j,0]
                if first >= numCoreSpecies: edgeSpeciesRates[first-numCoreSpecies] += reactionRate
                second = ip[j,1]
                if second != -1:
                    if second >= numCoreSpecies: edgeSpeciesRates[second-numCoreSpecies] += reactionRate
                    third = ip[j,2]
                    if third != -1:
                        if third >= numCoreSpecies: edgeSpeciesRates[third-numCoreSpecies] += reactionRate

        for j in range(inet.shape[0]):
            k = knet[j]
            if inet[j,1] == -1: # only one reactant
                reactionRate = k * y[inet[j,0]]
            elif inet[j,2] == -1: # only two reactants
                reactionRate = k * y[inet[j,0]] * y[inet[j,1]]
            else: # three reactants!! (really?)
                reactionRate = k * y[inet[j,0]] * y[inet[j,1]] * y[inet[j,2]]
            networkLeakRates[j] = reactionRate

        self.coreSpeciesConcentrations = coreSpeciesConcentrations
        self.coreSpeciesRates = coreSpeciesRates
        self.coreReactionRates = coreReactionRates
        self.edgeSpeciesRates = edgeSpeciesRates
        self.edgeReactionRates = edgeReactionRates
        self.networkLeakRates = networkLeakRates

        res = coreSpeciesRates - dydt
        return res, 0
