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
from scipy import sparse
from pydas import DASSL

class SimpleReactor(DASSL):
    """
    A reaction system consisting of a homogeneous, isothermal, isobaric batch
    reactor. These assumptions allow for a number of optimizations that enable
    this solver to complete very rapidly, even for large kinetic models.
    """

    def __init__(self, T, P, C0):
        self.T = T
        self.P = P
        self.C0 = C0

        # The reaction and species rates at the current time (in mol/m^3*s)
        self.reactionRates = None
        self.speciesRates = None

        # These are helper variables used within the solver
        self.reactantIndices = None
        self.productIndices = None
        self.forwardRateCoefficients = None
        self.reverseRateCoefficients = None


    def initialize(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """
        numCoreSpecies = len(coreSpecies)
        numCoreReactions = len(coreReactions)
        numEdgeSpecies = len(edgeSpecies)
        numEdgeReactions = len(edgeReactions)

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
        reactantIndices = -numpy.ones((numCoreReactions + numEdgeReactions, 3), numpy.int32 )
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
        
        self.reactantIndices = reactantIndices
        self.productIndices = productIndices
        self.forwardRateCoefficients = forwardRateCoefficients
        self.reverseRateCoefficients = reverseRateCoefficients

        self.reactionRates = numpy.zeros((numCoreReactions + numEdgeReactions), numpy.float64)
        self.speciesRates = numpy.zeros((numCoreSpecies + numEdgeSpecies), numpy.float64)

        # Set initial conditions
        t0 = 0.0
        y0 = numpy.zeros((numCoreSpecies), numpy.float64)
        for spec, C in self.C0.iteritems():
            y0[speciesIndex[spec]] = C
        
        # Initialize the model
        dydt0 = - self.residual(t0, y0, numpy.zeros((numCoreSpecies), numpy.float64))[0]
        DASSL.initialize(self, t0, y0, dydt0)

    def residual(self, t, y, dydt):
        """
        Return the residual function for the governing DAE system for the
        simple reaction system.
        """
        N = y.shape[0]
        res = numpy.zeros(N, numpy.float64)

        ir = self.reactantIndices
        ip = self.productIndices
        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients

        reactionRates = numpy.zeros_like(self.reactionRates)
        speciesRates = numpy.zeros_like(self.speciesRates)

        for j in range(ir.shape[0]):
            kf = self.forwardRateCoefficients[j]
            if ir[j,0] >= N or ir[j,1] >= N or ir[j,2] >= N:
                reactionRate = 0.0
            elif ir[j,1] == -1: # only one reactant
                reactionRate = kf * y[ir[j,0]]
            elif ir[j,2] == -1: # only two reactants
                reactionRate = kf * y[ir[j,0]] * y[ir[j,1]]
            else: # three reactants!! (really?)
                reactionRate = kf * y[ir[j,0]] * y[ir[j,1]] * y[ir[j,2]]
            kr = self.reverseRateCoefficients[j]
            if ip[j,0] >= N or ip[j,1] >= N or ip[j,2] >= N:
                pass
            elif ip[j,1] == -1: # only one reactant
                reactionRate -= kr * y[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                reactionRate -= kr * y[ip[j,0]] * y[ip[j,1]]
            else: # three reactants!! (really?)
                reactionRate -= kr * y[ip[j,0]] * y[ip[j,1]] * y[ip[j,2]]

            # Set the total reaction rate
            reactionRates[j] = reactionRate

            # Add/substract the total reaction rate from each species rate
            speciesRates[ir[j,0]] -= reactionRate
            second = ir[j,1]
            if second != -1:
                speciesRates[second] -= reactionRate
                third = ir[j,2]
                if third != -1:
                    speciesRates[third] -= reactionRate
            speciesRates[ip[j,0]] += reactionRate
            second = ip[j,1]
            if second != -1:
                speciesRates[second] += reactionRate
                third = ip[j,2]
                if third != -1:
                    speciesRates[third] += reactionRate

        self.reactionRates = reactionRates
        self.speciesRates = speciesRates

        res = speciesRates[0:N] - dydt
        return res, 0
