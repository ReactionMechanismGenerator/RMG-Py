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

include "settings.pxi"
if DASPK == 1:
    from pydas.daspk cimport DASPK as DASx
else:
    from pydas.dassl cimport DASSL as DASx
    
from base cimport ReactionSystem
cimport cython

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

cdef class LiquidReactor(ReactionSystem):
    """
    A reaction system consisting of a homogeneous, isothermal, constant volume batch
    reactor. These assumptions allow for a number of optimizations that enable
    this solver to complete very rapidly, even for large kinetic models.
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity P    
    cdef public double V
    cdef public bint constantVolume
    cdef public dict initialConcentrations
    cdef public list sensitiveSpecies
    cdef public double sensitivityThreshold

    cdef public numpy.ndarray reactantIndices
    cdef public numpy.ndarray productIndices
    cdef public numpy.ndarray networkIndices
    cdef public numpy.ndarray forwardRateCoefficients
    cdef public numpy.ndarray reverseRateCoefficients
    cdef public numpy.ndarray equilibriumConstants
    cdef public numpy.ndarray networkLeakCoefficients
    cdef public numpy.ndarray jacobianMatrix

    def __init__(self, T, initialConcentrations, termination, sensitiveSpecies=None, sensitivityThreshold=1e-3):
        ReactionSystem.__init__(self, termination)
        self.T = Quantity(T)
        self.P = Quantity(100000.,'kPa') # Arbitrary high pressure (1000 Bar) to get reactions in the high-pressure limit!
        self.initialConcentrations = initialConcentrations # should be passed in SI
        self.V = 0 # will be set from initialConcentrations in initializeModel
        self.constantVolume = True
      
        self.sensitiveSpecies = sensitiveSpecies
        self.sensitivityThreshold = sensitivityThreshold
        
        # These are helper variables used within the solver
        self.reactantIndices = None
        self.productIndices = None
        self.networkIndices = None
        self.forwardRateCoefficients = None
        self.reverseRateCoefficients = None
        self.equilibriumConstants = None
        self.networkLeakCoefficients = None
        self.jacobianMatrix = None
        
    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initialConcentrations dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initialConcentrations = {}
        for label, moleFrac in self.initialConcentrations.iteritems():
            initialConcentrations[speciesDict[label]] = moleFrac
        self.initialConcentrations = initialConcentrations

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list pdepNetworks=None, atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions, pdepNetworks, atol, rtol, sensitivity, sens_atol, sens_rtol)

        cdef int i, j, l, index
        cdef double V
        cdef numpy.ndarray[numpy.float64_t, ndim=1] forwardRateCoefficients, reverseRateCoefficients, equilibriumConstants

        # Generate reactant and product indices
        forwardRateCoefficients = numpy.zeros((numCoreReactions + numEdgeReactions), numpy.float64)
        reverseRateCoefficients = numpy.zeros_like(forwardRateCoefficients)
        equilibriumConstants = numpy.zeros_like(forwardRateCoefficients)
        self.generate_reactant_product_indices()

        # Generate forward and reverse rate coefficients k(T,P)
        self.generate_rate_coefficients(coreReactions, edgeReactions)
        
        for rxnList in [coreReactions, edgeReactions]:
            for rxn in rxnList:
                j = reactionIndex[rxn]
                forwardRateCoefficients[j] = rxn.getRateCoefficient(self.T.value_si, self.P.value_si)
                if rxn.reversible:
                    equilibriumConstants[j] = rxn.getEquilibriumConstant(self.T.value_si)
                    reverseRateCoefficients[j] = forwardRateCoefficients[j] / equilibriumConstants[j]

        ReactionSystem.compute_network_variables(pdepNetworks)
        
        # Set initial conditions
        self.set_initial_conditions()
        
        # Initialize the model
        DASx.initialize(self, t0, y0, dydt0, senpar, atol_array, rtol_array)

    def generate_rate_coefficients(self, list coreReactions, list edgeReactions):
        """
        Populates the forwardRateCoefficients, reverseRateCoefficients and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reacion system.
        """
        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = reactionIndex[rxn]
            forwardRateCoefficients[j] = rxn.getRateCoefficient(self.T.value_si, self.P.value_si)
            if rxn.reversible:
                equilibriumConstants[j] = rxn.getEquilibriumConstant(self.T.value_si)
                reverseRateCoefficients[j] = forwardRateCoefficients[j] / equilibriumConstants[j]

    def set_initial_conditions(self):
        ReactionSystem.set_initial_conditions()

        V = 1.0 / numpy.sum(self.coreSpeciesConcentrations)
        self.V = V  #: volume (m3) required to contain one mole total of core species at start

        for spec, conc in self.initialConcentrations.iteritems():
            self.coreSpeciesConcentrations[self.speciesIndex[spec]] = conc
        
            for j in range(self.numCoreSpecies):
                self.y0[j] = self.coreSpeciesConcentrations[j] * V


        self.dydt0 = - self.residual(self.t0, self.y0, numpy.zeros(neq, numpy.float64), self.senpar)[0]

    @cython.boundscheck(False)
    def residual(self, double t, numpy.ndarray[numpy.float64_t, ndim=1] y, numpy.ndarray[numpy.float64_t, ndim=1] dydt, numpy.ndarray[numpy.float64_t, ndim=1] senpar = numpy.zeros(1, numpy.float64)):

        """
        Return the residual function for the governing DAE system for the
        simple reaction system.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip, inet
        cdef numpy.ndarray[numpy.float64_t, ndim=1] res, kf, kr, knet, delta, equilibriumConstants
        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks
        cdef int i, j, z, first, second, third
        cdef double k, V, reactionRate
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesConcentrations, coreSpeciesRates, coreReactionRates, edgeSpeciesRates, edgeReactionRates, networkLeakRates
        cdef numpy.ndarray[numpy.float64_t, ndim=1] C
        cdef numpy.ndarray[numpy.float64_t, ndim=2] jacobian, dgdk

        ir = self.reactantIndices
        ip = self.productIndices
        equilibriumConstants = self.equilibriumConstants

        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients
        
        inet = self.networkIndices
        knet = self.networkLeakCoefficients

        numCoreSpecies = len(self.coreSpeciesRates)
        numCoreReactions = len(self.coreReactionRates)
        numEdgeSpecies = len(self.edgeSpeciesRates)
        numEdgeReactions = len(self.edgeReactionRates)
        numPdepNetworks = len(self.networkLeakRates)
        
        
        res = numpy.zeros(numCoreSpecies, numpy.float64)

        coreSpeciesConcentrations = numpy.zeros_like(self.coreSpeciesConcentrations)
        coreSpeciesRates = numpy.zeros_like(self.coreSpeciesRates)
        coreReactionRates = numpy.zeros_like(self.coreReactionRates)
        edgeSpeciesRates = numpy.zeros_like(self.edgeSpeciesRates)
        edgeReactionRates = numpy.zeros_like(self.edgeReactionRates)
        networkLeakRates = numpy.zeros_like(self.networkLeakRates)

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        V =  self.V # constant volume reactor

        for j in range(numCoreSpecies):
            C[j] = y[j] / V
            coreSpeciesConcentrations[j] = C[j]
        
        for j in range(ir.shape[0]):
            k = kf[j]
            if ir[j,0] >= numCoreSpecies or ir[j,1] >= numCoreSpecies or ir[j,2] >= numCoreSpecies:
                reactionRate = 0.0
            elif ir[j,1] == -1: # only one reactant
                reactionRate = k * C[ir[j,0]]
            elif ir[j,2] == -1: # only two reactants
                reactionRate = k * C[ir[j,0]] * C[ir[j,1]]
            else: # three reactants!! (really?)
                reactionRate = k * C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]]
            k = kr[j]
            if ip[j,0] >= numCoreSpecies or ip[j,1] >= numCoreSpecies or ip[j,2] >= numCoreSpecies:
                pass
            elif ip[j,1] == -1: # only one reactant
                reactionRate -= k * C[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                reactionRate -= k * C[ip[j,0]] * C[ip[j,1]]
            else: # three reactants!! (really?)
                reactionRate -= k * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]]

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
                reactionRate = k * C[inet[j,0]]
            elif inet[j,2] == -1: # only two reactants
                reactionRate = k * C[inet[j,0]] * C[inet[j,1]]
            else: # three reactants!! (really?)
                reactionRate = k * C[inet[j,0]] * C[inet[j,1]] * C[inet[j,2]]
            networkLeakRates[j] = reactionRate

        self.coreSpeciesConcentrations = coreSpeciesConcentrations
        self.coreSpeciesRates = coreSpeciesRates
        self.coreReactionRates = coreReactionRates
        self.edgeSpeciesRates = edgeSpeciesRates
        self.edgeReactionRates = edgeReactionRates
        self.networkLeakRates = networkLeakRates

        res = coreSpeciesRates * V 
        
        
        if self.sensitivity:
            delta = numpy.zeros(len(y), numpy.float64)
            delta[:numCoreSpecies] = res
            if self.jacobianMatrix is None:
                jacobian = self.jacobian(t,y,dydt,0,senpar)
            else:
                jacobian = self.jacobianMatrix
            dgdk = self.computeRateDerivative()
            for j in range(numCoreReactions+numCoreSpecies):
                for i in range(numCoreSpecies):
                    for z in range(numCoreSpecies):
                        delta[(j+1)*numCoreSpecies + i] += jacobian[i,z]*y[(j+1)*numCoreSpecies + z] 
                    delta[(j+1)*numCoreSpecies + i] += dgdk[i,j]

        else:
            delta = res
        delta = delta - dydt
        
        # Return DELTA, IRES.  IRES is set to 1 in order to tell DASPK to evaluate the sensitivity residuals
        return delta, 1
        
    @cython.boundscheck(False)
    def computeRateDerivative(self):
        """
        Returns derivative vector df/dk_j where dy/dt = f(y, t, k) and
        k_j is the rate parameter for the jth core reaction.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip
        cdef numpy.ndarray[numpy.float64_t, ndim=1] kf, kr, C, deriv
        cdef numpy.ndarray[numpy.float64_t, ndim=2] rateDeriv
        cdef double fderiv, rderiv, flux, V
        cdef int j, numCoreReactions, numCoreSpecies
        
        cdef double RT_inverse, gderiv
        
        ir = self.reactantIndices
        ip = self.productIndices
        
        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients
        
        numCoreReactions = len(self.coreReactionRates)
        numCoreSpecies = len(self.coreSpeciesConcentrations)      
        
        # Get constant volume of reactor
        RT_inverse = 1/(constants.R * self.T.value_si)
        V = self.V

        C = self.coreSpeciesConcentrations
        
        rateDeriv = numpy.zeros((numCoreSpecies,numCoreReactions+numCoreSpecies), numpy.float64)
        
        for j in range(numCoreReactions):
            if ir[j,1] == -1: # only one reactant
                fderiv = C[ir[j,0]]
            elif ir[j,2] == -1: # only two reactants
                fderiv = C[ir[j,0]] * C[ir[j,1]]                             
            else: # three reactants!! (really?)
                fderiv = C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]]          
                
            if ip[j,1] == -1: # only one reactant
                rderiv = kr[j] / kf [j] * C[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                rderiv = kr[j] / kf [j] * C[ip[j,0]] * C[ip[j,1]]
            else: # three reactants!! (really?)
                rderiv = kr[j] / kf [j] * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]]
    
            flux = fderiv - rderiv
            gderiv = rderiv * kf[j] * RT_inverse
            
            deriv = numpy.zeros(numCoreSpecies, numpy.float64) # derivative for reaction j with respect to dG_species i

            deriv[ir[j,0]] += gderiv
            if ir[j,1] != -1: # only two reactants
                deriv[ir[j,1]] += gderiv
                if ir[j,2] != -1: # three reactants!! (really?)
                    deriv[ir[j,2]] += gderiv
            
            deriv[ip[j,0]] -= gderiv
            if ip[j,1] != -1: # only two reactants
                deriv[ip[j,1]] -= gderiv
                if ip[j,2] != -1: # three reactants!! (really?)
                    deriv[ip[j,2]] -= gderiv
            
            rateDeriv[ir[j,0], j] -= flux
            rateDeriv[ir[j,0], numCoreReactions:numCoreReactions+numCoreSpecies] -= deriv
            if ir[j,1] != -1:
                rateDeriv[ir[j,1], j] -= flux
                rateDeriv[ir[j,1], numCoreReactions:numCoreReactions+numCoreSpecies] -= deriv
                if ir[j,2] != -1:
                    rateDeriv[ir[j,2], j] -= flux
                    rateDeriv[ir[j,2], numCoreReactions:numCoreReactions+numCoreSpecies] -= deriv
                
            rateDeriv[ip[j,0], j] += flux
            rateDeriv[ip[j,0], numCoreReactions:numCoreReactions+numCoreSpecies] += deriv
            if ip[j,1] != -1:
                rateDeriv[ip[j,1], j] += flux
                rateDeriv[ip[j,1], numCoreReactions:numCoreReactions+numCoreSpecies] += deriv
                if ip[j,2] != -1:
                    rateDeriv[ip[j,2], j] += flux          
                    rateDeriv[ip[j,2], numCoreReactions:numCoreReactions+numCoreSpecies] += deriv
                        
        rateDeriv = V * rateDeriv
                
        return rateDeriv