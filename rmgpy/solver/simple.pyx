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

cdef class SimpleReactor(ReactionSystem):
    """
    A reaction system consisting of a homogeneous, isothermal, isobaric batch
    reactor. These assumptions allow for a number of optimizations that enable
    this solver to complete very rapidly, even for large kinetic models.
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity P
    cdef public double V
    cdef public bint constantVolume
    cdef public dict initialMoleFractions
    cdef public list sensitiveSpecies
    cdef public double sensitivityThreshold
    cdef public list pdepColliderKinetics
    
    cdef public numpy.ndarray reactantIndices
    cdef public numpy.ndarray productIndices
    cdef public numpy.ndarray networkIndices
    cdef public numpy.ndarray forwardRateCoefficients
    cdef public numpy.ndarray reverseRateCoefficients
    cdef public numpy.ndarray equilibriumConstants
    cdef public numpy.ndarray networkLeakCoefficients
    cdef public numpy.ndarray jacobianMatrix
    cdef public numpy.ndarray pdepColliderReactionIndices
    cdef public numpy.ndarray colliderEfficiencies

    def __init__(self, T, P, initialMoleFractions, termination, sensitiveSpecies=None, sensitivityThreshold=1e-3):
        ReactionSystem.__init__(self, termination)
        self.T = Quantity(T)
        self.P = Quantity(P)
        self.initialMoleFractions = initialMoleFractions

        self.V = 0 # will be set in initializeModel
        self.constantVolume = False
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
        self.pdepColliderReactionIndices = None
        self.pdepColliderKinetics = None
        self.colliderEfficiencies = None
        

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (self.__class__, 
            (self.T, self.P, self.initialMoleFractions, self.termination, self.sensitiveSpecies, self.sensitivityThreshold, ))


    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initialMoleFractions dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initialMoleFractions = {}
        for label, moleFrac in self.initialMoleFractions.iteritems():
            initialMoleFractions[speciesDict[label]] = moleFrac
        self.initialMoleFractions = initialMoleFractions

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list pdepNetworks=None, atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions, pdepNetworks, atol, rtol, sensitivity, sens_atol, sens_rtol)

        cdef int i, j, l, index
        cdef double V, T, P, Peff

        cdef numpy.ndarray[numpy.float64_t, ndim=1] y0, y0_coreSpecies
        
        # Set initial conditions
        self.set_initial_conditions()
        
        set_colliders(coreReactions, edgeReactions)

        # Generate forward and reverse rate coefficients k(T,P)
        forwardRateCoefficients = numpy.zeros((numCoreReactions + numEdgeReactions), numpy.float64)
        reverseRateCoefficients = numpy.zeros_like(forwardRateCoefficients)
        equilibriumConstants = numpy.zeros_like(forwardRateCoefficients)

        for rxnList in [coreReactions, edgeReactions]:
            for rxn in rxnList:
                j = self.reactionIndex[rxn]
                for i in range(pdepColliderReactionIndices.shape[0]):
                    if j == pdepColliderReactionIndices[i]:
                        # Calculate effective pressure
                        Peff = P *numpy.sum(colliderEfficiencies[i]*y0_coreSpecies / numpy.sum(y0_coreSpecies))
                        forwardRateCoefficients[j] = rxn.getRateCoefficient(T, Peff)
                else:                    
                    forwardRateCoefficients[j] = rxn.getRateCoefficient(T, P)
                if rxn.reversible:
                    equilibriumConstants[j] = rxn.getEquilibriumConstant(T)
                    reverseRateCoefficients[j] = forwardRateCoefficients[j] / equilibriumConstants[j]

        ReactionSystem.compute_network_variables(pdepNetworks)

        self.reactantIndices = reactantIndices
        self.productIndices = productIndices
        self.forwardRateCoefficients = forwardRateCoefficients
        self.reverseRateCoefficients = reverseRateCoefficients
        self.equilibriumConstants = equilibriumConstants
        
        # Initialize the model
        DASx.initialize(self, self.t0, self.y0, dydt0, self.senpar, self.atol_array, self.rtol_array)

    def generate_rate_coefficients(self, list coreReactions, list edgeReactions):
        """
        Populates the forwardRateCoefficients, reverseRateCoefficients and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reacion system.
        """
        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            for i in xrange(pdepColliderReactionIndices.shape[0]):
                if j == pdepColliderReactionIndices[i]:
                    # Calculate effective pressure
                    Peff = P *numpy.sum(colliderEfficiencies[i]*y0_coreSpecies / numpy.sum(y0_coreSpecies))
                    forwardRateCoefficients[j] = rxn.getRateCoefficient(self.T.value_si, Peff)
            else:                    
                forwardRateCoefficients[j] = rxn.getRateCoefficient(self.T.value_si, self.P.value_si)
            if rxn.reversible:
                equilibriumConstants[j] = rxn.getEquilibriumConstant(self.T.value_si)
                reverseRateCoefficients[j] = forwardRateCoefficients[j] / equilibriumConstants[j]


    def set_colliders(self, coreReactions, edgeReactions):
        """
        Store collider efficiencies and reaction indices for pdep reactions that have specific collider efficiencies
        """
        pdepColliderReactionIndices = []
        pdepColliderKinetics = []
        colliderEfficiencies = []

        for rxn in itertools.chain(coreReactions, edgeReactions):
            if rxn.kinetics.isPressureDependent():
                if rxn.kinetics.efficiencies:
                    j = self.reactionIndex[rxn]
                    pdepColliderReactionIndices.append(j)
                    pdepColliderKinetics.append(rxn.kinetics)
                    colliderEfficiencies.append(rxn.kinetics.getEffectiveColliderEfficiencies(coreSpecies))
        
        self.pdepColliderReactionIndices = numpy.array(pdepColliderReactionIndices, numpy.int)
        self.colliderEfficiencies = numpy.array(colliderEfficiencies, numpy.float64)


    def set_initial_conditions(self):
        ReactionSystem.set_initial_conditions()

        for spec, moleFrac in self.initialMoleFractions.iteritems():
            self.y0[self.speciesIndex[spec]] = moleFrac
        
        # Use ideal gas law to compute volume
        self.V = constants.R * self.T.value_si * numpy.sum(self.y0[:self.numCoreSpecies]) / self.P.value_si# volume in m^3
        for j in range(numCoreSpecies):
            self.coreSpeciesConcentrations[j] = self.y0[j] / self.V

        dydt0 = - self.residual(self.t0, self.y0, numpy.zeros(neq, numpy.float64), self.senpar)[0]

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
        cdef double k, V, reactionRate, T, P, Peff
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesConcentrations, coreSpeciesRates, coreReactionRates, edgeSpeciesRates, edgeReactionRates, networkLeakRates
        cdef numpy.ndarray[numpy.float64_t, ndim=1] C, y_coreSpecies
        cdef numpy.ndarray[numpy.float64_t, ndim=2] jacobian, dgdk, colliderEfficiencies
        cdef numpy.ndarray[numpy.int_t, ndim=1] pdepColliderReactionIndices
        cdef list pdepColliderKinetics

        ir = self.reactantIndices
        ip = self.productIndices
        
        numCoreSpecies = len(self.coreSpeciesRates)
        numCoreReactions = len(self.coreReactionRates)
        numEdgeSpecies = len(self.edgeSpeciesRates)
        numEdgeReactions = len(self.edgeReactionRates)
        numPdepNetworks = len(self.networkLeakRates)
        
        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients
        
        y_coreSpecies = y[:numCoreSpecies]
        
        # Recalculate any forward and reverse rate coefficients that involve pdep collision efficiencies
        if self.pdepColliderReactionIndices.shape[0] != 0:
            T = self.T.value_si
            P = self.P.value_si
            equilibriumConstants = self.equilibriumConstants
            pdepColliderReactionIndices = self.pdepColliderReactionIndices
            pdepColliderKinetics = self.pdepColliderKinetics
            colliderEfficiencies = self.colliderEfficiencies
            for i in range(pdepColliderReactionIndices.shape[0]):
                # Calculate effective pressure
                Peff = P*numpy.sum(colliderEfficiencies[i]*y_coreSpecies / numpy.sum(y_coreSpecies))
                j = pdepColliderReactionIndices[i]
                kf[j] = pdepColliderKinetics[i].getRateCoefficient(T, Peff)
                kr[j] = kf[j] / equilibriumConstants[j]
            
        inet = self.networkIndices
        knet = self.networkLeakCoefficients
        
        
        res = numpy.zeros(numCoreSpecies, numpy.float64)

        coreSpeciesConcentrations = numpy.zeros_like(self.coreSpeciesConcentrations)
        coreSpeciesRates = numpy.zeros_like(self.coreSpeciesRates)
        coreReactionRates = numpy.zeros_like(self.coreReactionRates)
        edgeSpeciesRates = numpy.zeros_like(self.edgeSpeciesRates)
        edgeReactionRates = numpy.zeros_like(self.edgeReactionRates)
        networkLeakRates = numpy.zeros_like(self.networkLeakRates)

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        
        # Use ideal gas law to compute volume
        V = constants.R * self.T.value_si * numpy.sum(y_coreSpecies) / self.P.value_si
        self.V = V

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
    def jacobian(self, double t, numpy.ndarray[numpy.float64_t, ndim=1] y, numpy.ndarray[numpy.float64_t, ndim=1] dydt, double cj, numpy.ndarray[numpy.float64_t, ndim=1] senpar = numpy.zeros(1, numpy.float64)):
        """
        Return the analytical Jacobian for the reaction system.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip
        cdef numpy.ndarray[numpy.float64_t, ndim=1] kf, kr, C
        cdef numpy.ndarray[numpy.float64_t, ndim=2] pd
        cdef int numCoreReactions, numCoreSpecies, i, j
        cdef double k, V, Ctot, deriv, corr
        
        ir = self.reactantIndices
        ip = self.productIndices

        kf = self.forwardRateCoefficients
        kr = self.reverseRateCoefficients
        numCoreReactions = len(self.coreReactionRates)
        numCoreSpecies = len(self.coreSpeciesConcentrations)
        
        pd = -cj * numpy.identity(numCoreSpecies, numpy.float64)
        
        V = constants.R * self.T.value_si * numpy.sum(y[:numCoreSpecies]) / self.P.value_si
        
        Ctot = self.P.value_si /(constants.R * self.T.value_si)

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        for j in range(numCoreSpecies):
            C[j] = y[j] / V

        for j in range(numCoreReactions):
           
            k = kf[j]
            if ir[j,1] == -1: # only one reactant
                deriv = k
                pd[ir[j,0], ir[j,0]] -= deriv
                
                pd[ip[j,0], ir[j,0]] += deriv                
                if ip[j,1] != -1:
                    pd[ip[j,1], ir[j,0]] += deriv
                    if ip[j,2] != -1:
                        pd[ip[j,2], ir[j,0]] += deriv
                
                                
            elif ir[j,2] == -1: # only two reactants
                corr = - k * C[ir[j,0]] * C[ir[j,1]] / Ctot
                if ir[j,0] == ir[j,1]:  # reactants are the same
                    deriv = 2 * k * C[ir[j,0]]
                    pd[ir[j,0], ir[j,0]] -= 2 * deriv
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                    
                    pd[ip[j,0], ir[j,0]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr    
                    
                else:
                    # Derivative with respect to reactant 1
                    deriv = k * C[ir[j, 1]]
                    pd[ir[j,0], ir[j,0]] -= deriv                    
                    pd[ir[j,1], ir[j,0]] -= deriv                        
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv
                    
                    # Derivative with respect to reactant 2
                    deriv = k * C[ir[j, 0]] 
                    pd[ir[j,0], ir[j,1]] -= deriv                    
                    pd[ir[j,1], ir[j,1]] -= deriv                                           
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= corr     
                            
                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr               
                    
                    
            else: # three reactants!! (really?)
                corr = - 2* k * C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]] / Ctot
                if (ir[j,0] == ir[j,1] & ir[j,0] == ir[j,2]):
                    deriv = 3 * k * C[ir[j,0]] * C[ir[j,0]] 
                    pd[ir[j,0], ir[j,0]] -= 3 * deriv                                                           
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 3 * corr
                    
                    pd[ip[j,0], ir[j,0]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr        
                    
                elif ir[j,0] == ir[j,1]:
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ir[j,0]] * C[ir[j,2]]
                    pd[ir[j,0], ir[j,0]] -= 2 * deriv                  
                    pd[ir[j,2], ir[j,0]] -= deriv    
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv
                    
                    # derivative with respect to reactant 3
                    deriv = k * C[ir[j,0]] * C[ir[j,0]] 
                    pd[ir[j,0], ir[j,2]] -= 2 * deriv                  
                    pd[ir[j,2], ir[j,2]] -= deriv                                                                                           
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                        pd[ir[j,2], i] -= corr
                        
                    pd[ip[j,0], ir[j,2]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr    
                    
                    
                elif ir[j,1] == ir[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = k * C[ir[j,1]] * C[ir[j,1]] 
                    pd[ir[j,0], ir[j,0]] -= deriv                    
                    pd[ir[j,1], ir[j,0]] -= 2 * deriv
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv  
                    # derivative with respect to reactant 2
                    deriv = 2 * k * C[ir[j,0]] * C[ir[j,1]]
                    pd[ir[j,0], ir[j,1]] -= deriv                    
                    pd[ir[j,1], ir[j,1]] -= 2 * deriv                                                                                                         
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= 2 * corr

                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr     
                
                elif ir[j,0] == ir[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ir[j,0]] * C[ir[j,1]]
                    pd[ir[j,0], ir[j,0]] -= 2 * deriv                  
                    pd[ir[j,1], ir[j,0]] -= deriv    
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv
                    # derivative with respect to reactant 2
                    deriv = k * C[ir[j,0]] * C[ir[j,0]] 
                    pd[ir[j,0], ir[j,1]] -= 2 * deriv                    
                    pd[ir[j,1], ir[j,1]] -= deriv                                                                                                         
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                        pd[ir[j,1], i] -= corr

                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr     
                                
                else:
                    # derivative with respect to reactant 1
                    deriv = k * C[ir[j,1]] * C[ir[j,2]]
                    pd[ir[j,0], ir[j,0]] -= deriv                    
                    pd[ir[j,1], ir[j,0]] -= deriv
                    pd[ir[j,2], ir[j,0]] -= deriv
                    
                    pd[ip[j,0], ir[j,0]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv     
                                    
                    # derivative with respect to reactant 2
                    deriv = k * C[ir[j,0]] * C[ir[j,2]]
                    pd[ir[j,0], ir[j,1]] -= deriv                    
                    pd[ir[j,1], ir[j,1]] -= deriv   
                    pd[ir[j,2], ir[j,1]] -= deriv
                    
                    pd[ip[j,0], ir[j,1]] += deriv       
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv 
                                 
                    # derivative with respect to reactant 3
                    deriv = k * C[ir[j,0]] * C[ir[j,1]]             
                    pd[ir[j,0], ir[j,2]] -= deriv                    
                    pd[ir[j,1], ir[j,2]] -= deriv   
                    pd[ir[j,2], ir[j,2]] -= deriv
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= corr
                        pd[ir[j,2], i] -= corr
                        
                    pd[ip[j,0], ir[j,2]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ip[j,2], i] += corr     
                    
            
            
            k = kr[j]         
            if ip[j,1] == -1: # only one reactant
                deriv = k
                pd[ip[j,0], ip[j,0]] -= deriv
                
                pd[ir[j,0], ip[j,0]] += deriv                
                if ir[j,1] != -1:
                    pd[ir[j,1], ip[j,0]] += deriv
                    if ir[j,2] != -1:
                        pd[ir[j,2], ip[j,0]] += deriv
                
                                
            elif ip[j,2] == -1: # only two reactants
                corr = -k * C[ip[j,0]] * C[ip[j,1]] / Ctot
                if ip[j,0] == ip[j,1]:
                    deriv = 2 * k * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv                 
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        
                    pd[ir[j,0], ip[j,0]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr   
                    
                else:
                    # Derivative with respect to reactant 1
                    deriv = k * C[ip[j, 1]]
                    pd[ip[j,0], ip[j,0]] -= deriv                    
                    pd[ip[j,1], ip[j,0]] -= deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    
                    # Derivative with respect to reactant 2
                    deriv = k * C[ip[j, 0]] 
                    pd[ip[j,0], ip[j,1]] -= deriv                    
                    pd[ip[j,1], ip[j,1]] -= deriv              
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= corr
                      
                    pd[ir[j,0], ip[j,1]] += deriv                
                    for i in range(numCoreSpecies):
                         pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr              
                    
                    
            else: # three reactants!! (really?)
                corr = - 2 * k * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]] / Ctot
                if (ip[j,0] == ip[j,1] & ip[j,0] == ip[j,2]):
                    deriv = 3 * k * C[ip[j,0]] * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,0]] -= 3 * deriv          
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 3 * corr
                    
                    pd[ir[j,0], ip[j,0]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr       
                    
                elif ip[j,0] == ip[j,1]:
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,2]] 
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv                    
                    pd[ip[j,2], ip[j,0]] -= deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    # derivative with respect to reactant 3
                    deriv = k * C[ip[j,0]] * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,2]] -= 2 * deriv                    
                    pd[ip[j,2], ip[j,2]] -= deriv                       
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        pd[ip[j,2], i] -= corr

                    pd[ir[j,0], ip[j,2]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr     
                    
                    
                elif ip[j,1] == ip[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = k * C[ip[j,1]] * C[ip[j,1]] 
                    pd[ip[j,0], ip[j,0]] -= deriv                    
                    pd[ip[j,1], ip[j,0]] -= 2 * deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv                 
                    
                    # derivative with respect to reactant 2
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,1]] 
                    pd[ip[j,0], ip[j,1]] -= deriv                    
                    pd[ip[j,1], ip[j,1]] -= 2 * deriv   
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= 2 * corr
                        
                    pd[ir[j,0], ip[j,1]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr                    
                                
                elif ip[j,0] == ip[j,2]:                    
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,1]]
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv                  
                    pd[ip[j,1], ip[j,0]] -= deriv    
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    # derivative with respect to reactant 2
                    deriv = k * C[ip[j,0]] * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,1]] -= 2 * deriv                    
                    pd[ip[j,1], ip[j,1]] -= deriv                                                                                                         
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        pd[ip[j,1], i] -= corr

                    pd[ir[j,0], ip[j,1]] += deriv                       
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr    
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv                                               
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr    
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv                                          
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr     
                                
                else:
                    # derivative with respect to reactant 1
                    deriv = k * C[ip[j,1]] * C[ip[j,2]] 
                    pd[ip[j,0], ip[j,0]] -= deriv                    
                    pd[ip[j,1], ip[j,0]] -= deriv
                    pd[ip[j,2], ip[j,0]] -= deriv
                    
                    pd[ir[j,0], ip[j,0]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv     
                                    
                    # derivative with respect to reactant 2
                    deriv = k * C[ip[j,0]] * C[ip[j,2]] 
                    pd[ip[j,0], ip[j,1]] -= deriv                    
                    pd[ip[j,1], ip[j,1]] -= deriv   
                    pd[ip[j,2], ip[j,1]] -= deriv
                    
                    pd[ir[j,0], ip[j,1]] += deriv       
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv 
                                 
                    # derivative with respect to reactant 3
                    deriv = k * C[ip[j,0]] * C[ip[j,1]] 
                    pd[ip[j,0], ip[j,2]] -= deriv                    
                    pd[ip[j,1], ip[j,2]] -= deriv   
                    pd[ip[j,2], ip[j,2]] -= deriv 
                    for i in range(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= corr
                        pd[ip[j,2], i] -= corr
                    
                    pd[ir[j,0], ip[j,2]] += deriv                
                    for i in range(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv          
                        for i in range(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv  
                            for i in range(numCoreSpecies):
                                pd[ir[j,2], i] += corr  

        self.jacobianMatrix = pd + cj * numpy.identity(numCoreSpecies, numpy.float64)
        return pd
    
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
        
        # Use stored volume, since this function is only called from residual function. 
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
                rderiv = kr[j] / kf[j] * C[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                rderiv = kr[j] / kf[j] * C[ip[j,0]] * C[ip[j,1]]
            else: # three reactants!! (really?)
                rderiv = kr[j] / kf[j] * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]]
            
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