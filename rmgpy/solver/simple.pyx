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
Contains the :class:`SimpleReactor` class, providing a reaction system
consisting of a homogeneous, isothermal, isobaric batch reactor.
"""

import numpy, logging
cimport numpy

import itertools
    
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

    # collider variables

    """
    pdepColliderKinetics:
    an array that contains a reference to the kinetics object of the reaction
    that has pressure dependent kinetics.
    """
    cdef public list pdepColliderKinetics

    """
    colliderEfficiencies:
    an array consisting of array elements, each element corresponding to a reaction.
    Each element is an array with each position in the array corresponding to the collider efficiency
    of the core species. The collider efficiency is set to 1 if the species was not found in the list
    of colliders.
    """
    cdef public numpy.ndarray colliderEfficiencies

    """
    pdepColliderReactionIndices: 
    array that contains the indices of those reactions that 
    have pressure dependent kinetics. E.g. [4, 10, 2, 123]
    """
    cdef public numpy.ndarray pdepColliderReactionIndices

    """
    pdepSpecificColliderKinetics:
    an array that contains a reference to the kinetics object of the reaction
    that has pressure dependent kinetics with a specific species as a third body collider.
    """
    cdef public list pdepSpecificColliderKinetics

    """
    specificColliderSpecies:
    a list that contains object references to species which are specific third body colliders
    in the respective reactions in pdepSpecificColliderReactionIndices.
    """
    cdef public list specificColliderSpecies

    """
    pdepSpecificColliderReactionIndices:
    an array that contains the indices of reactions that have
    a specifcCollider attribyte. E.g. [16, 155, 90]
    """
    cdef public numpy.ndarray pdepSpecificColliderReactionIndices
    
    cdef public dict sensConditions
    
    cdef public list Trange
    cdef public list Prange
    cdef public int nSimsTerm

    def __init__(self, T, P, initialMoleFractions, nSimsTerm=None, termination=None, sensitiveSpecies=None, sensitivityThreshold=1e-3,sensConditions=None):
        ReactionSystem.__init__(self, termination, sensitiveSpecies, sensitivityThreshold)
        
        
        if type(T) != list:
            self.T = Quantity(T)
        else:
            self.Trange = [Quantity(t) for t in T]
            
        if type(P) != list:
            self.P = Quantity(P)
        else:
            self.Prange = [Quantity(p) for p in P]
        
        self.initialMoleFractions = initialMoleFractions

        self.V = 0 # will be set in initializeModel
        self.constantVolume = False

        self.pdepColliderReactionIndices = None
        self.pdepColliderKinetics = None
        self.colliderEfficiencies = None
        self.pdepSpecificColliderReactionIndices = None
        self.pdepSpecificColliderKinetics = None
        self.specificColliderSpecies = None
        self.sensConditions = sensConditions
        self.nSimsTerm = nSimsTerm

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (self.__class__, 
            (self.T, self.P, self.initialMoleFractions, self.nSimsTerm, self.termination))


    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initialMoleFractions dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initialMoleFractions = {}
        for label, moleFrac in self.initialMoleFractions.iteritems():
            initialMoleFractions[speciesDict[label]] = moleFrac
        self.initialMoleFractions = initialMoleFractions

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list surfaceSpecies=None,
                          list surfaceReactions=None, list pdepNetworks=None, atol=1e-16, rtol=1e-8, sensitivity=False, 
                          sens_atol=1e-6, sens_rtol=1e-4, filterReactions=False, dict conditions=None):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """
        
        if surfaceSpecies is None:
            surfaceSpecies = []
        if surfaceReactions is None:
            surfaceReactions = []
        
        
        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies=coreSpecies, coreReactions=coreReactions, edgeSpecies=edgeSpecies, 
                                       edgeReactions=edgeReactions, surfaceSpecies=surfaceSpecies, surfaceReactions=surfaceReactions,
                                       pdepNetworks=pdepNetworks, atol=atol, rtol=rtol, sensitivity=sensitivity, sens_atol=sens_atol, 
                                       sens_rtol=sens_rtol,filterReactions=filterReactions,conditions=conditions)
        
        # Set initial conditions
        self.set_initial_conditions()

        # Compute reaction thresholds if reaction filtering is turned on
        if filterReactions:
            ReactionSystem.set_initial_reaction_thresholds(self)
        
        self.set_colliders(coreReactions, edgeReactions, coreSpecies)
        
        ReactionSystem.compute_network_variables(self, pdepNetworks)

        # Generate forward and reverse rate coefficients k(T,P)
        self.generate_rate_coefficients(coreReactions, edgeReactions)
        
        ReactionSystem.set_initial_derivative(self)
        # Initialize the model
        ReactionSystem.initialize_solver(self)

    def calculate_effective_pressure(self, rxn):
        """
        Computes the effective pressure for a reaction as:

        .. math:: P_{eff} = P * \\sum_i \\frac{y_i * eff_i}{\\sum_j y_j}

        with:
            - P the pressure of the reactor,
            - y the array of initial moles of the core species

        or as:

        .. math:: P_{eff} = \\frac{P * y_{specificCollider}}{\\sum_j y_j}

        if a specificCollider is mentioned.
        """

        y0_coreSpecies = self.y0[:self.numCoreSpecies]
        sum_core_species = numpy.sum(y0_coreSpecies)
        
        j = self.reactionIndex[rxn]
        for i in xrange(self.pdepColliderReactionIndices.shape[0]):
            if j == self.pdepColliderReactionIndices[i]:
                # Calculate effective pressure
                if rxn.specificCollider is None:
                    Peff = self.P.value_si * numpy.sum(self.colliderEfficiencies[i]*y0_coreSpecies / sum_core_species)
                else:
                    logging.debug("Calculating Peff using {0} as a specificCollider".format(rxn.specificCollider))
                    Peff = self.P.value_si * self.y0[self.speciesIndex[rxn.specificCollider]] / sum_core_species
                return Peff
        return self.P.value_si

    def generate_rate_coefficients(self, coreReactions, edgeReactions):
        """
        Populates the forward rate coefficients (kf), reverse rate coefficients (kb)
        and equilibrium constants (Keq) arrays with the values computed at the temperature
        and (effective) pressure of the reaction system.
        """

        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            Peff = self.calculate_effective_pressure(rxn)
            self.kf[j] = rxn.getRateCoefficient(self.T.value_si, Peff)

            if rxn.reversible:
                self.Keq[j] = rxn.getEquilibriumConstant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]


    def set_colliders(self, coreReactions, edgeReactions, coreSpecies):
        """
        Store collider efficiencies and reaction indices for pdep reactions that have collider efficiencies,
        and store specific collider indices
        """
        pdepColliderReactionIndices = []
        self.pdepColliderKinetics = []
        colliderEfficiencies = []
        pdepSpecificColliderReactionIndices = []
        self.pdepSpecificColliderKinetics = []
        self.specificColliderSpecies = []

        for rxn in itertools.chain(coreReactions, edgeReactions):
            if rxn.kinetics.isPressureDependent():
                if rxn.kinetics.efficiencies:
                    j = self.reactionIndex[rxn]
                    pdepColliderReactionIndices.append(j)
                    self.pdepColliderKinetics.append(rxn.kinetics)
                    colliderEfficiencies.append(rxn.kinetics.getEffectiveColliderEfficiencies(coreSpecies))
                if rxn.specificCollider:
                    pdepSpecificColliderReactionIndices.append(self.reactionIndex[rxn])
                    self.pdepSpecificColliderKinetics.append(rxn.kinetics)
                    self.specificColliderSpecies.append(rxn.specificCollider)

        self.pdepColliderReactionIndices = numpy.array(pdepColliderReactionIndices, numpy.int)
        self.colliderEfficiencies = numpy.array(colliderEfficiencies, numpy.float64)
        self.pdepSpecificColliderReactionIndices = numpy.array(pdepSpecificColliderReactionIndices, numpy.int)


    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the 
        current reactor model.

        The volume is set to the value derived from the ideal gas law, using the 
        user-defined pressure, temperature, and the number of moles of initial species.

        The species moles array (y0) is set to the values stored in the
        initial mole fractions dictionary.

        The initial species concentration is computed and stored in the
        coreSpeciesConcentrations array.

        """

        ReactionSystem.set_initial_conditions(self)

        for spec, moleFrac in self.initialMoleFractions.iteritems():
            i = self.get_species_index(spec)
            self.y0[i] = moleFrac
        
        # Use ideal gas law to compute volume
        self.V = constants.R * self.T.value_si * numpy.sum(self.y0[:self.numCoreSpecies]) / self.P.value_si# volume in m^3
        for j in xrange(self.numCoreSpecies):
            self.coreSpeciesConcentrations[j] = self.y0[j] / self.V

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
        cdef double k, V, reactionRate, revReactionRate, T, P, Peff
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesConcentrations, coreSpeciesRates, coreReactionRates, edgeSpeciesRates, edgeReactionRates, networkLeakRates, coreSpeciesConsumptionRates, coreSpeciesProductionRates
        cdef numpy.ndarray[numpy.float64_t, ndim=1] C, y_coreSpecies
        cdef numpy.ndarray[numpy.float64_t, ndim=2] jacobian, dgdk, colliderEfficiencies
        cdef numpy.ndarray[numpy.int_t, ndim=1] pdepColliderReactionIndices, pdepSpecificColliderReactionIndices
        cdef list pdepColliderKinetics, pdepSpecificColliderKinetics

        ir = self.reactantIndices
        ip = self.productIndices
        
        numCoreSpecies = len(self.coreSpeciesRates)
        numCoreReactions = len(self.coreReactionRates)
        numEdgeSpecies = len(self.edgeSpeciesRates)
        numEdgeReactions = len(self.edgeReactionRates)
        numPdepNetworks = len(self.networkLeakRates)
        kf = self.kf
        kr = self.kb
        
        y_coreSpecies = y[:numCoreSpecies]
        
        # Recalculate any forward and reverse rate coefficients that involve pdep collision efficiencies
        if self.pdepColliderReactionIndices.shape[0] != 0:
            T = self.T.value_si
            P = self.P.value_si
            equilibriumConstants = self.Keq
            pdepColliderReactionIndices = self.pdepColliderReactionIndices
            pdepColliderKinetics = self.pdepColliderKinetics
            colliderEfficiencies = self.colliderEfficiencies
            for i in xrange(pdepColliderReactionIndices.shape[0]):
                # Calculate effective pressure
                Peff = P*numpy.sum(colliderEfficiencies[i]*y_coreSpecies / numpy.sum(y_coreSpecies))
                j = pdepColliderReactionIndices[i]
                kf[j] = pdepColliderKinetics[i].getRateCoefficient(T, Peff)
                kr[j] = kf[j] / equilibriumConstants[j]
        if self.pdepSpecificColliderReactionIndices.shape[0] != 0:
            T = self.T.value_si
            P = self.P.value_si
            equilibriumConstants = self.Keq
            pdepSpecificColliderReactionIndices = self.pdepSpecificColliderReactionIndices
            pdepSpecificColliderKinetics = self.pdepSpecificColliderKinetics
            specificColliderSpecies = self.specificColliderSpecies
            for i in xrange(pdepSpecificColliderReactionIndices.shape[0]):
                # Calculate effective pressure
                Peff = P * y[self.speciesIndex[specificColliderSpecies[i]]] / numpy.sum(y_coreSpecies)
                j = pdepSpecificColliderReactionIndices[i]
                kf[j] = pdepSpecificColliderKinetics[i].getRateCoefficient(T, Peff)
                kr[j] = kf[j] / equilibriumConstants[j]
            
        inet = self.networkIndices
        knet = self.networkLeakCoefficients
        
        
        res = numpy.zeros(numCoreSpecies, numpy.float64)

        coreSpeciesConcentrations = numpy.zeros_like(self.coreSpeciesConcentrations)
        coreSpeciesRates = numpy.zeros_like(self.coreSpeciesRates)
        coreReactionRates = numpy.zeros_like(self.coreReactionRates)
        coreSpeciesConsumptionRates = numpy.zeros_like(self.coreSpeciesConsumptionRates)
        coreSpeciesProductionRates = numpy.zeros_like(self.coreSpeciesProductionRates)
        edgeSpeciesRates = numpy.zeros_like(self.edgeSpeciesRates)
        edgeReactionRates = numpy.zeros_like(self.edgeReactionRates)
        networkLeakRates = numpy.zeros_like(self.networkLeakRates)

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        
        # Use ideal gas law to compute volume
        V = constants.R * self.T.value_si * numpy.sum(y_coreSpecies) / self.P.value_si
        self.V = V

        for j in xrange(numCoreSpecies):
            C[j] = y[j] / V
            coreSpeciesConcentrations[j] = C[j]
        
        for j in xrange(ir.shape[0]):
            k = kf[j]
            if ir[j,0] >= numCoreSpecies or ir[j,1] >= numCoreSpecies or ir[j,2] >= numCoreSpecies:
                fReactionRate = 0.0
            elif ir[j,1] == -1: # only one reactant
                fReactionRate = k * C[ir[j,0]]
            elif ir[j,2] == -1: # only two reactants
                fReactionRate = k * C[ir[j,0]] * C[ir[j,1]]
            else: # three reactants!! (really?)
                fReactionRate = k * C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]]
            k = kr[j]
            if ip[j,0] >= numCoreSpecies or ip[j,1] >= numCoreSpecies or ip[j,2] >= numCoreSpecies:
                revReactionRate = 0.0
            elif ip[j,1] == -1: # only one reactant
                revReactionRate = k * C[ip[j,0]]
            elif ip[j,2] == -1: # only two reactants
                revReactionRate = k * C[ip[j,0]] * C[ip[j,1]]
            else: # three reactants!! (really?)
                revReactionRate = k * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]]
                
            reactionRate = fReactionRate-revReactionRate
            
            # Set the reaction and species rates
            if j < numCoreReactions:
                # The reaction is a core reaction
                coreReactionRates[j] = reactionRate

                # Add/substract the total reaction rate from each species rate
                # Since it's a core reaction we know that all of its reactants
                # and products are core species
                first = ir[j,0]
                coreSpeciesRates[first] -= reactionRate
                coreSpeciesConsumptionRates[first] += fReactionRate
                coreSpeciesProductionRates[first] += revReactionRate
                second = ir[j,1]
                if second != -1:
                    coreSpeciesRates[second] -= reactionRate
                    coreSpeciesConsumptionRates[second] += fReactionRate
                    coreSpeciesProductionRates[second] += revReactionRate
                    third = ir[j,2]
                    if third != -1:
                        coreSpeciesRates[third] -= reactionRate
                        coreSpeciesConsumptionRates[third] += fReactionRate
                        coreSpeciesProductionRates[third] += revReactionRate
                first = ip[j,0]
                coreSpeciesRates[first] += reactionRate
                coreSpeciesProductionRates[first] += fReactionRate
                coreSpeciesConsumptionRates[first] += revReactionRate
                second = ip[j,1]
                if second != -1:
                    coreSpeciesRates[second] += reactionRate
                    coreSpeciesProductionRates[second] += fReactionRate
                    coreSpeciesConsumptionRates[second] += revReactionRate
                    third = ip[j,2]
                    if third != -1:
                        coreSpeciesRates[third] += reactionRate
                        coreSpeciesProductionRates[third] += fReactionRate
                        coreSpeciesConsumptionRates[third] += revReactionRate

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

        for j in xrange(inet.shape[0]):
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
        self.coreSpeciesProductionRates = coreSpeciesProductionRates
        self.coreSpeciesConsumptionRates = coreSpeciesConsumptionRates
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
            dgdk = ReactionSystem.computeRateDerivative(self)
            for j in xrange(numCoreReactions+numCoreSpecies):
                for i in xrange(numCoreSpecies):
                    for z in xrange(numCoreSpecies):
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

        kf = self.kf
        kr = self.kb
        numCoreReactions = len(self.coreReactionRates)
        numCoreSpecies = len(self.coreSpeciesConcentrations)
        
        pd = -cj * numpy.identity(numCoreSpecies, numpy.float64)
        
        V = constants.R * self.T.value_si * numpy.sum(y[:numCoreSpecies]) / self.P.value_si
        
        Ctot = self.P.value_si /(constants.R * self.T.value_si)

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        for j in xrange(numCoreSpecies):
            C[j] = y[j] / V

        for j in xrange(numCoreReactions):
           
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
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                    
                    pd[ip[j,0], ir[j,0]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv                                          
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= corr     
                            
                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in xrange(numCoreSpecies):
                                pd[ip[j,2], i] += corr               
                    
                    
            else: # three reactants!! (really?)
                corr = - 2* k * C[ir[j,0]] * C[ir[j,1]] * C[ir[j,2]] / Ctot
                if (ir[j,0] == ir[j,1] & ir[j,0] == ir[j,2]):
                    deriv = 3 * k * C[ir[j,0]] * C[ir[j,0]] 
                    pd[ir[j,0], ir[j,0]] -= 3 * deriv                                                           
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] -= 3 * corr
                    
                    pd[ip[j,0], ir[j,0]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv                                          
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                        pd[ir[j,2], i] -= corr
                        
                    pd[ip[j,0], ir[j,2]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv                                          
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= 2 * corr

                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] -= 2 * corr
                        pd[ir[j,1], i] -= corr

                    pd[ip[j,0], ir[j,1]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv                                          
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] -= corr
                        pd[ir[j,1], i] -= corr
                        pd[ir[j,2], i] -= corr
                        
                    pd[ip[j,0], ir[j,2]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] += corr    
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ip[j,1], i] += corr    
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv                                          
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        
                    pd[ir[j,0], ip[j,0]] += deriv                
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv          
                        for i in xrange(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv  
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= corr
                      
                    pd[ir[j,0], ip[j,1]] += deriv                
                    for i in xrange(numCoreSpecies):
                         pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv          
                        for i in xrange(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv  
                            for i in xrange(numCoreSpecies):
                                pd[ir[j,2], i] += corr              
                    
                    
            else: # three reactants!! (really?)
                corr = - 2 * k * C[ip[j,0]] * C[ip[j,1]] * C[ip[j,2]] / Ctot
                if (ip[j,0] == ip[j,1] & ip[j,0] == ip[j,2]):
                    deriv = 3 * k * C[ip[j,0]] * C[ip[j,0]] 
                    pd[ip[j,0], ip[j,0]] -= 3 * deriv          
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] -= 3 * corr
                    
                    pd[ir[j,0], ip[j,0]] += deriv                
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv          
                        for i in xrange(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv  
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        pd[ip[j,2], i] -= corr

                    pd[ir[j,0], ip[j,2]] += deriv                
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv          
                        for i in xrange(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv  
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= 2 * corr
                        
                    pd[ir[j,0], ip[j,1]] += deriv                
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv          
                        for i in xrange(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv  
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] -= 2 * corr
                        pd[ip[j,1], i] -= corr

                    pd[ir[j,0], ip[j,1]] += deriv                       
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] += corr    
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv                                               
                        for i in xrange(numCoreSpecies):
                            pd[ir[j,1], i] += corr    
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv                                          
                            for i in xrange(numCoreSpecies):
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
                    for i in xrange(numCoreSpecies):
                        pd[ip[j,0], i] -= corr
                        pd[ip[j,1], i] -= corr
                        pd[ip[j,2], i] -= corr
                    
                    pd[ir[j,0], ip[j,2]] += deriv                
                    for i in xrange(numCoreSpecies):
                        pd[ir[j,0], i] += corr   
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv          
                        for i in xrange(numCoreSpecies):
                            pd[ir[j,1], i] += corr   
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv  
                            for i in xrange(numCoreSpecies):
                                pd[ir[j,2], i] += corr  

        self.jacobianMatrix = pd + cj * numpy.identity(numCoreSpecies, numpy.float64)
        return pd
