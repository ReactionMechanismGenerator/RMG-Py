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
Contains the :class:`LiquidReactor` class, providing a reaction system
consisting of a homogeneous, isothermal, isobaric batch reactor.
"""

import numpy
cimport numpy

import itertools

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
    cdef public list constSPCNames
    cdef public list constSPCIndices
    cdef public dict initialConcentrations
    cdef public list Trange
    cdef public int nSimsTerm
    cdef public dict sensConditions
    
    def __init__(self, T, initialConcentrations, nSimsTerm=None, termination=None, sensitiveSpecies=None, sensitivityThreshold=1e-3, sensConditions=None, constSPCNames=None):
        
        ReactionSystem.__init__(self, termination, sensitiveSpecies, sensitivityThreshold)
        
        if type(T) != list:
            self.T = Quantity(T)
        else:
            self.Trange = [Quantity(t) for t in T]
        
        self.P = Quantity(100000.,'kPa') # Arbitrary high pressure (1000 Bar) to get reactions in the high-pressure limit!
        self.initialConcentrations = initialConcentrations # should be passed in SI
        self.V = 0 # will be set from initialConcentrations in initializeModel
        self.constantVolume = True
        #Constant concentration attributes
        self.constSPCIndices=None
        self.constSPCNames = constSPCNames #store index of constant species 
        self.sensConditions = sensConditions
        self.nSimsTerm = nSimsTerm
        
    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initialConcentrations dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initialConcentrations = {}
        for label, moleFrac in self.initialConcentrations.iteritems():
            if label == 'T':
                continue
            initialConcentrations[speciesDict[label]] = moleFrac
        self.initialConcentrations = initialConcentrations
    
    def get_constSPCIndices (self, coreSpecies):
        "Allow to identify constant Species position in solver"
        for spc in self.constSPCNames:
            if self.constSPCIndices is None: #initialize once the list if constant SPC declared
                self.constSPCIndices=[]
            for iter in coreSpecies: #Need to identify the species object corresponding to the the string written in the input file
                if iter.label == spc:
                    self.constSPCIndices.append(coreSpecies.index(iter))#get 
  
    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list surfaceSpecies=None,
                          list surfaceReactions=None, list pdepNetworks=None, atol=1e-16, rtol=1e-8, sensitivity=False, 
                          sens_atol=1e-6, sens_rtol=1e-4, filterReactions=False, dict conditions=None):
        """
        Initialize a simulation of the liquid reactor using the provided kinetic
        model.
        """
        if surfaceSpecies is None:
            surfaceSpecies = []
        if surfaceReactions is None:
            surfaceReactions = []
                    
        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions, surfaceSpecies, surfaceReactions, 
                                       pdepNetworks, atol, rtol, sensitivity, sens_atol, sens_rtol, filterReactions, conditions)

        # Set initial conditions
        self.set_initial_conditions()

        # Compute reaction thresholds if reaction filtering is turned on
        if filterReactions:
            ReactionSystem.set_initial_reaction_thresholds(self)

        # Generate forward and reverse rate coefficients k(T,P)
        self.generate_rate_coefficients(coreReactions, edgeReactions)

        ReactionSystem.compute_network_variables(self, pdepNetworks)
        
        ReactionSystem.set_initial_derivative(self)

        # Initialize the model
        ReactionSystem.initialize_solver(self)

    def generate_rate_coefficients(self, coreReactions, edgeReactions):
        """
        Populates the forwardRateCoefficients, reverseRateCoefficients and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reacion system.
        """
        
        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            self.kf[j] = rxn.getRateCoefficient(self.T.value_si, self.P.value_si)
            if rxn.reversible:
                self.Keq[j] = rxn.getEquilibriumConstant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]

    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the 
        current reactor model.

        The volume is set to the value in m3 required to contain 
        one mole total of core species at start.

        The coreSpeciesConcentrations array is set to the values stored in the
        initial concentrations dictionary.

        The initial number of moles of a species j is computed and stored in the
        y0 instance attribute.

        """
        ReactionSystem.set_initial_conditions(self)

        for spec, conc in self.initialConcentrations.iteritems():
            i = self.get_species_index(spec)
            self.coreSpeciesConcentrations[i] = conc
        
        V = 1.0 / numpy.sum(self.coreSpeciesConcentrations)
        self.V = V 
        
        for j in xrange(self.numCoreSpecies):
            self.y0[j] = self.coreSpeciesConcentrations[j] * V
       

    @cython.boundscheck(False)
    def residual(self, double t, numpy.ndarray[numpy.float64_t, ndim=1] y, numpy.ndarray[numpy.float64_t, ndim=1] dydt, numpy.ndarray[numpy.float64_t, ndim=1] senpar = numpy.zeros(1, numpy.float64)):

        """
        Return the residual function for the governing DAE system for the
        liquid reaction system.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip, inet
        cdef numpy.ndarray[numpy.float64_t, ndim=1] res, kf, kr, knet, delta, equilibriumConstants
        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks
        cdef int i, j, z, first, second, third
        cdef double k, V, reactionRate
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesConcentrations, coreSpeciesRates, coreReactionRates, edgeSpeciesRates, edgeReactionRates, networkLeakRates, coreSpeciesConsumptionRates, coreSpeciesProductionRates
        cdef numpy.ndarray[numpy.float64_t, ndim=1] C
        cdef numpy.ndarray[numpy.float64_t, ndim=2] jacobian, dgdk

        ir = self.reactantIndices
        ip = self.productIndices
        equilibriumConstants = self.Keq

        kf = self.kf
        kr = self.kb
        
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
        coreSpeciesConsumptionRates = numpy.zeros_like(self.coreSpeciesConsumptionRates)
        coreSpeciesProductionRates = numpy.zeros_like(self.coreSpeciesProductionRates)
        edgeSpeciesRates = numpy.zeros_like(self.edgeSpeciesRates)
        edgeReactionRates = numpy.zeros_like(self.edgeReactionRates)
        networkLeakRates = numpy.zeros_like(self.networkLeakRates)
        

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        V =  self.V # constant volume reactor

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


        #chatelak: Same as in Java, coreSpecies rate = 0 if declared as constatn 
        if self.constSPCIndices is not None:
            for spcIndice in self.constSPCIndices:
                coreSpeciesRates[spcIndice] = 0


        self.coreSpeciesConcentrations = coreSpeciesConcentrations
        self.coreSpeciesRates = coreSpeciesRates
        self.coreReactionRates = coreReactionRates
        self.coreSpeciesProductionRates = coreSpeciesProductionRates
        self.coreSpeciesConsumptionRates = coreSpeciesConsumptionRates
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

        V = self.V  # volume is constant

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
                if ir[j,0] == ir[j,1]:  # reactants are the same
                    deriv = 2 * k * C[ir[j,0]]
                    pd[ir[j,0], ir[j,0]] -= 2 * deriv

                    pd[ip[j,0], ir[j,0]] += deriv
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv

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

                    pd[ip[j,0], ir[j,1]] += deriv
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv


            else: # three reactants!! (really?)
                if (ir[j,0] == ir[j,1] & ir[j,0] == ir[j,2]):
                    deriv = 3 * k * C[ir[j,0]] * C[ir[j,0]]
                    pd[ir[j,0], ir[j,0]] -= 3 * deriv

                    pd[ip[j,0], ir[j,0]] += deriv
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,0]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,0]] += deriv

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

                    pd[ip[j,0], ir[j,2]] += deriv
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv


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

                    pd[ip[j,0], ir[j,1]] += deriv
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv

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

                    pd[ip[j,0], ir[j,1]] += deriv
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,1]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,1]] += deriv

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

                    pd[ip[j,0], ir[j,2]] += deriv
                    if ip[j,1] != -1:
                        pd[ip[j,1], ir[j,2]] += deriv
                        if ip[j,2] != -1:
                            pd[ip[j,2], ir[j,2]] += deriv



            k = kr[j]
            if ip[j,1] == -1: # only one product
                deriv = k
                pd[ip[j,0], ip[j,0]] -= deriv

                pd[ir[j,0], ip[j,0]] += deriv
                if ir[j,1] != -1:
                    pd[ir[j,1], ip[j,0]] += deriv
                    if ir[j,2] != -1:
                        pd[ir[j,2], ip[j,0]] += deriv


            elif ip[j,2] == -1: # only two products
                if ip[j,0] == ip[j,1]:
                    deriv = 2 * k * C[ip[j,0]]
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv

                    pd[ir[j,0], ip[j,0]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv

                else:
                    # Derivative with respect to product 1
                    deriv = k * C[ip[j, 1]]
                    pd[ip[j,0], ip[j,0]] -= deriv
                    pd[ip[j,1], ip[j,0]] -= deriv

                    pd[ir[j,0], ip[j,0]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv

                    # Derivative with respect to product 2
                    deriv = k * C[ip[j, 0]]
                    pd[ip[j,0], ip[j,1]] -= deriv
                    pd[ip[j,1], ip[j,1]] -= deriv

                    pd[ir[j,0], ip[j,1]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv


            else: # three products!! (really?)
                if (ip[j,0] == ip[j,1] & ip[j,0] == ip[j,2]):
                    deriv = 3 * k * C[ip[j,0]] * C[ip[j,0]]
                    pd[ip[j,0], ip[j,0]] -= 3 * deriv

                    pd[ir[j,0], ip[j,0]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv

                elif ip[j,0] == ip[j,1]:
                    # derivative with respect to product 1
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,2]]
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv
                    pd[ip[j,2], ip[j,0]] -= deriv

                    pd[ir[j,0], ip[j,0]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    # derivative with respect to product 3
                    deriv = k * C[ip[j,0]] * C[ip[j,0]]
                    pd[ip[j,0], ip[j,2]] -= 2 * deriv
                    pd[ip[j,2], ip[j,2]] -= deriv

                    pd[ir[j,0], ip[j,2]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv


                elif ip[j,1] == ip[j,2]:
                    # derivative with respect to product 1
                    deriv = k * C[ip[j,1]] * C[ip[j,1]]
                    pd[ip[j,0], ip[j,0]] -= deriv
                    pd[ip[j,1], ip[j,0]] -= 2 * deriv

                    pd[ir[j,0], ip[j,0]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv

                    # derivative with respect to product 2
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,1]]
                    pd[ip[j,0], ip[j,1]] -= deriv
                    pd[ip[j,1], ip[j,1]] -= 2 * deriv

                    pd[ir[j,0], ip[j,1]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv


                elif ip[j,0] == ip[j,2]:
                    # derivative with respect to product 1
                    deriv = 2 * k * C[ip[j,0]] * C[ip[j,1]]
                    pd[ip[j,0], ip[j,0]] -= 2 * deriv
                    pd[ip[j,1], ip[j,0]] -= deriv

                    pd[ir[j,0], ip[j,0]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv
                    # derivative with respect to product 2
                    deriv = k * C[ip[j,0]] * C[ip[j,0]]
                    pd[ip[j,0], ip[j,1]] -= 2 * deriv
                    pd[ip[j,1], ip[j,1]] -= deriv

                    pd[ir[j,0], ip[j,1]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv

                else:
                    # derivative with respect to product 1
                    deriv = k * C[ip[j,1]] * C[ip[j,2]]
                    pd[ip[j,0], ip[j,0]] -= deriv
                    pd[ip[j,1], ip[j,0]] -= deriv
                    pd[ip[j,2], ip[j,0]] -= deriv

                    pd[ir[j,0], ip[j,0]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,0]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,0]] += deriv

                    # derivative with respect to product 2
                    deriv = k * C[ip[j,0]] * C[ip[j,2]]
                    pd[ip[j,0], ip[j,1]] -= deriv
                    pd[ip[j,1], ip[j,1]] -= deriv
                    pd[ip[j,2], ip[j,1]] -= deriv

                    pd[ir[j,0], ip[j,1]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,1]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,1]] += deriv

                    # derivative with respect to product 3
                    deriv = k * C[ip[j,0]] * C[ip[j,1]]
                    pd[ip[j,0], ip[j,2]] -= deriv
                    pd[ip[j,1], ip[j,2]] -= deriv
                    pd[ip[j,2], ip[j,2]] -= deriv

                    pd[ir[j,0], ip[j,2]] += deriv
                    if ir[j,1] != -1:
                        pd[ir[j,1], ip[j,2]] += deriv
                        if ir[j,2] != -1:
                            pd[ir[j,2], ip[j,2]] += deriv

        self.jacobianMatrix = pd + cj * numpy.identity(numCoreSpecies, numpy.float64)
        return pd
