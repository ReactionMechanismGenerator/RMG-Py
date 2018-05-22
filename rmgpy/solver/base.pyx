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
Contains the :class:`ReactionSystem` class, a base class for all RMG reaction
systems.
"""

cdef extern from "math.h":
    double sqrt(double)

import numpy
cimport numpy
import rmgpy.constants as constants
cimport rmgpy.constants as constants

include "settings.pxi"
if DASPK == 1:
    from pydas.daspk cimport DASPK as DASx
    from pydas.daspk import DASPKError as DASxError
else:
    from pydas.dassl cimport DASSL as DASx
    from pydas.daspk import DASSLError as DASxError
    
import cython
import logging
import csv
import itertools
from cpython cimport bool
from rmgpy.quantity import Quantity
from rmgpy.chemkin import getSpeciesIdentifier
from rmgpy.reaction import Reaction
from rmgpy.species import Species

################################################################################

cdef class ReactionSystem(DASx):
    """
    A base class for all RMG reaction systems.
    """

    def __init__(self, termination=None, sensitiveSpecies=None, sensitivityThreshold=1e-3):
        DASx.__init__(self)

        # reactor state variables:
        self.t0 = 0.0
        self.y0 = None
        self.dydt0 = None

        #  variables that determine the dimensions of arrays and matrices:
        self.numCoreSpecies = -1
        self.numCoreReactions = -1
        self.numEdgeSpecies = -1
        self.numEdgeReactions = -1
        self.numPdepNetworks = -1
        
        """
        The number of differential equations that will 
        be solved simulatenously by the solver.

        """
        self.neq = -1

        # variables that store stoichiometry data
        """
        speciesIndex is a dictionary with species as keys, and 
        values equal to the index integer that is generated
        when enumerating the core and edge species .
        """
        self.speciesIndex = {}

        """
        reactionIndex is a dictionary with reactions as keys, and 
        values equal to the index integer that is generated
        when enumerating the core and edge reactions .
        """
        self.reactionIndex = {}


        """
        A matrix for the reactants and products.

        The matrix has dimensions n x 3, with :
        - n the sum of the number of core and edge reactions,
        - 3 the maximum number of molecules allowed in either the reactant or
            product side of a reaction.

        The matrix element (j,l), with
        - j the index of the reaction and 
        - l the index of the lth reactant/product in the reaction contains the index
            of the molecule in the speciesIndex dictionary.
        """
        self.reactantIndices = None
        self.productIndices = None


        self.networkIndices = None

        # matrices that cache kinetic and rate data
        self.kf = None # forward rate coefficients
        self.kb = None # reverse rate coefficients
        self.Keq = None # equilibrium constants
        self.networkLeakCoefficients = None
        self.jacobianMatrix = None
        
        self.coreSpeciesConcentrations = None
        
        # The reaction and species rates at the current time (in mol/m^3*s)
        self.coreSpeciesRates = None
        self.coreReactionRates = None
        self.coreSpeciesProductionRates = None
        self.coreSpeciesConsumptionRates = None
        self.edgeSpeciesRates = None
        self.edgeReactionRates = None

        self.networkLeakRates = None
        
        #surface indices
        self.surfaceSpeciesIndices = None
        self.surfaceReactionIndices = None
        self.validLayeringIndices = None
        
        # variables that cache maximum rate (ratio) data
        self.maxEdgeSpeciesRateRatios = None
        self.maxNetworkLeakRateRatios = None
        
        #for managing prunable edge species
        self.prunableSpecies = []
        self.prunableNetworks = []
        self.prunableSpeciesIndices = None
        self.prunableNetworkIndices = None

        # sensitivity variables
        self.sensmethod = 2 # sensmethod = 1 for staggered corrector sensitivities, 0 (simultaneous corrector), 2 (staggered direct)
        self.sensitivityCoefficients = None    
        self.sensitiveSpecies = sensitiveSpecies
        self.sensitivityThreshold = sensitivityThreshold
        self.senpar = None

        # tolerance settings

        """
        Arrays containing the absolute and relative tolerances.
        The dimension of the arrays correspond to the number of 
        differential equations in the system.
        """
        self.atol_array = None
        self.rtol_array = None

        self.termination = termination or []
        
        
        # reaction filtration, unimolecularThreshold is a vector with length of number of core species
        # bimolecularThreshold is a square matrix with length of number of core species
        # A value of 1 in the matrix indicates the species is above the threshold to react or participate in those reactions
        self.unimolecularThreshold = None
        self.bimolecularThreshold = None

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (self.__class__, (self.termination,))

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list surfaceSpecies=None,
                          list surfaceReactions=None, list pdepNetworks=None, atol=1e-16, rtol=1e-8, sensitivity=False, 
                          sens_atol=1e-6, sens_rtol=1e-4, filterReactions=False, dict conditions=None):
        """
        Initialize a simulation of the reaction system using the provided
        kinetic model. You will probably want to create your own version of this
        method in the derived class; don't forget to also call the base class
        version, too.
        """
        if surfaceSpecies is None:
            surfaceSpecies = []
        if surfaceReactions is None:
            surfaceReactions = []
            
        if conditions:
            isConc = hasattr(self,'initialConcentrations')
            keys = conditions.keys()
            if 'T' in keys and hasattr(self,'T'):
                self.T = Quantity(conditions['T'],'K')
            if 'P' in keys and hasattr(self,'P'):
                self.P = Quantity(conditions['P'],'Pa')
            for k in keys:
                if isConc:
                    if k in self.initialConcentrations.keys():
                        self.initialConcentrations[k] = conditions[k] #already in SI units
                else:
                    if k in self.initialMoleFractions.keys():
                        self.initialMoleFractions[k] = conditions[k]              
                    
        self.numCoreSpecies = len(coreSpecies)
        self.numCoreReactions = len(coreReactions)
        self.numEdgeSpecies = len(edgeSpecies)
        self.numEdgeReactions = len(edgeReactions)

        self.initiate_tolerances(atol, rtol, sensitivity, sens_atol, sens_rtol)

        pdepNetworks = pdepNetworks or []
        self.numPdepNetworks = len(pdepNetworks)

        self.kf = numpy.zeros((self.numCoreReactions + self.numEdgeReactions), numpy.float64)
        self.kb = numpy.zeros_like(self.kf)
        self.Keq = numpy.zeros_like(self.kf)

        self.generate_species_indices(coreSpecies, edgeSpecies)
        self.generate_reaction_indices(coreReactions, edgeReactions)
        self.generate_reactant_product_indices(coreReactions, edgeReactions)

        self.coreSpeciesConcentrations = numpy.zeros((self.numCoreSpecies), numpy.float64)
        self.coreSpeciesProductionRates = numpy.zeros((self.numCoreSpecies), numpy.float64)
        self.coreSpeciesConsumptionRates = numpy.zeros((self.numCoreSpecies), numpy.float64)
        self.coreReactionRates = numpy.zeros((self.numCoreReactions), numpy.float64)
        self.edgeReactionRates = numpy.zeros((self.numEdgeReactions), numpy.float64)
        self.coreSpeciesRates = numpy.zeros((self.numCoreSpecies), numpy.float64)
        self.edgeSpeciesRates = numpy.zeros((self.numEdgeSpecies), numpy.float64)
        self.networkLeakRates = numpy.zeros((self.numPdepNetworks), numpy.float64)
        self.maxEdgeSpeciesRateRatios = numpy.zeros((len(self.prunableSpecies)), numpy.float64)
        self.maxNetworkLeakRateRatios = numpy.zeros((len(self.prunableNetworks)), numpy.float64)
        self.sensitivityCoefficients = numpy.zeros((self.numCoreSpecies, self.numCoreReactions), numpy.float64)
        self.unimolecularThreshold = numpy.zeros((self.numCoreSpecies), bool)
        self.bimolecularThreshold = numpy.zeros((self.numCoreSpecies, self.numCoreSpecies), bool)

        surfaceSpecies,surfaceReactions = self.initialize_surface(coreSpecies,coreReactions,surfaceSpecies,surfaceReactions)
        
        self.set_prunable_indices(edgeSpecies, pdepNetworks)
        
    def initialize_solver(self):
        DASx.initialize(self, self.t0, self.y0, self.dydt0, self.senpar, self.atol_array, self.rtol_array)
    
    def set_prunable_indices(self,edgeSpecies,pdepNetworks):
        cdef object spc
        cdef list temp
        temp = []
        for i,spc in enumerate(self.prunableSpecies):
            try:
                temp.append(edgeSpecies.index(spc))
            except ValueError:
                self.maxEdgeSpeciesRateRatios[i] = numpy.inf #avoid pruning of species that have been moved to core
        
        self.prunableSpeciesIndices = numpy.array(temp)
        
        temp = []
        for i,spc in enumerate(self.prunableNetworks):
            try:
                temp.append(pdepNetworks.index(spc))
            except:
                self.maxNetworkLeakRateRatios[i] = numpy.inf #avoid pruning of lost networks
                
        self.prunableNetworkIndices = numpy.array(temp)
        
    @cython.boundscheck(False)
    cpdef initialize_surface(self,list coreSpecies,list coreReactions,list surfaceSpecies,list surfaceReactions):
        """
        removes surfaceSpecies and surfaceReactions from  until they are self consistent: 
            1) every reaction has one species in the surface
            2) every species participates in a surface reaction
        """
        cdef numpy.ndarray[numpy.int_t,ndim=2] productIndices,reactantIndices
        cdef list surfaceSpeciesIndices, surfaceReactionIndices, removeInds
        cdef set possibleSpeciesIndices
        cdef int i,j
        cdef bool notInSurface
        cdef object obj
        
        logging.info('initializing surface ...')
        
        productIndices = self.productIndices
        reactantIndices = self.reactantIndices
        
        surfaceSpeciesIndices = []
        surfaceReactionIndices = []
        possibleSpeciesIndices = set()
        removeInds = []
        
        for obj in surfaceSpecies:
            surfaceSpeciesIndices.append(coreSpecies.index(obj))
            
        for obj in surfaceReactions:
            surfaceReactionIndices.append(coreReactions.index(obj))
       
        for i in surfaceReactionIndices: #remove surface reactions whose species have been moved to the bulk core
            notInSurface = True
            for j in productIndices[i]:
                possibleSpeciesIndices.add(j)
                if j in surfaceSpeciesIndices:
                    notInSurface = False
            for j in reactantIndices[i]:
                possibleSpeciesIndices.add(j)
                if j in surfaceSpeciesIndices:
                    notInSurface = False
            if notInSurface:
                logging.info('removing disconnected reaction from surface: {0}'.format(str(coreReactions[i])))
                surfaceReactions.remove(coreReactions[i])
                surfaceReactionIndices.remove(i)
        
        possibleSpeciesIndices -= {-1} #remove the -1 indexes that indicate there is no third/second reactant/product
        
        for i in surfaceSpeciesIndices: #remove species without reactions in the surface
            if not(i in possibleSpeciesIndices):
                logging.info('removing disconnected species from surface: {0}'.format(coreSpecies[i].label))
                removeInds.append(i)
        
        for i in removeInds:
            surfaceSpecies.remove(coreSpecies[i])
            surfaceSpeciesIndices.remove(i)
        
        self.surfaceSpeciesIndices = numpy.array(surfaceSpeciesIndices,dtype=numpy.int)
        self.surfaceReactionIndices = numpy.array(surfaceReactionIndices,dtype=numpy.int)
        
        self.validLayeringIndices = self.getLayeringIndices()
        
        surfaceSpecies = [coreSpecies[i] for i in surfaceSpeciesIndices]
        surfaceReactions = [coreReactions[i] for i in surfaceReactionIndices]
        
        logging.info('surface initialization complete')

        return surfaceSpecies,surfaceReactions
        
    def initiate_tolerances(self, atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4):
        """
        Computes the number of differential equations and initializes the tolerance arrays.        
        """
        # Compute number of equations    
        if sensitivity:    
            # Set DASPK sensitivity analysis to ON
            self.sensitivity = True
            # Compute number of variables
            self.neq = self.numCoreSpecies * (self.numCoreReactions + self.numCoreSpecies + 1)
            
            self.atol_array = numpy.ones(self.neq, numpy.float64) * sens_atol
            self.atol_array[:self.numCoreSpecies] = atol
            
            self.rtol_array = numpy.ones(self.neq, numpy.float64) * sens_rtol
            self.rtol_array[:self.numCoreSpecies] = rtol
            
            self.senpar = numpy.zeros(self.numCoreReactions + self.numCoreSpecies, numpy.float64)
            
        else:
            self.neq = self.numCoreSpecies
            
            self.atol_array = numpy.ones(self.neq , numpy.float64) * atol
            self.rtol_array = numpy.ones(self.neq , numpy.float64) * rtol
            
            self.senpar = numpy.zeros(self.numCoreReactions, numpy.float64)

    def get_species_index(self, spc):
        """
        Retrieves the index that is associated with the parameter species
        from the species index dictionary.
        """
        return self.speciesIndex[spc]

    def generate_reactant_product_indices(self, coreReactions, edgeReactions):
        """
        Creates a matrix for the reactants and products.
        """

        self.reactantIndices = -numpy.ones((self.numCoreReactions + self.numEdgeReactions, 3), numpy.int )
        self.productIndices = -numpy.ones_like(self.reactantIndices)

        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            for l, spec in enumerate(rxn.reactants):
                i = self.get_species_index(spec)
                self.reactantIndices[j,l] = i
            for l, spec in enumerate(rxn.products):
                i = self.get_species_index(spec)
                self.productIndices[j,l] = i

    def generate_species_indices(self, coreSpecies, edgeSpecies):
        """
        Assign an index to each species (core first, then edge) and 
        store the (species, index) pair in a dictionary.
        """
        
        for index, spec in enumerate(itertools.chain(coreSpecies, edgeSpecies)):
            self.speciesIndex[spec] = index

    def generate_reaction_indices(self, coreReactions, edgeReactions):
        """
        Assign an index to each reaction (core first, then edge) and 
        store the (reaction, index) pair in a dictionary.
        """
        
        for index, rxn in enumerate(itertools.chain(coreReactions, edgeReactions)):
            self.reactionIndex[rxn] = index

    def set_initial_conditions(self):
        """
        Sets the common initial conditions of the rate equations that 
        represent the reaction system.

        - Sets the initial time of the reaction system to 0
        - Initializes the species moles to a n x 1 array with zeros
        """

        self.t0 = 0.0            

        self.y0 = numpy.zeros(self.neq, numpy.float64)
    
    def set_initial_reaction_thresholds(self):
        
        # Set unimolecular and bimolecular thresholds as true for any concentrations greater than 0
        numCoreSpecies = len(self.coreSpeciesConcentrations)
        for i in xrange(numCoreSpecies):
            if self.coreSpeciesConcentrations[i] > 0:
                self.unimolecularThreshold[i] = True
        for i in xrange(numCoreSpecies):
            for j in xrange(i, numCoreSpecies):
                if self.coreSpeciesConcentrations[i] > 0 and self.coreSpeciesConcentrations[j] > 0:
                    self.bimolecularThreshold[i,j] = True

    def set_initial_derivative(self):
        """
        Sets the derivative of the species moles with respect to the independent variable (time)
        equal to the residual.
        """
        self.dydt0 = - self.residual(self.t0, self.y0, numpy.zeros(self.neq, numpy.float64), self.senpar)[0]

    def compute_network_variables(self, pdepNetworks=None):
        """
        Initialize the arrays containing network information:

        - NetworkLeakCoefficients is a n x 1 array with 
            n the number of pressure-dependent networks.
        - NetworkIndices is a n x 3 matrix with 
            n the number of pressure-dependent networks and 
            3 the maximum number of molecules allowed in either the reactant or
            product side of a reaction. 
        """

        pdepNetworks = pdepNetworks or []

        self.networkIndices = -numpy.ones((self.numPdepNetworks, 3), numpy.int )
        self.networkLeakCoefficients = numpy.zeros((self.numPdepNetworks), numpy.float64)

        for j, network in enumerate(pdepNetworks):
            self.networkLeakCoefficients[j] = network.getLeakCoefficient(self.T.value_si, self.P.value_si)
            for l, spec in enumerate(network.source):
                i = self.get_species_index(spec)
                self.networkIndices[j,l] = i
   
    @cython.boundscheck(False)                               
    cpdef getLayeringIndices(self):
        """
        determines the edge reaction indices that indicate reactions that are valid for movement from
        edge to surface based on the layering constraint
        """
        
        cdef numpy.ndarray[numpy.int_t, ndim=1] surfSpeciesIndices
        cdef int i,j,k,index, numCoreSpecies, numCoreReactions, numEdgeReactions
        cdef list validIndices
        
        surfSpeciesIndices = self.surfaceSpeciesIndices
        numCoreSpecies = self.numCoreSpecies
        numCoreReactions = self.numCoreReactions
        numEdgeReactions = self.numEdgeReactions
        productIndices= self.productIndices
        reactantIndices = self.reactantIndices
        
        validIndices = []
        
        for index in xrange(numEdgeReactions):
            for j in productIndices[index+numCoreReactions]:
                if j in surfSpeciesIndices or j >= numCoreSpecies:
                    break
            else:
                validIndices.append(index)
                continue
            for j in reactantIndices[index+numCoreReactions]:
                if j in surfSpeciesIndices or j >= numCoreSpecies:
                    break
            else:
                validIndices.append(index)
        
        return numpy.array(validIndices)
    
    cpdef addReactionsToSurface(self,list newSurfaceReactions,list newSurfaceReactionInds,list surfaceSpecies,list surfaceReactions,list edgeSpecies):
        """
        moves new surface reactions to the surface
        done after the while loop before the simulate call ends
        """
        cdef object srxn
        cdef int k,sind,numCoreSpecies,i,numCoreReactions
        cdef numpy.ndarray[numpy.int_t, ndim=2] productIndices, reactantIndices
        
        productIndices = self.productIndices
        reactantIndices = self.reactantIndices
        numCoreSpecies = self.numCoreSpecies
        numCoreReactions = self.numCoreReactions
        
        for k in xrange(len(newSurfaceReactions)):
            srxn = newSurfaceReactions[k]
            sind = newSurfaceReactionInds[k]
            surfaceReactions.append(srxn) #add to surface trackers
                        
            for i in productIndices[sind+numCoreReactions]:
                if i >= numCoreSpecies:
                    surfaceSpecies.append(edgeSpecies[i-numCoreSpecies])
            for i in reactantIndices[sind+numCoreReactions]:
                if i >= numCoreSpecies:
                    surfaceSpecies.append(edgeSpecies[i-numCoreSpecies])
                    
        return surfaceSpecies,surfaceReactions

    @cython.boundscheck(False)
    cpdef simulate(self, list coreSpecies, list coreReactions, list edgeSpecies, 
        list edgeReactions,list surfaceSpecies, list surfaceReactions,
        list pdepNetworks=None, bool prune=False, bool sensitivity=False, list sensWorksheet=None, object modelSettings=None,
        object simulatorSettings=None, dict conditions=None):
        """
        Simulate the reaction system with the provided reaction model,
        consisting of lists of core species, core reactions, edge species, and
        edge reactions. As the simulation proceeds the system is monitored for
        validity. If the model becomes invalid (e.g. due to an excessively
        large edge flux), the simulation is interrupted and the object causing
        the model to be invalid is returned. If the simulation completes to
        the desired termination criteria and the model remains valid throughout,
        ``None`` is returned.
        """
        cdef double toleranceKeepInEdge,toleranceMoveToCore,toleranceMoveEdgeReactionToCore,toleranceInterruptSimulation
        cdef double toleranceMoveEdgeReactionToCoreInterrupt,toleranceMoveEdgeReactionToSurface
        cdef double toleranceMoveSurfaceSpeciesToCore,toleranceMoveSurfaceReactionToCore
        cdef double toleranceMoveEdgeReactionToSurfaceInterrupt
        cdef bool ignoreOverallFluxCriterion, filterReactions
        cdef double absoluteTolerance, relativeTolerance, sensitivityAbsoluteTolerance, sensitivityRelativeTolerance
        cdef dict speciesIndex
        cdef list row, tempSurfaceObjects
        cdef list sortedInds, tempNewObjects, tempNewObjectInds, tempNewObjectVals, tempInds
        cdef int index, spcIndex, maxSpeciesIndex, maxNetworkIndex
        cdef int numCoreSpecies, numEdgeSpecies, numPdepNetworks, numCoreReactions
        cdef double stepTime, charRate, maxSpeciesRate, maxNetworkRate, maxEdgeReactionAccum, stdan
        cdef numpy.ndarray[numpy.float64_t, ndim=1] y0 #: Vector containing the number of moles of each species
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesRates, edgeSpeciesRates, networkLeakRates, coreSpeciesProductionRates, coreSpeciesConsumptionRates, totalDivAccumNums
        cdef numpy.ndarray[numpy.float64_t, ndim=1] maxCoreSpeciesRates, maxEdgeSpeciesRates, maxNetworkLeakRates,maxEdgeSpeciesRateRatios, maxNetworkLeakRateRatios
        cdef bint terminated
        cdef object maxSpecies, maxNetwork
        cdef int i, j, k
        cdef numpy.float64_t maxSurfaceDifLnAccumNum, maxSurfaceSpeciesRate, conversion
        cdef int maxSurfaceAccumReactionIndex, maxSurfaceSpeciesIndex
        cdef object maxSurfaceAccumReaction, maxSurfaceSpecies
        cdef numpy.ndarray[numpy.float64_t,ndim=1] surfaceSpeciesProduction, surfaceSpeciesConsumption
        cdef numpy.ndarray[numpy.float64_t,ndim=1] surfaceTotalDivAccumNums, surfaceSpeciesRateRatios
        cdef numpy.ndarray[numpy.float64_t, ndim=1] forwardRateCoefficients, coreSpeciesConcentrations
        cdef double  prevTime, totalMoles, c, volume, RTP, unimolecularThresholdVal, bimolecularThresholdVal, maxCharRate
        cdef bool useDynamicsTemp, firstTime, useDynamics, terminateAtMaxObjects, schanged
        cdef numpy.ndarray[numpy.float64_t, ndim=1] edgeReactionRates
        cdef double reactionRate, production, consumption
        cdef numpy.ndarray[numpy.int_t,ndim=1] surfaceSpeciesIndices, surfaceReactionIndices
        # cython declations for sensitivity analysis
        cdef numpy.ndarray[numpy.int_t, ndim=1] sensSpeciesIndices
        cdef numpy.ndarray[numpy.float64_t, ndim=1] moleSens, dVdk, normSens
        cdef list time_array, normSens_array, newSurfaceReactions, newSurfaceReactionInds, newObjects, newObjectInds
        
        zeroProduction = False
        zeroConsumption = False
        pdepNetworks = pdepNetworks or []
        
        schanged = False
        
        numCoreSpecies = len(coreSpecies)
        numEdgeSpecies = len(edgeSpecies)
        numPdepNetworks = len(pdepNetworks)
        numCoreReactions = len(coreReactions)

        assert set(coreReactions) >= set(surfaceReactions), 'given surface reactions are not a subset of core reactions'
        assert set(coreSpecies) >= set(surfaceSpecies), 'given surface species are not a subset of core species' 

        toleranceKeepInEdge = modelSettings.fluxToleranceKeepInEdge if prune else 0
        toleranceMoveToCore = modelSettings.fluxToleranceMoveToCore
        toleranceMoveEdgeReactionToCore = modelSettings.toleranceMoveEdgeReactionToCore
        toleranceInterruptSimulation = modelSettings.fluxToleranceInterrupt
        toleranceMoveEdgeReactionToCoreInterrupt= modelSettings.toleranceMoveEdgeReactionToCore
        toleranceMoveEdgeReactionToSurface = modelSettings.toleranceMoveEdgeReactionToSurface
        toleranceMoveSurfaceSpeciesToCore = modelSettings.toleranceMoveSurfaceSpeciesToCore
        toleranceMoveSurfaceReactionToCore = modelSettings.toleranceMoveSurfaceReactionToCore
        toleranceMoveEdgeReactionToSurfaceInterrupt = modelSettings.toleranceMoveEdgeReactionToSurfaceInterrupt
        ignoreOverallFluxCriterion=modelSettings.ignoreOverallFluxCriterion
        absoluteTolerance = simulatorSettings.atol
        relativeTolerance = simulatorSettings.rtol
        sensitivityAbsoluteTolerance = simulatorSettings.sens_atol
        sensitivityRelativeTolerance = simulatorSettings.sens_rtol
        filterReactions = modelSettings.filterReactions
        maxNumObjsPerIter = modelSettings.maxNumObjsPerIter

        #if not pruning always terminate at max objects, otherwise only do so if terminateAtMaxObjects=True
        terminateAtMaxObjects = True if not prune else modelSettings.terminateAtMaxObjects 

        dynamicsTimeScale = modelSettings.dynamicsTimeScale
        
        useDynamics = not (toleranceMoveEdgeReactionToCore == numpy.inf and toleranceMoveEdgeReactionToSurface == numpy.inf)
        
        speciesIndex = {}
        for index, spec in enumerate(coreSpecies):
            speciesIndex[spec] = index
        
        self.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions, surfaceSpecies, surfaceReactions, 
                             pdepNetworks, absoluteTolerance, relativeTolerance, sensitivity, sensitivityAbsoluteTolerance, 
                             sensitivityRelativeTolerance, filterReactions,conditions)
        
        prunableSpeciesIndices = self.prunableSpeciesIndices
        prunableNetworkIndices = self.prunableNetworkIndices
        
        surfaceSpeciesIndices = self.surfaceSpeciesIndices
        surfaceReactionIndices = self.surfaceReactionIndices
        
        totalDivAccumNums = None #the product of the ratios between accumulation numbers with and without a given reaction for products and reactants

        invalidObjects = []
        newSurfaceReactions = []
        newSurfaceReactionInds = []
        terminated = False
        maxSpeciesIndex = -1
        maxSpecies = None
        maxSpeciesRate = 0.0
        maxNetworkIndex = -1
        maxNetwork = None
        maxNetworkRate = 0.0
        iteration = 0
        conversion = 0.0
        maxCharRate = 0.0
        
        maxEdgeSpeciesRateRatios = self.maxEdgeSpeciesRateRatios
        maxNetworkLeakRateRatios = self.maxNetworkLeakRateRatios
        forwardRateCoefficients = self.kf
        unimolecularThreshold = self.unimolecularThreshold
        bimolecularThreshold = self.bimolecularThreshold
        
        # Copy the initial conditions to use in evaluating conversions
        y0 = self.y.copy()
        
        # a list with the time, Volume, mole fractions of core species
        self.snapshots = []

        if sensitivity:
            time_array = []
            normSens_array = [[] for spec in self.sensitiveSpecies]    
            RTP = constants.R * self.T.value_si / self.P.value_si
            # identify sensitive species indices
            sensSpeciesIndices = numpy.array([speciesIndex[spec] for spec in self.sensitiveSpecies], numpy.int)  # index within coreSpecies list of the sensitive species
                
        
        stepTime = 1e-12
        prevTime = self.t

        firstTime = True
        
        while not terminated:
            # Integrate forward in time by one time step
            
            if not firstTime:
                try:
                    self.step(stepTime)
                except DASxError as e:
                    logging.error("Trying to step from time {} to {} resulted in a solver (DASPK) error".format(prevTime, stepTime))
                    
                    logging.info('Resurrecting Model...')
                    
                    conversion = 0.0
                    for term in self.termination:
                        if isinstance(term, TerminationConversion):
                            index = speciesIndex[term.species]
                            conversion = 1-(y_coreSpecies[index] / y0[index])

                    if invalidObjects == []:
                        #species flux criterion
                        if len(edgeSpeciesRateRatios) > 0:
                            ind = numpy.argmax(edgeSpeciesRateRatios)
                            obj = edgeSpecies[ind]
                            logging.info('At time {0:10.4e} s, species {1} at rate ratio {2} was added to model core in model resurrection process'.format(self.t, obj,edgeSpeciesRates[ind]))
                            invalidObjects.append(obj)
                        
                        if totalDivAccumNums and len(totalDivAccumNums) > 0: #if dynamics data available
                            ind = numpy.argmax(totalDivAccumNums)
                            obj = edgeReactions[ind]
                            logging.info('At time {0:10.4e} s, Reaction {1} at dynamics number {2} was added to model core in model resurrection process'.format(self.t, obj,totalDivAccumNums[ind]))
                            invalidObjects.append(obj)
                        
                        if pdepNetworks != [] and networkLeakRateRatios != []:
                            ind = numpy.argmax(networkLeakRateRatios)
                            obj = pdepNetworks[ind]
                            logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} at network leak rate {2} was sent for exploring during model resurrection process'.format(self.t, obj.index, networkLeakRateRatios[ind]))
                            invalidObjects.append(obj)
                    
                    if invalidObjects != []:
                        return False,True,invalidObjects,surfaceSpecies,surfaceReactions,self.t,conversion
                    else:
                        logging.error('Model Resurrection has failed')
                        logging.error("Core species names: {!r}".format([getSpeciesIdentifier(s) for s in coreSpecies]))
                        logging.error("Core species moles: {!r}".format(self.y[:numCoreSpecies]))
                        logging.error("Volume: {!r}".format(self.V))
                        logging.error("Core species net rates: {!r}".format(self.coreSpeciesRates))
                        logging.error("Edge species net rates: {!r}".format(self.edgeSpeciesRates))
                        logging.error("Network leak rates: {!r}".format(self.networkLeakRates))
                        raise ValueError('invalidObjects could not be filled during resurrection process')
            
            y_coreSpecies = self.y[:numCoreSpecies]
            totalMoles = numpy.sum(y_coreSpecies)
            if sensitivity:
                time_array.append(self.t)
                moleSens = self.y[numCoreSpecies:]#   
                volume = self.V
                
                dVdk = numpy.zeros(numCoreReactions + numCoreSpecies, numpy.float64)
                if not self.constantVolume:
                    for j in xrange(numCoreReactions + numCoreSpecies):
                        dVdk[j] = numpy.sum(moleSens[j*numCoreSpecies:(j+1)*numCoreSpecies])*RTP   # Contains [ dV_dk and dV_dG ]
                for i in xrange(len(self.sensitiveSpecies)):
                    normSens = numpy.zeros(numCoreReactions + numCoreSpecies, numpy.float64)
                    c = self.coreSpeciesConcentrations[sensSpeciesIndices[i]]
                    if c != 0:                        
                        for j in xrange(numCoreReactions):
                            normSens[j] = 1/volume*(moleSens[j*numCoreSpecies+sensSpeciesIndices[i]]-c*dVdk[j])*forwardRateCoefficients[j]/c
                        for j in xrange(numCoreReactions,numCoreReactions+numCoreSpecies):
                            normSens[j] = 1/volume*(moleSens[j*numCoreSpecies+sensSpeciesIndices[i]]-c*dVdk[j])/c*4184   # no normalization against dG, converstion to kcal/mol units
                    normSens_array[i].append(normSens)


            snapshot = [self.t, self.V]
            snapshot.extend(y_coreSpecies / numpy.sum(y_coreSpecies))
            self.snapshots.append(snapshot)            

            # Get the characteristic flux
            charRate = sqrt(numpy.sum(self.coreSpeciesRates * self.coreSpeciesRates))
            
            if charRate > maxCharRate:
                maxCharRate = charRate
                
            coreSpeciesRates = numpy.abs(self.coreSpeciesRates)
            edgeReactionRates = self.edgeReactionRates
            coreSpeciesConsumptionRates = self.coreSpeciesConsumptionRates
            coreSpeciesProductionRates = self.coreSpeciesProductionRates
            edgeSpeciesRates = numpy.abs(self.edgeSpeciesRates)
            networkLeakRates = numpy.abs(self.networkLeakRates)
            edgeSpeciesRateRatios = numpy.abs(self.edgeSpeciesRates/charRate)
            networkLeakRateRatios = numpy.abs(self.networkLeakRates/charRate)
            numEdgeReactions = self.numEdgeReactions
            coreReactionRates = self.coreReactionRates
            productIndices = self.productIndices
            reactantIndices = self.reactantIndices
            coreSpeciesConcentrations = self.coreSpeciesConcentrations
            
            # Update the maximum species rate and maximum network leak rate arrays
            for i,index in enumerate(prunableSpeciesIndices):
                if maxEdgeSpeciesRateRatios[i] < edgeSpeciesRateRatios[index]:
                    maxEdgeSpeciesRateRatios[i] = edgeSpeciesRateRatios[index]
            for i,index in enumerate(prunableNetworkIndices):
                if maxNetworkLeakRateRatios[i] < networkLeakRateRatios[index]:
                    maxNetworkLeakRateRatios[i] = networkLeakRateRatios[index]
            
            if charRate == 0 and len(edgeSpeciesRates)>0: #this deals with the case when there is no flux in the system
                maxSpeciesIndex = numpy.argmax(edgeSpeciesRates)
                maxSpecies = edgeSpecies[maxSpeciesIndex]
                maxSpeciesRate = edgeSpeciesRates[maxSpeciesIndex]
                logging.info('At time {0:10.4e} s, species {1} was added to model core to avoid singularity'.format(self.t, maxSpecies))
                invalidObjects.append(maxSpecies)
                break
            
            if useDynamics and not firstTime and self.t >= dynamicsTimeScale:
                #######################################################
                # Calculation of dynamics criterion for edge reactions#
                #######################################################

                totalDivAccumNums = numpy.ones(numEdgeReactions)
                for index in xrange(numEdgeReactions):
                    reactionRate = edgeReactionRates[index]
                    for spcIndex in self.reactantIndices[index+numCoreReactions,:]:
                        if spcIndex != -1 and spcIndex<numCoreSpecies:
                            consumption = coreSpeciesConsumptionRates[spcIndex]
                            if consumption != 0: #if consumption = 0 ignore species
                                totalDivAccumNums[index] *= (reactionRate+consumption)/consumption
                            
                    for spcIndex in self.productIndices[index+numCoreReactions,:]:
                        if spcIndex != -1 and spcIndex<numCoreSpecies:
                            production = coreSpeciesProductionRates[spcIndex]
                            if production != 0: #if production = 0 ignore species
                                totalDivAccumNums[index] *= (reactionRate+production)/production
                    
                totalDivLnAccumNums = numpy.log(totalDivAccumNums)
                
                surfaceSpeciesIndices = self.surfaceSpeciesIndices
                surfaceReactionIndices = self.surfaceReactionIndices
                
                ##########################################################
                # Calculation of dynamics criterion for surface reactions#
                ##########################################################
                
                surfaceTotalDivAccumNums = numpy.ones(len(surfaceReactionIndices))
                
                for i in xrange(len(surfaceReactionIndices)):
                    index = surfaceReactionIndices[i]
                    reactionRate = coreReactionRates[index]
                    for spcIndex in reactantIndices[index,:]:
                        if spcIndex != -1 and spcIndex<numCoreSpecies and not(spcIndex in surfaceSpeciesIndices):
                            consumption = coreSpeciesConsumptionRates[spcIndex]
                            if consumption != 0: #if consumption=0 ignore species
                                if abs(abs(consumption) - abs(reactionRate)) < absoluteTolerance:
                                    surfaceTotalDivAccumNums[i] = numpy.inf
                                elif reactionRate > 0:
                                    surfaceTotalDivAccumNums[i] *= consumption/(consumption-reactionRate)
                                else:
                                    surfaceTotalDivAccumNums[i] *= (consumption-reactionRate)/consumption
                            
                    for spcIndex in productIndices[index,:]:
                        if spcIndex != -1 and spcIndex<numCoreSpecies and not(spcIndex in surfaceSpeciesIndices):
                            production = coreSpeciesProductionRates[spcIndex]
                            if production != 0: #if production = 0 ignore species
                                if abs(abs(production) - abs(reactionRate)) < absoluteTolerance:
                                    surfaceTotalDivAccumNums[i] = numpy.inf
                                elif reactionRate > 0:
                                    surfaceTotalDivAccumNums[i] *= production/(production-reactionRate)
                                else:
                                    surfaceTotalDivAccumNums[i] *= (production-reactionRate)/production

                surfaceTotalDivAccumNums = numpy.log(surfaceTotalDivAccumNums)
                
                ###############################################
                # Move objects from surface to core on-the-fly#
                ###############################################
                
                surfaceObjects = []
                surfaceObjectIndices = []
                
                #Determination of reactions moving from surface to core on-the-fly
                
                for ind,stdan in enumerate(surfaceTotalDivAccumNums): 
                    if stdan > toleranceMoveSurfaceReactionToCore:
                        sind = surfaceReactionIndices[ind]
                        surfaceObjectIndices.append(sind)
                        surfaceObjects.append(coreReactions[sind])
                    
                #Determination of species moving from surface to core on-the-fly

                surfaceSpeciesRateRatios = numpy.zeros(len(surfaceSpeciesIndices))
                surfaceSpeciesProduction = coreSpeciesProductionRates[surfaceSpeciesIndices]
                surfaceSpeciesConsumption = coreSpeciesConsumptionRates[surfaceSpeciesIndices]
                    
                for i in xrange(len(surfaceSpeciesIndices)):
                    RR = max(abs(surfaceSpeciesProduction[i]),abs(surfaceSpeciesConsumption[i]))/charRate
                    surfaceSpeciesRateRatios[i] = RR

                    if RR > toleranceMoveSurfaceSpeciesToCore:
                        sind = surfaceSpeciesIndices[i]
                        surfaceObjectIndices.append(sind)
                        surfaceObjects.append(coreSpecies[sind])
                
                #process objects to be moved from surface to core
                
                if len(surfaceObjects) > 0:
                    schanged = True
                    for ind,obj in enumerate(surfaceObjects):
                        if isinstance(obj,Reaction):
                            logging.info('Moving reaction {0} from surface to core'.format(obj))
                            surfaceReactions.remove(obj)
                        elif isinstance(obj,Species):
                            logging.info('Moving species {0} from surface to core'.format(obj))
                            surfaceSpecies.remove(obj)
                        else:
                            raise ValueError
                    surfaceSpecies,surfaceReactions = self.initialize_surface(coreSpecies,coreReactions,surfaceSpecies,surfaceReactions)
                    logging.info('Surface now has {0} Species and {1} Reactions'.format(len(self.surfaceSpeciesIndices),len(self.surfaceReactionIndices)))
                    
            if filterReactions:
                # Calculate unimolecular and bimolecular thresholds for reaction
                # Set the maximum unimolecular rate to be kB*T/h
                unimolecularThresholdVal = toleranceMoveToCore * charRate / (2.08366122e10 * self.T.value_si)   
                # Set the maximum bimolecular rate to be 1e7 m^3/mol*s, or 1e13 cm^3/mol*s
                bimolecularThresholdVal = toleranceMoveToCore * charRate / 1e7 
                for i in xrange(numCoreSpecies):
                    if not unimolecularThreshold[i]:
                        # Check if core species concentration has gone above threshold for unimolecular reaction
                        if coreSpeciesConcentrations[i] > unimolecularThresholdVal:
                            unimolecularThreshold[i] = True
                for i in xrange(numCoreSpecies):
                    for j in xrange(i, numCoreSpecies):
                        if not bimolecularThreshold[i,j]:
                            if coreSpeciesConcentrations[i]*coreSpeciesConcentrations[j] > bimolecularThresholdVal:
                                bimolecularThreshold[i,j] = True
            
            
            ###############################################################################
            # Movement from edge to core or surface processing and interrupt determination#
            ###############################################################################
            
            newObjectInds = []
            newObjects = []
            newObjectVals = []
            
            tempNewObjects = []
            tempNewObjectInds = []
            tempNewObjectVals = []
            
            newSurfaceRxnInds = []
            interrupt = False
            
            #movement of species to core based on rate ratios
            
            if not ignoreOverallFluxCriterion:
                for ind,obj in enumerate(edgeSpecies):
                    RR = edgeSpeciesRateRatios[ind]
                    if RR > toleranceMoveToCore:
                        if not(obj in newObjects or obj in invalidObjects):
                            tempNewObjects.append(edgeSpecies[ind])
                            tempNewObjectInds.append(ind)
                            tempNewObjectVals.append(RR)
                    if RR > toleranceInterruptSimulation:
                        logging.info('At time {0:10.4e} s, species {1} at {2} exceeded the minimum rate for simulation interruption of {3}'.format(self.t, obj, RR, toleranceInterruptSimulation))
                        interrupt = True
                
                
                sortedInds = numpy.argsort(numpy.array(tempNewObjectVals)).tolist()[::-1]
                
                newObjects.extend([tempNewObjects[q] for q in sortedInds])
                newObjectInds.extend([tempNewObjectInds[q] for q in sortedInds])
                newObjectVals.extend([tempNewObjectVals[q] for q in sortedInds])
                
                tempNewObjects = []
                tempNewObjectInds = []
                tempNewObjectVals = []
                              
            if useDynamics and not firstTime and self.t >= dynamicsTimeScale:     
                #movement of reactions to core/surface based on dynamics number  
                validLayeringIndices = self.validLayeringIndices
                tempSurfaceObjects = []
                
                for ind,obj in enumerate(edgeReactions):
                    dlnaccum = totalDivLnAccumNums[ind]
                    if dlnaccum > toleranceMoveEdgeReactionToCore:
                        if not(obj in newObjects or obj in invalidObjects):
                            tempNewObjects.append(edgeReactions[ind])
                            tempNewObjectInds.append(ind)
                            tempNewObjectVals.append(dlnaccum)
                    elif dlnaccum > toleranceMoveEdgeReactionToSurface and ind in validLayeringIndices:
                        if not(obj in newObjects or obj in invalidObjects):
                            tempNewObjects.append(edgeReactions[ind])
                            tempNewObjectInds.append(ind)
                            tempNewObjectVals.append(dlnaccum)
                            tempSurfaceObjects.append(edgeReactions[ind])
                    if dlnaccum > toleranceMoveEdgeReactionToCoreInterrupt:
                        logging.info('At time {0:10.4e} s, Reaction {1} at {2} exceeded the minimum difference in total log(accumulation number) for simulation interruption of {3}'.format(self.t, obj,dlnaccum,toleranceMoveEdgeReactionToCoreInterrupt))
                        interrupt = True
                
                sortedInds = numpy.argsort(numpy.array(tempNewObjectVals)).tolist()[::-1]
                
                newObjects.extend([tempNewObjects[q] for q in sortedInds])
                newObjectInds.extend([tempNewObjectInds[q] for q in sortedInds])
                newObjectVals.extend([tempNewObjectVals[q] for q in sortedInds])
                
                newSurfaceRxnInds = [newObjects.index(obj) for obj in tempSurfaceObjects]
                
                tempNewObjects = []
                tempNewObjectInds = []
                tempNewObjectVals = []
                
            #Determination of pdepNetworks in need of exploring
            
            if pdepNetworks:
                for ind,obj in enumerate(pdepNetworks):
                    LR = networkLeakRateRatios[ind]
                    if LR > toleranceMoveToCore:
                        if not(obj in newObjects or obj in invalidObjects):
                            tempNewObjects.append(pdepNetworks[ind])
                            tempNewObjectInds.append(ind)
                            tempNewObjectVals.append(LR)
                    if LR > toleranceInterruptSimulation:
                        logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} at {2} exceeded the minimum rate for simulation interruption of {3}'.format(self.t, obj.index,LR,toleranceInterruptSimulation))
                        interrupt = True
                
                sortedInds = numpy.argsort(numpy.array(tempNewObjectVals)).tolist()[::-1]
                
                newObjects.extend([tempNewObjects[q] for q in sortedInds])
                newObjectInds.extend([tempNewObjectInds[q] for q in sortedInds])
                newObjectVals.extend([tempNewObjectVals[q] for q in sortedInds])
                
                tempNewObjects = []
                tempNewObjectInds = []
                tempNewObjectVals = []
            
            ###########################
            #Overall Object Processing#
            ###########################
            
            #remove excess objects
            if len(invalidObjects) + len(newObjects) > maxNumObjsPerIter:
                logging.info('Exceeded max number of objects...removing excess objects')
                num = maxNumObjsPerIter - len(invalidObjects)
                newObjects = newObjects[:num]
                newObjectInds = newObjectInds[:num]
                newObjectVals = newObjectVals[:num]
            
            if terminateAtMaxObjects and len(invalidObjects)+len(newObjects) >= maxNumObjsPerIter:
                logging.info('Reached max number of objects...preparing to terminate')
                interrupt = True
            
            if newObjects != []:
                for i,obj in enumerate(newObjects): #log everything and add reactions to surface
                    val = newObjectVals[i]
                    ind = newObjectInds[i]
                    if isinstance(obj,Species):
                        logging.info('At time {0:10.4e} s, species {1} at rate ratio {2} exceeded the minimum rate for moving to model core of {3}'.format(self.t, obj,val,toleranceMoveToCore))
                    elif isinstance(obj,Reaction):
                        if i in newSurfaceRxnInds:
                            logging.info('At time {0:10.4e} s, Reaction {1} at {2} exceeded the minimum difference in total log(accumulation number) for moving to model surface of {3}'.format(self.t, obj, val,toleranceMoveEdgeReactionToSurface))
                            newSurfaceReactions.append(obj)
                            newSurfaceReactionInds.append(ind)
                        else:
                            logging.info('At time {0:10.4e} s, Reaction {1} at {2} exceeded the minimum difference in total log(accumulation number) for moving to model core of {3}'.format(self.t, obj, val,toleranceMoveEdgeReactionToCore))
                    else: 
                        logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} at {2} exceeded the minimum rate for exploring of {3}'.format(self.t, obj.index, val,toleranceMoveToCore))
    
                
                invalidObjects += newObjects
                
            if schanged: #reinitialize surface
                surfaceSpecies,surfaceReactions = self.initialize_surface(coreSpecies,coreReactions,surfaceSpecies,surfaceReactions)
                schanged = False

            if firstTime: #turn off firstTime
                firstTime = False
                
            if interrupt: #breaks while loop terminating iterations
                logging.info('terminating simulation due to interrupt...')
                break
            
            # Finish simulation if any of the termination criteria are satisfied
            for term in self.termination:
                if isinstance(term, TerminationTime):
                    if self.t > term.time.value_si:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination time.'.format(term.time.value_si))
                        self.logConversions(speciesIndex, y0)
                        break
                elif isinstance(term, TerminationConversion):
                    index = speciesIndex[term.species]
                    conversion = 1-(y_coreSpecies[index] / y0[index])
                    if 1 - (y_coreSpecies[index] / y0[index]) > term.conversion:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination conversion: {1:f} of {2}'.format(self.t,term.conversion,term.species))
                        self.logConversions(speciesIndex, y0)
                        break
                elif isinstance(term, TerminationRateRatio):
                    if maxCharRate != 0.0 and charRate/maxCharRate < term.ratio:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination RateRatio: {1}'.format(self.t,charRate/maxCharRate))
                        self.logConversions(speciesIndex, y0)
                        
            # Increment destination step time if necessary
            if self.t >= 0.9999 * stepTime:
                stepTime *= 10.0
        
        #change surface species and reactions based on what will be added to the surface
        surfaceSpecies,surfaceReactions=self.addReactionsToSurface(newSurfaceReactions,newSurfaceReactionInds,surfaceSpecies,surfaceReactions,edgeSpecies)
        
        # notify reaction system listeners
        self.notify()

        if sensitivity:   
            for i in xrange(len(self.sensitiveSpecies)):
                with open(sensWorksheet[i], 'wb') as outfile:
                    worksheet = csv.writer(outfile)
                    reactionsAboveThreshold = []
                    for j in xrange(numCoreReactions + numCoreSpecies):
                        for k in xrange(len(time_array)):
                            if abs(normSens_array[i][k][j]) > self.sensitivityThreshold:
                                reactionsAboveThreshold.append(j)
                                break
                    species_name = getSpeciesIdentifier(self.sensitiveSpecies[i])
                    headers = ['Time (s)']
                    headers.extend(['dln[{0}]/dln[k{1}]: {2}'.format(species_name, j+1, coreReactions[j].toChemkin(kinetics=False)) if j < numCoreReactions 
                                    else 'dln[{0}]/dG[{1}]'.format(species_name, getSpeciesIdentifier(coreSpecies[j-numCoreReactions])) for j in reactionsAboveThreshold])
                    worksheet.writerow(headers)               
                
                    for k in xrange(len(time_array)):
                        row = [time_array[k]]
                        row.extend([normSens_array[i][k][j] for j in reactionsAboveThreshold])       
                        worksheet.writerow(row)  
        
        self.maxEdgeSpeciesRateRatios = maxEdgeSpeciesRateRatios
        self.maxNetworkLeakRateRatios = maxNetworkLeakRateRatios
        
        self.unimolecularThreshold = unimolecularThreshold
        self.bimolecularThreshold = bimolecularThreshold

        # Return the invalid object (if the simulation was invalid) or None
        # (if the simulation was valid)
        return terminated, False, invalidObjects, surfaceSpecies, surfaceReactions, self.t, conversion

    cpdef logRates(self, double charRate, object species, double speciesRate, double maxDifLnAccumNum, object network, double networkRate):
        """
        Log information about the current maximum species and network rates.
        """
        logging.info('    Characteristic rate: {0:10.4e} mol/m^3*s'.format(charRate))
        logging.info('    Dynamics number:  {0:10.4e}'.format(maxDifLnAccumNum))
        if charRate == 0.0:
            logging.info('    {0} rate: {1:10.4e} mol/m^3*s'.format(species, speciesRate))
            if network is not None:
                logging.info('    PDepNetwork #{0:d} leak rate: {1:10.4e} mol/m^3*s'.format(network.index, networkRate))
        else:
            logging.info('    {0} rate: {1:10.4e} mol/m^3*s ({2:.4g})'.format(species, speciesRate, speciesRate / charRate))
            if network is not None:
                logging.info('    PDepNetwork #{0:d} leak rate: {1:10.4e} mol/m^3*s ({2:.4g})'.format(network.index, networkRate, networkRate / charRate))
                
    cpdef logConversions(self, speciesIndex, y0):
        """
        Log information about the current conversion values.
        """
        for term in self.termination:
            if isinstance(term, TerminationConversion):
                index = speciesIndex[term.species]
                X = 1 - (self.y[index] / y0[index])
                logging.info('    {0} conversion: {1:<10.4g}'.format(term.species, X))

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
        
        kf = self.kf
        kr = self.kb    
        
        numCoreReactions = len(self.coreReactionRates)
        numCoreSpecies = len(self.coreSpeciesConcentrations)      
        
        # Use stored volume, since this function is only called from residual function. 
        RT_inverse = 1/(constants.R * self.T.value_si)
        V = self.V

        C = self.coreSpeciesConcentrations

        rateDeriv = numpy.zeros((numCoreSpecies,numCoreReactions+numCoreSpecies), numpy.float64)
        
        for j in xrange(numCoreReactions):
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
        
################################################################################

class TerminationTime:
    """
    Represent a time at which the simulation should be terminated. This class
    has one attribute: the termination `time` in seconds.
    """

    def __init__(self, time=(0.0,'s')):
        self.time = Quantity(time)

################################################################################

class TerminationConversion:
    """
    Represent a conversion at which the simulation should be terminated. This
    class has two attributes: the `species` to monitor and the fractional
    `conversion` at which to terminate.
    """

    def __init__(self, spec=None, conv=0.0):
        self.species = spec
        self.conversion = conv
        
class TerminationRateRatio:
    """
    Represent a fraction of the maximum characteristic rate of the simulation
    at which the simulation should be terminated.  This class has one attribute
    the ratio between the current and maximum characteristic rates at which
    to terminate
    """
    
    def __init__(self, ratio=0.01):
        self.ratio = ratio
        
