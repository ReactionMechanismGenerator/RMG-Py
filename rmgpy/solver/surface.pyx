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

import itertools
import logging

from base cimport ReactionSystem
cimport cython

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

cdef class SurfaceReactor(ReactionSystem):
    """
    A reaction system consisting of a heterogeneous, isothermal, constant volume batch
    reactor. 
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity initialP
    cdef public double V
    cdef public bint constantVolume

    cdef public dict initialGasMoleFractions
    cdef public dict initialSurfaceCoverages
    cdef public ScalarQuantity surfaceVolumeRatio
    cdef public ScalarQuantity surfaceSiteDensity
    cdef public numpy.ndarray surfaceReactions
    cdef public numpy.ndarray surfaceSpecies

    def __init__(self, 
                 T, 
                 initialP,
                 initialGasMoleFractions,
                 initialSurfaceCoverages,
                 surfaceVolumeRatio,
                 surfaceSiteDensity,
                 termination,
                 sensitiveSpecies=None,
                 sensitivityThreshold=1e-3):
        ReactionSystem.__init__(self,
                                termination, 
                                sensitiveSpecies, 
                                sensitivityThreshold)

        self.T = Quantity(T)
        self.initialP = Quantity(initialP)
        self.initialGasMoleFractions = initialGasMoleFractions
        self.initialSurfaceCoverages = initialSurfaceCoverages
        self.surfaceVolumeRatio = Quantity(surfaceVolumeRatio)
        self.surfaceSiteDensity = Quantity(surfaceSiteDensity)
        self.V = 0 # will be set from ideal gas law in initializeModel
        self.constantVolume = True
        
    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initialGasMoleFractions and initialSurfaceCoverages dictionaries
        from species names into species objects,
        using the given dictionary of species.
        """
        initialGasMoleFractions = {}
        for label, moleFrac in self.initialGasMoleFractions.iteritems():
            initialGasMoleFractions[speciesDict[label]] = moleFrac
        self.initialGasMoleFractions = initialGasMoleFractions
        initialSurfaceCoverages = {}
        for label, surfaceCoverage in self.initialSurfaceCoverages.iteritems():
            initialSurfaceCoverages[speciesDict[label]] = surfaceCoverage
        self.initialSurfaceCoverages = initialSurfaceCoverages

    cpdef initializeModel(self,
                          list coreSpecies,
                          list coreReactions,
                          list edgeSpecies,
                          list edgeReactions,
                          list pdepNetworks=None,
                          atol=1e-16,
                          rtol=1e-8,
                          sensitivity=False,
                          sens_atol=1e-6,
                          sens_rtol=1e-4,
                          filterReactions=False):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions, pdepNetworks, atol, rtol, sensitivity, sens_atol, sens_rtol)

        cdef numpy.ndarray[numpy.int_t, ndim=1] speciesOnSurface, reactionsOnSurface
        cdef int index
        #: 1 if it's on a surface, 0 if it's in the gas phase
        reactionsOnSurface = numpy.zeros((self.numCoreReactions + self.numEdgeReactions), numpy.int)
        speciesOnSurface = numpy.zeros((self.numCoreSpecies), numpy.int)
        for spec, index in self.speciesIndex.iteritems():
            if index >= self.numCoreSpecies:
                continue
            if spec.containsSurfaceSite():
                surfaceSpecies[index] = 1
        for rxn, index in self.reactionIndex.iteritems():
            if rxn.isSurfaceReaction():
                surfaceReactions[index] = 1
        self.speciesOnSurface = speciesOnSurface
        self.reactionsOnSurface = reactionsOnSurface
        
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
        Populates the kf, kb and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reacion system.
        """
        
        cdef double P, surfaceVolumeRatioSI
        
        surfaceVolumeRatioSI = self.surfaceVolumeRatio.value_si
        # ToDo: Pressure should come from ideal gas law?
        P = self.initialP.value_si

        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            
            # ToDo: getRateCoefficient should also depend on surface coverages vector
            assert not rxn.kinetics.isPressureDependent(), "Pressure may be varying."
            if rxn.isSurfaceReaction():
                """
                Be careful! From here on kf and kb will now be in Volume units,
                even for surface reactions (which you may expect to be in Area units).
                This is to avoid repeatedly multiplying a bunch of things inside every 
                loop of the ODE solver.
                """
                self.kf[j] = (surfaceVolumeRatioSI *
                              rxn.getSurfaceRateCoefficient(self.T.value_si,
                                                            self.surfaceSiteDensity.value_si
                                                            ))
            else:
                self.kf[j] = rxn.getRateCoefficient(self.T.value_si, P)
            if rxn.reversible:
                # ToDo: getEquilibriumConstant should be coverage dependent
                self.Keq[j] = rxn.getEquilibriumConstant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]

    def log_initial_conditions(self, number=None):
        """
        Log to the console some information about this reaction system.
        
        Should correspond to the calculations done in set_initial_conditions.
        """
        logging.info("\nSurface reaction system {}".format(number if number is not None else ""))
        logging.info("Gas phase mole fractions:")
        totalGasMoles = 0
        for spec, moleFrac in self.initialGasMoleFractions.iteritems():
            logging.info("  {0:20s} {1:.5g}".format(spec, moleFrac))
            totalGasMoles += moleFrac
        logging.info("Total gas phase:          {:.3g} moles".format(totalGasMoles))
        logging.info("Pressure:                 {:.3g} Pa".format(self.initialP.value_si))
        logging.info("Temperature:              {} K".format(self.T.value_si))
        V = constants.R * self.T.value_si * totalGasMoles / self.initialP.value_si
        logging.info("Reactor volume:           {:.3g} m3".format(V))
        surfaceVolumeRatio_si = self.surfaceVolumeRatio.value_si # 1/m
        logging.info("Surface/volume ratio:     {:.3g} m2/m3".format(surfaceVolumeRatio_si))
        logging.info("Surface site density:     {:.3g} mol/m2".format(self.surfaceSiteDensity.value_si))
        totalSurfaceSites = V * surfaceVolumeRatio_si * self.surfaceSiteDensity.value_si # total surface sites in reactor
        logging.info("Surface sites in reactor: {:.3g} moles".format(totalSurfaceSites))
        logging.info("Initial surface coverages (and amounts):")
        for spec, coverage in self.initialSurfaceCoverages.iteritems():
            logging.info("  {:18s} {:.5g} = {:.5g} moles".format(spec, coverage, totalSurfaceSites*coverage ))
        
        
    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the 
        current reactor model.

        The volume is set to the value in m3 required to contain 
        one mole total of gas phase core species at start.

        The total surface sites are calculated from surfaceVolumeRatio and surfaceSiteDensity
        allowing initialSurfaceCoverages to determine the number of moles of surface species.
        The number of moles of gas phase species is taken from initialGasMoleFractions.
        
        The coreSpeciesConcentrations array is then determined, in mol/m3 for gas phase
        and mol/m2 for surface species.

        The initial number of moles of a species j in the reactor is computed and stored in the
        y0 instance attribute.

        """
        ### 
        ### WARNING -- When updating this method, please be sure 
        ###            to also update the log_initial_conditions above
        ###            which unfortunately is maintained separately.
        ### 
        cdef double V, P, surfaceVolumeRatio_si
        
        ReactionSystem.set_initial_conditions(self)
        # self.y0 tracks number of moles of each core species

        # add only the gas phase species first
        self.y0 *= 0.
        for spec, moleFrac in self.initialGasMoleFractions.iteritems():
            i = self.get_species_index(spec)
            self.y0[i] = moleFrac # moles in reactor

        # Use ideal gas law to compute reactor volume
        V = constants.R * self.T.value_si * numpy.sum(self.y0[:self.numCoreSpecies]) / self.initialP.value_si
        self.V = V # volume in m^3  (assume reactor volume is gas phase volume, i.e catalyst takes no space)

        surfaceVolumeRatio_si = self.surfaceVolumeRatio.value_si # 1/m
        totalSurfaceSites = V * surfaceVolumeRatio_si * self.surfaceSiteDensity.value_si # total surface sites in reactor

        for spec, coverage in self.initialSurfaceCoverages.iteritems():
            i = self.get_species_index(spec)
            self.y0[i] = totalSurfaceSites * coverage # moles in reactor
        
        for j, isSurfaceSpecies in enumerate(self.surfaceSpecies): # should only go up to core species
            if isSurfaceSpecies:
                self.coreSpeciesConcentrations[j] = self.y0[j] / V / surfaceVolumeRatio_si # moles per m2 of surface
            else:
                self.coreSpeciesConcentrations[j] = self.y0[j] / V # moles per m3 of gas
        
    def compute_network_variables(self, pdepNetworks=None):
        # ToDo: this should allow pressure to vary?
        # for now, just call the base class version.
        ReactionSystem.compute_network_variables(self, pdepNetworks)

    @cython.boundscheck(False)
    def residual(self,
                 double t,
                 numpy.ndarray[numpy.float64_t, ndim=1] N,
                 numpy.ndarray[numpy.float64_t, ndim=1] dNdt,
                 numpy.ndarray[numpy.float64_t, ndim=1] senpar = numpy.zeros(1, numpy.float64)
                 ):

        """
        Return the residual function for the governing DAE system for the
        simple reaction system.
        """
        cdef numpy.ndarray[numpy.int_t, ndim=2] ir, ip, inet
        cdef numpy.ndarray[numpy.int_t, ndim=1] surfaceReactions, surfaceSpecies
        cdef numpy.ndarray[numpy.float64_t, ndim=1] res, kf, kr, knet, delta, equilibriumConstants
        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks
        cdef int i, j, z, first, second, third
        cdef double k, V, reactionRate, surfaceVolumeRatio_si
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesConcentrations, coreSpeciesRates, coreReactionRates, edgeSpeciesRates, edgeReactionRates, networkLeakRates
        cdef numpy.ndarray[numpy.float64_t, ndim=1] C
        cdef numpy.ndarray[numpy.float64_t, ndim=2] jacobian, dgdk

        ir = self.reactantIndices
        ip = self.productIndices
        equilibriumConstants = self.Keq

        kf = self.kf # are already 'per m3 of reactor' even for surface reactions
        kr = self.kb # are already 'per m3 of reactor' even for surface reactions
        
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
        
        surfaceReactions = self.surfaceReactions
        surfaceSpecies = self.surfaceSpecies
        surfaceVolumeRatio_si = self.surfaceVolumeRatio.value_si

        C = numpy.zeros_like(self.coreSpeciesConcentrations)
        V =  self.V # constant volume reactor

        for j in xrange(numCoreSpecies):
            if surfaceSpecies[j]:
                C[j] = (N[j] / V) / surfaceVolumeRatio_si
            else:
                C[j] = N[j] / V
            #: surface species are in mol/m2, gas phase are in mol/m3
            coreSpeciesConcentrations[j] = C[j]
        
        for j in xrange(ir.shape[0]):
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

            "reactionRate is now in mol/m3/s"
            # Set the reaction and species rates
            if j < numCoreReactions:
                # The reaction is a core reaction
                coreReactionRates[j] = reactionRate

                # Add/subtract the total reaction rate from each species rate
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
        self.coreReactionRates = coreReactionRates
        self.edgeSpeciesRates = edgeSpeciesRates
        self.edgeReactionRates = edgeReactionRates
        self.networkLeakRates = networkLeakRates

        res = coreSpeciesRates * V 
        # mol/s
        
        if self.sensitivity and False:
            delta = numpy.zeros(len(N), numpy.float64)
            delta[:numCoreSpecies] = res
            if self.jacobianMatrix is None:
                jacobian = self.jacobian(t,N,dNdt,0,senpar)
            else:
                jacobian = self.jacobianMatrix
            dgdk = ReactionSystem.computeRateDerivative(self)
            for j in xrange(numCoreReactions+numCoreSpecies):
                for i in xrange(numCoreSpecies):
                    for z in xrange(numCoreSpecies):
                        delta[(j+1)*numCoreSpecies + i] += jacobian[i,z]*N[(j+1)*numCoreSpecies + z] 
                    delta[(j+1)*numCoreSpecies + i] += dgdk[i,j]
        else:
            delta = res
        delta = delta - dNdt
        
        # Return DELTA, IRES.  IRES is set to 1 in order to tell DASPK to evaluate the sensitivity residuals
        return delta, 1