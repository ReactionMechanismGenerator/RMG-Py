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
Contains the :class:`ReactionSystemr` class, a base class for all RMG reaction
systems.
"""

cdef extern from "math.h":
    double sqrt(double)

import numpy
cimport numpy
from pydas cimport DASSL

import logging

################################################################################

cdef class ReactionSystem(DASSL):
    """
    A base class for all RMG reaction systems.
    """

    def __init__(self):
        DASSL.__init__(self)
        # The reaction and species rates at the current time (in mol/m^3*s)
        self.coreSpeciesRates = None
        self.coreReactionRates = None
        self.edgeSpeciesRates = None
        self.edgeReactionRates = None
        self.networkLeakRates = None

    cpdef simulate(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions,
        double toleranceKeepInEdge, double toleranceMoveToCore, double toleranceInterruptSimulation,
        list termination, list pdepNetworks=None):
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

        cdef dict speciesIndex
        cdef int index, maxSpeciesIndex, maxNetworkIndex
        cdef double stepTime, charRate, maxSpeciesRate, maxNetworkRate
        cdef numpy.ndarray[numpy.float64_t, ndim=1] y0
        cdef bint terminated
        cdef object maxSpecies, maxNetwork
        
        pdepNetworks = pdepNetworks or []

        speciesIndex = {}
        for index, spec in enumerate(coreSpecies):
            speciesIndex[spec] = index
        
        self.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions, pdepNetworks)

        invalidObject = None
        terminated = False
        maxSpeciesIndex = -1
        maxSpecies = None
        maxSpeciesRate = 0.0
        maxNetworkIndex = -1
        maxNetwork = None
        maxNetworkRate = 0.0

        # Copy the initial conditions to use in evaluating conversions
        y0 = self.y.copy()

        stepTime = 1e-12
        while not terminated:
            # Integrate forward in time by one time step
            self.step(stepTime)

            # Get the characteristic flux
            charRate = sqrt(numpy.sum(self.coreSpeciesRates * self.coreSpeciesRates))

            # Get the edge species with the highest flux
            maxSpeciesIndex = numpy.argmax(self.edgeSpeciesRates)
            maxSpecies = edgeSpecies[maxSpeciesIndex]
            maxSpeciesRate = self.edgeSpeciesRates[maxSpeciesIndex]
            if pdepNetworks:
                maxNetworkIndex = numpy.argmax(self.networkLeakRates)
                maxNetwork = pdepNetworks[maxNetworkIndex]
                maxNetworkRate = self.networkLeakRates[maxNetworkIndex]

            # Interrupt simulation if that flux exceeds the characteristic rate times a tolerance
            if maxSpeciesRate > toleranceMoveToCore * charRate and not invalidObject:
                logging.info('At time %10.4e s, species %s exceeded the minimum rate for moving to model core' % (self.t, maxSpecies))
                self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                self.logConversions(termination, speciesIndex, y0)
                invalidObject = maxSpecies
            if maxSpeciesRate > toleranceInterruptSimulation * charRate:
                logging.info('At time %10.4e s, species %s exceeded the minimum rate for simulation interruption' % (self.t, maxSpecies))
                self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                self.logConversions(termination, speciesIndex, y0)
                break

            # If pressure dependence, also check the network leak fluxes
            if pdepNetworks:
                if maxNetworkRate > toleranceMoveToCore * charRate and not invalidObject:
                    logging.info('At time %10.4e s, PDepNetwork #%i exceeded the minimum rate for exploring' % (self.t, maxNetwork.index))
                    self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                    self.logConversions(termination, speciesIndex, y0)
                    invalidObject = maxNetwork
                if maxNetworkRate > toleranceInterruptSimulation * charRate:
                    logging.info('At time %10.4e s, PDepNetwork #%i exceeded the minimum rate for simulation interruption' % (self.t, maxNetwork.index))
                    self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                    self.logConversions(termination, speciesIndex, y0)
                    break

            # Finish simulation if any of the termination criteria are satisfied
            for term in termination:
                if isinstance(term, TerminationTime):
                    if self.t > term.time:
                        terminated = True
                        break
                elif isinstance(term, TerminationConversion):
                    index = speciesIndex[term.species]
                    if (y0[index] - self.y[index]) / y0[index] > term.conversion:
                        terminated = True
                        break

            # Increment destination step time if necessary
            if self.t >= 0.9999 * stepTime:
                stepTime *= 10.0
                
        # Return the invalid object (if the simulation was invalid) or None
        # (if the simulation was valid)
        return invalidObject

    cpdef logRates(self, double charRate, object species, double speciesRate, object network, double networkRate):
        """
        Log information about the current maximum species and network rates.
        """
        logging.info('    Characteristic rate: %10.4e mol/m^3*s' % (charRate))
        if charRate == 0.0:
            logging.info('    %s rate: %10.4e mol/m^3*s' % (species, speciesRate))
            if network is not None:
                logging.info('    PDepNetwork #%i leak rate: %10.4e mol/m^3*s' % (network.index, networkRate))
        else:
            logging.info('    %s rate: %10.4e mol/m^3*s (%.4g)' % (species, speciesRate, speciesRate / charRate))
            if network is not None:
                logging.info('    PDepNetwork #%i leak rate: %10.4e mol/m^3*s (%.4g)' % (network.index, networkRate, networkRate / charRate))

    cpdef logConversions(self, termination, speciesIndex, y0):
        """
        Log information about the current conversion values.
        """
        for term in termination:
            if isinstance(term, TerminationConversion):
                index = speciesIndex[term.species]
                X = (y0[index] - self.y[index]) / y0[index]
                logging.info('    %s conversion: %-10.4g' % (term.species, X))

################################################################################

class TerminationTime:
    """
    Represent a time at which the simulation should be terminated. This class
    has one attribute: the termination `time` in seconds.
    """

    def __init__(self, time=0.0):
        self.time = time

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
