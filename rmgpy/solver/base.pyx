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

import cython
import logging

from rmgpy.quantity cimport Quantity

################################################################################

cdef class ReactionSystem(DASSL):
    """
    A base class for all RMG reaction systems.
    """

    def __init__(self, termination=None):
        DASSL.__init__(self)
        # The reaction and species rates at the current time (in mol/m^3*s)
        self.coreSpeciesConcentrations = None
        self.coreSpeciesRates = None
        self.coreReactionRates = None
        self.edgeSpeciesRates = None
        self.edgeReactionRates = None
        self.networkLeakRates = None
        self.maxCoreSpeciesRates = None
        self.maxEdgeSpeciesRates = None
        self.maxNetworkLeakRates = None
        self.termination = termination or []
    
    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions, list pdepNetworks=None, atol=1e-16, rtol=1e-8):
        """
        Initialize a simulation of the reaction system using the provided
        kinetic model. You will probably want to create your own version of this
        method in the derived class; don't forget to also call the base class
        version, too.
        """
        cdef int numCoreSpecies, numCoreReactions, numEdgeSpecies, numEdgeReactions, numPdepNetworks

        pdepNetworks = pdepNetworks or []

        numCoreSpecies = len(coreSpecies)
        numCoreReactions = len(coreReactions)
        numEdgeSpecies = len(edgeSpecies)
        numEdgeReactions = len(edgeReactions)
        numPdepNetworks = len(pdepNetworks)

        self.coreSpeciesConcentrations = numpy.zeros((numCoreSpecies), numpy.float64)
        self.coreReactionRates = numpy.zeros((numCoreReactions), numpy.float64)
        self.edgeReactionRates = numpy.zeros((numEdgeReactions), numpy.float64)
        self.coreSpeciesRates = numpy.zeros((numCoreSpecies), numpy.float64)
        self.edgeSpeciesRates = numpy.zeros((numEdgeSpecies), numpy.float64)
        self.networkLeakRates = numpy.zeros((numPdepNetworks), numpy.float64)
        self.maxCoreSpeciesRates = numpy.zeros((numCoreSpecies), numpy.float64)
        self.maxEdgeSpeciesRates = numpy.zeros((numEdgeSpecies), numpy.float64)
        self.maxNetworkLeakRates = numpy.zeros((numPdepNetworks), numpy.float64)

    
    cpdef writeWorksheetHeader(self, worksheet):
        """
        Write some descriptive information about the reaction system to the
        first two rows of the given `worksheet`.
        """
        import xlwt
        style0 = xlwt.easyxf('font: bold on')
        worksheet.write(0, 0, 'Reaction System', style0)

    @cython.boundscheck(False)
    cpdef simulate(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions,
        double toleranceKeepInEdge, double toleranceMoveToCore, double toleranceInterruptSimulation,
        list pdepNetworks=None, worksheet=None, absoluteTolerance=1e-16, relativeTolerance=1e-8):
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
        cdef int numCoreSpecies, numEdgeSpecies, numPdepNetworks
        cdef double stepTime, charRate, maxSpeciesRate, maxNetworkRate
        cdef numpy.ndarray[numpy.float64_t, ndim=1] y0
        cdef numpy.ndarray[numpy.float64_t, ndim=1] coreSpeciesRates, edgeSpeciesRates, networkLeakRates
        cdef numpy.ndarray[numpy.float64_t, ndim=1] maxCoreSpeciesRates, maxEdgeSpeciesRates, maxNetworkLeakRates
        cdef bint terminated
        cdef object maxSpecies, maxNetwork
        cdef int iteration, i
        
        pdepNetworks = pdepNetworks or []

        numCoreSpecies = len(coreSpecies)
        numEdgeSpecies = len(edgeSpecies)
        numPdepNetworks = len(pdepNetworks)

        speciesIndex = {}
        for index, spec in enumerate(coreSpecies):
            speciesIndex[spec] = index
        
        self.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions, pdepNetworks, absoluteTolerance, relativeTolerance)

        invalidObject = None
        terminated = False
        maxSpeciesIndex = -1
        maxSpecies = None
        maxSpeciesRate = 0.0
        maxNetworkIndex = -1
        maxNetwork = None
        maxNetworkRate = 0.0
        iteration = 0

        maxCoreSpeciesRates = self.maxCoreSpeciesRates
        maxEdgeSpeciesRates = self.maxEdgeSpeciesRates
        maxNetworkLeakRates = self.maxNetworkLeakRates
        
        # Copy the initial conditions to use in evaluating conversions
        y0 = self.y.copy()

        if worksheet:
            import xlwt
            self.writeWorksheetHeader(worksheet)
            style0 = xlwt.easyxf('font: bold on')
            style1 = xlwt.easyxf(num_format_str='0.000E+00')
            worksheet.write(3, 0, 'Time (s)', style0)
            worksheet.write(3, 1, 'Concentrations (mol/m^3)', style0)
            for i in range(numCoreSpecies):
                worksheet.write(4, i+1, str(coreSpecies[i]), style0)
            
        stepTime = 1e-12
        while not terminated:
            # Integrate forward in time by one time step
            self.step(stepTime)

            iteration += 1

            if worksheet:
                worksheet.write(iteration+4, 0, self.t, style1)
                for i in range(numCoreSpecies):
                    worksheet.write(iteration+4, i+1, self.coreSpeciesConcentrations[i], style1)

            # Get the characteristic flux
            charRate = sqrt(numpy.sum(self.coreSpeciesRates * self.coreSpeciesRates))

            coreSpeciesRates = numpy.abs(self.coreSpeciesRates)
            edgeSpeciesRates = numpy.abs(self.edgeSpeciesRates)
            networkLeakRates = numpy.abs(self.networkLeakRates)

            # Update the maximum species rate and maximum network leak rate arrays
            for index in range(numCoreSpecies):
                if maxCoreSpeciesRates[index] < coreSpeciesRates[index]:
                    maxCoreSpeciesRates[index] = coreSpeciesRates[index]
            for index in range(numEdgeSpecies):
                if maxEdgeSpeciesRates[index] < edgeSpeciesRates[index]:
                    maxEdgeSpeciesRates[index] = edgeSpeciesRates[index]
            for index in range(numPdepNetworks):
                if maxNetworkLeakRates[index] < networkLeakRates[index]:
                    maxNetworkLeakRates[index] = networkLeakRates[index]

            # Get the edge species with the highest flux
            if numEdgeSpecies > 0:
                maxSpeciesIndex = numpy.argmax(edgeSpeciesRates)
                maxSpecies = edgeSpecies[maxSpeciesIndex]
                maxSpeciesRate = edgeSpeciesRates[maxSpeciesIndex]
            else:
                maxSpeciesIndex = -1
                maxSpecies = None
                maxSpeciesRate = 0.0
            if pdepNetworks:
                maxNetworkIndex = numpy.argmax(networkLeakRates)
                maxNetwork = pdepNetworks[maxNetworkIndex]
                maxNetworkRate = networkLeakRates[maxNetworkIndex]

            # Interrupt simulation if that flux exceeds the characteristic rate times a tolerance
            if maxSpeciesRate > toleranceMoveToCore * charRate and not invalidObject:
                logging.info('At time {0:10.4e} s, species {1} exceeded the minimum rate for moving to model core'.format(self.t, maxSpecies))
                self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                self.logConversions(speciesIndex, y0)
                invalidObject = maxSpecies
            if maxSpeciesRate > toleranceInterruptSimulation * charRate:
                logging.info('At time {0:10.4e} s, species {1} exceeded the minimum rate for simulation interruption'.format(self.t, maxSpecies))
                self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                self.logConversions(speciesIndex, y0)
                break

            # If pressure dependence, also check the network leak fluxes
            if pdepNetworks:
                if maxNetworkRate > toleranceMoveToCore * charRate and not invalidObject:
                    logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} exceeded the minimum rate for exploring'.format(self.t, maxNetwork.index))
                    self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                    self.logConversions(speciesIndex, y0)
                    invalidObject = maxNetwork
                if maxNetworkRate > toleranceInterruptSimulation * charRate:
                    logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} exceeded the minimum rate for simulation interruption'.format(self.t, maxNetwork.index))
                    self.logRates(charRate, maxSpecies, maxSpeciesRate, maxNetwork, maxNetworkRate)
                    self.logConversions(speciesIndex, y0)
                    break

            # Finish simulation if any of the termination criteria are satisfied
            for term in self.termination:
                if isinstance(term, TerminationTime):
                    if self.t > term.time.value:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination time.'.format(term.time.value))
                        self.logConversions(speciesIndex, y0)
                        break
                elif isinstance(term, TerminationConversion):
                    index = speciesIndex[term.species]
                    if (y0[index] - self.y[index]) / y0[index] > term.conversion:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination conversion: {1:f} of {2}'.format(self.t,term.conversion,term.species))
                        self.logConversions(speciesIndex, y0)
                        break

            # Increment destination step time if necessary
            if self.t >= 0.9999 * stepTime:
                stepTime *= 10.0

        self.maxCoreSpeciesRates = maxCoreSpeciesRates
        self.maxEdgeSpeciesRates = maxEdgeSpeciesRates
        self.maxNetworkLeakRates = maxNetworkLeakRates

        # Return the invalid object (if the simulation was invalid) or None
        # (if the simulation was valid)
        return terminated, invalidObject

    cpdef logRates(self, double charRate, object species, double speciesRate, object network, double networkRate):
        """
        Log information about the current maximum species and network rates.
        """
        logging.info('    Characteristic rate: {0:10.4e} mol/m^3*s'.format(charRate))
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
                X = (y0[index] - self.y[index]) / y0[index]
                logging.info('    {0} conversion: {1:<10.4g}'.format(term.species, X))

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
