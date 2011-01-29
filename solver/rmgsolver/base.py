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

import math
import numpy
from pydas import DASSL

################################################################################

class ReactionSystem(DASSL):
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

    def solve(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions,
        toleranceKeepInEdge, toleranceMoveToCore, toleranceInterruptSimulation,
        termination):
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

        speciesIndex = {}
        for index, spec in enumerate(coreSpecies):
            speciesIndex[spec] = index
        
        self.initialize(coreSpecies, coreReactions, edgeSpecies, edgeReactions)

        invalidObject = None
        terminated = False

        # Copy the initial conditions to use in evaluating conversions
        y0 = self.y.copy()

        while not terminated:
            # Integrate forward in time by one time step
            self.step(1.0)
            
            # Get the characteristic flux
            charRate = math.sqrt(numpy.sum(self.coreSpeciesRates * self.coreSpeciesRates))

            # Get the edge species with the highest flux
            maxIndex = numpy.argmax(self.edgeSpeciesRates)

            # Interrupt simulation if that flux exceeds the characteristic rate times a tolerance
            if self.edgeSpeciesRates[maxIndex] > toleranceMoveToCore * charRate and not invalidObject:
                invalidObject = edgeSpecies[maxIndex]
            if self.edgeSpeciesRates[maxIndex] > toleranceInterruptSimulation * charRate:
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

        # Return the invalid object (if the simulation was invalid) or None
        # (if the simulation was valid)
        return invalidObject

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
