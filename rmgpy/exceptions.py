#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
This module contains classes which extend Exception for usage in the RMG module
"""

class ActionError(Exception):
    """
    An exception class for errors that occur while applying reaction recipe
    actions. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

class AtomTypeError(Exception):
    """
    An exception to be raised when an error occurs while working with atom
    types. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

class ChemicallySignificantEigenvaluesError(Exception):
    """
    An exception raised when the chemically significant eigenvalue method is 
    unsuccessful for any reason. 
    Pass a string describing the cause of the exceptional behavior.
    """
    pass

class ChemkinError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin files. Pass a
    string describing the circumstances that caused the exceptional behavior.
    """
    pass

class CollisionError(Exception):
    """
    An exception class for when RMG is unable to calculate collision efficiencies
    for the single exponential down pressure dependent solver. 
    Pass a string describing the circumstances that caused the exceptional behavior.
    """
    pass

class DatabaseError(Exception):
    """
    A exception that occurs when working with an RMG database. Pass a string
    giving specifics about the exceptional behavior.
    """
    pass

class DependencyError(Exception):
    """
    An exception that occurs when an error is encountered with a dependency.
    Pass a string describing the circumstances that caused the exception.
    """
    pass

class ElementError(Exception):
    """
    An exception class for errors that occur while working with elements.
    Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

class ForbiddenStructureException(Exception):
    """
    An exception passed when RMG encounters a forbidden structure. These are usually
    caught and the reaction that created it is ignored.
    """
    pass

class ILPSolutionError(Exception):
    """
    An exception to be raised when solving an integer linear programming problem if a solution
    could not be found or the solution is not valid. Can pass a string to indicate the reason
    that the solution is invalid.
    """
    pass

class ImplicitBenzeneError(Exception):
    """
    An exception class when encountering a group with too many implicit benzene
    atoms. These groups are hard to create sample molecules and hard for users
    to interpret. Pass a string describing the limitation.
    """
    pass

class InchiException(Exception):
    """
    An exception used when encountering a non-valid Inchi expression are encountered.
    Pass a string describing the error.
    """
    pass

class InputError(Exception):
    """
    An exception raised when parsing an input file for any module in RMG:
    mechanism generation, cantherm, conformer creation, etc.
    Pass a string describing the error.
    """
    pass

class InvalidActionError(Exception):
    """
    An exception to be raised when an invalid action is encountered in a
    reaction recipe.
    """
    pass

class InvalidAdjacencyListError(Exception):
    """
    An exception used to indicate that an RMG-style adjacency list is invalid.
    Pass a string describing the reason the adjacency list is invalid
    """
    pass

class KekulizationError(Exception):
    """
    An exception to be raised when encountering an error while kekulizing an aromatic molecule.
    Can pass a string to indicate the reason for failure.
    """
    pass

class KineticsError(Exception):
    """
    An exception class for problems with kinetics. This can be used when finding 
    degeneracy in reaction generation, modifying KineticsData objects, or finding the
    kinetics of reactions. Unable Pass a string describing the problem.
    """
    pass

class ModifiedStrongCollisionError(Exception): 
    """
    An exception raised when the modified strong collision method is unsuccessful
    for any reason. Pass a string describing the cause of the exceptional 
    behavior.
    """
    pass

class NegativeBarrierException(Exception):
    """ This Exception occurs when the energy barrier for a hindered Rotor is negative.
    This can occur if the scan or fourier fit is poor. """
    
    pass

class NetworkError(Exception): 
    """
    Raised when an error occurs while working with a pressure-dependent reaction network
    """
    pass

class OutputError(Exception):
    """
    This exception is raised whenever an error occurs while saving output
    information. Pass a string describing the circumstances of the exceptional
    behavior.
    """
    pass

class PressureDependenceError(Exception):
    """
    An exception class to use when an error involving pressure dependence is
    encountered. Pass a string describing the circumstances of the exceptional
    behavior.
    """
    pass

class QuantityError(Exception):
    """
    An exception to be raised when an error occurs while working with physical
    quantities in RMG. Pass a string describing the circumstances of the
    exceptional behavior.
    """
    pass

class ReactionError(Exception):
    """
    An exception class for exceptional behavior involving :class:`Reaction`
    objects. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

class ReactionPairsError(Exception):
    """
    An exception to be raised when an error occurs while working with reaction
    pairs.
    """
    pass

class ReservoirStateError(Exception): 
    """
    An exception raised when the reservoir state method is unsuccessful for
    any reason. Pass a string describing the cause of the exceptional behavior.
    """
    pass

class SettingsError(Exception):
    """
    An exception raised when dealing with settings.
    """
    pass

class SpeciesError(Exception):
    """
    An exception class for exceptional behavior that occurs while working with
    chemical species. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

class StatmechError(Exception):
    """
    An exception used when an error occurs in estimating Statmech.
    """

class StatmechFitError(StatmechError):
    """
    An exception used when attempting to fit molecular degrees of freedom to
    heat capacity data. Pass a string describing the circumstances of the
    exceptional behavior.
    """
    pass

class UnexpectedChargeError(Exception):
    """
    An exception class when encountering a group/molecule with unexpected charge
    Curently in RMG, we never expect to see -2/+2 or greater magnitude charge,
    we only except +1/-1 charges on nitrogen, oxygen, sulfur or specifically
    carbon monoxide/monosulfide.

    Attributes:
    `graph` is the molecule or group object with the unexpected charge
    """
    def __init__(self, graph):
        self.graph = graph

class VF2Error(Exception):
    """
    An exception raised if an error occurs within the VF2 graph isomorphism
    algorithm. Pass a string describing the error.
    """
    pass

class CoreError(Exception):
    """
    An exception raised if there is a problem within the model core
    """
    pass

class ResonanceError(Exception):
    """
    An exception class for when RMG is unable to generate resonance structures.
    """
    pass

################## move classes that extend off previous exceptions here

class InvalidMicrocanonicalRateError(NetworkError):
    """
    Used in pressure dependence when the k(E) calculation does not give 
    the correct kf(T) or Kc(T)
    """
    def __init__(self,message, k_ratio=1.0, Keq_ratio=1.0):
        self.message = message
        self.k_ratio = k_ratio
        self.Keq_ratio = Keq_ratio
    def badness(self):
        """
        How bad is the error?

        Returns the max of the absolute logarithmic errors of kf and Kc
        """
        import math
        return max(abs(math.log10(self.k_ratio)), abs(math.log10(self.Keq_ratio)))

class UndeterminableKineticsError(ReactionError):
    """
    An exception raised when attempts to estimate appropriate kinetic parameters
    for a chemical reaction are unsuccessful.
    """
    def __init__(self, reaction, message=''):
        new_message = 'Kinetics could not be determined. '+message
        ReactionError.__init__(self,reaction,new_message)
