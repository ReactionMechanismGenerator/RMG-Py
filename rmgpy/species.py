#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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
This module contains classes and functions for working with chemical species.

From the `IUPAC Compendium of Chemical Terminology 
<http://dx.doi.org/10.1351/goldbook>`_, a chemical species is "an 
ensemble of chemically identical molecular entities that can explore the same 
set of molecular energy levels on the time scale of the experiment". This
definition is purposefully vague to allow the user flexibility in application.

In RMG Py, a chemical species -- a local minimum on a potential energy surface
-- is represented in memory as a :class:`Species` object. This module also
contains the :class:`TransitionState` class for representing chemical reaction
transition states (first-order saddle points on a potential energy surface).
"""

from quantity import Quantity, constants
from molecule import Molecule

################################################################################

class SpeciesError(Exception):
    """
    An exception class for exceptional behavior that occurs while working with
    chemical species. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

################################################################################

class LennardJones:
    """
    A set of Lennard-Jones collision parameters. The Lennard-Jones parameters
    :math:`\\sigma` and :math:`\\epsilon` correspond to the potential

    .. math:: V(r) = 4 \\epsilon \\left[ \\left( \\frac{\\sigma}{r} \\right)^{12} - \\left( \\frac{\\sigma}{r} \\right)^{6} \\right]

    where the first term represents repulsion of overlapping orbitals and the
    second represents attraction due to van der Waals forces. The attributes
    are:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `sigma`         :class:`Quantity`   Distance at which the inter-particle potential is zero
    `epsilon`       :class:`Quantity`   Depth of the potential well
    =============== =================== ========================================
    
    """

    def __init__(self, sigma=0.0, epsilon=0.0):
        self.sigma = Quantity(sigma)
        self.epsilon = Quantity(epsilon)
        if self.epsilon.units == 'K':
            # We also accept K as valid units for epsilon
            # Let's convert it to proper energy units
            self.epsilon.value *= constants.kB
            self.epsilon.uncertainty *= constants.kB
            self.epsilon.units = 'J'

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        return 'LennardJones(sigma={0!r}, epsilon={1!r})'.format(self.sigma, self.epsilon)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (LennardJones, (self.sigma, self.epsilon))

################################################################################

class Species:
    """
    A chemical species, representing a local minimum on a potential energy
    surface. The attributes are:

    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `index`             :class:`int`            A unique nonnegative integer index
    `label`             :class:`str`            A descriptive string label
    `thermo`            :class:`ThermoModel`    The thermodynamics model for the species
    `states`            :class:`StatesModel`    The molecular degrees of freedom model for the species
    `molecule`          ``list``                The :class:`Molecule` objects describing the molecular structure
    `E0`                :class:`Quantity`       The ground-state energy
    `lennardJones`      :class:`LennardJones`   A set of Lennard-Jones collision parameters
    `molecularWeight`   :class:`Quantity`       The molecular weight of the species
    `reactive`          ``bool``                ``True`` if the species participates in reactions, ``False`` if not
    =================== ======================= ================================

    """

    def __init__(self, index=-1, label='', thermo=None, states=None, molecule=None, E0=None, lennardJones=None, molecularWeight=None, reactive=True):
        self.index = index
        self.label = label
        self.thermo = thermo
        self.states = states
        self.molecule = molecule or []
        self.E0 = Quantity(E0)
        self.lennardJones = lennardJones
        self.reactive = reactive
        self.molecularWeight = Quantity(molecularWeight)

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'Species('
        if self.index != -1: string += 'index={0:d}, '.format(self.index)
        if self.label != -1: string += 'label="{0}", '.format(self.label)
        if self.thermo is not None: string += 'thermo={0!r}, '.format(self.thermo)
        if self.states is not None: string += 'states={0!r}, '.format(self.states)
        if len(self.molecule) > 0: string += 'molecule=[{0!r}], '.format(self.molecule[0])
        if self.E0 is not None: string += 'E0={0!r}, '.format(self.E0)
        if self.lennardJones is not None: string += 'lennardJones={0!r}, '.format(self.lennardJones)
        if not self.reactive: string += 'reactive={0}, '.format(self.reactive)
        if self.molecularWeight is not None: string += 'molecularWeight={0!r}, '.format(self.molecularWeight)
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the species, in the form 'label(id)'.
        """
        if self.index == -1: return self.label
        else: return '{0}({1:d})'.format(self.label, self.index)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Species, (self.index, self.label, self.thermo, self.states, self.molecule, self.E0, self.lennardJones, self.molecularWeight, self.reactive))

    def generateResonanceIsomers(self):
        """
        Generate all of the resonance isomers of this species. The isomers are
        stored as a list in the `molecule` attribute. If the length of
        `molecule` is already greater than one, it is assumed that all of the
        resonance isomers have already been generated.
        """
        if len(self.molecule) == 1:
            self.molecule = self.molecule[0].generateResonanceIsomers()
    
    def isIsomorphic(self, other):
        """
        Return ``True`` if the species is isomorphic to `other`, which can be
        either a :class:`Molecule` object or a :class:`Species` object.
        """
        if isinstance(other, Molecule):
            for molecule in self.molecule:
                if molecule.isIsomorphic(other):
                    return True
        elif isinstance(other, Species):
            for molecule1 in self.molecule:
                for molecule2 in other.molecule:
                    if molecule1.isIsomorphic(molecule2):
                        return True
        else:
            raise ValueError('Unexpected value "{0!r}" for other parameter; should be a Molecule or Species object.'.format(other))
        return False
    
    def toAdjacencyList(self):
        """
        Return a string containing each of the molecules' adjacency lists.
        """
        output = '\n\n'.join([m.toAdjacencyList(label=self.label, removeH=True) for m in self.molecule])
        return output
            
    
################################################################################

class TransitionState:
    """
    A chemical transition state, representing a first-order saddle point on a
    potential energy surface. The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `label`         :class:`str`            A descriptive string label
    `states`        :class:`StatesModel`    The molecular degrees of freedom model for the species
    `E0`            :class:`Quantity`       The ground-state energy in J/mol
    `frequency`     :class:`Quantity`       The negative frequency of the first-order saddle point in cm^-1
    `degeneracy`    ``int``                 The reaction path degeneracy
    =============== ======================= ====================================

    """

    def __init__(self, label='', states=None, E0=None, frequency=None, degeneracy=1):
        self.label = label
        self.states = states
        self.E0 = Quantity(E0)
        self.frequency = Quantity(frequency)
        self.degeneracy = degeneracy

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'TransitionState('
        if self.label != '': string += 'label="{0}", '.format(self.label)
        if self.states is not None: string += 'states={0!r}, '.format(self.states)
        if self.E0 is not None: string += 'E0={0!r}, '.format(self.E0)
        if self.frequency is not None: string += 'frequency={0!r}, '.format(self.frequency)
        if self.degeneracy != 1: string += 'degeneracy={0}, '.format(self.degeneracy)
        string = string[:-2] + ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TransitionState, (self.label, self.states, self.E0, self.frequency, self.degeneracy))
