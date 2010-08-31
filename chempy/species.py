#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

In ChemPy, a chemical species is called a Species object and is represented in
memory as an instance of the :class:`Species` class.
"""

################################################################################

class LennardJones:
    """
    A set of Lennard-Jones collision parameters. The Lennard-Jones parameters
    :math:`\\sigma` and :math:`\\epsilon` correspond to the potential

    .. math:: V(r) = 4 \\epsilon \\left[ \\left( \\frac{\\sigma}{r} \\right)^{12} - \\left( \\frac{\\sigma}{r} \\right)^{6} \\right]

    where the first term represents repulsion of overlapping orbitals and the
    second represents attraction due to van der Waals forces.

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `sigma`         ``double``      Distance at which the inter-particle potential is zero
    `epsilon`       ``double``      Depth of the potential well in J
    =============== =============== ============================================

    """

    def __init__(self, sigma=0.0, epsilon=0.0):
        self.sigma = sigma
        self.epsilon = epsilon

################################################################################

class Species:
    """
    A chemical species.

    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `index`             :class:`int`            A unique nonnegative integer index
    `label`             :class:`str`            A descriptive string label
    `thermo`            :class:`ThermoModel`    The thermodynamics model for the species
    `states`            :class:`StatesModel`    The molecular degrees of freedom model for the species
    `molecule`          ``list``                The :class:`Molecule` objects describing the molecular structure
    `geometry`          :class:`Geometry`       The 3D geometry of the molecule
    `E0`                ``double``              The ground-state energy in J/mol
    `lennardJones`      :class:`LennardJones`   A set of Lennard-Jones collision parameters
    `molecularWeight`   ``double``              The molecular weight of the species in kg/mol
    `reactive`          ``bool``                ``True`` if the species participates in reactions, ``False`` if not
    =================== ======================= ================================

    """

    def __init__(self, index=-1, label='', thermo=None, states=None, molecule=None, geometry=None, E0=0.0, lennardJones=None, molecularWeight=0.0, reactive=True):
        self.index = index
        self.label = label
        self.thermo = thermo
        self.states = states
        self.molecule = molecule or []
        self.geometry = geometry
        self.E0 = E0
        self.lennardJones = lennardJones
        self.reactive = reactive
        self.molecularWeight = molecularWeight

    def __repr__(self):
        """
        Return a string representation of the species, suitable for console output.
        """
        return "<Species %i '%s'>" % (self.index, self.label)

    def __str__(self):
        """
        Return a string representation of the species, in the form 'label(id)'.
        """
        if self.index == -1: return '%s' % (self.label)
        else: return '%s(%i)' % (self.label, self.index)

    def generateResonanceIsomers(self):
        """
        Generate all of the resonance isomers of this species. The isomers are
        stored as a list in the `molecule` attribute. If the length of
        `molecule` is already greater than one, it is assumed that all of the
        resonance isomers have already been generated.
        """

        if len(self.molecule) != 1:
            return

        # Radicals
        if sum([atom.radicalElectrons for atom in self.molecule[0].atoms]) > 0:
            # Iterate over resonance isomers
            index = 0
            while index < len(self.molecule):
                isomer = self.molecule[index]
                newIsomers = isomer.getAdjacentResonanceIsomers()
                for newIsomer in newIsomers:
                    # Append to isomer list if unique
                    found = False
                    for isom in self.molecule:
                        if isom.isIsomorphic(newIsomer): found = True
                    if not found:
                        self.molecule.append(newIsomer)
                        newIsomer.updateAtomTypes()
                # Move to next resonance isomer
                index += 1

################################################################################

class TransitionState:
    """
    A chemical transition state, representing a first-order saddle point on a
    potential energy surface.

    =============== =========================== ================================
    Attribute       Type                        Description
    =============== =========================== ================================
    `label`         :class:`str`                A descriptive string label
    `states`        :class:`StatesModel`        The molecular degrees of freedom model for the species
    `geometry`      :class:`Geometry`           The 3D geometry of the molecule
    `E0`            ``double``                  The ground-state energy in J/mol
    `frequency`     ``double``                  The negative frequency of the first-order saddle point in cm^-1
    =============== =========================== ================================

    """

    def __init__(self, label='', states=None, geometry=None, E0=0.0, frequency=0.0):
        self.label = label
        self.states = states
        self.geometry = geometry
        self.E0 = E0
        self.frequency = frequency

    def __repr__(self):
        """
        Return a string representation of the species, suitable for console output.
        """
        return "<TransitionState '%s'>" % (self.label)

