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
Each atom in a molecular configuration has three spatial dimensions in which it
can move. Thus, a molecular configuration consisting of :math:`N` atoms has
:math:`3N` degrees of freedom. We can distinguish between those modes that
involve movement of atoms relative to the molecular center of mass (called
*internal* modes) and those that do not (called *external* modes). Of the
external degrees of freedom, three involve translation of the entire molecular
configuration, while either three (for a nonlinear molecule) or two (for a
linear molecule) involve rotation of the entire molecular configuration
around the center of mass. The remaining :math:`3N-6` (nonlinear) or
:math:`3N-5` (linear) degrees of freedom are the internal modes, and can be
divided into those that involve vibrational motions (symmetric and asymmetric
stretches, bends, etc.) and those that involve torsional rotation around single
bonds between nonterminal heavy atoms.

The mathematical description of these degrees of freedom falls under the purview
of quantum chemistry, and involves the solution of the time-independent
Schrodinger equation:

    .. math:: \\hat{H} \\psi = E \\psi

where :math:`\\hat{H}` is the Hamiltonian, :math:`\\hat{H}` is the wavefunction,
and :math:`E` is the energy. The exact form of the Hamiltonian varies depending
on the degree of freedom you are modeling. Since this is a quantum system, the
energy can only take on discrete values. Once the allowed energy levels are
known, the partition function :math:`Q(\\beta)` can be computed using the
summation

    .. math:: Q(\\beta) = \\sum_i g_i e^{-\\beta E_i}

where :math:`g_i` is the degeneracy of energy level :math:`i` (i.e. the number
of energy states at that energy level) and
:math:`\\beta \\equiv (k_\\mathrm{B} T)^{-1}`.

The partition function is an immensely useful quantity, as all sorts of
thermodynamic parameters can be evaluated using the partition function:

    .. math:: A = - k_\\mathrm{B} T \\ln Q

    .. math:: U = - \\frac{\\partial \\ln Q}{\\partial \\beta}

    .. math:: S = \\frac{\\partial}{\\partial T} \\left( k_\\mathrm{B} T \\ln Q \\right)

    .. math:: C_\\mathrm{v} = \\frac{1}{k_\\mathrm{B} T} \\frac{\\partial^2 \\ln Q}{\\partial \\beta^2}

Above, :math:`A`, :math:`U`, :math:`S`, and :math:`C_\\mathrm{v}` are the
Helmholtz free energy, internal energy, entropy, and constant-volume heat
capacity, respectively.

The partition function for a molecular configuration is the product of the
partition functions for each invidual degree of freedom:

    .. math:: Q = Q_\\mathrm{trans} Q_\\mathrm{rot} Q_\\mathrm{vib} Q_\\mathrm{tors} Q_\\mathrm{elec}

This means that the contributions to each thermodynamic quantity from each
molecular degree of freedom are additive.

This module contains models for various molecular degrees of freedom. All such
models derive from the :class:`Mode` base class. A list of molecular degrees of
freedom can be stored in a :class:`StatesModel` object.
"""

################################################################################

import math
import cython
import numpy

import constants
from exception import InvalidStatesModelError

################################################################################

class Mode:
    pass

################################################################################

class Translation(Mode):
    """
    A representation of translational modes. The `dimension` attribute is an
    integer representing the dimension of the translation (2 or 3) and the
    `mass` is the molar mass of the molecule in kg/mol.
    """

    def __init__(self, mass=0.0, volume=1.0, dimension=3):
        self.mass = mass
        self.volume = volume
        self.dimension = dimension

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the 
        object.
        """
        return 'Translation(mass=%g, volume=%g, dimension=%g)' % (self.mass, self.volume, self.dimension)

    def getPartitionFunction(self, Tlist):
        """
        Return the value of the partition function at the specified temperatures
        `Tlist` in K. The formula is

        .. math:: q_\\mathrm{trans}(T, V) = \\left( \\frac{2 \\pi m k_\\mathrm{B} T}{h^2} \\right)^{d/2} V

        where :math:`T` is temperature, :math:`V` is volume, :math:`m` is mass,
        :math:`d` is dimensionality, :math:`k_\\mathrm{B}` is the Boltzmann
        constant, and :math:`h` is the Planck constant.
        """
        cython.declare(qt=cython.double)
        qt = ((2 * constants.pi * self.mass / constants.Na) / (constants.h * constants.h))**(self.dimension/2.0) * self.volume
        return qt * (constants.kB * Tlist)**(self.dimension/2.0)
    
    def getHeatCapacity(self, Tlist):
        """
        Return the contribution to the heat capacity due to translation in
        J/mol*K at the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{C_\\mathrm{v}^\\mathrm{trans}(T)}{R} = \\frac{d}{2}

        where :math:`T` is temperature,	:math:`d` is dimensionality, and
        :math:`R` is the gas law constant.
        """
        return 0.5 * self.dimension * numpy.ones_like(Tlist) * constants.R
    
    def getEnthalpy(self, Tlist):
        """
        Return the contribution to the enthalpy due to translation in J/mol
        at the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{H^\\mathrm{trans}(T)}{RT} = \\frac{d}{2}

        where :math:`T` is temperature, :math:`d` is dimensionality, and
        :math:`R` is the gas law constant.
        """
        return 0.5 * self.dimension * constants.R * Tlist
    
    def getEntropy(self, Tlist):
        """
        Return the contribution to the entropy due to translation in J/mol*K
        at the specified temperatures `Tlist` in K. The formula	is

        .. math:: \\frac{S^\\mathrm{trans}(T)}{R} = \\frac{d}{2} \\ln \\left( \\frac{2 \\pi m k_\\mathrm{B} T}{h^2} \\right) + \\ln eV
        
        where :math:`T` is temperature, :math:`m` is mass, :math:`d` is
        dimensionality, :math:`V` is volume, :math:`k_\\mathrm{B}` is the
        Boltzmann constant, and :math:`R` is the gas law constant.
        """
        return (numpy.log(self.getPartitionFunction(Tlist)) + 0.5 * self.dimension + 1.0) * constants.R
    
    def getDensityOfStates(self, Elist):
        """
        Return the density of states at the specified energlies `Elist` in J/mol
        above the ground state. The formula is

        .. math:: \\rho(E) = \\left( \\frac{2 \\pi m}{h^2} \\right)^{d/2} \\frac{E^{d/2-1}}{(d/2-1)!} V

        where :math:`E` is energy, :math:`m` is mass, :math:`d` is
        dimensionality, :math:`k_\\mathrm{B}` is the Boltzmann constant, and
        :math:`R` is the gas law constant.
        """
        cython.declare(rho=numpy.ndarray, qt=cython.double)
        rho = numpy.zeros_like(Elist)
        qt = ((2 * constants.pi * self.mass / constants.Na / constants.Na) / (constants.h * constants.h))**(self.dimension/2.0) * self.volume
        if self.dimension == 2: # Dimension is 2
            rho = qt * numpy.ones_like(Elist)
        elif self.dimension == 3: # Dimension is 3
            rho = qt * numpy.sqrt(Elist / math.pi) * 0.5
        else:
            raise InvalidStatesModelError('Unexpected dimensionality "%i" encountered while evaluating translation density of states.')
        return rho

################################################################################

class RigidRotor(Mode):
    """
    A rigid rotor approximation of (external) rotational modes. The `linear`
    attribute is :data:`True` if the associated molecule is linear, and
    :data:`False` if nonlinear. For a linear molecule, `inertia` stores a
    list with one moment of inertia in kg*m^2. For a nonlinear molecule,
    `frequencies` stores a list of the three moments of inertia, even if two or
    three are equal, in kg*m^2. The symmetry number of the rotation is stored
    in the `symmetry` attribute.
    """

    def __init__(self, linear=False, inertia=None, symmetry=1):
        self.linear = linear
        self.inertia = inertia or []
        self.symmetry = symmetry

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the 
        object.
        """
        inertia = ', '.join(['%g' % i for i in self.inertia])
        return 'RigidRotor(linear=%s, inertia=[%s], symmetry=%s)' % (self.linear, inertia, self.symmetry)

    def getPartitionFunction(self, Tlist):
        """
        Return the value of the partition function at the specified temperatures
        `Tlist` in K. The formula is

        .. math:: q_\\mathrm{rot}(T) = \\frac{8 \\pi^2 I k_\\mathrm{B} T}{\\sigma h^2}

        for linear rotors and

        .. math:: q_\\mathrm{rot}(T) = \\frac{\\sqrt{\\pi}}{\\sigma} \\left( \\frac{8 \\pi^2 k_\\mathrm{B} T}{h^2} \\right)^{3/2} \\sqrt{I_\\mathrm{A} I_\\mathrm{B} I_\\mathrm{C}}

        for nonlinear rotors. Above, :math:`T` is temperature, :math:`\\sigma`
        is the symmetry number, :math:`I` is the moment of inertia, 
        :math:`k_\\mathrm{B}` is the Boltzmann constant, and :math:`h` is the
        Planck constant.
        """
        cython.declare(theta=cython.double, inertia=cython.double)
        if self.linear:
            theta = constants.h * constants.h / (8 * constants.pi * constants.pi * self.inertia[0] * constants.kB)
            return Tlist / theta / self.symmetry
        else:
            theta = 1.0
            for inertia in self.inertia:
                theta *= constants.h * constants.h / (8 * constants.pi * constants.pi * inertia * constants.kB)
            return numpy.sqrt(constants.pi * Tlist**len(self.inertia) / theta) / self.symmetry
        
    def getHeatCapacity(self, Tlist):
        """
        Return the contribution to the heat capacity due to rigid rotation
        in J/mol*K at the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{C_\\mathrm{v}^\\mathrm{rot}(T)}{R} = 1

        if linear and

        .. math:: \\frac{C_\\mathrm{v}^\\mathrm{rot}(T)}{R} = \\frac{3}{2}

        if nonlinear, where :math:`T` is temperature and :math:`R` is the gas
        law constant.
        """
        if self.linear:
            return numpy.ones_like(Tlist) * constants.R
        else:
            return 1.5 * numpy.ones_like(Tlist) * constants.R

    def getEnthalpy(self, Tlist):
        """
        Return the contribution to the enthalpy due to rigid rotation in J/mol
        at the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{H^\\mathrm{rot}(T)}{RT} = 1

        for linear rotors and

        .. math:: \\frac{H^\\mathrm{rot}(T)}{RT} = \\frac{3}{2}

        for nonlinear rotors, where :math:`T` is temperature and :math:`R` is
        the gas law constant.
        """
        if self.linear:
            return constants.R * Tlist
        else:
            return 1.5 * constants.R * Tlist
    
    def getEntropy(self, Tlist):
        """
        Return the contribution to the entropy due to rigid rotation in J/mol*K
        at the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{S^\\mathrm{rot}(T)}{R} = \\ln Q^\\mathrm{rot} + 1
        
        for linear rotors and

        .. math:: \\frac{S^\\mathrm{rot}(T)}{R} = \\ln Q^\\mathrm{rot} + \\frac{3}{2}

        for nonlinear rotors, where :math:`Q^\\mathrm{rot}` is the partition 
        function for a rigid rotor and :math:`R` is the gas law constant.
        """
        if self.linear:
            return (numpy.log(self.getPartitionFunction(Tlist)) + 1.0) * constants.R
        else:
            return (numpy.log(self.getPartitionFunction(Tlist)) + 1.5) * constants.R
    
    def getDensityOfStates(self, Elist):
        """
        Return the density of states at the specified energlies `Elist` in J/mol
        above the ground state in mol/J. The formula is

        .. math:: \\rho(E) = \\frac{8 \\pi^2 I}{\\sigma h^2}

        for linear rotors and

        .. math:: \\rho(E) = \\frac{\\sqrt{\\pi}}{\\sigma} \\left( \\frac{8 \\pi^2}{h^2} \\right)^{3/2} \\sqrt{I_\\mathrm{A} I_\\mathrm{B} I_\\mathrm{C}} \\frac{E^{1/2}}{\\frac{1}{2}!}

        for nonlinear rotors. Above, :math:`E` is energy, :math:`\\sigma`
        is the symmetry number, :math:`I` is the moment of inertia,
        :math:`k_\\mathrm{B}` is the Boltzmann constant, and :math:`h` is the
        Planck constant.
        """
        cython.declare(theta=cython.double, inertia=cython.double)
        if self.linear:
            theta = constants.h * constants.h / (8 * constants.pi * constants.pi * self.inertia[0]) * constants.Na
            return numpy.ones_like(Elist) / theta / self.symmetry
        else:
            theta = 1.0
            for inertia in self.inertia:
                theta *= constants.h * constants.h / (8 * constants.pi * constants.pi * inertia) * constants.Na
            return 2.0 * numpy.sqrt(Elist / theta) / self.symmetry
        
################################################################################

class HinderedRotor(Mode):
    """
    A one-dimensional hindered rotor using the simple potential function

    .. math:: V(\\phi) = \\frac{1}{2} V_0 \\left[1 - \\cos \\left( \\sigma \\phi \\right) \\right]

    where :math:`V_0` is the height of the potential barrier and
    :math:`\\sigma` is the number of minima or maxima in one revolution of
    angle :math:`\\phi`, equivalent to the symmetry number of that rotor. The
    hindered rotor is therefore described by three quantities: the moment of
    `inertia` in kg*m^2, the `barrier` height in J/mol, and the `symmetry`
    number.
    """

    def __init__(self, inertia=None, barrier=None, symmetry=None):
        self.inertia = inertia
        self.barrier = barrier
        self.symmetry = symmetry

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the 
        object.
        """
        return 'HinderedRotor(inertia=%g, barrier=%g, symmetry=%g)' % (self.inertia, self.barrier, self.symmetry)

    def getPartitionFunction(self, Tlist):
        """
        Return the value of the partition function at the specified temperatures
        `Tlist` in K. The formula makes use of the Pitzer-Gwynn approximation:

        .. math:: q_\\mathrm{hind}(T) = \\frac{q_\\mathrm{vib}^\\mathrm{quant}(T)}{q_\\mathrm{vib}^\\mathrm{class}(T)} q_\\mathrm{hind}^\\mathrm{class}(T)

        Substituting in for the right-hand side partition functions gives

        .. math:: q_\\mathrm{hind}(T) = \\frac{h \\nu}{k_\\mathrm{B} T} \\frac{1}{1 - \\exp \\left(- h \\nu / k_\\mathrm{B} T \\right)} \\left( \\frac{2 \\pi I k_\\mathrm{B} T}{h^2} \\right)^{1/2} \\frac{2 \\pi}{\\sigma} \\exp \\left( -\\frac{V_0}{2 k_\\mathrm{B} T} \\right) I_0 \\left( \\frac{V_0}{2 k_\\mathrm{B} T} \\right)

        where
        
        .. math:: \\nu = \\frac{\\sigma}{2 \\pi} \\sqrt{\\frac{V_0}{2 I}}

        :math:`T` is temperature, :math:`V_0` is the barrier height,
        :math:`I` is the moment of inertia, :math:`\\sigma` is the symmetry
        number, :math:`k_\\mathrm{B}` is the Boltzmann constant, and :math:`h`
        is the Planck constant. :math:`I_0(x)` is the modified Bessel function
        of order zero for argument :math:`x`.
        """
        cython.declare(frequency=cython.double, x=numpy.ndarray, z=numpy.ndarray)
        frequency = self.getFrequency() * constants.c * 100
        x = constants.h * frequency / (constants.kB * Tlist)
        z = 0.5 * self.barrier / (constants.R * Tlist)
        return x / (1 - numpy.exp(-x)) * numpy.sqrt(2 * math.pi * self.inertia * constants.kB * Tlist / constants.h / constants.h) * (2 * math.pi / self.symmetry) * numpy.exp(-z) * besseli0(z)
    
    def getHeatCapacity(self, Tlist):
        """
        Return the contribution to the heat capacity due to hindered rotation
        in J/mol*K at the specified temperatures `Tlist`
        in K. The formula is

        .. math:: \\frac{C_\\mathrm{v}^\\mathrm{hind}(T)}{R} = \\frac{C_\\mathrm{v}^\\mathrm{vib}(T)}{R} -\\frac{1}{2} + \\zeta^2 - \\left[ \\zeta \\frac{I_1(\\zeta)}{I_0(\\zeta)} \\right]^2 - \\zeta \\frac{I_1(\\zeta)}{I_0(\\zeta)}

        where :math:`\\zeta \\equiv V_0 / 2 k_\\mathrm{B} T`,
        :math:`T` is temperature, :math:`V_0` is the barrier height,
        :math:`k_\\mathrm{B}` is the Boltzmann constant, and :math:`R` is the
        gas law constant.
        """
        cython.declare(frequency=cython.double, x=numpy.ndarray, z=numpy.ndarray)
        cython.declare(exp_x=numpy.ndarray, one_minus_exp_x=numpy.ndarray, BB=numpy.ndarray)
        frequency = self.getFrequency() * constants.c * 100
        x = constants.h * frequency / (constants.kB * Tlist)
        z = 0.5 * self.barrier / (constants.R * Tlist)
        exp_x = numpy.exp(x)
        one_minus_exp_x = 1.0 - exp_x
        BB = besseli1(z) / besseli0(z)
        return (x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - BB - z * BB * BB)) * constants.R
        
    def getEnthalpy(self, Tlist):
        """
        Return the contribution to the heat capacity due to hindered rotation
        in J/mol at the specified temperatures `Tlist` in K. This is calculated
        numerically from the partition function.
        """
        Tlist_low = Tlist * 0.999
        Tlist_high = Tlist * 1.001
        return (Tlist *
            (numpy.log(self.getPartitionFunction(Tlist_high)) -
            numpy.log(self.getPartitionFunction(Tlist_low))) /
            (Tlist_high - Tlist_low) + 1.0) * constants.R * Tlist

    def getEntropy(self, Tlist):
        """
        Return the contribution to the heat capacity due to hindered rotation
        in J/mol*K at the specified temperatures `Tlist` in K. This is
        calculated numerically from the partition function.
        """
        Tlist_low = Tlist * 0.999
        Tlist_high = Tlist * 1.001
        return (numpy.log(self.getPartitionFunction(Tlist_high)) +
            Tlist * (numpy.log(self.getPartitionFunction(Tlist_high)) -
            numpy.log(self.getPartitionFunction(Tlist_low))) /
            (Tlist_high - Tlist_low)) * constants.R

    def getDensityOfStates(self, Elist):
        """
        Return the density of states at the specified energlies `Elist` in J/mol
        above the ground state. The formula is

        .. math:: \\rho(E) = \\frac{2 q_\\mathrm{1f}}{\\pi^{3/2} V_0^{1/2}} \\mathcal{K}(E / V_0) \\hspace{20pt} E < V_0

        and

        .. math:: \\rho(E) = \\frac{2 q_\\mathrm{1f}}{\\pi^{3/2} E^{1/2}} \\mathcal{K}(V_0 / E) \\hspace{20pt} E > V_0

        where

        .. math:: q_\\mathrm{1f} = \\frac{\\pi^{1/2}}{\\sigma} \\left( \\frac{8 \\pi^2 I}{h^2} \\right)^{1/2}

        :math:`E` is energy, :math:`V_0` is barrier height, and
        :math:`\\mathcal{K}(x)` is the complete elliptic integral of the first
        kind.
        """
        cython.declare(rho=numpy.ndarray, q1f=cython.double, pre=cython.double, V0=cython.double, i=cython.int)
        rho = numpy.zeros_like(Elist)
        q1f = math.sqrt(8 * math.pi * math.pi * math.pi * self.inertia / constants.h / constants.h / constants.Na) / self.symmetry
        V0 = self.barrier
        pre = 2.0 * q1f / math.sqrt(math.pi * math.pi * math.pi * V0)
        # The following is only valid in the classical limit
        # Note that cellipk(1) = infinity, so we must skip that value
        for i in range(len(Elist)):
            if Elist[i] / V0 < 1:
                rho[i] = pre * cellipk(Elist[i] / V0)
            elif Elist[i] / V0 > 1:
                rho[i] = pre * math.sqrt(V0 / Elist[i]) * cellipk(V0 / Elist[i])
        return rho

    def getFrequency(self):
        """
        Return the frequency of vibration corresponding to the limit of
        harmonic oscillation. For the cosine potential, the formula is

        .. math:: \\nu = \\frac{\\sigma}{2 \\pi} \\sqrt{\\frac{V_0}{2 I}}

        where :math:`\\sigma` is the symmetry number, :math:`V_0` the barrier
        height, and :math:`I` the reduced moment of inertia of the rotor. The
        units of the returned frequency are cm^-1.
        """
        return self.symmetry / 2.0 / math.pi * math.sqrt(self.barrier / constants.Na / 2 / self.inertia) / (constants.c * 100)

def besseli0(xlist):
    """
    Return the value of the zeroth-order modified Bessel function at `x`.
    """
    import scipy.special
    return scipy.special.i0(xlist)

def besseli1(xlist):
    """
    Return the value of the first-order modified Bessel function at `x`.
    """
    import scipy.special
    return scipy.special.i1(xlist)

def cellipk(x):
    """
    Return the value of the complete elliptic integral of the first kind at `x`.
    """
    import scipy.special
    return scipy.special.ellipk(x)

################################################################################

class HarmonicOscillator(Mode):
    """
    A representation of a set of vibrational modes as one-dimensional quantum 
    harmonic oscillator. The oscillators are defined by their `frequencies` in 
    cm^-1.
    """

    def __init__(self, frequencies=None):
        self.frequencies = frequencies or []
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the 
        object.
        """
        frequencies = ', '.join(['%g' % freq for freq in self.frequencies])
        return 'HarmonicOscillator(frequencies=[%s])' % (frequencies)

    def getPartitionFunction(self, Tlist):
        """
        Return the value of the partition function at the specified temperatures
        `Tlist` in K. The formula is

        .. math:: q_\\mathrm{vib}(T) = \\prod_i \\frac{1}{1 - e^{-\\xi_i}}
        
        where :math:`\\xi_i \\equiv h \\nu_i / k_\\mathrm{B} T`,
        :math:`T` is temperature, :math:`\\nu_i` is the frequency of vibration
        :math:`i`, :math:`k_\\mathrm{B}` is the Boltzmann constant, :math:`h` 
        is the Planck constant, and :math:`R` is the gas law constant. Note 
        that we have chosen our zero of energy to be at the zero-point energy
        of the molecule, *not* the bottom of the potential well.
        """
        cython.declare(Q=numpy.ndarray, freq=cython.double)
        Q = numpy.ones_like(Tlist)
        for freq in self.frequencies:
            Q = Q / (1 - numpy.exp(-freq / (0.695039 * Tlist)))  # kB = 0.695039 cm^-1/K
        return Q
    
    def getHeatCapacity(self, Tlist):
        """
        Return the contribution to the heat capacity due to vibration
        in J/mol*K at the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{C_\\mathrm{v}^\\mathrm{vib}(T)}{R} = \\sum_i \\xi_i^2 \\frac{e^{\\xi_i}}{\\left( 1 - e^{\\xi_i} \\right)^2}

        where :math:`\\xi_i \\equiv h \\nu_i / k_\\mathrm{B} T`,
        :math:`T` is temperature, :math:`\\nu_i` is the frequency of vibration
        :math:`i`, :math:`k_\\mathrm{B}` is the Boltzmann constant, :math:`h` 
        is the Planck constant, and :math:`R` is the gas law constant. 
        """
        cython.declare(Cv=numpy.ndarray, freq=cython.double)
        cython.declare(x=numpy.ndarray, exp_x=numpy.ndarray, one_minus_exp_x=numpy.ndarray)
        Cv = numpy.zeros_like(Tlist)
        for freq in self.frequencies:
            x = freq / (0.695039 * Tlist)	# kB = 0.695039 cm^-1/K
            exp_x = numpy.exp(x)
            one_minus_exp_x = 1.0 - exp_x
            Cv = Cv + x * x * exp_x / one_minus_exp_x / one_minus_exp_x
        return Cv * constants.R

    def getEnthalpy(self, Tlist):
        """
        Return the contribution to the enthalpy due to vibration in J/mol at
        the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{H^\\mathrm{vib}(T)}{RT} = \\sum_i \\frac{\\xi_i}{e^{\\xi_i} - 1}
        
        where :math:`\\xi_i \\equiv h \\nu_i / k_\\mathrm{B} T`,
        :math:`T` is temperature, :math:`\\nu_i` is the frequency of vibration
        :math:`i`, :math:`k_\\mathrm{B}` is the Boltzmann constant, :math:`h` 
        is the Planck constant, and :math:`R` is the gas law constant. 
        """
        cython.declare(H=numpy.ndarray, freq=cython.double)
        cython.declare(x=numpy.ndarray, exp_x=numpy.ndarray)
        H = numpy.zeros_like(Tlist)
        for freq in self.frequencies:
            x = freq / (0.695039 * Tlist)	# kB = 0.695039 cm^-1/K
            exp_x = numpy.exp(x)
            H = H + x / (exp_x - 1)
        return H * constants.R * Tlist

    def getEntropy(self, Tlist):
        """
        Return the contribution to the entropy due to vibration in J/mol*K at
        the specified temperatures `Tlist` in K. The formula is

        .. math:: \\frac{S^\\mathrm{vib}(T)}{R} = \\sum_i \\left[ - \\ln \\left(1 - e^{-\\xi_i} \\right) + \\frac{\\xi_i}{e^{\\xi_i} - 1} \\right]
        
        where :math:`\\xi_i \\equiv h \\nu_i / k_\\mathrm{B} T`,
        :math:`T` is temperature, :math:`\\nu_i` is the frequency of vibration
        :math:`i`, :math:`k_\\mathrm{B}` is the Boltzmann constant, :math:`h` 
        is the Planck constant, and :math:`R` is the gas law constant. 
        """
        cython.declare(S=numpy.ndarray, freq=cython.double)
        cython.declare(x=numpy.ndarray, exp_x=numpy.ndarray)
        S = numpy.log(self.getPartitionFunction(Tlist))
        for freq in self.frequencies:
            x = freq / (0.695039 * Tlist)	# kB = 0.695039 cm^-1/K
            exp_x = numpy.exp(x)
            S = S + x / (exp_x - 1)
        return S * constants.R

    def getDensityOfStates(self, Elist, rho0=None):
        """
        Return the density of states at the specified energies `Elist` in J/mol
        above the ground state. The Beyer-Swinehart method is used to 
        efficiently convolve the vibrational density of states into the
        density of states of other modes. To be accurate, this requires a small
        (:math:`1-10 \\ \\mathrm{cm^{-1}}` or so) energy spacing.
        """
        cython.declare(rho=numpy.ndarray, freq=cython.double)
        cython.declare(dE=cython.double, nE=cython.int, dn=cython.int, n=cython.int)
        if rho0 is not None:
            rho = rho0
        else:
            rho = numpy.zeros_like(Elist)
        dE = Elist[1] - Elist[0]
        nE = len(Elist)
        for freq in self.frequencies:
            dn = int(freq * constants.h * constants.c * 100 * constants.Na / dE)
            for n in range(dn+1, nE):
                rho[n] = rho[n] + rho[n-dn]
        return rho

################################################################################

class StatesModel:
    """
    A set of molecular degrees of freedom data for a given molecule, comprising
    the results of a quantum chemistry calculation.

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `modes`             ``list``            A list of the degrees of freedom
    `spinMultiplicity`  ``int``             The spin multiplicity of the molecule
    =================== =================== ====================================

    """

    def __init__(self, modes=None, spinMultiplicity=1):
        self.modes = modes or []
        self.spinMultiplicity = spinMultiplicity

    def getHeatCapacity(self, Tlist):
        """
        Return the constant-pressure heat capacity in J/mol*K at the specified
        temperatures `Tlist` in K.
        """
        cython.declare(Cp=numpy.ndarray)
        Cp = numpy.ones_like(Tlist) * constants.R
        for mode in self.modes:
            Cp += mode.getHeatCapacity(Tlist)
        return Cp

    def getEnthalpy(self, Tlist):
        """
        Return the enthalpy in J/mol at the specified temperatures `Tlist` in K.
        """
        cython.declare(H=numpy.ndarray)
        H = constants.R * Tlist
        for mode in self.modes:
            H += mode.getEnthalpy(Tlist)
        return H

    def getEntropy(self, Tlist):
        """
        Return the entropy in J/mol*K at the specified temperatures `Tlist` in
        K.
        """
        cython.declare(S=numpy.ndarray)
        S = numpy.zeros_like(Tlist)
        for mode in self.modes:
            S += mode.getEntropy(Tlist)
        return S

    def getPartitionFunction(self, Tlist):
        """
        Return the the partition function at the specified temperatures
        `Tlist` in K. An active K-rotor is automatically included if there are
        no external rotational modes.
        """
        cython.declare(Q=numpy.ndarray, Trot=cython.double)
        Q = numpy.ones_like(Tlist)
        # Active K-rotor
        rotors = [mode for mode in self.modes if isinstance(mode, RigidRotor)]
        if len(rotors) == 0:
            Trot = 1.0 / constants.R / 3.141592654
            Q *= numpy.sqrt(Tlist / Trot)
        # Other modes
        for mode in self.modes:
            Q *= mode.getPartitionFunction(Tlist)
        return Q * self.spinMultiplicity

    def getDensityOfStates(self, Elist):
        """
        Return the value of the density of states in mol/J at the specified
        energies `Elist` in J/mol above the ground state. An active K-rotor is
        automatically included if there are no external rotational modes.
        """
        cython.declare(rho=numpy.ndarray, i=cython.int, E=cython.double)
        rho = numpy.zeros_like(Elist)
        # Active K-rotor
        rotors = [mode for mode in self.modes if isinstance(mode, RigidRotor)]
        if len(rotors) == 0:
            rho0 = numpy.zeros_like(Elist)
            for i, E in enumerate(Elist):
                if E > 0: rho0[i] = 1.0 / math.sqrt(1.0 * E)
            rho = convolve(rho, rho0, Elist)
        # Other non-vibrational modes
        for mode in self.modes:
            if not isinstance(mode, HarmonicOscillator):
                rho = convolve(rho, mode.getDensityOfStates(Elist), Elist)
        # Vibrational modes
        for mode in self.modes:
            if isinstance(mode, HarmonicOscillator):
                rho = mode.getDensityOfStates(Elist, rho)
        return rho * self.spinMultiplicity

def convolve(rho1, rho2, Elist):
    """
    Convolutes two density of states arrays `rho1` and `rho2` with corresponding
    energies `Elist` together using the equation

    .. math:: \\rho(E) = \\int_0^E \\rho_1(x) \\rho_2(E-x) \\, dx

    The units of the parameters do not matter so long as they are consistent.
    """

    cython.declare(rho=numpy.ndarray, found1=cython.bint, found2=cython.bint)
    cython.declare(dE=cython.double, nE=cython.int, i=cython.int, j=cython.int)
    rho = numpy.zeros_like(Elist)

    found1 = rho1.any(); found2 = rho2.any()
    if not found1 and not found2:
        pass
    elif found1 and not found2:
        rho = rho1
    elif not found1 and found2:
        rho = rho2
    else:
        dE = Elist[1] - Elist[0]
        nE = len(Elist)
        for i in range(nE):
            for j in range(i+1):
                rho[i] += rho2[i-j] * rho1[i] * dE

    return rho
