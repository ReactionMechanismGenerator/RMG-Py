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
This module contains the molecular degrees of freedom model that are available
in ChemPy. All such models derive from the :class:`Mode` base class. A list of
molecular degrees of freedom can be stored in a :class:`StatesModel` object.
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
		return 'Translation(mass=%s, volume=%s, dimension=%s)' % (self.mass, self.volume, self.dimension)

	def getPartitionFunction(self, Tlist):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K. The formula is

		.. math:: q_\\mathrm{trans}(T, V) = \\left( \\frac{2 \\pi m k_\\mathrm{B} T}{h^2} \\right)^{d/2} V = q_\\mathrm{t} \\left( k_\\mathrm{B} T \\right)^{d/2}

		where :math:`T` is temperature, :math:`V` is volume, :math:`m` is mass,
		:math:`d` is dimensionality, :math:`k_\\mathrm{B}` is Boltzmann's
		cone-27stant, and :math:`h` is Planck's constant.
		"""
		Q = numpy.zeros_like(Tlist)
		qt = ((2 * constants.pi * self.mass / constants.Na) / (constants.h * constants.h))**(self.dimension/2.0) * self.volume
		Q = qt * (constants.kB * Tlist)**(self.dimension/2.0)
		return Q

	def getHeatCapacity(self, Tlist):
		"""
		Return the contribution to the heat capacity due to translation in 
		J/mol*K. The formula is

		.. math:: C_\\mathrm{v}^\\mathrm{trans}(T) = \\frac{d}{2} R

		where :math:`T` is temperature, :math:`V` is volume,
		:math:`d` is dimensionality, and :math:`R` is the gas law constant.
		"""
		return 0.5 * constants.R * self.dimension * numpy.ones_like(Tlist)
	
	def getEnthalpy(self, Tlist):
		"""
		Return the contribution to the enthalpy due to translation in 
		J/mol. The formula is

		.. math:: H^\\mathrm{trans}(T) = \\frac{d}{2} RT

		where :math:`T` is temperature, :math:`V` is volume,
		:math:`d` is dimensionality, and :math:`R` is the gas law constant.
		"""
		return 0.5 * self.dimension * Tlist
	
	def getEntropy(self, Tlist):
		"""
		Return the contribution to the entropy due to translation in 
		J/mol*K. The formula is

		.. math:: S^\\mathrm{trans}(T) = R \\ln \\left[ q_\\mathrm{t} \\left( k_\\mathrm{B} T \\right)^{d/2} e^{d/2} \\right]
		
		where :math:`T` is temperature, :math:`m` is mass, :math:`V` is volume,
		:math:`d` is dimensionality, and :math:`R` is the gas law constant.
		"""
		S = numpy.zeros_like(Tlist)
		qt = ((2 * constants.pi * self.mass / constants.Na) / (constants.h * constants.h))**(self.dimension/2.0) * self.volume
		S = constants.R * numpy.log( qt * (math.e * constants.R * Tlist)**(self.dimension/2.0) )
		return S
	
	def getDensityOfStates(self, Elist):
		"""
		Return the density of states at the specified energlies `Elist` in J/mol
		above the ground state. The formula is

		.. math:: \\rho(E) = \\frac{q_\\mathrm{t} E^{d/2-1}}{(d/2-1)!}

		where :math:`E` is energy, :math:`d` is dimensionality, and
		:math:`q_\\mathrm{t}` is defined in the equation for the partition
		function.
		"""
		rho = numpy.zeros_like(Elist)
		qt = ((2 * constants.pi * self.mass / constants.Na / constants.Na) / (constants.h * constants.h))**(self.dimension/2.0) * self.volume
		if self.dimension == 2: # Dimension is 2
			rho = qt
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
	:data:`False` if nonlinear. For a linear molecule, `frequencies` stores a
	list with one frequency, that of the rotation, in cm^-1. For a nonlinear
	molecule, `frequencies` stores a list of the three frequencies of rotation,
	even if two or three are equal, in cm^-1. Symmetry number corrections are
	*not* applied by this class.
	"""

	def __init__(self, linear=False, frequencies=None):
		self.linear = linear
		self.frequencies = frequencies or []

	def __repr__(self):
		return 'RigidRotor(linear=%s, frequencies=%s)' % (self.linear, self.frequencies)

	def getPartitionFunction(self, Tlist):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K. The formula is

		.. math:: q_\\mathrm{rot}(T) = \\frac{k_\\mathrm{B} T}{\\sigma h c \\tilde{\\omega}} = q_\\mathrm{r} k_\\mathrm{B} T

		for linear rotors and

		.. math:: q_\\mathrm{rot}(T) = \\frac{\\sqrt{\\pi}}{\\sigma} \\left[ \\frac{ \\left( k_\\mathrm{B} T \\right)^3 }{\\left( hc \\right)^3 \\tilde{\\omega_\\mathrm{A}} \\tilde{\\omega_\\mathrm{B}} \\tilde{\\omega_\\mathrm{C}} } \\right]^{1/2} = q_\\mathrm{r} \\left( k_\\mathrm{B} T \\right)^{3/2}

		for nonlinear rotors. Above, :math:`T` is temperature,
		:math:`\\tilde{\\omega}` is rotational frequency in cm^-1, :math:`c` is
		the speed of light,	:math:`k_\\mathrm{B}` is Boltzmann's constant, and
		:math:`h` is Planck's constant. :math:`\\sigma` is a placeholder for
		the symmetry number.
		"""
		Q = numpy.zeros_like(Tlist)
		if self.linear:
			theta = constants.h * constants.c * 100 * self.frequencies[0] / constants.kB
			Q = Tlist / theta
		else:
			theta = 1.0
			for freq in self.frequencies:
				theta = theta * constants.h * constants.c * 100 * freq / constants.kB
			Q = numpy.sqrt(constants.pi * Tlist**len(self.frequencies) / theta)
		return Q
		
	def getHeatCapacity(self, Tlist):
		"""
		Return the contribution to the heat capacity due to rigid rotation in
		J/mol*K. The formula is

		.. math:: C_\\mathrm{v}^\\mathrm{rot}(T) = \\frac{3}{2} R

		if nonlinear and

		.. math:: C_\\mathrm{v}^\\mathrm{rot}(T) = R

		if linear, where :math:`T` is temperature and :math:`R` is the gas law
		constant.
		"""
		Cv = constants.R * numpy.ones_like(Tlist)
		if not self.linear:
			Cv = 1.5 * Cv
		return Cv

	def getEnthalpy(self, Tlist):
		"""
		Return the contribution to the enthalpy due to rigid rotation in 
		J/mol. The formula is

		.. math:: H^\\mathrm{rot}(T) = \\frac{3}{2} RT

		where :math:`T` is temperature and :math:`R` is the gas law constant.
		"""
		return 1.5 * Tlist
	
	def getEntropy(self, Tlist):
		"""
		Return the contribution to the entropy due to rigid rotation in 
		J/mol*K. The formula is

		.. math:: S^\\mathrm{rot}(T) = R \\left( \\ln Q^\\mathrm{rot} + \\frac{3}{2} \\right)
		
		where :math:`Q^\\mathrm{rot}` is the partition function for a rigid 
		rotor and :math:`R` is the gas law constant.
		"""
		return constants.R * (numpy.log(self.getPartitionFunction(Tlist)) + 1.5)
	
	def getDensityOfStates(self, Elist):
		"""
		Return the density of states at the specified energlies `Elist` in J/mol
		above the ground state in mol/J. The formula is

		.. math:: \\rho(E) = q_\\mathrm{r}

		for linear rotors and

		.. math:: \\rho(E) = \\frac{q_\\mathrm{r} E^{1/2}}{\\frac{1}{2}!}

		for nonlinear rotors. Above,  :math:`E` is energy and
		:math:`q_\\mathrm{r}` is defined in the	equation for the partition
		function.
		"""
		rho = numpy.zeros_like(Elist)
		if self.linear:
			theta = constants.h * constants.c * 100 * self.frequencies[0] * constants.Na
			rho = 1.0 / theta
		else:
			theta = 1.0
			for freq in self.frequencies:
				theta = theta * constants.h * constants.c * 100 * freq * constants.Na
			rho = 2.0 * numpy.sqrt(Elist / theta)
		return rho

################################################################################

class HinderedRotor(Mode):
	"""
	A one-dimensional hindered rotor approximation of (internal) rotational
	modes. This class implements the Pitzer model of one-dimensional hindered
	rotation, which utilizes a hindered potential function

	.. math:: V(\\phi) = \\frac{1}{2} V_0 \\left[1 - \\cos \\left( \\sigma \\phi \\right) \\right]

	where :math:`V_0` is the height of the potential barrier and :math:`\\sigma`
	the number of minima or maxima in one revolution of angle :math:`\\phi`. The
	hindered rotor is therefore described by two quantities: the `frequency` of
	rotation in cm^-1, and the `barrier` height in cm^-1.
	"""

	def __init__(self, frequency=0.0, barrier=0.0):
		self.frequency = frequency
		self.barrier = barrier
	
	def __repr__(self):
		return '%s.HinderedRotor(%s, %s, %s)' % (self.__module__, self.frequency, self.barrier)

	def getPartitionFunction(self, Tlist):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K. The formula is

		.. math:: q_\\mathrm{hind}(T) = \\frac{1}{\\sigma} \\left( \\frac{\\pi k_\\mathrm{B} T}{h c \\tilde{\\omega}} \\right)^{1/2} \\exp \\left( -\\frac{V_0}{2 k_\\mathrm{B} T} \\right) I_0 \\left( \\frac{V_0}{2 k_\\mathrm{B} T} \\right) = q_\\mathrm{1f} \\left( k_\\mathrm{B} T \\right)^{1/2} \\exp \\left( -\\frac{V_0}{2 k_\\mathrm{B} T} \\right) I_0 \\left( \\frac{V_0}{2 k_\\mathrm{B} T} \\right)

		where :math:`T` is temperature, :math:`V_0` is the barrier height,
		:math:`\\tilde{\\omega}` is rotational frequency in cm^-1, :math:`c` is
		the speed of light,	:math:`k_\\mathrm{B}` is Boltzmann's constant, and
		:math:`h` is Planck's constant. :math:`\\sigma` is a placeholder for
		the symmetry number. :math:`I_0(x)` is the modified Bessel function of
		order zero.
		"""
		Q = []
		for T in Tlist:
			x = constants.h * constants.c * 100 * self.frequency / (constants.kB * T)
			z = 0.5 * constants.h * constants.c * 100 * self.barrier / (constants.kB * T)
			# The following is only valid in the classical limit
			Q.append(math.sqrt(2.0 * constants.pi * z) / x * math.exp(-z) * besseli0(z))
		return Q
	
	def getHeatCapacity(self, Tlist):
		"""
		Return the contribution to the heat capacity due to hindered rotation.
		The formula is

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{hind}(T)}{R} = \\frac{1}{2} + \\xi^2 - \\left[ \\xi \\frac{I_1(\\xi)}{I_0(\\xi)} \\right]^2 - \\xi \\frac{I_1(\\xi)}{I_0(\\xi)}

		where

		.. math:: \\xi \\equiv \\frac{V_0}{2 k_\\mathrm{B} T}

		:math:`T` is temperature, :math:`V_0` is barrier height,
		:math:`k_\\mathrm{B}` is Boltzmann's constant, and :math:`R` is the gas
		law constant.
		"""
		Cv = []
		for T in Tlist:
			x = constants.h * constants.c * 100 * self.frequency / (constants.kB * T)
			z = 0.5 * constants.h * constants.c * 100 * self.barrier / (constants.kB * T)
			
			exp_x = math.exp(x)
			one_minus_exp_x = 1.0 - exp_x
			BB = besseli1(z) / besseli0(z)
		
			Cv.append(x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - BB - z * BB * BB))
		
		return Cv
		
	def getDensityOfStates(self, Elist):
		"""
		Return the density of states at the specified energlies `Elist` in J/mol
		above the ground state. The formula is

		.. math:: \\rho(E) = \\frac{2 q_\\mathrm{1f}}{\\pi^{3/2} V_0^{1/2}} \\mathcal{K}(E / V_0) \\hspace{20pt} E < V_0

		and

		.. math:: \\rho(E) = \\frac{2 q_\\mathrm{1f}}{\\pi^{3/2} E^{1/2}} \\mathcal{K}(V_0 / E) \\hspace{20pt} E > V_0

		where :math:`E` is energy, :math:`V_0` is barrier height, and
		:math:`q_\\mathrm{1f}` is defined in the equation for the partition
		function. :math:`\\mathcal{K}(x)` is the complete elliptic integral of the first
		kind.
		"""
		
		# The following is only valid in the classical limit
		rho = []
		pre = 2.0 / constants.pi / (constants.h * constants.c * 100 * self.frequency * constants.Na)
		V0 = constants.h * constants.c * 100 * self.barrier * constants.Na
		for E in Elist:
			if E / V0 <= 1:
				rho.append(pre * cellipk(E / V0))
			else:
				rho.append(pre * math.sqrt(V0 / E) * cellipk(V0 / E))
		return rho

def besseli0(x):
	"""
	Return the value of the zeroth-order modified Bessel function at `x`.
	"""
	if math.abs(x) < 3.75:
		Y = x * x / 3.75 / 3.75
		return 1.0+Y*(3.5156229+Y*(3.0899424+Y*(1.2067429+Y*(0.2659732+Y*(0.0360768+Y*0.0045813)))))
	else:
		AX = math.abs(x)
		Y = 3.75 / AX
		BX = math.exp(AX) / math.sqrt(AX)
		AX = 0.39894228+Y*(0.01328592+Y*(0.00225319+Y*(-0.00157565+Y*(0.00916281+Y*(-0.02057706+Y*(0.02635537+Y*(-0.01647633+Y*0.00392377)))))))
		return AX * BX

def besseli1(x):
	"""
	Return the value of the first-order modified Bessel function at `x`.
	"""
	if math.abs(x) < 3.75:
		Y = x * x / 3.75 / 3.75
		return 0.5+Y*(0.87890594+Y*(0.51498869+Y*(0.15084934+Y*(0.02658733+Y*(0.00301532+Y*0.00032411)))))
	else:
		AX = math.abs(x)
		Y = 3.75 / AX
		BX = math.exp(AX) / math.sqrt(AX)
		AX = 0.39894228+Y*(-0.03988024+Y*(-0.00362018+Y*(0.00163801+Y*(-0.01031555+Y*(0.02282967+Y*(-0.02895312+Y*(0.01787654+Y*-0.00420059)))))))
		return AX * BX

def cellipk(x):
	"""
	Return the value of the complete elliptic integral of the first kind at `x`.
	"""
	if x < 0 or x > 1: return 0.0
	A = 1.0
	B = math.sqrt(1.0 - x)
	for n in range(100):
		A0 = A
		B0 = B
		A = (A0 + B0) / 2
		B = math.sqrt(A0 * B0)
		if math.abs(A - B) < 1.0e-12:
			break
			
	return 3.141592654 / 2.0 / A

################################################################################

class HarmonicOscillator(Mode):
	"""
	A representation of a vibrational mode as a one-dimensional quantum harmonic
	oscillator. The oscillator is defined by its `frequency` in cm^-1.
	"""

	def __init__(self, frequency=0.0):
		self.frequency = frequency
	
	def __repr__(self):
		return '%s.HarmonicOscillator(%s)' % (self.__module__, self.frequency)

	def getPartitionFunction(self, Tlist):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K. The formula is

		.. math:: q_\\mathrm{vib}(T) = \\frac{1.0}{1 - e^{- \\beta h c \\tilde{\\omega}}}

		where :math:`T` is temperature, :math:`\\beta \\equiv (k_\\mathrm{B} T)^{-1}`,
		:math:`\\tilde{\\omega}` is rotational frequency in cm^-1, :math:`c` is
		the speed of light,	:math:`k_\\mathrm{B}` is Boltzmann's constant, and
		:math:`h` is Planck's constant. Note that we have chosen our zero of
		energy to be at the zero-point energy of the molecule, *not* the bottom
		of the potential well.
		"""
		Q = []
		for T in Tlist:
			xi = self.frequency / (0.695039 * T)	# kB = 0.695039 cm^-1/K
			Q.append(1.0 / (1 - math.exp(-xi)))
		return Q
	
	def getHeatCapacity(self, Tlist):
		"""
		Return the contribution to the heat capacity due to hindered rotation.
		The formula is

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{vib}(T)}{R} = \\xi^2 \\frac{e^\\xi}{\\left( 1 - e^\\xi \\right)^2}

		where

		.. math:: \\xi \\equiv \\frac{h \\nu}{k_\\mathrm{B} T}

		:math:`T` is temperature, :math:`\\nu` is the vibration frequency,
		:math:`k_\\mathrm{B}` is Boltzmann's constant, and :math:`R` is the gas
		law constant.
		"""
		Cv = []
		for T in Tlist:
			x = self.frequency / (0.695039 * T)	# kB = 0.695039 cm^-1/K
			exp_x = math.exp(x)
			one_minus_exp_x = 1.0 - exp_x
			Cv.append(x * x * exp_x / one_minus_exp_x / one_minus_exp_x)
		return Cv

	def getDensityOfStates(self, Elist):
		"""
		Return the density of states at the specified energies `Elist` in J/mol
		above the ground state. The formula is

		.. math:: \\rho(E) = ?

		where :math:`E` is energy. Note that the Beyer-Swinehart algorithm
		provides a far more efficient method of convolving vibrational modes
		into a density of states expression, so this function should not be
		called for that purpose.
		"""
		pass

################################################################################

class StatesModel:
	"""
	A set of spectroscopic data for a given molecule. The degrees of freedom of
	the molecule are stored as a list in the `modes` attribute.
	"""

	def __init__(self, modes=None, symmetry=1):
		self.modes = modes or []
		self.symmetry = symmetry

	def getHeatCapacity(self, Tlist):
		"""
		Return the value of the heat capacity at the specified temperatures
		`Tlist` in K. The heat capacity returned is divided by the Boltzmann
		constant so as to be dimensionless.
		"""
		Cp = numpy.ones((len(Tlist)), numpy.float64)
		for mode in self.modes:
			Cp += mode.getHeatCapacity(Tlist)
		return Cp

	def getPartitionFunction(self, Tlist):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K.
		"""
		Q = numpy.ones((len(Tlist)), numpy.float64) / self.symmetry
		# Active K-rotor
		rotors = [mode for mode in self.modes if isinstance(mode, RigidRotor)]
		if len(rotors) == 0:
			Trot = 1.0 / constants.R / 3.141592654
			Q0 = [math.sqrt(T / Trot) for T in Tlist]
			for i in range(len(Tlist)):
				Q[i] *= Q0[i]
		# Other modes
		for mode in self.modes:
			Q0 = mode.getPartitionFunction(Tlist)
			for i in range(len(Tlist)):
				Q[i] *= Q0[i]
		return Q

	def getDensityOfStates(self, Elist, linear):
		"""
		Return the value of the density of states in mol/J at the specified
		energies `Elist` in J/mol above the ground state.
		"""
		
		import states
		
		# Prepare inputs for density of states function
		vib = numpy.array([mode.frequency for mode in self.modes if isinstance(mode, HarmonicOscillator)])
		rot = numpy.array([mode.frequencies for mode in self.modes if isinstance(mode, RigidRotor)])
		hind = numpy.array([[mode.frequency, mode.barrier] for mode in self.modes if isinstance(mode, HinderedRotor)])
		if len(hind) == 0: hind = numpy.zeros([0,2],numpy.float64)
		linear = 1 if linear else 0
		symm = self.symmetry

		# Calculate the density of states
		densStates, msg = states.densityofstates(Elist, vib, rot, hind, symm, linear)
		msg = msg.strip()
		if msg != '':
			raise Exception('Error while calculating the density of states for species %s: %s' % (self, msg))

		# Return result
		return densStates

	
