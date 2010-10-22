#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

import numpy as np
import math

import rmg.constants as constants

import _modes

################################################################################

class Translation:
	"""
	A representation of translational modes. The `dimension` attribute is an
	integer representing the dimension of the translation (2 or 3) and the
	`mass` is the molar mass of the molecule in kg/mol.
	"""

	def __init__(self, mass=0.0, dimension=3):
		self.mass = mass
		self.dimension = dimension

	def __repr__(self):
		return '%s.Translation(%s, %s)' % (self.__module__, self.mass, self.dimension)

	def getPartitionFunction(self, Tlist, V=1.0):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K. The formula is

		.. math:: q_\\mathrm{trans}(T, V) = \\left( \\frac{2 \\pi m k_\\mathrm{B} T}{h^2} \\right)^{d/2} V = q_\\mathrm{t} \\left( k_\\mathrm{B} T \\right)^{d/2}

		where :math:`T` is temperature, :math:`V` is volume, :math:`m` is mass,
		:math:`d` is dimensionality, :math:`k_\\mathrm{B}` is Boltzmann's
		constant, and :math:`h` is Planck's constant.
		"""
		return _modes.translation_partitionfunction(Tlist, self.mass, self.dimension, V)

	def getHeatCapacity(self, Tlist, V=1.0):
		"""
		Return the contribution to the heat capacity due to translation. The
		formula is

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{trans}(T)}{R} = \\frac{d}{2} R

		where :math:`T` is temperature, :math:`V` is volume,
		:math:`d` is dimensionality, and :math:`R` is the gas law constant.
		"""
		return _modes.translation_heatcapacity(Tlist, self.mass, self.dimension, V)

	def getDensityOfStates(self, Elist, V=1.0):
		"""
		Return the density of states at the specified energlies `Elist` in J/mol
		above the ground state. The formula is

		.. math:: \\rho(E) = \\frac{q_\\mathrm{t} E^{d/2-1}}{(d/2-1)!}

		where :math:`E` is energy, :math:`d` is dimensionality, and
		:math:`q_\\mathrm{t}` is defined in the equation for the partition
		function.
		"""
		return _modes.translation_densityofstates(Elist, self.mass, self.dimension, V)

################################################################################

class RigidRotor:
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
		return '%s.RigidRotor(%s, %s)' % (self.__module__, self.linear, self.frequencies)

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
		return _modes.freerotor_partitionfunction(Tlist, self.frequencies, 1 if self.linear else 0)

	def getHeatCapacity(self, Tlist):
		"""
		Return the contribution to the heat capacity due to rigid rotation. The
		formula is

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{rot}(T)}{R} = \\frac{3}{2}

		if nonlinear and

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{rot}(T)}{R} = 1

		if linear, where :math:`T` is temperature and :math:`R` is the gas law
		constant.
		"""
		return _modes.freerotor_heatcapacity(Tlist, self.frequencies, 1 if self.linear else 0)

	def getDensityOfStates(self, Elist):
		"""
		Return the density of states at the specified energlies `Elist` in J/mol
		above the ground state. The formula is

		.. math:: \\rho(E) = q_\\mathrm{r}

		for linear rotors and

		.. math:: \\rho(E) = \\frac{q_\\mathrm{r} E^{1/2}}{\\frac{1}{2}!}

		for nonlinear rotors. Above,  :math:`E` is energy and
		:math:`q_\\mathrm{r}` is defined in the	equation for the partition
		function.
		"""
		return _modes.freerotor_densityofstates(Elist, self.frequencies, 1 if self.linear else 0)

	def fromXML(self, document, rootElement):
		"""
		Convert a <rigidRotor> element from a standard RMG-style XML
		input file into a RigidRotor object. `document` is an
		:class:`io.XML` class representing the XML DOM tree, and `rootElement`
		is the <rigidRotor> element in that tree.
		"""

		# Read <frequencies> element
		self.frequencies = document.getChildQuantity(rootElement, 'frequencies', required=True)
		self.frequencies = [float(f) for f in self.frequency]

		# Read <linear> attribute
		self.linear = str(document.getAttribute(rootElement, 'linear', required=False, default='no')).lower()
		self.linear = (self.linear == 'yes' or self.linear == 'y' or self.linear == 'true')

	def toXML(self, document, rootElement):
		"""
		Add a <rigidRotor> element as a child of `rootElement` using
		RMG-style XML. `document` is an :class:`io.XML` class representing the
		XML DOM tree.
		"""

		rigidRotorElement = document.createElement('rigidRotor', rootElement)
		linear = 'yes' if self.linear else 'no'
		document.createAttribute('linear', rigidRotorElement, linear)
		document.createQuantity('frequencies', rigidRotorElement, self.frequencies, 'cm^-1')

################################################################################

class HinderedRotor:
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

	def __init__(self, frequency=0.0, barrier=0.0, degeneracy=1):
		self.frequency = frequency
		self.barrier = barrier
		self.degeneracy = degeneracy

	def __repr__(self):
		return '%s.HinderedRotor(%s, %s, %s)' % (self.__module__, self.frequency, self.barrier, self.degeneracy)

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
		return _modes.hinderedrotor_partitionfunction(Tlist, self.frequency, self.barrier) ** self.degeneracy

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
		return _modes.hinderedrotor_heatcapacity(Tlist, self.frequency, self.barrier) * self.degeneracy

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
		return _modes.hinderedrotor_densityofstates(Elist, self.frequency, self.barrier)

	def fromXML(self, document, rootElement, frequencyScaleFactor=1.0):
		"""
		Convert a <hinderedRotor> element from a standard RMG-style XML
		input file into a HinderedRotor object. `document` is an
		:class:`io.XML` class representing the XML DOM tree, and `rootElement`
		is the <hinderedRotor> element in that tree.
		"""

		# Read <frequency> element
		self.frequency = document.getChildQuantity(rootElement, 'frequency', required=True)
		self.frequency = float(self.frequency) * frequencyScaleFactor

		# Read <barrier> element
		self.barrier = document.getChildQuantity(rootElement, 'barrier', required=True)
		self.barrier = float(self.barrier)

		# Read <degeneracy> element
		self.degeneracy = int(document.getChildElementText(rootElement, 'degeneracy', required=False, default='1'))

	def toXML(self, document, rootElement):
		"""
		Add a <hinderedRotor> element as a child of `rootElement` using
		RMG-style XML. `document` is an :class:`io.XML` class representing the
		XML DOM tree.
		"""
		hinderedRotorElement = document.createElement('hinderedRotor', rootElement)
		document.createQuantity('frequency', hinderedRotorElement, self.frequency, 'cm^-1')
		document.createQuantity('barrier', hinderedRotorElement, self.barrier, 'cm^-1')
		document.createTextElement('degeneracy', hinderedRotorElement, str(self.degeneracy))

################################################################################

class HarmonicOscillator:
	"""
	A representation of a vibrational mode as a one-dimensional quantum harmonic
	oscillator. The oscillator is defined by its `frequency` in cm^-1.
	"""

	def __init__(self, frequency=0.0, degeneracy=1):
		self.frequency = frequency
		self.degeneracy = degeneracy

	def __repr__(self):
		return '%s.HarmonicOscillator(%s, %s)' % (self.__module__, self.frequency, self.degeneracy)

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
		return _modes.harmonicoscillator_partitionfunction(Tlist, self.frequency) ** self.degeneracy

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
		return _modes.harmonicoscillator_heatcapacity(Tlist, self.frequency) * self.degeneracy

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

	def fromXML(self, document, rootElement, frequencyScaleFactor=1.0):
		"""
		Convert a <harmonicOscillator> element from a standard RMG-style XML
		input file into a HarmonicOscillator object. `document` is an
		:class:`io.XML` class representing the XML DOM tree, and `rootElement`
		is the <harmonicOscillator> element in that tree.
		"""

		# Read <frequency> element
		self.frequency = document.getChildQuantity(rootElement, 'frequency', required=True)
		self.frequency = float(self.frequency) * frequencyScaleFactor

		# Read <degeneracy> element
		self.degeneracy = int(document.getChildElementText(rootElement, 'degeneracy', required=False, default='1'))
	
	def toXML(self, document, rootElement):
		"""
		Add a <harmonicOscillator> element as a child of `rootElement` using
		RMG-style XML. `document` is an :class:`io.XML` class representing the
		XML DOM tree.
		"""
		harmonicOscillatorElement = document.createElement('harmonicOscillator', rootElement)
		document.createQuantity('frequency', harmonicOscillatorElement, self.frequency, 'cm^-1')
		document.createTextElement('degeneracy', harmonicOscillatorElement, str(self.degeneracy))

################################################################################

class SpectralData:
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
		Cp = np.ones((len(Tlist)), np.float64)
		for mode in self.modes:
			Cp += mode.getHeatCapacity(Tlist)
		return Cp

	def getPartitionFunction(self, Tlist):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K.
		"""
		Q = np.ones((len(Tlist)), np.float64) / self.symmetry
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
		vib = np.array([mode.frequency for mode in self.modes if isinstance(mode, HarmonicOscillator)])
		rot = np.array([mode.frequencies for mode in self.modes if isinstance(mode, RigidRotor)])
		hind = np.array([[mode.frequency, mode.barrier] for mode in self.modes if isinstance(mode, HinderedRotor)])
		if len(hind) == 0: hind = np.zeros([0,2],np.float64)
		linear = 1 if linear else 0
		symm = self.symmetry

		# Calculate the density of states
		densStates, msg = states.densityofstates(Elist, vib, rot, hind, symm, linear)
		msg = msg.strip()
		if msg != '':
			raise Exception('Error while calculating the density of states for species %s: %s' % (self, msg))

		# Return result
		return densStates

	def phi(self, beta, E):
		beta = float(beta)
		T = 1.0 / (constants.R * beta)
		Q = self.partitionFunction([T])[0]
		return math.log(Q) + beta * float(E)

	def getDensityOfStatesForst(self, Elist):

		import scipy.optimize

		# Initialize density of states array
		rho = [0.0]

		# Initial guess for first minimization
		x = 1e-5

		for E in Elist[1:]:

			# Find minimum of phi

			#                       func x0 arg xtol  ftol  maxi maxf  fullout disp retall  callback
			x = scipy.optimize.fmin(self.phi, x, (), 1e-8, 1e-8, 100, 1000, False, False, False, None)
			x = float(x)
			dx = 1e-4 * x

			# Apply first-order steepest descents approximation (accurate to 1-3%, smoother)
			f = self.phi(x)
			d2fdx2 = (self.phi(x+dx) - 2 * self.phi(x) + self.phi(x-dx)) / (dx**2)
			i1 = math.exp(f) / math.sqrt(2 * math.pi * d2fdx2)

			# Apply second-order steepest descents approximation (more accurate, less smooth)
			#d3fdx3 = (self.phi(x+1.5*dx) - 3 * self.phi(x+0.5*dx) + 3 * self.phi(x-0.5*dx) - self.phi(x-1.5*dx)) / (dx**3)
			#d4fdx4 = (self.phi(x+2*dx) - 4 * self.phi(x+dx) + 6 * self.phi(x) - 4 * self.phi(x-dx) + self.phi(x-2*dx)) / (dx**4)
			#i2 = i1 * (1 + d4fdx4 / 8 / (d2fdx2**2) - 5 * (d3fdx3**2) / 24 / (d2fdx2**3))

			rho.append(i1)

		return rho

	def fromXML(self, document, rootElement):
		"""
		Convert a <spectralData> element from a standard RMG-style XML input
		file into a SpectralData object. `document` is an :class:`io.XML` class
		representing the XML DOM tree, and `rootElement` is the <spectralData>
		element in that tree.
		"""

		# Read <frequencyScaleFactor> element
		frequencyScaleFactor = float(document.getChildQuantity(rootElement, 'frequencyScaleFactor', required=False, default=1.0))

		# Read <harmonicOscillator> elements
		hoElements = document.getChildElements(rootElement, 'harmonicOscillator', required=False)
		for hoElement in hoElements:
			mode = HarmonicOscillator()
			mode.fromXML(document, hoElement, frequencyScaleFactor)
			self.modes.append(mode)

		# Read <hinderedRotor> elements
		hrElements = document.getChildElements(rootElement, 'hinderedRotor', required=False)
		for hrElement in hrElements:
			mode = HinderedRotor()
			mode.fromXML(document, hrElement, frequencyScaleFactor)
			self.modes.append(mode)

		# Read <rigidRotor> elements
		rrElements = document.getChildElements(rootElement, 'rigidRotor', required=False)
		for rrElement in rrElements:
			mode = RigidRotor()
			mode.fromXML(document, rrElement)
			self.modes.append(mode)

		# Read <symmetryNumber> element
		self.symmetry = float(document.getChildQuantity(rootElement, 'symmetryNumber', required=False, default=1))

	def toXML(self, document, rootElement):
		"""
		Add a <spectralData> element as a child of `rootElement` using
		RMG-style XML. `document` is an :class:`io.XML` class representing the
		XML DOM tree.
		"""

		spectralDataElement = document.createElement('spectralData', rootElement)
		for mode in self.modes:
			mode.toXML(document, spectralDataElement)
		document.createTextElement('symmetryNumber', spectralDataElement, str(self.symmetry))

################################################################################
