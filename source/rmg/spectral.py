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

import logging

import data
import constants
import thermo

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

	def partitionFunction(self, Tlist, V=1.0):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K. The formula is

		.. math:: q_\\mathrm{trans}(T, V) = \\left( \\frac{2 \\pi m k_\\mathrm{B} T}{h^2} \\right)^{d/2} V = q_\\mathrm{t} \\left( k_\\mathrm{B} T \\right)^{d/2}

		where :math:`T` is temperature, :math:`V` is volume, :math:`m` is mass,
		:math:`d` is dimensionality, :math:`k_\\mathrm{B}` is Boltzmann's
		constant, and :math:`h` is Planck's constant.
		"""

		Q = np.zeros(len(Tlist), np.float64)
		for i in range(len(Tlist)):
			Q[i] = ((2 * math.pi * self.mass / constants.Na * constants.kB * Tlist[i]) / (constants.h * constants.h))**(self.dimension/2.0) * V
		return Q

	def heatCapacity(self, Tlist, V=1.0):
		"""
		Return the contribution to the heat capacity due to translation. The 
		formula is

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{trans}(T)}{R} = \\frac{d}{2} R
				
		where :math:`T` is temperature, :math:`V` is volume,
		:math:`d` is dimensionality, and :math:`R` is the gas law constant.
		"""
		return np.ones(len(Tlist), numpy.float64) * self.dimension / 2.0

	def densityOfStates(self, Elist, V=1.0):
		"""
		Return the density of states at the specified energlies `Elist` in J/mol
		above the ground state. The formula is

		.. math:: \\rho(E) = \\frac{q_\\mathrm{t} E^{d/2-1}}{(d/2-1)!}

		where :math:`E` is energy, :math:`d` is dimensionality, and
		:math:`q_\\mathrm{t}` is defined in the equation for the partition
		function.
		"""

		qt = ((2 * math.pi * self.mass) / (constants.h * constants.h))**(self.dimension/2.0) * V

		rho = np.zeros(len(Elist), np.float64)
		if self.dimension == 2:
			for i in range(len(Elist)):
				rho[i] = qt
		else:
			for i in range(len(Elist)):
				rho[i] = qt * math.sqrt(Elist[i]) / (math.sqrt(math.pi) / 2)		# qt * E^0.5 / 0.5!

		return rho


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

	def partitionFunction(self, Tlist):
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

		Q = np.zeros(len(Tlist), np.float64)
		if self.linear:
			freqProduct = constants.h * constants.c * 100.0 * self.frequencies[0]
			for i in range(len(Tlist)):
				Q[i] = constants.kB * Tlist[i] / freqProduct
		else:
			freqProduct = 1.0
			for freq in self.frequencies:
				freqProduct *= freq * constants.h * constants.c * 100.0
			for i in range(len(Tlist)):
				Q[i] = math.sqrt(math.pi) * (constants.kB * Tlist[i])**(len(self.frequencies)/2.0) / math.sqrt(freqProduct)

		return Q

	def heatCapacity(self, Tlist):
		"""
		Return the contribution to the heat capacity due to rigid rotation. The
		formula is

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{rot}(T)}{R} = \\frac{3}{2}

		if nonlinear and

		.. math:: \\frac{C_\\mathrm{v}^\\mathrm{rot}(T)}{R} = 1

		if linear, where :math:`T` is temperature and :math:`R` is the gas law
		constant.
		"""
		return np.ones(len(Tlist), numpy.float64) * (1.0 if self.linear else 1.5)

	def densityOfStates(self, Elist):
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

		rho = np.zeros(len(Elist), np.float64)
		if self.linear:
			freqProduct = constants.h * constants.c * 100.0 * self.frequencies[0] * constants.Na
			for i in range(len(Elist)):
				rho[i] = 1.0 / freqProduct
		else:
			freqProduct = 1.0
			for freq in self.frequencies:
				freqProduct *= freq * constants.h * constants.c * 100.0
			qr = math.sqrt(math.pi) / math.sqrt(freqProduct)
			for i in range(len(Elist)):
				rho[i] = qr  * math.sqrt(Elist[i]) / (math.sqrt(math.pi) / 2)

		return rho

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

	def __init__(self, frequency=0.0, barrier=0.0):
		self.frequency = frequency
		self.barrier = barrier

	def __repr__(self):
		return '%s.HinderedRotor(%s, %s)' % (self.__module__, self.frequency, self.barrier)

	def partitionFunction(self, Tlist):
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

		import scipy.special

		q1f = math.sqrt(math.pi / (constants.h * constants.c * 100.0 * self.frequency))

		Q = np.zeros(len(Tlist), np.float64)
		for i in range(len(Tlist)):
			b = self.barrier * constants.h * constants.c * 100.0 / (2 * constants.kB * Tlist[i])
			Q[i] = q1f * math.sqrt(constants.kB * Tlist[i]) * math.exp(-b) * scipy.special.i0(b)

		return Q

	def heatCapacity(self, Tlist):
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
		import scipy.special
		xi = self.barrier * constants.h * constants.c * 100.0 / (2.0 * constants.kB * Tlist)
		I0 = scipy.special.i0(xi)
		I1 = scipy.special.i1(xi)
		return (0.5 + xi*xi - xi*xi*I1*I1/I0/I0 - xi*I1/I0)

	def densityOfStates(self, Elist):
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

		import scipy.special

		q1f = math.sqrt(math.pi / (constants.h * constants.c * 100.0 * self.frequency * constants.Na))
		V0 = self.barrier * constants.h * constants.c * 100.0 * constants.Na

		rho = np.zeros(len(Elist), np.float64)
		for i in range(len(Elist)):
			E = Elist[i]
			if Elist[i] < V0:
				rho[i] = 2 * q1f / (math.pi**1.5 * math.sqrt(V0)) * scipy.special.ellipk(E / V0)
			else:
				rho[i] = 2 * q1f / (math.pi**1.5 * math.sqrt(E)) * scipy.special.ellipk(V0 / E)

		return rho

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

	def toXML(self, document, rootElement):
		"""
		Add a <hinderedRotor> element as a child of `rootElement` using
		RMG-style XML. `document` is an :class:`io.XML` class representing the
		XML DOM tree.
		"""
		hinderedRotorElement = document.createElement('hinderedRotor', rootElement)
		document.createQuantity('frequency', hinderedRotorElement, self.frequency, 'cm^-1')
		document.createQuantity('barrier', hinderedRotorElement, self.barrier, 'cm^-1')

################################################################################

class HarmonicOscillator:
	"""
	A representation of a vibrational mode as a one-dimensional quantum harmonic
	oscillator. The oscillator is defined by its `frequency` in cm^-1.
	"""

	def __init__(self, frequency=0.0):
		self.frequency = frequency

	def __repr__(self):
		return '%s.HarmonicOscillator(%s)' % (self.__module__, self.frequency)

	def partitionFunction(self, Tlist):
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

		Q = np.zeros(len(Tlist), np.float64)
		for i in range(len(Tlist)):
			exponent = constants.h * constants.c * 100.0 * self.frequency / (constants.kB * Tlist[i])
			Q[i] = 1.0 / (1 - math.exp(-exponent))
		return Q

	def heatCapacity(self, Tlist):
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
		from numpy import exp
		xi = constants.h * constants.c * 100.0 * self.frequency / (constants.kB * Tlist)
		exp_xi = exp(xi)
		return (xi*xi * exp_xi / (1 - exp_xi) / (1 - exp_xi))

	def densityOfStates(self, Elist):
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

	def toXML(self, document, rootElement):
		"""
		Add a <harmonicOscillator> element as a child of `rootElement` using
		RMG-style XML. `document` is an :class:`io.XML` class representing the
		XML DOM tree.
		"""
		harmonicOscillatorElement = document.createElement('harmonicOscillator', rootElement)
		document.createQuantity('frequency', harmonicOscillatorElement, self.frequency, 'cm^-1')

################################################################################

class SpectralData:
	"""
	A set of spectroscopic data for a given molecule. The degrees of freedom of
	the molecule are stored as a list in the `modes` attribute.
	"""

	def __init__(self, modes=None, symmetry=1):
		self.modes = modes or []
		self.symmetry = symmetry

	def partitionFunction(self, Tlist):
		"""
		Return the value of the partition function at the specified temperatures
		`Tlist` in K.
		"""
		Q = np.ones((len(Tlist)), np.float64) / self.symmetry
		# Active K-rotor
		rotors = [mode for mode in self.modes if isinstance(mode, RigidRotor)]
		if len(rotors) == 0:
			Trot = constants.h * constants.c * 100.0 * 1.0 / constants.kB
			Q0 = [math.sqrt(T / Trot) for T in Tlist]
			for i in range(len(Tlist)):
				Q[i] *= Q0[i]
		# Other modes
		for mode in self.modes:
			Q0 = mode.partitionFunction(Tlist)
			for i in range(len(Tlist)):
				Q[i] *= Q0[i]
		return Q

	def densityOfStates(self, Elist):
		"""
		Return the value of the density of states in mol/J at the specified
		energies `Elist` in J/mol above the ground state.
		"""

		rho = np.zeros((len(Elist)), np.float64)

		if len(self.modes) == 0:
			return rho

		# First convolve nonvibrational modes
		for mode in self.modes:
			if not isinstance(mode, HarmonicOscillator):
				rho0 = mode.densityOfStates(Elist)
				if not rho.any():
					rho = rho0
				else:
					rho = convolve(rho, rho0, Elist)
		# Use Beyer-Swinehart for vibrational modes
		if not rho.any():
			rho[0] = 1.0 / (Elist[1] - Elist[0])
		rho = beyerSwinehart(Elist, [mode.frequency for mode in self.modes if isinstance(mode, HarmonicOscillator)], rho)
		# Return result
		return smooth(rho)

	def phi(self, beta, E):
		beta = float(beta)
		T = 1.0 / (constants.R * beta)
		Q = self.partitionFunction([T])[0]
		return math.log(Q) + beta * float(E)

	def densityOfStatesForst(self, Elist):

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

class CharacteristicFrequency:
	"""
	Represent a characteristic frequency in the frequency database. The
	characteristic frequency has a real lower bound `lower`, a real upper bound
	`upper`, an integer `degeneracy`.
	"""

	def __init__(self, lower=0.0, upper=0.0, degeneracy=1):
		self.lower = lower
		self.upper = upper
		self.degeneracy = degeneracy

	def generateFrequencies(self, count=1):
		"""
		Generate a set of frequencies. The number of frequencies returned is
		self.degeneracy * count, and these are distributed linearly between
		`self.lower` and `self.upper`.
		"""
		from numpy import linspace
		number = self.degeneracy * count
		if number == 1:
			return [(self.lower + self.upper) / 2.0]
		else:
			return linspace(self.lower, self.upper, number, endpoint=True)

################################################################################

class FrequencyDatabase(data.Database):
	"""
	Represent an RMG frequency database.
	"""

	def __init__(self):
		"""
		Call the generic `data.Database.__init__()` method. This in turn creates
			* self.dictionary = Dictionary()
			* self.library = Library()
			* self.tree = Tree()
		"""
		data.Database.__init__(self)


	def load(self, dictstr, treestr, libstr):
		"""
		Load a group frequency database. The database is stored
		in three files: `dictstr` is the path to the dictionary, `treestr` to
		the tree, and `libstr` to the library. The tree is optional, and should
		be set to '' if not desired.
		"""

		# Load dictionary, library, and (optionally) tree
		data.Database.load(self, dictstr, treestr, libstr)

		# Convert data in library to ThermoData objects or lists of
		# [link, comment] pairs
		for label, item in self.library.iteritems():

			if item is None:
				pass
			elif not item.__class__ is tuple:
				raise data.InvalidDatabaseException('Frequencies library should be tuple at this point. Instead got %r'%data)
			else:
				index, item = item
				if not (item.__class__ is str or item.__class__ is unicode):
					raise data.InvalidDatabaseException('Frequencies library data format is unrecognized.')

				items = item.split()

				# First item is the symmetry correction 
				frequencies = [int(items.pop(0))]

				# Items should be a multiple of three (no comments allowed at the moment)
				if len(items) % 3 != 0:
					raise data.InvalidDatabaseException('Unexpected number of items encountered in frequencies library.')

				# Convert list of data into a list of characteristic frequencies
				count = len(items) / 3
				for i in range(count):
					frequency = CharacteristicFrequency(
						lower=float(items[3*i+1]),
						upper=float(items[3*i+2]),
						degeneracy=int(items[3*i]))
					frequencies.append(frequency)

				self.library[label] = frequencies

		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (dictstr))

frequencyDatabase = None

################################################################################

def loadFrequencyDatabase(databasePath):
	"""
	Create and load the frequencies databases.
	"""
	import os.path
	import logging

	databasePath = os.path.join(databasePath, 'frequencies')

	# Create and load thermo databases
	database = FrequencyDatabase()
	logging.debug('\tFrequencies database')
	database.load(
		dictstr=os.path.join(databasePath, 'Dictionary.txt'),
		treestr=os.path.join(databasePath, 'Tree.txt'),
		libstr=os.path.join(databasePath, 'Library.txt'))

	return database

################################################################################

def generateSpectralData(struct, thermoData):
	"""
	Generate the spectral data for a :class:`structure.Structure` object
	`struct` with corresponding :class:`thermo.ThermoGAData` object `thermo`.
	The group frequency method is used to do so; this method has two steps:

	1.	Search the structure for certain functional groups for which
		characteristic frequencies are known, and use those frequencies.
	2.	For any remaining degrees of freedom, fit the parameters such that they
		replicate the heat capacity.

	This method only fits the internal modes (i.e. vibrations and hindered
	rotors).
	"""

	# No spectral data for single atoms
	if len(struct.atoms()) < 2:
		return None

	# Determine linearity of structure
	linear = struct.isLinear()

	# Determine the number of internal modes for this structure
	numModes = 3 * len(struct.atoms()) - (5 if linear else 6)

	# Determine the number of vibrational and hindered rotor modes for this
	# structure
	numRotors = struct.countInternalRotors()
	numVibrations = numModes - numRotors
	#print 'For %s, I found %i internal rotors and %i vibrations for a total of %i modes' % (struct, numRotors, numVibrations, numModes)

	# For each group in library, find all subgraph isomorphisms
	groupCount = {}
	for node, data in frequencyDatabase.library.iteritems():
		ismatch, map12List, map21List = struct.findSubgraphIsomorphisms(frequencyDatabase.dictionary[node])
		if ismatch:
			count = len(map12List)
		else:
			count = 0
		if count % data[0] != 0:
			raise Exception('Incorrect number of matches of node "%s" while estimating frequencies of %s.' % (node, struct))
		groupCount[node] = count / data[0]

	# For debugging, print a list of the groups found
	#print 'Groups found:'
	#for node, count in groupCount.iteritems():
	#	if count != 0: print '\t', node, count

	# Get characteristic frequencies
	frequencies = []
	for node, count in groupCount.iteritems():
		for charFreq in frequencyDatabase.library[node][1:]:
			frequencies.extend(charFreq.generateFrequencies(count))

	# Check that we have the right number of degrees of freedom specified
	if len(frequencies) > numVibrations:
		# We have too many vibrational modes
		difference = len(frequencies) - numVibrations
		# First try to remove hindered rotor modes until the proper number of modes remains
		if numRotors >= difference:
			numRotors -= difference
			logging.warning('For structure %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i internal rotors to compensate.' % (struct, difference))
		# If that doesn't work, turn off functional groups until the problem is underspecified again
		else:
			groupsRemoved = 0
			freqsRemoved = 0
			freqCount = len(frequencies)
			while freqCount > numVibrations:
				minDegeneracy, minNode = min([(sum([charFreq.degeneracy for charFreq in frequencyDatabase.library[node][1:]]), node) for node in groupCount])
				if groupCount[minNode] > 1:
					groupCount[minNode] -= 1
				else:
					del groupCount[minNode]
				groupsRemoved += 1
				freqsRemoved += minDegeneracy
				freqCount -= minDegeneracy
			# Log warning
			logging.warning('For structure %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i groups (%i frequencies) to compensate.' % (struct, groupsRemoved, freqsRemoved))
			# Regenerate characteristic frequencies
			frequencies = []
			for node, count in groupCount.iteritems():
				for charFreq in frequencyDatabase.library[node][1:]:
					frequencies.extend(charFreq.generateFrequencies(count))

	# Create spectral data object with characteristic frequencies
	spectralData = SpectralData()
	for freq in frequencies:
		spectralData.modes.append(HarmonicOscillator(frequency=freq))

	# Subtract out contributions to heat capacity from the characteristic modes
	import numpy
	Cp = [thermoData.getHeatCapacity(T) for T in thermo.ThermoGAData.CpTlist]
	Cv = numpy.array(Cp) / constants.R
	Tlist = numpy.array(thermo.ThermoGAData.CpTlist)
	for mode in spectralData.modes:
		Cv -= mode.heatCapacity(Tlist)
	# Subtract out translational modes
	Cv -= 1.5
	# Subtract out external rotational modes
	Cv -= (1.5 if not linear else 1.0)
	# Subtract out PV term (Cp -> Cv)
	Cv -= 1.0
	# Check that all Cv values are still positive (should we do this?)
	#for C in Cv:
	#	if C <= 0.0: raise Exception('Remaining heat capacity is negative.')
	
	# Fit remaining frequencies and hindered rotors to the heat capacity data
	import spectralfit
	freeVibrations = numVibrations - len(frequencies)
	if freeVibrations > 0 and numRotors > 0:
		vib, hind = spectralfit.fitspectraldata(Cv, Tlist, freeVibrations, numRotors)
		for v in vib:
			spectralData.modes.append(HarmonicOscillator(frequency=v))
		for v, b in hind:
			spectralData.modes.append(HinderedRotor(frequency=v, barrier=b))
	elif freeVibrations > 0 and numRotors == 0:
		vib = spectralfit.fitspectraldatanorotors(Cv, Tlist, freeVibrations)
		for v in vib:
			spectralData.modes.append(HarmonicOscillator(frequency=v))
	elif freeVibrations == 0 and numRotors > 0:
		hind = spectralfit.fitspectraldatanooscillators(Cv, Tlist, numRotors)
		for v, b in hind:
			spectralData.modes.append(HinderedRotor(frequency=v, barrier=b))

	return spectralData

################################################################################
