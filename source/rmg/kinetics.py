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

"""
Provides classes for working with the kinetics of chemical reactions:

* :class:`Kinetics` - A base class from which kinetics classes are derived

* :class:`ArrheniusKinetics` - A representation of an Arrhenius kinetic model

* :class:`ArrheniusEPKinetics` - A representation of an Arrhenius kinetic model with Evans-Polanyi correction

"""

import math
import quantities as pq

import constants
import data

################################################################################

class Kinetics:
	"""
	Represent a set of kinetic data. The details of the form of the kinetic
	data are left to a derived class. The attributes are:

	===============  ===========================================================
	Attribute        Description
	===============  ===========================================================
	`Trange`         A list of the minimum and maximum valid temperatures in K
	`rank`           An integer rank denoting the degree of confidence in the
	                 data (1 = high, 5 = low, 0 = none)
	`comment`        Comments, including but not limited to the source of the
	                 data
	`numReactants`   The number of reactants (used to determine the units of the
	                 kinetics)
	===============  ===========================================================

	"""

	def __init__(self, Trange=None, rank=0, comment=''):
		self.Trange = Trange
		self.rank = 0
		self.comment = ''
		self.numReactants = None

	def isTemperatureInRange(self, T):
		"""
		Return :data:`True` if temperature `T` is within the valid temperature
		range specified by self.Trange and :data:`False` if not. If no
		temperature range is specified in self.Trange, the kinetic data is
		assumed to be valid at all temperatures.
		"""
		if self.Trange is not None:
			if T < self.Trange[0] or T > self.Trange[1]:
				return False
		return True

################################################################################

class ArrheniusKinetics(Kinetics):
	"""
	Represent a set of modified Arrhenius kinetics. The kinetic expression has
	the form

	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_\\mathrm{a}}{RT} \\right)

	The attributes are:

	=========  ===========================================================
	Attribute  Description
	=========  ===========================================================
	`A`        The preexponential factor in s^-1, m^3/mol*s, etc.
	`n`        The temperature exponent
	`Ea`       The activation energy in J/mol
	=========  ===========================================================

	"""

	def __init__(self, A=0.0, Ea=0.0, n=0.0):
		"""If calling without keyword arguments be careful of the order!"""
		# in fact, should (can?) we enforce keyword calling?
		Kinetics.__init__(self)
		self.A = A
		self.Ea = Ea
		self.n = n

	def __str__(self):
		return 'k(T) = %s * T ** %s * math.exp(-%s / constants.R / T)\t%s < T < %s' % (self.A, self.n, self.Ea, self.Trange[0], self.Trange[1])

	def __repr__(self):
		"""How it looks on the console"""
		return '<ArrheniusKinetics A=%.0e E=%.0fkJ/mol n=%.1f >'%(self.A,
			self.Ea/1000.0, self.n )

	def getRateConstant(self, T):
		"""
		Return the rate constant k(T) at temperature `T` by evaluating the
		Arrhenius expression.
		"""

		# Raise exception if T is outside of valid temperature range
		#if not self.isTemperatureInRange(T):
		#	raise Exception('Attempted to evaluate a rate constant at an invalid temperature.')

		return self.A * (T ** self.n) * math.exp(-self.Ea / constants.R / T)

	def getReverse(self, dHrxn, Keq, T):
		"""
		Generate the reverse of the current kinetics for a reaction with
		standard enthalpy of reaction `dHrxn` and equilibrium constant `Keq` at
		298 K, respectively, defined in the same direction that these kinetics
		are.
		"""

		kinetics = ArrheniusKinetics(self.A / Keq * math.exp(-dHrxn / constants.R / T), self.Ea - dHrxn, self.n)
		kinetics.Trange = self.Trange
		kinetics.rank = self.rank
		kinetics.comment = self.comment

		return kinetics

	def toXML(self, dom, root, numReactants):
		"""
		Generate the XML for these kinetics using the :data:`xml.dom.minidom`
		package. The `dom` and `root` parameters refer to the DOM and the
		point within the DOM to place this item.
		"""
		kinetics = dom.createElement('arrheniusKinetics')
		root.appendChild(kinetics)
		kinetics.setAttribute('Trange', '%s-%s K' % (self.Trange[0], self.Trange[1]))
		kinetics.setAttribute('rank', str(self.rank))
		kinetics.setAttribute('comment', self.comment)

		preexponential = dom.createElement('preexponential')
		kinetics.appendChild(preexponential)
		preexponentialUnits = None
		if len(self.reactants) == 1:
			preexponentialUnits = 's^-1'
		else:
			preexponentialUnits = 'm^%s/(mol^%s*s)' % ((numReactants-1)*3, numReactants-1)
		data.createXMLQuantity(dom, preexponential, self.A, preexponentialUnits)

		exponent = dom.createElement('exponent')
		kinetics.appendChild(exponent)
		data.createXMLQuantity(dom, exponent, self.n, '')

		activationEnergy = dom.createElement('activationEnergy')
		kinetics.appendChild(activationEnergy)
		data.createXMLQuantity(dom, activationEnergy, self.Ea, 'J/mol')

################################################################################

class ArrheniusEPKinetics(Kinetics):
	"""
	Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
	kinetic expression has the form

	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_0 + \\alpha \\Delta H_\\mathrm{rxn}}{RT} \\right)

	The attributes are:

	=========  ===========================================================
	Attribute  Description
	=========  ===========================================================
	`A`        The preexponential factor in s^-1, m^3/mol*s, etc.
	`n`        The temperature exponent
	`E0`       The activation energy at zero enthalpy of reaction in J/mol
	`alpha`    The linear dependence of activation energy on enthalpy of
	           reaction
	=========  ===========================================================

	"""

	def __init__(self, A=0.0, E0=0.0, n=0.0, alpha=0.0):
		"""If calling without keyword arguments be careful of the order!"""
		Kinetics.__init__(self)
		self.A = A
		self.E0 = E0
		self.n = n
		self.alpha = alpha

	def __str__(self):
		return 'k(T) = %s * T ** %s * math.exp(-(%s + %s * DHrxn) / constants.R / T)\t%s < T < %s' % (self.A, self.n, self.E0, self.alpha, self.Trange[0], self.Trange[1])

	def __repr__(self):
		"""How it looks on the console"""
		return '<ArrheniusEPKinetics A=%.0e E0=%.0fkJ/mol n=%.1f alpha=%.1g>'%(
			self.A, self.E0/1000.0, self.n, self.alpha)

	def getActivationEnergy(self, dHrxn):
		"""
		Return the activation energy using the enthalpy of reaction `dHrxn`.
		"""
		return self.E0 + self.alpha * dHrxn

	def getArrhenius(self, dHrxn):
		"""
		Return the Arrhenius form of k(T) at temperature `T` by correcting E0
		to Ea using the enthalpy of reaction `dHrxn`.
		"""

		Ea = self.getActivationEnergy(dHrxn)

		kinetics = ArrheniusKinetics(self.A, Ea, self.n)
		kinetics.Trange = self.Trange
		kinetics.rank = self.rank
		kinetics.comment = self.comment

		return kinetics

	def getRateConstant(self, T, dHrxn):
		"""
		Return the rate constant k(T) at temperature `T` by evaluating the
		Arrhenius expression. The reaction has an enthalpy of reaction `dHrxn`.
		"""

		# Raise exception if T is outside of valid temperature range
		#if not self.isTemperatureInRange(T):
		#	raise Exception('Attempted to evaluate a rate constant at an invalid temperature.')

		Ea = self.getActivationEnergy(dHrxn)

		return self.A * T ** self.n * math.exp(-Ea / constants.R / T)

	def fromDatabase(self, data, comment, numReactants):
		"""
		Process a list of numbers `data` and associated description `comment`
		generated while reading from a kinetics database. The `numReactants`
		parameter is used to interpret the units of the preexponential.
		"""

		if len(data) != 11:
			raise Exception('Invalid list of kinetic data; should be a list of numbers of length 11.')

		Tmin, Tmax, A, n, alpha, E0, dA, dn, dalpha, dE0, rank = data

		originalUnits = 's^-1'; desiredUnits = 's^-1'
		if numReactants == 2:
			originalUnits = 'cm^3/(mol*s)'
		elif numReactants > 2:
			originalUnits = 'cm^%s/(mol^%s*s)' % ((numReactants-1)*3, numReactants-1)

		self.Trange = pq.Quantity([Tmin, Tmax], 'K').simplified
		self.Trange = [float(self.Trange[0]), float(self.Trange[1])]

		self.A = float(pq.Quantity(A, originalUnits).simplified)
		self.E0 = float(pq.Quantity(E0, 'kcal/mol').simplified)
		self.n = float(pq.Quantity(n, '').simplified)
		self.alpha = float(pq.Quantity(alpha, '').simplified)

		self.rank = rank
		self.comment = comment
		self.numReactants = numReactants

	def toXML(self, dom, root):
		"""
		Generate the XML for these kinetics using the :data:`xml.dom.minidom`
		package. The `dom` and `root` parameters refer to the DOM and the
		point within the DOM to place this item.
		"""
		kinetics = dom.createElement('arrheniusEPKinetics')
		root.appendChild(kinetics)
		kinetics.setAttribute('Trange', '%s-%s K' % (self.Trange[0], self.Trange[1]))
		kinetics.setAttribute('rank', str(self.rank))
		kinetics.setAttribute('comment', self.comment)


		preexponential = dom.createElement('preexponential')
		kinetics.appendChild(preexponential)
		preexponentialUnits = None
		numReactants = self.numReactants
		if numReactants == 1:
			preexponentialUnits = 's^-1'
		else:
			preexponentialUnits = 'm^%s/(mol^%s*s)' % ((numReactants-1)*3, numReactants-1)
		data.createXMLQuantity(dom, preexponential, self.A, preexponentialUnits)

		exponent = dom.createElement('exponent')
		kinetics.appendChild(exponent)
		data.createXMLQuantity(dom, exponent, self.n, '')

		alpha = dom.createElement('evansPolanyiSlope')
		kinetics.appendChild(alpha)
		data.createXMLQuantity(dom, alpha, self.alpha, '')

		activationEnergy = dom.createElement('evansPolanyiIntercept')
		kinetics.appendChild(activationEnergy)
		data.createXMLQuantity(dom, activationEnergy, self.E0, 'J/mol')

################################################################################

if __name__ == '__main__':
	
	pass
