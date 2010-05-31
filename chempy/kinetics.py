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
This module contains the kinetics models that are available in ChemPy.
All such models derive from the :class:`KineticsModel` base class.
"""

################################################################################

import math
import cython

import constants
from exception import InvalidKineticsModelError

################################################################################

class KineticsModel:
	"""
	Represent a set of kinetic data. The details of the form of the kinetic
	data are left to a derived class. The attributes are:

	=============== =============== ============================================
	Attribute       Type            Description
	=============== =============== ============================================
	`Tmin`          :class:`float`  The minimum absolute temperature in K at which the model is valid
	`Tmax`          :class:`float`  The maximum absolute temperature in K at which the model is valid
	`Pmin`          :class:`float`  The minimum absolute pressure in Pa at which the model is valid
	`Pmax`          :class:`float`  The maximum absolute pressure in Pa at which the model is valid
	`numReactants`  :class:`int`    The number of reactants (used to determine the units of the kinetics)
	`comment`       :class:`str`    A string containing information about the model (e.g. its source)
	=============== =============== ============================================
	
	"""

	def __init__(self, Tmin=0.0, Tmax=1.0e10, Pmin=0.0, Pmax=1.0e100, numReactants=-1, comment=''):
		self.Tmin = Tmin
		self.Tmax = Tmax
		self.Pmin = Pmin
		self.Pmax = Pmax
		self.numReactants = numReactants
		self.comment = comment

	def isTemperatureValid(self, T):
		"""
		Return :data:`True` if temperature `T` is within the valid temperature
		range and :data:`False` if not. 
		"""
		return (self.Tmin <= T and T <= self.Tmax)

	def isPressureValid(self, P):
		"""
		Return :data:`True` if pressure `P` is within the valid pressure
		range, and :data:`False` if not.
		"""
		return (self.Pmin <= P and P <= self.Pmax)

################################################################################

class ArrheniusModel(KineticsModel):
	"""
	Represent a set of modified Arrhenius kinetics. The kinetic expression has
	the form

	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_\\mathrm{a}}{RT} \\right)

	The attributes are:

	=============== =============== ============================================
	Attribute       Type            Description
	=============== =============== ============================================
	`A`             :class:`float`  The preexponential factor in s^-1, m^3/mol*s, etc.
	`n`             :class:`float`  The temperature exponent
	`Ea`            :class:`float`  The activation energy in J/mol
	=============== =============== ============================================
	
	"""
	
	def __init__(self, A=0.0, Ea=0.0, n=0.0):
		KineticsModel.__init__(self)
		self.A = A
		self.Ea = Ea
		self.n = n
	
	def __str__(self):
		return 'k(T) = %s * T ** %s * math.exp(-%s / constants.R / T)\t%s < T < %s' % (self.A, self.n, self.Ea, self.Tmin, self.Tmax)
	
	def __repr__(self):
		return '<ArrheniusModel A=%.0e Ea=%.0f kJ/mol n=%.1f>' % (self.A,self.Ea/1000.0, self.n)
	
	def getRateConstant(self, T):
		"""
		Return the rate constant k(T) at temperature `T` in K by evaluating the
		Arrhenius expression.
		"""
		if not self.isTemperatureValid(T):
			raise InvalidKineticsModelError('Invalid temperature for rate coefficient evaluation from ArrheniusModel.')
		if cython.compiled:
			return self.A * pow(T, self.n) * exp(-self.Ea / constants.R / T)
		else:
			return self.A * (T ** self.n) * math.exp(-self.Ea / constants.R / T)

################################################################################

class ArrheniusEPModel(KineticsModel):
	"""
	Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
	kinetic expression has the form

	.. math:: k(T) = A T^n \\exp \\left( - \\frac{E_0 + \\alpha \\Delta H_\\mathrm{rxn}}{RT} \\right)

	The attributes are:

	=============== =============== ============================================
	Attribute       Type            Description
	=============== =============== ============================================
	`A`             :class:`float`  The preexponential factor in s^-1, m^3/mol*s, etc.
	`n`             :class:`float`  The temperature exponent
	`E0`            :class:`float`  The activation energy at zero enthalpy of reaction in J/mol
	`alpha`         :class:`float`  The linear dependence of activation energy on enthalpy of reaction
	=============== =============== ============================================
	
	"""

	def __init__(self, A=0.0, E0=0.0, n=0.0, alpha=0.0):
		KineticsModel.__init__(self)
		self.A = A
		self.E0 = E0
		self.n = n
		self.alpha = alpha

	def __str__(self):
		return 'k(T) = %s * T ** %s * math.exp(-(%s + %s * DHrxn) / constants.R / T)\t%s < T < %s' % (self.A, self.n, self.E0, self.alpha, self.Tmin, self.Tmax)
		
	def __repr__(self):
		return '<ArrheniusEPModel A=%.0e E0=%.0f kJ/mol n=%.1f alpha=%.1g>' % (self.A, self.E0/1000.0, self.n, self.alpha)
	
	def getActivationEnergy(self, dHrxn):
		"""
		Return the activation energy in J/mol using the enthalpy of reaction 
		`dHrxn` in J/mol.
		"""
		return self.E0 + self.alpha * dHrxn
	
	def getRateConstant(self, T, dHrxn):
		"""
		Return the rate constant k(T) at temperature `T` by evaluating the
		ArrheniusModel expression for a reaction having an enthalpy of reaction
		`dHrxn` in J/mol.
		"""
		if not self.isTemperatureValid(T):
			raise InvalidKineticsModelError('Invalid temperature for rate coefficient evaluation from ArrheniusEPModel.')
		Ea = cython.declare(cython.double)
		Ea = self.getActivationEnergy(dHrxn)
		if cython.compiled:
			return self.A * pow(T, self.n) * exp(-self.Ea / constants.R / T)
		else:
			return self.A * (T ** self.n) * math.exp(-self.Ea / constants.R / T)

################################################################################

class PDepArrheniusModel(KineticsModel):
	"""
	A kinetic model of a phenomenological rate coefficient k(T, P) using the
	expression

	.. math:: k(T,P) = A(P) T^{n(P)} \\exp \\left[ \\frac{-E_\\mathrm{a}(P)}{RT} \\right]

	where the modified Arrhenius parameters are stored at a variety of pressures
	and interpolated between on a logarithmic scale. The attributes are:

	=============== =============== ============================================
	Attribute       Type            Description
	=============== =============== ============================================
	`pressures`     :class:`list`   The list of pressures in Pa
	`arrhenius`     :class:`list`   The list of :class:`ArrheniusModel` objects at each pressure
	=============== =============== ============================================
	
	"""

	def __init__(self, pressures=None, arrhenius=None):
		KineticsModel.__init__(self)
		self.pressures = pressures or []
		self.arrhenius = arrhenius or []

	def __getAdjacentExpressions(self, P):
		"""
		Returns the pressures and ArrheniusModel expressions for the pressures that
		most closely bound the specified pressure `P` in Pa.
		"""
		cython.declare(Plow=cython.double, Phigh=cython.double)
		cython.declare(arrh=ArrheniusModel)
		cython.declare(i=cython.int, ilow=cython.int, ihigh=cython.int)
		
		if P in self.pressures:
			arrh = self.arrhenius[self.pressures.index(P)]
			return P, P, arrh, arrh
		else:
			ilow = 0; ihigh = -1; Plow = self.pressures[0]; Phigh = 0.0
			for i in range(1, len(self.pressures)):
				if self.pressures[i] <= P:
					ilow = i; Plow = P
				if self.pressures[i] > P and ihigh is None:
					ihigh = i; Phigh = P
			
			return Plow, Phigh, self.arrhenius[ilow], self.arrhenius[ihigh]
	
	def getRateConstant(self, T, P):
		"""
		Return the rate constant k(T, P) at temperature `T` in K and pressure 
		`P` in Pa by evaluating the pressure-dependent Arrhenius expression.
		"""
		if not self.isTemperatureValid(T):
			raise InvalidKineticsModelError('Invalid temperature for rate coefficient evaluation from PDepArrheniusModel.')
		if not self.isPressureValid(P):
			raise InvalidKineticsModelError('Invalid pressure for rate coefficient evaluation from PDepArrheniusModel.')
		
		cython.declare(Plow=cython.double, Phigh=cython.double)
		cython.declare(alow=ArrheniusModel, ahigh=ArrheniusModel)
		cython.declare(klow=cython.double, khigh=cython.double)
		
		Plow, Phigh, alow, ahigh = self.__getAdjacentExpressions(P)
		if Plow == Phigh: return alow.getRateConstant(T)
		klow = alow.getRateConstant(T)
		khigh = ahigh.getRateConstant(T)
		if cython.compiled:
			return pow(log10(P/Plow)/log10(Phigh/Plow)*log10(khigh/klow), 10)
		else:
			return 10**(math.log10(P/Plow)/math.log10(Phigh/Plow)*math.log10(khigh/klow))

################################################################################

class ChebyshevModel(KineticsModel):
	"""
	A kinetic model of a phenomenological rate coefficient k(T, P) using the
	expression
	
	.. math:: \\log k(T,P) = \\sum_{t=1}^{N_T} \\sum_{p=1}^{N_P} \\alpha_{tp} \\phi_t(\\tilde{T}) \\phi_p(\\tilde{P})
	
	where :math:`\\alpha_{tp}` is a constant, :math:`\\phi_n(x)` is the
	Chebyshev polynomial of degree :math:`n` evaluated at :math:`x`, and
	
	.. math:: \\tilde{T} \\equiv \\frac{2T^{-1} - T_\\mathrm{min}^{-1} - T_\\mathrm{max}^{-1}}{T_\\mathrm{max}^{-1} - T_\\mathrm{min}^{-1}}
	
	.. math:: \\tilde{P} \\equiv \\frac{2 \\log P - \\log P_\\mathrm{min} - \\log P_\\mathrm{max}}{\\log P_\\mathrm{max} - \\log P_\\mathrm{min}}
	
	are reduced temperature and reduced pressures designed to map the ranges
	:math:`(T_\\mathrm{min}, T_\\mathrm{max})` and
	:math:`(P_\\mathrm{min}, P_\\mathrm{max})` to :math:`(-1, 1)`.
	The attributes are:
	
	=============== =============== ============================================
	Attribute       Type            Description
	=============== =============== ============================================
	`coeffs`        :class:`list`   Matrix of Chebyshev coefficients
	`degreeT`       :class:`int`    The number of terms in the inverse temperature direction
	`degreeP`       :class:`int`    The number of terms in the log pressure direction
	=============== =============== ============================================
	
	"""

	def __init__(self, Tmin=0.0, Tmax=0.0, Pmin=0.0, Pmax=0.0, coeffs=None):
		KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax)
		self.coeffs = coeffs
		self.degreeT = 0
		self.degreeP = 0

	def __chebyshev(self, n, x):
		if cython.compiled:
			return cos(n * acos(x))
		else:
			return math.cos(n * math.acos(x))

	def __getReducedTemperature(self, T):
		return (2.0/T - 1.0/self.Tmin - 1.0/self.Tmax) / (1.0/self.Tmax - 1.0/self.Tmin)
	
	def __getReducedPressure(self, P):
		if cython.compiled:
			return (2.0*log10(P) - log10(self.Pmin) - log10(self.Pmax)) / (log10(self.Pmax) - log10(self.Pmin))
		else:
			return (2.0*math.log(P) - math.log(self.Pmin) - math.log(self.Pmax)) / (math.log(self.Pmax) - math.log(self.Pmin))
	
	def getRateConstant(self, T, P):
		"""
		Return the rate constant k(T, P) in SI units at temperature `T` in K and
		pressure `P` in Pa by evaluating the Chebyshev expression.
		"""
		
		if not self.isTemperatureValid(T):
			raise InvalidKineticsModelError('Invalid temperature for rate coefficient evaluation from ChebyshevModel.')
		if not self.isPressureValid(P):
			raise InvalidKineticsModelError('Invalid pressure for rate coefficient evaluation from ChebyshevModel.')
		
		cython.declare(Tred=cython.double, Pred=cython.double, k=cython.double)
		cython.declare(t=cython.int, p=cython.int)
		
		Tred = self.__getReducedTemperature(T)
		Pred = self.__getReducedPressure(P)
		
		k = 0.0
		for t in range(self.degreeT):
			for p in range(self.degreeP):
				k += self.coeffs[t,p] * self.__chebyshev(t, Tred) * self.__chebyshev(p, Pred)
		
		if cython.compiled:
			return pow(k, 10.0)
		else:
			return 10.0**k
