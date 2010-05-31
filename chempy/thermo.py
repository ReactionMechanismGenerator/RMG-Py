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
This module contains the thermodynamics models that are available in ChemPy.
All such models derive from the :class:`ThermoModel` base class.
"""

################################################################################

import math
import cython

import constants
from exception import InvalidThermoModelError

################################################################################

class ThermoModel:
	"""
	A base class for thermodynamics models, containing several attributes
	common to all models:
	
	=============== =============== ============================================
	Attribute       Type            Description
	=============== =============== ============================================
	`Tmin`          :class:`float`  The minimum temperature in K at which the model is valid
	`Tmax`          :class:`float`  The maximum temperature in K at which the model is valid
	`comment`       :class:`str`    A string containing information about the model (e.g. its source)
	=============== =============== ============================================
	
	"""
	
	def __init__(self, Tmin=0.0, Tmax=1.0e10, comment=''):
		self.Tmin = Tmin
		self.Tmax = Tmax
		self.comment = comment
	
	def isTemperatureValid(self, T):
		"""
		Return ``True`` if the temperature `T` in K is within the valid
		temperature range of the thermodynamic data, or ``False`` if not.
		"""
		return self.Tmin <= T and T <= self.Tmax

################################################################################

class WilhoitModel(ThermoModel):
	"""
	A thermodynamics model based on the Wilhoit equation for heat capacity,

	.. math::
		C_\\mathrm{p}(T) = C_\\mathrm{p}(0) + \\left[ C_\\mathrm{p}(\\infty) -
		C_\\mathrm{p}(0) \\right] y^2 \\left[ 1 + (y - 1) \\sum_{i=0}^3 a_i y^i \\right]

	where :math:`y \\equiv \\frac{T}{T + B}` is a scaled temperature that ranges
	from zero to one. (The characteristic temperature :math:`B` is chosen by
	default to be 500 K.) This formulation has the advantage of correctly
	reproducting the heat capacity behavior as :math:`T \\rightarrow 0` and
	:math:`T \\rightarrow \\infty`. The low-temperature limit 
	:math:`C_\\mathrm{p}(0)` is taken to be :math:`3.5R` for linear molecules
	and :math:`4R` for nonlinear molecules. The high-temperature limit 
	:math:`C_\\mathrm{p}(\\infty)` is taken to be 
	:math:`\\left[ 3 N_\\mathrm{atoms} - 1.5 \\right] R` for linear molecules and
	:math:`\\left[ 3 N_\\mathrm{atoms} - (2 + 0.5 N_\\mathrm{rotors}) \\right] R`
	for nonlinear molecules, for a molecule composed of :math:`N_\\mathrm{atoms}`
	atoms and :math:`N_\\mathrm{rotors}` internal rotors.
	
	The Wilhoit parameters are stored in the attributes `cp0`, `cpInf`, `a0`,
	`a1`, `a2`, `a3`, and `B`. There are also integration constants `H0` and
	`S0` that are needed to evaluate the enthalpy and entropy, respectively.
	"""

	
	def __init__(self, cp0, cpInf, a0, a1, a2, a3, H0, S0, comment='', B=500.0):
		"""Initialise the Wilhoit polynomial. (Tmin,Tmax) is set to (0,9999.9)"""
		ThermoModel.__init__(self, comment=comment)
		self.cp0 = cp0
		self.cpInf = cpInf
		self.B = B
		self.a0 = a0
		self.a1 = a1
		self.a2 = a2
		self.a3 = a3
		self.H0 = H0
		self.S0 = S0
	
	def getHeatCapacity(self, T):
		"""
		Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
		"""
		cython.declare(y=cython.double)
		y = T/(T+self.B)
		return self.cp0+(self.cpInf-self.cp0)*y*y*( 1 +
			(y-1)*(self.a0 + y*(self.a1 + y*(self.a2 + y*self.a3))) )
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K. The formula used
		is

		.. math::
			H(T) = H_0 +
			C_\\mathrm{p}(0) T + \\left[ C_\\mathrm{p}(\\infty) - C_\\mathrm{p}(0) \\right] T
			\\left\\{ \\left[ 2 + \\sum_{i=0}^3 a_i \\right]
			\\left[ \\frac{1}{2}y - 1 + \\left( \\frac{1}{y} - 1 \\right) \\ln \\frac{T}{y} \\right]
			+ y^2 \\sum_{i=0}^3 \\frac{y^i}{(i+2)(i+3)} \\sum_{j=0}^3 f_{ij} a_j
			\\right\\}

		where :math:`f_{ij} = 3 + j` if :math:`i = j`, :math:`f_{ij} = 1` if
		:math:`i > j`, and :math:`f_{ij} = 0` if :math:`i < j`.
		"""
		return self.H0 + self.__integral_T0(T)
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K. The formula used
		is

		.. math::
			S(T) = S_0 +
			C_\\mathrm{p}(0) \\ln T - \\left[ C_\\mathrm{p}(\\infty) - C_\\mathrm{p}(0) \\right]
			\\left[ \\ln y + \\left( 1 + y \\sum_{i=0}^3 \\frac{a_i y^i}{2+i} \\right) y
			\\right]

		"""
		return self.S0 + self.__integral_TM1(T)
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		return self.getEnthalpy(T) - T * self.getEntropy(T)
	
	#a faster version of the integral based on H from Yelvington's thesis; it differs from the original (see above) by a constant (dependent on parameters but independent of t)
	def __integral_T0(self, t):
		#output: the quantity Integrate[Cp(Wilhoit)/R, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(y=cython.double, y2=cython.double, logBplust=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		y = t/(t+B)
		y2 = y*y
		if cython.compiled:
			logBplust = log(B + t)
		else:
			logBplust = math.log(B + t)
		result = cp0*t - (cpInf-cp0)*t*(y2*((3*a0 + a1 + a2 + a3)/6. + (4*a1 + a2 + a3)*y/12. + (5*a2 + a3)*y2/20. + a3*y2*y/5.) + (2 + a0 + a1 + a2 + a3)*( y/2. - 1 + (1/y-1)*logBplust))
		return result
	
	#a faster version of the integral based on S from Yelvington's thesis; it differs from the original by a constant (dependent on parameters but independent of t)
	def __integral_TM1(self, t):
		#output: the quantity Integrate[Cp(Wilhoit)/R*t^-1, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(y=cython.double, logt=cython.double, logy=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		y = t/(t+B)
		if cython.compiled:
			logy = log(y); logt = log(t)
		else:
			logy = math.log(y); logt = math.log(t)
		result = cpInf*logt-(cpInf-cp0)*(logy+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5)))))
		return result

################################################################################

class NASAPolynomial(ThermoModel):
	"""
	A single NASA polynomial for thermodynamic data. The `coeffs` attribute
	stores the seven polynomial coefficients
	:math:`\\mathbf{a} = \\left[a_1\\ a_2\\ a_3\\ a_4\\ a_5\\ a_6\\ a_7 \\right]`
	from which the relevant thermodynamic parameters are evaluated via the
	expressions
	
	.. math:: \\frac{C_\\mathrm{p}(T)}{R} = a_1 + a_2 T + a_3 T^2 + a_4 T^3 + a_5 T^4
	
	.. math:: \\frac{H(T)}{RT} = a_1 + \\frac{1}{2} a_2 T + \\frac{1}{3} a_3 T^2 + \\frac{1}{4} a_4 T^3 + \\frac{1}{5} a_5 T^4 + \\frac{a_6}{T}
	
	.. math:: \\frac{S(T)}{R} = a_1 \\ln T + a_2 T + \\frac{1}{2} a_3 T^2 + \\frac{1}{3} a_4 T^3 + \\frac{1}{4} a_5 T^4 + a_7
	
	The above was adapted from `this page <http://www.me.berkeley.edu/gri-mech/data/nasa_plnm.html>`_.
	"""
	
	def __init__(self, Tmin=0.0, Tmax=0.0, coeffs=None, comment=''):
		ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
		coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
		
	def getHeatCapacity(self, T):
		"""
		Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise InvalidThermoModelError('Invalid temperature for heat capacity estimation from NASA polynomial.')
		# Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
		return (self.c0 + T*(self.c1 + T*(self.c2 + T*(self.c3 + self.c4*T)))) * constants.R
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise InvalidThermoModelError('Invalid temperature for enthalpy estimation from NASA polynomial.')
		T2 = cython.declare(cython.double)
		T4 = cython.declare(cython.double)
		T2 = T*T
		T4 = T2*T2
		# H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
		return ( self.c0*T + self.c1*T2/2 + self.c2*T2*T/3 + self.c3*T4/4 +
						  self.c4*T4*T/5 + self.c5 ) * constants.R
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise InvalidThermoModelError('Invalid temperature for entropy estimation from NASA polynomial.')
		# S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
		T2 = cython.declare(cython.double)
		T4 = cython.declare(cython.double)
		T2 = T*T
		T4 = T2*T2
		if cython.compiled: # we've imported log from math.h in the pxd file
			return ( self.c0*log(T) + self.c1*T + self.c2*T2/2 +
					self.c3*T2*T/3 + self.c4*T4/4 + self.c6 ) * constants.R
		else:
			return ( self.c0*math.log(T) + self.c1*T + self.c2*T2/2 +
					self.c3*T2*T/3 + self.c4*T4/4 + self.c6 ) * constants.R
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise InvalidThermoModelError('Invalid temperature for free energy estimation from NASA polynomial.')
		return self.getEnthalpy(T) - T * self.getEntropy(T)

################################################################################

class NASAModel(ThermoModel):
	"""
	A set of thermodynamic parameters given by NASA polynomials. This class
	stores a list of :class:`NASAPolynomial` objects in the `polynomials`
	attribute. When evaluating a thermodynamic quantity, a polynomial that
	contains the desired temperature within its valid range will be used.
	"""
	
	def __init__(self, polynomials=None, Tmin=0.0, Tmax=0.0, comment=''):
		ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
		self.polynomials = polynomials or []
	
	def getHeatCapacity(self, T):
		"""
		Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
		"""
		poly = cython.declare(NASAPolynomial)
		poly = self.__selectPolynomialForTemperature(T)
		return poly.getHeatCapacity(T)
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		poly = cython.declare(NASAPolynomial)
		poly = self.__selectPolynomialForTemperature(T)
		return poly.getEnthalpy(T)
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		poly = cython.declare(NASAPolynomial)
		poly = self.__selectPolynomialForTemperature(T)
		return poly.getEntropy(T)
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		poly = cython.declare(NASAPolynomial)
		poly = self.__selectPolynomialForTemperature(T)
		return poly.getFreeEnergy(T)
	
	def __selectPolynomialForTemperature(self, T):
		poly = cython.declare(NASAPolynomial)
		for poly in self.polynomials:
			if poly.isTemperatureValid(T): return poly
		else:
			raise InvalidThermoModelError("No polynomial found for T=%s" % T)

################################################################################