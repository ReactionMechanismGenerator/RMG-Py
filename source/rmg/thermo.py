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
A module for working with the thermodynamics of chemical species. This module
seeks to provide functionality for answering the question, "Given a species,
what are its thermodynamics?"

This module can be compiled using Cython to a shared library, which provides a
significant speed boost to running in pure Python mode.
"""

import quantities as pq
import constants
import data
import math
import scipy
from scipy import linalg
from scipy import optimize
from scipy import integrate
import cython
import log as logging

import ctml_writer

################################################################################

# this will be replaced with something in io.py
forbiddenStructures = None

class ThermoData:
	"""
	A base class for all forms of thermodynamic data used by RMG. The common
	attributes are:
	
	========= ==============================================================
	Attribute Meaning
	========= ==============================================================
	`Trange`  a list of length 2 containing the min and max temperature in K
	`comment` a string describing the source of the data
	========= ==============================================================
	"""
	
	def __init__(self, Trange=None, comment=''):
		Trange = Trange or [0.0, 0.0]
		self.Trange=Trange
		self.Tmin = Trange[0]
		self.Tmax = Trange[1]
		self.comment = comment
	
	def isTemperatureValid(self, T):
		"""
		Return :data:`True` if the temperature `T` in K is within the valid
		temperature range of the thermodynamic data, or :data:`False` if not.
		
		If  Tmax == Tmin == 0 then returns :data:`True`.
		(This case is for when the range is undefined. Tmax and Tmin must be of
		type :data:`float` because of Cython declaration.)
		"""
		if self.Tmax == 0 and self.Tmin == 0:
			return True
		else:
			return self.Tmin <= T and T <= self.Tmax

	def fromXML(self, document, rootElement):
		"""
		Convert a <thermoData> element from a standard RMG-style XML input file
		into a ThermoData object. `document` is an :class:`io.XML` class
		representing the XML DOM tree, and `rootElement` is the <thermoData>
		element in that tree.
		"""

		# Read comment attribute
		self.comment = document.getAttribute(rootElement, 'comment', required=False, default='')

		# Temperature range not currently read
		self.Trange = [0.0, 0.0]
		self.Tmin = 0.0
		self.Tmax = 0.0

	def getGroundStateEnergy(self):
		"""
		Calculate and return the ground-state energy using this thermo model.
		This is done by calculating the enthalpy at 0 K. The returned energy
		has units of J/mol.
		"""
		# We may not be able to evaluate the Cp(T) model at exactly 0 K
		# Instead, we just evaluate it very close to 0 K
		return self.getEnthalpy(0.001)

################################################################################

class ThermoGAData(ThermoData):
	"""
	A set of thermodynamic parameters as determined from Benson's group
	additivity data. The attributes are:
	
	========= ========================================================
	Attribute Meaning
	========= ========================================================
	`H298`    the standard enthalpy of formation at 298 K in J/mol
	`S298`    the standard entropy of formation at 298 K in J/mol*K
	`Cp`      the standard heat capacity in J/mol*K at 300, 400, 500, 600, 800, 1000, and 1500 K
	========= ========================================================
	"""
	# I think putting a comment with a colon just before the thing is defined
	# makes it show up in the documentation with autodoc. (or is that just in constants.py?)
	#: The list of temperatures at which heat capacity is defined = [300,400,500,600,800,1000,1500]
	CpTlist = list(pq.Quantity([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K').simplified)
	CpTlist = [float(T) for T in CpTlist]
	# refer to it as ThermoGAData.CpTlist, even when calling it from methods within this Class.
	
	def __init__(self, H298=0.0, S298=0.0, Cp=None, comment='', index=''):
		"""Initialize a set of group additivity thermodynamic data."""
		ThermoData.__init__(self, Trange=(298.0, 2500.0), comment=comment)
		self.H298 = H298
		self.S298 = S298
		self.Cp = Cp or [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
		self.index = index
	
	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (ThermoGAData, (self.H298, self.S298, self.Cp, self.comment, self.index))

	# how do we cythonize 'special' methods?
	# it's so confusing. This breaks pure python mode:
	#@cython.locals(self=ThermoGAData, other=ThermoGAData, new=ThermoGAData, i=cython.int)
	def __add__(self, other):
		"""
		Add two sets of thermodynamic data together. All parameters are
		considered additive. Returns a new :class:`ThermoGAData` object that is
		the sum of the two sets of thermodynamic data.
		"""
		cython.declare(i=int, new=ThermoGAData)
		new = ThermoGAData()
		new.H298 = self.H298 + other.H298
		new.S298 = self.S298 + other.S298
		new.Cp = [self.Cp[i] + other.Cp[i] for i in range(len(self.Cp))]
		if self.comment == '': new.comment = other.comment
		elif other.comment == '': new.comment = self.comment
		else: new.comment = self.comment + '+ ' + other.comment
		new.index = self.index + '+' + other.index
		return new
	
	def __repr__(self):
		string = 'ThermoGAData(H298=%s, S298=%s, Cp=%s, index="%s")'%(self.H298,self.S298,self.Cp,self.index)
		return string
	
	def __str__(self):
		"""
		Return a string summarizing the thermodynamic data.
		"""
		string = ''
		string += 'Enthalpy of formation: %s J/mol\n' % (self.H298)
		string += 'Entropy of formation: %s J/mol*K\n' % (self.S298)
		string += 'Heat capacity (J/mol*K): '
		for T, Cp in zip(ThermoGAData.CpTlist, self.Cp):
			string += '%s(%sK) ' % (Cp,T)
		string += '\n'
		string += 'Index: %s\tComment: %s' % (self.index, self.comment)
		return string
	
	def equals(self, other):
		"""
		Equality comparison.
		"""
		if self.Trange != other.Trange: return False
		if self.comment != other.comment: return False
		if self.H298 != other.H298: return False
		if self.S298 != other.S298: return False
		if len(self.Cp) != len(other.Cp): return False
		for i in range(len(self.Cp)):
			if self.Cp[i] != other.Cp[i]: return False
		return True

	def getHeatCapacity(self, T):
		"""
		Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from group additivity.')
		if T < 300.0:
			return self.Cp[0]
		elif T > ThermoGAData.CpTlist[-1]:
			return self.Cp[-1]
		else:
			cython.declare(Tmin=cython.double, Tmax=cython.double, Cpmin=cython.double, Cpmax=cython.double)
			for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
					ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
				if Tmin <= T and T <= Tmax:
					return (Cpmax - Cpmin) * ((T - Tmin) / (Tmax - Tmin)) + Cpmin
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""	
		
		cython.declare(H=cython.double, slope=cython.double, intercept=cython.double,
		     Tmin=cython.double, Tmax=cython.double, Cpmin=cython.double, Cpmax=cython.double)
		
		H = self.H298
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for enthalpy estimation from group additivity.')
		for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
				ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
			if T > Tmin:
				slope = (Cpmax - Cpmin) / (Tmax - Tmin)
				intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
				if T < Tmax:	H += 0.5 * slope * (T*T - Tmin*Tmin) + intercept * (T - Tmin)
				else:			H += 0.5 * slope * (Tmax*Tmax - Tmin*Tmin) + intercept * (Tmax - Tmin)
		if T > ThermoGAData.CpTlist[-1]:
			H += self.Cp[-1] * (T - ThermoGAData.CpTlist[-1])
		return H
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		
		S = self.S298
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for entropy estimation from group additivity.')
		for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
				ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
			if T > Tmin:
				slope = (Cpmax - Cpmin) / (Tmax - Tmin)
				intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
				if T < Tmax:	S += slope * (T - Tmin) + intercept * math.log(T/Tmin)
				else:			S += slope * (Tmax - Tmin) + intercept * math.log(Tmax/Tmin)
		if T > ThermoGAData.CpTlist[-1]:
			S += self.Cp[-1] * math.log(T / ThermoGAData.CpTlist[-1])
		return S
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
			
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from group additivity.')
		G = self.getEnthalpy(T) - T * self.getEntropy(T)
		return G
	
	def fromDatabase(self, data, comment):
		"""
		Process a list of numbers `data` and associated description `comment`
		generated while reading from a thermodynamic database.
		"""
		
		if len(data) != 12:
			raise Exception('Invalid list of thermo data; should be a list of numbers of length 12.')
		
		H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, \
			dH, dS, dCp = data
		
		self.H298 = float(pq.Quantity(H298, 'kcal/mol').simplified)
		self.S298 = float(pq.Quantity(S298, 'cal/(mol*K)').simplified)
		self.Cp = list(pq.Quantity([Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500], 'cal/(mol*K)').simplified)
		for i in range(len(self.Cp)): self.Cp[i] = float(self.Cp[i])
		self.comment = comment
	
	def fromXML(self, document, rootElement):
		"""
		Convert a <thermoData> element from a standard RMG-style XML input file
		into a ThermoGAData object. `document` is an :class:`io.XML` class
		representing the XML DOM tree, and `rootElement` is the <thermoData>
		element in that tree.
		"""

		# 'comment' attribute parsed by base class
		ThermoData.fromXML(self, document, rootElement)

		# Read <enthalpyOfFormation> element
		H298 = document.getChildQuantity(rootElement, 'enthalpyOfFormation', required=True)
		self.H298 = float(H298.simplified)

		# Read <entropyOfFormation> element
		S298 = document.getChildQuantity(rootElement, 'entropyOfFormation', required=True)
		self.S298 = float(S298.simplified)

		# Read <heatCapacities> element
		Cp = document.getChildQuantity(rootElement, 'heatCapacities', required=True)
		self.Cp = [float(C.simplified) for C in Cp]

	def toXML(self, document, rootElement):
		"""
		Create a <thermoData> element as a child of `rootElement` in the XML DOM
		tree represented by `document`, an :class:`io.XML` class. The format
		matches the format of the :meth:`ThermoGAData.fromXML()` function.
		"""
		
		# Create <thermoData> element
		thermoDataElement = document.createElement('thermoData', rootElement)
		document.createAttribute('format', thermoDataElement, 'group additivity')

		document.createQuantity('enthalpyOfFormation', thermoDataElement, self.H298/1000.0, 'kJ/mol')
		document.createQuantity('entropyOfFormation', thermoDataElement, self.S298, 'J/(mol*K)')
		document.createQuantity('heatCapacities', thermoDataElement, self.Cp, 'J/(mol*K)')

	
################################################################################

class ThermoNASAPolynomial(ThermoData):
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
	
	def __init__(self, T_range=None, coeffs=None, comment=''):
		ThermoData.__init__(self, Trange=T_range, comment=comment)
		coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		
		self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
		
	def __repr__(self):
		return "ThermoNASAPolynomial(%r,%r,'%s')"%(self.Trange,(self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6),self.comment)
	def __reduce__(self):
		return (ThermoNASAPolynomial,(self.Trange,(self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6),self.comment))

	def equals(self, other):
		"""
		Equality comparison.
		"""
		if self.Trange != other.Trange: return False
		if self.comment != other.comment: return False
		if self.c0 != other.c0: return False
		if self.c1 != other.c1: return False
		if self.c2 != other.c2: return False
		if self.c3 != other.c3: return False
		if self.c4 != other.c4: return False
		if self.c5 != other.c5: return False
		if self.c6 != other.c6: return False
		return True

	def getHeatCapacity(self, T):
		"""
		Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from NASA polynomial.')
		T2 = cython.declare(cython.double)
		HeatCapacityOverR = cython.declare(cython.double)
		
		T2 = T*T
		# Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
		# HeatCapacityOverR = self.c0 + self.c1*T + self.c2*T*T + self.c3*T*T*T + self.c4*T*T*T*T
		HeatCapacityOverR = self.c0 + T*(self.c1 + T*(self.c2 + T*(self.c3 + self.c4*T)))
		return HeatCapacityOverR * constants.R
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for enthalpy estimation from NASA polynomial.')
		T2 = cython.declare(cython.double)
		T4 = cython.declare(cython.double)
		EnthalpyOverR = cython.declare(cython.double)
		
		T2 = T*T
		T4 = T2*T2
		# H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
		EnthalpyOverR = ( self.c0*T + self.c1*T2/2 + self.c2*T2*T/3 + self.c3*T4/4 +
						  self.c4*T4*T/5 + self.c5 )
		return EnthalpyOverR * constants.R
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for entropy estimation from NASA polynomial.')
		# S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
		T2 = cython.declare(cython.double)
		T4 = cython.declare(cython.double)
		EntropyOverR = cython.declare(cython.double)
		
		T2 = T*T
		T4 = T2*T2
		if cython.compiled: # we've imported log from math.h in the pxd file
			EntropyOverR = ( self.c0*log(T) + self.c1*T + self.c2*T2/2 +
					self.c3*T2*T/3 + self.c4*T4/4 + self.c6 )
		else:
			EntropyOverR = ( self.c0*math.log(T) + self.c1*T + self.c2*T2/2 +
					self.c3*T2*T/3 + self.c4*T4/4 + self.c6 )
		return EntropyOverR * constants.R
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from NASA polynomial.')
		return self.getEnthalpy(T) - T * self.getEntropy(T)
	
	def toXML(self, dom, root):
		"""Append to xml `dom` at node `root`"""
### prime-like:
#   <polynomial>
#      <validRange>
#         <bound kind="lower" property="temperature" units="K">300.000</bound>
#         <bound kind="upper" property="temperature" units="K">1000</bound>
#      </validRange>
#      <coefficient id="1" label="a1">4.0733</coefficient>
#      <coefficient id="2" label="a2">0.011308</coefficient>
#      <coefficient id="3" label="a3">-1.6565e-005</coefficient>
#      <coefficient id="4" label="a4">1.1784e-008</coefficient>
#      <coefficient id="5" label="a5">-3.3006e-012</coefficient>
#      <coefficient id="6" label="a6">-19054.23</coefficient>
#      <coefficient id="7" label="a7">4.4083</coefficient>
#   </polynomial>
		pass
		
	def toCantera(self):
		"""Return a Cantera ctml_writer instance"""
		return ctml_writer.NASA([self.Tmin,self.Tmax], [self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6])
		
	
	def integral2_T0(self, t):
		#input: NASA parameters for Cp/R, c1, c2, c3, c4, c5 (either low or high temp parameters), temperature t (in kiloKelvin; an endpoint of the low or high temp range
		#output: the quantity Integrate[(Cp(NASA)/R)^2, t'] evaluated at t'=t 
		#can speed further by precomputing and storing e.g. thigh^2, tlow^2, etc.
		cython.declare(c1=cython.double, c2=cython.double,c3=cython.double,c4=cython.double,c5=cython.double)
		cython.declare(result=cython.double)
		c1, c2, c3, c4, c5 = self.c0, self.c1, self.c2, self.c3, self.c4
		result = c1*c1*t + c1*c2*t*t + (2*c1*c3+c2*c2)/3*t*t*t + (c1*c4+c2*c3)/2*t*t*t*t + (2*c1*c5 + 2*c2*c4 + c3*c3)/5*t*t*t*t*t + (c2*c5 + c3*c4)/3*t*t*t*t*t*t + (2*c3*c5 + c4*c4)/7*t*t*t*t*t*t*t + c4*c5/4*t*t*t*t*t*t*t*t + c5*c5/9*t*t*t*t*t*t*t*t*t
		return result

	def integral2_TM1(self, t):
		#input: NASA parameters for Cp/R, c1, c2, c3, c4, c5 (either low or high temp parameters), temperature t (in kiloKelvin; an endpoint of the low or high temp range
		#output: the quantity Integrate[(Cp(NASA)/R)^2*t^-1, t'] evaluated at t'=t 
		#can speed further by precomputing and storing e.g. thigh^2, tlow^2, etc.
		cython.declare(c1=cython.double, c2=cython.double,c3=cython.double,c4=cython.double,c5=cython.double)
		cython.declare(result=cython.double)
		c1, c2, c3, c4, c5 = self.c0, self.c1, self.c2, self.c3, self.c4
		if cython.compiled: # we've imported log from math.h in the pxd file
			result = c1*c1*log(t)
		else:
			result = c1*c1*math.log(t)
		result = result + 2*c1*c2*t + (2*c1*c3+c2*c2)/2*t*t + 2*(c1*c4+c2*c3)/3*t*t*t + (2*c1*c5 + 2*c2*c4 + c3*c3)/4*t*t*t*t + 2*(c2*c5 + c3*c4)/5*t*t*t*t*t + (2*c3*c5 + c4*c4)/6*t*t*t*t*t*t + 2*c4*c5/7*t*t*t*t*t*t*t + c5*c5/8*t*t*t*t*t*t*t*t
		return result

################################################################################

class ThermoNASAData(ThermoData):
	"""
	A set of thermodynamic parameters given by NASA polynomials. This class
	stores a list of :class:`ThermoNASAPolynomial` objects in the `polynomials`
	attribute. When evaluating a thermodynamic quantity, a polynomial that
	contains the desired temperature within its valid range will be used.
	"""
	
	def __init__(self, polynomials=None, comment='', Trange=None):
		ThermoData.__init__(self, Trange=Trange, comment=comment)
		self.polynomials = polynomials or []
	def __reduce__(self):
		return (ThermoNASAData,(self.polynomials,self.comment,self.Trange))
		
	def equals(self, other):
		"""
		Equality comparison.
		"""
		if self.Trange != other.Trange: return False
		if self.comment != other.comment: return False
		for i in range(len(self.polynomials)):
			if not self.polynomials[i].equals(other.polynomials[i]): return False
		return True

	def addPolynomial(self, polynomial):
		if not isinstance(polynomial,ThermoNASAPolynomial):
			raise TypeError("Polynomial attribute should be instance of ThermoNASAPolynomial")
		self.polynomials.append(polynomial)
	
	def selectPolynomialForTemperature(self, T):
		for poly in self.polynomials:
			if poly.isTemperatureValid(T): break
		else:
			raise data.TemperatureOutOfRangeException("No polynomial found for T=%s" % T)
		return poly
	
	def __repr__(self):
		return "ThermoNASAData(%s, '%s')"%(repr(self.polynomials),self.comment)
	
	def toXML(self, dom, root):
### prime-like: 
# <thermodynamicPolynomials type="nasa7">
#   <referenceState>
#      <Tref units="K">298.15</Tref>
#      <Pref units="Pa">100000</Pref>
#   </referenceState>
#   <dfH units="J/mol">-145177.5731</dfH>
#   <polynomial> FROM NASA_polynomials </polynomial>
# </thermodynamicPolynomials>
		pass
	
	def toCantera(self):
		"""Return a Cantera ctml_writer instance"""
		return tuple([poly.toCantera() for poly in self.polynomials])
	
	def getHeatCapacity(self, T):
		"""
		Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
		return poly.getHeatCapacity(T)
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
		return poly.getEnthalpy(T)
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
		return poly.getEntropy(T)
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		poly = self.selectPolynomialForTemperature(T)
		return poly.getFreeEnergy(T)
	
	def toCHEMKIN(self):
		"""Return the latter ~half of the first line of a CHEMKIN thermo line 
		and the other three full lines for a case with two polynomials; 
		note, this is a quick and dirty implementation; there is an extra 
		zero in each exponent that must be manually deleted"""
		#potential fix to avoid need to manually delete zeroes:
		#replace "+0" and "-0" substrings with "+" and "-", respectively
		
		# nb. for me it works as written and produces numbers like 3.00000000E+00 4.00000000E-11 
		# what is the problem?  --Richard
		
		low = self.polynomials[0]
		high = self.polynomials[1]
		line1 = "G  %8.3F  %8.3F  %8.3F    1\n"%(self.Tmin,self.Tmax, low.Tmax)
		line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n"%(high.c0,high.c1,high.c2,high.c3,high.c4)
		line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n"%(high.c5,high.c6,low.c0,low.c1,low.c2)
		line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4"%(low.c3,low.c4,low.c5,low.c6)
		return line1 + line2 + line3 + line4

	def rmsErr(self, thermoGAdata):
		"""
		Calculate the RMS error between the NASA polynomial and training data points

		input: thermoGAdata
		output: value is in non-dimensional units (/R);
		"""
		t=ThermoGAData.CpTlist
		cp=thermoGAdata.Cp
		R = constants.R
		m = len(t)
		assert (len(cp)==m), 'cp and t are different lengths'
		rms = 0.0
		for i in range(m):
			err = cp[i]-self.getHeatCapacity(t[i])
			rms += err*err
		rms = rms/m
		rms = math.sqrt(rms)/R

		return rms

################################################################################

#: The default temperature in K
ThermoWilhoitDataB = 500.0

class ThermoWilhoitData(ThermoData):
	"""
	A set of thermodynamic parameters given by Wilhoit polynomials, which have
	the form

	.. math::
		C_\\mathrm{p}(T) = C_\\mathrm{p}(0) + \\left[ C_\\mathrm{p}(\\infty) -
		C_\\mathrm{p}(0) \\right] y^2 \\left[ 1 + (y - 1) \\sum_{i=0}^3 a_i y^i \\right]

	where :math:`y \\equiv \\frac{T}{T + B}` is a scaled temperature that ranges
	from zero to one. The characteristic temperature :math:`B` is chosen by
	default to be 500 K. This formulation has the advantage of correctly
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

	
	def __init__(self, cp0, cpInf, a0, a1, a2, a3, H0, S0, comment='', B=ThermoWilhoitDataB):
		"""Initialise the Wilhoit polynomial. Trange is set to (0,9999.9)"""
		Trange = (0,9999.9) # Wilhoit valid over all temperatures
		ThermoData.__init__(self, Trange=Trange, comment=comment)
		self.cp0 = cp0
		self.cpInf = cpInf
		self.B = B
		self.a0 = a0
		self.a1 = a1
		self.a2 = a2
		self.a3 = a3
		self.H0 = H0
		self.S0 = S0
	
	def __repr__(self):
		return "ThermoWilhoitData(%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,'%s',B=%.4g)"%(self.cp0, self.cpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.comment, self.B)
	
	def __reduce__(self):
		return (ThermoWilhoitData,(self.cp0, self.cpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.comment, self.B))

	def equals(self, other):
		"""
		Equality comparison.
		"""
		if self.Trange != other.Trange: return False
		if self.comment != other.comment: return False
		if self.cp0 != other.cp0: return False
		if self.cpInf != other.cpInf: return False
		if self.B != other.B: return False
		if self.a0 != other.a0: return False
		if self.a1 != other.a1: return False
		if self.a2 != other.a2: return False
		if self.a3 != other.a3: return False
		if self.H0 != other.H0: return False
		if self.S0 != other.S0: return False
		return True

	def toXML(self, dom, root):
		pass
	
	def getHeatCapacity(self, T):
		"""
		Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
		"""
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
		return self.H0 + self.integral_T0(T)
	
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
		return self.S0 + self.integral_TM1(T)
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		return self.getEnthalpy(T) - T * self.getEntropy(T)
	
	def rmsErrWilhoit(self,t,cp):
		#calculate the RMS error between the Wilhoit form and training data points; result will have same units as cp inputs; cp, cp0, and cpInf should agree in units (e.g. Cp/R); units of B and t should be consistent, based, for example on kK or K 
		m = len(t)
		rms = 0.0
		for i in range(m):
			err = cp[i]-self.getHeatCapacity(t[i])
			rms += err*err
		rms = rms/m
		rms = math.sqrt(rms)
		
		return rms
	
	#a faster version of the integral based on H from Yelvington's thesis; it differs from the original (see above) by a constant (dependent on parameters but independent of t)
	def integral_T0(self, t):
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
	def integral_TM1(self, t):
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
	
	def integral_T1(self, t):
		#output: the quantity Integrate[Cp(Wilhoit)/R*t, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(logBplust=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		if cython.compiled:
			logBplust = log(B + t)
		else:
			logBplust = math.log(B + t)
		result = ( (2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t + (cpInf*t**2)/2. + (a3*B**7*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 6*a3)*B**6*(cp0 - cpInf))/(4.*(B + t)**4) -
			((a1 + 5*(a2 + 3*a3))*B**5*(cp0 - cpInf))/(3.*(B + t)**3) + ((a0 + 4*a1 + 10*(a2 + 2*a3))*B**4*(cp0 - cpInf))/(2.*(B + t)**2) -
			((1 + 3*a0 + 6*a1 + 10*a2 + 15*a3)*B**3*(cp0 - cpInf))/(B + t) - (3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*logBplust)
		return result
	
	def integral_T2(self, t):
		#output: the quantity Integrate[Cp(Wilhoit)/R*t^2, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(logBplust=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		if cython.compiled:
			logBplust = log(B + t)
		else:
			logBplust = math.log(B + t)
		result = ( -((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*t) + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**2)/2. + (cpInf*t**3)/3. + (a3*B**8*(cp0 - cpInf))/(5.*(B + t)**5) -
			((a2 + 7*a3)*B**7*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 6*a2 + 21*a3)*B**6*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 5*(a1 + 3*a2 + 7*a3))*B**5*(cp0 - cpInf))/(2.*(B + t)**2) +
			((1 + 4*a0 + 10*a1 + 20*a2 + 35*a3)*B**4*(cp0 - cpInf))/(B + t) + (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*logBplust)
		return result

	def integral_T3(self, t):
		#output: the quantity Integrate[Cp(Wilhoit)/R*t^3, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(logBplust=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		if cython.compiled:
			logBplust = log(B + t)
		else:
			logBplust = math.log(B + t)
		result = ( (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t + ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**2)/2. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**3)/3. +
			(cpInf*t**4)/4. + (a3*B**9*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 8*a3)*B**8*(cp0 - cpInf))/(4.*(B + t)**4) - ((a1 + 7*(a2 + 4*a3))*B**7*(cp0 - cpInf))/(3.*(B + t)**3) +
			((a0 + 6*a1 + 21*a2 + 56*a3)*B**6*(cp0 - cpInf))/(2.*(B + t)**2) - ((1 + 5*a0 + 15*a1 + 35*a2 + 70*a3)*B**5*(cp0 - cpInf))/(B + t) -
			(5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*logBplust)
		return result

	def integral_T4(self, t):
		#output: the quantity Integrate[Cp(Wilhoit)/R*t^4, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(logBplust=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		if cython.compiled:
			logBplust = log(B + t)
		else:
			logBplust = math.log(B + t)
		result = ( -((5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*t) + ((4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t**2)/2. +
			((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**3)/3. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**4)/4. + (cpInf*t**5)/5. + (a3*B**10*(cp0 - cpInf))/(5.*(B + t)**5) -
			((a2 + 9*a3)*B**9*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 8*a2 + 36*a3)*B**8*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 7*(a1 + 4*(a2 + 3*a3)))*B**7*(cp0 - cpInf))/(2.*(B + t)**2) +
			((1 + 6*a0 + 21*a1 + 56*a2 + 126*a3)*B**6*(cp0 - cpInf))/(B + t) + (6 + 15*a0 + 35*a1 + 70*a2 + 126*a3)*B**5*(cp0 - cpInf)*logBplust)
		return result

	def integral2_T0(self, t):
		#output: the quantity Integrate[(Cp(Wilhoit)/R)^2, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(logBplust=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		if cython.compiled:
			logBplust = log(B + t)
		else:
			logBplust = math.log(B + t)
		result = (cpInf**2*t - (a3**2*B**12*(cp0 - cpInf)**2)/(11.*(B + t)**11) + (a3*(a2 + 5*a3)*B**11*(cp0 - cpInf)**2)/(5.*(B + t)**10) -
			((a2**2 + 18*a2*a3 + a3*(2*a1 + 45*a3))*B**10*(cp0 - cpInf)**2)/(9.*(B + t)**9) + ((4*a2**2 + 36*a2*a3 + a1*(a2 + 8*a3) + a3*(a0 + 60*a3))*B**9*(cp0 - cpInf)**2)/(4.*(B + t)**8) -
			((a1**2 + 14*a1*(a2 + 4*a3) + 2*(14*a2**2 + a3 + 84*a2*a3 + 105*a3**2 + a0*(a2 + 7*a3)))*B**8*(cp0 - cpInf)**2)/(7.*(B + t)**7) +
			((3*a1**2 + a2 + 28*a2**2 + 7*a3 + 126*a2*a3 + 126*a3**2 + 7*a1*(3*a2 + 8*a3) + a0*(a1 + 6*a2 + 21*a3))*B**7*(cp0 - cpInf)**2)/(3.*(B + t)**6) -
			(B**6*(cp0 - cpInf)*(a0**2*(cp0 - cpInf) + 15*a1**2*(cp0 - cpInf) + 10*a0*(a1 + 3*a2 + 7*a3)*(cp0 - cpInf) + 2*a1*(1 + 35*a2 + 70*a3)*(cp0 - cpInf) +
			 2*(35*a2**2*(cp0 - cpInf) + 6*a2*(1 + 21*a3)*(cp0 - cpInf) + a3*(5*(4 + 21*a3)*cp0 - 21*(cpInf + 5*a3*cpInf)))))/(5.*(B + t)**5) +
			(B**5*(cp0 - cpInf)*(14*a2*cp0 + 28*a2**2*cp0 + 30*a3*cp0 + 84*a2*a3*cp0 + 60*a3**2*cp0 + 2*a0**2*(cp0 - cpInf) + 10*a1**2*(cp0 - cpInf) +
			 a0*(1 + 10*a1 + 20*a2 + 35*a3)*(cp0 - cpInf) + a1*(5 + 35*a2 + 56*a3)*(cp0 - cpInf) - 15*a2*cpInf - 28*a2**2*cpInf - 35*a3*cpInf - 84*a2*a3*cpInf - 60*a3**2*cpInf))/
			 (2.*(B + t)**4) - (B**4*(cp0 - cpInf)*((1 + 6*a0**2 + 15*a1**2 + 32*a2 + 28*a2**2 + 50*a3 + 72*a2*a3 + 45*a3**2 + 2*a1*(9 + 21*a2 + 28*a3) + a0*(8 + 20*a1 + 30*a2 + 42*a3))*cp0 -
			 (1 + 6*a0**2 + 15*a1**2 + 40*a2 + 28*a2**2 + 70*a3 + 72*a2*a3 + 45*a3**2 + a0*(8 + 20*a1 + 30*a2 + 42*a3) + a1*(20 + 42*a2 + 56*a3))*cpInf))/(3.*(B + t)**3) +
			(B**3*(cp0 - cpInf)*((2 + 2*a0**2 + 3*a1**2 + 9*a2 + 4*a2**2 + 11*a3 + 9*a2*a3 + 5*a3**2 + a0*(5 + 5*a1 + 6*a2 + 7*a3) + a1*(7 + 7*a2 + 8*a3))*cp0 -
			 (2 + 2*a0**2 + 3*a1**2 + 15*a2 + 4*a2**2 + 21*a3 + 9*a2*a3 + 5*a3**2 + a0*(6 + 5*a1 + 6*a2 + 7*a3) + a1*(10 + 7*a2 + 8*a3))*cpInf))/(B + t)**2 -
			(B**2*((2 + a0 + a1 + a2 + a3)**2*cp0**2 - 2*(5 + a0**2 + a1**2 + 8*a2 + a2**2 + 9*a3 + 2*a2*a3 + a3**2 + 2*a0*(3 + a1 + a2 + a3) + a1*(7 + 2*a2 + 2*a3))*cp0*cpInf +
			 (6 + a0**2 + a1**2 + 12*a2 + a2**2 + 14*a3 + 2*a2*a3 + a3**2 + 2*a1*(5 + a2 + a3) + 2*a0*(4 + a1 + a2 + a3))*cpInf**2))/(B + t) +
			2*(2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*cpInf*logBplust)
		return result

	def integral2_TM1(self, t):
		#output: the quantity Integrate[(Cp(Wilhoit)/R)^2*t^-1, t'] evaluated at t'=t
		cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
		cython.declare(logBplust=cython.double, logt=cython.double, result=cython.double)
		cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
		if cython.compiled:
			logBplust = log(B + t); logt = log(t)
		else:
			logBplust = math.log(B + t); logt = math.log(t)
		result = ( (a3**2*B**11*(cp0 - cpInf)**2)/(11.*(B + t)**11) - (a3*(2*a2 + 9*a3)*B**10*(cp0 - cpInf)**2)/(10.*(B + t)**10) +
			((a2**2 + 16*a2*a3 + 2*a3*(a1 + 18*a3))*B**9*(cp0 - cpInf)**2)/(9.*(B + t)**9) -
			((7*a2**2 + 56*a2*a3 + 2*a1*(a2 + 7*a3) + 2*a3*(a0 + 42*a3))*B**8*(cp0 - cpInf)**2)/(8.*(B + t)**8) +
			((a1**2 + 21*a2**2 + 2*a3 + 112*a2*a3 + 126*a3**2 + 2*a0*(a2 + 6*a3) + 6*a1*(2*a2 + 7*a3))*B**7*(cp0 - cpInf)**2)/(7.*(B + t)**7) -
			((5*a1**2 + 2*a2 + 30*a1*a2 + 35*a2**2 + 12*a3 + 70*a1*a3 + 140*a2*a3 + 126*a3**2 + 2*a0*(a1 + 5*(a2 + 3*a3)))*B**6*(cp0 - cpInf)**2)/(6.*(B + t)**6) +
			(B**5*(cp0 - cpInf)*(10*a2*cp0 + 35*a2**2*cp0 + 28*a3*cp0 + 112*a2*a3*cp0 + 84*a3**2*cp0 + a0**2*(cp0 - cpInf) + 10*a1**2*(cp0 - cpInf) + 2*a1*(1 + 20*a2 + 35*a3)*(cp0 - cpInf) +
			4*a0*(2*a1 + 5*(a2 + 2*a3))*(cp0 - cpInf) - 10*a2*cpInf - 35*a2**2*cpInf - 30*a3*cpInf - 112*a2*a3*cpInf - 84*a3**2*cpInf))/(5.*(B + t)**5) -
			(B**4*(cp0 - cpInf)*(18*a2*cp0 + 21*a2**2*cp0 + 32*a3*cp0 + 56*a2*a3*cp0 + 36*a3**2*cp0 + 3*a0**2*(cp0 - cpInf) + 10*a1**2*(cp0 - cpInf) +
			2*a0*(1 + 6*a1 + 10*a2 + 15*a3)*(cp0 - cpInf) + 2*a1*(4 + 15*a2 + 21*a3)*(cp0 - cpInf) - 20*a2*cpInf - 21*a2**2*cpInf - 40*a3*cpInf - 56*a2*a3*cpInf - 36*a3**2*cpInf))/
			(4.*(B + t)**4) + (B**3*(cp0 - cpInf)*((1 + 3*a0**2 + 5*a1**2 + 14*a2 + 7*a2**2 + 18*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(5 + 6*a2 + 7*a3))*cp0 -
			(1 + 3*a0**2 + 5*a1**2 + 20*a2 + 7*a2**2 + 30*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(6 + 6*a2 + 7*a3))*cpInf))/(3.*(B + t)**3) -
			(B**2*((3 + a0**2 + a1**2 + 4*a2 + a2**2 + 4*a3 + 2*a2*a3 + a3**2 + 2*a1*(2 + a2 + a3) + 2*a0*(2 + a1 + a2 + a3))*cp0**2 -
			2*(3 + a0**2 + a1**2 + 7*a2 + a2**2 + 8*a3 + 2*a2*a3 + a3**2 + 2*a1*(3 + a2 + a3) + a0*(5 + 2*a1 + 2*a2 + 2*a3))*cp0*cpInf +
			(3 + a0**2 + a1**2 + 10*a2 + a2**2 + 12*a3 + 2*a2*a3 + a3**2 + 2*a1*(4 + a2 + a3) + 2*a0*(3 + a1 + a2 + a3))*cpInf**2))/(2.*(B + t)**2) +
			(B*(cp0 - cpInf)*(cp0 - (3 + 2*a0 + 2*a1 + 2*a2 + 2*a3)*cpInf))/(B + t) + cp0**2*logt + (-cp0**2 + cpInf**2)*logBplust)
		return result

################################################################################

def convertGAtoWilhoit(GAthermo, atoms, rotors, linear, fixedB=1, Bmin=300.0, Bmax=6000.0):
	"""Convert a Group Additivity thermo instance into a Wilhoit thermo instance.
	
	Takes a `ThermoGAData` instance of themochemical data, and some extra information 
	about the molecule used to calculate high- and low-temperature limits of Cp.
	These are the number of atoms (integer `atoms`), the number of rotors (integer `rotors`)
	and whether the molecule is linear (boolean `linear`)
	Returns a `ThermoWilhoitData` instance.
	
	cf. Paul Yelvington's thesis, p. 185-186
	"""
	
	# get info from incoming group additivity thermo
	H298 = GAthermo.H298
	S298 = GAthermo.S298
	Cp_list = GAthermo.Cp
	T_list = ThermoGAData.CpTlist  # usually [300, 400, 500, 600, 800, 1000, 1500] but why assume?
	R = constants.R
	B = ThermoWilhoitDataB # Constant (if fixed=1), set once in the class def.
	
	# convert from K to kK
	T_list = [t/1000. for t in T_list] 
	B = B/1000.
	Bmin=Bmin/1000.
	Bmax=Bmax/1000.
	
	Cp_list = [x/R for x in Cp_list] # convert to Cp/R
	
	(cp0, cpInf) = CpLimits(atoms, rotors, linear) # determine the heat capacity limits (non-dimensional)
	
	if (cp0==cpInf):
		a0=0.0
		a1=0.0
		a2=0.0
		a3=0.0
		resid = 0.0
	elif(fixedB == 1):
		(a0, a1, a2, a3, resid) = GA2Wilhoit(B, T_list, Cp_list, cp0, cpInf)		
	else:
		(a0, a1, a2, a3, B, resid) = GA2Wilhoit_BOpt(T_list, Cp_list, cp0, cpInf, Bmin, Bmax)
	m = len(T_list)
	err = math.sqrt(resid/m) # gmagoon 1/19/10: this is a (probably) faster alternative to using rmsErrWilhoit, and it fits better within a scheme where we modify B

	# scale everything back
	T_list = [t*1000. for t in T_list]
	B = B*1000.
	Cp_list = [x*R for x in Cp_list]
	
	# cp0 and cpInf should be in units of J/mol-K
	cp0 = cp0*R
	cpInf = cpInf*R
	
	# output comment
	comment = ''
	
	# first set H0 = S0 = 0, then calculate what they should be
	# by referring to H298, S298
	H0 = 0
	S0 = 0
	# create Wilhoit instance
	WilhoitThermo = ThermoWilhoitData( cp0, cpInf, a0, a1, a2, a3, H0, S0, B=B, comment=comment)
	# calculate correct I, J (integration constants for H, S, respectively)
	H0 = H298 - WilhoitThermo.getEnthalpy(298.15)
	S0 = S298 - WilhoitThermo.getEntropy(298.15)
	# update Wilhoit instance with correct I,J
	WilhoitThermo.H0 = H0
	WilhoitThermo.S0 = S0

	# calculate the correct err for the monoatomic case; there seems to be a bug in linalg.lstsq() where resid is incorrectly returned as [] when A matrix is all zeroes (this is also the reason for the check above for cp0==cpInf, where we set resid = 0)
	if(cp0==cpInf):
		err = WilhoitThermo.rmsErrWilhoit(T_list, Cp_list)/R #rms Error (J/mol-K units until it is divided by R)
	WilhoitThermo.comment = WilhoitThermo.comment + 'Wilhoit function fitted to GA data with Cp0=%2g and Cp_inf=%2g. RMS error = %.3f*R. '%(cp0,cpInf,err) + GAthermo.comment

	#print a warning if the rms fit is worse that 0.25*R
	if (err>0.25):
		logging.warning("Poor GA-to-Wilhoit fit quality: RMS error = %.3f*R" % err)
	
	return WilhoitThermo

def GA2Wilhoit(B, T_list, Cp_list, cp0, cpInf):
	#input: B (in kiloKelvin), GA temperature and Cp_list (non-dimensionalized), Wilhoit parameters, Cp0/R and CpInf/R
	#output: Wilhoit parameters a0-a3, and the sum of squared errors between Wilhoit and GA data

	#create matrices for linear least squares problem
	m = len(T_list)	 # probably m=7
	# A = mx4
	# b = mx1
	# x = 4x1
	A = scipy.zeros([m,4])
	b = scipy.zeros([m])
	for i in range(m):
		y = T_list[i]/(T_list[i]+B)
		A[i,0] = (cpInf-cp0) * y*y*(y-1)
		A[i,1] = A[i,0] * y
		A[i,2] = A[i,1] * y
		A[i,3] = A[i,2] * y
		b[i] = Cp_list[i]-cp0 - y*y*(cpInf-cp0)
		
	#solve least squares problem A*x = b; http://docs.scipy.org/doc/scipy/reference/tutorial/linalg.html#solving-linear-least-squares-problems-and-pseudo-inverses
	x,resid,rank,sigma = linalg.lstsq(A,b, overwrite_a=1, overwrite_b=1)
	a0 = x[0]
	a1 = x[1]
	a2 = x[2]
	a3 = x[3]	
	
	return a0, a1, a2, a3, resid
	
def GA2Wilhoit_BOpt(T_list, Cp_list, cp0, cpInf, Bmin, Bmax):
	#input: GA temperature and Cp_list (scaled/non-dimensionalized), Wilhoit parameters, Cp0/R and CpInf/R, and maximum and minimum bounds for B (in kK)
	#output: Wilhoit parameters, including optimized B value (in kK), and the sum of squared errors between Wilhoit and GA data (dimensionless)
	B = optimize.fminbound(BOpt_objFun, Bmin, Bmax, args=(T_list, Cp_list, cp0, cpInf))
	(a0, a1, a2, a3, resid) = GA2Wilhoit(B[0], T_list, Cp_list, cp0, cpInf)
	return a0, a1, a2, a3, B[0], resid

def BOpt_objFun(B, T_list, Cp_list, cp0, cpInf):
	#input: B (in kiloKelvin), GA temperature and Cp_list (scaled/non-dimensionalized), Wilhoit parameters, Cp0/R and CpInf/R
	#output: the sum of squared errors between Wilhoit and GA data (dimensionless)
	(a0, a1, a2, a3, resid) = GA2Wilhoit(B, T_list, Cp_list, cp0, cpInf)
	return resid


def CpLimits(atoms, rotors, linear):
	"""Calculate the zero and infinity limits for heat capacity.
	
	Input: number of atoms, number of rotors, linearity (`True` if molecule is linear)
	Output: Cp(0 K)/R, Cp(infinity)/R
	"""
	#(based off of lsfp_wilh1.f in GATPFit)
	
	if (atoms == 1): # monoatomic
		cp0 = 2.5
		cpInf = 2.5
	elif linear: # linear
		cp0	 = 3.5
		cpInf = 3*atoms - 1.5
	else: # nonlinear
		cp0	 = 4.0
		cpInf = 3*atoms - (2 + 0.5*rotors)
	return cp0, cpInf

################################################################################
def convertWilhoitToNASA(Wilhoit, fixed=1, weighting=1, tint=1000.0, Tmin = 298.0, Tmax=6000.0):
	"""Convert a Wilhoit thermo instance into a NASA polynomial thermo instance.
	
	Takes a `ThermoWilhoitData` instance of themochemical data.
	Returns a `ThermoNASAData` instance containing two `ThermoNASAPolynomial` 
	polynomials
	"""
	
	# Scale the temperatures to kK
	Tmin = Tmin/1000
	tint = tint/1000
	Tmax = Tmax/1000

	# Make copy of Wilhoit data so we don't modify the original
	wilhoit_scaled = ThermoWilhoitData(Wilhoit.cp0, Wilhoit.cpInf, Wilhoit.a0, Wilhoit.a1, Wilhoit.a2, Wilhoit.a3, Wilhoit.H0, Wilhoit.S0, Wilhoit.comment)
	# Rescale Wilhoit parameters
	wilhoit_scaled.cp0 /= constants.R
	wilhoit_scaled.cpInf /= constants.R
	wilhoit_scaled.B /= 1000.
	
	#if we are using fixed tint, do not allow tint to float
	if(fixed == 1):
		nasa_low, nasa_high = Wilhoit2NASA(wilhoit_scaled, Tmin, Tmax, tint, weighting) 
	else:
		nasa_low, nasa_high, tint = Wilhoit2NASA_TintOpt(wilhoit_scaled, Tmin, Tmax, weighting)
	iseUnw = TintOpt_objFun(tint, wilhoit_scaled, Tmin, Tmax, 0) #the scaled, unweighted ISE (integral of squared error)
	rmsUnw = math.sqrt(iseUnw/(Tmax-Tmin))
	rmsStr = '(Unweighted) RMS error = %.3f*R;'%(rmsUnw)
	if(weighting == 1):
		iseWei= TintOpt_objFun(tint, wilhoit_scaled, Tmin, Tmax, weighting) #the scaled, weighted ISE
		rmsWei = math.sqrt(iseWei/math.log(Tmax/Tmin))
		rmsStr = 'Weighted RMS error = %.3f*R;'%(rmsWei)+rmsStr

	#print a warning if the rms fit is worse that 0.25*R
	if(rmsUnw > 0.25 or rmsWei > 0.25):
		logging.warning("Poor Wilhoit-to-NASA fit quality: RMS error = %.3f*R" % (rmsWei if weighting == 1 else rmsUnw))
		
	#restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
	tint=tint*1000.
	Tmin = Tmin*1000
	Tmax = Tmax*1000
	
	nasa_low.c1 /= 1000.
	nasa_low.c2 /= 1000000.
	nasa_low.c3 /= 1000000000.
	nasa_low.c4 /= 1000000000000.
	
	nasa_high.c1 /= 1000.
	nasa_high.c2 /= 1000000.
	nasa_high.c3 /= 1000000000.
	nasa_high.c4 /= 1000000000000.
	
	# output comment
	comment = 'NASA function fitted to Wilhoit function. ' + rmsStr + Wilhoit.comment
	nasa_low.Trange = (Tmin,tint); nasa_low.Tmin = Tmin; nasa_low.Tmax = tint
	nasa_low.comment = 'Low temperature range polynomial'
	nasa_high.Trange = (tint,Tmax); nasa_high.Tmin = tint; nasa_high.Tmax = Tmax
	nasa_high.comment = 'High temperature range polynomial'
	
	#for the low polynomial, we want the results to match the Wilhoit value at 298.15K
	#low polynomial enthalpy:
	Hlow = (Wilhoit.getEnthalpy(298.15) - nasa_low.getEnthalpy(298.15))/constants.R
	###polynomial_low.coeffs[5] = (Wilhoit.getEnthalpy(298.15) - polynomial_low.getEnthalpy(298.15))/constants.R
	#low polynomial entropy:
	Slow = (Wilhoit.getEntropy(298.15) - nasa_low.getEntropy(298.15))/constants.R
	###polynomial_low.coeffs[6] = (Wilhoit.getEntropy(298.15) - polynomial_low.getEntropy(298.15))/constants.R
	
	# update last two coefficients
	nasa_low.c5 = Hlow
	nasa_low.c6 = Slow
	
	#for the high polynomial, we want the results to match the low polynomial value at tint
	#high polynomial enthalpy:
	Hhigh = (nasa_low.getEnthalpy(tint) - nasa_high.getEnthalpy(tint))/constants.R
	#high polynomial entropy:
	Shigh = (nasa_low.getEntropy(tint) - nasa_high.getEntropy(tint))/constants.R
	
	# update last two coefficients
	#polynomial_high.coeffs = (b6,b7,b8,b9,b10,Hhigh,Shigh)
	nasa_high.c5 = Hhigh
	nasa_high.c6 = Shigh
	
	NASAthermo = ThermoNASAData( Trange=(Tmin,Tmax), polynomials=[nasa_low,nasa_high], comment=comment)
	return NASAthermo

################################################################################

def Wilhoit2NASA(wilhoit, tmin, tmax, tint, weighting, contCons=3):
	"""
	input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, 
	       Tmin (minimum temperature (in kiloKelvin), 
	       Tmax (maximum temperature (in kiloKelvin), 
	       Tint (intermediate temperature, in kiloKelvin)
	       weighting (boolean: should the fit be weighted by 1/T?)
	       contCons: a measure of the continutity constraints on the fitted NASA polynomials; possible values are:
		    5: constrain Cp, dCp/dT, d2Cp/dT2, d3Cp/dT3, and d4Cp/dT4 to be continuous at tint; note: this effectively constrains all the coefficients to be equal and should be equivalent to fitting only one polynomial (rather than two)
		    4: constrain Cp, dCp/dT, d2Cp/dT2, and d3Cp/dT3 to be continuous at tint
		    3 (default): constrain Cp, dCp/dT, and d2Cp/dT2 to be continuous at tint
		    2: constrain Cp and dCp/dT to be continuous at tint
		    1: constrain Cp to be continous at tint
		    0: no constraints on continuity of Cp(T) at tint
		    note: 5th (and higher) derivatives of NASA Cp(T) are zero and hence will automatically be continuous at tint by the form of the Cp(T) function
	output: NASA polynomials (nasa_low, nasa_high) with scaled parameters
	"""
	#construct (typically 13*13) symmetric A matrix (in A*x = b); other elements will be zero
	A = scipy.zeros([10+contCons,10+contCons])
	b = scipy.zeros([10+contCons])

	if weighting:
		A[0,0] = 2*math.log(tint/tmin)
		A[0,1] = 2*(tint - tmin)
		A[0,2] = tint*tint - tmin*tmin
		A[0,3] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
		A[0,4] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
		A[1,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
		A[2,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
		A[3,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
		A[4,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
	else:
		A[0,0] = 2*(tint - tmin)
		A[0,1] = tint*tint - tmin*tmin
		A[0,2] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
		A[0,3] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
		A[0,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
		A[1,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
		A[2,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
		A[3,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
		A[4,4] = 2.*(tint*tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/9
	A[1,1] = A[0,2]
	A[1,2] = A[0,3]
	A[1,3] = A[0,4]
	A[2,2] = A[0,4]
	A[2,3] = A[1,4]
	A[3,3] = A[2,4]

	if weighting:
		A[5,5] = 2*math.log(tmax/tint)
		A[5,6] = 2*(tmax - tint)
		A[5,7] = tmax*tmax - tint*tint
		A[5,8] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
		A[5,9] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
		A[6,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
		A[7,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
		A[8,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
		A[9,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
	else:
		A[5,5] = 2*(tmax - tint)
		A[5,6] = tmax*tmax - tint*tint
		A[5,7] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
		A[5,8] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
		A[5,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
		A[6,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
		A[7,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
		A[8,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
		A[9,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint*tint)/9
	A[6,6] = A[5,7]
	A[6,7] = A[5,8]
	A[6,8] = A[5,9]
	A[7,7] = A[5,9]
	A[7,8] = A[6,9]
	A[8,8] = A[7,9]

	if(contCons > 0):#set non-zero elements in the 11th column for Cp(T) continuity contraint
		A[0,10] = 1.
		A[1,10] = tint
		A[2,10] = tint*tint
		A[3,10] = A[2,10]*tint
		A[4,10] = A[3,10]*tint
		A[5,10] = -A[0,10]
		A[6,10] = -A[1,10]
		A[7,10] = -A[2,10]
		A[8,10] = -A[3,10]
		A[9,10] = -A[4,10]
		if(contCons > 1): #set non-zero elements in the 12th column for dCp/dT continuity constraint
			A[1,11] = 1.
			A[2,11] = 2*tint
			A[3,11] = 3*A[2,10]
			A[4,11] = 4*A[3,10]
			A[6,11] = -A[1,11]
			A[7,11] = -A[2,11]
			A[8,11] = -A[3,11]
			A[9,11] = -A[4,11]
			if(contCons > 2): #set non-zero elements in the 13th column for d2Cp/dT2 continuity constraint
				A[2,12] = 2.
				A[3,12] = 6*tint
				A[4,12] = 12*A[2,10]
				A[7,12] = -A[2,12]
				A[8,12] = -A[3,12]
				A[9,12] = -A[4,12]
				if(contCons > 3): #set non-zero elements in the 14th column for d3Cp/dT3 continuity constraint
					A[3,13] = 6
					A[4,13] = 24*tint
					A[8,13] = -A[3,13]
					A[9,13] = -A[4,13]
					if(contCons > 4): #set non-zero elements in the 15th column for d4Cp/dT4 continuity constraint
						A[4,14] = 24
						A[9,14] = -A[4,14]

	# make the matrix symmetric
	for i in range(1,10+contCons):
		for j in range(0, i):
			A[i,j] = A[j,i]

	#construct b vector
	w0int = wilhoit.integral_T0(tint)
	w1int = wilhoit.integral_T1(tint)
	w2int = wilhoit.integral_T2(tint)
	w3int = wilhoit.integral_T3(tint)
	w0min = wilhoit.integral_T0(tmin)
	w1min = wilhoit.integral_T1(tmin)
	w2min = wilhoit.integral_T2(tmin)
	w3min = wilhoit.integral_T3(tmin)
	w0max = wilhoit.integral_T0(tmax)
	w1max = wilhoit.integral_T1(tmax)
	w2max = wilhoit.integral_T2(tmax)
	w3max = wilhoit.integral_T3(tmax)
	if weighting:
		wM1int = wilhoit.integral_TM1(tint)
		wM1min = wilhoit.integral_TM1(tmin)
		wM1max = wilhoit.integral_TM1(tmax)
	else:
		w4int = wilhoit.integral_T4(tint)
		w4min = wilhoit.integral_T4(tmin)
		w4max = wilhoit.integral_T4(tmax)

	if weighting:
		b[0] = 2*(wM1int - wM1min)
		b[1] = 2*(w0int - w0min)
		b[2] = 2*(w1int - w1min)
		b[3] = 2*(w2int - w2min)
		b[4] = 2*(w3int - w3min)
		b[5] = 2*(wM1max - wM1int)
		b[6] = 2*(w0max - w0int)
		b[7] = 2*(w1max - w1int)
		b[8] = 2*(w2max - w2int)
		b[9] = 2*(w3max - w3int)
	else:
		b[0] = 2*(w0int - w0min)
		b[1] = 2*(w1int - w1min)
		b[2] = 2*(w2int - w2min)
		b[3] = 2*(w3int - w3min)
		b[4] = 2*(w4int - w4min)
		b[5] = 2*(w0max - w0int)
		b[6] = 2*(w1max - w1int)
		b[7] = 2*(w2max - w2int)
		b[8] = 2*(w3max - w3int)
		b[9] = 2*(w4max - w4int)

	# solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A
	# matrix is not required; not including it should give same result, except
	# Lagrange multipliers will differ by a factor of two)
	x = linalg.solve(A,b,overwrite_a=1,overwrite_b=1)

	nasa_low = ThermoNASAPolynomial(T_range=(0,0), coeffs=[x[0], x[1], x[2], x[3], x[4], 0.0, 0.0], comment='')
	nasa_high = ThermoNASAPolynomial(T_range=(0,0), coeffs=[x[5], x[6], x[7], x[8], x[9], 0.0, 0.0], comment='')

	return nasa_low, nasa_high
	
def Wilhoit2NASA_TintOpt(wilhoit, tmin, tmax, weighting):
	#input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
	#1. vary Tint, bounded by tmin and tmax, to minimize TintOpt_objFun
	#cf. http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html and http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html#scipy.optimize.fminbound)
	tint = optimize.fminbound(TintOpt_objFun, tmin, tmax, args=(wilhoit, tmin, tmax, weighting))
	#note that we have not used any guess when using this minimization routine
	#2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
	(nasa1, nasa2) = Wilhoit2NASA(wilhoit, tmin, tmax, tint[0] ,weighting)
	return nasa1, nasa2, tint[0]

def TintOpt_objFun(tint, wilhoit, tmin, tmax, weighting):
	#input: Tint (intermediate temperature, in kiloKelvin); Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	if (weighting == 1):
		result = TintOpt_objFun_W(tint, wilhoit, tmin, tmax)
	else:
		result = TintOpt_objFun_NW(tint, wilhoit, tmin, tmax)

	# numerical errors could accumulate to give a slightly negative result
	# this is unphysical (it's the integral of a *squared* error) so we
	# set it to zero to avoid later problems when we try find the square root.
	if result<0:
		logging.error("Greg thought he fixed the numerical problem, but apparently it is still an issue; please e-mail him with the following results:")
		logging.error(tint)
		logging.error(wilhoit)
		logging.error(tmin)
		logging.error(tmax)
		logging.error(weighting)
		logging.error(result)
		result = 0

	return result

def TintOpt_objFun_NW(tint, wilhoit, tmin, tmax):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.
	
	input: Tint (intermediate temperature, in kiloKelvin)
			Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, 
			Tmin (minimum temperature (in kiloKelvin), 
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Wilhoit2NASA(wilhoit,tmin,tmax,tint, 0)
	b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
	b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

	q0=wilhoit.integral_T0(tint)
	q1=wilhoit.integral_T1(tint)
	q2=wilhoit.integral_T2(tint)
	q3=wilhoit.integral_T3(tint)
	q4=wilhoit.integral_T4(tint)
	result = (wilhoit.integral2_T0(tmax) - wilhoit.integral2_T0(tmin) +
				 nasa_low.integral2_T0(tint)-nasa_low.integral2_T0(tmin) + nasa_high.integral2_T0(tmax) - nasa_high.integral2_T0(tint)
				 - 2* (b6*(wilhoit.integral_T0(tmax)-q0)+b1*(q0-wilhoit.integral_T0(tmin))
				 +b7*(wilhoit.integral_T1(tmax) - q1) +b2*(q1 - wilhoit.integral_T1(tmin))
				 +b8*(wilhoit.integral_T2(tmax) - q2) +b3*(q2 - wilhoit.integral_T2(tmin))
				 +b9*(wilhoit.integral_T3(tmax) - q3) +b4*(q3 - wilhoit.integral_T3(tmin))
				 +b10*(wilhoit.integral_T4(tmax) - q4)+b5*(q4 - wilhoit.integral_T4(tmin))))

	return result

def TintOpt_objFun_W(tint, wilhoit, tmin, tmax):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.
	
	If fit is close to perfect, result may be slightly negative due to numerical errors in evaluating this integral.
	input: Tint (intermediate temperature, in kiloKelvin)
			Wilhoit parameters: Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, 
			Tmin (minimum temperature (in kiloKelvin), 
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[1/t*(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Wilhoit2NASA(wilhoit,tmin,tmax,tint, 1)
	b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
	b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

	qM1=wilhoit.integral_TM1(tint)
	q0=wilhoit.integral_T0(tint)
	q1=wilhoit.integral_T1(tint)
	q2=wilhoit.integral_T2(tint)
	q3=wilhoit.integral_T3(tint)
	result = (wilhoit.integral2_TM1(tmax) - wilhoit.integral2_TM1(tmin) +
				 nasa_low.integral2_TM1(tint)-nasa_low.integral2_TM1(tmin) + nasa_high.integral2_TM1(tmax) - nasa_high.integral2_TM1(tint)
				 - 2* (b6*(wilhoit.integral_TM1(tmax)-qM1)+b1*(qM1 - wilhoit.integral_TM1(tmin))
				 +b7*(wilhoit.integral_T0(tmax)-q0)+b2*(q0 - wilhoit.integral_T0(tmin))
				 +b8*(wilhoit.integral_T1(tmax)-q1)+b3*(q1 - wilhoit.integral_T1(tmin))
				 +b9*(wilhoit.integral_T2(tmax)-q2)+b4*(q2 - wilhoit.integral_T2(tmin))
				 +b10*(wilhoit.integral_T3(tmax)-q3)+b5*(q3 - wilhoit.integral_T3(tmin))))

	return result

####################################################################################################
#below are functions for conversion of general Cp to NASA polynomials
#because they use numerical integration, they are, in general, likely to be slower and less accurate than versions with analytical integrals for the starting Cp form (e.g. Wilhoit polynomials)
#therefore, this should only be used when no analytic alternatives are available
def convertCpToNASA(CpObject, H298, S298, fixed=1, weighting=0, tint=1000.0, Tmin = 298.0, Tmax=6000.0):
	"""Convert an arbitrary heat capacity function into a NASA polynomial thermo instance (using numerical integration)

	Takes:  CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
		H298: enthalpy at 298.15 K (in J/mol)
		S298: entropy at 298.15 K (in J/mol-K)
		fixed: 1 (default) to fix tint; 0 to allow it to float to get a better fit
		weighting: 0 (default) to not weight the fit by 1/T; 1 to weight by 1/T to emphasize good fit at lower temperatures
		Tint, Tmin, Tmax: intermediate, minimum, and maximum temperatures in Kelvin
	Returns a `ThermoNASAData` instance containing two `ThermoNASAPolynomial` polynomials
	"""

	# Scale the temperatures to kK
	Tmin = Tmin/1000
	tint = tint/1000
	Tmax = Tmax/1000

	#if we are using fixed tint, do not allow tint to float
	if(fixed == 1):
		nasa_low, nasa_high = Cp2NASA(CpObject, Tmin, Tmax, tint, weighting)
	else:
		nasa_low, nasa_high, tint = Cp2NASA_TintOpt(CpObject, Tmin, Tmax, weighting)
	iseUnw = Cp_TintOpt_objFun(tint, CpObject, Tmin, Tmax, 0) #the scaled, unweighted ISE (integral of squared error)
	rmsUnw = math.sqrt(iseUnw/(Tmax-Tmin))
	rmsStr = '(Unweighted) RMS error = %.3f*R;'%(rmsUnw)
	if(weighting == 1):
		iseWei= Cp_TintOpt_objFun(tint, CpObject, Tmin, Tmax, weighting) #the scaled, weighted ISE
		rmsWei = math.sqrt(iseWei/math.log(Tmax/Tmin))
		rmsStr = 'Weighted RMS error = %.3f*R;'%(rmsWei)+rmsStr

	#print a warning if the rms fit is worse that 0.25*R
	if(rmsUnw > 0.25 or rmsWei > 0.25):
		logging.warning("Poor Cp-to-NASA fit quality: RMS error = %.3f*R" % (rmsWei if weighting == 1 else rmsUnw))

	#restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
	tint=tint*1000.
	Tmin = Tmin*1000
	Tmax = Tmax*1000

	nasa_low.c1 /= 1000.
	nasa_low.c2 /= 1000000.
	nasa_low.c3 /= 1000000000.
	nasa_low.c4 /= 1000000000000.

	nasa_high.c1 /= 1000.
	nasa_high.c2 /= 1000000.
	nasa_high.c3 /= 1000000000.
	nasa_high.c4 /= 1000000000000.

	# output comment
	comment = 'Cp function fitted to NASA function. ' + rmsStr
	nasa_low.Trange = (Tmin,tint); nasa_low.Tmin = Tmin; nasa_low.Tmax = tint
	nasa_low.comment = 'Low temperature range polynomial'
	nasa_high.Trange = (tint,Tmax); nasa_high.Tmin = tint; nasa_high.Tmax = Tmax
	nasa_high.comment = 'High temperature range polynomial'

	#for the low polynomial, we want the results to match the given values at 298.15K
	#low polynomial enthalpy:
	Hlow = (H298 - nasa_low.getEnthalpy(298.15))/constants.R
	#low polynomial entropy:
	Slow = (S298 - nasa_low.getEntropy(298.15))/constants.R
	#***consider changing this to use getEnthalpy and getEntropy methods of thermoObject

	# update last two coefficients
	nasa_low.c5 = Hlow
	nasa_low.c6 = Slow

	#for the high polynomial, we want the results to match the low polynomial value at tint
	#high polynomial enthalpy:
	Hhigh = (nasa_low.getEnthalpy(tint) - nasa_high.getEnthalpy(tint))/constants.R
	#high polynomial entropy:
	Shigh = (nasa_low.getEntropy(tint) - nasa_high.getEntropy(tint))/constants.R

	# update last two coefficients
	#polynomial_high.coeffs = (b6,b7,b8,b9,b10,Hhigh,Shigh)
	nasa_high.c5 = Hhigh
	nasa_high.c6 = Shigh

	NASAthermo = ThermoNASAData( Trange=(Tmin,Tmax), polynomials=[nasa_low,nasa_high], comment=comment)
	return NASAthermo

################################################################################

def Cp2NASA(CpObject, tmin, tmax, tint, weighting, contCons=3):
	"""
	input: CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
	       Tmin (minimum temperature (in kiloKelvin),
	       Tmax (maximum temperature (in kiloKelvin),
	       Tint (intermediate temperature, in kiloKelvin)
	       weighting (boolean: should the fit be weighted by 1/T?)
	       contCons: a measure of the continutity constraints on the fitted NASA polynomials; possible values are:
		    5: constrain Cp, dCp/dT, d2Cp/dT2, d3Cp/dT3, and d4Cp/dT4 to be continuous at tint; note: this effectively constrains all the coefficients to be equal and should be equivalent to fitting only one polynomial (rather than two)
		    4: constrain Cp, dCp/dT, d2Cp/dT2, and d3Cp/dT3 to be continuous at tint
		    3 (default): constrain Cp, dCp/dT, and d2Cp/dT2 to be continuous at tint
		    2: constrain Cp and dCp/dT to be continuous at tint
		    1: constrain Cp to be continous at tint
		    0: no constraints on continuity of Cp(T) at tint
		    note: 5th (and higher) derivatives of NASA Cp(T) are zero and hence will automatically be continuous at tint by the form of the Cp(T) function
	output: NASA polynomials (nasa_low, nasa_high) with scaled parameters
	"""
	#construct (typically 13*13) symmetric A matrix (in A*x = b); other elements will be zero
	A = scipy.zeros([10+contCons,10+contCons])
	b = scipy.zeros([10+contCons])

	if weighting:
		A[0,0] = 2*math.log(tint/tmin)
		A[0,1] = 2*(tint - tmin)
		A[0,2] = tint*tint - tmin*tmin
		A[0,3] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
		A[0,4] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
		A[1,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
		A[2,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
		A[3,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
		A[4,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
	else:
		A[0,0] = 2*(tint - tmin)
		A[0,1] = tint*tint - tmin*tmin
		A[0,2] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
		A[0,3] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
		A[0,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
		A[1,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
		A[2,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
		A[3,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
		A[4,4] = 2.*(tint*tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/9
	A[1,1] = A[0,2]
	A[1,2] = A[0,3]
	A[1,3] = A[0,4]
	A[2,2] = A[0,4]
	A[2,3] = A[1,4]
	A[3,3] = A[2,4]

	if weighting:
		A[5,5] = 2*math.log(tmax/tint)
		A[5,6] = 2*(tmax - tint)
		A[5,7] = tmax*tmax - tint*tint
		A[5,8] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
		A[5,9] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
		A[6,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
		A[7,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
		A[8,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
		A[9,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
	else:
		A[5,5] = 2*(tmax - tint)
		A[5,6] = tmax*tmax - tint*tint
		A[5,7] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
		A[5,8] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
		A[5,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
		A[6,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
		A[7,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
		A[8,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
		A[9,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint*tint)/9
	A[6,6] = A[5,7]
	A[6,7] = A[5,8]
	A[6,8] = A[5,9]
	A[7,7] = A[5,9]
	A[7,8] = A[6,9]
	A[8,8] = A[7,9]

	if(contCons > 0):#set non-zero elements in the 11th column for Cp(T) continuity contraint
		A[0,10] = 1.
		A[1,10] = tint
		A[2,10] = tint*tint
		A[3,10] = A[2,10]*tint
		A[4,10] = A[3,10]*tint
		A[5,10] = -A[0,10]
		A[6,10] = -A[1,10]
		A[7,10] = -A[2,10]
		A[8,10] = -A[3,10]
		A[9,10] = -A[4,10]
		if(contCons > 1): #set non-zero elements in the 12th column for dCp/dT continuity constraint
			A[1,11] = 1.
			A[2,11] = 2*tint
			A[3,11] = 3*A[2,10]
			A[4,11] = 4*A[3,10]
			A[6,11] = -A[1,11]
			A[7,11] = -A[2,11]
			A[8,11] = -A[3,11]
			A[9,11] = -A[4,11]
			if(contCons > 2): #set non-zero elements in the 13th column for d2Cp/dT2 continuity constraint
				A[2,12] = 2.
				A[3,12] = 6*tint
				A[4,12] = 12*A[2,10]
				A[7,12] = -A[2,12]
				A[8,12] = -A[3,12]
				A[9,12] = -A[4,12]
				if(contCons > 3): #set non-zero elements in the 14th column for d3Cp/dT3 continuity constraint
					A[3,13] = 6
					A[4,13] = 24*tint
					A[8,13] = -A[3,13]
					A[9,13] = -A[4,13]
					if(contCons > 4): #set non-zero elements in the 15th column for d4Cp/dT4 continuity constraint
						A[4,14] = 24
						A[9,14] = -A[4,14]

	# make the matrix symmetric
	for i in range(1,10+contCons):
		for j in range(0, i):
			A[i,j] = A[j,i]

	#construct b vector
	w0low = Nintegral_T0(CpObject,tmin,tint)
	w1low = Nintegral_T1(CpObject,tmin,tint)
	w2low = Nintegral_T2(CpObject,tmin,tint)
	w3low = Nintegral_T3(CpObject,tmin,tint)
	w0high = Nintegral_T0(CpObject,tint,tmax)
	w1high = Nintegral_T1(CpObject,tint,tmax)
	w2high = Nintegral_T2(CpObject,tint,tmax)
	w3high = Nintegral_T3(CpObject,tint,tmax)
	if weighting:
		wM1low = Nintegral_TM1(CpObject,tmin,tint)
		wM1high = Nintegral_TM1(CpObject,tint,tmax)
	else:
		w4low = Nintegral_T4(CpObject,tmin,tint)
		w4high = Nintegral_T4(CpObject,tint,tmax)

	if weighting:
		b[0] = 2*wM1low
		b[1] = 2*w0low
		b[2] = 2*w1low
		b[3] = 2*w2low
		b[4] = 2*w3low
		b[5] = 2*wM1high
		b[6] = 2*w0high
		b[7] = 2*w1high
		b[8] = 2*w2high
		b[9] = 2*w3high
	else:
		b[0] = 2*w0low
		b[1] = 2*w1low
		b[2] = 2*w2low
		b[3] = 2*w3low
		b[4] = 2*w4low
		b[5] = 2*w0high
		b[6] = 2*w1high
		b[7] = 2*w2high
		b[8] = 2*w3high
		b[9] = 2*w4high

	# solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A
	# matrix is not required; not including it should give same result, except
	# Lagrange multipliers will differ by a factor of two)
	x = linalg.solve(A,b,overwrite_a=1,overwrite_b=1)

	nasa_low = ThermoNASAPolynomial(T_range=(0,0), coeffs=[x[0], x[1], x[2], x[3], x[4], 0.0, 0.0], comment='')
	nasa_high = ThermoNASAPolynomial(T_range=(0,0), coeffs=[x[5], x[6], x[7], x[8], x[9], 0.0, 0.0], comment='')

	return nasa_low, nasa_high

def Cp2NASA_TintOpt(CpObject, tmin, tmax, weighting):
	#input: CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
	#1. vary Tint, bounded by tmin and tmax, to minimize TintOpt_objFun
	#cf. http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html and http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html#scipy.optimize.fminbound)
	tint = optimize.fminbound(Cp_TintOpt_objFun, tmin, tmax, args=(CpObject, tmin, tmax, weighting))
	#note that we have not used any guess when using this minimization routine
	#2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
	(nasa1, nasa2) = Cp2NASA(CpObject, tmin, tmax, tint[0] ,weighting)
	return nasa1, nasa2, tint[0]

def Cp_TintOpt_objFun(tint, CpObject, tmin, tmax, weighting):
	#input: Tint (intermediate temperature, in kiloKelvin); CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	if (weighting == 1):
		result = Cp_TintOpt_objFun_W(tint, CpObject, tmin, tmax)
	else:
		result = Cp_TintOpt_objFun_NW(tint, CpObject, tmin, tmax)

	# numerical errors could accumulate to give a slightly negative result
	# this is unphysical (it's the integral of a *squared* error) so we
	# set it to zero to avoid later problems when we try find the square root.
	if result<0:
		logging.error("Numerical integral results suggest sum of squared errors is negative; please e-mail Greg with the following results:")
		logging.error(tint)
		logging.error(CpObject)
		logging.error(tmin)
		logging.error(tmax)
		logging.error(weighting)
		logging.error(result)
		result = 0

	return result

def Cp_TintOpt_objFun_NW(tint, CpObject, tmin, tmax):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.

	input: Tint (intermediate temperature, in kiloKelvin)
			CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
			Tmin (minimum temperature (in kiloKelvin),
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[(Cp/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Cp2NASA(CpObject,tmin,tmax,tint, 0)
	b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
	b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

	result = (Nintegral2_T0(CpObject,tmin,tmax) +
				 nasa_low.integral2_T0(tint)-nasa_low.integral2_T0(tmin) + nasa_high.integral2_T0(tmax) - nasa_high.integral2_T0(tint)
				 - 2* (b6*Nintegral_T0(CpObject,tint,tmax)+b1*Nintegral_T0(CpObject,tmin,tint)
				 +b7*Nintegral_T1(CpObject,tint,tmax) +b2*Nintegral_T1(CpObject,tmin,tint)
				 +b8*Nintegral_T2(CpObject,tint,tmax) +b3*Nintegral_T2(CpObject,tmin,tint)
				 +b9*Nintegral_T3(CpObject,tint,tmax) +b4*Nintegral_T3(CpObject,tmin,tint)
				 +b10*Nintegral_T4(CpObject,tint,tmax)+b5*Nintegral_T4(CpObject,tmin,tint)))

	return result

def Cp_TintOpt_objFun_W(tint, CpObject, tmin, tmax):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.

	If fit is close to perfect, result may be slightly negative due to numerical errors in evaluating this integral.
	input: Tint (intermediate temperature, in kiloKelvin)
			CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
			Tmin (minimum temperature (in kiloKelvin),
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[1/t*(Cp/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Cp2NASA(CpObject,tmin,tmax,tint, 1)
	b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
	b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

	result = (Nintegral2_TM1(CpObject,tmin,tmax) +
				 nasa_low.integral2_TM1(tint)-nasa_low.integral2_TM1(tmin) + nasa_high.integral2_TM1(tmax) - nasa_high.integral2_TM1(tint)
				 - 2* (b6*Nintegral_TM1(CpObject,tint,tmax)+b1*Nintegral_TM1(CpObject,tmin,tint)
				 +b7*Nintegral_T0(CpObject,tint,tmax) +b2*Nintegral_T0(CpObject,tmin,tint)
				 +b8*Nintegral_T1(CpObject,tint,tmax) +b3*Nintegral_T1(CpObject,tmin,tint)
				 +b9*Nintegral_T2(CpObject,tint,tmax) +b4*Nintegral_T2(CpObject,tmin,tint)
				 +b10*Nintegral_T3(CpObject,tint,tmax)+b5*Nintegral_T3(CpObject,tmin,tint)))

	return result

#the numerical integrals:

def Nintegral_T0(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,0,0)

def Nintegral_TM1(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,-1,0)

def Nintegral_T1(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,1,0)

def Nintegral_T2(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,2,0)

def Nintegral_T3(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,3,0)

def Nintegral_T4(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,4,0)

def Nintegral2_T0(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,0,1)

def Nintegral2_TM1(CpObject, tmin, tmax):
	#units of input and output are same as Nintegral
	return Nintegral(CpObject,tmin,tmax,-1,1)

def Nintegral(CpObject, tmin, tmax, n, squared):
	#inputs:CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
	#	tmin, tmax: limits of integration in kiloKelvin
	#	n: integeer exponent on t (see below), typically -1 to 4
	#	squared: 0 if integrating Cp/R(t)*t^n; 1 if integrating Cp/R(t)^2*t^n
	#output: a numerical approximation to the quantity Integrate[Cp/R(t)*t^n, {t, tmin, tmax}] or Integrate[Cp/R(t)^2*t^n, {t, tmin, tmax}], in units based on kiloKelvin

	return integrate.quad(integrand,tmin,tmax,args=(CpObject,n,squared))[0]

def integrand(t, CpObject , n, squared):
	#input requirements same as Nintegral above
	result = CpObject.getHeatCapacity(t*1000)/constants.R#note that we multiply t by 1000, since the Cp function uses Kelvin rather than kiloKelvin; also, we divide by R to get the dimensionless Cp/R
	if(squared):
		result = result*result
	if(n < 0):
		for i in range(0,abs(n)):#divide by t, |n| times
			result = result/t
	else:
		for i in range(0,n):#multiply by t, n times
			result = result*t
	return result
################################################################################

if __name__ == '__main__':
	pass
