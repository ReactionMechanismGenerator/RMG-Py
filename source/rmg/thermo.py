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
"""

import quantities as pq
import constants
import data
import math
import scipy
from scipy import linalg
from scipy import optimize

################################################################################

class ThermoData:
	"""
	A base class for all forms of thermodynamic data used by RMG. The common
	attributes are:
	
	========= ============================================================
	Attribute Meaning
	========= ============================================================
	`Trange`  a 2-tuple containin the minimum and maximum temperature in K
	`comment` a string describing the source of the data
	========= ============================================================
	"""
	
	def __init__(self, Trange=None, comment=''):
		self.Trange = Trange or (0.0, 0.0)
		self.comment = comment
	
	def isTemperatureValid(self, T):
		"""
		Return :data:`True` if the temperature `T` in K is within the valid
		temperature range of the thermodynamic data, or :data:`False` if not.
		"""
		if self.Trange is None:
			return True
		else:
			return self.Trange[0] <= T and T <= self.Trange[1]

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
	`Cp`      the standard heat capacity in J/mol*K at 300, 400, 500,
		      600, 800, 1000, and 1500 K
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
		self.__cache_T=None
		self.__cache = dict()
	
	def __add__(self, other):
		"""
		Add two sets of thermodynamic data together. All parameters are
		considered additive. Returns a new :class:`ThermoGAData` object that is
		the sum of the two sets of thermodynamic data.
		"""
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
	
	def getHeatCapacity(self, T):
		"""
		Return the heat capacity in J/mol*K at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from group additivity.')
		if T < 300.0:
			return self.Cp[0]
		elif T > ThermoGAData.CpTlist[-1]:
			return self.Cp[-1]
		else:
			for Tmin, Tmax, Cpmin, Cpmax in zip(ThermoGAData.CpTlist[:-1], \
					ThermoGAData.CpTlist[1:], self.Cp[:-1], self.Cp[1:]):
				if Tmin <= T and T <= Tmax:
					return (Cpmax - Cpmin) * ((T - Tmin) / (Tmax - Tmin)) + Cpmin
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""	
		if self.__cache_T==T:
			try: return self.__cache['H']
			except KeyError: pass # 'H' isn't cached yet
		else:  # temperature has changed
			self.__cache = dict()
			self.__cache_T = T
			
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
		self.__cache['H'] = H
		return H
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		import math
		
		if self.__cache_T==T:
			try: return self.__cache['S']
			except KeyError: pass # 'S' isn't cached yet
		else:  # temperature has changed
			self.__cache = dict()
			self.__cache_T = T
			
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
		self.__cache['S'] = S
		return S
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		if self.__cache_T==T:
			try: return self.__cache['G']
			except KeyError: pass # 'G' isn't cached yet
		else:  # temperature has changed
			self.__cache = dict()
			self.__cache_T = T
			
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from group additivity.')
		G = self.getEnthalpy(T) - T * self.getEntropy(T)
		self.__cache['G'] = G
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
	
	def toXML(self, dom, root):
		"""
		Generate an XML representation of the thermodynamic data using the
		:data:`xml.dom.minidom` package. The `dom` and `root` parameters
		represent the DOM and the element in the DOM used as the parent of
		the generated XML.
		"""
		
		thermo = dom.createElement('thermodynamics')
		thermo.setAttribute('index', self.index)
		thermo.setAttribute('comment', self.comment)
		root.appendChild(thermo)
		
		enthalpy = dom.createElement('enthalpyOfFormation')
		enthalpy.setAttribute('temperature', '298.0 K')
		thermo.appendChild(enthalpy)
		data.createXMLQuantity(dom, enthalpy, self.H298, 'J/mol')
		
		entropy = dom.createElement('entropyOfFormation')
		entropy.setAttribute('temperature', '298.0 K')
		thermo.appendChild(entropy)
		data.createXMLQuantity(dom, entropy, self.S298, 'J/(mol*K)')
		
		for i, Cp in enumerate(self.Cp):
			heatCapacity = dom.createElement('heatCapacity')
			heatCapacity.setAttribute('temperature', '%s K' % (ThermoGAData.CpTlist[i]) )
			thermo.appendChild(heatCapacity)
			data.createXMLQuantity(dom, heatCapacity, Cp, 'J/(mol*K)')
	

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
		self.coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
	
	def getHeatCapacity(self, T):
		"""
		Return the heat capacity in J/mol*K at temperature `T` in K.
		"""
		import constants
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for heat capacity estimation from NASA polynomial.')
		c = self.coeffs
		T2 = T*T
		# Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
#		HeatCapacityOverR = c[0] + c[1]*T + c[2]*T*T + c[3]*T*T*T + c[4]*T*T*T*T
		HeatCapacityOverR = c[0] + T*(c[1] + T*(c[2] + T*(c[3] + c[4]*T)))
		Cp = HeatCapacityOverR * constants.R
		return Cp
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		import constants
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for enthalpy estimation from NASA polynomial.')
		c = self.coeffs
		T2 = T*T
		T4 = T2*T2
		# H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
		EnthalpyOverR = ( c[0]*T + c[1]*T2/2 + c[2]*T2*T/3 + c[3]*T4/4 +
						  c[4]*T4*T/5 + c[5] )
		H = EnthalpyOverR * constants.R
		return H
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		import math
		import constants
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for entropy estimation from NASA polynomial.')
		# S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
		c = self.coeffs
		T2 = T*T
		T4 = T2*T2
		EntropyOverR = ( c[0]*math.log(T) + c[1]*T + c[2]*T2/2 +
					c[3]*T2*T/3 + c[4]*T4/4 + c[6] )
		S = EntropyOverR * constants.R
		return S
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		if not self.isTemperatureValid(T):
			raise data.TemperatureOutOfRangeException('Invalid temperature for free energy estimation from NASA polynomial.')
		return self.getEnthalpy(T) - T * self.getEntropy(T)
	
	def toXML(self, dom, root):
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
	
	def getHeatCapacity(self, T):
		"""
		Return the heat capacity in J/mol*K at temperature `T` in K.
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
	

################################################################################

class ThermoWilhoitData(ThermoData):
	"""
	A set of thermodynamic parameters given by Wilhoit polynomials. 
	
	"""
	B = 500 # Kelvin. Default temperature.
	
	def __init__(self, cp0, cpInf, a0, a1, a2, a3, I, J, comment='', ):
		"""Initialise the Wilhoit polynomial. Trange is set to (0,9999.9)"""
		Trange = (0,9999.9) # Wilhoit valid over all temperatures
		ThermoData.__init__(self, Trange=Trange, comment=comment)
		self.cp0 = cp0
		self.cpInf = cpInf
		self.B = ThermoWilhoitData.B
		self.a0 = a0
		self.a1 = a1
		self.a2 = a2
		self.a3 = a3
		self.I = I
		self.J = J
	
	def __repr__(self):
		return "ThermoWilhoitData(%.2g,%.2g,%.2g,%.2g,%.2g,%.2g,%.2g,%.2g,%.2g'%s')"%(self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3, self.I, self.J, self.comment)
	
	def toXML(self, dom, root):
		pass
	
	def getHeatCapacity(self, T):
		"""
		Return the heat capacity in J/mol*K at temperature `T` in K.
		"""
		y = T/(T+self.B)
		Cp = self.cp0+(self.cpInf-self.cp0)*y*y*( 1 + 
			(y-1)*(self.a0 + y*(self.a1 + y*(self.a2 + y*self.a3))) )
		return Cp
	
	def getEnthalpy(self, T):
		"""
		Return the enthalpy in J/mol at temperature `T` in K.
		"""
		#cp0 = self.cp0 * R
		#cpInf = self.cpInf * R
		cp0 = self.cp0
		cpInf = self.cpInf
		B = self.B
		a0 = self.a0
		a1 = self.a1
		a2 = self.a2
		a3 = self.a3
		I = self.I
		
		enthalpy = I + WilhoitInt0(cp0, cpInf, B, a0, a1, a2, a3, T)
		
		return enthalpy
	
	def getEntropy(self, T):
		"""
		Return the entropy in J/mol*K at temperature `T` in K.
		"""
		#cp0 = self.cp0 * R
		#cpInf = self.cpInf * R
		cp0 = self.cp0
		cpInf = self.cpInf
		B = self.B
		a0 = self.a0
		a1 = self.a1
		a2 = self.a2
		a3 = self.a3
		J = self.J
		
		entropy = J + WilhoitIntM1(cp0, cpInf, B, a0, a1, a2, a3, T)
		
		return entropy
	
	def getFreeEnergy(self, T):
		"""
		Return the Gibbs free energy in J/mol at temperature `T` in K.
		"""
		return self.getEnthalpy(T) - T*self.getEntropy(T)
	

###########
def convertGAtoWilhoit(GAthermo, atoms, rotors, linear):
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
	cp = GAthermo.Cp
	t = ThermoGAData.CpTlist  # usually [300, 400, 500, 600, 800, 1000, 1500] but why assume?
	R = constants.R
	
	#convert from K to kK
	t = [x/1000. for x in t] 
	
	B = ThermoWilhoitData.B # Constant, set once in the class def.
	B = B/1000.  
	
	cp = [x/R for x in cp] #convert to Cp/R
	
	(cp0, cpInf)=CpLimits(atoms, rotors, linear)#determine the heat capacity limits (non-dimensional)
	
	#create matrices for linear least squares problem
	m = len(t)
	#A = mx4
	#b = mx1
	#x = 4x1
	A = scipy.zeros([m,4])
	b = scipy.zeros([m,1])
	for i in range(m):
		y = t[i]/(t[i]+B)
		A[i,0] = y*y*(y-1)*(cpInf-cp0)
		A[i,1] = A[i,0]*y
		A[i,2] = A[i,1]*y
		A[i,3] = A[i,2]*y
		b[i] = cp[i]-cp0 - y*y*(cpInf-cp0)
	#solve least squares problem A*x = b; http://docs.scipy.org/doc/scipy/reference/tutorial/linalg.html#solving-linear-least-squares-problems-and-pseudo-inverses
	x,resid,rank,sigma = linalg.lstsq(A,b)
	a0 = x[0]
	a1 = x[1]
	a2 = x[2]
	a3 = x[3]
	
	# err = rmsErrWilhoit(t, cp, cp0, cpInf, a0, a1, a2, a3) #(optional); display rmsError (dimensionless units)
	
	# scale everything back
	#t = [x*1000. for x in t] #not needed
	# B = B*1000. # not needed
	#cp = [x*R for x in cp] #not needed
	#cp0 and cpInf will be in units of J/mol-K
	cp0 = cp0*R
	cpInf = cpInf*R
	
	# output comment; ****we could also include fitting accuracy ("err") in the output below
	comment = 'Fitted to GA data with Cp0=%2g and Cp_inf=%2g. '%(cp0,cpInf) + GAthermo.comment
	
	# first create an instance with I=J=0, then calculate what they should be 
	# by referring to H298, S298
	I = 0
	J = 0
	
	# create Wilhoit instances
	WilhoitThermo = ThermoWilhoitData( cp0, cpInf, a0, a1, a2, a3, I, J, comment=comment)

	# calculate I, J (for H, S, respectively); self.getX(T) should return results
	# with I, J = 0, thereby allowing easy solution of the additive parameters I, J 
	I = H298 - WilhoitThermo.getEnthalpy(298.15)
	J = S298 - WilhoitThermo.getEntropy(298.15)
	
	WilhoitThermo.I = I
	WilhoitThermo.J = J
	
	return WilhoitThermo

def CpLimits(atoms, rotors, linear):
	#(based off of lsfp_wilh1.f in GATPFit)
	#input: number of atoms, number of rotors, and linearity (False for non-linear, True for linear)
	#output: Cp0/R, CpInf/R
	#monatomic case
	if(atoms == 1):
		cp0 = 2.5
		cpInf = 2.5
	#non-linear case
	elif(not linear):
		cp0	 = 4.0
		cpInf=3.*atoms-(2.+0.5*rotors)
		#linear case
	else:
		cp0	 = 3.5
		cpInf=3.*atoms-1.5
	
	return cp0, cpInf

def rmsErrWilhoit(self,t):
	#calculate the RMS error between the Wilhoit form and training data points; result will have same units as cp inputs; cp, cp0, and cpInf should agree in units (e.g. Cp/R); units of B and t should be consistent, based, for example on kK or K 
	m = len(t)
	rms = 0.0
	for i in range(m):
		err = cp[i]-getHeatCapacity(self,t[i])
		rms += err*err
	rms = rms/m
	rms = math.sqrt(rms)
	
	return rms

#this function already exists in getHeatCapacity
#def Wilhoit_Cp(t, cp0, cpInf, B, a0, a1, a2, a3):
#	#calculate Cp/R based on Wilhoit form;  result will have same units as cp0 and cpInf (e.g. Cp/R); cp is Cp/R; units of B and t should be consistent, based, for example on kK or K
#	y = t/(t+B)
#	cp = cp0+(cpInf-cp0)*y*y*(1+(y-1)*(a0+a1*y+a2*y*y+a3*y*y*y))
#	
#	return cp

################################################################################
def convertWilhoitToNASA(Wilhoit):
	"""Convert a Wilhoit thermo instance into a NASA polynomial thermo instance.
	
	Takes a `ThermoWilhoitData` instance of themochemical data.
	Returns a `ThermoNASAData` instance containing two `ThermoNASAPolynomial` 
	polynomials
	"""
	
	# Temperature ranges for resulting polynomials
	Tmin = 298.0
	Tintg = 1000.0
	Tmax = 6000.0
	
	#for now, do not allow tint to float so tint = Tint(guess)
	fixed = 1
	
	#for now, use weighting of 1/T
	weighting = 1
	
	# get info from incoming Wilhoit thermo
	cp0 = Wilhoit.cp0
	cpInf = Wilhoit.cpInf
	B = Wilhoit.B
	a0 = Wilhoit.a0
	a1 = Wilhoit.a1
	a2 = Wilhoit.a2
	a3 = Wilhoit.a3
	R = constants.R
	
	#scale the values to kK and dimensionless (for heat capacity)
	cp0 = cp0/R
	cpInf = cpInf/R
	B = B/1000
	Tmin = Tmin/1000
	Tintg = Tintg/1000
	Tmax = Tmax/1000
	
	# do your clever maths here
	#if we are using fixed tint, set tint equal to Tintg and do not allow tint to float
	if(fixed == 1):
		(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA(cp0, cpInf, B, a0, a1, a2, a3, Tmin, Tmax, Tintg, weighting)
		err = TintOpt_objFun(Tintg, cp0, cpInf, B, a0, a1, a2, a3, Tmin, Tmax, weighting) #to print the objective function value
		tint = Tintg
	else:
		(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint) = Wilhoit2NASA_TintOpt(cp0, cpInf, B, a0, a1, a2, a3, Tmin, Tmax, Tintg, weighting)
		# rmsErr = rmsErrNASA(t, cp, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint) #this needs group data
		
	#restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
	tint=tint*1000.
	b2 = b2/1000.
	b7 = b7/1000.
	b3 = b3/1000000.
	b8 = b8/1000000.
	b4 = b4/1000000000.
	b9 = b9/1000000000.
	b5 = b5/1000000000000.
	b10= b10/1000000000000.
	
	#for now, set H and S parameters equal to zero; this will need to be fixed****
	Hlow = 0.0
	Slow = 0.0
	Hhigh = 0.0
	Shigh = 0.0
	
	coeffs_low = (Hlow,Slow,b1,b2,b3,b4,b5)
	coeffs_high = (Hhigh,Shigh,b6,b7,b8,b9,b10)
	
	# could we include fitting accuracy in the expression below?
	# output comment
	comment = 'Fitted to Wilhoit data. '+Wilhoit.comment
		
	# create ThermoNASAPolynomial instances
	polynomial_low = ThermoNASAPolynomial( T_range=(Tmin,tint), comment=comment, coeffs=coeffs_low)
	polynomial_high = ThermoNASAPolynomial( T_range=(tint,Tmax), comment=comment, coeffs=coeffs_high)
	
	NASAthermo = ThermoNASAData( Trange=(Tmin,Tmax), polynomials=[polynomial_low,polynomial_high], comment=comment)
	return NASAthermo

################################################################################
def rmsErrNASA(t, cp, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint):
	#calculate the RMS error between the NASA polynomial and training data points in non-dimensional units (/R); cp is Cp/R; units of tint, t, and bi should be consistent, based, for example on kK or K 
	m = len(t)
	rms = 0.0
	for i in range(m):
		err = cp[i]-NASA_CpR(t[i],b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint)
		rms += err*err
	rms = rms/m
	rms = math.sqrt(rms)
	
	return rms


def NASA_CpR(t, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint):
	#calculate Cp/R based on NASA polynomial; cp is Cp/R; units of tint, t, and bi should be consistent, based, for example on kK or K; does not take into account lower and upper temperature limits, and will extrapolate
	if(t < tint):
		cp = b1 + b2*t + b3*t*t + b4*t*t*t + b5*t*t*t*t
	else:
		cp = b6 + b7*t + b8*t*t + b9*t*t*t + b10*t*t*t*t
	
	return cp



def Wilhoit2NASA(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tint, weighting):
	#input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin), Tint (intermediate temperature, in kiloKelvin)
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters)
	if (weighting == 1):
		(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA_W(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tint)
	else:
		(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA_NW(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tint)
	return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10


def Wilhoit2NASA_NW(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tint):
	#this is the case with no weighting
	#input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin), Tint (intermediate temperature, in kiloKelvin)
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters)
	
	#construct 13*13 symmetric A matrix (in A*x = b); other elements will be zero
	A = scipy.zeros([13,13])
	b = scipy.zeros([13,1])
	A[0,0] = 2*(tint - tmin)
	A[0,1] = tint*tint - tmin*tmin
	A[0,2] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
	A[0,3] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
	A[0,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
	A[1,1] = A[0,2]
	A[1,2] = A[0,3]
	A[1,3] = A[0,4]
	A[1,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
	A[2,2] = A[0,4]
	A[2,3] = A[1,4]
	A[2,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
	A[3,3] = A[2,4]
	A[3,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
	A[4,4] = 2.*(tint*tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/9
	
	A[5,5] = 2*(tmax - tint)
	A[5,6] = tmax*tmax - tint*tint
	A[5,7] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
	A[5,8] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
	A[5,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
	A[6,6] = A[5,7]
	A[6,7] = A[5,8]
	A[6,8] = A[5,9]
	A[6,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
	A[7,7] = A[5,9]
	A[7,8] = A[6,9]
	A[7,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
	A[8,8] = A[7,9]
	A[8,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
	A[9,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint*tint)/9
	
	A[0,10] = 1.
	A[1,10] = tint
	A[1,11] = 1.
	A[2,10] = tint*tint
	A[2,11] = 2*tint
	A[2,12] = 2.
	A[3,10] = A[2,10]*tint
	A[3,11] = 3*A[2,10]
	A[3,12] = 6*tint
	A[4,10] = A[3,10]*tint
	A[4,11] = 4*A[3,10]
	A[4,12] = 12*A[2,10]
	
	A[5,10] = -A[0,10]
	A[6,10] = -A[1,10]
	A[6,11] = -A[1,11]
	A[7,10] = -A[2,10]
	A[7,11] = -A[2,11]
	A[7,12] = -A[2,12]
	A[8,10] = -A[3,10]
	A[8,11] = -A[3,11]
	A[8,12] = -A[3,12]
	A[9,10] = -A[4,10]
	A[9,11] = -A[4,11]
	A[9,12] = -A[4,12]
	
	# make the matrix symmetric
	for i in range(1,13):
		for j in range(0, i):
			A[i,j] = A[j,i]
			
	#construct b vector
	#store values at tint (this will avoid evaluating them twice)
	w0int = WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w1int = WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w2int = WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w3int = WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w4int = WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tint)

	b[0] = 2*(w0int - WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[1] = 2*(w1int - WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[2] = 2*(w2int - WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[3] = 2*(w3int - WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[4] = 2*(w4int - WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[5] = 2*(WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w0int)
	b[6] = 2*(WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w1int)
	b[7] = 2*(WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w2int)
	b[8] = 2*(WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w3int)
	b[9] = 2*(WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w4int)
	# b[10] = 0
	# b[11] = 0
	# b[12] = 0
	
	#solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A matrix is not required; not including it should give same result, except Lagrange multipliers will differ by a factor of two)
	#from linalg import solve
	#print A
	x = linalg.solve(A,b,overwrite_a=1,overwrite_b=1)
	b1 = x[0]
	b2 = x[1]
	b3 = x[2]
	b4 = x[3]
	b5 = x[4]
	b6 = x[5]
	b7 = x[6]
	b8 = x[7]
	b9 = x[8]
	b10 = x[9]
	
	return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10

def Wilhoit2NASA_W(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tint):
	#this is the case WITH weighting
	#input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin), Tint (intermediate temperature, in kiloKelvin)
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters)
	
	#construct 13*13 symmetric A matrix (in A*x = b); other elements will be zero
	A = scipy.zeros([13,13])
	b = scipy.zeros([13,1])
	
	A[0,0] = 2*math.log(tint/tmin)
	A[0,1] = 2*(tint - tmin)
	A[0,2] = tint*tint - tmin*tmin
	A[0,3] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
	A[0,4] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2	 
	A[1,1] = A[0,2]
	A[1,2] = A[0,3]
	A[1,3] = A[0,4]
	A[1,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
	A[2,2] = A[0,4]
	A[2,3] = A[1,4]
	A[2,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
	A[3,3] = A[2,4]
	A[3,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
	A[4,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
	
	A[5,5] = 2*math.log(tmax/tint)
	A[5,6] = 2*(tmax - tint)
	A[5,7] = tmax*tmax - tint*tint
	A[5,8] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
	A[5,9] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
	A[6,6] = A[5,7]
	A[6,7] = A[5,8]
	A[6,8] = A[5,9]
	A[6,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
	A[7,7] = A[5,9]
	A[7,8] = A[6,9]
	A[7,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
	A[8,8] = A[7,9]
	A[8,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
	A[9,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
	
	
	A[0,10] = 1.
	A[1,10] = tint
	A[1,11] = 1.
	A[2,10] = tint*tint
	A[2,11] = 2*tint
	A[2,12] = 2.
	A[3,10] = A[2,10]*tint
	A[3,11] = 3*A[2,10]
	A[3,12] = 6*tint
	A[4,10] = A[3,10]*tint
	A[4,11] = 4*A[3,10]
	A[4,12] = 12*A[2,10]
	
	A[5,10] = -A[0,10]
	A[6,10] = -A[1,10]
	A[6,11] = -A[1,11]
	A[7,10] = -A[2,10]
	A[7,11] = -A[2,11]
	A[7,12] = -A[2,12]
	A[8,10] = -A[3,10]
	A[8,11] = -A[3,11]
	A[8,12] = -A[3,12]
	A[9,10] = -A[4,10]
	A[9,11] = -A[4,11]
	A[9,12] = -A[4,12]
	
	# make the matrix symmetric
	for i in range(1,13):
		for j in range(0, i):
			A[i,j] = A[j,i]
	
	#construct b vector
	#store values at tint (this will avoid evaluating them twice)
	wM1int = WilhoitIntM1(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w0int = WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w1int = WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w2int = WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tint)
	w3int = WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tint)


	b[0] = 2*(wM1int - WilhoitIntM1(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[1] = 2*(w0int - WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[2] = 2*(w1int - WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[3] = 2*(w2int - WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[4] = 2*(w3int - WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmin))
	b[5] = 2*(WilhoitIntM1(cp0,cpInf,B,a0,a1,a2,a3,tmax) - wM1int)
	b[6] = 2*(WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w0int)
	b[7] = 2*(WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w1int)
	b[8] = 2*(WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w2int)
	b[9] = 2*(WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w3int)
	# b[10] = 0
	# b[11] = 0
	# b[12] = 0

	#solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A matrix is not required; not including it should give same result, except Lagrange multipliers will differ by a factor of two)
	#from linalg import solve
	#print A
	x = linalg.solve(A,b,overwrite_a=1,overwrite_b=1)
	b1 = x[0]
	b2 = x[1]
	b3 = x[2]
	b4 = x[3]
	b5 = x[4]
	b6 = x[5]
	b7 = x[6]
	b8 = x[7]
	b9 = x[8]
	b10 = x[9]
	
	return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
	
def Wilhoit2NASA_TintOpt(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tintg, weighting):
	#input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin), Tintg (guess intermediate temperature, in kiloKelvin)
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
	#1. vary Tint, using Tintg as a starting guess, to minimize TintOpt_objFun
	#from optimize import fminbound
	tint = optimize.fminbound(TintOpt_objFun, tmin, tmax, args=(cp0, cpInf,B,a0,a1,a2,a3,tmin,tmax,weighting))
	#note that we have not used the specified guess, tintg when using this minimization routine
	#2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
	(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA(cp0,cpInf,B,a0,a1,a2,a3,tmin,tmax,tint,weighting)
	return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint

def TintOpt_objFun(tint, cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax,weighting):
	#input: Tint (intermediate temperature, in kiloKelvin); Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	if (weighting == 1):
		result = TintOpt_objFun_W(tint, cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax)
	else:
		result = TintOpt_objFun_NW(tint, cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax)
	print tint
	print result
	return result

def TintOpt_objFun_NW(tint, cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax):
	#input: Tint (intermediate temperature, in kiloKelvin); Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA(cp0,cpInf,B,a0,a1,a2,a3,tmin,tmax,tint, 0)
	result = (Wilhoit2Int(cp0,cpInf,B,a0,a1,a2,a3,tmax) - Wilhoit2Int(cp0,cpInf,B,a0,a1,a2,a3,tmin) +
				 NASA2Int(b1,b2,b3,b4,b5,tint)-NASA2Int(b1,b2,b3,b4,b5,tmin) + NASA2Int(b6,b7,b8,b9,b10,tmax) - NASA2Int(b6,b7,b8,b9,b10,tint)
				 - 2* (b6*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b1-b6)*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tint) - b1*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b7*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b2-b7)*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tint) - b2*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b8*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b3-b8)*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tint) - b3*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b9*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b4-b9)*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tint) - b4*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b10*WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b5-b10)*WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tint) - b5*WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmin)))
	return result

def TintOpt_objFun_W(tint, cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax):
	#input: Tint (intermediate temperature, in kiloKelvin); Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA(cp0,cpInf,B,a0,a1,a2,a3,tmin,tmax,tint, 1)
	result = (Wilhoit2IntM1(cp0,cpInf,B,a0,a1,a2,a3,tmax) - Wilhoit2IntM1(cp0,cpInf,B,a0,a1,a2,a3,tmin) +
				 NASA2IntM1(b1,b2,b3,b4,b5,tint)-NASA2IntM1(b1,b2,b3,b4,b5,tmin) + NASA2IntM1(b6,b7,b8,b9,b10,tmax) - NASA2IntM1(b6,b7,b8,b9,b10,tint)
				 - 2* (b6*WilhoitIntM1(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b1-b6)*WilhoitIntM1(cp0,cpInf,B,a0,a1,a2,a3,tint) - b1*WilhoitIntM1(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b7*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b2-b7)*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tint) - b2*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b8*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b3-b8)*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tint) - b3*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b9*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b4-b9)*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tint) - b4*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmin)
				 +b10*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b5-b10)*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tint) - b5*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmin)))
	return result

#analytical integrals:

#input (for all functions WilhoitXIntN): Wilhoit parameters: Cp0, CpInf, B, a0, a1, a2, a3 and t; units of Cp0/CpInf and B/t should be consistent and determine the units of the output 

### These functions belong inside the ThermoWilhoitData class:

def WilhoitInt0Orig(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R, t'] evaluated at t'=t
	result = ( cpInf*t + (a3*B**6*(cp0 - cpInf))/(5.*(B + t)**5) - ((a2 + 5*a3)*B**5*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 4*a2 + 10*a3)*B**4*(cp0 - cpInf))/(3.*(B + t)**3) - 
		((a0 + 3*a1 + 6*a2 + 10*a3)*B**3*(cp0 - cpInf))/(2.*(B + t)**2) + ((1 + 2*a0 + 3*a1 + 4*a2 + 5*a3)*B**2*(cp0 - cpInf))/(B + t) + (2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*math.log(B + t))
	return result

#a faster version of the integral based on H from Yelvington's thesis; it differs from the original (see above) by a constant (dependent on parameters but independent of t)
def WilhoitInt0(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R, t'] evaluated at t'=t
	y = t/(t+B)
	result = cp0*t - (-cp0 + cpInf)*t*(y**2*((3*a0 + a1 + a2 + a3)/6. + ((4*a1 + a2 + a3)*y)/12. + ((5*a2 + a3)*y**2)/20. + (a3*y**3)/5.) + (2 + a0 + a1 + a2 + a3)*(-1 + y/2. + (-1 + 1/y)*math.log(B + t)))
	return result

def WilhoitInt1(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R*t, t'] evaluated at t'=t
	result = ( (2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t + (cpInf*t**2)/2. + (a3*B**7*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 6*a3)*B**6*(cp0 - cpInf))/(4.*(B + t)**4) - 
		((a1 + 5*(a2 + 3*a3))*B**5*(cp0 - cpInf))/(3.*(B + t)**3) + ((a0 + 4*a1 + 10*(a2 + 2*a3))*B**4*(cp0 - cpInf))/(2.*(B + t)**2) - 
		((1 + 3*a0 + 6*a1 + 10*a2 + 15*a3)*B**3*(cp0 - cpInf))/(B + t) - (3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*math.log(B + t))
	return result

def WilhoitInt2(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R*t^2, t'] evaluated at t'=t
	result = ( -((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*t) + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**2)/2. + (cpInf*t**3)/3. + (a3*B**8*(cp0 - cpInf))/(5.*(B + t)**5) - 
		((a2 + 7*a3)*B**7*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 6*a2 + 21*a3)*B**6*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 5*(a1 + 3*a2 + 7*a3))*B**5*(cp0 - cpInf))/(2.*(B + t)**2) + 
		((1 + 4*a0 + 10*a1 + 20*a2 + 35*a3)*B**4*(cp0 - cpInf))/(B + t) + (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*math.log(B + t))
	return result

def WilhoitInt3(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R*t^3, t'] evaluated at t'=t
	result = ( (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t + ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**2)/2. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**3)/3. + 
		(cpInf*t**4)/4. + (a3*B**9*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 8*a3)*B**8*(cp0 - cpInf))/(4.*(B + t)**4) - ((a1 + 7*(a2 + 4*a3))*B**7*(cp0 - cpInf))/(3.*(B + t)**3) + 
		((a0 + 6*a1 + 21*a2 + 56*a3)*B**6*(cp0 - cpInf))/(2.*(B + t)**2) - ((1 + 5*a0 + 15*a1 + 35*a2 + 70*a3)*B**5*(cp0 - cpInf))/(B + t) - 
		(5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*math.log(B + t))
	return result

def WilhoitInt4(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R*t^4, t'] evaluated at t'=t
	result = ( -((5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*t) + ((4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t**2)/2. + 
		((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**3)/3. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**4)/4. + (cpInf*t**5)/5. + (a3*B**10*(cp0 - cpInf))/(5.*(B + t)**5) - 
		((a2 + 9*a3)*B**9*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 8*a2 + 36*a3)*B**8*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 7*(a1 + 4*(a2 + 3*a3)))*B**7*(cp0 - cpInf))/(2.*(B + t)**2) + 
		((1 + 6*a0 + 21*a1 + 56*a2 + 126*a3)*B**6*(cp0 - cpInf))/(B + t) + (6 + 15*a0 + 35*a1 + 70*a2 + 126*a3)*B**5*(cp0 - cpInf)*math.log(B + t))
	return result

def Wilhoit2Int(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[(Cp(Wilhoit)/R)^2, t'] evaluated at t'=t
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
		2*(2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*cpInf*math.log(B + t))
	return result

def NASA2Int(c1,c2,c3,c4,c5,t) :
	#input: NASA parameters for Cp/R, c1, c2, c3, c4, c5 (either low or high temp parameters), temperature t (in kiloKelvin; an endpoint of the low or high temp range
	#output: the quantity Integrate[(Cp(NASA)/R)^2, t'] evaluated at t'=t 
	#can speed further by precomputing and storing e.g. thigh^2, tlow^2, etc.
	result = c1*c1*t + c1*c2*t*t + (2*c1*c3+c2*c2)/3*t*t*t + (c1*c4+c2*c3)/2*t*t*t*t + (2*c1*c5 + 2*c2*c4 + c3*c3)/5*t*t*t*t*t + (c2*c5 + c3*c4)/3*t*t*t*t*t*t + (2*c3*c5 + c4*c4)/7*t*t*t*t*t*t*t + c4*c5/4*t*t*t*t*t*t*t*t + c5*c5/9*t*t*t*t*t*t*t*t*t
	return result

def WilhoitIntM1Orig(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R*t^-1, t'] evaluated at t'=t
	result = (a3*B**5*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 4*a3)*B**4*(cp0 - cpInf))/(4.*(B + t)**4) - ((a1 + 3*a2 + 6*a3)*B**3*(cp0 - cpInf))/(3.*(B + t)**3) + ((a0 + 2*a1 + 3*a2 + 4*a3)*B**2*(cp0 - cpInf))/(2.*(B + t)**2) - ((1 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf))/(B + t) + cp0*math.log(t) + (-cp0 + cpInf)*math.log(B + t) 
	return result

#a faster version of the integral based on S from Yelvington's thesis; it differs from the original by a constant (dependent on parameters but independent of t)
def WilhoitIntM1(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[Cp(Wilhoit)/R*t^-1, t'] evaluated at t'=t
	y = t/(t+B)
	result= cpInf*math.log(t)-(cpInf-cp0)*(math.log(y)+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5))))) 
	return result

def Wilhoit2IntM1(cp0, cpInf, B, a0, a1, a2, a3, t):
	#output: the quantity Integrate[(Cp(Wilhoit)/R)^2*t^-1, t'] evaluated at t'=t
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
		(B*(cp0 - cpInf)*(cp0 - (3 + 2*a0 + 2*a1 + 2*a2 + 2*a3)*cpInf))/(B + t) + cp0**2*math.log(t) + (-cp0**2 + cpInf**2)*math.log(B + t))
	return result

def NASA2IntM1(c1,c2,c3,c4,c5,t) :
	#input: NASA parameters for Cp/R, c1, c2, c3, c4, c5 (either low or high temp parameters), temperature t (in kiloKelvin; an endpoint of the low or high temp range
	#output: the quantity Integrate[(Cp(NASA)/R)^2*t^-1, t'] evaluated at t'=t 
	#can speed further by precomputing and storing e.g. thigh^2, tlow^2, etc.
	result = c1*c1*math.log(t) + 2*c1*c2*t + (2*c1*c3+c2*c2)/2*t*t + 2*(c1*c4+c2*c3)/3*t*t*t + (2*c1*c5 + 2*c2*c4 + c3*c3)/4*t*t*t*t + 2*(c2*c5 + c3*c4)/5*t*t*t*t*t + (2*c3*c5 + c4*c4)/6*t*t*t*t*t*t + 2*c4*c5/7*t*t*t*t*t*t*t + c5*c5/8*t*t*t*t*t*t*t*t
	return result



################################################################################
class ThermoDatabase(data.Database):
	"""
	Represent an RMG thermodynamics database.
	"""
	
	def __init__(self):
		"""Call the generic `data.Database.__init__()` method.
		
		This in turn creates
			* self.dictionary = Dictionary()
			* self.library = Library()
			* self.tree = Tree()
		"""
		data.Database.__init__(self)
		

	def load(self, dictstr, treestr, libstr):
		"""
		Load a thermodynamics group additivity database. The database is stored
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
				raise data.InvalidDatabaseException('Kinetics library should be tuple at this point. Instead got %r'%data) 
			else: 
				index,item = item # break apart tuple, recover the 'index' - the beginning of the line in the library file.
				# Is't it dangerous having a local variable with the same name as a module?
				# what if we want to raise another data.InvalidDatabaseException() ?
				if not ( item.__class__ is str or item.__class__ is unicode) :
					raise data.InvalidDatabaseException('Kinetics library data format is unrecognized.')
				
				items = item.split()
				try:
					thermoData = []; comment = ''
					# First 12 entries are thermo data
					for i in range(0, 12):
						thermoData.append(float(items[i]))
					# Remaining entries are comment
					for i in range(12, len(items)):
						comment += items[i] + ' '
					
					thermoGAData = ThermoGAData()
					thermoGAData.fromDatabase(thermoData, comment)
					thermoGAData.index = index
					self.library[label] = thermoGAData
				except (ValueError, IndexError), e:
					# Split data into link string and comment string; store
					# as list of length 2
					link = items[0]
					comment = item[len(link)+1:].strip()
					self.library[label] = [link, comment]
		
		# Check for well-formedness
		if not self.isWellFormed():
			raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (dictstr))
		
		#self.library.removeLinks()

	def toXML(self):
		"""
		Return an XML representation of the thermo database.
		"""
		import xml.dom.minidom
		dom = xml.dom.minidom.parseString('<database type="thermodynamics"></database>')
		root = dom.documentElement
		
		data.Database.toXML(self, dom, root)
		
		return dom.toprettyxml()

	def getThermoData(self, structure, atom):
		"""
		Determine the group additivity thermodynamic data for the atom `atom`
		in the structure `structure`.
		"""
		
		node = self.descendTree(structure, atom, None)
		
		if node not in self.library:
			# No data present (e.g. bath gas)
			data = ThermoGAData()
		else:
			data = self.library[node]
		
		while data.__class__ != ThermoGAData and data is not None:
			if data[0].__class__ == str or data[0].__class__ == unicode:
				data = self.library[data[0]]

		# This code prints the hierarchy of the found node; useful for debugging
#		result = ''
#		while node is not None:
#			result = ' -> ' + node + result
#			node = self.tree.parent[node]
#		print result[4:]

		return data

	def contains(self, structure):
		"""
		Search the dictionary for the specified `structure`. If found, the label
		corresponding to the structure is returned. If not found, :data:`None`
		is returned.
		"""
		for label, struct in self.dictionary.iteritems():
			if structure.isIsomorphic(struct):
				return label
		return None

################################################################################

class ThermoDatabaseSet:
	"""
	A set of thermodynamics group additivity databases, consisting of a primary
	database of functional groups and a number of secondary databases to provide
	corrections for 1,5-interactions, gauche interactions, radicals, rings,
	and other functionality. There is also a primary library containing data for
	individual species.
	"""
	
	
	def __init__(self):
		self.groupDatabase = ThermoDatabase()
		self.int15Database = ThermoDatabase()
		self.gaucheDatabase = ThermoDatabase()
		self.otherDatabase = ThermoDatabase()
		self.radicalDatabase = ThermoDatabase()
		self.ringDatabase = ThermoDatabase()
		self.primaryDatabase = ThermoDatabase()
	
	def load(self, datapath):
		"""
		Load a set of thermodynamics group additivity databases from the general
		database specified at `datapath`.
		"""
		import logging
		
		logging.debug('\tThermodynamics databases:')
		
		self.groupDatabase.load(datapath + 'thermo/Group_Dictionary.txt', \
			datapath + 'thermo/Group_Tree.txt', \
			datapath + 'thermo/Group_Library.txt')
		logging.debug('\t\tFunctional groups')
		
		self.int15Database.load(datapath + 'thermo/15_Dictionary.txt', \
			datapath + 'thermo/15_Tree.txt', \
			datapath + 'thermo/15_Library.txt')
		logging.debug('\t\t1,5 interactions')
		
		self.gaucheDatabase.load(datapath + 'thermo/Gauche_Dictionary.txt', \
			datapath + 'thermo/Gauche_Tree.txt', \
			datapath + 'thermo/Gauche_Library.txt')
		logging.debug('\t\tGauche interactions')
		
		self.otherDatabase.load(datapath + 'thermo/Other_Dictionary.txt', \
			datapath + 'thermo/Other_Tree.txt', \
			datapath + 'thermo/Other_Library.txt')
		logging.debug('\t\tOther corrections')
		
		self.radicalDatabase.load(datapath + 'thermo/Radical_Dictionary.txt', \
			datapath + 'thermo/Radical_Tree.txt', \
			datapath + 'thermo/Radical_Library.txt')
		logging.debug('\t\tRadical corrections')
		
		self.ringDatabase.load(datapath + 'thermo/Ring_Dictionary.txt', \
			datapath + 'thermo/Ring_Tree.txt', \
			datapath + 'thermo/Ring_Library.txt')
		logging.debug('\t\tRing corrections')
		
		self.primaryDatabase.load(datapath + 'thermo/Primary_Dictionary.txt', \
			'', \
			datapath + 'thermo/Primary_Library.txt')
		logging.debug('\t\tPrimary thermo database')
	
	def saveXML(self, datapath):
		"""
		Save the loaded databases in the set to XML files.
		"""
		
		f = open(datapath + 'thermo/group.xml', 'w')
		f.write(self.groupDatabase.toXML())
		f.close()
		
		f = open(datapath + 'thermo/1,5-interactions.xml', 'w')
		f.write(self.int15Database.toXML())
		f.close()
		
		f = open(datapath + 'thermo/gauche-interactions.xml', 'w')
		f.write(self.gaucheDatabase.toXML())
		f.close()
		
		f = open(datapath + 'thermo/other.xml', 'w')
		f.write(self.otherDatabase.toXML())
		f.close()
		
		f = open(datapath + 'thermo/radicals.xml', 'w')
		f.write(self.radicalDatabase.toXML())
		f.close()
		
		f = open(datapath + 'thermo/ring.xml', 'w')
		f.write(self.ringDatabase.toXML())
		f.close()
		
		f = open(datapath + 'thermo/primary.xml', 'w')
		f.write(self.primaryDatabase.toXML())
		f.close()
	
	def getThermoData(self, struct):
		"""
		Determine the group additivity thermodynamic data for the structure
		`structure`.
		"""
		
		import chem
		import structure
		
		# First check to see if structure is in primary thermo library
		label = self.primaryDatabase.contains(struct)
		if label is not None:
			return self.primaryDatabase.library[label]
		
		thermoData = ThermoGAData()
		
		if struct.getRadicalCount() > 0:
			# Make a copy of the structure so we don't change the original
			saturatedStruct = struct.copy()
			# Saturate structure by replacing all radicals with bonds to
			# hydrogen atoms
			added = {}
			for atom in saturatedStruct.atoms():
				for i in range(0, atom.getFreeElectronCount()):
					H = chem.Atom('H', '0')
					bond = chem.Bond([atom, H], 'S')
					saturatedStruct.addAtom(H)
					saturatedStruct.addBond(bond)
					atom.decreaseFreeElectron()
					if atom not in added:
						added[atom] = []
					added[atom].append(bond)
			# Update the atom types of the saturated structure (not sure why
			# this is necessary, because saturating with H shouldn't be
			# changing atom types, but it doesn't hurt anything and is not
			# very expensive, so will do it anyway)
			saturatedStruct.simplifyAtomTypes()
			saturatedStruct.updateAtomTypes()
			# Get thermo estimate for saturated form of structure
			thermoData = self.getThermoData(saturatedStruct)
			# For each radical site, get radical correction
			# Only one radical site should be considered at a time; all others
			# should be saturated with hydrogen atoms
			for atom in added:
			
				# Remove the added hydrogen atoms and bond and restore the radical
				for bond in added[atom]:
					H = bond.atoms[1]
					saturatedStruct.removeBond(bond)
					saturatedStruct.removeAtom(H)
					atom.increaseFreeElectron()
				
				thermoData += self.radicalDatabase.getThermoData(saturatedStruct, {'*':atom})
				
				# Re-saturate
				for bond in added[atom]:
					H = bond.atoms[1]
					saturatedStruct.addAtom(H)
					saturatedStruct.addBond(bond)
					atom.decreaseFreeElectron()
			
			# Subtract the enthalpy of the added hydrogens
			thermoData_H = self.primaryDatabase.library['H']
			for bond in added[atom]:
				thermoData.H298 -= thermoData_H.H298
				#thermoData.S298 -= thermoData_H.S298
			
			# Correct the entropy for the symmetry number
		
		else:
			# Generate estimate of thermodynamics
			for atom in struct.atoms():
				# Iterate over heavy (non-hydrogen) atoms
				if atom.isNonHydrogen():
					# Get initial thermo estimate from main group database
					thermoData += self.groupDatabase.getThermoData(struct, {'*':atom})
					# Correct for gauche and 1,5- interactions
					thermoData += self.gaucheDatabase.getThermoData(struct, {'*':atom})
					thermoData += self.int15Database.getThermoData(struct, {'*':atom})
					thermoData += self.otherDatabase.getThermoData(struct, {'*':atom})
			
			# Do ring corrections separately because we only want to match
			# each ring one time; this doesn't work yet
			rings = struct.getSmallestSetOfSmallestRings()
			for ring in rings:
				
				# Make a temporary structure containing only the atoms in the ring
				ringStructure = structure.Structure()
				for atom in ring: ringStructure.addAtom(atom)
				for atom1 in ring:
					for atom2 in ring:
						if struct.hasBond(atom1, atom2):
							ringStructure.addBond(struct.getBond(atom1, atom2))
				
				# Get thermo correction for this ring
				thermoData += self.ringDatabase.getThermoData(ringStructure, {})
		
		return thermoData

thermoDatabase = None
forbiddenStructures = None

################################################################################

def getThermoData(struct):
	"""
	Get the thermodynamic data associated with `structure` by looking in the
	loaded thermodynamic database.
	"""
	import constants
	import math
		
	GAthermoData = thermoDatabase.getThermoData(struct)

	# Correct entropy for symmetry number
	struct.calculateSymmetryNumber()
	GAthermoData.S298 -= constants.R * math.log(struct.symmetryNumber)
	
	return GAthermoData  # return here because Wilhoit conversion not working yet
	
	# Convert to Wilhoit
	rotors = struct.calculateNumberOfRotors()
	atoms = len(struct.atoms())
	linear = struct.isLinear()
	WilhoitData = convertGAtoWilhoit(GAthermoData,atoms,rotors,linear)

	return WilhoitData

################################################################################

if __name__ == '__main__':

	pass
