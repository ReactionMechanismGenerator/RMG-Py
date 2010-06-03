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
This module contains classes and functions for working with chemical reactions.

From the `IUPAC Compendium of Chemical Terminology 
<http://dx.doi.org/10.1351/goldbook>`_, a chemical reaction is "a process that 
results in the interconversion of chemical species".

In ChemPy, a chemical reaction is called a Reaction object and is represented in
memory as an instance of the :class:`Reaction` class.
"""

import cython
import math

import constants
from exception import ChemPyError

################################################################################

class Reaction:
	"""
	A chemical reaction.
	
	=============== =========================== ================================
	Attribute       Type                        Description
	=============== =========================== ================================
	`index`         :class:`int`                A unique nonnegative integer index
	`reactants`     :class:`list`               The reactant species (as :class:`Species` objects)
	`products`      :class:`list`               The product species (as :class:`Species` objects)
	`kinetics`      :class:`KineticsModel`      The kinetics model to use for the reaction
	=============== =========================== ================================
	
	"""
	
	def __init__(self, index=-1, reactants=None, products=None, kinetics=None):
		self.index = index
		self.reactants = reactants
		self.products = products
		self.kinetics = kinetics
	
	def __repr__(self):
		"""
		Return a string representation of the reaction, suitable for console output.
		"""
		return "<Reaction %i '%s'>" % (self.index, str(self))
	
	def __str__(self):
		"""
		Return a string representation of the reaction, in the form 'A + B <=> C + D'.
		"""
		return ' <=> '.join([' + '.join([str(s) for s in self.reactants]), ' + '.join([str(s) for s in self.products])])

	def getEnthalpyOfReaction(self, T):
		"""
		Return the enthalpy of reaction in J/mol evaluated at temperature `T`
		in K.
		"""
		cython.declare(dHrxn=cython.double, reactant=Species, product=Species)
		dHrxn = -self.reactants[0].getEnthalpy(T)
		for reactant in self.reactants[1:]:
			dHrxn -= reactant.getEnthalpy(T)
		for product in self.products:
			dHrxn += product.getEnthalpy(T)
		return dHrxn

	def getEntropyOfReaction(self, T):
		"""
		Return the entropy of reaction in J/mol*K evaluated at temperature `T`
		in K.
		"""
		cython.declare(dSrxn=cython.double, reactant=Species, product=Species)
		dSrxn = -self.reactants[0].getEntropy(T)
		for reactant in self.reactants[1:]:
			dSrxn -= reactant.getEntropy(T)
		for product in self.products:
			dSrxn += product.getEntropy(T)
		return dSrxn

	def getFreeEnergyOfReaction(self, T):
		"""
		Return the Gibbs free energy of reaction in J/mol evaluated at
		temperature `T` in K.
		"""
		cython.declare(dGrxn=cython.double, reactant=Species, product=Species)
		dGrxn = -self.reactants[0].getFreeEnergy(T)
		for reactant in self.reactants[1:]:
			dGrxn -= reactant.getFreeEnergy(T)
		for product in self.products:
			dGrxn += product.getFreeEnergy(T)
		return dGrxn

	def getEquilibriumConstant(self, T, type='Kc'):
		"""
		Return the equilibrium constant for the reaction at the specified
		temperature `T` in K. The `type` parameter lets	you specify the
		quantities used in the equilibrium constant: ``Ka`` for	activities,
		``Kc`` for concentrations (default), or ``Kp`` for pressures. Note that
		this function currently assumes an ideal gas mixture.
		"""
		cython.declare(dGrxn=cython.double, K=cython.double, C0=cython.double, P0=cython.double)
		# Use free energy of reaction to calculate Ka
		dGrxn = self.getFreeEnergyOfReaction(T)
		K = math.exp(-dGrxn / constants.R / T)
		# Convert Ka to Kc or Kp if specified
		if type == 'Kc':
			# Convert from Ka to Kc; C0 is the reference concentration
			C0 = 1e5 / constants.R / T
			K *= C0 ** (len(self.products) - len(self.reactants))
		elif type == 'Kp':
			# Convert from Ka to Kp; P0 is the reference pressure
			P0 = 1e5
			K *= P0 ** (len(self.products) - len(self.reactants))
		elif type != 'Ka' or type != '':
			raise ChemPyError('Invalid type "%s" passed to Reaction.getEquilibriumConstant(); should be "Ka", "Kc", or "Kp".')
		return K

	def getRateConstant(self, T, P):
		"""
		Return the value of the rate constant at the temperature `T` in K and
		pressure `P` in Pa. The units of the returned rate constant depend on
		the number of reactants, but are some combination of moles, m^3, and s.
		"""
		return self.kinetics.getRateConstant(T, P)

	def getStoichiometricCoefficient(self, spec):
		"""
		Return the stoichiometric coefficient of species `spec` in the reaction.
		The stoichiometric coefficient is increased by one for each time `spec`
		appears as a product and decreased by one for each time `spec` appears
		as a reactant.
		"""
		cython.declare(stoich=cython.int, reactant=Species, product=Species)
		stoich = 0
		for reactant in self.reactants:
			if reactant is spec: stoich -= 1
		for product in self.products:
			if product is spec: stoich += 1
		return stoich

	def getRate(self, T, P, conc):
		"""
		Return the net rate of reaction in mol/m^3*s at temperature `T` in K
		and pressure `P` in Pa. The parameter `conc` is a dict with species as
		keys and concentrations as values. A reactant not found in the `conc`
		map is treated as having zero concentration.
		"""

		cython.declare(rateConstant=cython.double, equilibriumConstant=cython.double)
		cython.declare(forward=cython.double, reverse=cython.double)
		cython.declare(reactant=Species, product=Species)

		# Evaluate rate constant
		rateConstant = self.getRateConstant(T, P)
		if self.thirdBody: rateConstant *= sum(conc.values())

		# Evaluate equilibrium constant
		equilibriumConstant = self.getEquilibriumConstant(T)

		# Evaluate forward concentration product
		forward = 1.0
		for reactant in self.reactants:
			if reactant in conc:
				forward = forward * conc[reactant]
			else:
				forward = 0.0
				break

		# Evaluate reverse concentration product
		reverse = 1.0
		for product in self.products:
			if product in conc:
				reverse = reverse * conc[product]
			else:
				reverse = 0.0
				break

		# Return rate
		return rateConstant * (forward - reverse / equilibriumConstant)

