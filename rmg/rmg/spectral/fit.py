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
Contains functions for fitting of molecular degrees of freedom from
macroscopic properties, particularly the heat capacity. The primary function
of interest is :meth:`fitSpectralDataToHeatCapacity`, which does exactly this.
This function utilizes a number of helper functions that are called based on
the numbers and types of modes that should be fitted. The helper functions
are divided into those that setup the initial guess and bounds (named
``setupCase?vib?rot``) and those that convert the fitted results into
mode objects (named ``postprocessCase?vib?rot``).
"""

import math
import numpy

import rmg.log as logging

import _fit

################################################################################

# This block contains a number of global constants that control various
# aspects of the optimization (bounds, iterations, etc.)

# The lower bound for harmonic oscillator frequencies in cm^-1
hoFreqLowerBound = 180.0
# The upper bound for harmonic oscillator frequencies in cm^-1
hoFreqUpperBound = 4000.0

# The lower bound for hindered rotor frequencies in cm^-1
hrFreqLowerBound = 180.0
# The upper bound for hindered rotor frequencies in cm^-1
hrFreqUpperBound = 4000.0

# The lower bound for hindered rotor barrier heights in cm^-1
hrBarrLowerBound = 10.0
# The upper bound for hindered rotor barrier heights in cm^-1
hrBarrUpperBound = 10000.0

# The maximum number of iterations for the optimization solver to use
maxIter = 1000

################################################################################

class SpectralFitException(Exception):

	"""
	An exception used when attempting to fit molecular degrees of freedom to
	heat capacity data. Use the `msg` attribute to describe the cause of the
	exception.
	"""

	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return self.msg

################################################################################

def fitSpectralDataToHeatCapacity(struct, Tlist, Cvlist, Nvib, Nrot):
	"""
	For a given set of dimensionless heat capacity data `Cvlist` corresponding
	to temperature list `Tlist` in K, fit `Nvib` harmonic oscillator and `Nrot`
	hindered internal rotor modes. The fitting is done by calling into Fortran
	in order to use the DQED nonlinear constrained optimization code. This
	function returns two lists: the first contains lists of frequency-degeneracy
	pairs representing fitted harmonic oscillator modes, and the second contains
	lists of frequency-barrier-degeneracy triples representing fitted hindered
	rotor modes.
	"""

	# You must specify at least 7 heat capacity points to use in the fitting;
	# you can specify as many as you like above that minimum
	if len(Tlist) < 7:
		raise Exception('Unable to fit spectral data; you need to specify at least 7 heat capacity points.')

	# The number of optimization variables available is constrained to be less
	# than the number of heat capacity points
	# This is also capped to a (somewhat arbitrarily chosen) maximum of 16
	maxVariables = len(Tlist) - 1
	if maxVariables > 16:
		maxVariables = 16

	# Setup the initial guess and the bounds for the solver variables
	# The format of the variables depends on the numbers of oscillators and 
	# rotors being fitted:
	#	-	For low values of Nvib and Nrot, we can fit the individual  
	#		parameters directly
	#	-	For high values of Nvib and/or Nrot we are limited by the number of
	#		temperatures we are fitting at, and so we can only fit 
	#		pseudo-oscillators and/or pseudo-rotors
	if Nvib <= 0 and Nrot <= 0:
		return [], []
	elif Nvib + 2 * Nrot <= maxVariables:
		x0, bl, bu, ind = setupCaseDirect(Nvib, Nrot)
	elif Nvib + 2 <= maxVariables:
		x0, bl, bu, ind = setupCasePseudoRot(Nvib, Nrot)
	else:
		x0, bl, bu, ind = setupCasePseudo(Nvib, Nrot)

	# Set parameters that are not needed by the solver but are needed to 
	# evaluate the objective function and its Jacobian
	# These are stored in a Fortran 90 module called params
	_fit.setparams(
		p_tlist = numpy.array(Tlist),
		p_cvlist = numpy.array(Cvlist),
		p_mcon = 0,
		p_mequa = len(Tlist),
		p_nvars = len(x0),
		p_nvib = Nvib,
		p_nrot = Nrot
		)

	# Execute the optimization, passing the initial guess and bounds and other
	# solver options
	x, igo = _fit.fitmodes(
		x0 = numpy.array(x0),
		bl = numpy.array(bl),
		bu = numpy.array(bu),
		ind = numpy.array(ind),
		maxiter = maxIter,
		)
	
	# Clean up the temporary variables stored via _fit.setparams() earlier
	_fit.cleanup()

	if not numpy.isfinite(x).all():
		raise SpectralFitException('Returned solution vector is nonsensical: x = %s.' % (x))
	if igo == 8:
		logging.warning('Maximum number of iterations reached when fitting spectral data for %s.' % str(struct))
		#raise SpectralFitException('Maximum number of iterations reached; solution may be invalid.\nI was trying to fit %s oscillators and %s rotors.' % (Nvib, Nrot))
	
	# Convert the solution vector into sets of oscillators and rotors
	# The procedure for doing this depends on the content of the solution
	# vector, which itself depends on the number of oscillators and rotors
	# being fitted
	if Nvib + 2 * Nrot <= maxVariables:
		vib, rot = postprocessCaseDirect(Nvib, Nrot, x)
	elif Nvib + 2 <= maxVariables:
		vib, rot = postprocessCasePseudoRot(Nvib, Nrot, x)
	else:
		vib, rot = postprocessCasePseudo(Nvib, Nrot, x)
	
	return vib, rot

################################################################################

def setupCaseDirect(Nvib, Nrot):
	"""
	Fit a known number of harmonic oscillator and hindered rotor modes to
	the heat capacity data. This case is to be called when the total number of
	variables is less than the number of data points used in the fit. For this
	case, we can fit the oscillator frequencies and rotor frequency-barrier
	pairs directly.
	"""

	# Initialize the variables that are set by this function
	count = Nvib + 2 * Nrot
	ind = numpy.zeros(count, numpy.float64)		# Indicates type of bounds to apply (1 = lower, 2 = upper, 3 = both, 4 = neither)
	lb = numpy.zeros(count, numpy.float64)		# Lower bounds
	ub = numpy.zeros(count, numpy.float64)		# Upper bounds
	x0 = numpy.zeros(count, numpy.float64)		# Initial guess

	# The first Nvib variables correspond to real harmonic oscillator frequencies
	for i in range(Nvib):
		ind[i] = 3
		lb[i]  = hoFreqLowerBound
		ub[i]  = hoFreqUpperBound

	# The remaining 2 * Nrot variables correspond to real hindered rotor frequencies and barrier heights
	for i in range(Nrot):
		ind[Nvib+2*i] = 3
		lb[Nvib+2*i]  = hrFreqLowerBound
		ub[Nvib+2*i]  = hrFreqUpperBound
		ind[Nvib+2*i+1] = 3
		lb[Nvib+2*i+1]  = hrBarrLowerBound
		ub[Nvib+2*i+1]  = hrBarrUpperBound

	# Initial guesses within each mode type must be distinct or else the 
	# optimization will fail
	
	# Initial guess for harmonic oscillators
	if Nvib > 0:
		x0[0] = 200.0
		x0[1:Nvib] = numpy.linspace(800.0, 1600.0, Nvib-1)
	# Initial guess for hindered rotors
	if Nrot > 0:
		x0[Nvib] = 100.0
		x0[Nvib+1] = 100.0
		for i in range(1, Nrot):
			x0[Nvib+2*i] = x0[Nvib+2*i-2] + 20.0
			x0[Nvib+2*i+1] = x0[Nvib+2*i-1] + 100.0

	return x0, lb, ub, ind

def postprocessCaseDirect(Nvib, Nrot, x):
	"""
	Convert the optimization solution vector `x` to a set of harmonic
	oscillators and a set of hindered rotors. This case is to be called when
	the total number of variables is less than the number of data points used
	in the fit. For this case, we have fitted the oscillator frequencies and
	rotor frequency-barrier pairs directly.
	"""

	vib = []
	for i in range(Nvib):
		vib.append([x[i], 1])
	rot = []
	for i in range(Nrot):
		rot.append([x[Nvib+2*i], x[Nvib+2*i+1], 1])

	return vib, rot

################################################################################

def setupCasePseudo(Nvib, Nrot):
	"""
	Fit an arbitrary number of harmonic oscillator and hindered rotor modes to
	the heat capacity data. This case will fit two pseudo-oscillators and two
	pseudo-rotors such that their degeneracies sum to `Nvib` and `Nrot`,
	respectively.
	"""

	# Initialize the variables that are set by this function
	ind = numpy.zeros(6, numpy.float64)		# Indicates type of bounds to apply (1 = lower, 2 = upper, 3 = both, 4 = neither)
	lb = numpy.zeros(6, numpy.float64)		# Lower bounds
	ub = numpy.zeros(6, numpy.float64)		# Upper bounds
	x0 = numpy.zeros(6, numpy.float64)		# Initial guess

	# x[0] corresponds to the first harmonic oscillator (real) frequency
	ind[0] = 3
	lb[0]  = hoFreqLowerBound
	ub[0]  = hoFreqUpperBound
	
	# x[1] corresponds to the degeneracy of the second harmonic oscillator
	ind[1] = 3
	lb[1]  = 1.0
	ub[1]  = float(Nvib - 2)

	# x[2] corresponds to the second harmonic oscillator pseudo-frequency
	ind[2] = 3
	lb[2]  = hoFreqLowerBound
	ub[2]  = hoFreqUpperBound

	# x[3] corresponds to the third harmonic oscillator pseudo-frequency
	ind[3] = 3
	lb[3]  = hoFreqLowerBound
	ub[3]  = hoFreqUpperBound

	# x[4] corresponds to the hindered rotor pseudo-frequency
	ind[4] = 3
	lb[4]  = hrFreqLowerBound
	ub[4]  = hrFreqUpperBound

	# x[5] corresponds to the hindered rotor pseudo-barrier
	ind[5] = 3
	lb[5]  = hrBarrLowerBound
	ub[5]  = hrBarrUpperBound

	# Initial guess
	x0[0] = 300.0
	x0[1] = float(math.floor((Nvib - 1) / 2.0))
	x0[2] = 800.0
	x0[3] = 1600.0
	x0[4] = 100.0
	x0[5] = 300.0

	return x0, lb, ub, ind

def postprocessCasePseudo(Nvib, Nrot, x):
	"""
	Convert the optimization solution vector `x` to a set of harmonic
	oscillators and a set of hindered rotors for the case of an arbitrary
	number of oscillator and rotor modes `Nvib` and `Nrot`, respectively.
	"""

	Nvib2 = int(round(x[1]))
	Nvib3 = Nvib - Nvib2 - 1
	if Nvib2 < 0 or Nvib2 > Nvib-1 or Nvib3 < 0 or Nvib3 > Nvib-1:
		raise SpectralFitException('Invalid degeneracies %s and %s fitted for pseudo-frequencies.' % (Nvib2, Nvib3))

	vib = []
	vib.append([x[0], 1])
	if Nvib2 > 0: vib.append([x[2], Nvib2])
	if Nvib3 > 0: vib.append([x[3], Nvib3])
	rot = []
	if Nrot > 0: rot.append([x[4], x[5], Nrot])

	return vib, rot

################################################################################

def setupCasePseudoRot(Nvib, Nrot):
	"""
	Fit an arbitrary number of harmonic oscillator and hindered rotor modes to
	the heat capacity data. This case will fit `Nvib` real oscillators and two
	pseudo-rotors such that their degeneracies sum to `Nrot`.
	"""

	# Initialize the variables that are set by this function
	ind = numpy.zeros(Nvib+2, numpy.float64)		# Indicates type of bounds to apply (1 = lower, 2 = upper, 3 = both, 4 = neither)
	lb = numpy.zeros(Nvib+2, numpy.float64)		# Lower bounds
	ub = numpy.zeros(Nvib+2, numpy.float64)		# Upper bounds
	x0 = numpy.zeros(Nvib+2, numpy.float64)		# Initial guess

	# x[0] corresponds to the first harmonic oscillator (real) frequency
	for i in range(Nvib):
		ind[i] = 3
		lb[i]  = hoFreqLowerBound
		ub[i]  = hoFreqUpperBound

	# x[Nvib] corresponds to the hindered rotor pseudo-frequency
	ind[Nvib] = 3
	lb[Nvib]  = hrFreqLowerBound
	ub[Nvib]  = hrFreqUpperBound

	# x[Nvib+1] corresponds to the hindered rotor pseudo-barrier
	ind[Nvib+1] = 3
	lb[Nvib+1]  = hrBarrLowerBound
	ub[Nvib+1]  = hrBarrUpperBound

	# Initial guess for harmonic oscillators
	if Nvib > 0:
		x0[0] = 200.0
		x0[1:Nvib] = numpy.linspace(800.0, 1600.0, Nvib-1)
	# Initial guess for hindered pseudo-rotors
	x0[-2] = 100.0
	x0[-1] = 300.0

	return x0, lb, ub, ind

def postprocessCasePseudoRot(Nvib, Nrot, x):
	"""
	Convert the optimization solution vector `x` to a set of harmonic
	oscillators and a set of hindered rotors for the case of real oscillator
	and pseudo rotor modes `Nvib` and `Nrot`, respectively.
	"""

	vib = []
	for i in range(Nvib):
		vib.append([x[i], 1])
	rot = []
	if Nrot > 0: rot.append([x[-2], x[-1], Nrot])

	return vib, rot

################################################################################
