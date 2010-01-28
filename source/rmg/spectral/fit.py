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

import math
import numpy

import _fit

################################################################################

hoFreqLowerBound = 180.0
hoFreqUpperBound = 4000.0

hrFreqLowerBound = 180.0
hrFreqUpperBound = 4000.0

hrBarrLowerBound = 10.0
hrBarrUpperBound = 10000.0

hrFreqLow = 100.0
hrFreqMid = 150.0
hrFreqHigh = 300.0

maxIter = 1000

################################################################################

def fitSpectralDataToHeatCapacity(Tlist, Cvlist, Nvib, Nrot):

	x0, bl, bu, ind = setupCaseNvibNrot(Nvib, Nrot)
	
	_fit.params.setparams(
		p_tlist = numpy.array(Tlist),
		p_cvlist = numpy.array(Cvlist),
		p_mcon = 0,
		p_mequa = len(Tlist),
		p_nvars = len(x0),
		p_nvib = Nvib,
		p_nrot = Nrot,
		p_nu_low = hrFreqLow,
		p_nu_mid = hrFreqMid,
		p_nu_high = hrFreqHigh,
		)

	x, igo = _fit.fitmodes(
		x0 = numpy.array(x0),
		bl = numpy.array(bl),
		bu = numpy.array(bu),
		ind = numpy.array(ind),
		maxiter = maxIter,
		)
		
	print 'x =', list(x)
	print 'igo =', igo

	vib, rot = postprocessCaseNvibNrot(Nvib, Nrot, x)
	
	return vib, rot

################################################################################

def setupCaseNvibNrot(Nvib, Nrot):
	"""
	Fit an arbitrary number of harmonic oscillator and hindered rotor modes to
	the heat capacity data.
	"""

	ind = numpy.zeros(6, numpy.float64)		# Indicates type of bounds to apply (1 = lower, 2 = upper, 3 = both, 4 = neither)
	lb = numpy.zeros(6, numpy.float64)		# Lower bounds
	ub = numpy.zeros(6, numpy.float64)		# Upper bounds
	x0 = numpy.zeros(6, numpy.float64)		# Initial guess

	# x[0] corresponds to the degeneracy of the first harmonic oscillator
	ind[0] = 3
	lb[0]  = 0.0
	ub[0]  = float(Nvib)

	# x[1] corresponds to the first harmonic oscillator pseudo-frequency
	ind[1] = 3
	lb[1]  = hoFreqLowerBound
	ub[1]  = hoFreqUpperBound

	# x[2] corresponds to the second harmonic oscillator pseudo-frequency
	ind[2] = 3
	lb[2]  = hoFreqLowerBound
	ub[2]  = hoFreqUpperBound

	# x[3] corresponds to the degeneracy of the first hindered rotor
	ind[3] = 3
	lb[3]  = 0.0
	ub[3]  = float(Nrot)

	# x[4] corresponds to the first hindered rotor pseudo-barrier
	ind[4] = 3
	lb[4]  = hrBarrLowerBound
	ub[4]  = hrBarrUpperBound

	# x[5] corresponds to the second hindered rotor pseudo-barrier
	ind[5] = 3
	lb[5]  = hrBarrLowerBound
	ub[5]  = hrBarrUpperBound

	# Initial guess
	x0[0] = float(math.floor(Nvib / 2.0))
	x0[1] = 300.0
	x0[2] = 3000.0
	x0[3] = float(math.floor(Nrot / 2.0))
	x0[4] = 500.0
	x0[5] = 4000.0

	return x0, lb, ub, ind

def postprocessCaseNvibNrot(Nvib, Nrot, x):

	Nvib1 = int(round(x[0]))
	Nvib2 = Nvib - Nvib1
	Nrot1 = int(round(x[3]))
	Nrot2 = Nrot - Nrot1

	vib = []
	if Nvib1 > 0: vib.append([x[1], Nvib1])
	if Nvib2 > 0: vib.append([x[2], Nvib2])
	rot = []
	if Nrot1 > 0: rot.append([hrFreqLow, x[4], Nrot1])
	if Nrot2 > 0: rot.append([hrFreqHigh, x[5], Nrot2])

	return vib, rot

################################################################################
