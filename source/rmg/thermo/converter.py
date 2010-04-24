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

import math
import scipy
from scipy import linalg, optimize, integrate

import rmg.constants
import rmg.log as logging

from model import *

################################################################################

def convertGAtoWilhoit(GAthermo, atoms, rotors, linear, fixedB=1, Bmin=300.0, Bmax=6000.0):
	"""Convert a Group Additivity thermo instance into a Wilhoit thermo instance.
	
	Takes a `ThermoGAModel` instance of themochemical data, and some extra information
	about the molecule used to calculate high- and low-temperature limits of Cp.
	These are the number of atoms (integer `atoms`), the number of rotors (integer `rotors`)
	and whether the molecule is linear (boolean `linear`)
	Returns a `WilhoitModel` instance.
	
	cf. Paul Yelvington's thesis, p. 185-186
	"""
	
	# get info from incoming group additivity thermo
	H298 = GAthermo.H298
	S298 = GAthermo.S298
	Cp_list = GAthermo.Cp
	T_list = ThermoGAModel.CpTlist  # usually [300, 400, 500, 600, 800, 1000, 1500] but why assume?
	R = constants.R
	B = WilhoitModelB # Constant (if fixed=1), set once in the class def.
	
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
	WilhoitThermo = WilhoitModel( cp0, cpInf, a0, a1, a2, a3, H0, S0, B=B, comment=comment)
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

def convertWilhoitToNASA(Wilhoit, fixed=1, weighting=1, tint=1000.0, Tmin = 298.0, Tmax=6000.0, contCons=3):
	"""Convert a Wilhoit thermo instance into a NASA polynomial thermo instance.
	
	Takes: a `WilhoitModel` instance of themochemical data.
		fixed: 1 (default) to fix tint; 0 to allow it to float to get a better fit
		weighting: 0 to not weight the fit by 1/T; 1 (default) to weight by 1/T to emphasize good fit at lower temperatures
		tint, Tmin, Tmax: intermediate, minimum, and maximum temperatures in Kelvin
		contCons: a measure of the continutity constraints on the fitted NASA polynomials; possible values are:
			    5: constrain Cp, dCp/dT, d2Cp/dT2, d3Cp/dT3, and d4Cp/dT4 to be continuous at tint; note: this effectively constrains all the coefficients to be equal and should be equivalent to fitting only one polynomial (rather than two)
			    4: constrain Cp, dCp/dT, d2Cp/dT2, and d3Cp/dT3 to be continuous at tint
			    3 (default): constrain Cp, dCp/dT, and d2Cp/dT2 to be continuous at tint
			    2: constrain Cp and dCp/dT to be continuous at tint
			    1: constrain Cp to be continous at tint
			    0: no constraints on continuity of Cp(T) at tint
			    note: 5th (and higher) derivatives of NASA Cp(T) are zero and hence will automatically be continuous at tint by the form of the Cp(T) function
	Returns a `NASAModel` instance containing two `NASAPolynomial`
	polynomials
	"""
	
	# Scale the temperatures to kK
	Tmin = Tmin/1000
	tint = tint/1000
	Tmax = Tmax/1000

	# Make copy of Wilhoit data so we don't modify the original
	wilhoit_scaled = WilhoitModel(Wilhoit.cp0, Wilhoit.cpInf, Wilhoit.a0, Wilhoit.a1, Wilhoit.a2, Wilhoit.a3, Wilhoit.H0, Wilhoit.S0, Wilhoit.comment, B=Wilhoit.B)
	# Rescale Wilhoit parameters
	wilhoit_scaled.cp0 /= constants.R
	wilhoit_scaled.cpInf /= constants.R
	wilhoit_scaled.B /= 1000.
	
	#if we are using fixed tint, do not allow tint to float
	if(fixed == 1):
		nasa_low, nasa_high = Wilhoit2NASA(wilhoit_scaled, Tmin, Tmax, tint, weighting, contCons)
	else:
		nasa_low, nasa_high, tint = Wilhoit2NASA_TintOpt(wilhoit_scaled, Tmin, Tmax, weighting, contCons)
	iseUnw = TintOpt_objFun(tint, wilhoit_scaled, Tmin, Tmax, 0, contCons) #the scaled, unweighted ISE (integral of squared error)
	rmsUnw = math.sqrt(iseUnw/(Tmax-Tmin))
	rmsStr = '(Unweighted) RMS error = %.3f*R;'%(rmsUnw)
	if(weighting == 1):
		iseWei= TintOpt_objFun(tint, wilhoit_scaled, Tmin, Tmax, weighting, contCons) #the scaled, weighted ISE
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
	nasa_low.Tmin = Tmin; nasa_low.Tmax = tint
	nasa_low.comment = 'Low temperature range polynomial'
	nasa_high.Tmin = tint; nasa_high.Tmax = Tmax
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
	
	NASAthermo = NASAModel(Tmin=Tmin, Tmax=Tmax, polynomials=[nasa_low,nasa_high], comment=comment)
	return NASAthermo

def Wilhoit2NASA(wilhoit, tmin, tmax, tint, weighting, contCons):
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

	nasa_low = NASAPolynomial(Tmin=0, Tmax=0, coeffs=[x[0], x[1], x[2], x[3], x[4], 0.0, 0.0], comment='')
	nasa_high = NASAPolynomial(Tmin=0, Tmax=0, coeffs=[x[5], x[6], x[7], x[8], x[9], 0.0, 0.0], comment='')

	return nasa_low, nasa_high
	
def Wilhoit2NASA_TintOpt(wilhoit, tmin, tmax, weighting, contCons):
	#input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
	#1. vary Tint, bounded by tmin and tmax, to minimize TintOpt_objFun
	#cf. http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html and http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html#scipy.optimize.fminbound)
	tint = optimize.fminbound(TintOpt_objFun, tmin, tmax, args=(wilhoit, tmin, tmax, weighting, contCons))
	#note that we have not used any guess when using this minimization routine
	#2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
	(nasa1, nasa2) = Wilhoit2NASA(wilhoit, tmin, tmax, tint[0] ,weighting, contCons)
	return nasa1, nasa2, tint[0]

def TintOpt_objFun(tint, wilhoit, tmin, tmax, weighting, contCons):
	#input: Tint (intermediate temperature, in kiloKelvin); Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	if (weighting == 1):
		result = TintOpt_objFun_W(tint, wilhoit, tmin, tmax, contCons)
	else:
		result = TintOpt_objFun_NW(tint, wilhoit, tmin, tmax, contCons)

	# numerical errors could accumulate to give a slightly negative result
	# this is unphysical (it's the integral of a *squared* error) so we
	# set it to zero to avoid later problems when we try find the square root.
	if result < 0:
		if result<-1E-13:
			logging.error("Greg thought he fixed the numerical problem, but apparently it is still an issue; please e-mail him with the following results:")
			logging.error(tint)
			logging.error(wilhoit)
			logging.error(tmin)
			logging.error(tmax)
			logging.error(weighting)
			logging.error(result)
		logging.info("Negative ISE of %f reset to zero."%(result))
		result = 0

	return result

def TintOpt_objFun_NW(tint, wilhoit, tmin, tmax, contCons):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.
	
	input: Tint (intermediate temperature, in kiloKelvin)
			Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, 
			Tmin (minimum temperature (in kiloKelvin), 
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Wilhoit2NASA(wilhoit,tmin,tmax,tint, 0, contCons)
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

def TintOpt_objFun_W(tint, wilhoit, tmin, tmax, contCons):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.
	
	If fit is close to perfect, result may be slightly negative due to numerical errors in evaluating this integral.
	input: Tint (intermediate temperature, in kiloKelvin)
			Wilhoit parameters: Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, 
			Tmin (minimum temperature (in kiloKelvin), 
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[1/t*(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Wilhoit2NASA(wilhoit,tmin,tmax,tint, 1, contCons)
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
def convertCpToNASA(CpObject, H298, S298, fixed=1, weighting=0, tint=1000.0, Tmin = 298.0, Tmax=6000.0, contCons=3):
	"""Convert an arbitrary heat capacity function into a NASA polynomial thermo instance (using numerical integration)

	Takes:  CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
		H298: enthalpy at 298.15 K (in J/mol)
		S298: entropy at 298.15 K (in J/mol-K)
		fixed: 1 (default) to fix tint; 0 to allow it to float to get a better fit
		weighting: 0 (default) to not weight the fit by 1/T; 1 to weight by 1/T to emphasize good fit at lower temperatures
		tint, Tmin, Tmax: intermediate, minimum, and maximum temperatures in Kelvin
		contCons: a measure of the continutity constraints on the fitted NASA polynomials; possible values are:
			    5: constrain Cp, dCp/dT, d2Cp/dT2, d3Cp/dT3, and d4Cp/dT4 to be continuous at tint; note: this effectively constrains all the coefficients to be equal and should be equivalent to fitting only one polynomial (rather than two)
			    4: constrain Cp, dCp/dT, d2Cp/dT2, and d3Cp/dT3 to be continuous at tint
			    3 (default): constrain Cp, dCp/dT, and d2Cp/dT2 to be continuous at tint
			    2: constrain Cp and dCp/dT to be continuous at tint
			    1: constrain Cp to be continous at tint
			    0: no constraints on continuity of Cp(T) at tint
			    note: 5th (and higher) derivatives of NASA Cp(T) are zero and hence will automatically be continuous at tint by the form of the Cp(T) function
	Returns a `NASAModel` instance containing two `NASAPolynomial` polynomials
	"""

	# Scale the temperatures to kK
	Tmin = Tmin/1000
	tint = tint/1000
	Tmax = Tmax/1000

	#if we are using fixed tint, do not allow tint to float
	if(fixed == 1):
		nasa_low, nasa_high = Cp2NASA(CpObject, Tmin, Tmax, tint, weighting, contCons)
	else:
		nasa_low, nasa_high, tint = Cp2NASA_TintOpt(CpObject, Tmin, Tmax, weighting, contCons)
	iseUnw = Cp_TintOpt_objFun(tint, CpObject, Tmin, Tmax, 0, contCons) #the scaled, unweighted ISE (integral of squared error)
	rmsUnw = math.sqrt(iseUnw/(Tmax-Tmin))
	rmsStr = '(Unweighted) RMS error = %.3f*R;'%(rmsUnw)
	if(weighting == 1):
		iseWei= Cp_TintOpt_objFun(tint, CpObject, Tmin, Tmax, weighting, contCons) #the scaled, weighted ISE
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
	nasa_low.Tmin = Tmin; nasa_low.Tmax = tint
	nasa_low.comment = 'Low temperature range polynomial'
	nasa_high.Tmin = tint; nasa_high.Tmax = Tmax
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

	NASAthermo = NASAModel(Tmin=Tmin, Tmax=Tmax, polynomials=[nasa_low,nasa_high], comment=comment)
	return NASAthermo

def Cp2NASA(CpObject, tmin, tmax, tint, weighting, contCons):
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

	nasa_low = NASAPolynomial(Tmin=0, Tmax=0, coeffs=[x[0], x[1], x[2], x[3], x[4], 0.0, 0.0], comment='')
	nasa_high = NASAPolynomial(Tmin=0, Tmax=0, coeffs=[x[5], x[6], x[7], x[8], x[9], 0.0, 0.0], comment='')

	return nasa_low, nasa_high

def Cp2NASA_TintOpt(CpObject, tmin, tmax, weighting, contCons):
	#input: CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
	#output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
	#1. vary Tint, bounded by tmin and tmax, to minimize TintOpt_objFun
	#cf. http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html and http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html#scipy.optimize.fminbound)
	tint = optimize.fminbound(Cp_TintOpt_objFun, tmin, tmax, args=(CpObject, tmin, tmax, weighting, contCons))
	#note that we have not used any guess when using this minimization routine
	#2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
	(nasa1, nasa2) = Cp2NASA(CpObject, tmin, tmax, tint[0] ,weighting, contCons)
	return nasa1, nasa2, tint[0]

def Cp_TintOpt_objFun(tint, CpObject, tmin, tmax, weighting, contCons):
	#input: Tint (intermediate temperature, in kiloKelvin); CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	if (weighting == 1):
		result = Cp_TintOpt_objFun_W(tint, CpObject, tmin, tmax, contCons)
	else:
		result = Cp_TintOpt_objFun_NW(tint, CpObject, tmin, tmax, contCons)

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

def Cp_TintOpt_objFun_NW(tint, CpObject, tmin, tmax, contCons):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.

	input: Tint (intermediate temperature, in kiloKelvin)
			CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
			Tmin (minimum temperature (in kiloKelvin),
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[(Cp/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Cp2NASA(CpObject,tmin,tmax,tint, 0, contCons)
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

def Cp_TintOpt_objFun_W(tint, CpObject, tmin, tmax, contCons):
	"""
	Evaluate the objective function - the integral of the square of the error in the fit.

	If fit is close to perfect, result may be slightly negative due to numerical errors in evaluating this integral.
	input: Tint (intermediate temperature, in kiloKelvin)
			CpObject: an object with method "getHeatCapacity(self,T) that will return Cp in J/mol-K with argument T in K
			Tmin (minimum temperature (in kiloKelvin),
			Tmax (maximum temperature (in kiloKelvin)
	output: the quantity Integrate[1/t*(Cp/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	"""
	nasa_low, nasa_high = Cp2NASA(CpObject,tmin,tmax,tint, 1, contCons)
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
