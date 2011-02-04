#!/usr/bin/python
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
Contains functions for converting between some of the thermodynamics models
given in the :mod:`chempy.thermo` module. The two primary functions are:

* :func:`convertGAtoWilhoit()` - converts a :class:`ThermoGAModel` to a :class:`WilhoitModel`

* :func:`convertWilhoitToNASA()` - converts a :class:`WilhoitModel` to a :class:`NASAModel`

"""

import math
import numpy
import logging
import cython
from scipy import zeros, linalg, optimize, integrate

import rmgpy.chem.constants as constants

from rmgpy.chem.thermo import ThermoGAModel, WilhoitModel, NASAPolynomial, NASAModel

################################################################################

def convertGAtoWilhoit(GAthermo, atoms, rotors, linear, B0=500.0, constantB=False):
    """
    Convert a :class:`ThermoGAModel` object `GAthermo` to a
    :class:`WilhoitModel` object. You must specify the number of `atoms`,
    internal `rotors` and the linearity `linear` of the molecule so that the
    proper limits of heat capacity at zero and infinite temperature can be
    determined. You can also specify an initial guess of the scaling temperature
    `B0` to use, and whether or not to allow that parameter to vary
    (`constantB`). Returns the fitted :class:`WilhoitModel` object.
    """
    freq = 3 * atoms - (5 if linear else 6) - rotors
    wilhoit = WilhoitModel()
    if constantB:
        wilhoit.fitToDataForConstantB(GAthermo.Tdata, GAthermo.Cpdata, linear, freq, rotors, GAthermo.H298, GAthermo.S298, B0)
    else:
        wilhoit.fitToData(GAthermo.Tdata, GAthermo.Cpdata, linear, freq, rotors, GAthermo.H298, GAthermo.S298, B0)
    return wilhoit

################################################################################

def convertWilhoitToNASA(wilhoit, Tmin, Tmax, Tint, fixedTint=False, weighting=True, continuity=3):
    """
    Convert a :class:`WilhoitModel` object `Wilhoit` to a :class:`NASAModel` 
    object. You must specify the minimum and maximum temperatures of the fit
    `Tmin` and `Tmax`, as well as the intermediate temperature `Tint` to use
    as the bridge between the two fitted polynomials. The remaining parameters
    can be used to modify the fitting algorithm used:
    
    * `fixedTint` - ``False`` to allow `Tint` to vary in order to improve the fit, or ``True`` to keep it fixed

    * `weighting` - ``True`` to weight the fit by :math:`T^{-1}` to emphasize good fit at lower temperatures, or ``False`` to not use weighting

    * `continuity` - The number of continuity constraints to enforce at `Tint`:

        - 0: no constraints on continuity of :math:`C_\\mathrm{p}(T)` at `Tint`

        - 1: constrain :math:`C_\\mathrm{p}(T)` to be continous at `Tint`

        - 2: constrain :math:`C_\\mathrm{p}(T)` and :math:`\\frac{d C_\\mathrm{p}}{dT}` to be continuous at `Tint`

        - 3: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, and :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}` to be continuous at `Tint`

        - 4: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, and :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}` to be continuous at `Tint`

        - 5: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}`, and :math:`\\frac{d^4 C_\\mathrm{p}}{dT^4}` to be continuous at `Tint`
        
    Note that values of `continuity` of 5 or higher effectively constrain all
    the coefficients to be equal and should be equivalent to fitting only one
    polynomial (rather than two).

    Returns the fitted :class:`NASAModel` object containing the two fitted
    :class:`NASAPolynomial` objects.
    """

    # Scale the temperatures to kK
    Tmin /= 1000.
    Tint /= 1000.
    Tmax /= 1000.
    
    # Make copy of Wilhoit data so we don't modify the original
    wilhoit_scaled = WilhoitModel(wilhoit.cp0, wilhoit.cpInf, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3, wilhoit.H0, wilhoit.S0, wilhoit.B, Tmin=wilhoit.Tmin, Tmax=wilhoit.Tmax, comment=wilhoit.comment)
    # Rescale Wilhoit parameters
    wilhoit_scaled.cp0 /= constants.R
    wilhoit_scaled.cpInf /= constants.R
    wilhoit_scaled.B /= 1000.

    #if we are using fixed Tint, do not allow Tint to float
    if fixedTint:
        nasa_low, nasa_high = Wilhoit2NASA(wilhoit_scaled, Tmin, Tmax, Tint, weighting, continuity)
    else:
        nasa_low, nasa_high, Tint = Wilhoit2NASA_TintOpt(wilhoit_scaled, Tmin, Tmax, weighting, continuity)
    iseUnw = TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, 0, continuity) #the scaled, unweighted ISE (integral of squared error)
    rmsUnw = math.sqrt(iseUnw/(Tmax-Tmin))
    rmsStr = '(Unweighted) RMS error = %.3f*R;'%(rmsUnw)
    if(weighting == 1):
        iseWei= TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, weighting, continuity) #the scaled, weighted ISE
        rmsWei = math.sqrt(iseWei/math.log(Tmax/Tmin))
        rmsStr = 'Weighted RMS error = %.3f*R;'%(rmsWei)+rmsStr

    #print a warning if the rms fit is worse that 0.25*R
    if(rmsUnw > 0.25 or rmsWei > 0.25):
        logging.warning("Poor Wilhoit-to-NASA fit quality: RMS error = %.3f*R" % (rmsWei if weighting == 1 else rmsUnw))

    #restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
    Tint *= 1000.
    Tmin *= 1000.
    Tmax *= 1000.

    nasa_low.c1 /= 1000.
    nasa_low.c2 /= 1000000.
    nasa_low.c3 /= 1000000000.
    nasa_low.c4 /= 1000000000000.

    nasa_high.c1 /= 1000.
    nasa_high.c2 /= 1000000.
    nasa_high.c3 /= 1000000000.
    nasa_high.c4 /= 1000000000000.

    # output comment
    comment = 'NASA function fitted to Wilhoit function. ' + rmsStr + wilhoit.comment
    nasa_low.Tmin = Tmin; nasa_low.Tmax = Tint
    nasa_low.comment = 'Low temperature range polynomial'
    nasa_high.Tmin = Tint; nasa_high.Tmax = Tmax
    nasa_high.comment = 'High temperature range polynomial'

    #for the low polynomial, we want the results to match the Wilhoit value at 298.15K
    #low polynomial enthalpy:
    Hlow = (wilhoit.getEnthalpy(298.15) - nasa_low.getEnthalpy(298.15))/constants.R
    #low polynomial entropy:
    Slow = (wilhoit.getEntropy(298.15) - nasa_low.getEntropy(298.15))/constants.R

    # update last two coefficients
    nasa_low.c5 = Hlow
    nasa_low.c6 = Slow

    #for the high polynomial, we want the results to match the low polynomial value at tint
    #high polynomial enthalpy:
    Hhigh = (nasa_low.getEnthalpy(Tint) - nasa_high.getEnthalpy(Tint))/constants.R
    #high polynomial entropy:
    Shigh = (nasa_low.getEntropy(Tint) - nasa_high.getEntropy(Tint))/constants.R

    # update last two coefficients
    #polynomial_high.coeffs = (b6,b7,b8,b9,b10,Hhigh,Shigh)
    nasa_high.c5 = Hhigh
    nasa_high.c6 = Shigh

    return NASAModel(Tmin=Tmin, Tmax=Tmax, polynomials=[nasa_low,nasa_high], comment=comment)

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
    A = zeros([10+contCons,10+contCons])
    b = zeros([10+contCons])

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
    w0int = Wilhoit_integral_T0(wilhoit, tint)
    w1int = Wilhoit_integral_T1(wilhoit, tint)
    w2int = Wilhoit_integral_T2(wilhoit, tint)
    w3int = Wilhoit_integral_T3(wilhoit, tint)
    w0min = Wilhoit_integral_T0(wilhoit, tmin)
    w1min = Wilhoit_integral_T1(wilhoit, tmin)
    w2min = Wilhoit_integral_T2(wilhoit, tmin)
    w3min = Wilhoit_integral_T3(wilhoit, tmin)
    w0max = Wilhoit_integral_T0(wilhoit, tmax)
    w1max = Wilhoit_integral_T1(wilhoit, tmax)
    w2max = Wilhoit_integral_T2(wilhoit, tmax)
    w3max = Wilhoit_integral_T3(wilhoit, tmax)
    if weighting:
        wM1int = Wilhoit_integral_TM1(wilhoit, tint)
        wM1min = Wilhoit_integral_TM1(wilhoit, tmin)
        wM1max = Wilhoit_integral_TM1(wilhoit, tmax)
    else:
        w4int = Wilhoit_integral_T4(wilhoit, tint)
        w4min = Wilhoit_integral_T4(wilhoit, tmin)
        w4max = Wilhoit_integral_T4(wilhoit, tmax)

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
    (nasa1, nasa2) = Wilhoit2NASA(wilhoit, tmin, tmax, tint, weighting, contCons)
    return nasa1, nasa2, tint

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

    q0=Wilhoit_integral_T0(wilhoit, tint)
    q1=Wilhoit_integral_T1(wilhoit, tint)
    q2=Wilhoit_integral_T2(wilhoit, tint)
    q3=Wilhoit_integral_T3(wilhoit, tint)
    q4=Wilhoit_integral_T4(wilhoit, tint)
    result = (Wilhoit_integral2_T0(wilhoit, tmax) - Wilhoit_integral2_T0(wilhoit, tmin) +
                 NASAPolynomial_integral2_T0(nasa_low, tint) - NASAPolynomial_integral2_T0(nasa_low, tmin) +
                 NASAPolynomial_integral2_T0(nasa_high, tmax) - NASAPolynomial_integral2_T0(nasa_high, tint)
                 - 2* (b6*(Wilhoit_integral_T0(wilhoit, tmax)-q0)+b1*(q0-Wilhoit_integral_T0(wilhoit, tmin))
                 +b7*(Wilhoit_integral_T1(wilhoit, tmax) - q1) +b2*(q1 - Wilhoit_integral_T1(wilhoit, tmin))
                 +b8*(Wilhoit_integral_T2(wilhoit, tmax) - q2) +b3*(q2 - Wilhoit_integral_T2(wilhoit, tmin))
                 +b9*(Wilhoit_integral_T3(wilhoit, tmax) - q3) +b4*(q3 - Wilhoit_integral_T3(wilhoit, tmin))
                 +b10*(Wilhoit_integral_T4(wilhoit, tmax) - q4)+b5*(q4 - Wilhoit_integral_T4(wilhoit, tmin))))

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

    qM1=Wilhoit_integral_TM1(wilhoit, tint)
    q0=Wilhoit_integral_T0(wilhoit, tint)
    q1=Wilhoit_integral_T1(wilhoit, tint)
    q2=Wilhoit_integral_T2(wilhoit, tint)
    q3=Wilhoit_integral_T3(wilhoit, tint)
    result = (Wilhoit_integral2_TM1(wilhoit, tmax) - Wilhoit_integral2_TM1(wilhoit, tmin) +
                 NASAPolynomial_integral2_TM1(nasa_low, tint) - NASAPolynomial_integral2_TM1(nasa_low, tmin) +
                 NASAPolynomial_integral2_TM1(nasa_high, tmax) - NASAPolynomial_integral2_TM1(nasa_high, tint)
                 - 2* (b6*(Wilhoit_integral_TM1(wilhoit, tmax)-qM1)+b1*(qM1 - Wilhoit_integral_TM1(wilhoit, tmin))
                 +b7*(Wilhoit_integral_T0(wilhoit, tmax)-q0)+b2*(q0 - Wilhoit_integral_T0(wilhoit, tmin))
                 +b8*(Wilhoit_integral_T1(wilhoit, tmax)-q1)+b3*(q1 - Wilhoit_integral_T1(wilhoit, tmin))
                 +b9*(Wilhoit_integral_T2(wilhoit, tmax)-q2)+b4*(q2 - Wilhoit_integral_T2(wilhoit, tmin))
                 +b10*(Wilhoit_integral_T3(wilhoit, tmax)-q3)+b5*(q3 - Wilhoit_integral_T3(wilhoit, tmin))))

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
    else:
        rmsWei = 0.0

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
    A = zeros([10+contCons,10+contCons])
    b = zeros([10+contCons])

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
    (nasa1, nasa2) = Cp2NASA(CpObject, tmin, tmax, tint, weighting, contCons)
    return nasa1, nasa2, tint

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

################################################################################

#a faster version of the integral based on H from Yelvington's thesis; it differs from the original (see above) by a constant (dependent on parameters but independent of t)
def Wilhoit_integral_T0(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(y=cython.double, y2=cython.double, logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    y = t/(t+B)
    y2 = y*y
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = cp0*t - (cpInf-cp0)*t*(y2*((3*a0 + a1 + a2 + a3)/6. + (4*a1 + a2 + a3)*y/12. + (5*a2 + a3)*y2/20. + a3*y2*y/5.) + (2 + a0 + a1 + a2 + a3)*( y/2. - 1 + (1/y-1)*logBplust))
    return result

#a faster version of the integral based on S from Yelvington's thesis; it differs from the original by a constant (dependent on parameters but independent of t)
def Wilhoit_integral_TM1(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^-1, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(y=cython.double, logt=cython.double, logy=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    y = t/(t+B)
    if cython.compiled:
        logy = log(y); logt = log(t)
    else:
        logy = math.log(y); logt = math.log(t)
    result = cpInf*logt-(cpInf-cp0)*(logy+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5)))))
    return result

def Wilhoit_integral_T1(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( (2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t + (cpInf*t**2)/2. + (a3*B**7*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 6*a3)*B**6*(cp0 - cpInf))/(4.*(B + t)**4) -
        ((a1 + 5*(a2 + 3*a3))*B**5*(cp0 - cpInf))/(3.*(B + t)**3) + ((a0 + 4*a1 + 10*(a2 + 2*a3))*B**4*(cp0 - cpInf))/(2.*(B + t)**2) -
        ((1 + 3*a0 + 6*a1 + 10*a2 + 15*a3)*B**3*(cp0 - cpInf))/(B + t) - (3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral_T2(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^2, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( -((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*t) + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**2)/2. + (cpInf*t**3)/3. + (a3*B**8*(cp0 - cpInf))/(5.*(B + t)**5) -
        ((a2 + 7*a3)*B**7*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 6*a2 + 21*a3)*B**6*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 5*(a1 + 3*a2 + 7*a3))*B**5*(cp0 - cpInf))/(2.*(B + t)**2) +
        ((1 + 4*a0 + 10*a1 + 20*a2 + 35*a3)*B**4*(cp0 - cpInf))/(B + t) + (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral_T3(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^3, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t + ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**2)/2. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**3)/3. +
        (cpInf*t**4)/4. + (a3*B**9*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 8*a3)*B**8*(cp0 - cpInf))/(4.*(B + t)**4) - ((a1 + 7*(a2 + 4*a3))*B**7*(cp0 - cpInf))/(3.*(B + t)**3) +
        ((a0 + 6*a1 + 21*a2 + 56*a3)*B**6*(cp0 - cpInf))/(2.*(B + t)**2) - ((1 + 5*a0 + 15*a1 + 35*a2 + 70*a3)*B**5*(cp0 - cpInf))/(B + t) -
        (5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral_T4(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^4, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( -((5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*t) + ((4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t**2)/2. +
        ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**3)/3. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**4)/4. + (cpInf*t**5)/5. + (a3*B**10*(cp0 - cpInf))/(5.*(B + t)**5) -
        ((a2 + 9*a3)*B**9*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 8*a2 + 36*a3)*B**8*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 7*(a1 + 4*(a2 + 3*a3)))*B**7*(cp0 - cpInf))/(2.*(B + t)**2) +
        ((1 + 6*a0 + 21*a1 + 56*a2 + 126*a3)*B**6*(cp0 - cpInf))/(B + t) + (6 + 15*a0 + 35*a1 + 70*a2 + 126*a3)*B**5*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral2_T0(wilhoit, t):
    #output: the quantity Integrate[(Cp(Wilhoit)/R)^2, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
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

def Wilhoit_integral2_TM1(wilhoit, t):
    #output: the quantity Integrate[(Cp(Wilhoit)/R)^2*t^-1, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, logt=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0, wilhoit.cpInf, wilhoit.B, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
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

def NASAPolynomial_integral2_T0(polynomial, T):
    #output: the quantity Integrate[(Cp(NASAPolynomial)/R)^2, t'] evaluated at t'=t
    cython.declare(c0=cython.double, c1=cython.double, c2=cython.double, c3=cython.double, c4=cython.double)
    cython.declare(T2=cython.double, T4=cython.double, T8=cython.double)
    c0, c1, c2, c3, c4 = polynomial.c0, polynomial.c1, polynomial.c2, polynomial.c3, polynomial.c4
    T2=T*T; T4=T2*T2; T8=T4*T4
    result = (
        c0*c0*T + c0*c1*T2 + 2./3.*c0*c2*T2*T + 0.5*c0*c3*T4 + 0.4*c0*c4*T4*T +
        c1*c1*T2*T/3. + 0.5*c1*c2*T4 + 0.4*c1*c3*T4*T + c1*c4*T4*T2/3. +
        0.2*c2*c2*T4*T + c2*c3*T4*T2/3. + 2./7.*c2*c4*T4*T2*T +
        c3*c3*T4*T2*T/7. + 0.25*c3*c4*T8 +
        c4*c4*T8*T/9.
    )
    return result

def NASAPolynomial_integral2_TM1(polynomial, T):
    #output: the quantity Integrate[(Cp(NASAPolynomial)/R)^2*t^-1, t'] evaluated at t'=t
    cython.declare(c0=cython.double, c1=cython.double, c2=cython.double, c3=cython.double, c4=cython.double)
    cython.declare(T2=cython.double, T4=cython.double, logT=cython.double)
    c0, c1, c2, c3, c4 = polynomial.c0, polynomial.c1, polynomial.c2, polynomial.c3, polynomial.c4
    T2=T*T; T4=T2*T2
    if cython.compiled:
        logT = log(T)
    else:
        logT = math.log(T)
    result = (
        c0*c0*logT + 2*c0*c1*T + c0*c2*T2 + 2./3.*c0*c3*T2*T + 0.5*c0*c4*T4 +
        0.5*c1*c1*T2 + 2./3.*c1*c2*T2*T + 0.5*c1*c3*T4 + 0.4*c1*c4*T4*T +
        0.25*c2*c2*T4 + 0.4*c2*c3*T4*T + c2*c4*T4*T2/3. +
        c3*c3*T4*T2/6. + 2./7.*c3*c4*T4*T2*T +
        c4*c4*T4*T4/8.
    )
    return result

################################################################################

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
