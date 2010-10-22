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

from chempy.thermo cimport ThermoGAModel, WilhoitModel, NASAPolynomial, NASAModel

cdef extern from "math.h":
    double log(double)


################################################################################

cpdef WilhoitModel convertGAtoWilhoit(ThermoGAModel GAthermo, int atoms, int rotors, bint linear, double B0=?, bint constantB=?)

cpdef NASAModel convertWilhoitToNASA(WilhoitModel wilhoit, double Tmin, double Tmax, double Tint, bint fixedTint=?, bint weighting=?, int continuity=?)

cpdef Wilhoit2NASA(WilhoitModel wilhoit, double tmin, double tmax, double tint, bint weighting, int contCons)

cpdef Wilhoit2NASA_TintOpt(WilhoitModel wilhoit, double tmin, double tmax, bint weighting, int contCons)

cpdef TintOpt_objFun(tint, WilhoitModel wilhoit, double tmin, double tmax, bint weighting, int contCons)

cpdef TintOpt_objFun_NW(double tint, WilhoitModel wilhoit, double tmin, double tmax, int contCons)

cpdef TintOpt_objFun_W(double tint, WilhoitModel wilhoit, double tmin, double tmax, int contCons)

cpdef convertCpToNASA(CpObject, double H298, double S298, int fixed=?, bint weighting=?, double tint=?, double Tmin=?, double Tmax=?, int contCons=?)

cpdef Cp2NASA(CpObject, double tmin, double tmax, double tint, bint weighting, int contCons)

cpdef Cp2NASA_TintOpt(CpObject, double tmin, double tmax, bint weighting, int contCons)

cpdef Cp_TintOpt_objFun(double tint, CpObject, double tmin, double tmax, bint weighting, int contCons)

cpdef Cp_TintOpt_objFun_NW(double tint, CpObject, double tmin, double tmax, int contCons)

cpdef Cp_TintOpt_objFun_W(double tint, CpObject, double tmin, double tmax, int contCons)

################################################################################

cpdef double Wilhoit_integral_T0(WilhoitModel wilhoit, double t)

cpdef double Wilhoit_integral_TM1(WilhoitModel wilhoit, double t)

cpdef double Wilhoit_integral_T1(WilhoitModel wilhoit, double t)

cpdef double Wilhoit_integral_T2(WilhoitModel wilhoit, double t)

cpdef double Wilhoit_integral_T3(WilhoitModel wilhoit, double t)

cpdef double Wilhoit_integral_T4(WilhoitModel wilhoit, double t)

cpdef double Wilhoit_integral2_T0(WilhoitModel wilhoit, double t)

cpdef double Wilhoit_integral2_TM1(WilhoitModel wilhoit, double t)

################################################################################

cpdef double NASAPolynomial_integral2_T0(NASAPolynomial polynomial, double T)

cpdef double NASAPolynomial_integral2_TM1(NASAPolynomial polynomial, double T)

################################################################################

cpdef Nintegral_T0(CpObject, double tmin, double tmax)

cpdef Nintegral_TM1(CpObject, double tmin, double tmax)

cpdef Nintegral_T1(CpObject, double tmin, double tmax)

cpdef Nintegral_T2(CpObject, double tmin, double tmax)

cpdef Nintegral_T3(CpObject, double tmin, double tmax)

cpdef Nintegral_T4(CpObject, double tmin, double tmax)

cpdef Nintegral2_T0(CpObject, double tmin, double tmax)

cpdef Nintegral2_TM1(CpObject, double tmin, double tmax)

cpdef Nintegral(CpObject, double tmin, double tmax, int n, int squared)

cpdef integrand(double t, CpObject, int n, int squared)
