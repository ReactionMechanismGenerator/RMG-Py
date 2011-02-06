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
import scipy.special

import rmgpy.chem.constants as constants
from pydqed import DQED

################################################################################

def harmonicOscillator_heatCapacity(T, freq):
    """
    Return the heat capacity in J/mol*K at the given set of temperatures `Tlist`
    in K for the harmonic oscillator with a frequency `freq` in cm^-1.
    """
    x = freq / (0.695039 * T)        # kB = 0.695039 cm^-1/K
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    return x * x * exp_x / one_minus_exp_x / one_minus_exp_x

def harmonicOscillator_d_heatCapacity_d_freq(T, freq):
    """
    Return the first derivative of the heat capacity with respect to the
    harmonic oscillator frequency in J/mol*K/cm^-1 at the given set of
    temperatures `Tlist` in K, evaluated at the frequency `freq` in cm^-1.
    """
    x = freq / (0.695039 * T)        # kB = 0.695039 cm^-1/K
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    return x * exp_x / one_minus_exp_x / one_minus_exp_x * (2.0 + x + 2.0 * x * exp_x / one_minus_exp_x) * x / freq

def hinderedRotor_heatCapacity(T, freq, barr):
    """
    Return the heat capacity in J/mol*K at the given set of temperatures `Tlist`
    in K for the 1D hindered rotor with a frequency `freq` in cm^-1 and a
    barrier height `barr` in cm^-1.
    """
    x = constants.h * constants.c * 100. * freq / constants.kB / T
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    z = 0.5 * constants.h * constants.c * 100. * barr / constants.kB / T
    BB = scipy.special.i1(z) / scipy.special.i0(z)
    return x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - BB - z * BB * BB)

def hinderedRotor_d_heatCapacity_d_freq(T, freq, barr):
    """
    Return the first derivative of the heat capacity with respect to the
    hindered rotor frequency in J/mol*K/cm^-1 at the given set of temperatures
    `Tlist` in K, evaluated at the frequency `freq` in cm^-1 and a barrier
    height `barr` in cm^-1.
    """
    x = constants.h * constants.c * 100. * freq / constants.kB / T
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    return x * exp_x / one_minus_exp_x / one_minus_exp_x * (2 + x + 2 * x * exp_x / one_minus_exp_x) * x / freq

def hinderedRotor_d_heatCapacity_d_barr(T, freq, barr):
    """
    Return the first derivative of the heat capacity with respect to the
    hindered rotor frequency in J/mol*K/cm^-1 at the given set of temperatures
    `Tlist` in K, evaluated at the frequency `freq` in cm^-1 and a barrier
    height `barr` in cm^-1.
    """
    z = 0.5 * constants.h * constants.c * 100. * barr / constants.kB / T
    BB = scipy.special.i1(z) / scipy.special.i0(z)
    return z * (1 - 2 * z * BB + BB * BB + 2 * z * BB * BB * BB) * z / barr

################################################################################

class DirectFit(DQED):
    """
    Class for fitting vibrational frequencies and hindered rotor
    frequency-barrier pairs for the case when there are few enough oscillators
    and rotors that their values can be fit directly.
    """

    def __init__(self, Tdata, Cvdata, Nvib, Nrot):
        self.Tdata = Tdata
        self.Cvdata = Cvdata
        self.Nvib = Nvib
        self.Nrot = Nrot

    def evaluate(self, x):
        Neq = self.Neq; Nvars = self.Nvars; Ncons = self.Ncons
        f = numpy.zeros((Neq), numpy.float64)
        J = numpy.zeros((Neq, Nvars), numpy.float64)
        fcons = numpy.zeros((Ncons), numpy.float64)
        Jcons = numpy.zeros((Ncons, Nvars), numpy.float64)

        Nvib = self.Nvib
        Nrot = self.Nrot

        for i in range(len(self.Tdata)):
            # Residual
            for n in range(Nvib):
                f[i] += harmonicOscillator_heatCapacity(self.Tdata[i], x[n])
            for n in range(Nrot):
                f[i] += hinderedRotor_heatCapacity(self.Tdata[i], x[Nvib+2*n], x[Nvib+2*n+1])
            f[i] -= self.Cvdata[i]
            # Jacobian
            for n in range(Nvib):
                J[i,n         ] = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata[i], x[n])
            for n in range(Nrot):
                J[i,Nvib+2*n  ] = hinderedRotor_d_heatCapacity_d_freq(self.Tdata[i], x[Nvib+2*n], x[Nvib+2*n+1])
                J[i,Nvib+2*n+1] = hinderedRotor_d_heatCapacity_d_barr(self.Tdata[i], x[Nvib+2*n], x[Nvib+2*n+1])
        
        return f, J, fcons, Jcons

class PseudoRotorFit(DQED):
    """
    Class for fitting vibrational frequencies and hindered rotor
    frequency-barrier pairs for the case when there are too many oscillators
    and rotors for their values can be fit directly, and where collapsing the
    rotors into a single pseudo-rotor allows for fitting the vibrational
    frequencies directly.
    """

    def __init__(self, Tdata, Cvdata, Nvib, Nrot):
        self.Tdata = Tdata
        self.Cvdata = Cvdata
        self.Nvib = Nvib
        self.Nrot = Nrot

    def evaluate(self, x):
        Neq = self.Neq; Nvars = self.Nvars; Ncons = self.Ncons
        f = numpy.zeros((Neq), numpy.float64)
        J = numpy.zeros((Neq, Nvars), numpy.float64)
        fcons = numpy.zeros((Ncons), numpy.float64)
        Jcons = numpy.zeros((Ncons, Nvars), numpy.float64)

        Nvib = self.Nvib
        Nrot = self.Nrot

        Cv = numpy.zeros((len(self.Tdata), Nvib+1), numpy.float64)
        dCv = numpy.zeros((len(self.Tdata), Nvib+2), numpy.float64)
    
        for i in range(len(self.Tdata)):
            for j in range(Nvib):
                Cv[i,j] = harmonicOscillator_heatCapacity(self.Tdata[i], x[j])
                dCv[i,j] = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata[i], x[j])
            Cv[i,Nvib] = hinderedRotor_heatCapacity(self.Tdata[i], x[Nvib], x[Nvib+1])
            dCv[i,Nvib] = hinderedRotor_d_heatCapacity_d_freq(self.Tdata[i], x[Nvib], x[Nvib+1])
            dCv[i,Nvib+1] = hinderedRotor_d_heatCapacity_d_barr(self.Tdata[i], x[Nvib], x[Nvib+1])

        for i in range(len(self.Tdata)):
            # Residual
            for j in range(Nvib):
                f[i] += Cv[i,j]
            f[i] += Nrot * Cv[i,Nvib]
            f[i] -= self.Cvdata[i]
            # Jacobian
            for j in range(Nvib):
                J[i,j] = 2.0 * f[i] * dCv[i,j]
            J[i,Nvib] = 2.0 * f[i] * Nrot * dCv[i,Nvib]
            J[i,Nvib+1] = 2.0 * f[i] * Nrot * dCv[i,Nvib+1]

        return f, J, fcons, Jcons

class PseudoFit(DQED):
    """
    Class for fitting vibrational frequencies and hindered rotor
    frequency-barrier pairs for the case when there are too many oscillators
    and rotors for their values can be fit directly, and where we must collapse
    both the vibrations and hindered rotations into "pseudo-oscillators" and
    "pseudo-rotors".
    """

    def __init__(self, Tdata, Cvdata, Nvib, Nrot):
        self.Tdata = Tdata
        self.Cvdata = Cvdata
        self.Nvib = Nvib
        self.Nrot = Nrot

    def evaluate(self, x):
        Neq = self.Neq; Nvars = self.Nvars; Ncons = self.Ncons
        f = numpy.zeros((Neq), numpy.float64)
        J = numpy.zeros((Neq, Nvars), numpy.float64)
        fcons = numpy.zeros((Ncons), numpy.float64)
        Jcons = numpy.zeros((Ncons, Nvars), numpy.float64)

        Nvib = self.Nvib
        Nrot = self.Nrot

        for i in range(len(self.Tdata)):
            Cv1 = harmonicOscillator_heatCapacity(self.Tdata[i], x[0])
            Cv2 = harmonicOscillator_heatCapacity(self.Tdata[i], x[2])
            Cv3 = harmonicOscillator_heatCapacity(self.Tdata[i], x[3])
            Cv4 = hinderedRotor_heatCapacity(self.Tdata[i], x[4], x[5])
            dCv1 = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata[i], x[0])
            dCv2 = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata[i], x[2])
            dCv3 = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata[i], x[3])
            dCv4 = hinderedRotor_d_heatCapacity_d_freq(self.Tdata[i], x[4], x[5])
            dCv5 = hinderedRotor_d_heatCapacity_d_barr(self.Tdata[i], x[4], x[5])

            # Residual
            f[i] = Cv1 + x[1] * Cv2 + (Nvib - x[1] - 1) * Cv3 + Nrot * Cv4 - self.Cvdata[i]

            # Jacobian
            J[i,0] = 2.0 * f[i] * dCv1
            J[i,1] = 2.0 * f[i] * (Cv2 - Cv3)
            J[i,2] = 2.0 * f[i] * x[1] * dCv2
            J[i,3] = 2.0 * f[i] * ((Nvib - x[1] - 1) * dCv3)
            J[i,4] = 2.0 * f[i] * Nrot * dCv4
            J[i,5] = 2.0 * f[i] * Nrot * dCv5

        return f, J, fcons, Jcons

################################################################################

def fitModes(mode, x0, bounds, maxIter, Tdata, Cvdata, Nvib, Nrot):
    """
    Fit a set of vibrational frequencies and hindered rotor frequency-barrier
    pairs. The `mode` of the fitting is ``'direct'``, ``'pseudo'``, or
    ``'pseudo-rotor'``. `x0` is the initial guess, `bounds` are the lower and
    upper bounds for each variable, `maxIter` is the maximum number of
    iterations, `Tdata` are the temperatures at which heat capacity data is
    known, `Cvdata` are the corresponding heat capacities in J/mol*K, `Nvib`
    is the total number of vibrational modes being fit, and `Nrot` is the total
    number of hindered rotors being fit. Returns the best fit parameters and
    the status flag from DQED.
    """

    if mode == 'direct':
        fit = DirectFit(Tdata, Cvdata, Nvib, Nrot)
    elif mode == 'pseudo':
        fit = PseudoFit(Tdata, Cvdata, Nvib, Nrot)
    elif mode == 'pseudo-rotors':
        fit = PseudoRotorFit(Tdata, Cvdata, Nvib, Nrot)

    fit.initialize(Neq=len(Tdata), Nvars=len(x0), Ncons=0, bounds=bounds, maxIter=maxIter)
    x, igo = fit.solve(x0)

    return x, igo
    