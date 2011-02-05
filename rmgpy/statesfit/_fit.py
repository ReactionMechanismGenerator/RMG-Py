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

def harmonicOscillator_heatCapacity(Tlist, freq):
    """
    Return the heat capacity in J/mol*K at the given set of temperatures `Tlist`
    in K for the harmonic oscillator with a frequency `freq` in cm^-1.
    """
    nT = len(Tlist)
    Cvlist = numpy.zeros_like(Tlist)
    for i in range(nT):
        x = freq / (0.695039 * Tlist[i])        # kB = 0.695039 cm^-1/K
        exp_x = math.exp(x)
        one_minus_exp_x = 1.0 - exp_x
        Cvlist[i] = x * x * exp_x / one_minus_exp_x / one_minus_exp_x
    return Cvlist

def harmonicOscillator_d_heatCapacity_d_freq(Tlist, freq):
    """
    Return the first derivative of the heat capacity with respect to the
    harmonic oscillator frequency in J/mol*K/cm^-1 at the given set of
    temperatures `Tlist` in K, evaluated at the frequency `freq` in cm^-1.
    """
    nT = len(Tlist)
    dCvlist = numpy.zeros_like(Tlist)

    for i in range(nT):
        x = freq / (0.695039 * Tlist[i])        # kB = 0.695039 cm^-1/K
        exp_x = math.exp(x)
        one_minus_exp_x = 1.0 - exp_x
        dCvlist[i] = x * exp_x / one_minus_exp_x / one_minus_exp_x * (2.0 + x + 2.0 * x * exp_x / one_minus_exp_x) * x / freq
    return dCvlist

def hinderedRotor_heatCapacity(Tlist, freq, barr):
    """
    Return the heat capacity in J/mol*K at the given set of temperatures `Tlist`
    in K for the 1D hindered rotor with a frequency `freq` in cm^-1 and a
    barrier height `barr` in cm^-1.
    """
    nT = len(Tlist)
    Cvlist = numpy.zeros_like(Tlist)
    for i in range(nT):
        x = constants.h * constants.c * 100. * freq / constants.kB / Tlist[i]
        exp_x = math.exp(x)
        one_minus_exp_x = 1.0 - exp_x
        z = 0.5 * constants.h * constants.c * 100. * barr / constants.kB / Tlist[i]
        BB = scipy.special.i1(z) / scipy.special.i0(z)
        Cvlist[i] = x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - BB - z * BB * BB)
    return Cvlist

def hinderedRotor_d_heatCapacity_d_freq(Tlist, freq, barr):
    """
    Return the first derivative of the heat capacity with respect to the
    hindered rotor frequency in J/mol*K/cm^-1 at the given set of temperatures
    `Tlist` in K, evaluated at the frequency `freq` in cm^-1 and a barrier
    height `barr` in cm^-1.
    """
    nT = len(Tlist)
    dCvlist = numpy.zeros_like(Tlist)
    for i in range(nT):
        x = constants.h * constants.c * 100. * freq / constants.kB / Tlist[i]
        exp_x = math.exp(x)
        one_minus_exp_x = 1.0 - exp_x
        dCvlist[i] = x * exp_x / one_minus_exp_x / one_minus_exp_x * (2 + x + 2 * x * exp_x / one_minus_exp_x) * x / freq
    return dCvlist

def hinderedRotor_d_heatCapacity_d_barr(Tlist, freq, barr):
    """
    Return the first derivative of the heat capacity with respect to the
    hindered rotor frequency in J/mol*K/cm^-1 at the given set of temperatures
    `Tlist` in K, evaluated at the frequency `freq` in cm^-1 and a barrier
    height `barr` in cm^-1.
    """
    nT = len(Tlist)
    dCvlist = numpy.zeros_like(Tlist)
    for i in range(nT):
        z = 0.5 * constants.h * constants.c * 100. * barr / constants.kB / Tlist[i]
        BB = scipy.special.i1(z) / scipy.special.i0(z)
        dCvlist[i] = z * (1 - 2 * z * BB + BB * BB + 2 * z * BB * BB * BB) * z / barr
    return dCvlist

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

        # Residual
        for n in range(Nvib):
            f += harmonicOscillator_heatCapacity(self.Tdata, x[n])
        for n in range(Nrot):
            f += hinderedRotor_heatCapacity(self.Tdata, x[Nvib+2*n], x[Nvib+2*n+1])
        f -= self.Cvdata

        # Jacobian
        for n in range(Nvib):
            J[:,n         ] = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata, x[n])
        for n in range(Nrot):
            J[:,Nvib+2*n  ] = hinderedRotor_d_heatCapacity_d_freq(self.Tdata, x[Nvib+2*n], x[Nvib+2*n+1])
            J[:,Nvib+2*n+1] = hinderedRotor_d_heatCapacity_d_barr(self.Tdata, x[Nvib+2*n], x[Nvib+2*n+1])
        
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

        Cv = numpy.zeros((Nvib+1, len(self.Tdata)), numpy.float64)
        dCv = numpy.zeros((Nvib+2, len(self.Tdata)), numpy.float64)
    
        for i in range(Nvib):
            Cv[i,:] = harmonicOscillator_heatCapacity(self.Tdata, x[i])
            dCv[i,:] = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata, x[i])
        Cv[Nvib,:] = hinderedRotor_heatCapacity(self.Tdata, x[Nvib], x[Nvib+1])
        dCv[Nvib,:] = hinderedRotor_d_heatCapacity_d_freq(self.Tdata, x[Nvib], x[Nvib+1])
        dCv[Nvib+1,:] = hinderedRotor_d_heatCapacity_d_barr(self.Tdata, x[Nvib], x[Nvib+1])

        # Residual
        for i in range(Nvib):
            f += Cv[i,:]
        f += Nrot * Cv[Nvib,:]
        f -= self.Cvdata

        # Jacobian
        for i in range(Nvib):
            J[:,i] = 2.0 * f * dCv[i,:]
        J[:,Nvib] = 2.0 * f * Nrot * dCv[Nvib,:]
        J[:,Nvib+1] = 2.0 * f * Nrot * dCv[Nvib+1,:]

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

        Cv1 = harmonicOscillator_heatCapacity(self.Tdata, x[0])
        Cv2 = harmonicOscillator_heatCapacity(self.Tdata, x[2])
        Cv3 = harmonicOscillator_heatCapacity(self.Tdata, x[3])
        Cv4 = hinderedRotor_heatCapacity(self.Tdata, x[4], x[5])
        dCv1 = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata, x[0])
        dCv2 = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata, x[2])
        dCv3 = harmonicOscillator_d_heatCapacity_d_freq(self.Tdata, x[3])
        dCv4 = hinderedRotor_d_heatCapacity_d_freq(self.Tdata, x[4], x[5])
        dCv5 = hinderedRotor_d_heatCapacity_d_barr(self.Tdata, x[4], x[5])

        # Residual
        f = Cv1 + x[1] * Cv2 + (Nvib - x[1] - 1) * Cv3 + Nrot * Cv4 - self.Cvdata

        # Jacobian
        J[:,0] = 2.0 * f * dCv1
        J[:,1] = 2.0 * f * (Cv2 - Cv3)
        J[:,2] = 2.0 * f * x[1] * dCv2
        J[:,3] = 2.0 * f * ((Nvib - x[1] - 1) * dCv3)
        J[:,4] = 2.0 * f * Nrot * dCv4
        J[:,5] = 2.0 * f * Nrot * dCv5

        return f, J, fcons, Jcons

################################################################################

def fitModes(mode, x0, bl, bu, ind, maxIter, Tdata, Cvdata, Nvib, Nrot):

    bounds = [(l,u) for l, u in zip(bl, bu)]

    if mode == 'direct':
        fit = DirectFit(Tdata, Cvdata, Nvib, Nrot)
    elif mode == 'pseudo':
        fit = PseudoFit(Tdata, Cvdata, Nvib, Nrot)
    elif mode == 'pseudo-rotors':
        fit = PseudoRotorFit(Tdata, Cvdata, Nvib, Nrot)

    fit.initialize(Neq=len(Tdata), Nvars=Nvib+2*Nrot, Ncons=0, bounds=bounds, maxIter=maxIter)
    x, igo = fit.solve(x0)

    return x, igo
    