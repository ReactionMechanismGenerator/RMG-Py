#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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
Contains functions for fitting of molecular degrees of freedom from
macroscopic properties, particularly the heat capacity.
"""

import math
import numpy
import scipy.special
import logging

import rmgpy.constants as constants
from rmgpy.statmech import HarmonicOscillator, HinderedRotor
from pydqed import DQED

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
maxIter = 200

################################################################################

class StatmechFitError(Exception):
    """
    An exception used when attempting to fit molecular degrees of freedom to
    heat capacity data. Pass a string describing the circumstances of the
    exceptional behavior.
    """
    pass

################################################################################

def fitStatmechToHeatCapacity(Tlist, Cvlist, Nvib, Nrot, molecule=None):
    """
    For a given set of dimensionless heat capacity data `Cvlist` corresponding
    to temperature list `Tlist` in K, fit `Nvib` harmonic oscillator and `Nrot`
    hindered internal rotor modes. External and other previously-known modes
    should have already been removed from `Cvlist` prior to calling this
    function. You must provide at least 7 values for `Cvlist`.

    This function returns a list containing the fitted vibrational frequencies
    in a :class:`HarmonicOscillator` object and the fitted 1D hindered rotors
    in :class:`HinderedRotor` objects.
    """

    # You must specify at least 7 heat capacity points to use in the fitting;
    # you can specify as many as you like above that minimum
    if len(Tlist) < 7:
        raise StatmechFitError('You must specify at least 7 heat capacity points to fitStatmechToHeatCapacity().')
    if len(Tlist) != len(Cvlist):
        raise StatmechFitError('The number of heat capacity points ({0:d}) does not match the number of temperatures provided ({1:d}).'.format(len(Cvlist), len(Tlist)))

    # The number of optimization variables available is constrained to be less
    # than the number of heat capacity points
    # This is also capped to a (somewhat arbitrarily chosen) maximum of 16
    maxVariables = len(Tlist) - 1
    if maxVariables > 16: maxVariables = 16

    # The type of variables fitted depends on the values of Nvib and Nrot and
    # the number of heat capacity points provided
    # For low values of Nvib and Nrot, we can fit the individual
    # parameters directly
    # For high values of Nvib and/or Nrot we are limited by the number of
    # temperatures we are fitting at, and so we can only fit
    # pseudo-oscillators and/or pseudo-rotors
    vib = []; hind = []
    if Nvib <= 0 and Nrot <= 0:
        pass
    elif Nvib + 2 * Nrot <= maxVariables:
        vib, hind = fitStatmechDirect(Tlist, Cvlist, Nvib, Nrot, molecule)
    elif Nvib + 2 <= maxVariables:
        vib, hind = fitStatmechPseudoRotors(Tlist, Cvlist, Nvib, Nrot, molecule)
    else:
        vib, hind = fitStatmechPseudo(Tlist, Cvlist, Nvib, Nrot, molecule)

    modes = []
    if Nvib > 0:
        vib.sort()
        ho = HarmonicOscillator(frequencies=(vib[:],"cm^-1"))
        modes.append(ho)
    for i in range(Nrot):
        freq = hind[i][0]
        barr = hind[i][1]
        inertia = (barr*constants.c*100.0*constants.h) / (8 * math.pi * math.pi * (freq*constants.c*100.0)**2)
        barrier = barr*constants.c*100.0*constants.h*constants.Na
        hr = HinderedRotor(inertia=(inertia*constants.Na*1e23,"amu*angstrom^2"), barrier=(barrier/1000.,"kJ/mol"), symmetry=1, semiclassical=False, quantum=False)
        modes.append(hr)

    # Return the fitted modes
    return modes

################################################################################

def fitStatmechDirect(Tlist, Cvlist, Nvib, Nrot, molecule=None):
    """
    Fit `Nvib` harmonic oscillator and `Nrot` hindered internal rotor modes to
    the provided dimensionless heat capacities `Cvlist` at temperatures `Tlist`
    in K. This method assumes that there are enough heat capacity points
    provided that the vibrational frequencies and hindered rotation frequency-
    barrier pairs can be fit directly.
    """

    # Construct the lower and upper bounds for each variable
    bounds = []
    # Bounds for harmonic oscillator frequencies
    for i in range(Nvib):
        bounds.append((hoFreqLowerBound, hoFreqUpperBound))
    # Bounds for hindered rotor frequencies and barrier heights
    for i in range(Nrot):
        bounds.append((hrFreqLowerBound, hrFreqUpperBound))
        bounds.append((hrBarrLowerBound, hrBarrUpperBound))

    # Construct the initial guess
    # Initial guesses within each mode type must be distinct or else the
    # optimization will fail
    x0 = numpy.zeros(Nvib + 2*Nrot, numpy.float64)
    # Initial guess for harmonic oscillator frequencies
    if Nvib > 0:
        x0[0] = 200.0
        x0[1:Nvib] = numpy.linspace(800.0, 1600.0, Nvib-1)
    # Initial guess for hindered rotor frequencies and barrier heights
    if Nrot > 0:
        x0[Nvib] = 100.0
        x0[Nvib+1] = 100.0
        for i in range(1, Nrot):
            x0[Nvib+2*i] = x0[Nvib+2*i-2] + 20.0
            x0[Nvib+2*i+1] = x0[Nvib+2*i-1] + 100.0

    # Execute the optimization
    fit = DirectFit(Tlist, Cvlist, Nvib, Nrot)
    fit.initialize(Neq=len(Tlist), Nvars=len(x0), Ncons=0, bounds=bounds, maxIter=maxIter)
    x, igo = fit.solve(x0)

    # Check that the results of the optimization are valid
    if not numpy.isfinite(x).all():
        raise StatmechFitError('Returned solution vector is nonsensical: x = {0}.'.format(x))
    if igo == 8:
        logging.warning('Maximum number of iterations reached when fitting spectral data for {0}.'.format(molecule.toSMILES()))

    # Postprocess optimization results
    vib = list(x[0:Nvib])
    hind = []
    for i in range(Nrot):
        hind.append((x[Nvib+2*i], x[Nvib+2*i+1]))

    return vib, hind

################################################################################

def fitStatmechPseudoRotors(Tlist, Cvlist, Nvib, Nrot, molecule=None):
    """
    Fit `Nvib` harmonic oscillator and `Nrot` hindered internal rotor modes to
    the provided dimensionless heat capacities `Cvlist` at temperatures `Tlist`
    in K. This method assumes that there are enough heat capacity points
    provided that the vibrational frequencies can be fit directly, but the
    hindered rotors must be combined into a single "pseudo-rotor".
    """

    # Construct the lower and upper bounds for each variable
    bounds = []
    # Bounds for harmonic oscillator frequencies
    for i in range(Nvib):
        bounds.append((hoFreqLowerBound, hoFreqUpperBound))
    # Bounds for pseudo-hindered rotor frequency and barrier height
    bounds.append((hrFreqLowerBound, hrFreqUpperBound))
    bounds.append((hrBarrLowerBound, hrBarrUpperBound))

    # Construct the initial guess
    # Initial guesses within each mode type must be distinct or else the
    # optimization will fail
    x0 = numpy.zeros(Nvib + 2, numpy.float64)
    # Initial guess for harmonic oscillator frequencies
    if Nvib > 0:
        x0[0] = 200.0
        x0[1:Nvib] = numpy.linspace(800.0, 1600.0, Nvib-1)
    # Initial guess for hindered rotor frequencies and barrier heights
    x0[Nvib] = 100.0
    x0[Nvib+1] = 300.0

    # Execute the optimization
    fit = PseudoRotorFit(Tlist, Cvlist, Nvib, Nrot)
    fit.initialize(Neq=len(Tlist), Nvars=len(x0), Ncons=0, bounds=bounds, maxIter=maxIter)
    x, igo = fit.solve(x0)

    # Check that the results of the optimization are valid
    if not numpy.isfinite(x).all():
        raise StatmechFitError('Returned solution vector is nonsensical: x = {0}.'.format(x))
    if igo == 8:
        logging.warning('Maximum number of iterations reached when fitting spectral data for {0}.'.format(molecule.toSMILES()))

    # Postprocess optimization results
    vib = list(x[0:Nvib])
    hind = []
    for i in range(Nrot):
        hind.append((x[Nvib], x[Nvib+1]))

    return vib, hind

################################################################################

def fitStatmechPseudo(Tlist, Cvlist, Nvib, Nrot, molecule=None):
    """
    Fit `Nvib` harmonic oscillator and `Nrot` hindered internal rotor modes to
    the provided dimensionless heat capacities `Cvlist` at temperatures `Tlist`
    in K. This method assumes that there are relatively few heat capacity points
    provided, so the vibrations must be combined into one real vibration and
    two "pseudo-vibrations" and the hindered rotors must be combined into a
    single "pseudo-rotor".
    """


    # Construct the lower and upper bounds for each variable
    bounds = []
    # x[0] corresponds to the first harmonic oscillator (real) frequency
    bounds.append((hoFreqLowerBound, hoFreqUpperBound))
    # x[1] corresponds to the degeneracy of the second harmonic oscillator
    bounds.append((1.0, float(Nvib - 2)))
    # x[2] corresponds to the second harmonic oscillator pseudo-frequency
    bounds.append((hoFreqLowerBound, hoFreqUpperBound))
    # x[3] corresponds to the third harmonic oscillator pseudo-frequency
    bounds.append((hoFreqLowerBound, hoFreqUpperBound))
    # x[4] corresponds to the hindered rotor pseudo-frequency
    bounds.append((hrFreqLowerBound, hrFreqUpperBound))
    # x[5] corresponds to the hindered rotor pseudo-barrier
    bounds.append((hrBarrLowerBound, hrBarrUpperBound))

    # Construct the initial guess
    x0 = numpy.zeros(6, numpy.float64)      # Initial guess
    x0[0] = 300.0
    x0[1] = float(math.floor((Nvib - 1) / 2.0))
    x0[2] = 800.0
    x0[3] = 1600.0
    x0[4] = 100.0
    x0[5] = 300.0

    # Execute the optimization
    fit = PseudoFit(Tlist, Cvlist, Nvib, Nrot)
    fit.initialize(Neq=len(Tlist), Nvars=len(x0), Ncons=0, bounds=bounds, maxIter=maxIter)
    x, igo = fit.solve(x0)

    # Check that the results of the optimization are valid
    if not numpy.isfinite(x).all():
        raise StatmechFitError('Returned solution vector is nonsensical: x = {0}.'.format(x))
    if igo == 8:
        logging.warning('Maximum number of iterations reached when fitting spectral data for {0}.'.format(molecule.toSMILES()))

    # Postprocess optimization results
    Nvib2 = int(round(x[1]))
    Nvib3 = Nvib - Nvib2 - 1
    if Nvib2 < 0 or Nvib2 > Nvib-1 or Nvib3 < 0 or Nvib3 > Nvib-1:
        raise StatmechFitError('Invalid degeneracies {0} and {1} fitted for pseudo-frequencies.'.format(Nvib2, Nvib3))

    vib = [x[0]]
    for i in range(Nvib2): vib.append(x[2])
    for i in range(Nvib3): vib.append(x[3])
    hind = []
    for i in range(Nrot):
        hind.append((x[4], x[5]))

    return vib, hind

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
