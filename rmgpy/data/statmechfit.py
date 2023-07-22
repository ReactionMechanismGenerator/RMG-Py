#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Contains functions for fitting of molecular degrees of freedom from
macroscopic properties, particularly the heat capacity.
"""
import logging
import math

import numpy as np
import scipy.special
from pydqed import DQED

import rmgpy.constants as constants
from rmgpy.exceptions import StatmechFitError
from rmgpy.statmech import HarmonicOscillator, HinderedRotor

################################################################################

# This block contains a number of global constants that control various
# aspects of the optimization (bounds, iterations, etc.)

# The lower bound for harmonic oscillator frequencies in cm^-1
ho_freq_lower_bound = 180.0
# The upper bound for harmonic oscillator frequencies in cm^-1
ho_freq_upper_bound = 4000.0

# The lower bound for hindered rotor frequencies in cm^-1
hr_freq_lower_bound = 180.0
# The upper bound for hindered rotor frequencies in cm^-1
hr_freq_upper_bound = 4000.0

# The lower bound for hindered rotor barrier heights in cm^-1
hr_barr_lower_bound = 10.0
# The upper bound for hindered rotor barrier heights in cm^-1
hr_barr_upper_bound = 10000.0

# The maximum number of iterations for the optimization solver to use
max_iter = 200


################################################################################


def fit_statmech_to_heat_capacity(Tlist, Cvlist, n_vib, n_rot, molecule=None):
    """
    For a given set of dimensionless heat capacity data `Cvlist` corresponding
    to temperature list `Tlist` in K, fit `n_vib` harmonic oscillator and `n_rot`
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
        raise StatmechFitError("You must specify at least 7 heat capacity points to fit_statmech_to_heat_capacity().")
    if len(Tlist) != len(Cvlist):
        raise StatmechFitError(
            "The number of heat capacity points ({0:d}) does not match the number of temperatures "
            "provided ({1:d}).".format(len(Cvlist), len(Tlist))
        )

    # The number of optimization variables available is constrained to be less
    # than the number of heat capacity points
    # This is also capped to a (somewhat arbitrarily chosen) maximum of 16
    max_variables = len(Tlist) - 1
    if max_variables > 16:
        max_variables = 16

    # The type of variables fitted depends on the values of n_vib and n_rot and
    # the number of heat capacity points provided
    # For low values of n_vib and n_rot, we can fit the individual
    # parameters directly
    # For high values of n_vib and/or n_rot we are limited by the number of
    # temperatures we are fitting at, and so we can only fit
    # pseudo-oscillators and/or pseudo-rotors
    vib = []
    hind = []
    if n_vib <= 0 and n_rot <= 0:
        pass
    elif n_vib + 2 * n_rot <= max_variables:
        vib, hind = fit_statmech_direct(Tlist, Cvlist, n_vib, n_rot, molecule)
    elif n_vib + 2 <= max_variables:
        vib, hind = fit_statmech_pseudo_rotors(Tlist, Cvlist, n_vib, n_rot, molecule)
    else:
        vib, hind = fit_statmech_pseudo(Tlist, Cvlist, n_vib, n_rot, molecule)

    modes = []
    if n_vib > 0:
        vib.sort()
        ho = HarmonicOscillator(frequencies=(vib[:], "cm^-1"))
        modes.append(ho)
    for i in range(n_rot):
        freq = hind[i][0]
        barr = hind[i][1]
        inertia = (barr * constants.c * 100.0 * constants.h) / (8 * math.pi * math.pi * (freq * constants.c * 100.0) ** 2)
        barrier = barr * constants.c * 100.0 * constants.h * constants.Na
        hr = HinderedRotor(
            inertia=(inertia * constants.Na * 1e23, "amu*angstrom^2"),
            barrier=(barrier / 1000.0, "kJ/mol"),
            symmetry=1,
            semiclassical=False,
            quantum=False,
        )
        modes.append(hr)

    # Return the fitted modes
    return modes


################################################################################


def fit_statmech_direct(Tlist, Cvlist, n_vib, n_rot, molecule=None):
    """
    Fit `n_vib` harmonic oscillator and `n_rot` hindered internal rotor modes to
    the provided dimensionless heat capacities `Cvlist` at temperatures `Tlist`
    in K. This method assumes that there are enough heat capacity points
    provided that the vibrational frequencies and hindered rotation frequency-
    barrier pairs can be fit directly.
    """

    # Construct the lower and upper bounds for each variable
    bounds = []
    # Bounds for harmonic oscillator frequencies
    for i in range(n_vib):
        bounds.append((ho_freq_lower_bound, ho_freq_upper_bound))
    # Bounds for hindered rotor frequencies and barrier heights
    for i in range(n_rot):
        bounds.append((hr_freq_lower_bound, hr_freq_upper_bound))
        bounds.append((hr_barr_lower_bound, hr_barr_upper_bound))

    # Construct the initial guess
    # Initial guesses within each mode type must be distinct or else the
    # optimization will fail
    x0 = np.zeros(n_vib + 2 * n_rot, float)
    # Initial guess for harmonic oscillator frequencies
    if n_vib > 0:
        x0[0] = 200.0
        x0[1:n_vib] = np.linspace(800.0, 1600.0, n_vib - 1)
    # Initial guess for hindered rotor frequencies and barrier heights
    if n_rot > 0:
        x0[n_vib] = 100.0
        x0[n_vib + 1] = 100.0
        for i in range(1, n_rot):
            x0[n_vib + 2 * i] = x0[n_vib + 2 * i - 2] + 20.0
            x0[n_vib + 2 * i + 1] = x0[n_vib + 2 * i - 1] + 100.0

    # Execute the optimization
    fit = DirectFit(Tlist, Cvlist, n_vib, n_rot)
    fit.initialize(Neq=len(Tlist), Nvars=len(x0), Ncons=0, bounds=bounds, maxIter=max_iter)
    x, igo = fit.solve(x0)

    # Check that the results of the optimization are valid
    if not np.isfinite(x).all():
        raise StatmechFitError("Returned solution vector is nonsensical: x = {0}.".format(x))
    if igo == 8:
        logging.warning("Maximum number of iterations reached when fitting spectral data for " "{0}.".format(molecule.to_smiles()))
    elif igo > 8:
        logging.warning("A solver error occured when fitting spectral data for {0}.".format(molecule.to_smiles()))
    logging.debug("Fitting remaining heat capacity to {0} vibrations and {1} rotations".format(n_vib, n_rot))
    logging.debug("The residuals for heat capacity values is {}".format(fit.evaluate(x)[0]))

    # Postprocess optimization results
    vib = list(x[0:n_vib])
    hind = []
    for i in range(n_rot):
        hind.append((x[n_vib + 2 * i], x[n_vib + 2 * i + 1]))

    return vib, hind


################################################################################


def fit_statmech_pseudo_rotors(Tlist, Cvlist, n_vib, n_rot, molecule=None):
    """
    Fit `n_vib` harmonic oscillator and `n_rot` hindered internal rotor modes to
    the provided dimensionless heat capacities `Cvlist` at temperatures `Tlist`
    in K. This method assumes that there are enough heat capacity points
    provided that the vibrational frequencies can be fit directly, but the
    hindered rotors must be combined into a single "pseudo-rotor".
    """

    # Construct the lower and upper bounds for each variable
    bounds = []
    # Bounds for harmonic oscillator frequencies
    for i in range(n_vib):
        bounds.append((ho_freq_lower_bound, ho_freq_upper_bound))
    # Bounds for pseudo-hindered rotor frequency and barrier height
    bounds.append((hr_freq_lower_bound, hr_freq_upper_bound))
    bounds.append((hr_barr_lower_bound, hr_barr_upper_bound))

    # Construct the initial guess
    # Initial guesses within each mode type must be distinct or else the
    # optimization will fail
    x0 = np.zeros(n_vib + 2, float)
    # Initial guess for harmonic oscillator frequencies
    if n_vib > 0:
        x0[0] = 200.0
        x0[1:n_vib] = np.linspace(800.0, 1600.0, n_vib - 1)
    # Initial guess for hindered rotor frequencies and barrier heights
    x0[n_vib] = 100.0
    x0[n_vib + 1] = 300.0

    # Execute the optimization
    fit = PseudoRotorFit(Tlist, Cvlist, n_vib, n_rot)
    fit.initialize(Neq=len(Tlist), Nvars=len(x0), Ncons=0, bounds=bounds, maxIter=max_iter)
    x, igo = fit.solve(x0)

    # Check that the results of the optimization are valid
    if not np.isfinite(x).all():
        raise StatmechFitError("Returned solution vector is nonsensical: x = {0}.".format(x))
    if igo == 8:
        logging.warning("Maximum number of iterations reached when fitting spectral data for " "{0}.".format(molecule.to_smiles()))
    if igo > 8:
        logging.warning("A solver error occured when fitting spectral data for {0}.".format(molecule.to_smiles()))
    logging.debug("Fitting remaining heat capacity to {0} vibrations and {1} rotations".format(n_vib, n_rot))
    logging.debug("The residuals for heat capacity values is {}".format(fit.evaluate(x)[0]))

    # Postprocess optimization results
    vib = list(x[0:n_vib])
    hind = []
    for i in range(n_rot):
        hind.append((x[n_vib], x[n_vib + 1]))

    return vib, hind


################################################################################


def fit_statmech_pseudo(Tlist, Cvlist, n_vib, n_rot, molecule=None):
    """
    Fit `n_vib` harmonic oscillator and `n_rot` hindered internal rotor modes to
    the provided dimensionless heat capacities `Cvlist` at temperatures `Tlist`
    in K. This method assumes that there are relatively few heat capacity points
    provided, so the vibrations must be combined into one real vibration and
    two "pseudo-vibrations" and the hindered rotors must be combined into a
    single "pseudo-rotor".
    """

    # Construct the lower and upper bounds for each variable
    bounds = []
    # x[0] corresponds to the first harmonic oscillator (real) frequency
    bounds.append((ho_freq_lower_bound, ho_freq_upper_bound))
    # x[1] corresponds to the degeneracy of the second harmonic oscillator
    bounds.append((1.0, float(n_vib - 2)))
    # x[2] corresponds to the second harmonic oscillator pseudo-frequency
    bounds.append((ho_freq_lower_bound, ho_freq_upper_bound))
    # x[3] corresponds to the third harmonic oscillator pseudo-frequency
    bounds.append((ho_freq_lower_bound, ho_freq_upper_bound))
    # x[4] corresponds to the hindered rotor pseudo-frequency
    bounds.append((hr_freq_lower_bound, hr_freq_upper_bound))
    # x[5] corresponds to the hindered rotor pseudo-barrier
    bounds.append((hr_barr_lower_bound, hr_barr_upper_bound))

    # Construct the initial guess
    x0 = np.zeros(6, float)  # Initial guess
    x0[0] = 300.0
    x0[1] = float(math.floor((n_vib - 1) / 2.0))
    x0[2] = 800.0
    x0[3] = 1600.0
    x0[4] = 100.0
    x0[5] = 300.0

    # Execute the optimization
    fit = PseudoFit(Tlist, Cvlist, n_vib, n_rot)
    fit.initialize(Neq=len(Tlist), Nvars=len(x0), Ncons=0, bounds=bounds, maxIter=max_iter)
    x, igo = fit.solve(x0)

    # Check that the results of the optimization are valid
    if not np.isfinite(x).all():
        raise StatmechFitError("Returned solution vector is nonsensical: x = {0}.".format(x))
    if igo == 8:
        logging.warning("Maximum number of iterations reached when fitting spectral data for " "{0}.".format(molecule.to_smiles()))
    if igo > 8:
        logging.warning("A solver error occured when fitting spectral data for {0}.".format(molecule.to_smiles()))
    logging.debug("Fitting remaining heat capacity to {0} vibrations and {1} rotations".format(n_vib, n_rot))
    logging.debug("The residuals for heat capacity values is {}".format(fit.evaluate(x)[0]))

    # Postprocess optimization results
    n_vib2 = int(round(x[1]))
    n_vib_3 = n_vib - n_vib2 - 1
    if n_vib2 < 0 or n_vib2 > n_vib - 1 or n_vib_3 < 0 or n_vib_3 > n_vib - 1:
        raise StatmechFitError("Invalid degeneracies {0} and {1} fitted for " "pseudo-frequencies.".format(n_vib2, n_vib_3))

    vib = [x[0]]
    for i in range(n_vib2):
        vib.append(x[2])
    for i in range(n_vib_3):
        vib.append(x[3])
    hind = []
    for i in range(n_rot):
        hind.append((x[4], x[5]))

    return vib, hind


################################################################################


def harmonic_oscillator_heat_capacity(T, freq):
    """
    Return the heat capacity in J/mol*K at the given set of temperatures `Tlist`
    in K for the harmonic oscillator with a frequency `freq` in cm^-1.
    """
    x = freq / (0.695039 * T)  # kB = 0.695039 cm^-1/K
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    return x * x * exp_x / one_minus_exp_x / one_minus_exp_x


def harmonic_oscillator_d_heat_capacity_d_freq(T, freq):
    """
    Return the first derivative of the heat capacity with respect to the
    harmonic oscillator frequency in J/mol*K/cm^-1 at the given set of
    temperatures `Tlist` in K, evaluated at the frequency `freq` in cm^-1.
    """
    x = freq / (0.695039 * T)  # kB = 0.695039 cm^-1/K
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    return x * exp_x / one_minus_exp_x / one_minus_exp_x * (2.0 + x + 2.0 * x * exp_x / one_minus_exp_x) * x / freq


def hindered_rotor_heat_capacity(T, freq, barr):
    """
    Return the heat capacity in J/mol*K at the given set of temperatures `Tlist`
    in K for the 1D hindered rotor with a frequency `freq` in cm^-1 and a
    barrier height `barr` in cm^-1.
    """
    x = constants.h * constants.c * 100.0 * freq / constants.kB / T
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    z = 0.5 * constants.h * constants.c * 100.0 * barr / constants.kB / T
    bb = scipy.special.i1(z) / scipy.special.i0(z)
    return x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - bb - z * bb * bb)


def hindered_rotor_d_heat_capacity_d_freq(T, freq, barr):
    """
    Return the first derivative of the heat capacity with respect to the
    hindered rotor frequency in J/mol*K/cm^-1 at the given set of temperatures
    `Tlist` in K, evaluated at the frequency `freq` in cm^-1 and a barrier
    height `barr` in cm^-1.
    """
    x = constants.h * constants.c * 100.0 * freq / constants.kB / T
    exp_x = math.exp(x)
    one_minus_exp_x = 1.0 - exp_x
    return x * exp_x / one_minus_exp_x / one_minus_exp_x * (2 + x + 2 * x * exp_x / one_minus_exp_x) * x / freq


def hindered_rotor_d_heat_capacity_d_barr(T, freq, barr):
    """
    Return the first derivative of the heat capacity with respect to the
    hindered rotor frequency in J/mol*K/cm^-1 at the given set of temperatures
    `Tlist` in K, evaluated at the frequency `freq` in cm^-1 and a barrier
    height `barr` in cm^-1.
    """
    z = 0.5 * constants.h * constants.c * 100.0 * barr / constants.kB / T
    bb = scipy.special.i1(z) / scipy.special.i0(z)
    return z * (1 - 2 * z * bb + bb * bb + 2 * z * bb * bb * bb) * z / barr


################################################################################


class DirectFit(DQED):
    """
    Class for fitting vibrational frequencies and hindered rotor
    frequency-barrier pairs for the case when there are few enough oscillators
    and rotors that their values can be fit directly.
    """

    def __init__(self, Tdata, Cvdata, n_vib, n_rot):
        self.Tdata = Tdata
        self.Cvdata = Cvdata
        self.n_vib = n_vib
        self.n_rot = n_rot

    def evaluate(self, x):
        n_eq = self.Neq
        n_vars = self.Nvars
        n_cons = self.Ncons
        f = np.zeros((n_eq), float)
        J = np.zeros((n_eq, n_vars), float)
        f_cons = np.zeros((n_cons), float)
        J_cons = np.zeros((n_cons, n_vars), float)

        n_vib = self.n_vib
        n_rot = self.n_rot

        for i in range(len(self.Tdata)):
            # Residual
            for n in range(n_vib):
                f[i] += harmonic_oscillator_heat_capacity(self.Tdata[i], x[n])
            for n in range(n_rot):
                f[i] += hindered_rotor_heat_capacity(self.Tdata[i], x[n_vib + 2 * n], x[n_vib + 2 * n + 1])
            f[i] -= self.Cvdata[i]
            # Jacobian
            for n in range(n_vib):
                J[i, n] = harmonic_oscillator_d_heat_capacity_d_freq(self.Tdata[i], x[n])
            for n in range(n_rot):
                J[i, n_vib + 2 * n] = hindered_rotor_d_heat_capacity_d_freq(self.Tdata[i], x[n_vib + 2 * n], x[n_vib + 2 * n + 1])
                J[i, n_vib + 2 * n + 1] = hindered_rotor_d_heat_capacity_d_barr(self.Tdata[i], x[n_vib + 2 * n], x[n_vib + 2 * n + 1])

        return f, J, f_cons, J_cons


class PseudoRotorFit(DQED):
    """
    Class for fitting vibrational frequencies and hindered rotor
    frequency-barrier pairs for the case when there are too many oscillators
    and rotors for their values can be fit directly, and where collapsing the
    rotors into a single pseudo-rotor allows for fitting the vibrational
    frequencies directly.
    """

    def __init__(self, Tdata, Cvdata, n_vib, n_rot):
        self.Tdata = Tdata
        self.Cvdata = Cvdata
        self.n_vib = n_vib
        self.n_rot = n_rot

    def evaluate(self, x):
        n_eq = self.Neq
        n_vars = self.Nvars
        n_cons = self.Ncons
        f = np.zeros((n_eq), float)
        J = np.zeros((n_eq, n_vars), float)
        f_cons = np.zeros((n_cons), float)
        J_cons = np.zeros((n_cons, n_vars), float)

        n_vib = self.n_vib
        n_rot = self.n_rot

        cv = np.zeros((len(self.Tdata), n_vib + 1), float)
        d_cv = np.zeros((len(self.Tdata), n_vib + 2), float)

        for i in range(len(self.Tdata)):
            for j in range(n_vib):
                cv[i, j] = harmonic_oscillator_heat_capacity(self.Tdata[i], x[j])
                d_cv[i, j] = harmonic_oscillator_d_heat_capacity_d_freq(self.Tdata[i], x[j])
            cv[i, n_vib] = hindered_rotor_heat_capacity(self.Tdata[i], x[n_vib], x[n_vib + 1])
            d_cv[i, n_vib] = hindered_rotor_d_heat_capacity_d_freq(self.Tdata[i], x[n_vib], x[n_vib + 1])
            d_cv[i, n_vib + 1] = hindered_rotor_d_heat_capacity_d_barr(self.Tdata[i], x[n_vib], x[n_vib + 1])

        for i in range(len(self.Tdata)):
            # Residual
            for j in range(n_vib):
                f[i] += cv[i, j]
            f[i] += n_rot * cv[i, n_vib]
            f[i] -= self.Cvdata[i]
            # Jacobian
            for j in range(n_vib):
                J[i, j] = 2.0 * f[i] * d_cv[i, j]
            J[i, n_vib] = 2.0 * f[i] * n_rot * d_cv[i, n_vib]
            J[i, n_vib + 1] = 2.0 * f[i] * n_rot * d_cv[i, n_vib + 1]

        return f, J, f_cons, J_cons


class PseudoFit(DQED):
    """
    Class for fitting vibrational frequencies and hindered rotor
    frequency-barrier pairs for the case when there are too many oscillators
    and rotors for their values can be fit directly, and where we must collapse
    both the vibrations and hindered rotations into "pseudo-oscillators" and
    "pseudo-rotors".
    """

    def __init__(self, Tdata, Cvdata, n_vib, n_rot):
        self.Tdata = Tdata
        self.Cvdata = Cvdata
        self.n_vib = n_vib
        self.n_rot = n_rot

    def evaluate(self, x):
        n_eq = self.Neq
        n_vars = self.Nvars
        n_cons = self.Ncons
        f = np.zeros((n_eq), float)
        J = np.zeros((n_eq, n_vars), float)
        f_cons = np.zeros((n_cons), float)
        J_cons = np.zeros((n_cons, n_vars), float)

        n_vib = self.n_vib
        n_rot = self.n_rot

        for i in range(len(self.Tdata)):
            cv1 = harmonic_oscillator_heat_capacity(self.Tdata[i], x[0])
            cv2 = harmonic_oscillator_heat_capacity(self.Tdata[i], x[2])
            cv3 = harmonic_oscillator_heat_capacity(self.Tdata[i], x[3])
            cv4 = hindered_rotor_heat_capacity(self.Tdata[i], x[4], x[5])
            d_cv1 = harmonic_oscillator_d_heat_capacity_d_freq(self.Tdata[i], x[0])
            d_cv2 = harmonic_oscillator_d_heat_capacity_d_freq(self.Tdata[i], x[2])
            d_cv3 = harmonic_oscillator_d_heat_capacity_d_freq(self.Tdata[i], x[3])
            d_cv4 = hindered_rotor_d_heat_capacity_d_freq(self.Tdata[i], x[4], x[5])
            d_cv5 = hindered_rotor_d_heat_capacity_d_barr(self.Tdata[i], x[4], x[5])

            # Residual
            f[i] = cv1 + x[1] * cv2 + (n_vib - x[1] - 1) * cv3 + n_rot * cv4 - self.Cvdata[i]

            # Jacobian
            J[i, 0] = 2.0 * f[i] * d_cv1
            J[i, 1] = 2.0 * f[i] * (cv2 - cv3)
            J[i, 2] = 2.0 * f[i] * x[1] * d_cv2
            J[i, 3] = 2.0 * f[i] * ((n_vib - x[1] - 1) * d_cv3)
            J[i, 4] = 2.0 * f[i] * n_rot * d_cv4
            J[i, 5] = 2.0 * f[i] * n_rot * d_cv5

        return f, J, f_cons, J_cons
