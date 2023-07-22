# cython: embedsignature=True, cdivision=True

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
This module contains classes that represent various models of vibrational
motion. For most molecular systems, a quantum treatment of vibrational motion
is required since typical vibrational energies are of similar magnitude to
:math:`k_\\mathrm{B} T`.
"""

import numpy as np
cimport numpy as np
from libc.math cimport log, exp
from scipy.special import factorial

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity
import rmgpy.statmech.schrodinger as schrodinger
cimport rmgpy.statmech.schrodinger as schrodinger

################################################################################

cdef class Vibration(Mode):
    """
    A base class for all vibrational degrees of freedom. The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================

    In the majority of chemical applications, the vibrational energy levels are
    widely spaced compared to :math:`k_\\mathrm{B} T`, which makes a quantum
    mechanical treatment required.
    """

    def __init__(self, quantum=True):
        Mode.__init__(self, quantum)

################################################################################

cdef class HarmonicOscillator(Vibration):
    """
    A statistical mechanical model of a set of one-dimensional independent
    harmonic oscillators. The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `frequencies`   The vibrational frequencies of the oscillators
    `quantum`       ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    =============== ============================================================
    
    In the majority of chemical applications, the energy levels of the
    harmonic oscillator are of similar magnitude to :math:`k_\\mathrm{B} T`,
    requiring a quantum mechanical treatment. Fortunately, the harmonic
    oscillator has an analytical quantum mechanical solution.
    """

    def __init__(self, frequencies=None, quantum=True):
        Vibration.__init__(self, quantum)
        self.frequencies = quantity.Frequency(frequencies)
        self.quantum = quantum

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        HarmonicOscillator object.
        """
        result = 'HarmonicOscillator(frequencies={0!r}'.format(self.frequencies)
        if not self.quantum:
            result += ', quantum=False'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a HarmonicOscillator object.
        """
        return (HarmonicOscillator, (self.frequencies, self.quantum))

    property frequencies:
        """The vibrational frequencies of the oscillators."""
        def __get__(self):
            return self._frequencies
        def __set__(self, value):
            self._frequencies = quantity.Frequency(value)

    cpdef double get_partition_function(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        cdef double Q = 1.0, beta = 1.0 / (constants.kB * T), freq
        cdef int i
        cdef np.ndarray[float_t, ndim=1] frequencies
        frequencies = self._frequencies.value_si
        if self.quantum:
            for i in range(frequencies.shape[0]):
                freq = frequencies[i] * constants.c * 100.
                Q *= 1.0 / (1 - exp(-beta * constants.h * freq))
        else:
            for i in range(frequencies.shape[0]):
                freq = frequencies[i] * constants.c * 100.
                Q *= 1.0 / (beta * constants.h * freq)
        return Q

    cpdef double get_heat_capacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Cv = 0.0, beta = 1.0 / (constants.kB * T), freq, x, exp_x
        cdef int i
        cdef np.ndarray[float_t, ndim=1] frequencies
        frequencies = self._frequencies.value_si
        if self.quantum:
            for i in range(frequencies.shape[0]):
                freq = frequencies[i] * constants.c * 100.
                x = beta * constants.h * freq
                if x > 500.0:
                    # exp(x) approaches infinity, thus x^2 exp(x)/(1-exp(x))^2 tends to zero
                    Cv += 0.0
                else:
                    exp_x = exp(x)
                    Cv += x * x * exp_x / (1 - exp_x) / (1 - exp_x)
        else:
            Cv = frequencies.shape[0]
        return Cv * constants.R

    cpdef double get_enthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double H = 0.0, beta = 1.0 / (constants.kB * T), freq, x
        cdef int i
        cdef np.ndarray[float_t, ndim=1] frequencies
        frequencies = self._frequencies.value_si
        if self.quantum:
            for i in range(frequencies.shape[0]):
                freq = frequencies[i] * constants.c * 100.
                x = beta * constants.h * freq
                H += x / (exp(x) - 1)
        else:
            H = frequencies.shape[0]
        return H * constants.R * T

    cpdef double get_entropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double S = 0.0, beta = 1.0 / (constants.kB * T), freq, x
        cdef int i
        cdef np.ndarray[float_t, ndim=1] frequencies
        S = log(self.get_partition_function(T))
        frequencies = self._frequencies.value_si
        if self.quantum:
            for i in range(frequencies.shape[0]):
                freq = frequencies[i] * constants.c * 100.
                x = beta * constants.h * freq
                S += x / (exp(x) - 1)
        else:
            S += frequencies.shape[0]
        return S * constants.R

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. If an initial sum of states 
        `sum_states_0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double freq
        cdef int i, Nfreq
        cdef np.ndarray[float_t, ndim=1] frequencies
        cdef np.ndarray sum_states
        frequencies = self._frequencies.value_si
        if self.quantum:
            if sum_states_0 is None:
                sum_states = np.ones_like(e_list)
            else:
                sum_states = sum_states_0
            for i in range(frequencies.shape[0]):
                freq = frequencies[i] * constants.c * 100.
                sum_states = schrodinger.convolve_bs(e_list, sum_states, constants.h * freq * constants.Na, 1)
        elif sum_states_0 is not None:
            sum_states = schrodinger.convolve(sum_states_0, self.get_density_of_states(e_list))
        else:
            Nfreq = frequencies.shape[0]
            sum_states = e_list ** Nfreq / factorial(Nfreq)
            for i in range(Nfreq):
                freq = frequencies[i] * constants.c * 100.
                sum_states /= constants.h * freq * constants.Na
        return sum_states

    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state. If an initial density
        of states `dens_states_0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double freq, dE
        cdef int i, Nfreq
        cdef np.ndarray[float_t, ndim=1] frequencies
        cdef np.ndarray dens_states
        frequencies = self._frequencies.value_si
        if self.quantum:
            if dens_states_0 is None:
                dens_states = np.zeros_like(e_list)
                dens_states[0] = 1.0
            else:
                dens_states = dens_states_0
            for i in range(frequencies.shape[0]):
                freq = frequencies[i] * constants.c * 100.
                dens_states = schrodinger.convolve_bs(e_list, dens_states, constants.h * freq * constants.Na, 1)
        else:
            Nfreq = frequencies.shape[0]
            dE = e_list[1] - e_list[0]
            dens_states = e_list ** (Nfreq - 1) / factorial(Nfreq - 1) * dE
            for i in range(Nfreq):
                freq = frequencies[i] * constants.c * 100.
                dens_states /= constants.h * freq * constants.Na
            if dens_states_0 is not None:
                dens_states = schrodinger.convolve(dens_states_0, dens_states)
        return dens_states
