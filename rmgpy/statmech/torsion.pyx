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
This module contains classes that represent various models of torsional
motion. Torsional modes have both vibrational and rotational character, and
therefore should be treated semiclassically or quantum mechanically.
"""
import logging

import cython
import numpy as np
cimport numpy as np
import scipy.linalg
from libc.math cimport log, exp, sqrt, sin, cos
from scipy.special import i0, i1, ellipk, ellipe

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity
import rmgpy.statmech.schrodinger as schrodinger
cimport rmgpy.statmech.schrodinger as schrodinger
from rmgpy.exceptions import NegativeBarrierException
from rmgpy.rmgobject import recursive_make_object

# Prior to numpy 1.14, `numpy.linalg.lstsq` does not accept None as a value
RCOND = -1 if int(np.__version__.split('.')[1]) < 14 else None

################################################################################

cdef class Torsion(Mode):
    """
    A base class for all torsional degrees of freedom. The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================

    In the majority of chemical applications, the torsional energy levels are
    mostly close together compared to :math:`k_\\mathrm{B} T`, which makes the
    torsional motion well-approximated by a semiclassical treatment. 
    """

    def __init__(self, symmetry=1, quantum=False):
        Mode.__init__(self, quantum)
        self.symmetry = symmetry

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Rotation object.
        """
        result = 'Torsion(symmetry={0:d}'.format(self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a Rotation object.
        """
        return (Torsion, (self.symmetry, self.quantum))

    cpdef make_object(self, dict data, dict class_dict):
        kwargs = recursive_make_object(data, class_dict, make_final_object=False)
        if ('inertia' in kwargs) and ('rotationalConstant' in kwargs):  # Only one of these can be specified
            del kwargs['inertia']
        self.__init__(**kwargs)

################################################################################

cdef class HinderedRotor(Torsion):
    """
    A statistical mechanical model of a one-dimensional hindered rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `fourier`                The :math:`2 x N` array of Fourier series coefficients
    `barrier`                The barrier height of the cosine potential
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    `semiclassical`          ``True`` to use the semiclassical correction, ``False`` otherwise
    ======================== ===================================================

    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.
    """

    def __init__(self, inertia=None, symmetry=1, barrier=None, fourier=None, rotationalConstant=None, quantum=True,
                 semiclassical=True, frequency=None, energies=None):
        Torsion.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.barrier = barrier
        self.fourier = fourier
        self.semiclassical = False if quantum else semiclassical
        self.quantum = quantum
        self.frequency = frequency if frequency is not None else 0.0
        self.energies = energies

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        HinderedRotor object.
        """
        result = 'HinderedRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.fourier is not None:
            result += ', fourier={0!r}'.format(self.fourier)
        if self.barrier is not None:
            result += ', barrier={0!r}'.format(self.barrier)
        if self.quantum:
            result += ', quantum=True'
        if not self.semiclassical:
            result += ', semiclassical=False'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a HinderedRotor object.
        """
        return (HinderedRotor, (self.inertia, self.symmetry, self.barrier, self.fourier, None, self.quantum,
                                self.semiclassical))

    property inertia:
        """The moment of inertia of the rotor."""
        def __get__(self):
            return self._inertia
        def __set__(self, value):
            self._inertia = quantity.Inertia(value)

    property rotationalConstant:
        """The rotational constant of the rotor."""
        def __get__(self):
            cdef double I = self._inertia.value_si
            cdef double B = constants.h / (8 * constants.pi * constants.pi * I) / (constants.c * 100.)
            return quantity.Quantity(B, "cm^-1")
        def __set__(self, B):
            cdef double I
            B = quantity.Frequency(B)
            I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
            self._inertia = quantity.ScalarQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")

    property fourier:
        """The :math:`2 x N` array of Fourier series coefficients."""
        def __get__(self):
            return self._fourier
        def __set__(self, value):
            self._fourier = quantity.Energy(value)

    property barrier:
        """The barrier height of the cosine potential."""
        def __get__(self):
            return self._barrier
        def __set__(self, value):
            self._barrier = quantity.Energy(value)

    cdef double get_rotational_constant_energy(self):
        """
        Return the value of the rotational constant in J/mol.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia.value_si) * constants.Na

    cpdef double get_frequency(self) except -1:
        """
        Return the frequency of vibration in cm^-1 corresponding to the limit of
        harmonic oscillation.
        """
        cdef np.ndarray[np.float64_t, ndim=2] fourier
        cdef double V0, I, frequency
        cdef int k
        I = self._inertia.value_si
        if self.frequency != 0:
            return self.frequency
        elif self._fourier is not None:
            fourier = self._fourier.value_si
            V0 = 0.0
            for k in range(fourier.shape[1]):
                V0 -= fourier[0, k] * (k + 1) * (k + 1)
            V0 /= constants.Na
            if V0 < 0:
                raise NegativeBarrierException("Hindered rotor barrier height is less than 0 \n     Try running Arkane"
                                               "in verbose mode, -v, to identify which rotor caused the error\n"
                                               "Also, try changing the Hindered rotor fit to 'cosine'")
            frequency = 1.0 / (2. * constants.pi) * sqrt(V0 / I)
        else:
            V0 = self._barrier.value_si / constants.Na
            if V0 < 0:
                raise Exception('V0 barrier height is less than 0')
            frequency = self.symmetry / (2. * constants.pi) * sqrt(V0 / (2. * I))
        self.frequency = frequency / (constants.c * 100.)
        return self.frequency

    cpdef double get_level_energy(self, int J) except -1:
        """
        Return the energy of level `J` in J.
        """
        return self.energies[J] if J < self.energies.shape[0] else 0.0

    cpdef int get_level_degeneracy(self, int J) except -1:
        """
        Return the degeneracy of level `J`.
        """
        return 1

    cpdef np.ndarray solve_schrodinger_equation(self, int n_basis=401):
        """
        Solves the one-dimensional time-independent Schrodinger equation to 
        determine the energy levels of a one-dimensional hindered rotor with a
        Fourier series potential using `n_basis` basis functions. For the
        purposes of this function it is usually sufficient to use 401 basis
        functions (the default). Returns the energy eigenvalues of the
        Hamiltonian matrix in J/mol.
        """
        # Populate Hamiltonian matrix (banded in lower triangular form)
        H = self.get_hamiltonian(n_basis)
        # The overlap matrix is the identity matrix, i.e. this is a standard
        # eigenvalue problem

        # Find the eigenvalues and eigenvectors of the Hamiltonian matrix
        E = scipy.linalg.eig_banded(H, lower=True, eigvals_only=True, overwrite_a_band=True)
        # Don't consider zero-point energy here
        self.energies = E - np.min(E)

        # Return the eigenvalues
        return self.energies

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef np.ndarray get_hamiltonian(self, int n_basis):
        """
        Return the to the Hamiltonian matrix for the hindered rotor for the
        given number of basis functions `n_basis`. The Hamiltonian matrix is
        returned in banded lower triangular form and with units of J/mol.
        """
        cdef int M, m, col, n
        cdef np.ndarray[np.float64_t, ndim=2] coeffs
        cdef double V0

        # The number of terms to use is 2*M + 1, ranging from -M to M inclusive
        if n_basis % 2 == 0:
            M = n_basis / 2
        else:
            M = (n_basis - 1) / 2

        if self._fourier is not None:
            coeffs = self._fourier.value_si
            V0 = -np.sum(coeffs[0, :])
        else:
            coeffs = np.zeros((2, self.symmetry), np.float64)
            V0 = 0.5 * self._barrier.value_si
            coeffs[0, self.symmetry - 1] = -V0

        # Populate Hamiltonian matrix (banded in lower triangular form)
        H = np.zeros((coeffs.shape[1] + 1, 2 * M + 1), np.complex64)
        B = self.get_rotational_constant_energy()
        col = 0
        for m in range(-M, M + 1):
            H[0, col] = B * m * m + V0
            for n in range(coeffs.shape[1]):
                H[n + 1, col] = 0.5 * coeffs[0, n] - 0.5j * coeffs[1, n]
            col += 1

        return H

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double get_potential(self, double phi) except -100000000:
        """
        Return the value of the hindered rotor potential :math:`V(\\phi)`
        in J/mol at the angle `phi` in radians.
        """
        cdef np.ndarray[np.float64_t, ndim=2] fourier
        cdef double V = 0.0
        cdef int k
        if self._fourier is not None:
            fourier = self._fourier.value_si
            for k in range(fourier.shape[1]):
                V += fourier[0, k] * (cos((k + 1) * phi) - 1.0) + fourier[1, k] * sin((k + 1) * phi)
        else:
            V = 0.5 * self._barrier.value_si * (1 - cos(self.symmetry * phi))
        return V

    cpdef double get_partition_function(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        cdef double frequency, x, z, Q
        cdef double beta = 1. / (constants.R * T), V, phi, dphi
        cdef int k

        frequency = self.get_frequency() * constants.c * 100
        x = constants.h * frequency / (constants.kB * T)

        if self.quantum:
            if self.energies is None: self.solve_schrodinger_equation()
            return np.sum(np.exp(-self.energies / constants.R / T)) / self.symmetry
        elif self._fourier is not None:
            # Fourier series data found, so use it
            # Numerically evaluate the configuration integral
            dphi = constants.pi / 32.
            Q = 0.0
            for phi in np.arange(0, 2 * constants.pi, dphi):
                Q += exp(-beta * self.get_potential(phi)) * dphi
            Q *= sqrt(constants.R * T / (4 * constants.pi * self.get_rotational_constant_energy())) / self.symmetry
        else:
            # No Fourier data, so use the cosine potential data
            z = 0.5 * self._barrier.value_si / (constants.R * T)
            if z > 100:
                Q = sqrt(0.5 / z / constants.pi)
            else:
                Q = exp(-z) * i0(z)
            Q *= sqrt(constants.pi * constants.R * T / self.get_rotational_constant_energy()) / self.symmetry

        # Semiclassical correction
        if self.semiclassical:
            Q *= x / (1 - exp(-x))

        return Q

    cpdef double get_heat_capacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double frequency, x, z, exp_x, one_minus_exp_x, BB, Cv
        cdef double Tlow, Thigh, logQlow, logQhigh, logQ

        if self.quantum:
            if self.energies is None: self.solve_schrodinger_equation()
            E = self.energies
            e_kT = np.exp(-E / constants.R / T)
            Cv = (np.sum(E*E*e_kT) * np.sum(e_kT) - np.sum(E*e_kT)**2) / (constants.R*constants.R*T*T * np.sum(e_kT)**2)
        elif self._fourier is not None:
            # Fourier series data found, so use it
            Tlow = T * 0.999
            Thigh = T * 1.001
            logQlow = log(self.get_partition_function(Tlow))
            logQhigh = log(self.get_partition_function(Thigh))
            logQ = log(self.get_partition_function(T))
            Cv = T * T * (logQhigh - 2 * logQ + logQlow) / ((Thigh - T) * (T - Tlow)) + 2 * T * (logQhigh - logQlow) / (Thigh - Tlow)
        else:
            # No Fourier data, so use the cosine potential data
            frequency = self.get_frequency() * constants.c * 100
            x = constants.h * frequency / (constants.kB * T)
            z = 0.5 * self._barrier.value_si / (constants.R * T)
            exp_x = exp(x)
            one_minus_exp_x = 1.0 - exp_x
            BB = i1(z) / i0(z)
            Cv = (x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - BB - z * BB * BB))
        return Cv * constants.R

    cpdef double get_enthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Tlow, Thigh

        if self.quantum:
            if self.energies is None: self.solve_schrodinger_equation()
            E = self.energies
            e_kT = np.exp(-E / constants.R / T)
            return np.sum(E * e_kT) / np.sum(e_kT)
        else:
            # No Fourier data, so use the cosine potential data
            Tlow = T * 0.999
            Thigh = T * 1.001
            return (T *
                    (log(self.get_partition_function(Thigh)) -
                     log(self.get_partition_function(Tlow))) /
                    (Thigh - Tlow)) * constants.R * T

    cpdef double get_entropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Tlow, Thigh

        if self.quantum:
            if self.energies is None: self.solve_schrodinger_equation()
            E = self.energies
            e_kT = np.exp(-E / constants.R / T)
            return np.log(self.get_partition_function(T)) * constants.R + np.sum(E * e_kT) / (T * np.sum(e_kT))
        else:
            # No Fourier data, so use the cosine potential data
            Tlow = T * 0.999
            Thigh = T * 1.001
            return (log(self.get_partition_function(T)) +
                    T * (log(self.get_partition_function(Thigh)) -
                         log(self.get_partition_function(Tlow))) /
                    (Thigh - Tlow)) * constants.R

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. If an initial sum of states 
        `sum_states_0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef np.ndarray[np.float64_t, ndim=1] sum_states, _e_list = e_list
        cdef double q1f, pre, V0
        cdef int i

        if sum_states_0 is not None:
            return schrodinger.convolve(sum_states_0, self.get_density_of_states(e_list))
        elif self.quantum:
            if self.energies is None: self.solve_schrodinger_equation()
            sum_states = schrodinger.get_sum_of_states(e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                       sum_states_0) / self.symmetry
        elif self.semiclassical:
            raise NotImplementedError
        elif self._fourier is not None:
            # Fourier series data found, so use it
            raise NotImplementedError
        else:
            # No Fourier data, so use the cosine potential data
            sum_states = np.zeros_like(_e_list)
            q1f = sqrt(constants.pi / self.get_rotational_constant_energy()) / self.symmetry
            V0 = self._barrier.value_si
            pre = 4.0 * q1f * sqrt(V0 / (constants.pi * constants.pi * constants.pi))
            # The following is only valid in the classical limit
            # Note that ellipk(1) = infinity, so we must skip that value
            for i in range(_e_list.shape[0]):
                if _e_list[i] < V0:
                    sum_states[i] = pre * (ellipe(_e_list[i] / V0) - (1 - _e_list[i] / V0) * ellipk(_e_list[i] / V0))
                elif _e_list[i] > V0:
                    sum_states[i] = pre * sqrt(_e_list[i] / V0) * ellipe(V0 / _e_list[i])
        return sum_states

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state. If an initial density
        of states `dens_states_0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef np.ndarray[np.float64_t, ndim=1] dens_states, _e_list = e_list
        cdef double q1f, pre, V0
        cdef int i

        if self.quantum:
            if self.energies is None: self.solve_schrodinger_equation()
            dens_states = schrodinger.get_density_of_states(_e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                            dens_states_0) / self.symmetry
        elif self.semiclassical:
            raise NotImplementedError
        elif self._fourier is not None:
            # Fourier series data found, so use it
            raise NotImplementedError
        else:
            # No Fourier data, so use the cosine potential data
            dens_states = np.zeros_like(e_list)
            q1f = sqrt(constants.pi / self.get_rotational_constant_energy()) / self.symmetry
            V0 = self._barrier.value_si
            pre = 2.0 * q1f / sqrt(constants.pi * constants.pi * constants.pi * V0)
            # The following is only valid in the classical limit
            # Note that ellipk(1) = infinity, so we must skip that value
            for i in range(_e_list.shape[0]):
                if _e_list[i] < V0:
                    dens_states[i] = pre * ellipk(_e_list[i] / V0)
                elif _e_list[i] > V0:
                    dens_states[i] = pre * sqrt(V0 / _e_list[i]) * ellipk(V0 / _e_list[i])
            dens_states *= _e_list[1] - _e_list[0]
            if dens_states_0 is not None:
                dens_states = schrodinger.convolve(dens_states_0, dens_states)
        return dens_states

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef fit_fourier_potential_to_data(self, np.ndarray angle, np.ndarray V):
        """
        Fit the given angles in radians and corresponding potential energies in
        J/mol to the Fourier series potential. For best results, the angle
        should begin at zero and end at :math:`2 \pi`, with the minimum energy
        conformation having a potential of zero be placed at zero angle.
        """
        cdef np.ndarray[np.float64_t, ndim=2] A, fourier
        cdef np.ndarray[np.float64_t, ndim=1] b, fit
        cdef double phi, value, V0
        cdef int N, i, m, numterms, maxterms
        numterms = 6
        cdef bint negative_barrier
        negative_barrier = True
        # numterms is actually half the number of terms. It is called numterms
        # because it is the number of terms of either the cosine or sine fit

        maxterms = np.floor(len(angle) / 3.0)
        while negative_barrier and numterms <= maxterms:
            # Fit Fourier series potential
            N = V.shape[0]
            # A: [1, cos(phi), ..., cos(M * phi), sin(phi), ..., sin(M * phi)]
            A = np.zeros((N + 1, 2 * numterms - 1), np.float64)
            A[:-1, 0] = 1
            for m in range(1, numterms):
                A[:-1, m] = np.cos(m * angle)
                A[:-1, numterms + m - 1] = np.sin(m * angle)
            # This row forces dV/dangle = 0 at angle = 0
            A[N, numterms:] = np.arange(1., numterms)
            b = np.concatenate((V, np.array([0.])))

            x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)
            fit = np.dot(A, x)
            x *= 0.001
            # This checks if there are any negative values in the Fourier fit.
            negative_barrier = False
            V0 = 0.0
            self.fourier = ([x[1:numterms], x[numterms:2 * numterms - 1]], "kJ/mol")
            fourier = self._fourier.value_si
            for k in range(fourier.shape[1]):
                V0 -= fourier[0, k] * (k + 1) * (k + 1)
            if V0 < 0:
                negative_barrier = True
                logging.warning("Fourier fit for hindered rotor gave a negative barrier when fit with {0} terms, "
                                "retrying with {1} terms...".format(2 * numterms, 2 * numterms + 4))
                numterms = numterms + 2
        if V0 < 0:
            logging.error("Fourier fit for hindered rotor gave a negative barrier on final try with "
                          "{0} terms".format(numterms * 2))

        self.barrier = None
        return self

    cpdef fit_cosine_potential_to_data(self, np.ndarray angle, np.ndarray V):
        """
        Fit the given angles in radians and corresponding potential energies in
        J/mol to the cosine potential. For best results, the angle should 
        begin at zero and end at :math:`2 \pi`, with the minimum energy
        conformation having a potential of zero be placed at zero angle. The
        fit is attempted at several possible values of the symmetry number in
        order to determine which one is correct.
        """
        cdef double barrier, barr, num, den
        cdef int symmetry, symm

        # We fit at integral symmetry numbers in the range [1, 9]
        # The best fit will have the maximum barrier height
        symmetry = 0
        barrier = 0.0
        for symm in range(1, 10):
            num = np.sum(V * (1 - np.cos(symm * angle)))
            den = np.sum((1 - np.cos(symm * angle)) ** 2)
            barr = 2 * num / den
            if barr > barrier:
                symmetry = symm
                barrier = barr

        self.fourier = None
        self.barrier = (barrier * 0.001, "kJ/mol")
        self.symmetry = symmetry

        return self

cdef class FreeRotor(Torsion):
    """
    A statistical mechanical model of a one-dimensional hindered rotor.  
    Based on Pfaendtner et al. 2007.  
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    ======================== ===================================================

    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.
    """
    def __init__(self, inertia=None, symmetry=1, rotationalConstant=None):
        Torsion.__init__(self, symmetry, False)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        FreeRotor object.
        """
        result = 'FreeRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a FreeRotor object.
        """
        return (FreeRotor, (self.inertia, self.symmetry))

    property inertia:
        """The moment of inertia of the rotor."""
        def __get__(self):
            return self._inertia
        def __set__(self, value):
            self._inertia = quantity.Inertia(value)

    property rotationalConstant:
        """The rotational constant of the rotor."""
        def __get__(self):
            cdef double I = self._inertia.value_si
            cdef double B = constants.h / (8 * constants.pi * constants.pi * I) / (constants.c * 100.)
            return quantity.Quantity(B, "cm^-1")
        def __set__(self, B):
            cdef double I
            B = quantity.Frequency(B)
            I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
            self._inertia = quantity.ScalarQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")

    cdef double get_rotational_constant_energy(self):
        """
        Return the value of the rotational constant in J/mol.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia.value_si) * constants.Na

    cpdef double get_partition_function(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        return np.sqrt(8 * np.pi ** 3 * constants.kB * T * self._inertia.value_si) / (self.symmetry * constants.h)

    cpdef double get_heat_capacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        return constants.R / 2.0

    cpdef double get_enthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        return constants.R * T / 2.0

    cpdef double get_entropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Q
        Q = self.get_partition_function(T)
        return constants.R * (np.log(Q) + .5)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. 
        formula from
        Forst 1995 Journal of Computational Chemistry, Vol. 17, No. 8 954-961 (1996)
        """
        if sum_states_0:
            raise NotImplementedError

        A = constants.hbar / (2.0 * self.inertia.value_si)
        return 2.0 / self.symmetry * np.sqrt(e_list / A)
