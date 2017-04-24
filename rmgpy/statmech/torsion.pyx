# cython: embedsignature=True, cdivision=True

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains classes that represent various models of torsional
motion. Torsional modes have both vibrational and rotational character, and
therefore should be treated semiclassically or quantum mechanically.
"""

import math
import numpy
import cython
import scipy.linalg
from scipy.misc import factorial
from scipy.special import i0, i1, ellipk, ellipe
from libc.math cimport log, exp, sqrt, sin, cos

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity
cimport rmgpy.statmech.schrodinger as schrodinger 
import rmgpy.statmech.schrodinger as schrodinger 

import logging

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

################################################################################

class NegativeBarrierException(Exception):
    """This Exception occurs when the energy barrier for a hindered Rotor is negative.
    This can occur if the scan or fourier fit is poor. """
    
    pass

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
    
    def __init__(self, inertia=None, symmetry=1, barrier=None, fourier=None, rotationalConstant=None, quantum=False, semiclassical=True):
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
        return (HinderedRotor, (self.inertia, self.symmetry, self.barrier, self.fourier, None, self.quantum, self.semiclassical))
    
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
            return quantity.Quantity(B,"cm^-1")
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
    
    cdef double getRotationalConstantEnergy(self):
        """
        Return the value of the rotational constant in J/mol.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia.value_si) * constants.Na

    cpdef double getFrequency(self) except -1:
        """
        Return the frequency of vibration in cm^-1 corresponding to the limit of
        harmonic oscillation.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=2] fourier
        cdef double V0, I, frequency
        cdef int k
        I = self._inertia.value_si
        if self.frequency != 0:
            return self.frequency
        elif self._fourier is not None:
            fourier = self._fourier.value_si
            V0 = 0.0
            for k in range(fourier.shape[1]):
                V0 -= fourier[0,k] * (k+1) * (k+1)
            V0 /= constants.Na
            if V0 < 0:
                raise NegativeBarrierException(" Hindered rotor barrier height is less than 0 \n     Try running cantherm in verbose mode, -v,  to identify which hindered rotor \n    Try changing the Hindered rotor fit to 'cosine'")
#                 raise Exception("Hindered rotor barrier height is less than 0 \nTry running cantherm in verbose mode, -v,  to identify which hindered rotor \nTry changing the Hindered rotor fit to 'cosine'")
            frequency = 1.0 / (2. * constants.pi) * sqrt(V0 / I)
        else:
            V0 = self._barrier.value_si / constants.Na
            if V0 < 0:
                raise Exception('V0 barrier height is less than 0')
            frequency = self.symmetry / (2. * constants.pi) * sqrt(V0 / (2. * I))
        self.frequency = frequency / (constants.c * 100.)
        return self.frequency

    cpdef double getLevelEnergy(self, int J) except -1:
        """
        Return the energy of level `J` in J.
        """
        return self.energies[J] if J < self.energies.shape[0] else 0.0
    
    cpdef int getLevelDegeneracy(self, int J) except -1:
        """
        Return the degeneracy of level `J`.
        """
        return 1
    
    cpdef numpy.ndarray solveSchrodingerEquation(self, int Nbasis=401):
        """
        Solves the one-dimensional time-independent Schrodinger equation to 
        determine the energy levels of a one-dimensional hindered rotor with a
        Fourier series potential using `Nbasis` basis functions. For the
        purposes of this function it is usually sufficient to use 401 basis
        functions (the default). Returns the energy eigenvalues of the
        Hamiltonian matrix in J/mol.
        """
        # Populate Hamiltonian matrix (banded in lower triangular form)
        H = self.getHamiltonian(Nbasis)
        # The overlap matrix is the identity matrix, i.e. this is a standard
        # eigenvalue problem
        
        # Find the eigenvalues and eigenvectors of the Hamiltonian matrix
        E = scipy.linalg.eig_banded(H, lower=True, eigvals_only=True, overwrite_a_band=True)
        # Don't consider zero-point energy here
        self.energies = E - numpy.min(E)
        
        # Return the eigenvalues
        return self.energies

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef numpy.ndarray getHamiltonian(self, int Nbasis):
        """
        Return the to the Hamiltonian matrix for the hindered rotor for the
        given number of basis functions `Nbasis`. The Hamiltonian matrix is
        returned in banded lower triangular form and with units of J/mol.
        """
        cdef int M, m, col, n
        cdef numpy.ndarray[numpy.float64_t,ndim=2] coeffs
        cdef double V0
    
        # The number of terms to use is 2*M + 1, ranging from -M to M inclusive
        if Nbasis % 2 == 0:
            M = Nbasis / 2
        else:
            M = (Nbasis - 1) / 2
        
        if self._fourier is not None:
            coeffs = self._fourier.value_si
            V0 = -numpy.sum(coeffs[0,:])
        else:
            coeffs = numpy.zeros((2,self.symmetry), numpy.float64)
            V0 = 0.5 * self._barrier.value_si
            coeffs[0,self.symmetry-1] = -V0
            
        # Populate Hamiltonian matrix (banded in lower triangular form)
        H = numpy.zeros((coeffs.shape[1]+1,2*M+1), numpy.complex64)
        B = self.getRotationalConstantEnergy()
        col = 0
        for m in range(-M, M+1):
            H[0,col] = B * m * m + V0
            for n in range(coeffs.shape[1]):
                H[n+1,col] = 0.5 * coeffs[0,n] - 0.5j * coeffs[1,n]
            col += 1

        return H

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double getPotential(self, double phi) except -100000000:
        """
        Return the value of the hindered rotor potential :math:`V(\\phi)`
        in J/mol at the angle `phi` in radians.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=2] fourier
        cdef double V = 0.0
        cdef int k
        if self._fourier is not None:
            fourier = self._fourier.value_si
            for k in range(fourier.shape[1]):
                V += fourier[0,k] * (cos((k+1) * phi) - 1.0) + fourier[1,k] * sin((k+1) * phi)
        else:
            V = 0.5 * self._barrier.value_si * (1 - cos(self.symmetry * phi))
        return V
    
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        cdef double frequency, x, z, Q
        cdef double beta = 1. / (constants.R * T), V, phi, dphi
        cdef int k
        
        frequency = self.getFrequency() * constants.c * 100
        x = constants.h * frequency / (constants.kB * T)
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            return numpy.sum(numpy.exp(-self.energies / constants.R / T)) / self.symmetry
        elif self._fourier is not None:
            # Fourier series data found, so use it
            # Numerically evaluate the configuration integral
            dphi = constants.pi/32.
            Q = 0.0
            for phi in numpy.arange(0, 2*constants.pi, dphi):
                Q += exp(-beta * self.getPotential(phi)) * dphi
            Q *= sqrt(constants.R * T / (4 * constants.pi * self.getRotationalConstantEnergy())) / self.symmetry
        else:
            # No Fourier data, so use the cosine potential data
            z = 0.5 * self._barrier.value_si / (constants.R * T)
            if z > 100:
                Q = sqrt(0.5 / z / constants.pi)
            else:
                Q = exp(-z) * i0(z)
            Q *= sqrt(constants.pi * constants.R * T / self.getRotationalConstantEnergy()) / self.symmetry
        
        # Semiclassical correction
        if self.semiclassical:
            Q *= x / (1 - exp(-x))
        
        return Q 
  
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double frequency, x, z, exp_x, one_minus_exp_x, BB, Cv
        cdef double Tlow, Thigh, logQlow, logQhigh, logQ
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            E = self.energies
            e_kT = numpy.exp(-E / constants.R / T)
            Cv = (numpy.sum(E*E*e_kT) * numpy.sum(e_kT) - numpy.sum(E*e_kT)**2) / (constants.R*constants.R*T*T * numpy.sum(e_kT)**2)
        elif self._fourier is not None:
            # Fourier series data found, so use it
            Tlow = T * 0.999
            Thigh = T * 1.001
            logQlow = log(self.getPartitionFunction(Tlow))
            logQhigh = log(self.getPartitionFunction(Thigh))
            logQ = log(self.getPartitionFunction(T))
            Cv = T * T * (logQhigh - 2 * logQ + logQlow) / ((Thigh - T) * (T - Tlow)) + 2 * T * (logQhigh - logQlow) / (Thigh - Tlow)
        else:
            # No Fourier data, so use the cosine potential data
            frequency = self.getFrequency() * constants.c * 100
            x = constants.h * frequency / (constants.kB * T)
            z = 0.5 * self._barrier.value_si / (constants.R * T)
            exp_x = exp(x)
            one_minus_exp_x = 1.0 - exp_x
            BB = i1(z) / i0(z)
            Cv = (x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - BB - z * BB * BB))
        return Cv * constants.R
    
    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Tlow, Thigh
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            E = self.energies
            e_kT = numpy.exp(-E / constants.R / T)
            return numpy.sum(E*e_kT) / numpy.sum(e_kT)
        else:
            # No Fourier data, so use the cosine potential data
            Tlow = T * 0.999
            Thigh = T * 1.001
            return (T *
                (log(self.getPartitionFunction(Thigh)) -
                log(self.getPartitionFunction(Tlow))) /
                (Thigh - Tlow)) * constants.R * T

    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Tlow, Thigh
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            E = self.energies
            e_kT = numpy.exp(-E / constants.R / T)
            return numpy.log(self.getPartitionFunction(T)) * constants.R + numpy.sum(E*e_kT) / (T * numpy.sum(e_kT))
        else:
            # No Fourier data, so use the cosine potential data
            Tlow = T * 0.999
            Thigh = T * 1.001
            return (log(self.getPartitionFunction(T)) +
                T * (log(self.getPartitionFunction(Thigh)) -
                log(self.getPartitionFunction(Tlow))) /
                (Thigh - Tlow)) * constants.R 

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in J/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] sumStates, _Elist = Elist
        cdef double q1f, pre, V0
        cdef int i
            
        if sumStates0 is not None:
            return schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist))
        elif self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            sumStates = schrodinger.getSumOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, sumStates0) / self.symmetry
        elif self.semiclassical:
            raise NotImplementedError
        elif self._fourier is not None:
            # Fourier series data found, so use it
            raise NotImplementedError
        else:
            # No Fourier data, so use the cosine potential data
            sumStates = numpy.zeros_like(_Elist)
            q1f = sqrt(constants.pi / self.getRotationalConstantEnergy()) / self.symmetry
            V0 = self._barrier.value_si
            pre = 4.0 * q1f * sqrt(V0 / (constants.pi * constants.pi * constants.pi))
            # The following is only valid in the classical limit
            # Note that ellipk(1) = infinity, so we must skip that value
            for i in range(_Elist.shape[0]):
                if _Elist[i] < V0:
                    sumStates[i] = pre * (ellipe(_Elist[i] / V0) - (1 - _Elist[i] / V0) * ellipk(_Elist[i] / V0))
                elif _Elist[i] > V0:
                    sumStates[i] = pre * sqrt(_Elist[i] / V0) * ellipe(V0 / _Elist[i])
        return sumStates
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in J/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] densStates, _Elist = Elist
        cdef double q1f, pre, V0
        cdef int i
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            densStates = schrodinger.getDensityOfStates(_Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, densStates0) / self.symmetry
        elif self.semiclassical:
            raise NotImplementedError
        elif self._fourier is not None:
            # Fourier series data found, so use it
            raise NotImplementedError
        else:
            # No Fourier data, so use the cosine potential data
            densStates = numpy.zeros_like(Elist)
            q1f = sqrt(constants.pi / self.getRotationalConstantEnergy()) / self.symmetry
            V0 = self._barrier.value_si
            pre = 2.0 * q1f / sqrt(constants.pi * constants.pi * constants.pi * V0)
            # The following is only valid in the classical limit
            # Note that ellipk(1) = infinity, so we must skip that value
            for i in range(_Elist.shape[0]):
                if _Elist[i] < V0:
                    densStates[i] = pre * ellipk(_Elist[i] / V0)
                elif _Elist[i] > V0:
                    densStates[i] = pre * sqrt(V0 / _Elist[i]) * ellipk(V0 / _Elist[i])
            densStates *= _Elist[1] - _Elist[0]
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef fitFourierPotentialToData(self, numpy.ndarray angle, numpy.ndarray V):
        """
        Fit the given angles in radians and corresponding potential energies in
        J/mol to the Fourier series potential. For best results, the angle
        should begin at zero and end at :math:`2 \pi`, with the minimum energy
        conformation having a potential of zero be placed at zero angle.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=2] A
        cdef numpy.ndarray[numpy.float64_t,ndim=1] b
        cdef double phi
        cdef int N, i, m
        
        # Fit Fourier series potential
        N = V.shape[0]
        A = numpy.zeros((N+1,12), numpy.float64)
        b = numpy.zeros(N+1, numpy.float64)
        for i in range(N):
            phi = angle[i]
            for m in range(6):
                A[i,m] = cos(m * phi)
                A[i,6+m] = sin(m * phi)
                b[i] = V[i]
        # This row forces dV/dangle = 0 at angle = 0
        for m in range(6):
            A[N,m+6] = m
        x, residues, rank, s = numpy.linalg.lstsq(A, b)
        x *= 0.001
        
        self.fourier = ([x[1:6], x[7:12]], "kJ/mol")
        self.barrier = None
        
        return self

    cpdef fitCosinePotentialToData(self, numpy.ndarray angle, numpy.ndarray V):
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
        symmetry = 0; barrier = 0.0
        for symm in range(1, 10):
            num = numpy.sum(V * (1 - numpy.cos(symm * angle)))
            den = numpy.sum((1 - numpy.cos(symm * angle))**2)
            barr = 2 * num / den
            if barr > barrier:
                symmetry = symm
                barrier = barr

        self.fourier = None
        self.barrier = (barrier*0.001,"kJ/mol")
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
    def __init__(self,inertia=None,symmetry=1,rotationalConstant=None):
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
            return quantity.Quantity(B,"cm^-1")
        def __set__(self, B):
            cdef double I
            B = quantity.Frequency(B)
            I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
            self._inertia = quantity.ScalarQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")
    
    cdef double getRotationalConstantEnergy(self):
        """
        Return the value of the rotational constant in J/mol.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia.value_si) * constants.Na
    
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        return numpy.sqrt(8*numpy.pi**3*constants.kb*T*self._inertia.value_si)/(self.symmetry*constants.h)
        
  
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        return constants.R/2.0
       
    
    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        return constants.R*T/2.0
        

    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Q
        Q = self.getPartitionFunction(T)
        return constants.R*(numpy.log(Q)+.5)
        
        


    

