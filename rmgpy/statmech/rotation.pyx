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
This module contains classes that represent various models of rotational
motion. For most molecular systems, a classical treatment of rotational motion
is sufficient since typical rotational energies are much smaller in magnitude
than :math:`k_\\mathrm{B} T`.
"""

import math
import numpy
from libc.math cimport log, sqrt

cimport rmgpy.constants as constants
import rmgpy.statmech.schrodinger as schrodinger 
cimport rmgpy.statmech.schrodinger as schrodinger 
import rmgpy.quantity as quantity

################################################################################

cdef class Rotation(Mode):
    """
    A base class for all rotational degrees of freedom. The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================

    In the majority of chemical applications, the rotational energy levels are
    very close together compared to :math:`k_\\mathrm{B} T`, which makes the
    rotational motion well-approximated by the classical limit at all relevant
    temperatures. Therefore, the classical model is used by default. 
    """

    def __init__(self, symmetry=1, quantum=False):
        Mode.__init__(self, quantum)
        self.symmetry = symmetry

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Rotation object.
        """
        result = 'Rotation(symmetry={0:d}'.format(self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a Rotation object.
        """
        return (Rotation, (self.symmetry, self.quantum))

################################################################################

cdef inline double getRotationalConstantEnergy(double inertia) except -1:
    """
    Return the value of the rotational constant in energy units (J/mol)
    corresponding to the given moment of `inertia` in kg*m^2.
    """
    return constants.hbar * constants.hbar / (2. * inertia) * constants.Na

################################################################################

cdef class LinearRotor(Rotation):
    """
    A statistical mechanical model of a two-dimensional (linear) rigid rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.

    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """
    
    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        LinearRotor object.
        """
        result = 'LinearRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a LinearRotor object.
        """
        return (LinearRotor, (self.inertia, self.symmetry, self.quantum, None))

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

    cpdef double getLevelEnergy(self, int J) except -1:
        """
        Return the energy of level `J` in kJ/mol.
        """
        return getRotationalConstantEnergy(self.inertia.value_si) * J * (J + 1)
    
    cpdef int getLevelDegeneracy(self, int J) except -1:
        """
        Return the degeneracy of level `J`.
        """
        return 2 * J + 1
    
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        cdef double B, Q
        if self.quantum:
            Q = schrodinger.getPartitionFunction(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) / self.symmetry
        else:
            B = getRotationalConstantEnergy(self._inertia.value_si)
            Q = constants.R * T / B / self.symmetry
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            Cv = schrodinger.getHeatCapacity(T, self.getLevelEnergy, self.getLevelDegeneracy, 0)
        else:
            Cv = 1.0
        return Cv * constants.R

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            H = schrodinger.getEnthalpy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0)
        else:
            H = 1.0
        return H * constants.R * T
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            S = schrodinger.getEntropy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) - log(self.symmetry)
        else:
            S = log(self.getPartitionFunction(T)) + 1.0
        return S * constants.R

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in J/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double B
        cdef numpy.ndarray sumStates
        if self.quantum:
            sumStates = schrodinger.getSumOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, sumStates0) / self.symmetry
        elif sumStates0 is not None:
            sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
        else:
            B = getRotationalConstantEnergy(self._inertia.value_si)
            sumStates = Elist / B / self.symmetry
        return sumStates
            
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in J/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double B, dE
        cdef numpy.ndarray densStates
        if self.quantum:
            densStates = schrodinger.getDensityOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, densStates0) / self.symmetry
        else:
            B = getRotationalConstantEnergy(self._inertia.value_si)
            dE = Elist[1] - Elist[0]
            densStates = numpy.ones_like(Elist) * dE / B / self.symmetry
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates

################################################################################

cdef class NonlinearRotor(Rotation):
    """
    A statistical mechanical model of an N-dimensional nonlinear rigid rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moments of inertia of the rotor
    `rotationalConstant`     The rotational constants of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================

    Note that the moments of inertia and the rotational constants are simply two
    ways of representing the same quantity; only one set of these can be 
    specified independently.
    
    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. In the
    current implementation, the quantum mechanical model has not been 
    implemented, and a :class:`NotImplementedError` will be raised if you try
    to use it. 
    """
    
    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        NonlinearRotor object.
        """
        result = 'NonlinearRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a NonlinearRotor object.
        """
        return (NonlinearRotor, (self.inertia, self.symmetry, self.quantum, None))

    property inertia:
        """The moments of inertia of the rotor."""
        def __get__(self):
            return self._inertia
        def __set__(self, value):
            self._inertia = quantity.Inertia(value)
    
    property rotationalConstant:
        """The rotational constant of the rotor."""
        def __get__(self):
            cdef numpy.ndarray I = self._inertia.value_si
            cdef numpy.ndarray B = constants.h / (8 * constants.pi * constants.pi * I) / (constants.c * 100.)
            return quantity.Quantity(B,"cm^-1")
        def __set__(self, B):
            cdef numpy.ndarray I
            B = quantity.Frequency(B)
            I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
            self._inertia = quantity.ArrayQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")

    cdef numpy.ndarray getRotationalConstantEnergy(self):
        """
        Return the values of the rotational constants in kJ/mol.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia.value_si) * constants.Na

    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        cdef double Q, theta = 1.0
        cdef int i
        cdef numpy.ndarray[numpy.float64_t,ndim=1] rotationalConstants
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            rotationalConstants = self.getRotationalConstantEnergy()
            for i in range(rotationalConstants.shape[0]):
                theta *= rotationalConstants[i] / constants.R
            Q = sqrt(constants.pi * T**rotationalConstants.shape[0] / theta) / self.symmetry
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            Cv = 0.5 * self._inertia.value_si.shape[0]
        return Cv * constants.R

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            H = 0.5 * self._inertia.value_si.shape[0]
        return H * constants.R * T
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            S = numpy.log(self.getPartitionFunction(T)) + 0.5 * self._inertia.value_si.shape[0]
        return S * constants.R
    
    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in J/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double theta = 1.0
        cdef int i
        cdef numpy.ndarray[numpy.float64_t,ndim=1] rotationalConstants
        cdef numpy.ndarray sumStates
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            if sumStates0 is not None:
                sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
            else:
                rotationalConstants = self.getRotationalConstantEnergy()
                assert rotationalConstants.shape[0] == 3
                for i in range(rotationalConstants.shape[0]):
                    theta *= rotationalConstants[i]
                sumStates = 4.0/3.0 * Elist * numpy.sqrt(Elist / theta) / self.symmetry
        return sumStates
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in J/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double theta = 1.0, inertia, dE
        cdef int i
        cdef numpy.ndarray[numpy.float64_t,ndim=1] rotationalConstants
        cdef numpy.ndarray densStates
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            dE = Elist[1] - Elist[0]
            rotationalConstants = self.getRotationalConstantEnergy()
            assert rotationalConstants.shape[0] == 3
            for i in range(rotationalConstants.shape[0]):
                theta *= rotationalConstants[i]
            densStates = 2.0 * numpy.sqrt(Elist / theta) / self.symmetry * dE
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates

################################################################################

cdef class KRotor(Rotation):
    """
    A statistical mechanical model of an active K-rotor (a one-dimensional
    rigid rotor). The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor in amu*angstrom^2
    `rotationalConstant`     The rotational constant of the rotor in cm^-1
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.
    
    In the majority of chemical applications, the energies involved in the
    K-rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """
    
    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KRotor object.
        """
        result = 'KRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a KRotor object.
        """
        return (KRotor, (self.inertia, self.symmetry, self.quantum, None))

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

    cpdef double getLevelEnergy(self, int J) except -1:
        """
        Return the energy of level `J` in kJ/mol.
        """
        return getRotationalConstantEnergy(self.inertia.value_si) * J * J
    
    cpdef int getLevelDegeneracy(self, int J) except -1:
        """
        Return the degeneracy of level `J`.
        """
        return 1 if J == 0 else 2
        
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        cdef double Q, B
        if self.quantum:
            Q = schrodinger.getPartitionFunction(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) / self.symmetry
        else:
            B = getRotationalConstantEnergy(self._inertia.value_si)
            Q = sqrt(constants.pi * constants.R * T / B) / self.symmetry
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            Cv = schrodinger.getHeatCapacity(T, self.getLevelEnergy, self.getLevelDegeneracy, 0)
        else:
            Cv = 0.5
        return Cv * constants.R
    
    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            H = schrodinger.getEnthalpy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0)
        else:
            H = 0.5
        return H * constants.R * T
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            S = schrodinger.getEntropy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) - log(self.symmetry)
        else:
            S = log(self.getPartitionFunction(T)) + 0.5
        return S * constants.R
    
    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in J/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double B
        cdef numpy.ndarray sumStates
        if self.quantum:
            sumStates = schrodinger.getSumOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, sumStates0) / self.symmetry
        elif sumStates0 is not None:
            sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
        else:
            B = getRotationalConstantEnergy(self._inertia.value_si)
            sumStates = 2.0 * numpy.sqrt(Elist / B) / self.symmetry
        return sumStates
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in J/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double B, dE
        cdef int r
        cdef numpy.ndarray densStates
        if self.quantum:
            densStates = schrodinger.getDensityOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, densStates0) / self.symmetry
        else:
            densStates = numpy.zeros(Elist.shape[0], dtype=numpy.float)
            B = getRotationalConstantEnergy(self._inertia.value_si)
            dE = Elist[1] - Elist[0]
            for r in range(Elist.shape[0]):
                if Elist[r] == 0: continue
                densStates[r] = 1.0 / sqrt(B * Elist[r]) * dE / self.symmetry
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates

################################################################################

cdef class SphericalTopRotor(Rotation):
    """
    A statistical mechanical model of a three-dimensional rigid rotor with a
    single rotational constant: a spherical top. The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.

    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """
    
    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        SphericalTopRotor object.
        """
        result = 'SphericalTopRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a SphericalTopRotor object.
        """
        return (SphericalTopRotor, (self.inertia, self.symmetry, self.quantum, None))
            
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

    cpdef double getLevelEnergy(self, int J) except -1:
        """
        Return the energy of level `J` in kJ/mol.
        """
        return getRotationalConstantEnergy(self.inertia.value_si) * J * (J + 1)
    
    cpdef int getLevelDegeneracy(self, int J) except -1:
        """
        Return the degeneracy of level `J`.
        """
        return (2 * J + 1) * (2 * J + 1)
    
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        cdef double B, Q, theta
        if self.quantum:
            Q = schrodinger.getPartitionFunction(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) / self.symmetry
        else:
            B = getRotationalConstantEnergy(self._inertia.value_si)
            theta = constants.R * T / B
            Q = sqrt(theta * theta * theta * constants.pi) / self.symmetry
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            Cv = schrodinger.getHeatCapacity(T, self.getLevelEnergy, self.getLevelDegeneracy, 0)
        else:
            Cv = 1.5
        return Cv * constants.R

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            H = schrodinger.getEnthalpy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0)
        else:
            H = 1.5
        return H * constants.R * T
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            S = schrodinger.getEntropy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) - log(self.symmetry)
        else:
            S = log(self.getPartitionFunction(T)) + 1.5
        return S * constants.R

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in J/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double B, theta
        cdef numpy.ndarray sumStates
        if self.quantum:
            sumStates = schrodinger.getSumOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, sumStates0) / self.symmetry
        elif sumStates0 is not None:
            sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
        else:
            if sumStates0 is not None:
                sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
            else:
                B = getRotationalConstantEnergy(self._inertia.value_si)
                theta = B * B * B
                sumStates = 4.0/3.0 * Elist * numpy.sqrt(Elist / theta) / self.symmetry
        return sumStates
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in J/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double theta, inertia, dE
        cdef numpy.ndarray densStates
        if self.quantum:
            densStates = schrodinger.getDensityOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, densStates0) / self.symmetry
        else:
            dE = Elist[1] - Elist[0]
            B = getRotationalConstantEnergy(self._inertia.value_si)
            theta = B * B * B
            densStates = 2.0 * numpy.sqrt(Elist / theta) / self.symmetry * dE
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates
