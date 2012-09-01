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
This module provides the :class:`Conformer` class, used for working with
statistical mechanical models of a molecular system involving multiple degrees
of freedom.
"""

import numpy
import cython

from libc.math cimport log, exp, sqrt

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity

from rmgpy.statmech.translation cimport *
from rmgpy.statmech.rotation cimport *
from rmgpy.statmech.vibration cimport *
from rmgpy.statmech.torsion cimport *

################################################################################

cdef class Conformer:
    """
    A representation of an individual molecular conformation. The attributes 
    are:
    
    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `E0`                The ground-state energy (including zero-point energy) of the conformer
    `modes`             A list of the molecular degrees of freedom
    `spinMultiplicity`  The degeneracy of the electronic ground state
    `opticalIsomers`    The number of optical isomers
    `number`            An array of atomic numbers of each atom in the conformer
    `mass`              An array of masses of each atom in the conformer
    `coordinates`       An array of 3D coordinates of each atom in the conformer
    =================== ========================================================
    
    Note that the `spinMultiplicity` reflects the electronic mode of the
    molecular system.    
    """
    
    def __init__(self, E0=None, modes=None, spinMultiplicity=1, opticalIsomers=1, number=None, mass=None, coordinates=None):
        self.E0 = E0
        self.modes = modes or []
        self.spinMultiplicity = spinMultiplicity
        self.opticalIsomers = opticalIsomers
        self.number = number
        self.mass = mass
        self.coordinates = coordinates

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Conformer object.
        """
        result = 'Conformer(E0={0!r}, modes={0!r}'.format(self.E0, self.modes)
        if self.spinMultiplicity != 1:
            result += ', spinMultiplicity={0:d}'.format(self.spinMultiplicity)
        if self.opticalIsomers != 1:
            result += ', opticalIsomers={0:d}'.format(self.opticalIsomers)
        result += ')'
        return result
    
    def __reduce__(self):
        """
        A helper function used when pickling a Conformer object.
        """
        return (Conformer, (self.E0, self.modes, self.spinMultiplicity, self.opticalIsomers, self.number, self.mass, self.coordinates))
    
    property E0:
        """The ground-state energy (including zero-point energy) of the conformer."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    property number:
        """An array of atomic numbers of each atom in the conformer."""
        def __get__(self):
            return self._number
        def __set__(self, value):
            self._number = quantity.Dimensionless(value)

    property mass:
        """An array of masses of each atom in the conformer."""
        def __get__(self):
            return self._mass
        def __set__(self, value):
            self._mass = quantity.Mass(value)

    property coordinates:
        """An array of 3D coordinates of each atom in the conformer."""
        def __get__(self):
            return self._coordinates
        def __set__(self, value):
            self._coordinates = quantity.Length(value)

    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the partition function :math:`Q(T)` for the system at the
        specified temperature `T` in K.
        """
        cdef double Q = 1.0
        cdef Mode mode
        for mode in self.modes:
            Q *= mode.getPartitionFunction(T)
        return Q * self.spinMultiplicity * self.opticalIsomers

    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the dimensionless heat capacity :math:`C_\\mathrm{v}(T)/R` for 
        the system at the specified temperature `T` in K.
        """
        cdef double Cp = 0.0
        cdef Mode mode
        for mode in self.modes:
            Cp += mode.getHeatCapacity(T)
        return Cp

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the dimensionless enthalpy :math:`H(T)/RT` for the system at the
        specified temperature `T` in K.
        """
        cdef double H = 0.0
        cdef Mode mode
        for mode in self.modes:
            H += mode.getEnthalpy(T)
        return H
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the dimensionless entropy :math:`S(T)/R` for the system at the
        specified temperature `T` in K.
        """
        cdef double S = log(self.spinMultiplicity * self.opticalIsomers) * constants.R
        cdef Mode mode
        for mode in self.modes:
            S += mode.getEntropy(T)
        return S

    cpdef double getFreeEnergy(self, double T) except 100000000:
        """
        Return the dimensionless Gibbs free energy :math:`G(T)/RT` for the
        system at the specified temperature `T` in K.
        """
        return self.getEnthalpy(T) - self.getEntropy(T)
        
    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state.
        """
        cdef numpy.ndarray sumStates = None
        cdef Mode mode
        for mode in self.modes:
            sumStates = mode.getSumOfStates(Elist, sumStates)        
        return sumStates * self.spinMultiplicity * self.opticalIsomers
        
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` above the ground state.
        """
        cdef numpy.ndarray densStates = None
        cdef Mode mode
        for mode in self.modes:
            densStates = mode.getDensityOfStates(Elist, densStates)
        return densStates * self.spinMultiplicity * self.opticalIsomers

    cpdef double getTotalMass(self, atoms=None) except -1:
        """
        Calculate and return the total mass of the atoms in the conformer in 
        kg. If a list `atoms` of atoms is specified, only those atoms will
        be used to calculate the center of mass. Otherwise, all atoms will be
        used.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] mass
        cdef int atom
        mass = self._mass.value_si
        if atoms is None: atoms = range(1, mass.shape[0]+1)
        return sum([mass[atom-1] for atom in atoms])

    cpdef numpy.ndarray getCenterOfMass(self, atoms=None):
        """
        Calculate and return the [three-dimensional] position of the center of
        mass of the conformer in m. If a list `atoms` of atoms is specified, 
        only those atoms will be used to calculate the center of mass. 
        Otherwise, all atoms will be used.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] mass
        cdef numpy.ndarray[numpy.float64_t,ndim=2] coordinates
        cdef double totalMass
        cdef int atom
        
        mass = self._mass.value_si
        coordinates = self._coordinates.value_si
        
        if atoms is None: atoms = range(1, mass.shape[0]+1)
        center = numpy.zeros(3); totalMass = 0.0
        for atom in atoms:
            center += mass[atom-1] * coordinates[atom-1,:]
            totalMass += mass[atom-1]
        center /= totalMass
        
        return center

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef numpy.ndarray getMomentOfInertiaTensor(self):
        """
        Calculate and return the moment of inertia tensor for the conformer 
        in kg*m^2. If the coordinates are not at the center of mass, they are
        temporarily shifted there for the purposes of this calculation.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] mass
        cdef numpy.ndarray[numpy.float64_t,ndim=2] coordinates
        cdef numpy.ndarray[numpy.float64_t,ndim=1] coord, centerOfMass
        cdef numpy.ndarray[numpy.float64_t,ndim=2] I
        cdef double m
        cdef int atom
        
        mass = self._mass.value_si
        coordinates = self._coordinates.value_si
        
        I = numpy.zeros((3,3), numpy.float64)
        centerOfMass = self.getCenterOfMass()
        
        atoms = range(1, mass.shape[0]+1)
        for atom in atoms:
            m = mass[atom-1]
            coord = coordinates[atom-1,:] - centerOfMass
            I[0,0] += m * (coord[1] * coord[1] + coord[2] * coord[2])
            I[1,1] += m * (coord[0] * coord[0] + coord[2] * coord[2])
            I[2,2] += m * (coord[0] * coord[0] + coord[1] * coord[1])
            I[0,1] -= m * coord[0] * coord[1]
            I[0,2] -= m * coord[0] * coord[2]
            I[1,2] -= m * coord[1] * coord[2]
        I[1,0] = I[0,1]
        I[2,0] = I[0,2]
        I[2,1] = I[1,2]
        
        return I
    
    cpdef getPrincipalMomentsOfInertia(self):
        """
        Calculate and return the principal moments of inertia and corresponding 
        principal axes for the conformer. The moments of inertia are in 
        kg*m^2, while the principal axes have unit length.
        """
        I0 = self.getMomentOfInertiaTensor()
        # Since I0 is real and symmetric, diagonalization is always possible
        I, V = numpy.linalg.eig(I0)
        return I, V
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double getInternalReducedMomentOfInertia(self, pivots, top1) except -1:
        """
        Calculate and return the reduced moment of inertia for an internal
        torsional rotation around the axis defined by the two atoms in 
        `pivots`. The list `top1` contains the atoms that should be considered
        as part of the rotating top; this list should contain the pivot atom
        connecting the top to the rest of the molecule.    The procedure used is
        that of Pitzer [1]_, which is described as :math:`I^{(2,3)}` by East
        and Radom [2]_. In this procedure, the molecule is divided into two
        tops: those at either end of the hindered rotor bond. The moment of
        inertia of each top is evaluated using an axis passing through the
        center of mass of both tops. Finally, the reduced moment of inertia is
        evaluated from the moment of inertia of each top via the formula

        .. math:: \\frac{1}{I^{(2,3)}} = \\frac{1}{I_1} + \\frac{1}{I_2}
        
        .. [1] Pitzer, K. S. *J. Chem. Phys.* **14**, p. 239-243 (1946).
        
        .. [2] East, A. L. L. and Radom, L. *J. Chem. Phys.* **106**, p. 6655-6674 (1997).
        
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] mass
        cdef numpy.ndarray[numpy.float64_t,ndim=2] coordinates
        cdef numpy.ndarray[numpy.float64_t,ndim=1] top1CenterOfMass, top2CenterOfMass, axis
        cdef double I1, I2
        cdef int Natoms, atom
        
        mass = self._mass.value_si
        coordinates = self._coordinates.value_si
        
        # The total number of atoms in the geometry
        Natoms = mass.shape[0]

        # Check that exactly one pivot atom is in the specified top
        if pivots[0] not in top1 and pivots[1] not in top1:
            raise ValueError('No pivot atom included in top; you must specify which pivot atom belongs with the specified top.')
        elif pivots[0] in top1 and pivots[1] in top1:
            raise ValueError('Both pivot atoms included in top; you must specify only one pivot atom that belongs with the specified top.')

        # Enumerate atoms in other top
        top2 = [i+1 for i in range(Natoms) if i+1 not in top1]
        
        # Determine centers of mass of each top
        top1CenterOfMass = self.getCenterOfMass(top1)
        top2CenterOfMass = self.getCenterOfMass(top2)
        
        # Determine axis of rotation
        axis = (top1CenterOfMass - top2CenterOfMass)
        axis /= numpy.linalg.norm(axis)
        
        # Determine moments of inertia of each top
        I1 = 0.0
        for atom in top1:
            r1 = coordinates[atom-1,:] - top1CenterOfMass
            r1 -= numpy.dot(r1, axis) * axis
            I1 += mass[atom-1] * numpy.linalg.norm(r1)**2
        I2 = 0.0
        for atom in top2:
            r2 = coordinates[atom-1,:] - top2CenterOfMass
            r2 -= numpy.dot(r2, axis) * axis
            I2 += mass[atom-1] * numpy.linalg.norm(r2)**2
        
        return 1.0 / (1.0 / I1 + 1.0 / I2)
