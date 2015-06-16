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
import logging

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
        result = 'Conformer(E0={0!r}, modes={1!r}'.format(self.E0, self.modes)
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
            logging.debug('        Calculating Partition Function for ' + mode.__class__.__name__)
            Q *= mode.getPartitionFunction(T)
        return Q * self.spinMultiplicity * self.opticalIsomers

    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the system at the specified
        temperature `T` in K.
        """
        cdef double Cp = 0.0
        cdef Mode mode
        for mode in self.modes:
            Cp += mode.getHeatCapacity(T)
        return Cp

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the system at the specified
        temperature `T` in K.
        """
        cdef double H = 0.0
        cdef Mode mode
        for mode in self.modes:
            H += mode.getEnthalpy(T)
        return H
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the system at the specified
        temperature `T` in K.
        """
        cdef double S = log(self.spinMultiplicity * self.opticalIsomers) * constants.R
        cdef Mode mode
        for mode in self.modes:
            S += mode.getEntropy(T)
        return S

    cpdef double getFreeEnergy(self, double T) except 100000000:
        """
        Return the Gibbs free energy in J/mol for the system at the specified
        temperature `T` in K.
        """
        return self.getEnthalpy(T) - T * self.getEntropy(T)
        
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

    cpdef getNumberDegreesOfFreedom(self):
        """
        Return the number of degrees of freedom in a species object, which should be 3N,
        and raises an exception if it is not.
        """
        cdef int N, Natoms, expectedDegreesOfFreedom
        cdef Mode mode
        N = 0 
        for mode in self.modes:
            if isinstance(mode, HinderedRotor):
                N += 1
            elif hasattr(mode,'frequencies'):
                N += len(mode.frequencies.value) # found the harmonic frequencies
            elif isinstance(mode, IdealGasTranslation):
                # found the translational degrees of freedom
                N += 3
            elif type(mode) == NonlinearRotor: 
                N += 3  
            elif type(mode) == LinearRotor:
                N += 2
            elif type(mode) == KRotor:
                N += 1
            elif type(mode) == SphericalTopRotor:
                N += 3
            else:
                raise TypeError("Mode type {0!r} not supported".format(mode))
        
        if self.mass:
            Natoms =  len(self.mass.value)
            # what the total number of degrees of freedom for the species should be
            expectedDegreesOfFreedom = Natoms * 3
            if N != expectedDegreesOfFreedom:
                raise ValueError('The total degrees of molecular freedom for this species should be {0}'.format(expectedDegreesOfFreedom))          

        return N
    
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef getSymmetricTopRotors(self):
        """
        Return objects representing the external J-rotor and K-rotor under the
        symmetric top approximation. For nonlinear molecules, the J-rotor is
        a 2D rigid rotor with a rotational constant :math:`B` determined as the
        geometric mean of the two most similar rotational constants. The
        K-rotor is a 1D rigid rotor with a rotational constant :math:`A-B`
        determined by the difference between the remaining molecular rotational
        constant and the J-rotor rotational constant.
        """
        cdef double A, B
        cdef numpy.ndarray[numpy.float64_t,ndim=1] Blist

        Jrotor = None; Krotor = None
        for mode in self.modes:
            if isinstance(mode, LinearRotor):
                Jrotor = mode
                Krotor = None
            elif isinstance(mode, NonlinearRotor):
                Blist = numpy.array(sorted(mode.rotationalConstant.value_si))
                assert len(Blist) == 3
                if Blist[1] / Blist[0] < Blist[2] / Blist[1]:
                    B = sqrt(Blist[1] * Blist[0])
                    A = Blist[2]
                else:
                    B = sqrt(Blist[1] * Blist[2])
                    A = Blist[0]
                Jrotor = LinearRotor(rotationalConstant=(B,"cm^-1"), symmetry=1)
                Krotor = KRotor(rotationalConstant=(A-B,"cm^-1"), symmetry=mode.symmetry)

        return Jrotor, Krotor
    
    cpdef list getActiveModes(self, bint activeJRotor=False, bint activeKRotor=True):
        """
        Return a list of the active molecular degrees of freedom of the
        molecular system.
        """
        cdef list modes = []
        
        for mode in self.modes:
            if isinstance(mode, IdealGasTranslation):
                continue
            elif isinstance(mode, LinearRotor):
                if activeJRotor: 
                    modes.append(mode)
            elif isinstance(mode, NonlinearRotor):
                if activeJRotor and activeKRotor: 
                    modes.append(mode)
                elif not activeJRotor and activeKRotor:
                    Jrotor, Krotor = self.getSymmetricTopRotors()
                    modes.append(Krotor)
                else:
                    continue
            else:
                modes.append(mode)
        return modes

################################################################################

cpdef double phi(double beta, int k, double E, logQ) except -10000000:
    """
    Evaluate the value of the objective function used in the method of
    steepest descents to compute the sum and/or density of states from the
    partition function.
    
    :param beta: The value of :math:`\\left( k_\\mathrm{B} T \\right)^{-1}`
                 in mol/J to evaluate the objective function at
    :param k:    0 if computing the density of states, 1 if computing the sum 
                 of states
    :param E:    The energy at which to compute the sum or density of states in
                 J/mol
    :param logQ: A callable object that accepts the current temperature in K as
                 a parameter and returns the natural logarithm of the partition
                 function at that temperature
    :returns:    The value of the objective function to minimize for the
                 method of steepest descents
    """
    cdef double T
    T = 1.0 / (constants.R * beta)
    return logQ(T) - k * log(beta) + beta * E

@cython.boundscheck(False)
@cython.wraparound(False)
def getDensityOfStatesForst(numpy.ndarray[numpy.float64_t,ndim=1] Elist, logQ, int order=1):
    """
    Return the density of states :math:`\\rho(E) \\ dE` and sum of states
    :math:`N(E)` at the specified total energies `Elist` in J/mol above the
    ground state. The parameter `logQ` should be a callable object that accepts
    the current temperature in K as a parameter and returns the natural
    logarithm of the partition function at that temperature. The optional
    `order` parameter indicates the order of the steepest descents 
    approximation to apply (1 or 2); the first-order approximation is smoother,
    faster to compute, and generally accurate enough for most applications.
    """

    cdef numpy.ndarray[numpy.float64_t,ndim=1] densStates, sumStates
    cdef double x, dx, v, E, dE
    cdef int i, k
    
    if order != 1 and order != 2:
        raise ValueError('Invalid value {0} for order parameter; valid values are 1 or 2.'.format(order))
    
    import scipy.optimize
    dE = Elist[1] - Elist[0]
    
    densStates = numpy.zeros_like(Elist)
    sumStates = numpy.zeros_like(Elist)
    
    # Initial guess for first minimization
    x = 1e-4
    
    # Use method of steepest descents to compute sum of states
    k = 1
    
    # Iterate over energies
    for i in range(1, Elist.shape[0]):
        E = Elist[i]
        
        # Find minimum of phi  func x0  arg     xtol  ftol maxi  maxf fullout  disp retall  callback
        try:
            x = scipy.optimize.fmin(phi, x, (k, E, logQ), 1e-8, 1e-8, 100, 1000, False, False, False, None)
        except ValueError:
            break
        x = float(x)
        dx = 1e-2 * x
        
        # Evaluate derivatives needed for steepest descents approximation numerically
        d2fdx2 = (phi(x+dx, k, E, logQ) - 2 * phi(x, k, E, logQ) + phi(x-dx, k, E, logQ)) / (dx*dx)
        
        # Apply first-order steepest descents approximation (accurate to 1-3%, smoother)
        v = phi(x, k, E, logQ)
        if k == 1:
            sumStates[i] = exp(v) / sqrt(2 * constants.pi * d2fdx2)
        else:
            densStates[i] = exp(v) / sqrt(2 * constants.pi * d2fdx2)
        
        if order == 2:
            # Apply second-order steepest descents approximation (more accurate, less smooth)
            d3fdx3 = (phi(x+1.5*dx, k, E, logQ) - 3 * phi(x+0.5*dx, k, E, logQ) + 3 * phi(x-0.5*dx, k, E, logQ) - phi(x-1.5*dx, k, E, logQ)) / (dx**3)
            d4fdx4 = (phi(x+2*dx, k, E, logQ) - 4 * phi(x+dx, k, E, logQ) + 6 * phi(x, k, E, logQ) - 4 * phi(x-dx, k, E, logQ) + phi(x-2*dx, k, E, logQ)) / (dx**4)
            if k == 1:
                sumStates[i] *= 1 + d4fdx4 / 8. / (d2fdx2**2) - 5. * (d3fdx3**2) / 24. / (d2fdx2**3)
            else:
                densStates[i] *= 1 + d4fdx4 / 8. / (d2fdx2**2) - 5. * (d3fdx3**2) / 24. / (d2fdx2**3)
    
        if k == 1:
            d3fdx3 = (phi(x+1.5*dx, k, E, logQ) - 3 * phi(x+0.5*dx, k, E, logQ) + 3 * phi(x-0.5*dx, k, E, logQ) - phi(x-1.5*dx, k, E, logQ)) / (dx**3)
            densStates[i] = sumStates[i] * (x + d3fdx3 / (2. * d2fdx2 * d2fdx2))
        
    return densStates * dE, sumStates
