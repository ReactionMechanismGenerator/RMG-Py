#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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
Contains the :class:`Geometry` class for working with the three-dimensional
geometry of molecules and evaluating properties based on the geometry
information, e.g. moments of inertia.
"""

import numpy
import cython

import rmgpy.constants as constants
from rmgpy.quantity import Quantity

################################################################################

class GeometryError(Exception):
    """
    An exception class for errors that occur while working with molecular
    geometries. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

################################################################################

class Geometry:
    """
    The three-dimensional geometry of a molecular configuration. The attributes
    are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `coordinates`   :class:`numpy.ndarray`  An N x 3 array containing the 3D coordinates of each atom
    `number`        :class:`numpy.ndarray`  An array containing the integer atomic number of each atom
    `mass`          :class:`numpy.ndarray`  An array containing the atomic mass in kg/mol of each atom
    =============== ======================= ====================================

    The integer index of each atom is consistent across all three attributes.
    """
    
    def __init__(self, coordinates, number, mass):
        self.coordinates = Quantity(coordinates).value_si
        self.number = numpy.array(number)
        self.mass = Quantity(mass).value_si
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        coordinates = '(['
        for i in range(self.coordinates.shape[0]):
            if i > 0: coordinates += ', '
            coordinates += '[{0}]'.format(','.join(['{0:g}'.format(self.coordinates[i,j]) for j in range(self.coordinates.shape[1])]))
        coordinates += '],"m")'
        number = '[{0}]'.format(','.join(['{0:d}'.format(n) for n in self.number]))
        mass = '([{0}],"g/mol")'.format(','.join(['{0:g}'.format(m * 1000.) for m in self.mass]))
        return 'Geometry(coordinates={0}, number={1}, mass={2})'.format(coordinates, number, mass)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Geometry, (self.coordinates, self.number, self.mass))

    def getTotalMass(self, atoms=None):
        """
        Calculate and return the total mass of the atoms in the geometry in 
        kg/mol. If a list `atoms` of atoms is specified, only those atoms will
        be used to calculate the center of mass. Otherwise, all atoms will be
        used.
        """
        if atoms is None: atoms = range(len(self.mass))
        return sum([self.mass[atom] for atom in atoms])

    def getCenterOfMass(self, atoms=None):
        """
        Calculate and return the [three-dimensional] position of the center of
        mass of the current geometry. If a list `atoms` of atoms is specified,
        only those atoms will be used to calculate the center of mass. 
        Otherwise, all atoms will be used.
        """

        cython.declare(center=numpy.ndarray, mass=cython.double, atom=cython.int)

        if atoms is None: atoms = range(len(self.mass))
        center = numpy.zeros(3, numpy.float64); mass = 0.0
        for atom in atoms:
            center += self.mass[atom] * self.coordinates[atom]
            mass += self.mass[atom]
        center /= mass
        return center

    def getMomentOfInertiaTensor(self):
        """
        Calculate and return the moment of inertia tensor for the current 
        geometry in kg*m^2. If the coordinates are not at the center of mass,
        they are temporarily shifted there for the purposes of this calculation.
        """
        
        cython.declare(I=numpy.ndarray, mass=cython.double, atom=cython.int)
        cython.declare(coord0=numpy.ndarray, coord=numpy.ndarray, centerOfMass=numpy.ndarray)

        I = numpy.zeros((3,3), numpy.float64)
        centerOfMass = self.getCenterOfMass()
        for atom, coord0 in enumerate(self.coordinates):
            mass = self.mass[atom] / constants.Na
            coord = coord0 - centerOfMass
            I[0,0] += mass * (coord[1] * coord[1] + coord[2] * coord[2])
            I[1,1] += mass * (coord[0] * coord[0] + coord[2] * coord[2])
            I[2,2] += mass * (coord[0] * coord[0] + coord[1] * coord[1])
            I[0,1] -= mass * coord[0] * coord[1]
            I[0,2] -= mass * coord[0] * coord[2]
            I[1,2] -= mass * coord[1] * coord[2]
        I[1,0] = I[0,1]
        I[2,0] = I[0,2]
        I[2,1] = I[1,2]
        
        return I
    
    def getPrincipalMomentsOfInertia(self):
        """
        Calculate and return the principal moments of inertia and corresponding 
        principal axes for the current geometry. The moments of inertia are in
        kg*m^2, while the principal axes have unit length.
        """
        I0 = self.getMomentOfInertiaTensor()
        # Since I0 is real and symmetric, diagonalization is always possible
        I, V = numpy.linalg.eig(I0)
        return I, V
    
    def getInternalReducedMomentOfInertia(self, pivots, top1):
        """
        Calculate and return the reduced moment of inertia for an internal
        torsional rotation around the axis defined by the two atoms in 
        `pivots`. The list `top1` contains the atoms that should be considered
        as part of the rotating top; this list should contain the pivot atom
        connecting the top to the rest of the molecule.	The procedure used is
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

        cython.declare(Natoms=cython.int, top2=list, top1CenterOfMass=numpy.ndarray, top2CenterOfMass=numpy.ndarray)
        cython.declare(axis=numpy.ndarray, I1=cython.double, I2=cython.double, atom=cython.int, i=cython.int)

        # The total number of atoms in the geometry
        Natoms = len(self.mass)

        # Check that exactly one pivot atom is in the specified top
        if pivots[0] not in top1 and pivots[1] not in top1:
            raise GeometryError('No pivot atom included in top; you must specify which pivot atom belongs with the specified top.')
        elif pivots[0] in top1 and pivots[1] in top1:
            raise GeometryError('Both pivot atoms included in top; you must specify only one pivot atom that belongs with the specified top.')

        # Determine atoms in other top
        top2 = []
        for i in range(Natoms):
            if i not in top1: top2.append(i)
        
        # Determine centers of mass of each top
        top1CenterOfMass = self.getCenterOfMass(top1)
        top2CenterOfMass = self.getCenterOfMass(top2)
        
        # Determine axis of rotation
        axis = (top1CenterOfMass - top2CenterOfMass)
        axis /= numpy.linalg.norm(axis)
        
        # Determine moments of inertia of each top
        I1 = 0.0
        for atom in top1:
            r1 = self.coordinates[atom,:] - top1CenterOfMass
            r1 -= numpy.dot(r1, axis) * axis
            I1 += self.mass[atom] / constants.Na * numpy.linalg.norm(r1)**2
        I2 = 0.0
        for atom in top2:
            r2 = self.coordinates[atom,:] - top2CenterOfMass
            r2 -= numpy.dot(r2, axis) * axis
            I2 += self.mass[atom] / constants.Na * numpy.linalg.norm(r2)**2
        
        return 1.0 / (1.0 / I1 + 1.0 / I2)
