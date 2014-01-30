#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest

from rmgpy.cantherm.geometry import Geometry
import rmgpy.constants as constants

################################################################################

class GeometryTest(unittest.TestCase):

    def testEthaneInternalReducedMomentOfInertia(self):
        """
        Uses an optimum geometry for ethane (CC) to test that the
        proper moments of inertia for its internal hindered rotor is 
        calculated.
        """
        
        # Masses should be in kg/mol
        mass = numpy.array([12.0, 1.0, 1.0, 1.0, 12.0, 1.0, 1.0, 1.0], numpy.float64) * 0.001

        # Atomic numbers
        number = numpy.array([6, 1, 1, 1, 6, 1, 1, 1], numpy.int)

        # Coordinates should be in m
        position = numpy.zeros((8,3), numpy.float64)
        position[0,:] = numpy.array([ 0.001294,  0.002015,  0.000152]) * 1e-10
        position[1,:] = numpy.array([ 0.397758,  0.629904, -0.805418]) * 1e-10
        position[2,:] = numpy.array([-0.646436,  0.631287,  0.620549]) * 1e-10
        position[3,:] = numpy.array([ 0.847832, -0.312615,  0.620435]) * 1e-10
        position[4,:] = numpy.array([-0.760734, -1.204707, -0.557036]) * 1e-10
        position[5,:] = numpy.array([-1.15728 , -1.832718,  0.248402]) * 1e-10
        position[6,:] = numpy.array([-1.607276, -0.890277, -1.177452]) * 1e-10
        position[7,:] = numpy.array([-0.11271 , -1.833701, -1.177357]) * 1e-10
        
        geometry = Geometry(position, number, mass)
        
        pivots = [0, 4]
        top = [0, 1, 2, 3]
        
        # Returned moment of inertia is in kg*m^2; convert to amu*A^2
        I = geometry.getInternalReducedMomentOfInertia(pivots, top) * 1e23 * constants.Na
        self.assertAlmostEqual(I / 1.5595197928, 1.0, 2)
        
    def testButanolInternalReducedMomentOfInertia(self):
        """
        Uses an optimum geometry for s-butanol (CCC(O)C) to test that the
        proper moments of inertia for its internal hindered rotors are 
        calculated.
        """
        
        # Masses should be in kg/mol
        mass = numpy.array([12.0107, 1.00794, 1.00794, 1.00794, 12.0107, 1.00794, 1.00794, 12.0107, 1.00794, 12.0107, 1.00794, 1.00794, 1.00794, 15.9994, 1.00794], numpy.float64) * 0.001
        
        # Atomic numbers
        number = numpy.array([6, 1, 1, 1, 6, 1, 1, 6, 1, 6, 1, 1, 1, 8, 1], numpy.int)

        # Coordinates should be in m
        position = numpy.zeros((15,3), numpy.float64)
        position[0,:] = numpy.array([-2.066968, -0.048470, -0.104326]) * 1e-10
        position[1,:] = numpy.array([-2.078133,  1.009166,  0.165745]) * 1e-10
        position[2,:] = numpy.array([-2.241129, -0.116565, -1.182661]) * 1e-10
        position[3,:] = numpy.array([-2.901122, -0.543098,  0.400010]) * 1e-10
        position[4,:] = numpy.array([-0.729030, -0.686020,  0.276105]) * 1e-10
        position[5,:] = numpy.array([-0.614195, -0.690327,  1.369198]) * 1e-10
        position[6,:] = numpy.array([-0.710268, -1.736876, -0.035668]) * 1e-10
        position[7,:] = numpy.array([ 0.482521,  0.031583, -0.332519]) * 1e-10
        position[8,:] = numpy.array([ 0.358535,  0.069368, -1.420087]) * 1e-10
        position[9,:] = numpy.array([ 1.803404, -0.663583, -0.006474]) * 1e-10
        position[10,:] = numpy.array([ 1.825001, -1.684006, -0.400007]) * 1e-10
        position[11,:] = numpy.array([ 2.638619, -0.106886, -0.436450]) * 1e-10
        position[12,:] = numpy.array([ 1.953652, -0.720890,  1.077945]) * 1e-10
        position[13,:] = numpy.array([ 0.521504,  1.410171,  0.056819]) * 1e-10
        position[14,:] = numpy.array([ 0.657443,  1.437685,  1.010704]) * 1e-10
        
        geometry = Geometry(position, number, mass)
        
        pivots = [0, 4]
        top = [0, 1, 2, 3]
        I = geometry.getInternalReducedMomentOfInertia(pivots, top) * 1e23 * constants.Na
        self.assertAlmostEqual(I / 2.73090431938, 1.0, 3)
        
        pivots = [4, 7]
        top = [4, 5, 6, 0, 1, 2, 3]
        I = geometry.getInternalReducedMomentOfInertia(pivots, top) * 1e23 * constants.Na
        self.assertAlmostEqual(I / 12.1318136515, 1.0, 3)
        
        pivots = [13, 7]
        top = [13, 14]
        I = geometry.getInternalReducedMomentOfInertia(pivots, top) * 1e23 * constants.Na
        self.assertAlmostEqual(I / 0.853678578741, 1.0, 3)
        
        pivots = [9, 7]
        top = [9, 10, 11, 12]
        I = geometry.getInternalReducedMomentOfInertia(pivots, top) * 1e23 * constants.Na
        self.assertAlmostEqual(I / 2.97944840397, 1.0, 3)

    def testPickle(self):
        """
        Test that a Geometry object can be successfully pickled and unpickled
        with no loss of information.
        """

        # Masses should be in kg/mol
        mass = numpy.array([12.0, 1.0, 1.0, 1.0, 12.0, 1.0, 1.0, 1.0], numpy.float64) * 0.001
        # Atomic numbers
        number = numpy.array([6, 1, 1, 1, 6, 1, 1, 1], numpy.int)
        # Coordinates should be in m
        position = numpy.zeros((8,3), numpy.float64)
        position[0,:] = numpy.array([ 0.001294,  0.002015,  0.000152]) * 1e-10
        position[1,:] = numpy.array([ 0.397758,  0.629904, -0.805418]) * 1e-10
        position[2,:] = numpy.array([-0.646436,  0.631287,  0.620549]) * 1e-10
        position[3,:] = numpy.array([ 0.847832, -0.312615,  0.620435]) * 1e-10
        position[4,:] = numpy.array([-0.760734, -1.204707, -0.557036]) * 1e-10
        position[5,:] = numpy.array([-1.15728 , -1.832718,  0.248402]) * 1e-10
        position[6,:] = numpy.array([-1.607276, -0.890277, -1.177452]) * 1e-10
        position[7,:] = numpy.array([-0.11271 , -1.833701, -1.177357]) * 1e-10

        g0 = Geometry(position, number, mass)

        import cPickle
        g = cPickle.loads(cPickle.dumps(g0,-1))

        Natoms = len(g.number)
        self.assertEqual(len(g0.number), len(g.number))
        for i in range(Natoms):
            for j in range(3):
                self.assertEqual(g0.coordinates[i,j], g.coordinates[i,j])
            self.assertEqual(g0.number[i], g.number[i])
            self.assertEqual(g0.mass[i], g.mass[i])

    def testOutput(self):
        """
        Test that a Geometry object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        # Masses should be in kg/mol
        mass = numpy.array([12.0, 1.0, 1.0, 1.0, 12.0, 1.0, 1.0, 1.0], numpy.float64) * 0.001
        # Atomic numbers
        number = numpy.array([6, 1, 1, 1, 6, 1, 1, 1], numpy.int)
        # Coordinates should be in m
        position = numpy.zeros((8,3), numpy.float64)
        position[0,:] = numpy.array([ 0.001294,  0.002015,  0.000152]) * 1e-10
        position[1,:] = numpy.array([ 0.397758,  0.629904, -0.805418]) * 1e-10
        position[2,:] = numpy.array([-0.646436,  0.631287,  0.620549]) * 1e-10
        position[3,:] = numpy.array([ 0.847832, -0.312615,  0.620435]) * 1e-10
        position[4,:] = numpy.array([-0.760734, -1.204707, -0.557036]) * 1e-10
        position[5,:] = numpy.array([-1.15728 , -1.832718,  0.248402]) * 1e-10
        position[6,:] = numpy.array([-1.607276, -0.890277, -1.177452]) * 1e-10
        position[7,:] = numpy.array([-0.11271 , -1.833701, -1.177357]) * 1e-10

        g0 = Geometry(position, number, mass)
        exec('g = %r' % g0)
        
        Natoms = len(g.number)
        self.assertEqual(len(g0.number), len(g.number))
        for i in range(Natoms):
            for j in range(3):
                self.assertAlmostEqual(g0.coordinates[i,j], g.coordinates[i,j], 6)
            self.assertEqual(g0.number[i], g.number[i])
            self.assertAlmostEqual(g0.mass[i], g.mass[i], 6)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
