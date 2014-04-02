#!/usr/bin/env python
# encoding: utf-8

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
This script contains unit tests of the :mod:`rmgpy.statmech.torsion` module.
"""

import unittest
import math
import numpy

from rmgpy.statmech.torsion import HinderedRotor
import rmgpy.constants as constants

################################################################################

class TestHinderedRotor(unittest.TestCase):
    """
    Contains unit tests of the HinderedRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 1.56764
        self.symmetry = 3
        self.barrier = 11.373 
        self.quantum = True
        self.mode = HinderedRotor(
            inertia = (self.inertia,"amu*angstrom^2"), 
            symmetry = self.symmetry,
            barrier = (self.barrier,"kJ/mol"),
            fourier = ([ [4.58375, 0.841648, -5702.71, 6.02657, 4.7446], [0.726951, -0.677255, 0.207032, 0.553307, -0.503303] ],"J/mol"),
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the HinderedRotor.rotationalConstant property.
        """
        Bexp = 10.7535
        Bact = self.mode.rotationalConstant.value_si
        self.assertAlmostEqual(Bexp, Bact, 4)
        
    def test_setRotationalConstant(self):
        """
        Test setting the HinderedRotor.rotationalConstant property.
        """
        B = self.mode.rotationalConstant
        B.value_si *= 2
        self.mode.rotationalConstant = B
        Iexp = 0.5 * self.inertia
        Iact = self.mode.inertia.value_si * constants.Na * 1e23
        self.assertAlmostEqual(Iexp, Iact, 4)
    
    def test_getPotential_cosine(self):
        """
        Test the HinderedRotor.getPotential() method for a cosine potential.
        """
        self.mode.fourier = None
        phi = numpy.arange(0.0, 2 * constants.pi + 0.0001, constants.pi / 24.)
        V = numpy.zeros_like(phi)
        for i in range(phi.shape[0]):
            V[i] = self.mode.getPotential(phi[i])
    
    def test_getPotential_fourier(self):
        """
        Test the HinderedRotor.getPotential() method for a Fourier series
        potential.
        """
        phi = numpy.arange(0.0, 2 * constants.pi + 0.0001, constants.pi / 24.)
        V = numpy.zeros_like(phi)
        for i in range(phi.shape[0]):
            V[i] = self.mode.getPotential(phi[i])

    def test_getPartitionFunction_classical_cosine(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a cosine
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([0.741953, 1.30465, 2.68553, 3.88146, 4.91235])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getPartitionFunction_classical_fourier(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a Fourier
        series potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([0.745526, 1.30751, 2.68722, 3.88258, 4.91315])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getPartitionFunction_quantum_cosine(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a cosine
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([1.39947, 1.94793, 3.30171, 4.45856, 5.45188])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getPartitionFunction_quantum_fourier(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a Fourier
        series potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([1.39364, 1.94182, 3.29509, 4.45205, 5.44563])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=5e-4*Qexp)

    def test_getHeatCapacity_classical_cosine(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a cosine
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.01741, 0.951141, 0.681919, 0.589263, 0.552028]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
    
    def test_getHeatCapacity_classical_fourier(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.17682, 1.01369, 0.698588, 0.596797, 0.556293]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
        
    def test_getHeatCapacity_quantum_cosine(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a cosine
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.01271, 0.945341, 0.684451, 0.591949, 0.554087]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
     
    def test_getHeatCapacity_quantum_fourier(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.01263, 0.946618, 0.685345, 0.592427, 0.554374]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-3*Cvexp)
 
    def test_getEnthalpy_classical_cosine(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.09556, 1.09949, 0.962738, 0.854617, 0.784333]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
     
    def test_getEnthalpy_classical_fourier(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a Fourier series 
        potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.08882, 1.09584, 0.961543, 0.854054, 0.784009]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)

    def test_getEnthalpy_quantum_cosine(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.545814, 0.727200, 0.760918, 0.717496, 0.680767]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEnthalpy_quantum_fourier(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a Fourier series 
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.548251, 0.728974, 0.762396, 0.718702, 0.681764]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-3*Hexp)

    def test_getEntropy_classical_cosine(self):
        """
        Test the HinderedRotor.getEntropy() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.797089, 1.36543, 1.95062, 2.21083, 2.37608]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)

    def test_getEntropy_classical_fourier(self):
        """
        Test the HinderedRotor.getEntropy() method using a Fourier series 
        potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.795154, 1.36396, 1.95005, 2.21055, 2.37592]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)

    def test_getEntropy_quantum_cosine(self):
        """
        Test the HinderedRotor.getEntropy() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.881906, 1.39397, 1.95536, 2.21232, 2.37673]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)
    
    def test_getEntropy_quantum_fourier(self):
        """
        Test the HinderedRotor.getEntropy() method using a Fourier series 
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.880170, 1.39260, 1.95483, 2.21207, 2.37658]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-3*Sexp)

    def test_getSumOfStates_classical_cosine(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getSumOfStates_classical_fourier(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        try:
            sumStates = self.mode.getSumOfStates(Elist)
            self.fail('NotImplementedError not raised by HinderedRotor.getSumOfStates()')
        except NotImplementedError:
            pass
        
    def test_getSumOfStates_quantum_cosine(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getSumOfStates_quantum_fourier(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getDensityOfStates_classical_cosine(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a classical
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_getDensityOfStates_classical_fourier(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a Fourier 
        series potential in the classical limit.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        try:
            densStates = self.mode.getDensityOfStates(Elist)
            self.fail('NotImplementedError not raised by HinderedRotor.getDensityOfStates()')
        except NotImplementedError:
            pass
        
    def test_getDensityOfStates_quantum_cosine(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a classical
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_getDensityOfStates_quantum_fourier(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a Fourier 
        series potential in the quantum limit.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_repr(self):
        """
        Test that a HinderedRotor object can be reconstructed from its repr()
        output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units, 6)
        self.assertEqual(self.mode.fourier.value.shape, mode.fourier.value.shape)
        for A0, A in zip(self.mode.fourier.value.flat, mode.fourier.value.flat):
            self.assertAlmostEqual(A0 / A, 1.0, 6)
        self.assertEqual(self.mode.fourier.units, mode.fourier.units)
        self.assertAlmostEqual(self.mode.barrier.value, mode.barrier.value, 6)
        self.assertEqual(self.mode.barrier.units, mode.barrier.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
    def test_pickle(self):
        """
        Test that a HinderedRotor object can be pickled and unpickled with no
        loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode,-1))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units, 6)
        self.assertEqual(self.mode.fourier.value.shape, mode.fourier.value.shape)
        for A0, A in zip(self.mode.fourier.value.flat, mode.fourier.value.flat):
            self.assertAlmostEqual(A0 / A, 1.0, 6)
        self.assertEqual(self.mode.fourier.units, mode.fourier.units)
        self.assertAlmostEqual(self.mode.barrier.value, mode.barrier.value, 6)
        self.assertEqual(self.mode.barrier.units, mode.barrier.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
