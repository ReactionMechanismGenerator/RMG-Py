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
This script contains unit tests of the :mod:`rmgpy.statmech.rotation` module.
"""

import unittest
import math
import numpy

from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
import rmgpy.constants as constants

################################################################################

class TestLinearRotor(unittest.TestCase):
    """
    Contains unit tests of the LinearRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75
        self.symmetry = 2
        self.quantum = False
        self.mode = LinearRotor(
            inertia = (self.inertia,"amu*angstrom^2"), 
            symmetry = self.symmetry, 
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the LinearRotor.rotationalConstant property.
        """
        Bexp = 1.434692
        Bact = self.mode.rotationalConstant.value_si
        self.assertAlmostEqual(Bexp, Bact, 4)
        
    def test_setRotationalConstant(self):
        """
        Test setting the LinearRotor.rotationalConstant property.
        """
        B = self.mode.rotationalConstant
        B.value_si *= 2
        self.mode.rotationalConstant = B
        Iexp = 0.5 * self.inertia
        Iact = self.mode.inertia.value_si * constants.Na * 1e23
        self.assertAlmostEqual(Iexp, Iact, 4)
        
    def test_getLevelEnergy(self):
        """
        Test the LinearRotor.getLevelEnergy() method.
        """
        B = self.mode.rotationalConstant.value_si * constants.h * constants.c * 100.
        B *= constants.Na
        for J in range(0, 100):
            Eexp = B * J * (J + 1)
            Eact = self.mode.getLevelEnergy(J)
            if J == 0:
                self.assertEqual(Eact, 0)
            else:
                self.assertAlmostEqual(Eexp, Eact, delta=1e-4*Eexp)
    
    def test_getLevelDegeneracy(self):
        """
        Test the LinearRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            gexp = 2 * J + 1
            gact = self.mode.getLevelDegeneracy(J)
            self.assertEqual(gexp, gact)
    
    def test_getPartitionFunction_classical(self):
        """
        Test the LinearRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([72.6691, 121.115, 242.230, 363.346, 484.461])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getPartitionFunction_quantum(self):
        """
        Test the LinearRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([72.8360, 121.282, 242.391, 363.512, 484.627])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getHeatCapacity_classical(self):
        """
        Test the LinearRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1, 1, 1, 1, 1]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
            
    def test_getHeatCapacity_quantum(self):
        """
        Test the LinearRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1, 1, 1, 1, 1]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
       
    def test_getEnthalpy_classical(self):
        """
        Test the LinearRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1, 1, 1, 1, 1]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEnthalpy_quantum(self):
        """
        Test the LinearRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.997705, 0.998624, 0.999312, 0.999541, 0.999656]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)

    def test_getEntropy_classical(self):
        """
        Test the LinearRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([5.28592, 5.79674, 6.48989, 6.89535, 7.18304]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)
    
    def test_getEntropy_quantum(self):
        """
        Test the LinearRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([5.28592, 5.79674, 6.48989, 6.89535, 7.18304]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)

    def test_getSumOfStates_classical(self):
        """
        Test the LinearRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 2000*11.96, 1.0*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        sumStates = self.mode.getSumOfStates(Elist)
        for n in range(1, len(Elist)):
            self.assertAlmostEqual(numpy.sum(densStates[0:n]) / sumStates[n], 1.0, 3)

    def test_getSumOfStates_quantum(self):
        """
        Test the LinearRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 4000.*11.96, 2.0*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        sumStates = self.mode.getSumOfStates(Elist)
        for n in range(1, len(Elist)):
            self.assertAlmostEqual(numpy.sum(densStates[0:n+1]) / sumStates[n], 1.0, 3)

    def test_getDensityOfStates_classical(self):
        """
        Test the LinearRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,400,500])
        Elist = numpy.arange(0, 4000.*11.96, 1.0*11.96)
        for T in Tlist:
            densStates = self.mode.getDensityOfStates(Elist)
            Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
            Qexp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_getDensityOfStates_quantum(self):
        """
        Test the LinearRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,400,500])
        Elist = numpy.arange(0, 4000.*11.96, 2.0*11.96)
        for T in Tlist:
            densStates = self.mode.getDensityOfStates(Elist)
            Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
            Qexp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_repr(self):
        """
        Test that a LinearRotor object can be reconstructed from its repr()
        output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
    def test_pickle(self):
        """
        Test that a LinearRotor object can be pickled and unpickled with no
        loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode,-1))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

################################################################################

class TestNonlinearRotor(unittest.TestCase):
    """
    Contains unit tests of the NonlinearRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = numpy.array([3.415, 16.65, 20.07])
        self.symmetry = 4
        self.quantum = False
        self.mode = NonlinearRotor(
            inertia = (self.inertia,"amu*angstrom^2"), 
            symmetry = self.symmetry, 
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the NonlinearRotor.rotationalConstant property.
        """
        Bexp = numpy.array([4.93635, 1.0125, 0.839942])
        Bact = self.mode.rotationalConstant.value_si
        for B0, B in zip(Bexp, Bact):
            self.assertAlmostEqual(B0, B, 4)
        
    def test_setRotationalConstant(self):
        """
        Test setting the NonlinearRotor.rotationalConstant property.
        """
        B = self.mode.rotationalConstant
        B.value_si *= 2
        self.mode.rotationalConstant = B
        Iexp = 0.5 * self.inertia
        Iact = self.mode.inertia.value_si * constants.Na * 1e23
        for I0, I in zip(Iexp, Iact):
            self.assertAlmostEqual(I0, I, 4)
        
    def test_getPartitionFunction_classical(self):
        """
        Test the NonlinearRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([651.162, 1401.08, 3962.84, 7280.21, 11208.6])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getHeatCapacity_classical(self):
        """
        Test the NonlinearRotor.getHeatCapacity() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
    
    def test_getEnthalpy_classical(self):
        """
        Test the NonlinearRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
   
    def test_getEntropy_classical(self):
        """
        Test the NonlinearRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([7.97876, 8.74500, 9.78472, 10.3929, 10.8244]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)
    
    def test_getSumOfStates_classical(self):
        """
        Test the NonlinearRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 1000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getDensityOfStates_classical(self):
        """
        Test the NonlinearRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 1000*11.96, 1*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_repr(self):
        """
        Test that a NonlinearRotor object can be reconstructed from its
        repr() output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertEqual(self.mode.inertia.value.shape, mode.inertia.value.shape)
        for I0, I in zip(self.mode.inertia.value, mode.inertia.value):
            self.assertAlmostEqual(I0, I, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
    def test_pickle(self):
        """
        Test that a NonlinearRotor object can be pickled and unpickled with
        no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode,-1))
        self.assertEqual(self.mode.inertia.value.shape, mode.inertia.value.shape)
        for I0, I in zip(self.mode.inertia.value, mode.inertia.value):
            self.assertAlmostEqual(I0, I, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

################################################################################

class TestKRotor(unittest.TestCase):
    """
    Contains unit tests of the KRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75
        self.symmetry = 2
        self.quantum = False
        self.mode = KRotor(
            inertia = (self.inertia,"amu*angstrom^2"),
            symmetry = self.symmetry, 
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the KRotor.rotationalConstant property.
        """
        Bexp = 1.434692
        Bact = self.mode.rotationalConstant.value_si
        self.assertAlmostEqual(Bexp, Bact, 4)
        
    def test_setRotationalConstant(self):
        """
        Test setting the KRotor.rotationalConstant property.
        """
        B = self.mode.rotationalConstant
        B.value_si *= 2
        self.mode.rotationalConstant = B
        Iexp = 0.5 * self.inertia
        Iact = self.mode.inertia.value_si * constants.Na * 1e23
        self.assertAlmostEqual(Iexp, Iact, 4)
        
    def test_getLevelEnergy(self):
        """
        Test the KRotor.getLevelEnergy() method.
        """
        B = self.mode.rotationalConstant.value_si * constants.h * constants.c * 100.
        B *= constants.Na
        for J in range(0, 100):
            Eexp = float(B * J * J)
            Eact = float(self.mode.getLevelEnergy(J))
            if J == 0:
                self.assertEqual(Eact, 0)
            else:
                self.assertAlmostEqual(Eexp, Eact, delta=1e-4*Eexp)
    
    def test_getLevelDegeneracy(self):
        """
        Test the KRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            gexp = 1 if J == 0 else 2
            gact = self.mode.getLevelDegeneracy(J)
            self.assertEqual(gexp, gact, '{0} != {1}'.format(gact, gexp))
    
    def test_getPartitionFunction_classical(self):
        """
        Test the KRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([10.6839, 13.7929, 19.5060, 23.8899, 27.5857])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getPartitionFunction_quantum(self):
        """
        Test the KRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([10.6839, 13.7929, 19.5060, 23.8899, 27.5857])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getHeatCapacity_classical(self):
        """
        Test the KRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
    
    def test_getHeatCapacity_quantum(self):
        """
        Test the KRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
       
    def test_getEnthalpy_classical(self):
        """
        Test the KRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEnthalpy_quantum(self):
        """
        Test the KRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)

    def test_getEntropy_classical(self):
        """
        Test the KRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([2.86874, 3.12415, 3.47072, 3.67346, 3.81730]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)
    
    def test_getEntropy_quantum(self):
        """
        Test the KRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([2.86874, 3.12415, 3.47072, 3.67346, 3.81730]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)

    def test_getSumOfStates_classical(self):
        """
        Test the KRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 1000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.75 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.3333, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))

    def test_getSumOfStates_quantum(self):
        """
        Test the KRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 1000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))

    def test_getDensityOfStates_classical(self):
        """
        Test the KRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 3000*11.96, 0.05*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 500
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_getDensityOfStates_quantum(self):
        """
        Test the KRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 4000*11.96, 2*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 500
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_repr(self):
        """
        Test that a KRotor object can be reconstructed from its repr() output
        with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
    def test_pickle(self):
        """
        Test that a KRotor object can be pickled and unpickled with no loss
        of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode,-1))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

################################################################################

class TestSphericalTopRotor(unittest.TestCase):
    """
    Contains unit tests of the SphericalTopRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75
        self.symmetry = 2
        self.quantum = False
        self.mode = SphericalTopRotor(
            inertia = (self.inertia,"amu*angstrom^2"), 
            symmetry = self.symmetry, 
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the SphericalTopRotor.rotationalConstant property.
        """
        Bexp = 1.434692
        Bact = self.mode.rotationalConstant.value_si
        self.assertAlmostEqual(Bexp, Bact, 4)
        
    def test_setRotationalConstant(self):
        """
        Test setting the SphericalTopRotor.rotationalConstant property.
        """
        B = self.mode.rotationalConstant
        B.value_si *= 2
        self.mode.rotationalConstant = B
        Iexp = 0.5 * self.inertia
        Iact = self.mode.inertia.value_si * constants.Na * 1e23
        self.assertAlmostEqual(Iexp, Iact, 4)
        
    def test_getLevelEnergy(self):
        """
        Test the SphericalTopRotor.getLevelEnergy() method.
        """
        B = self.mode.rotationalConstant.value_si * constants.h * constants.c * 100.
        B *= constants.Na
        for J in range(0, 100):
            Eexp = B * J * (J + 1)
            Eact = self.mode.getLevelEnergy(J)
            if J == 0:
                self.assertEqual(Eact, 0)
            else:
                self.assertAlmostEqual(Eexp, Eact, delta=1e-4*Eexp)
    
    def test_getLevelDegeneracy(self):
        """
        Test the SphericalTopRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            gexp = (2 * J + 1)**2
            gact = self.mode.getLevelDegeneracy(J)
            self.assertEqual(gexp, gact, '{0} != {1}'.format(gact, gexp))
    
    def test_getPartitionFunction_classical(self):
        """
        Test the SphericalTopRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([1552.74, 3340.97, 9449.69, 17360.2, 26727.8])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getPartitionFunction_quantum(self):
        """
        Test the SphericalTopRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([1555.42, 3344.42, 9454.57, 17366.2, 26734.7])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getHeatCapacity_classical(self):
        """
        Test the SphericalTopRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
    
    def test_getHeatCapacity_quantum(self):
        """
        Test the SphericalTopRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
       
    def test_getEnthalpy_classical(self):
        """
        Test the SphericalTopRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEnthalpy_quantum(self):
        """
        Test the SphericalTopRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.49828, 1.49897, 1.49948, 1.49966, 1.49974]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)

    def test_getEntropy_classical(self):
        """
        Test the SphericalTopRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([8.84778, 9.61402, 10.6537, 11.2619, 11.6935]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)
    
    def test_getEntropy_quantum(self):
        """
        Test the SphericalTopRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([8.84778, 9.61402, 10.6537, 11.2619, 11.6935]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)

    def test_getSumOfStates_classical(self):
        """
        Test the SphericalTopRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 2000*11.96, 1.0*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        sumStates = self.mode.getSumOfStates(Elist)
        for n in range(20, len(Elist)):
            self.assertAlmostEqual(numpy.sum(densStates[0:n+1]) / sumStates[n], 1.0, 1)

    def test_getSumOfStates_quantum(self):
        """
        Test the SphericalTopRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 2000*11.96, 1.0*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        sumStates = self.mode.getSumOfStates(Elist)
        for n in range(1, len(Elist)):
            self.assertAlmostEqual(numpy.sum(densStates[0:n+1]) / sumStates[n], 1.0, 3)

    def test_getDensityOfStates_classical(self):
        """
        Test the SphericalTopRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,400,500])
        Elist = numpy.arange(0, 2000*11.96, 1.0*11.96)
        for T in Tlist:
            densStates = self.mode.getDensityOfStates(Elist)
            Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
            Qexp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_getDensityOfStates_quantum(self):
        """
        Test the SphericalTopRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,400,500])
        Elist = numpy.arange(0, 4000*11.96, 2.0*11.96)
        for T in Tlist:
            densStates = self.mode.getDensityOfStates(Elist)
            Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
            Qexp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-2*Qexp)

    def test_repr(self):
        """
        Test that a SphericalTopRotor object can be reconstructed from its
        repr() output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
    def test_pickle(self):
        """
        Test that a SphericalTopRotor object can be pickled and unpickled
        with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode,-1))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
