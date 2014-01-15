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
This script contains unit tests of the :mod:`rmgpy.statmech.vibration` module.
"""

import unittest
import math
import numpy

from rmgpy.statmech.vibration import HarmonicOscillator
import rmgpy.constants as constants

################################################################################

class TestHarmonicOscillator(unittest.TestCase):
    """
    Contains unit tests of the HarmonicOscillator class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.frequencies = numpy.array([500, 1000, 2000])
        self.quantum = True
        self.mode = HarmonicOscillator(
            frequencies = (self.frequencies,"cm^-1"),
            quantum = self.quantum,
        )
        
    def test_getPartitionFunction_classical(self):
        """
        Test the HarmonicOscillator.getPartitionFunction() method for a set of
        classical oscillators.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([0.00906536, 0.04196925, 0.335754, 1.13316978, 2.68603])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getPartitionFunction_quantum(self):
        """
        Test the HarmonicOscillator.getPartitionFunction() method for a set of
        quantum oscillators.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([1.10923, 1.39358, 2.70819, 4.98825, 8.459780])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getHeatCapacity_classical(self):
        """
        Test the HarmonicOscillator.getHeatCapacity() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([3, 3, 3, 3, 3]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
    
    def test_getHeatCapacity_quantum(self):
        """
        Test the HarmonicOscillator.getHeatCapacity() method using a set of
        quantum oscillators.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([0.832004, 1.47271, 2.32513, 2.65024, 2.79124]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
       
    def test_getEnthalpy_classical(self):
        """
        Test the HarmonicOscillator.getEnthalpy() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([3, 3, 3, 3, 3]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEnthalpy_quantum(self):
        """
        Test the HarmonicOscillator.getEnthalpy() method using a set of quantum
        oscillators.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.280395, 0.637310, 1.30209, 1.70542, 1.96142]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)

    def test_getEntropy_classical(self):
        """
        Test the HarmonicOscillator.getEntropy() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([-1.70329, -0.170818, 1.90862, 3.12502, 3.98807]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*abs(Sexp))
        
    def test_getEntropy_quantum(self):
        """
        Test the HarmonicOscillator.getEntropy() method using a set of quantum
        oscillators.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.384065, 0.969182, 2.29837, 3.31251, 4.09675]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)

    def test_getSumOfStates_classical(self):
        """
        Test the HarmonicOscillator.getSumOfStates() method using a set of 
        classical oscillators.
        """
        self.mode.quantum = False
        self.mode.frequencies = ([500, 1000],"cm^-1")
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))
        
    def test_getSumOfStates_quantum(self):
        """
        Test the HarmonicOscillator.getSumOfStates() method using a set of 
        quantum oscillators.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(1, len(Elist)):
            if sumStates[n-1] == 0:
                self.assertTrue(numpy.sum(densStates[0:n]) == 0, '{0} != {1}'.format(numpy.sum(densStates[0:n]), 0))
            else:
                self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getDensityOfStates_classical(self):
        """
        Test the HarmonicOscillator.getDensityOfStates() method using a set of 
        classical oscillators.
        """
        self.mode.quantum = False
        factor = constants.h * constants.c * 100. * constants.Na # cm^-1 to J/mol
        Elist = numpy.arange(0, 10000*factor, 1*factor)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)

    def test_getDensityOfStates_quantum(self):
        """
        Test the HarmonicOscillator.getDensityOfStates() method using a set of 
        quantum oscillators.
        """
        self.mode.quantum = True
        factor = constants.h * constants.c * 100. * constants.Na # cm^-1 to J/mol
        Elist = numpy.arange(0, 10000*factor, 1*factor)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(len(Elist)):
            if densStates[n] != 0:
                # The peaks should occur near a multiple of 500 cm^-1
                E = float(Elist[n]) / factor
                self.assertTrue(E % 500 < 5 or E % 500 > 495)

    def test_repr(self):
        """
        Test that a HarmonicOscillator object can be reconstructed from its
        repr() output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertEqual(self.mode.frequencies.value.shape, mode.frequencies.value.shape)
        for freq0, freq in zip(self.mode.frequencies.value, mode.frequencies.value):
            self.assertAlmostEqual(freq0, freq, 6)
        self.assertEqual(self.mode.frequencies.units, mode.frequencies.units)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
    def test_pickle(self):
        """
        Test that a HarmonicOscillator object can be pickled and unpickled
        with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode,-1))
        self.assertEqual(self.mode.frequencies.value.shape, mode.frequencies.value.shape)
        for freq0, freq in zip(self.mode.frequencies.value, mode.frequencies.value):
            self.assertAlmostEqual(freq0, freq, 6)
        self.assertEqual(self.mode.frequencies.units, mode.frequencies.units)
        self.assertEqual(self.mode.quantum, mode.quantum)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
