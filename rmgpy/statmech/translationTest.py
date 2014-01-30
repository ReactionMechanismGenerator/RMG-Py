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
This script contains unit tests of the :mod:`rmgpy.statmech.translation` module.
"""

import unittest
import math
import numpy

from rmgpy.statmech.translation import IdealGasTranslation
import rmgpy.constants as constants

################################################################################

class TestIdealGasTranslation(unittest.TestCase):
    """
    Contains unit tests of the IdealGasTranslation class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.mass = 32.0
        self.quantum = False
        self.mode = IdealGasTranslation(
            mass = (self.mass,"amu"), 
            quantum = self.quantum, 
        )
        
    def test_getPartitionFunction_classical(self):
        """
        Test the IdealGasTranslation.getPartitionFunction() method for a
        classical translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([7.22597e+06, 2.59130e+07, 1.46586e+08, 4.03944e+08, 8.29217e+08])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)
            
    def test_getHeatCapacity_classical(self):
        """
        Test the IdealGasTranslation.getHeatCapacity() method using a classical
        translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([2.5, 2.5, 2.5, 2.5, 2.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, delta=1e-4*Cvexp)
    
    def test_getEnthalpy_classical(self):
        """
        Test the IdealGasTranslation.getEnthalpy() method using a classical
        translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([2.5, 2.5, 2.5, 2.5, 2.5]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEntropy_classical(self):
        """
        Test the IdealGasTranslation.getEntropy() method using a classical
        translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([18.2932, 19.5703, 21.3031, 22.3168, 23.0360]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, delta=1e-4*Sexp)
    
    def test_getSumOfStates_classical(self):
        """
        Test the IdealGasTranslation.getSumOfStates() method using a classical
        translator.
        """
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getDensityOfStates_classical(self):
        """
        Test the IdealGasTranslation.getDensityOfStates() method using a 
        classical translator.
        """
        Elist = numpy.arange(0, 10000*11.96, 1*11.96)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-6*Qexp)

    def test_repr(self):
        """
        Test that an IdealGasTranslation object can be reconstructed from its
        repr() output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.mass.value, mode.mass.value, 6)
        self.assertEqual(self.mode.mass.units, mode.mass.units)
        self.assertEqual(self.mode.quantum, mode.quantum)
        
    def test_pickle(self):
        """
        Test that a IdealGasTranslation object can be pickled and unpickled
        with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode,-1))
        self.assertAlmostEqual(self.mode.mass.value, mode.mass.value, 6)
        self.assertEqual(self.mode.mass.units, mode.mass.units)
        self.assertEqual(self.mode.quantum, mode.quantum)
