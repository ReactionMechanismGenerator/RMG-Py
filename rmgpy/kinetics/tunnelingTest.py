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
This script contains unit tests of the :mod:`rmgpy.kinetics.tunneling` module.
"""

import unittest
import numpy

from rmgpy.kinetics.tunneling import Wigner, Eckart

################################################################################

class TestWigner(unittest.TestCase):
    """
    Contains unit tests of the :class:`Wigner` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.frequency = -2017.96
        self.tunneling = Wigner(
            frequency = (self.frequency,"cm^-1"),
        )
    
    def test_frequency(self):
        """
        Test that the Wigner frequency property was properly set.
        """
        self.assertAlmostEqual(self.tunneling.frequency.value_si, self.frequency, 4)

    def test_calculateTunnelingFactor(self):
        """
        Test the Wigner.calculateTunnelingFactor() method.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        kexplist = numpy.array([4.90263, 2.40495, 1.35124, 1.15611, 1.08781])
        for T, kexp in zip(Tlist, kexplist):
            kact = self.tunneling.calculateTunnelingFactor(T)
            self.assertAlmostEqual(kexp, kact, 4)

    def test_pickle(self):
        """
        Test that a Wigner object can be successfully pickled and unpickled
        with no loss of information.
        """
        import cPickle
        tunneling = cPickle.loads(cPickle.dumps(self.tunneling,-1))
        self.assertAlmostEqual(self.tunneling.frequency.value, tunneling.frequency.value, 2)
        self.assertEqual(self.tunneling.frequency.units, tunneling.frequency.units)

    def test_repr(self):
        """
        Test that a Wigner object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        tunneling = None
        exec('tunneling = {0!r}'.format(self.tunneling))
        self.assertAlmostEqual(self.tunneling.frequency.value, tunneling.frequency.value, 2)
        self.assertEqual(self.tunneling.frequency.units, tunneling.frequency.units)

################################################################################

class TestEckart(unittest.TestCase):
    """
    Contains unit tests of the :class:`Eckart` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.frequency = -2017.96
        self.E0_reac = -295.563
        self.E0_TS = -12.7411
        self.E0_prod = (-10.2664) + (-253.48)
        self.tunneling = Eckart(
            frequency = (self.frequency,"cm^-1"),
            E0_reac = (self.E0_reac,"kJ/mol"),
            E0_TS = (self.E0_TS,"kJ/mol"),
            E0_prod = (self.E0_prod,"kJ/mol"),
        )

    def test_frequency(self):
        """
        Test that the Eckart frequency property was properly set.
        """
        self.assertAlmostEqual(self.tunneling.frequency.value_si, self.frequency, 4)

    def test_E0_reac(self):
        """
        Test that the Eckart E0_reac property was properly set.
        """
        self.assertAlmostEqual(self.tunneling.E0_reac.value_si*0.001, self.E0_reac, 4)
        
    def test_E0_TS(self):
        """
        Test that the Eckart E0_TS property was properly set.
        """
        self.assertAlmostEqual(self.tunneling.E0_TS.value_si*0.001, self.E0_TS, 4)
        
    def test_E0_prod(self):
        """
        Test that the Eckart E0_prod property was properly set.
        """
        self.assertAlmostEqual(self.tunneling.E0_prod.value_si*0.001, self.E0_prod, 4)
    
    def test_calculateTunnelingFactor(self):
        """
        Test the Eckart.calculateTunnelingFactor() method.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        kexplist = numpy.array([1623051., 7.69349, 1.46551, 1.18111, 1.09858])
        for T, kexp in zip(Tlist, kexplist):
            kact = self.tunneling.calculateTunnelingFactor(T)
            self.assertAlmostEqual(kexp, kact, delta=1e-3*kexp)

    def test_pickle(self):
        """
        Test that an Eckart object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        tunneling = cPickle.loads(cPickle.dumps(self.tunneling,-1))
        self.assertAlmostEqual(self.tunneling.frequency.value, tunneling.frequency.value, 2)
        self.assertEqual(self.tunneling.frequency.units, tunneling.frequency.units)
        self.assertAlmostEqual(self.tunneling.E0_reac.value, tunneling.E0_reac.value, 3)
        self.assertEqual(self.tunneling.E0_reac.units, tunneling.E0_reac.units)
        self.assertAlmostEqual(self.tunneling.E0_TS.value, tunneling.E0_TS.value, 3)
        self.assertEqual(self.tunneling.E0_TS.units, tunneling.E0_TS.units)
        self.assertAlmostEqual(self.tunneling.E0_prod.value, tunneling.E0_prod.value, 3)
        self.assertEqual(self.tunneling.E0_prod.units, tunneling.E0_prod.units)
        
    def test_repr(self):
        """
        Test that an Eckart object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        tunneling = None
        exec('tunneling = {0!r}'.format(self.tunneling))
        self.assertAlmostEqual(self.tunneling.frequency.value, tunneling.frequency.value, 2)
        self.assertEqual(self.tunneling.frequency.units, tunneling.frequency.units)
        self.assertAlmostEqual(self.tunneling.E0_reac.value, tunneling.E0_reac.value, 3)
        self.assertEqual(self.tunneling.E0_reac.units, tunneling.E0_reac.units)
        self.assertAlmostEqual(self.tunneling.E0_TS.value, tunneling.E0_TS.value, 3)
        self.assertEqual(self.tunneling.E0_TS.units, tunneling.E0_TS.units)
        self.assertAlmostEqual(self.tunneling.E0_prod.value, tunneling.E0_prod.value, 3)
        self.assertEqual(self.tunneling.E0_prod.units, tunneling.E0_prod.units)
