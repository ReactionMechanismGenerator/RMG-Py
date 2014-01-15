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
This script contains unit tests of the :mod:`rmgpy.kinetics.arrhenius` module.
"""

import unittest
import numpy

from rmgpy.kinetics.kineticsdata import KineticsData, PDepKineticsData
import rmgpy.constants as constants

################################################################################

class TestKineticsData(unittest.TestCase):
    """
    Contains unit tests of the :class:`KineticsData` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Tdata = numpy.array([300,400,500,600,700,800,900,1000,1500,2000], numpy.float64)
        self.kdata = numpy.array([4.73e-19,3.93e-17,6.51e-16,4.60e-15,2.03e-14,6.28e-14,1.58e-13,3.31e-13,3.72e-12,1.49e-11], numpy.float64)
        self.Tmin = 300.
        self.Tmax = 3000.
        self.comment = 'H + CH4 <=> H2 + CH3 (RPMD)'
        self.kinetics = KineticsData(
            Tdata = (self.Tdata,"K"),
            kdata = (self.kdata,"cm^3/(molecule*s)"),
            Tmin = (self.Tmin,"K"),
            Tmax = (self.Tmax,"K"),
            comment = self.comment,
        )
    
    def test_Tdata(self):
        """
        Test that the KineticsData Tdata property was properly set.
        """
        self.assertEqual(self.kinetics.Tdata.value_si.shape, self.Tdata.shape)
        for T, T0 in zip(self.kinetics.Tdata.value_si, self.Tdata):
            self.assertAlmostEqual(T, T0, 4)
        
    def test_kdata(self):
        """
        Test that the KineticsData kdata property was properly set.
        """
        self.assertEqual(self.kinetics.kdata.value_si.shape, self.kdata.shape)
        for k, k0 in zip(self.kinetics.kdata.value_si, self.kdata):
            k0 *= constants.Na * 1e-6
            self.assertAlmostEqual(k, k0, delta=1e-6*k0)
        
    def test_Tmin(self):
        """
        Test that the KineticsData Tmin property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Tmin.value_si, self.Tmin, 6)
        
    def test_Tmax(self):
        """
        Test that the KineticsData Tmax property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Tmax.value_si, self.Tmax, 6)
        
    def test_comment(self):
        """
        Test that the KineticsData comment property was properly set.
        """
        self.assertEqual(self.kinetics.comment, self.comment)
        
    def test_isTemperatureValid(self):
        """
        Test the KineticsData.isTemperatureValid() method.
        """
        Tdata = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        validdata = numpy.array([False,True,True,True,True,True,True,True,True,True], numpy.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
                
    def test_getRateCoefficient(self):
        """
        Test the KineticsData.getRateCoefficient() method.
        """
        Tlist = numpy.array([300,400,600,800,1000,1200,1400,1600,1800,2000])
        kexplist = numpy.array([2.84847e-01, 2.36670e+01, 2.77019e+03, 3.78191e+04, 1.99333e+05, 5.24644e+05, 1.38086e+06, 2.95680e+06, 5.15086e+06, 8.97299e+06])
        for T, kexp in zip(Tlist, kexplist):
            kact = self.kinetics.getRateCoefficient(T)
            self.assertAlmostEqual(kexp, kact, delta=1e-4*kexp)

    def test_pickle(self):
        """
        Test that a KineticsData object can be pickled and unpickled with no
        loss of information.
        """
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics,-1))
        self.assertEqual(self.kinetics.Tdata.value.shape, kinetics.Tdata.value.shape)
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            self.assertAlmostEqual(T, T0, 4)
        self.assertEqual(self.kinetics.Tdata.units, kinetics.Tdata.units)
        self.assertEqual(self.kinetics.kdata.value.shape, kinetics.kdata.value.shape)
        for k, k0 in zip(self.kinetics.kdata.value, kinetics.kdata.value):
            self.assertAlmostEqual(k, k0, 4)
        self.assertEqual(self.kinetics.kdata.units, kinetics.kdata.units)
        self.assertAlmostEqual(self.kinetics.Tmin.value, kinetics.Tmin.value, 4)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertAlmostEqual(self.kinetics.Tmax.value, kinetics.Tmax.value, 4)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
    
    def test_repr(self):
        """
        Test that a KineticsData object can be reconstructed from its repr()
        output with no loss of information.
        """
        kinetics = None
        exec('kinetics = {0!r}'.format(self.kinetics))
        self.assertEqual(self.kinetics.Tdata.value.shape, kinetics.Tdata.value.shape)
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            self.assertAlmostEqual(T, T0, 4)
        self.assertEqual(self.kinetics.Tdata.units, kinetics.Tdata.units)
        self.assertEqual(self.kinetics.kdata.value.shape, kinetics.kdata.value.shape)
        for k, k0 in zip(self.kinetics.kdata.value, kinetics.kdata.value):
            self.assertAlmostEqual(k, k0, 4)
        self.assertEqual(self.kinetics.kdata.units, kinetics.kdata.units)
        self.assertAlmostEqual(self.kinetics.Tmin.value, kinetics.Tmin.value, 4)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertAlmostEqual(self.kinetics.Tmax.value, kinetics.Tmax.value, 4)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

################################################################################

class TestPDepKineticsData(unittest.TestCase):
    """
    Contains unit tests of the :class:`PDepKineticsData` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Tdata = numpy.array([300,400,500,600,700,800,900,1000,1500,2000], numpy.float64)
        self.Pdata = numpy.array([1e-1,1e0,1e1], numpy.float64)
        self.kdata = numpy.array([
            [4.73e-21,3.93e-19,6.51e-18,4.60e-17,2.03e-16,6.28e-16,1.58e-15,3.31e-15,3.72e-14,1.49e-13],
            [4.73e-20,3.93e-18,6.51e-17,4.60e-16,2.03e-15,6.28e-15,1.58e-14,3.31e-14,3.72e-13,1.49e-12],
            [4.73e-19,3.93e-17,6.51e-16,4.60e-15,2.03e-14,6.28e-14,1.58e-13,3.31e-13,3.72e-12,1.49e-11],
            ], numpy.float64).T
        self.Tmin = 300.
        self.Tmax = 3000.
        self.Pmin = 1e-1
        self.Pmax = 1e1
        self.comment = ''
        self.kinetics = PDepKineticsData(
            Tdata = (self.Tdata,"K"),
            Pdata = (self.Pdata,"bar"),
            kdata = (self.kdata,"cm^3/(molecule*s)"),
            Tmin = (self.Tmin,"K"),
            Tmax = (self.Tmax,"K"),
            Pmin = (self.Pmin,"bar"),
            Pmax = (self.Pmax,"bar"),
            comment = self.comment,
        )
    
    def test_Tdata(self):
        """
        Test that the PDepKineticsData Tdata property was properly set.
        """
        self.assertEqual(self.kinetics.Tdata.value_si.shape, self.Tdata.shape)
        for T, T0 in zip(self.kinetics.Tdata.value_si, self.Tdata):
            self.assertAlmostEqual(T, T0, 4)
        
    def test_Pdata(self):
        """
        Test that the PDepKineticsData Pdata property was properly set.
        """
        self.assertEqual(self.kinetics.Pdata.value_si.shape, self.Pdata.shape)
        for P, P0 in zip(self.kinetics.Pdata.value_si, self.Pdata):
            self.assertAlmostEqual(P*1e-5, P0, 4)
        
    def test_kdata(self):
        """
        Test that the PDepKineticsData kdata property was properly set.
        """
        self.assertEqual(self.kinetics.kdata.value_si.shape, self.kdata.shape)
        for i in range(self.kdata.shape[0]):
            for j in range(self.kdata.shape[1]):
                k0 = self.kdata[i,j] * constants.Na * 1e-6
                k = self.kinetics.kdata.value_si[i,j]
                self.assertAlmostEqual(k, k0, delta=1e-6*k0)
        
    def test_Tmin(self):
        """
        Test that the PDepKineticsData Tmin property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Tmin.value_si, self.Tmin, 6)
        
    def test_Tmax(self):
        """
        Test that the PDepKineticsData Tmax property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Tmax.value_si, self.Tmax, 6)
        
    def test_Pmin(self):
        """
        Test that the PDepKineticsData Pmin property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Pmin.value_si*1e-5, self.Pmin, 6)
        
    def test_Pmax(self):
        """
        Test that the PDepKineticsData Pmax property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Pmax.value_si*1e-5, self.Pmax, 6)
        
    def test_comment(self):
        """
        Test that the PDepKineticsData comment property was properly set.
        """
        self.assertEqual(self.kinetics.comment, self.comment)
        
    def test_isTemperatureValid(self):
        """
        Test the PDepKineticsData.isTemperatureValid() method.
        """
        Tdata = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        validdata = numpy.array([False,True,True,True,True,True,True,True,True,True], numpy.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
                
    def test_isPressureValid(self):
        """
        Test the PDepKineticsData.isPressureValid() method.
        """
        Pdata = numpy.array([1e3,1e4,1e5,1e6,1e7])
        validdata = numpy.array([False,True,True,True,False], numpy.bool)
        for P, valid in zip(Pdata, validdata):
            valid0 = self.kinetics.isPressureValid(P)
            self.assertEqual(valid0, valid)
                
    def test_getRateCoefficient(self):
        """
        Test the PDepKineticsData.getRateCoefficient() method.
        """
        Tlist = numpy.array([300,400,600,800,1000,1200,1400,1600,1800,2000])
        Plist = numpy.array([1e4,1e5,1e6])
        kexplist = numpy.array([
            [2.84847e-03, 2.36670e-01, 2.77019e+01, 3.78191e+02, 1.99333e+03, 5.24644e+03, 1.38086e+04, 2.95680e+04, 5.15086e+04, 8.97299e+04],
            [2.84847e-02, 2.36670e+00, 2.77019e+02, 3.78191e+03, 1.99333e+04, 5.24644e+04, 1.38086e+05, 2.95680e+05, 5.15086e+05, 8.97299e+05],
            [2.84847e-01, 2.36670e+01, 2.77019e+03, 3.78191e+04, 1.99333e+05, 5.24644e+05, 1.38086e+06, 2.95680e+06, 5.15086e+06, 8.97299e+06],
        ]).T
        for i in range(Tlist.shape[0]):
            for j in range(Plist.shape[0]):
                kexp = kexplist[i,j]
                kact = self.kinetics.getRateCoefficient(Tlist[i], Plist[j])
                self.assertAlmostEqual(kexp, kact, delta=1e-4*kexp)

    def test_pickle(self):
        """
        Test that a PDepKineticsData object can be pickled and unpickled with no
        loss of information.
        """
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics,-1))
        self.assertEqual(self.kinetics.Tdata.value.shape, kinetics.Tdata.value.shape)
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            self.assertAlmostEqual(T, T0, 4)
        self.assertEqual(self.kinetics.Tdata.units, kinetics.Tdata.units)
        self.assertEqual(self.kinetics.Pdata.value.shape, kinetics.Pdata.value.shape)
        for P, P0 in zip(self.kinetics.Pdata.value, kinetics.Pdata.value):
            self.assertAlmostEqual(P, P0, 4)
        self.assertEqual(self.kinetics.Pdata.units, kinetics.Pdata.units)
        self.assertEqual(self.kinetics.kdata.value.shape, kinetics.kdata.value.shape)
        for i in range(self.kinetics.kdata.value.shape[0]):
            for j in range(self.kinetics.kdata.value.shape[1]):
                k0 = self.kinetics.kdata.value[i,j]
                k = kinetics.kdata.value[i,j]
                self.assertAlmostEqual(k, k0, delta=1e-6*k0)
        self.assertEqual(self.kinetics.kdata.units, kinetics.kdata.units)
        self.assertAlmostEqual(self.kinetics.Tmin.value, kinetics.Tmin.value, 4)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertAlmostEqual(self.kinetics.Tmax.value, kinetics.Tmax.value, 4)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
    
    def test_repr(self):
        """
        Test that a PDepKineticsData object can be reconstructed from its repr()
        output with no loss of information.
        """
        kinetics = None
        exec('kinetics = {0!r}'.format(self.kinetics))
        self.assertEqual(self.kinetics.Tdata.value.shape, kinetics.Tdata.value.shape)
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            self.assertAlmostEqual(T, T0, 4)
        self.assertEqual(self.kinetics.Tdata.units, kinetics.Tdata.units)
        self.assertEqual(self.kinetics.Pdata.value.shape, kinetics.Pdata.value.shape)
        for P, P0 in zip(self.kinetics.Pdata.value, kinetics.Pdata.value):
            self.assertAlmostEqual(P, P0, 4)
        self.assertEqual(self.kinetics.Pdata.units, kinetics.Pdata.units)
        self.assertEqual(self.kinetics.kdata.value.shape, kinetics.kdata.value.shape)
        for i in range(self.kinetics.kdata.value.shape[0]):
            for j in range(self.kinetics.kdata.value.shape[1]):
                k0 = self.kinetics.kdata.value[i,j]
                k = kinetics.kdata.value[i,j]
                self.assertAlmostEqual(k, k0, delta=1e-6*k0)
        self.assertEqual(self.kinetics.kdata.units, kinetics.kdata.units)
        self.assertAlmostEqual(self.kinetics.Tmin.value, kinetics.Tmin.value, 4)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertAlmostEqual(self.kinetics.Tmax.value, kinetics.Tmax.value, 4)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
