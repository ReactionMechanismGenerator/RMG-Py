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
This script contains unit tests of the :mod:`rmgpy.thermo.thermodata` module.
"""

import unittest
import numpy

from rmgpy.thermo.thermodata import ThermoData
import rmgpy.constants as constants

################################################################################

class TestThermoData(unittest.TestCase):
    """
    Contains unit tests of the :class:`ThermoData` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.H298 = -32.9725
        self.S298 = 27.5727
        self.Tdata = numpy.array([300,400,500,600,800,1000,1500])
        self.Cpdata = numpy.array([6.3827,7.80327,9.22175,10.5528,12.8323,14.6013,17.4089])
        self.Cp0 = 4.0
        self.CpInf = 21.5 
        self.Tmin = 100.
        self.Tmax = 3000.
        self.E0 = -782292.
        self.comment = 'C2H6'
        self.thermodata = ThermoData(
            Tdata = (self.Tdata,"K"),
            Cpdata = (self.Cpdata*constants.R,"J/(mol*K)"),
            H298 = (self.H298*0.001*constants.R*298.,"kJ/mol"),
            S298 = (self.S298*constants.R,"J/(mol*K)"),
            Cp0 = (self.Cp0*constants.R,"J/(mol*K)"),
            CpInf = (self.CpInf*constants.R,"J/(mol*K)"),
            Tmin = (self.Tmin,"K"),
            Tmax = (self.Tmax,"K"),
            E0 = (self.E0,'J/mol'),
            comment = self.comment,
        )
    
    def test_Tdata(self):
        """
        Test that the ThermoData Tdata property was properly set.
        """
        self.assertEqual(self.thermodata.Tdata.value_si.shape, self.Tdata.shape)
        for T, T0 in zip(self.thermodata.Tdata.value_si, self.Tdata):
            self.assertAlmostEqual(T, T0, 4)
    
    def test_Cpdata(self):
        """
        Test that the ThermoData Cpdata property was properly set.
        """
        self.assertEqual(self.thermodata.Cpdata.value_si.shape, self.Cpdata.shape)
        for Cp, Cp0 in zip(self.thermodata.Cpdata.value_si / constants.R, self.Cpdata):
            self.assertAlmostEqual(Cp, Cp0, 4)
    
    def test_H298(self):
        """
        Test that the ThermoData H298 property was properly set.
        """
        self.assertAlmostEqual(self.thermodata.H298.value_si / constants.R / 298., self.H298, 4)
    
    def test_S298(self):
        """
        Test that the ThermoData S298 property was properly set.
        """
        self.assertAlmostEqual(self.thermodata.S298.value_si / constants.R, self.S298, 4)
    
    def test_Cp0(self):
        """
        Test that the ThermoData Cp0 property was properly set.
        """
        self.assertAlmostEqual(self.thermodata.Cp0.value_si / constants.R, self.Cp0, 4)
    
    def test_CpInf(self):
        """
        Test that the ThermoData CpInf property was properly set.
        """
        self.assertAlmostEqual(self.thermodata.CpInf.value_si / constants.R, self.CpInf, 4)
        
    def test_Tmin(self):
        """
        Test that the ThermoData Tmin property was properly set.
        """
        self.assertAlmostEqual(self.thermodata.Tmin.value_si, self.Tmin, 6)
    
    def test_Tmax(self):
        """
        Test that the ThermoData Tmax property was properly set.
        """
        self.assertAlmostEqual(self.thermodata.Tmax.value_si, self.Tmax, 6)

    def test_E0(self):
        """
        Test that the ThermoData E0 property was properly set.
        """
        self.assertAlmostEqual(self.thermodata.E0.value_si, self.E0, 6)
    
    def test_Comment(self):
        """
        Test that the ThermoData comment property was properly set.
        """
        self.assertEqual(self.thermodata.comment, self.comment)
    
    def test_isTemperatureValid(self):
        """
        Test the ThermoData.isTemperatureValid() method.
        """
        Tdata = [200,400,600,800,1000,1200,1400,1600,1800,2000]
        validdata = [True,True,True,True,True,True,True,True,True,True]
        for T, valid in zip(Tdata, validdata):
            valid0 = self.thermodata.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
        
    def test_getHeatCapacity(self):
        """
        Test the ThermoData.getHeatCapacity() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        Cpexplist = numpy.array([4.96208, 7.80327, 10.5528, 12.8323, 14.6013, 15.7243, 16.8473, 17.9704, 19.0934, 20.2165]) * constants.R
        for T, Cpexp in zip(Tlist, Cpexplist):
            Cpact = self.thermodata.getHeatCapacity(T)
            self.assertAlmostEqual(Cpexp, Cpact, 2)

    def test_getEnthalpy(self):
        """
        Test the ThermoData.getEnthalpy() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        Hexplist = numpy.array([-51.9015, -22.7594, -12.1063, -6.15660, -2.18192, 0.708869, 2.93415, 4.74350, 6.27555, 7.61349]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.thermodata.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e0)
        
    def test_getEntropy(self):
        """
        Test the ThermoData.getEntropy() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        Sexplist = numpy.array([25.3347, 29.6460, 33.3386, 36.6867, 39.7402, 42.5016, 45.0098, 47.3328, 49.5142, 51.5841]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.thermodata.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, 3)
            
    def test_getFreeEnergy(self):
        """
        Test the ThermoData.getFreeEnergy() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        for T in Tlist:
            Gexp = self.thermodata.getEnthalpy(T) - T * self.thermodata.getEntropy(T)
            Gact = self.thermodata.getFreeEnergy(T)
            self.assertAlmostEqual(Gexp, Gact, 3)
    
    def test_pickle(self):
        """
        Test that a ThermoData object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thermodata = cPickle.loads(cPickle.dumps(self.thermodata))
        self.assertEqual(self.thermodata.Tdata.value.shape, thermodata.Tdata.value.shape)
        for T, T0 in zip(self.thermodata.Tdata.value, thermodata.Tdata.value):
            self.assertAlmostEqual(T, T0, 4)
        self.assertEqual(self.thermodata.Tdata.units, thermodata.Tdata.units)
        self.assertEqual(self.thermodata.Cpdata.value.shape, thermodata.Cpdata.value.shape)
        for Cp, Cp0 in zip(self.thermodata.Cpdata.value, thermodata.Cpdata.value):
            self.assertAlmostEqual(Cp, Cp0, 3)
        self.assertEqual(self.thermodata.Cpdata.units, thermodata.Cpdata.units)
        self.assertAlmostEqual(self.thermodata.H298.value, thermodata.H298.value, 4)
        self.assertEqual(self.thermodata.H298.units, thermodata.H298.units)
        self.assertAlmostEqual(self.thermodata.S298.value, thermodata.S298.value, 2)
        self.assertEqual(self.thermodata.S298.units, thermodata.S298.units)
        self.assertAlmostEqual(self.thermodata.Cp0.value, thermodata.Cp0.value, 4)
        self.assertEqual(self.thermodata.Cp0.units, thermodata.Cp0.units)
        self.assertAlmostEqual(self.thermodata.CpInf.value, thermodata.CpInf.value, 3)
        self.assertEqual(self.thermodata.CpInf.units, thermodata.CpInf.units)
        self.assertAlmostEqual(self.thermodata.Tmin.value, thermodata.Tmin.value, 4)
        self.assertEqual(self.thermodata.Tmin.units, thermodata.Tmin.units)
        self.assertAlmostEqual(self.thermodata.Tmax.value, thermodata.Tmax.value, 4)
        self.assertEqual(self.thermodata.Tmax.units, thermodata.Tmax.units)
        self.assertAlmostEqual(self.thermodata.E0.value, thermodata.E0.value, 4)
        self.assertEqual(self.thermodata.E0.units, thermodata.E0.units)
        self.assertEqual(self.thermodata.comment, thermodata.comment)
    
    def test_repr(self):
        """
        Test that a ThermoData object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        thermodata = None
        exec('thermodata = {0!r}'.format(self.thermodata))
        self.assertEqual(self.thermodata.Tdata.value.shape, thermodata.Tdata.value.shape)
        for T, T0 in zip(self.thermodata.Tdata.value, thermodata.Tdata.value):
            self.assertAlmostEqual(T, T0, 4)
        self.assertEqual(self.thermodata.Tdata.units, thermodata.Tdata.units)
        self.assertEqual(self.thermodata.Cpdata.value.shape, thermodata.Cpdata.value.shape)
        for Cp, Cp0 in zip(self.thermodata.Cpdata.value, thermodata.Cpdata.value):
            self.assertAlmostEqual(Cp, Cp0, 3)
        self.assertEqual(self.thermodata.Cpdata.units, thermodata.Cpdata.units)
        self.assertAlmostEqual(self.thermodata.H298.value, thermodata.H298.value, 4)
        self.assertEqual(self.thermodata.H298.units, thermodata.H298.units)
        self.assertAlmostEqual(self.thermodata.S298.value, thermodata.S298.value, 2)
        self.assertEqual(self.thermodata.S298.units, thermodata.S298.units)
        self.assertAlmostEqual(self.thermodata.Cp0.value, thermodata.Cp0.value, 4)
        self.assertEqual(self.thermodata.Cp0.units, thermodata.Cp0.units)
        self.assertAlmostEqual(self.thermodata.CpInf.value, thermodata.CpInf.value, 3)
        self.assertEqual(self.thermodata.CpInf.units, thermodata.CpInf.units)
        self.assertAlmostEqual(self.thermodata.Tmin.value, thermodata.Tmin.value, 4)
        self.assertEqual(self.thermodata.Tmin.units, thermodata.Tmin.units)
        self.assertAlmostEqual(self.thermodata.Tmax.value, thermodata.Tmax.value, 4)
        self.assertEqual(self.thermodata.Tmax.units, thermodata.Tmax.units)
        self.assertAlmostEqual(self.thermodata.E0.value, thermodata.E0.value, 4)
        self.assertEqual(self.thermodata.E0.units, thermodata.E0.units)
        self.assertEqual(self.thermodata.comment, thermodata.comment)
