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
This script contains unit tests of the :mod:`rmgpy.thermo` object conversion
methods.
"""

import unittest
import math
import numpy

from rmgpy.thermo import Wilhoit, NASA, NASAPolynomial, ThermoData
import rmgpy.constants as constants

################################################################################

class TestConverter(unittest.TestCase):
    """
    Contains unit tests of the thermodynamics model conversion functions.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """        
        self.wilhoit = Wilhoit(
            Cp0 = (4.0*constants.R,"J/(mol*K)"),
            CpInf = (21.5*constants.R,"J/(mol*K)"),
            a0 = 0.0977518,
            a1 = -16.3067,
            a2 = 26.2524,
            a3 = -12.6785,
            B = (1068.68,"K"),
            H0 = (-94.088*constants.R,"kJ/mol"),
            S0 = (-118.46*constants.R,"J/(mol*K)"),
            Tmin = (10,"K"),
            Tmax = (3000,"K"),
            comment = 'C2H6',
        )
        self.nasa = NASA(
            polynomials = [
                NASAPolynomial(coeffs=[4.03055,-0.00214171,4.90611e-05,-5.99027e-08,2.38945e-11,-11257.6,3.5613], Tmin=(10,"K"), Tmax=(650.73,"K")),
                NASAPolynomial(coeffs=[-0.307954,0.0245269,-1.2413e-05,3.07724e-09,-3.01467e-13,-10693,22.628], Tmin=(650.73,"K"), Tmax=(3000,"K")),
            ],
            Tmin = (10,"K"),
            Tmax = (3000,"K"),
            E0 = (-93.6077,'kJ/mol'),
            comment = 'C2H6',
        )
        self.thermodata = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],"K"),
            Cpdata = (numpy.array([6.38268,7.80327,9.22175,10.5528,12.8323,14.6013,17.40890])*constants.R,"J/(mol*K)"),
            H298 = (-81.7,"kJ/mol"),
            S298 = (27.5727*constants.R,"J/(mol*K)"),
            Cp0 = (4.0*constants.R,"J/(mol*K)"),
            CpInf = (21.5*constants.R,"J/(mol*K)"),
            Tmin = (10,"K"),
            Tmax = (3000,"K"),
            E0 = (-93.6077,'kJ/mol'),
            comment = 'C2H6',
        )
    
    def test_convert_Wilhoit_to_NASA(self):
        """
        Test the conversion of a Wilhoit model to a NASA model.
        """
        wilhoit = self.wilhoit
        nasa = wilhoit.toNASA(Tmin=10, Tmax=3000, Tint=1000)
        Tlist = numpy.arange(10, 3000, 10)
        for T in Tlist:
            Cp_wilhoit = wilhoit.getHeatCapacity(T)
            Cp_nasa = nasa.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_nasa, Cp_wilhoit, delta=1e0)
            H_wilhoit = wilhoit.getEnthalpy(T)
            H_nasa = nasa.getEnthalpy(T)
            self.assertAlmostEqual(H_nasa, H_wilhoit, delta=2e1)
            S_wilhoit = wilhoit.getEntropy(T)
            S_nasa = nasa.getEntropy(T)
            self.assertAlmostEqual(S_nasa, S_wilhoit, delta=1e0)
        self.assertAlmostEqual(wilhoit.E0.value_si, nasa.E0.value_si, delta=1e1)

    def test_convert_Wilhoit_to_ThermoData(self):
        """
        Test the conversion of a Wilhoit model to a ThermoData model.
        """
        wilhoit = self.wilhoit
        thermodata = wilhoit.toThermoData()
        Tlist = numpy.array([300,400,500,600,800,1000,1500])
        for T in Tlist:
            Cp_wilhoit = wilhoit.getHeatCapacity(T)
            Cp_thermodata = thermodata.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_thermodata, Cp_wilhoit, 4)
        T = 298
        H_wilhoit = wilhoit.getEnthalpy(T)
        H_thermodata = thermodata.getEnthalpy(T)
        self.assertAlmostEqual(H_thermodata, H_wilhoit, 4)
        S_wilhoit = wilhoit.getEntropy(T)
        S_thermodata = thermodata.getEntropy(T)
        self.assertAlmostEqual(S_thermodata, S_wilhoit, 4)
        self.assertAlmostEqual(wilhoit.E0.value_si, thermodata.E0.value_si,  delta=1e1)

    def test_convert_NASA_to_Wilhoit(self):
        """
        Test the conversion of a NASA model to a Wilhoit model.
        """
        nasa = self.nasa
        wilhoit = nasa.toWilhoit(Cp0=self.wilhoit.Cp0.value_si, CpInf=self.wilhoit.CpInf.value_si)
        Tlist = numpy.arange(10, 3000, 10)
        for T in Tlist:
            Cp_wilhoit = wilhoit.getHeatCapacity(T)
            Cp_nasa = nasa.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_nasa, Cp_wilhoit, delta=1e0)
            H_wilhoit = wilhoit.getEnthalpy(T)
            H_nasa = nasa.getEnthalpy(T)
            self.assertAlmostEqual(H_nasa, H_wilhoit, delta=2e1)
            S_wilhoit = wilhoit.getEntropy(T)
            S_nasa = nasa.getEntropy(T)
            self.assertAlmostEqual(S_nasa, S_wilhoit, delta=1e0)
        self.assertAlmostEqual(nasa.E0.value_si, wilhoit.E0.value_si, delta=2e1)
        
    def test_convert_NASA_to_ThermoData(self):
        """
        Test the conversion of a NASA model to a ThermoData model.
        """
        nasa = self.nasa
        thermodata = nasa.toThermoData(Cp0=self.thermodata.Cp0.value_si, CpInf=self.thermodata.CpInf.value_si)
        Tlist = numpy.array([300,400,500,600,800,1000,1500])
        for T in Tlist:
            Cp_thermodata = thermodata.getHeatCapacity(T)
            Cp_nasa = nasa.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_nasa, Cp_thermodata, 4)
        T = 298
        H_thermodata = thermodata.getEnthalpy(T)
        H_nasa = nasa.getEnthalpy(T)
        self.assertAlmostEqual(H_nasa, H_thermodata, 4)
        S_thermodata = thermodata.getEntropy(T)
        S_nasa = nasa.getEntropy(T)
        self.assertAlmostEqual(S_nasa, S_thermodata, 4)
        self.assertAlmostEqual(nasa.E0.value_si, thermodata.E0.value_si,  delta=1e1)
        
    def test_convert_ThermoData_to_Wilhoit(self):
        """
        Test the conversion of a ThermoData model to a Wilhoit model.
        """
        thermodata = self.thermodata
        wilhoit = thermodata.toWilhoit()
        Tlist = numpy.array([300,400,500,600,800,1000,1500])
        for T in Tlist:
            Cp_wilhoit = wilhoit.getHeatCapacity(T)
            Cp_thermodata = thermodata.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_thermodata, Cp_wilhoit, 3)
        T = 298
        H_wilhoit = wilhoit.getEnthalpy(T)
        H_thermodata = thermodata.getEnthalpy(T)
        self.assertAlmostEqual(H_thermodata, H_wilhoit, 3)
        S_wilhoit = wilhoit.getEntropy(T)
        S_thermodata = thermodata.getEntropy(T)
        self.assertAlmostEqual(S_thermodata, S_wilhoit, 3)
        self.assertAlmostEqual(thermodata.E0.value_si, wilhoit.E0.value_si, delta=1e1)

    def test_convert_ThermoData_to_NASA(self):
        """
        Test the conversion of a ThermoData model to a NASA model.
        """
        thermodata = self.thermodata
        nasa = thermodata.toNASA(Tmin=10, Tmax=3000, Tint=1000)
        Tlist = numpy.array([300,400,500,600,800,1000,1500])
        for T in Tlist:
            Cp_thermodata = thermodata.getHeatCapacity(T)
            Cp_nasa = nasa.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_nasa, Cp_thermodata, delta=1e0)
        T = 298
        H_thermodata = thermodata.getEnthalpy(T)
        H_nasa = nasa.getEnthalpy(T)
        self.assertAlmostEqual(H_nasa, H_thermodata, delta=1e0)
        S_thermodata = thermodata.getEntropy(T)
        S_nasa = nasa.getEntropy(T)
        self.assertAlmostEqual(S_nasa, S_thermodata, delta=1e0)
        self.assertAlmostEqual(thermodata.E0.value_si, nasa.E0.value_si, delta=1e1)
        
    def test_Wilhoit_NASA_Wilhoit(self):
        """
        Test round-trip conversion from Wilhoit to NASA and back
        """
        wilhoit1 = self.wilhoit
        nasa = wilhoit1.toNASA(Tmin=10, Tmax=3000, Tint=1000)
        wilhoit2 = nasa.toWilhoit(Cp0=self.wilhoit.Cp0.value_si, CpInf=self.wilhoit.CpInf.value_si)
        Tlist = numpy.arange(10, 3000, 10)
        for T in Tlist:
            Cp_1 = wilhoit1.getHeatCapacity(T)
            Cp_2 = wilhoit2.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_1, Cp_2, delta=1e0)
            H_1 = wilhoit1.getEnthalpy(T)
            H_2 = wilhoit2.getEnthalpy(T)
            self.assertAlmostEqual(H_1, H_2, delta=2e1)
            S_1 = wilhoit1.getEntropy(T)
            S_2 = wilhoit2.getEntropy(T)
            self.assertAlmostEqual(S_1, S_2, delta=1e0)
        self.assertAlmostEqual(wilhoit1.E0.value_si, wilhoit2.E0.value_si, delta=1e1)

    def test_Wilhoit_ThermoData_Wilhoit(self):
        """
        Test round-trip conversion from Wilhoit to ThermoData and back
        """
        wilhoit1 = self.wilhoit
        thermodata = wilhoit1.toThermoData()
        wilhoit2 = thermodata.toWilhoit()
        Tlist = numpy.arange(10, 3000, 10)
        for T in Tlist:
            Cp_1 = wilhoit1.getHeatCapacity(T)
            Cp_2 = wilhoit2.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_1, Cp_2, delta=1e0)
            H_1 = wilhoit1.getEnthalpy(T)
            H_2 = wilhoit2.getEnthalpy(T)
            self.assertAlmostEqual(H_1, H_2, delta=2e1)
            S_1 = wilhoit1.getEntropy(T)
            S_2 = wilhoit2.getEntropy(T)
            self.assertAlmostEqual(S_1, S_2, delta=1e0)
        self.assertAlmostEqual(wilhoit1.E0.value_si, wilhoit2.E0.value_si, delta=1e1)
            
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
