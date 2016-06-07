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
This script contains unit tests of the :mod:`rmgpy.thermo.nasa` module.
"""

import unittest
import numpy

from rmgpy.thermo.nasa import NASA, NASAPolynomial
import rmgpy.constants as constants

################################################################################

class TestNASA(unittest.TestCase):
    """
    Contains unit tests of the MultiNASA class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.coeffs_low = [4.03055,-0.00214171,4.90611e-05,-5.99027e-08,2.38945e-11,-11257.6,3.5613]
        self.coeffs_high = [-0.307954,0.0245269,-1.2413e-05,3.07724e-09,-3.01467e-13,-10693,22.628]
        self.Tmin = 300.
        self.Tmax = 3000.
        self.Tint = 650.73
        self.E0 = -782292. # J/mol.
        self.comment = "C2H6"
        self.nasa = NASA(
            polynomials = [
                NASAPolynomial(coeffs=self.coeffs_low, Tmin=(self.Tmin,"K"), Tmax=(self.Tint,"K")),
                NASAPolynomial(coeffs=self.coeffs_high, Tmin=(self.Tint,"K"), Tmax=(self.Tmax,"K")),
            ],
            Tmin = (self.Tmin,"K"),
            Tmax = (self.Tmax,"K"),
            E0 = (self.E0, "J/mol"),
            comment = self.comment,
        )
    
    def test_polyLow(self):
        """
        Test that the NASA low-temperature polynomial was properly set.
        """
        self.assertEqual(len(self.nasa.poly1.coeffs), len(self.coeffs_low))
        for coeff0, coeff in zip(self.nasa.poly1.coeffs, self.coeffs_low):
            self.assertAlmostEqual(coeff / coeff0, 1.0, 6)
        self.assertEqual(self.nasa.poly1.Tmin.value_si, self.Tmin)
        self.assertEqual(self.nasa.poly1.Tmax.value_si, self.Tint)
    
    def test_polyHigh(self):
        """
        Test that the NASA high-temperature polynomial was properly set.
        """
        self.assertEqual(len(self.nasa.poly2.coeffs), len(self.coeffs_high))
        for coeff0, coeff in zip(self.nasa.poly2.coeffs, self.coeffs_high):
            self.assertAlmostEqual(coeff / coeff0, 1.0, 6)
        self.assertEqual(self.nasa.poly2.Tmin.value_si, self.Tint)
        self.assertEqual(self.nasa.poly2.Tmax.value_si, self.Tmax)
    
    def test_Tmin(self):
        """
        Test that the NASA Tmin property was properly set.
        """
        self.assertAlmostEqual(self.nasa.Tmin.value_si / self.Tmin, 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.Tmin, self.Tmin))
    
    def test_Tmax(self):
        """
        Test that the NASA Tmax property was properly set.
        """
        self.assertAlmostEqual(self.nasa.Tmax.value_si / self.Tmax, 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.Tmax, self.Tmax))

    def test_E0(self):
        """
        Test that the NASA E0 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.E0.value_si / self.E0, 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.Tmax, self.Tmax))
        
    def test_Comment(self):
        """
        Test that the NASA comment property was properly set.
        """
        self.assertEqual(self.nasa.comment, self.comment)

    def test_isTemperatureValid(self):
        """
        Test the NASA.isTemperatureValid() method.
        """
        Tdata = [200,400,600,800,1000,1200,1400,1600,1800,2000]
        validdata = [False,True,True,True,True,True,True,True,True,True]
        for T, valid in zip(Tdata, validdata):
            valid0 = self.nasa.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
        
    def test_getHeatCapacity(self):
        """
        Test the NASA.getHeatCapacity() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        Cpexplist = numpy.array([7.80157, 10.5653, 12.8213, 14.5817, 15.9420, 16.9861, 17.78645, 18.4041, 18.8883]) * constants.R
        for T, Cpexp in zip(Tlist, Cpexplist):
            Cpact = self.nasa.getHeatCapacity(T)
            self.assertAlmostEqual(Cpexp / Cpact, 1.0, 4, '{0} != {1}'.format(Cpexp, Cpact))
        
    def test_getEnthalpy(self):
        """
        Test the NASA.getEnthalpy() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        Hexplist = numpy.array([-22.7613, -12.1027, -6.14236, -2.16615, 0.743456, 2.99256, 4.79397, 6.27334, 7.51156]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.nasa.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 3, '{0} != {1}'.format(Hexp, Hact))
        
    def test_getEntropy(self):
        """
        Test the NASA.getEntropy() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        Sexplist = numpy.array([29.6534, 33.3516, 36.7131, 39.7715, 42.5557, 45.0952, 47.4179, 49.5501, 51.5152]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.nasa.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getFreeEnergy(self):
        """
        Test the NASA.getFreeEnergy() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        for T in Tlist:
            Gexp = self.nasa.getEnthalpy(T) - T * self.nasa.getEntropy(T)
            Gact = self.nasa.getFreeEnergy(T)
            self.assertAlmostEqual(Gexp / Gact, 1.0, 4, '{0} != {1}'.format(Gexp, Gact))
    
    def test_pickle(self):
        """
        Test that a NASA object can be pickled and unpickled with no loss of
        information.
        """
        import cPickle
        nasa = cPickle.loads(cPickle.dumps(self.nasa))
        self.assertEqual(len(self.nasa.poly1.coeffs), len(nasa.poly1.coeffs))
        for coeff0, coeff in zip(self.nasa.poly1.coeffs, nasa.poly1.coeffs):
            self.assertAlmostEqual(coeff / coeff0, 1.0, 6)
        self.assertEqual(self.nasa.poly1.Tmin.value, nasa.poly1.Tmin.value)
        self.assertEqual(self.nasa.poly1.Tmin.units, nasa.poly1.Tmin.units)
        self.assertEqual(self.nasa.poly1.Tmax.value, nasa.poly1.Tmax.value)
        self.assertEqual(self.nasa.poly1.Tmax.units, nasa.poly1.Tmax.units)
        self.assertEqual(self.nasa.poly1.comment, nasa.poly1.comment)
        self.assertEqual(len(self.nasa.poly2.coeffs), len(nasa.poly2.coeffs))
        for coeff0, coeff in zip(self.nasa.poly2.coeffs, nasa.poly2.coeffs):
            self.assertAlmostEqual(coeff / coeff0, 1.0, 6)
        self.assertEqual(self.nasa.poly2.Tmin.value, nasa.poly2.Tmin.value)
        self.assertEqual(self.nasa.poly2.Tmin.units, nasa.poly2.Tmin.units)
        self.assertEqual(self.nasa.poly2.Tmax.value, nasa.poly2.Tmax.value)
        self.assertEqual(self.nasa.poly2.Tmax.units, nasa.poly2.Tmax.units)
        self.assertEqual(self.nasa.poly2.comment, nasa.poly2.comment)
        self.assertEqual(self.nasa.Tmin.value, nasa.Tmin.value)
        self.assertEqual(self.nasa.Tmin.units, nasa.Tmin.units)
        self.assertEqual(self.nasa.Tmax.value, nasa.Tmax.value)
        self.assertEqual(self.nasa.Tmax.units, nasa.Tmax.units)
        self.assertEqual(self.nasa.E0.value, nasa.E0.value)
        self.assertEqual(self.nasa.E0.units, nasa.E0.units)
        self.assertEqual(self.nasa.comment, nasa.comment)

    def test_repr(self):
        """
        Test that a NASA object can be reconstructed from its repr() output
        with no loss of information.
        """
        nasa = None
        exec('nasa = {0!r}'.format(self.nasa))
        self.assertEqual(len(self.nasa.poly1.coeffs), len(nasa.poly1.coeffs))
        for coeff0, coeff in zip(self.nasa.poly1.coeffs, nasa.poly1.coeffs):
            self.assertAlmostEqual(coeff / coeff0, 1.0, 6)
        self.assertEqual(self.nasa.poly1.Tmin.value, nasa.poly1.Tmin.value)
        self.assertEqual(self.nasa.poly1.Tmin.units, nasa.poly1.Tmin.units)
        self.assertEqual(self.nasa.poly1.Tmax.value, nasa.poly1.Tmax.value)
        self.assertEqual(self.nasa.poly1.Tmax.units, nasa.poly1.Tmax.units)
        self.assertEqual(self.nasa.poly1.comment, nasa.poly1.comment)
        self.assertEqual(len(self.nasa.poly2.coeffs), len(nasa.poly2.coeffs))
        for coeff0, coeff in zip(self.nasa.poly2.coeffs, nasa.poly2.coeffs):
            self.assertAlmostEqual(coeff / coeff0, 1.0, 6)
        self.assertEqual(self.nasa.poly2.Tmin.value, nasa.poly2.Tmin.value)
        self.assertEqual(self.nasa.poly2.Tmin.units, nasa.poly2.Tmin.units)
        self.assertEqual(self.nasa.poly2.Tmax.value, nasa.poly2.Tmax.value)
        self.assertEqual(self.nasa.poly2.Tmax.units, nasa.poly2.Tmax.units)
        self.assertEqual(self.nasa.poly2.comment, nasa.poly2.comment)
        self.assertEqual(self.nasa.Tmin.value, nasa.Tmin.value)
        self.assertEqual(self.nasa.Tmin.units, nasa.Tmin.units)
        self.assertEqual(self.nasa.Tmax.value, nasa.Tmax.value)
        self.assertEqual(self.nasa.Tmax.units, nasa.Tmax.units)
        self.assertEqual(self.nasa.E0.value, nasa.E0.value)
        self.assertEqual(self.nasa.E0.units, nasa.E0.units)
        self.assertEqual(self.nasa.comment, nasa.comment)
        
    def test_toCantera(self):
        """
        Test that conversion to a Cantera NasaPoly2 object works
        """
        nasapoly2 = self.nasa.toCantera()
        # NasaPoly2 units use J/kmol rather than J/mol
        self.assertAlmostEqual(self.nasa.getEnthalpy(900), nasapoly2.h(900)/1000, 1)
        self.assertAlmostEqual(self.nasa.getEntropy(700), nasapoly2.s(700)/1000, 1)
