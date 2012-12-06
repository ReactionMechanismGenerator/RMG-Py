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
This script contains unit tests of the :mod:`rmgpy.statmech.schrodinger` 
module.
"""

import unittest

import numpy
from rmgpy.statmech.schrodinger import getPartitionFunction, getHeatCapacity, getEnthalpy, getEntropy, getDensityOfStates
import rmgpy.constants as constants

################################################################################

class TestSchrodinger(unittest.TestCase):
    """
    Contains unit tests of the various methods of the :mod:`schrodinger`
    module. The solution to the Schrodinger equation used for these tests is
    that of a linear rigid rotor with a rotational constant of 1 cm^-1.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.B = 1.0 * 11.96
        self.energy = lambda J: self.B * J * (J + 1)
        self.degeneracy = lambda J: 2 * J + 1
        self.n0 = 0
        
    def test_getPartitionFunction(self):
        """
        Test the getPartitionFunction() method.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([208.8907, 347.9285, 695.5234, 1043.118, 1390.713])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = getPartitionFunction(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1} within 4 figures'.format(Qexp, Qact))
            
    def test_getHeatCapacity(self):
        """
        Test the getHeatCapacity() method.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1, 1, 1, 1, 1])
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = getHeatCapacity(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1} within 4 figures'.format(Cvexp, Cvact))
       
    def test_getEnthalpy(self):
        """
        Test the getEnthalpy() method.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.9984012, 0.9990409, 0.9995205, 0.9996803, 0.9997603])
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = getEnthalpy(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1} within 4 figures'.format(Hexp, Hact))

    def test_getEntropy(self):
        """
        Test the getEntropy() method.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([6.340212, 6.851038, 7.544185, 7.949650, 8.237332])
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = getEntropy(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1} within 4 figures'.format(Sexp, Sact))

#    def test_getSumOfStates(self):
#        """
#        Test the getSumOfStates() method.
#        """
#        Elist = numpy.arange(0, 10., 0.01)
#        densStates = getDensityOfStates(Elist, self.energy, self.degeneracy, self.n0)
#        sumStates = getSumOfStates(Elist, self.energy, self.degeneracy, self.n0)
#        for n in range(1, len(Elist)):
#            self.assertAlmostEqual(numpy.sum(densStates[0:n+1]) / sumStates[n], 1.0, 3)

    def test_getDensityOfStates(self):
        """
        Test the getDensityOfStates() method.
        """
        Tlist = numpy.array([300,400,500,600])
        Elist = numpy.arange(0, 40000., 20.)
        for T in Tlist:
            densStates = getDensityOfStates(Elist, self.energy, self.degeneracy, self.n0)
            Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
            Qexp = getPartitionFunction(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 2, '{0} != {1} within 2 figures'.format(Qexp, Qact))
            
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
