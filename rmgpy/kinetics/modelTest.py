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
This script contains unit tests of the :mod:`rmgpy.kinetics.model` module.
"""

import unittest

from rmgpy.kinetics.model import getReactionOrderFromRateCoefficientUnits, \
                                 getRateCoefficientUnitsFromReactionOrder
import rmgpy.constants as constants

################################################################################

class TestOrder(unittest.TestCase):
    """
    Contains unit tests of the functions for converting rate coefficient units
    to/from reaction orders.
    """
    
    def test_toOrder_zeroth(self):
        """
        Test the conversion of zeroth-order rate coefficient units to an integer
        reaction order.
        """
        self.assertEqual(0, getReactionOrderFromRateCoefficientUnits('mol/(m^3*s)'))
        self.assertEqual(0, getReactionOrderFromRateCoefficientUnits('mol/(cm^3*s)'))
        self.assertEqual(0, getReactionOrderFromRateCoefficientUnits('molecule/(m^3*s)'))
        self.assertEqual(0, getReactionOrderFromRateCoefficientUnits('molecule/(cm^3*s)'))
        
    def test_toOrder_first(self):
        """
        Test the conversion of first-order rate coefficient units to an integer
        reaction order.
        """
        self.assertEqual(1, getReactionOrderFromRateCoefficientUnits('s^-1'))
        
    def test_toOrder_second(self):
        """
        Test the conversion of second-order rate coefficient units to an integer
        reaction order.
        """
        self.assertEqual(2, getReactionOrderFromRateCoefficientUnits('m^3/(mol*s)'))
        self.assertEqual(2, getReactionOrderFromRateCoefficientUnits('cm^3/(mol*s)'))
        self.assertEqual(2, getReactionOrderFromRateCoefficientUnits('m^3/(molecule*s)'))
        self.assertEqual(2, getReactionOrderFromRateCoefficientUnits('cm^3/(molecule*s)'))
        
    def test_toOrder_third(self):
        """
        Test the conversion of third-order rate coefficient units to an integer
        reaction order.
        """
        self.assertEqual(3, getReactionOrderFromRateCoefficientUnits('m^6/(mol^2*s)'))
        self.assertEqual(3, getReactionOrderFromRateCoefficientUnits('cm^6/(mol^2*s)'))
        self.assertEqual(3, getReactionOrderFromRateCoefficientUnits('m^6/(molecule^2*s)'))
        self.assertEqual(3, getReactionOrderFromRateCoefficientUnits('cm^6/(molecule^2*s)'))
        
    def test_toUnits_zeroth(self):
        """
        Test the conversion of a reaction order of zero to rate coefficient
        units.
        """
        self.assertEqual('mol/(m^3*s)', getRateCoefficientUnitsFromReactionOrder(0))
        
    def test_toUnits_first(self):
        """
        Test the conversion of a reaction order of one to rate coefficient
        units.
        """
        self.assertEqual('s^-1', getRateCoefficientUnitsFromReactionOrder(1))
        
    def test_toUnits_second(self):
        """
        Test the conversion of a reaction order of two to rate coefficient
        units.
        """
        self.assertEqual('m^3/(mol*s)', getRateCoefficientUnitsFromReactionOrder(2))
        
    def test_toUnits_third(self):
        """
        Test the conversion of a reaction order of three to rate coefficient
        units.
        """
        self.assertEqual('m^6/(mol^2*s)', getRateCoefficientUnitsFromReactionOrder(3))
