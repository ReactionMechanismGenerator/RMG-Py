#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This script contains unit tests of the :mod:`rmgpy.kinetics.model` module.
"""

import unittest

from rmgpy.kinetics.model import get_reaction_order_from_rate_coefficient_units, \
                                 get_rate_coefficient_units_from_reaction_order

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
        self.assertEqual(0, get_reaction_order_from_rate_coefficient_units('mol/(m^3*s)'))
        self.assertEqual(0, get_reaction_order_from_rate_coefficient_units('mol/(cm^3*s)'))
        self.assertEqual(0, get_reaction_order_from_rate_coefficient_units('molecule/(m^3*s)'))
        self.assertEqual(0, get_reaction_order_from_rate_coefficient_units('molecule/(cm^3*s)'))
        
    def test_toOrder_first(self):
        """
        Test the conversion of first-order rate coefficient units to an integer
        reaction order.
        """
        self.assertEqual(1, get_reaction_order_from_rate_coefficient_units('s^-1'))
        
    def test_toOrder_second(self):
        """
        Test the conversion of second-order rate coefficient units to an integer
        reaction order.
        """
        self.assertEqual(2, get_reaction_order_from_rate_coefficient_units('m^3/(mol*s)'))
        self.assertEqual(2, get_reaction_order_from_rate_coefficient_units('cm^3/(mol*s)'))
        self.assertEqual(2, get_reaction_order_from_rate_coefficient_units('m^3/(molecule*s)'))
        self.assertEqual(2, get_reaction_order_from_rate_coefficient_units('cm^3/(molecule*s)'))
        
    def test_toOrder_third(self):
        """
        Test the conversion of third-order rate coefficient units to an integer
        reaction order.
        """
        self.assertEqual(3, get_reaction_order_from_rate_coefficient_units('m^6/(mol^2*s)'))
        self.assertEqual(3, get_reaction_order_from_rate_coefficient_units('cm^6/(mol^2*s)'))
        self.assertEqual(3, get_reaction_order_from_rate_coefficient_units('m^6/(molecule^2*s)'))
        self.assertEqual(3, get_reaction_order_from_rate_coefficient_units('cm^6/(molecule^2*s)'))
        
    def test_toUnits_zeroth(self):
        """
        Test the conversion of a reaction order of zero to rate coefficient
        units.
        """
        self.assertEqual('mol/(m^3*s)', get_rate_coefficient_units_from_reaction_order(0))
        
    def test_toUnits_first(self):
        """
        Test the conversion of a reaction order of one to rate coefficient
        units.
        """
        self.assertEqual('s^-1', get_rate_coefficient_units_from_reaction_order(1))
        
    def test_toUnits_second(self):
        """
        Test the conversion of a reaction order of two to rate coefficient
        units.
        """
        self.assertEqual('m^3/(mol*s)', get_rate_coefficient_units_from_reaction_order(2))
        
    def test_toUnits_third(self):
        """
        Test the conversion of a reaction order of three to rate coefficient
        units.
        """
        self.assertEqual('m^6/(mol^2*s)', get_rate_coefficient_units_from_reaction_order(3))
