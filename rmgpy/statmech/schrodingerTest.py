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
This script contains unit tests of the :mod:`rmgpy.statmech.schrodinger`
module.
"""

from __future__ import division

import unittest
import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech.schrodinger import getPartitionFunction, getHeatCapacity, getEnthalpy, getEntropy, \
    getDensityOfStates

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
        self.energy = lambda j: self.B * j * (j + 1)
        self.degeneracy = lambda j: 2 * j + 1
        self.n0 = 0

    def test_get_partition_function(self):
        """
        Test the getPartitionFunction() method.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([208.8907, 347.9285, 695.5234, 1043.118, 1390.713])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = getPartitionFunction(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(q_exp / q_act, 1.0, 4,
                                   '{0} != {1} within 4 figures'.format(q_exp, q_act))

    def test_get_heat_capacity(self):
        """
        Test the getHeatCapacity() method.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1, 1, 1, 1, 1])
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = getHeatCapacity(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(cv_exp / cv_act, 1.0, 4,
                                   '{0} != {1} within 4 figures'.format(cv_exp, cv_act))

    def test_get_enthalpy(self):
        """
        Test the getEnthalpy() method.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([0.9984012, 0.9990409, 0.9995205, 0.9996803, 0.9997603])
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = getEnthalpy(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(h_exp / h_act, 1.0, 4,
                                   '{0} != {1} within 4 figures'.format(h_exp, h_act))

    def test_get_entropy(self):
        """
        Test the getEntropy() method.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([6.340212, 6.851038, 7.544185, 7.949650, 8.237332])
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = getEntropy(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(s_exp / s_act, 1.0, 4,
                                   '{0} != {1} within 4 figures'.format(s_exp, s_act))

#    def test_get_sum_of_states(self):
#        """
#        Test the getSumOfStates() method.
#        """
#        e_list = np.arange(0, 10., 0.01)
#        dens_states = getDensityOfStates(e_list, self.energy, self.degeneracy, self.n0)
#        sum_states = getSumOfStates(e_list, self.energy, self.degeneracy, self.n0)
#        for n in range(1, len(e_list)):
#            self.assertAlmostEqual(np.sum(dens_states[0:n+1]) / sum_states[n], 1.0, 3)

    def test_get_density_of_states(self):
        """
        Test the getDensityOfStates() method.
        """
        t_list = np.array([300, 400, 500, 600])
        e_list = np.arange(0, 40000., 20.)
        for T in t_list:
            dens_states = getDensityOfStates(e_list, self.energy, self.degeneracy, self.n0)
            q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
            q_exp = getPartitionFunction(T, self.energy, self.degeneracy, self.n0)
            self.assertAlmostEqual(q_exp / q_act, 1.0, 2,
                                   '{0} != {1} within 2 figures'.format(q_exp, q_act))

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
