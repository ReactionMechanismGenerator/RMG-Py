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
This module contains unit tests of the :mod:`arkane.statmech` module.
"""

import os
import unittest

import numpy as np

from rmgpy.exceptions import InputError

from arkane import Arkane
from arkane.qchem import QChemLog
from arkane.statmech import StatMechJob, determine_rotor_symmetry, is_linear

################################################################################


class TestStatmech(unittest.TestCase):
    """
    Contains unit tests of the StatmechJob class.
    """

    @classmethod
    def setUp(cls):
        """A method that is run before each unit test in this class"""
        arkane = Arkane()
        cls.job_list = arkane.load_input_file(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                           'data', 'Benzyl', 'input.py'))

    def test_gaussian_log_file_error(self):
        """Test that the proper error is raised if gaussian geometry and frequency file paths are not the same"""
        job = self.job_list[-2]
        self.assertTrue(isinstance(job, StatMechJob))
        with self.assertRaises(InputError):
            job.load()

    def test_rotor_symmetry_determination(self):
        """
        Test that the correct symmetry number is determined for rotor potential scans.
        """
        path1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'NCC_NRotor.out')
        path2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'NCC_CRotor.out')
        scan_log1 = QChemLog(path1)
        scan_log2 = QChemLog(path2)
        v_list1, angle = scan_log1.load_scan_energies()
        v_list2, angle = scan_log2.load_scan_energies()
        symmetry1 = determine_rotor_symmetry(energies=v_list1, label='NCC', pivots=[])
        symmetry2 = determine_rotor_symmetry(energies=v_list2, label='NCC', pivots=[])
        self.assertEqual(symmetry1, 1)
        self.assertEqual(symmetry2, 3)

    def test_is_linear(self):
        """Test that we can determine the linearity of a molecule from it's coordinates"""
        xyz1 = np.array([
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.159076],
            [0.000000, 0.000000, -1.159076]])  # a trivial case
        xyz2 = np.array([
            [-0.06618943, -0.12360663, -0.07631983],
            [-0.79539707, 0.86755487, 1.02675668],
            [-0.68919931, 0.25421823, -1.34830853],
            [0.01546439, -1.54297548, 0.44580391],
            [1.94428095, 0.40772394, 1.03719428],
            [2.20318015, -0.14715186, -0.64755729],
            [1.59252246, 1.51178950, -0.33908352],
            [-0.87856890, -2.02453514, 0.38494433],
            [-1.34135876, 1.49608206, 0.53295071]])  # a non-linear multi-atom molecule
        xyz3 = np.array([
            [0.0000000000, 0.0000000000, 0.3146069129],
            [-1.0906813653, 0.0000000000, -0.1376405244],
            [1.0906813653, 0.0000000000, -0.1376405244]])  # NO2, a non-linear 3-atom molecule
        xyz4 = np.array([
            [0.0000000000, 0.0000000000, 0.1413439534],
            [-0.8031792912, 0.0000000000, -0.4947038368],
            [0.8031792912, 0.0000000000, -0.4947038368]])  # NH2, a non-linear 3-atom molecule
        xyz5 = np.array([
            [-0.5417345330, 0.8208150346, 0.0000000000],
            [0.9206183692, 1.6432038228, 0.0000000000],
            [-1.2739176462, 1.9692549926, 0.0000000000]])  # HSO, a non-linear 3-atom molecule
        xyz6 = np.array([
            [1.18784533, 0.98526702, 0.00000000],
            [0.04124533, 0.98526702, 0.00000000],
            [-1.02875467, 0.98526702, 0.00000000]])  # HCN, a linear 3-atom molecule
        xyz7 = np.array([
            [-4.02394116, 0.56169428, 0.00000000],
            [-5.09394116, 0.56169428, 0.00000000],
            [-2.82274116, 0.56169428, 0.00000000],
            [-1.75274116, 0.56169428, 0.00000000]])  # C2H2, a linear 4-atom molecule
        xyz8 = np.array([
            [-1.02600933, 2.12845307, 0.00000000],
            [-0.77966935, 0.95278385, 0.00000000],
            [-1.23666197, 3.17751246, 0.00000000],
            [-0.56023545, -0.09447399, 0.00000000]])  # C2H2, just 0.5 degree off from linearity, so NOT linear
        xyz9 = np.array([
            [-1.1998, 0.1610, 0.0275],
            [-1.4021, 0.6223, -0.8489],
            [-1.48302, 0.80682, -1.19946]])  # just 3 points in space on a straight line (not a physical molecule)
        xyz10 = np.array([
            [-1.1998, 0.1610, 0.0275]])  # mono-atomic species, non-linear
        xyz11 = np.array([
            [1.06026500, -0.07706800, 0.03372800],
            [3.37340700, -0.07706800, 0.03372800],
            [2.21683600, -0.07706800, 0.03372800]])  # CO2 at wb97xd/6-311+g(d,p), linear
        xyz12 = np.array([
            [1.05503600, -0.00335000, 0.09823600],
            [2.42816800, -0.00335000, 0.09823600],
            [-0.14726400, -0.00335000, 0.09823600],
            [3.63046800, -0.00335000, 0.09823600],
            [-1.21103500, -0.00335000, 0.09823600],
            [4.69423900, -0.00335000, 0.09823600]])  # C#CC#C at wb97xd/6-311+g(d,p), linear

        self.assertTrue(is_linear(xyz1))
        self.assertTrue(is_linear(xyz6))
        self.assertTrue(is_linear(xyz7))
        self.assertTrue(is_linear(xyz9))
        self.assertTrue(is_linear(xyz11))
        self.assertTrue(is_linear(xyz12))
        self.assertFalse(is_linear(xyz2))
        self.assertFalse(is_linear(xyz3))
        self.assertFalse(is_linear(xyz4))
        self.assertFalse(is_linear(xyz5))
        self.assertFalse(is_linear(xyz8))
        self.assertFalse(is_linear(xyz10))


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
