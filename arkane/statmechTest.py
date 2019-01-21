#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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
This script contains unit tests of the :mod:`arkane.main` module.
"""

import unittest
import os

from rmgpy.exceptions import InputError

from arkane import Arkane
from arkane.statmech import StatMechJob, determine_rotor_symmetry
from arkane.qchem import QChemLog

################################################################################


class TestStatmech(unittest.TestCase):
    """
    Contains unit tests of the StatmechJob class.
    """
    @classmethod
    def setUp(self):
        arkane = Arkane()
        self.job_list = arkane.loadInputFile(os.path.join(os.path.dirname(os.path.abspath(__file__)),
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
        v_list1, angle = scan_log1.loadScanEnergies()
        v_list2, angle = scan_log2.loadScanEnergies()
        symmetry1 = determine_rotor_symmetry(energies=v_list1, label='NCC', pivots=[])
        symmetry2 = determine_rotor_symmetry(energies=v_list2, label='NCC', pivots=[])
        self.assertEqual(symmetry1, 1)
        self.assertEqual(symmetry2, 3)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
