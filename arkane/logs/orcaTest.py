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
This module contains unit tests of the :mod:`arkane.logs.orca` module.
"""

import os
import unittest

from arkane.logs.orca import OrcaLog

################################################################################


class OrcaTest(unittest.TestCase):
    """
    Contains unit tests for the orca module, used for parsing Orca log files.
    """

    def test_number_of_atoms_from_orca_log(self):
        """
        Uses a Orca log files to test that
        number of atoms can be properly read.
        """
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_opt_freq_test.log'))
        self.assertEqual(log.get_number_of_atoms(), 3)
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_dlpno_test.log'))
        self.assertEqual(log.get_number_of_atoms(), 30)

    def test_read_coordinates_from_orca_log(self):
        """
        Uses a Orca log files to test that
        coordinate can be properly read.
        """
        log1 = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_opt_freq_test.log'))
        coord, number, mass = log1.load_geometry()
        self.assertEqual(len(coord), 3)
        log2 = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_dlpno_test.log'))
        coord, number, mass = log2.load_geometry()
        self.assertEqual(len(coord), 30)

    def test_energy_from_orca_log(self):
        """
        Uses a Orca log files to test that
        molecular energies can be properly read.
        """
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_opt_freq_test.log'))
        self.assertAlmostEqual(log.load_energy(), -200656222.56292167, delta=1e-3)
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_TS_test.log'))
        self.assertAlmostEqual(log.load_energy(), -88913623.10592663, delta=1e-3)
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_dlpno_test.log'))
        self.assertAlmostEqual(log.load_energy(), -1816470909.1153, delta=1e-3)

    def test_load_zero_point_energy_from_orca_log(self):
        """
        Uses a Orca log files to test that
        molecular zero point_energy can be properly read.
        """
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_opt_freq_test.log'))
        self.assertAlmostEqual(log.load_zero_point_energy(), 55502.673180815, delta=1e-3)
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_TS_test.log'))
        self.assertAlmostEqual(log.load_zero_point_energy(), 93500.08860598055, delta=1e-3)

    def test_load_negative_frequency_from_orca_log(self):
        """
        Uses a orca log file for npropyl to test that its
        negative frequency can be properly read.
        """
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_TS_test.log'))
        self.assertAlmostEqual(log.load_negative_frequency(), -503.24, delta=1e-1)

    def test_T1_diagnostic_from_orca_log(self):
        """
        Uses a Orca log file for npropyl to test that its
        T1_diagnostic of freedom can be properly read.
        """
        log = OrcaLog(os.path.join(os.path.dirname(__file__), 'data', 'Orca_dlpno_test.log'))
        self.assertAlmostEqual(log.get_T1_diagnostic(), 0.009872238, delta=1e-3)

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
