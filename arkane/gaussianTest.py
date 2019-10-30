#!/usr/bin/env python3

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
This module contains unit tests of the :mod:`arkane.gaussian` module.
"""

import os
import unittest

import numpy as np

import rmgpy.constants as constants
from external.wip import work_in_progress
from rmgpy.statmech import IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator, HinderedRotor

from arkane.gaussian import GaussianLog
from arkane.statmech import determine_qm_software

################################################################################


class GaussianTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.io.gaussian module, used for reading
    and writing Gaussian files.
    """

    @work_in_progress
    def test_load_ethylene_from_gaussian_log_cbsqb3(self):
        """
        Uses a Gaussian03 log file for ethylene (C2H4) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'ethylene.log'))
        conformer, unscaled_frequencies = log.load_conformer()
        e0 = log.load_energy()

        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)]) == 0)

        trans = [mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)][0]
        rot = [mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)][0]
        vib = [mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)][0]
        t_list = np.array([298.15], np.float64)
        self.assertAlmostEqual(trans.get_partition_function(t_list), 5.83338e6, delta=1e1)
        self.assertAlmostEqual(rot.get_partition_function(t_list), 2.59622e3, delta=1e-2)
        self.assertAlmostEqual(vib.get_partition_function(t_list), 1.0481e0, delta=1e-4)

        self.assertAlmostEqual(e0 / constants.Na / constants.E_h, -78.467452, 4)
        self.assertEqual(conformer.spin_multiplicity, 1)
        self.assertEqual(conformer.optical_isomers, 1)

    def test_load_oxygen_from_gaussian_log(self):
        """
        Uses a Gaussian03 log file for oxygen (O2) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'oxygen.log'))
        conformer, unscaled_frequencies = log.load_conformer()
        e0 = log.load_energy()

        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, LinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)]) == 0)

        trans = [mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)][0]
        rot = [mode for mode in conformer.modes if isinstance(mode, LinearRotor)][0]
        vib = [mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)][0]
        t_list = np.array([298.15], np.float64)
        self.assertAlmostEqual(trans.get_partition_function(t_list), 7.11169e6, delta=1e1)
        self.assertAlmostEqual(rot.get_partition_function(t_list), 7.13316e1, delta=1e-4)
        self.assertAlmostEqual(vib.get_partition_function(t_list), 1.00037e0, delta=1e-4)

        self.assertAlmostEqual(e0 / constants.Na / constants.E_h, -150.3784877, 4)
        self.assertEqual(conformer.spin_multiplicity, 3)
        self.assertEqual(conformer.optical_isomers, 1)

    @work_in_progress
    def test_load_ethylene_from_gaussian_log_g3(self):
        """
        Uses a Gaussian03 log file for ethylene (C2H4) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'ethylene_G3.log'))
        conformer, unscaled_frequencies = log.load_conformer()
        e0 = log.load_energy()

        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)]) == 0)

        trans = [mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)][0]
        rot = [mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)][0]
        vib = [mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)][0]
        t_list = np.array([298.15], np.float64)

        self.assertAlmostEqual(trans.get_partition_function(t_list), 5.83338e6, delta=1e1)
        self.assertAlmostEqual(rot.get_partition_function(t_list), 2.53410e3, delta=1e-2)
        self.assertAlmostEqual(vib.get_partition_function(t_list), 1.0304e0, delta=1e-4)

        self.assertAlmostEqual(e0 / constants.Na / constants.E_h, -78.562189, 4)
        self.assertEqual(conformer.spin_multiplicity, 1)
        self.assertEqual(conformer.optical_isomers, 1)

    def test_load_symmetry_and_optics(self):
        """
        Uses a Gaussian03 log file for oxygen (O2) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'oxygen.log'))
        optical, symmetry, _ = log.get_symmetry_properties()
        self.assertEqual(optical, 1)
        self.assertEqual(symmetry, 2)

        conf = log.load_conformer()[0]
        self.assertEqual(conf.optical_isomers, 1)
        found_rotor = False
        for mode in conf.modes:
            if isinstance(mode, LinearRotor):
                self.assertEqual(mode.symmetry, 2)
                found_rotor = True
        self.assertTrue(found_rotor)

    def test_load_scan_angle(self):
        """
        Ensures proper scan angle found in Gaussian scan job
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'isobutanolQOOH_scan.log'))
        self.assertAlmostEqual(log.load_scan_angle(), 10.0)

    def test_load_number_scans(self):
        """
        Ensures proper scan angle found in Gaussian scan job
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'isobutanolQOOH_scan.log'))
        self.assertAlmostEqual(log.load_number_scans(), 36)

    def test_determine_qm_software(self):
        """
        Ensures that determine_qm_software returns a GaussianLog object
        """
        log = determine_qm_software(os.path.join(os.path.dirname(__file__), 'data', 'oxygen.log'))
        self.assertIsInstance(log, GaussianLog)

    def test_gap_in_scan(self):
        """
        Ensures when an error occurs in the hindered rotors, proper distribution occurs
        """
        # load gaussian log with 10 degree separations, which has a failed optimization
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'isobutanolQOOH_scan.log'))
        with self.assertLogs(level=30):  # warnings only
            vlist, angles = log.load_scan_energies()
        self.assertAlmostEqual(angles[1], 10. * np.pi / 180)
        self.assertAlmostEqual(angles[-1], 2 * np.pi)

    def test_uncompleted_runs_throw_warning(self):
        """
        Ensures a warning is thrown when the gaussian scan file is not complete
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'alcohol_ether_scan.log'))
        with self.assertLogs(level=30):  # warnings only
            log.load_scan_energies()

    def test_rigid_scan_pivots(self):
        """
        Tests that load_scan_pivot_atoms returns correct value for rigid scans
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'frozen_scan.log'))
        pivots = log.load_scan_pivot_atoms()
        self.assertEqual(pivots[0], 13)
        self.assertEqual(pivots[1], 12)
        self.assertEqual(pivots[2], 1)
        self.assertEqual(pivots[3], 4)

    def test_rigid_scan_frozen(self):
        """
        Tests that load_scan_frozen_atoms returns 'rigid scan' for rigid scans
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'frozen_scan.log'))
        pivots = log.load_scan_frozen_atoms()
        self.assertEqual(pivots[0], 'rigid scan')

    def test_load_scan_angle_on_rigid(self):
        """
        Ensures proper scan angle found in Gaussian rigid scan job
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'frozen_scan.log'))
        self.assertAlmostEqual(log.load_scan_angle(), 10.0)

    def test_load_number_scans_on_rigid(self):
        """
        Ensures proper scan angle found in Gaussian rigid scan job
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'frozen_scan.log'))
        self.assertAlmostEqual(log.load_number_scans(), 36)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
