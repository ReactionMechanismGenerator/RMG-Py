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
This module contains unit tests of the :mod:`arkane.logs.gaussian` module.
"""

import os
import unittest

import numpy as np

import rmgpy.constants as constants
from external.wip import work_in_progress
from rmgpy.statmech import IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator, HinderedRotor

from arkane.logs.gaussian import GaussianLog
from arkane.statmech import determine_qm_software
from arkane.exceptions import LogError


################################################################################


class GaussianLogTest(unittest.TestCase):
    """
    Contains unit tests for the gaussian module, used for parsing Gaussian log files.
    """
    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        cls.data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'gaussian')

    @work_in_progress
    def test_load_ethylene_from_gaussian_log_cbsqb3(self):
        """
        Uses a Gaussian03 log file for ethylene (C2H4) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(self.data_path, 'ethylene.log'))
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

    def test_gaussian_energies(self):
        """
        test parsing double hydride, MP2, CCSD, CCSD(T) form Gaussian log
        """
        log_doublehybrid = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'B2PLYP.LOG'))
        log_mp2 = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'UMP2_C_ATOM.LOG'))
        log_ccsd = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'UCCSD_C_ATOM.LOG'))
        log_ccsdt = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'UCCSDT_C_ATOM.LOG'))
        log_qb3 = GaussianLog(os.path.join(os.path.dirname(__file__), '../examples/arkane/species/C2H5/', 'ethyl_cbsqb3.log'))

        self.assertAlmostEqual(log_doublehybrid.load_energy() / constants.Na / constants.E_h, -0.40217794572194e+02,
                               delta=1e-6)
        self.assertAlmostEqual(log_mp2.load_energy() / constants.Na / constants.E_h, -0.37504683723025e+02,
                               delta=1e-6)
        self.assertAlmostEqual(log_ccsd.load_energy() / constants.Na / constants.E_h, -37.517151426,
                               delta=1e-6)
        self.assertAlmostEqual(log_ccsdt.load_energy() / constants.Na / constants.E_h, -0.37517454469e+02,
                               delta=1e-6)
        self.assertAlmostEqual(log_qb3.load_energy() / constants.Na / constants.E_h, -79.029798,
                               delta=1e-6)

    def test_load_oxygen_from_gaussian_log(self):
        """
        Uses a Gaussian03 log file for oxygen (O2) to test that its
        molecular degrees of freedom can be properly read.
        """

        log = GaussianLog(os.path.join(self.data_path, 'oxygen.log'))
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

        log = GaussianLog(os.path.join(self.data_path, 'ethylene_G3.log'))
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

        log = GaussianLog(os.path.join(self.data_path, 'oxygen.log'))
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
        self.assertAlmostEqual(log._load_scan_angle(), 10.0)

    def test_load_number_scans(self):
        """
        Ensures proper scan angle found in Gaussian scan job
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'isobutanolQOOH_scan.log'))
        self.assertAlmostEqual(log._load_number_scans(), 36)

    def test_determine_qm_software(self):
        """
        Ensures that determine_qm_software returns a GaussianLog object
        """
        log = determine_qm_software(os.path.join(os.path.dirname(__file__), 'data', 'oxygen.log'))
        self.assertIsInstance(log, GaussianLog)

    def test_gaussian_log_error_termination(self):
        """
        Ensures that error termination gaussian log file raises an logError
        """
        file_path = os.path.join(os.path.dirname(__file__), 'data', 'error_termination.out')
        log = GaussianLog(file_path)
        with self.assertRaises(LogError) as log_error:
            log.load_conformer()
        self.assertTrue(f'The Gaussian job in {file_path} did not converge.' in str(log_error.exception))

    def test_load_scan_with_freq(self):
        """
        Ensures that the length of enegies with hr scans and freq calc is correct
        """
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'hr_scan_with_freq.log'))
        self.assertAlmostEqual(log._load_number_scans(), 36)
        self.assertAlmostEqual(log._load_scan_angle(), 10.0)
        vlist, _ = log.load_scan_energies()
        self.assertEqual(len(vlist), 37)

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
