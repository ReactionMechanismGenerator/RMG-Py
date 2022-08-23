#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains unit tests of the :mod:`arkane.ess.molpro` module.
"""

import os
import unittest

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, HarmonicOscillator, HinderedRotor

from arkane.ess.molpro import MolproLog
from arkane.exceptions import LogError

################################################################################


class MolproLogTest(unittest.TestCase):
    """
    Contains unit tests for the molpro module, used for parsing Molpro log files.
    """
    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        cls.data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'molpro')

    def test_check_for_errors(self):
        """
        Uses Molpro log files that failed due to insufficient memory and unrecognized basis set
        to test if errors are properly parsed.
        """
        with self.assertRaises(LogError):
            MolproLog(os.path.join(self.data_path, 'insufficient_memory.out'))
        with self.assertRaises(LogError):
            MolproLog(os.path.join(self.data_path, 'unrecognized_basis_set.out'))

    def test_number_of_atoms_from_molpro_log(self):
        """
        Uses a Molpro log file for ethylene_dz (C2H4) to test that the
        number of atoms can be properly read.
        """
        log = MolproLog(os.path.join(self.data_path, 'ethylene_f12_dz.out'))
        self.assertEqual(log.get_number_of_atoms(), 6)

    def test_load_dz_from_molpro_log_f12(self):
        """
        Uses a Molpro log file for ethylene_dz (C2H4) to test that F12a
        energy can be properly read.
        """
        log = MolproLog(os.path.join(self.data_path, 'ethylene_f12_dz.out'))
        e0 = log.load_energy()

        self.assertAlmostEqual(e0 / constants.Na / constants.E_h, -78.474353559604, 5)

    def test_load_qz_from_molpro_log_f12(self):
        """
        Uses a Molpro log file for ethylene_qz (C2H4) to test that F12b
        energy can be properly read.
        """
        log = MolproLog(os.path.join(self.data_path, 'ethylene_f12_qz.out'))
        e0 = log.load_energy()

        self.assertAlmostEqual(e0 / constants.Na / constants.E_h, -78.472682755635, 5)

    def test_load_rad_from_molpro_log_f12(self):
        """
        Uses a Molpro log file for OH (C2H4) to test that radical
        energy can be properly read.
        """
        log = MolproLog(os.path.join(self.data_path, 'OH_f12.out'))
        e0 = log.load_energy()

        self.assertAlmostEqual(e0 / constants.Na / constants.E_h, -75.663696424380, 5)

    def test_load_hosi_from_molpro_log(self):
        """
        Uses a molpro log file for HOSI to test that its
        molecular degrees of freedom can be properly read.
        """
        log = MolproLog(os.path.join(self.data_path, 'HOSI_ccsd_t1.out'))
        conformer, unscaled_frequencies = log.load_conformer(spin_multiplicity=1)
        e0 = log.load_energy()

        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)]) == 0)

        trans = [mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)][0]
        rot = [mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)][0]
        vib = [mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)][0]
        t_list = np.array([298.15], np.float64)

        self.assertAlmostEqual(trans.get_partition_function(t_list), 9.175364e7, delta=1e1)
        self.assertAlmostEqual(rot.get_partition_function(t_list), 1.00005557e5, delta=1e-2)
        self.assertAlmostEqual(vib.get_partition_function(t_list), 1.9734989e0, delta=1e-4)

        self.assertAlmostEqual(e0 / constants.Na / constants.E_h, -768.275662, 4)
        self.assertEqual(conformer.spin_multiplicity, 1)
        self.assertEqual(conformer.optical_isomers, 1)

    def test_load_non_f12_e0(self):
        """
        Load E0 for CCSD(T) (without F12) from a molpro output file
        """
        molpro_log = MolproLog(os.path.join(self.data_path, 'TS_CCSD(T)_no_F12_sp_molpro.out'))
        e0 = molpro_log.load_energy()
        self.assertAlmostEqual(e0, -301585968.58196217, places=7)

    def test_load_mrci_e0(self):
        """
        Load the MRCI and MRCI+Davidson energies from a molpro output file
        """
        mrci_log = MolproLog(os.path.join(self.data_path, 'molpro_mrci.out'))
        mrciq_log = MolproLog(os.path.join(self.data_path, 'molpro_mrci+q.out'))
        mrci_e0 = mrci_log.load_energy()
        mrciq_e0 = mrciq_log.load_energy()
        self.assertAlmostEqual(mrci_e0, -293217091.0381712, places=7)
        self.assertAlmostEqual(mrciq_e0, -293284017.3925107, places=7)

    def test_load_negative_frequency(self):
        """
        Load an imaginary frequency from a molpro output file
        """
        freq_log = MolproLog(os.path.join(self.data_path, 'molpro_TS.out'))
        imaginary_freq = freq_log.load_negative_frequency()
        self.assertEqual(imaginary_freq, -1997.98)

    def test_get_d1_diagnostic(self):
        """
        Ensure molpro can retrieve the T1 diagnostic from CCSD calculations
        """
        log = MolproLog(os.path.join(self.data_path, 'ethylene_f12_dz.out'))
        d1_diagnostic = log.get_D1_diagnostic()
        self.assertAlmostEqual(d1_diagnostic, 0.03369031)

    def test_get_t1_diagnostic(self):
        """
        Ensure molpro can retrieve the T1 diagnostic from CCSD calculations
        """
        log = MolproLog(os.path.join(self.data_path, 'ethylene_f12_dz.out'))
        t1_diagnostic = log.get_T1_diagnostic()
        self.assertAlmostEqual(t1_diagnostic, 0.01152184)

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
