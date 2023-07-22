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
This module contains unit tests of the :mod:`arkane.ess.qchem` module.
"""

import os
import unittest

from rmgpy.statmech import (
    IdealGasTranslation,
    LinearRotor,
    NonlinearRotor,
    HarmonicOscillator,
    HinderedRotor,
)

from arkane.ess.qchem import QChemLog
from arkane.exceptions import LogError

################################################################################


class QChemLogTest(unittest.TestCase):
    """
    Contains unit tests for the qchem module, used for parsing QChem log files.
    """

    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        cls.data_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "data", "qchem"
        )

    def test_check_for_errors(self):
        """
        Uses a QChem log file that reached to maximum number of optimization cycles
        to test if errors are properly parsed.
        """
        with self.assertRaises(LogError):
            QChemLog(os.path.join(self.data_path, "formyl_azide.out"))

    def test_number_of_atoms_from_qchem_log(self):
        """
        Uses QChem log files to test that
        number of atoms can be properly read.
        """
        log = QChemLog(os.path.join(self.data_path, "npropyl.out"))
        self.assertEqual(log.get_number_of_atoms(), 10)
        log = QChemLog(os.path.join(self.data_path, "co.out"))
        self.assertEqual(log.get_number_of_atoms(), 2)

    def test_energy_from_qchem_log(self):
        """
        Uses QChem log files to test that
        molecular energies can be properly read.
        """
        log = QChemLog(os.path.join(self.data_path, "npropyl.out"))
        self.assertAlmostEqual(log.load_energy(), -310896203.5432524, delta=1e-7)
        log = QChemLog(os.path.join(self.data_path, "co.out"))
        self.assertAlmostEqual(log.load_energy(), -297402545.0217114, delta=1e-7)
        log = QChemLog(os.path.join(self.data_path, "CH4_sp.out"))
        self.assertAlmostEqual(log.load_energy(), -106356735.53661588, delta=1e-7)

    def test_zero_point_energy_from_qchem_log(self):
        """
        Uses QChem log files to test that
        zero point energies can be properly read.
        """
        log = QChemLog(os.path.join(self.data_path, "npropyl.out"))
        self.assertAlmostEqual(log.load_zero_point_energy(), 228785.304, delta=1e-3)
        log = QChemLog(os.path.join(self.data_path, "co.out"))
        self.assertAlmostEqual(log.load_zero_point_energy(), 13476.664, delta=1e-3)
        log = QChemLog(
            os.path.join(self.data_path, "formyl_azide.out"), check_for_errors=False
        )
        self.assertAlmostEqual(log.load_zero_point_energy(), 83014.744, delta=1e-3)

    def test_load_vibrations_from_qchem_log(self):
        """
        Uses QChem log files to test that
        molecular energies can be properly read.
        """
        log = QChemLog(os.path.join(self.data_path, "npropyl.out"))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(len(conformer.modes[2]._frequencies.value), 24)
        self.assertEqual(conformer.modes[2]._frequencies.value[5], 881.79)
        log = QChemLog(os.path.join(self.data_path, "co.out"))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(len(conformer.modes[2]._frequencies.value), 1)
        self.assertEqual(conformer.modes[2]._frequencies.value, 2253.16)

    def test_load_npropyl_modes_from_qchem_log(self):
        """
        Uses a QChem log file for npropyl to test that its
        molecular modes can be properly read.
        """
        log = QChemLog(os.path.join(self.data_path, "npropyl.out"))
        conformer, unscaled_frequencies = log.load_conformer()

        self.assertTrue(
            len(
                [
                    mode
                    for mode in conformer.modes
                    if isinstance(mode, IdealGasTranslation)
                ]
            )
            == 1
        )
        self.assertTrue(
            len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)])
            == 1
        )
        self.assertTrue(
            len(
                [
                    mode
                    for mode in conformer.modes
                    if isinstance(mode, HarmonicOscillator)
                ]
            )
            == 1
        )
        self.assertTrue(
            len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)])
            == 0
        )

    def test_spin_multiplicity_from_qchem_log(self):
        """
        Uses QChem log files to test that
        molecular degrees of freedom can be properly read.
        """
        log = QChemLog(os.path.join(self.data_path, "npropyl.out"))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(conformer.spin_multiplicity, 2)
        log = QChemLog(os.path.join(self.data_path, "co.out"))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(conformer.spin_multiplicity, 1)

    def test_load_co_modes_from_qchem_log(self):
        """
        Uses a QChem log file for CO to test that its
        molecular degrees of freedom can be properly read.
        """
        log = QChemLog(os.path.join(self.data_path, "co.out"))
        conformer, unscaled_frequencies = log.load_conformer()
        E0 = log.load_energy()

        self.assertTrue(
            len(
                [
                    mode
                    for mode in conformer.modes
                    if isinstance(mode, IdealGasTranslation)
                ]
            )
            == 1
        )
        self.assertTrue(
            len([mode for mode in conformer.modes if isinstance(mode, LinearRotor)])
            == 1
        )
        self.assertTrue(
            len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)])
            == 0
        )
        self.assertTrue(
            len(
                [
                    mode
                    for mode in conformer.modes
                    if isinstance(mode, HarmonicOscillator)
                ]
            )
            == 1
        )
        self.assertTrue(
            len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)])
            == 0
        )

    def test_load_negative_frequency(self):
        """
        Load an imaginary frequency from a QChem output file.
        """
        log = QChemLog(os.path.join(self.data_path, "ts004630.log"))
        imaginary_freq = log.load_negative_frequency()
        self.assertEqual(imaginary_freq, -647.47)

        # verify that an error is raised if there are no negative frequencies
        with self.assertRaises(LogError):
            log = QChemLog(os.path.join(self.data_path, "npropyl.out"))
            imaginary_freq = log.load_negative_frequency()
