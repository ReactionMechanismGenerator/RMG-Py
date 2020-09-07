#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains unit tests of the :mod:`arkane.ess.psi4` module.
"""

import os
import unittest

from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, HarmonicOscillator, HinderedRotor

from arkane.ess.psi4 import Psi4Log

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
        cls.data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'psi4')

    def test_number_of_atoms_from_psi4_log(self):
        """
        Uses a QChem log files to test that
        number of atoms can be properly read.
        """
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq.out'))
        self.assertEqual(log.get_number_of_atoms(), 3)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_ts.out'))
        self.assertEqual(log.get_number_of_atoms(), 4)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_dft.out'))
        self.assertEqual(log.get_number_of_atoms(), 3)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_dft_ts.out'))
        self.assertEqual(log.get_number_of_atoms(), 4)

    def test_energy_from_psi4_log(self):
        """
        Uses a QChem log files to test that
        molecular energies can be properly read.
        """
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq.out'))
        self.assertAlmostEqual(log.load_energy(), -199599899.9822719, delta=1e-2)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_ts.out'))
        self.assertAlmostEqual(log.load_energy(), -395828407.5987777, delta=1e-2)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_dft.out'))
        self.assertAlmostEqual(log.load_energy(), -200640009.37231186, delta=1e-2)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_dft_ts.out'))
        self.assertAlmostEqual(log.load_energy(), -397841662.56434655, delta=1e-2)

    def test_load_vibrations_from_psi4_log(self):
        """
        Uses a QChem log files to test that
        molecular energies can be properly read.
        """
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq.out'))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(len(conformer.modes[2]._frequencies.value), 3)
        self.assertEqual(conformer.modes[2]._frequencies.value[2], 4261.7445)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_dft_ts.out'))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(len(conformer.modes[2]._frequencies.value), 5)
        self.assertEqual(conformer.modes[2]._frequencies.value[2], 1456.2449)

    def test_load_modes_from_psi4_log(self):
        """
        Uses a psi4 log file for opt_freq.out to test that its
        molecular modes can be properly read.
        """
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq.out'))
        conformer, unscaled_frequencies = log.load_conformer()

        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 1)
        self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1)
        self.assertEqual(len(unscaled_frequencies), 3)

    def test_spin_multiplicity_from_psi4_log(self):
        """
        Uses a QChem log file for opt_freq_dft_ts.out to test that its
        molecular degrees of freedom can be properly read.
        """
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq.out'))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(conformer.spin_multiplicity, 1)
        log = Psi4Log(os.path.join(self.data_path, 'opt_freq_dft_ts.out'))
        conformer, unscaled_frequencies = log.load_conformer()
        self.assertEqual(conformer.spin_multiplicity, 1)


################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
