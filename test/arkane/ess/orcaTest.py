#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains unit tests of the :mod:`arkane.ess.orca` module.
"""

import os

import numpy as np

from arkane.ess.orca import OrcaLog
from arkane.exceptions import LogError
from rmgpy.statmech import IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator
import pytest


class OrcaTest:
    """
    Contains unit tests for the orca module, used for parsing Orca log files.
    """

    @classmethod
    def setup_class(cls):
        """
        A method that is run before all unit tests in this class.
        """
        cls.data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "..", "..", "arkane", "data", "orca")

    def test_check_for_errors(self):
        """
        Uses an Orca log file with SCF error
        to test if errors are properly parsed.
        """
        with pytest.raises(LogError):
            OrcaLog(os.path.join(self.data_path, "Orca_scf_error.log"))

    def test_number_of_atoms_from_orca_log(self):
        """
        Uses Orca log files to test that
        number of atoms can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "N_MRCI.log"))
        assert log.get_number_of_atoms() == 1
        log = OrcaLog(os.path.join(self.data_path, "Orca_opt_freq_test.log"))
        assert log.get_number_of_atoms() == 3
        log = OrcaLog(os.path.join(self.data_path, "Orca_dlpno_test.log"))
        assert log.get_number_of_atoms() == 30

    def test_read_coordinates_from_orca_log(self):
        """
        Uses Orca log files to test that
        coordinate can be properly read.
        """
        log1 = OrcaLog(os.path.join(self.data_path, "Orca_opt_freq_test.log"))
        coord, number, mass = log1.load_geometry()
        assert len(coord) == 3
        log2 = OrcaLog(os.path.join(self.data_path, "Orca_dlpno_test.log"))
        coord, number, mass = log2.load_geometry()
        assert len(coord) == 30

    def test_energy_from_orca_log(self):
        """
        Uses Orca log files to test that
        molecular energies can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "Orca_opt_freq_test.log"))
        assert abs(log.load_energy() - -200656222.56292167) < 1e-3
        log = OrcaLog(os.path.join(self.data_path, "Orca_TS_test.log"))
        assert abs(log.load_energy() - -88913623.10592663) < 1e-3
        log = OrcaLog(os.path.join(self.data_path, "Orca_dlpno_test.log"))
        assert abs(log.load_energy() - -1816470909.1153) < 1e-3

    def test_load_zero_point_energy_from_orca_log(self):
        """
        Uses Orca log files to test that
        molecular zero point_energy can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "Orca_opt_freq_test.log"))
        assert abs(log.load_zero_point_energy() - 55502.673180815) < 1e-3
        log = OrcaLog(os.path.join(self.data_path, "Orca_TS_test.log"))
        assert abs(log.load_zero_point_energy() - 93500.08860598055) < 1e-3

    def test_load_negative_frequency_from_orca_log(self):
        """
        Uses an Orca log file for npropyl to test that its
        negative frequency can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "Orca_TS_test.log"))
        assert abs(log.load_negative_frequency() - -503.24) < 1e-1

    def test_T1_diagnostic_from_orca_log(self):
        """
        Uses an Orca log file for a DLPNO calculation to test that its
        T1_diagnostic of freedom can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "Orca_dlpno_test.log"))
        assert abs(log.get_T1_diagnostic() - 0.009872238) < 1e-3

    def test_load_vibrations_from_orca_log(self):
        """
        Test that molecular energies can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "Orca_opt_freq_test.log"))
        conformer, unscaled_frequencies = log.load_conformer()
        assert len(conformer.modes[2]._frequencies.value) == 3
        assert conformer.modes[2]._frequencies.value[1] == 3780.96
        log = OrcaLog(os.path.join(self.data_path, "Orca_TS_test.log"))
        conformer, unscaled_frequencies = log.load_conformer()
        assert len(conformer.modes[2]._frequencies.value) == 11
        assert conformer.modes[2]._frequencies.value[2] == 331.23
        log = OrcaLog(os.path.join(self.data_path, "freq_orca_5.0.4.log"))
        conformer, unscaled_frequencies = log.load_conformer()
        assert len(conformer.modes[2]._frequencies.value) == 6
        assert conformer.modes[2]._frequencies.value[0] == 228.25

    def test_load_modes_from_orca_log(self):
        """
        Test that molecular modes can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "Orca_opt_freq_test.log"))
        conformer, unscaled_frequencies = log.load_conformer()
        assert len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1
        assert len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 1
        assert len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1
        assert len(unscaled_frequencies) == 3

        log = OrcaLog(os.path.join(self.data_path, "N_MRCI.log"))
        conformer, unscaled_frequencies = log.load_conformer()
        assert len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1
        assert len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 0
        assert len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 0
        assert len(unscaled_frequencies) == 0

    def test_spin_multiplicity_from_orca_log(self):
        """
        Test that molecular degrees of freedom can be properly read.
        """
        log = OrcaLog(os.path.join(self.data_path, "Orca_opt_freq_test.log"))
        conformer, unscaled_frequencies = log.load_conformer()
        assert conformer.spin_multiplicity == 1
        log = OrcaLog(os.path.join(self.data_path, "Orca_TS_test.log"))
        conformer, unscaled_frequencies = log.load_conformer()
        assert conformer.spin_multiplicity == 1

    def test_linear_rotor_from_oblique_diatomic(self, tmp_path):
        """
        A diatomic should get a LinearRotor regardless of its orientation,
        even when the smallest principal moment is a floating-point residual.
        """
        oblique = np.array([[0.803286, 1.531929, -1.098384], [0.663152, 2.004019, -0.249782]])
        bond = float(np.linalg.norm(oblique[1] - oblique[0]))
        aligned = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, bond]])
        q_rot = []
        for name, coords in (("oblique", oblique), ("aligned", aligned)):
            atoms = "\n".join(
                f"  {symbol}  " + "  ".join(f"{x:.17g}" for x in position)
                for symbol, position in zip(("O", "H"), coords)
            )
            path = tmp_path / f"oh_{name}.log"
            path.write_text(
                "CARTESIAN COORDINATES (ANGSTROEM)\n"
                "---------------------------------\n"
                f"{atoms}\n\n"
                " Multiplicity           Mult            ....    2\n"
                " Mode    freq (cm**-1)   T**2         TX         TY         TZ\n"
                "---------------------------------------------------------\n"
                "   5:      3669.30   10.000000  (  0.000000   0.000000   1.000000)\n\n"
                "                             ****ORCA TERMINATED NORMALLY****\n"
            )
            log = OrcaLog(str(path))
            conformer, _ = log.load_conformer(symmetry=1, optical_isomers=1)
            rotors = [mode for mode in conformer.modes if isinstance(mode, (LinearRotor, NonlinearRotor))]
            assert len(rotors) == 1
            assert isinstance(rotors[0], LinearRotor), name
            q_rot.append(rotors[0].get_partition_function(340.3))
        assert q_rot[0] == pytest.approx(q_rot[1])
        assert q_rot[0] == pytest.approx(12.805036, rel=1e-5)
