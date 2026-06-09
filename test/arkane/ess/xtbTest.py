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

"""Tests for the xTB ESS adapter."""

import os

import numpy as np
import pytest

import rmgpy.constants as constants
from arkane.exceptions import LogError
from arkane.ess.factory import ess_factory
from arkane.ess.xtb import XTBLog

DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), '..', '..', 'arkane', 'data', 'xTB')


class TestFactory:

    def test_detects_freq(self):
        assert isinstance(ess_factory(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False), XTBLog)

    def test_detects_opt(self):
        assert isinstance(ess_factory(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False), XTBLog)

    def test_detects_co2(self):
        assert isinstance(ess_factory(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False), XTBLog)

    def test_detects_ts(self):
        assert isinstance(ess_factory(os.path.join(DATA_PATH, 'TS_NH2+N2H3_xtb.out'), check_for_errors=False), XTBLog)

    def test_detects_sp(self):
        assert isinstance(ess_factory(os.path.join(DATA_PATH, 'NCC_xTB.out'), check_for_errors=False), XTBLog)


class TestGetNumberOfAtoms:

    def test_freq_c3h7(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        assert log.get_number_of_atoms() == 10

    def test_opt_c3h7(self):
        log = XTBLog(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False)
        assert log.get_number_of_atoms() == 10

    def test_co2(self):
        """CO2 opt+freq output (Turbomol $coord format)."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        assert log.get_number_of_atoms() == 3

    def test_opt2_c3h11(self):
        log = XTBLog(os.path.join(DATA_PATH, 'xtb_opt_2.out'), check_for_errors=False)
        assert log.get_number_of_atoms() == 11


class TestGeometryV2000:

    def test_opt_c3h7_shape(self):
        log = XTBLog(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False)
        coord, number, mass = log.load_geometry()
        assert coord.shape == (10, 3)
        assert len(number) == 10

    def test_opt_c3h7_atoms(self):
        log = XTBLog(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False)
        _, number, _ = log.load_geometry()
        assert list(number[:3]) == [6, 6, 6]
        assert all(n == 1 for n in number[3:])

    def test_opt_c3h7_coords(self):
        log = XTBLog(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False)
        coord, _, _ = log.load_geometry()
        np.testing.assert_allclose(coord[0], [-1.2948, 0.0509, -0.1043], atol=1e-4)

    def test_opt2_v2000(self):
        """xtb_opt_2.out uses V2000 format with 11 atoms."""
        log = XTBLog(os.path.join(DATA_PATH, 'xtb_opt_2.out'), check_for_errors=False)
        coord, number, mass = log.load_geometry()
        assert coord.shape == (11, 3)
        assert sum(number == 6) == 3  # 3 carbons
        assert sum(number == 1) == 8  # 8 hydrogens


class TestGeometryTurbomol:

    def test_co2_shape(self):
        """CO2 opt output uses $coord (Bohr) format."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        coord, number, mass = log.load_geometry()
        assert coord.shape == (3, 3)

    def test_co2_atoms(self):
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        _, number, _ = log.load_geometry()
        assert list(number) == [8, 6, 8]  # O, C, O

    def test_co2_linearity(self):
        """CO2 should be linear — central C at ~origin, O atoms symmetric."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        coord, _, _ = log.load_geometry()
        # Central carbon should be near origin
        assert abs(coord[1, 0]) < 0.01
        # O atoms should be roughly symmetric about C
        np.testing.assert_allclose(coord[0, 0], -coord[2, 0], atol=0.01)

    def test_opt1_turbomol(self):
        """xtb_opt_1.out uses $coord Turbomol format with 8 atoms."""
        log = XTBLog(os.path.join(DATA_PATH, 'xtb_opt_1.out'), check_for_errors=False)
        coord, number, mass = log.load_geometry()
        assert coord.shape == (8, 3)
        assert sum(number == 6) == 3  # 3 C
        assert sum(number == 8) == 1  # 1 O
        assert sum(number == 1) == 4  # 4 H

    def test_opt1_units_are_angstrom(self):
        """Verify Bohr→Angstrom conversion: C-O bond ~1.14 Å, not ~2.16 Bohr."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        coord, _, _ = log.load_geometry()
        co_dist = np.linalg.norm(coord[0] - coord[1])
        assert 1.0 < co_dist < 1.3  # C=O bond is ~1.14 Å


class TestGeometryFreqOnly:

    def test_freq_only_raises(self):
        """Freq-only (--hess) output has no embedded geometry."""
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        with pytest.raises(LogError):
            log.load_geometry()


class TestEnergy:

    def test_freq_c3h7(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        energy = log.load_energy()
        expected = -9.907945789911 * constants.E_h * constants.Na
        assert abs(energy - expected) < abs(expected) * 1e-6

    def test_opt_c3h7(self):
        log = XTBLog(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False)
        energy = log.load_energy()
        expected = -9.907945789911 * constants.E_h * constants.Na
        assert abs(energy - expected) < abs(expected) * 1e-6

    def test_co2(self):
        """CO2: total energy -10.308452243026 Eh."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        energy = log.load_energy()
        expected = -10.308452243026 * constants.E_h * constants.Na
        assert abs(energy - expected) < abs(expected) * 1e-6

    def test_sp_ncc(self):
        """Single point NCC: -10.752192838993 Eh."""
        log = XTBLog(os.path.join(DATA_PATH, 'NCC_xTB.out'), check_for_errors=False)
        energy = log.load_energy()
        expected = -10.752192838993 * constants.E_h * constants.Na
        assert abs(energy - expected) < abs(expected) * 1e-6

    def test_ts(self):
        log = XTBLog(os.path.join(DATA_PATH, 'TS_NH2+N2H3_xtb.out'), check_for_errors=False)
        energy = log.load_energy()
        assert energy < 0  # should be negative

    def test_units_j_per_mol(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        energy = log.load_energy()
        assert energy < 0
        assert abs(energy) > 1e6  # millions of J/mol


class TestZPE:

    def test_freq_c3h7(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        zpe = log.load_zero_point_energy()
        expected = 0.086665171024 * constants.E_h * constants.Na
        assert abs(zpe - expected) < abs(expected) * 1e-6
        assert zpe > 0

    def test_ts_zpe(self):
        """TS_NH2+N2H3: ZPE = 0.056690417480 Eh."""
        log = XTBLog(os.path.join(DATA_PATH, 'TS_NH2+N2H3_xtb.out'), check_for_errors=False)
        zpe = log.load_zero_point_energy()
        expected = 0.056690417480 * constants.E_h * constants.Na
        assert abs(zpe - expected) < abs(expected) * 1e-6

    def test_opt_only_raises(self):
        log = XTBLog(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False)
        with pytest.raises(LogError):
            log.load_zero_point_energy()

    def test_sp_raises(self):
        log = XTBLog(os.path.join(DATA_PATH, 'NCC_xTB.out'), check_for_errors=False)
        with pytest.raises(LogError):
            log.load_zero_point_energy()


class TestFrequencies:

    def test_c3h7_count(self):
        """C3H7 (10 atoms): 3*10 - 6 = 24 real vibrational modes."""
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        freqs = log._load_frequencies()
        assert len(freqs) == 24

    def test_c3h7_values(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        freqs = log._load_frequencies()
        assert abs(freqs[0] - 63.60) < 0.1
        assert abs(freqs[-1] - 3104.88) < 0.1
        assert all(f > 0 for f in freqs)

    def test_co2_frequencies(self):
        """CO2 (3 atoms, linear): 3*3 - 5 = 4 real modes."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        freqs = log._load_frequencies()
        np.testing.assert_allclose(freqs, [600.70, 600.70, 1424.29, 2592.18], atol=0.01)

    def test_ts_frequencies(self):
        """TS should have real frequencies (imaginary excluded by _load_frequencies)."""
        log = XTBLog(os.path.join(DATA_PATH, 'TS_NH2+N2H3_xtb.out'), check_for_errors=False)
        freqs = log._load_frequencies()
        assert all(f > 0 for f in freqs)
        # TS has 8 atoms, 3*8 - 6 - 1(imaginary) = 17 real modes
        assert len(freqs) == 17

    def test_opt_only_no_freqs(self):
        """Optimization-only output without --hess should have no frequencies."""
        log = XTBLog(os.path.join(DATA_PATH, 'NCC_xTB.out'), check_for_errors=False)
        freqs = log._load_frequencies()
        assert freqs == []


class TestNegativeFrequency:

    def test_ts_imaginary(self):
        """TS_NH2+N2H3 should have one imaginary frequency at ~-781.89 cm^-1."""
        log = XTBLog(os.path.join(DATA_PATH, 'TS_NH2+N2H3_xtb.out'), check_for_errors=False)
        neg_freq = log.load_negative_frequency()
        assert abs(neg_freq - (-781.89)) < 0.1

    def test_stable_molecule_raises(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        with pytest.raises(LogError):
            log.load_negative_frequency()


class TestSpinMultiplicity:

    def test_c3h7_doublet(self):
        """C3H7 radical: spin=0.5, multiplicity=2."""
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        assert log._load_spin_multiplicity() == 2

    def test_co2_singlet(self):
        """CO2 closed-shell: should default to singlet."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        mult = log._load_spin_multiplicity()
        assert mult == 1


class TestConformer:

    def test_opt_c3h7(self):
        """Conformer from optimization output (has geometry, no freqs)."""
        log = XTBLog(os.path.join(DATA_PATH, 'output.out'), check_for_errors=False)
        conformer, unscaled_freqs = log.load_conformer(
            symmetry=1, spin_multiplicity=2, optical_isomers=1, label='C3H7')
        assert conformer.spin_multiplicity == 2
        assert conformer.optical_isomers == 1
        mode_types = [type(m).__name__ for m in conformer.modes]
        assert 'IdealGasTranslation' in mode_types
        assert 'NonlinearRotor' in mode_types

    def test_co2_linear_rotor(self):
        """CO2 is linear — conformer should have LinearRotor."""
        log = XTBLog(os.path.join(DATA_PATH, 'CO2_xtb.out'), check_for_errors=False)
        conformer, freqs = log.load_conformer(
            symmetry=2, spin_multiplicity=1, optical_isomers=1, label='CO2')
        mode_types = [type(m).__name__ for m in conformer.modes]
        assert 'LinearRotor' in mode_types
        assert 'HarmonicOscillator' in mode_types
        assert len(freqs) == 4  # CO2: 4 real modes

    def test_opt1_turbomol(self):
        """Conformer from Turbomol $coord geometry."""
        log = XTBLog(os.path.join(DATA_PATH, 'xtb_opt_1.out'), check_for_errors=False)
        conformer, freqs = log.load_conformer(
            symmetry=1, spin_multiplicity=1, optical_isomers=1, label='test')
        assert conformer is not None
        mode_types = [type(m).__name__ for m in conformer.modes]
        assert 'IdealGasTranslation' in mode_types

    def test_freq_only_raises(self):
        """Freq-only output can't build conformer (no geometry)."""
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        with pytest.raises(LogError):
            log.load_conformer(symmetry=1, spin_multiplicity=2, optical_isomers=1)


class TestUnsupported:

    def test_force_constant_matrix(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        assert log.load_force_constant_matrix() is None

    def test_scan_energies(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        with pytest.raises(NotImplementedError):
            log.load_scan_energies()

    def test_scan_pivot_atoms(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        with pytest.raises(NotImplementedError):
            log.load_scan_pivot_atoms()

    def test_scan_frozen_atoms(self):
        log = XTBLog(os.path.join(DATA_PATH, 'freq.out'), check_for_errors=False)
        with pytest.raises(NotImplementedError):
            log.load_scan_frozen_atoms()
