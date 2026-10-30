#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'), #
# to deal in the Software without restriction, including without limitation  #
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
Tests for the Arkane adapter of ARC's YAML outputs (arkane/ess/arc_yaml.py).
The fixtures are outputs of ARC's ASE adapter with the UMA-S 1.2 potential.
"""

import math
import os
import shutil
import tempfile

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.statmech import HarmonicOscillator, NonlinearRotor

from arkane.exceptions import LogError
from arkane.ess.arc_yaml import ARCYAMLLog
from arkane.ess.factory import ess_factory

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '..', '..', 'arkane', 'data', 'arc_yaml')


def wavenumbers(force_si, mass_amu):
    """Positive wavenumbers (cm^-1) of a force constant matrix in J/m^2."""
    weights = np.repeat(np.sqrt(np.asarray(mass_amu) * constants.amu), 3)
    eig = np.linalg.eigvalsh(force_si / np.outer(weights, weights))
    return np.sort(np.sqrt(eig[eig > 0]) / (2 * math.pi * constants.c * 100))


class TestARCYAMLLog:
    def test_factory_recognizes_yaml_by_extension_and_by_content(self):
        assert isinstance(ess_factory(os.path.join(DATA, 'H2O_freq_uma.yml')), ARCYAMLLog)
        with tempfile.TemporaryDirectory() as tmp:
            freq_out = os.path.join(tmp, 'freq.out')  # ARC copies its outputs to this name
            shutil.copyfile(os.path.join(DATA, 'H2O_freq_uma.yml'), freq_out)
            assert isinstance(ess_factory(freq_out), ARCYAMLLog)

    def test_load_geometry(self):
        coord, number, mass = ARCYAMLLog(os.path.join(DATA, 'H2O_freq_uma.yml')).load_geometry()
        assert coord.shape == (3, 3)
        assert sorted(number.tolist()) == [1, 1, 8]
        assert sum(mass) == pytest.approx(18.0106, abs=1e-3)

    def test_load_conformer(self):
        log = ARCYAMLLog(os.path.join(DATA, 'H2O_freq_uma.yml'))
        conformer, frequencies = log.load_conformer(symmetry=2, spin_multiplicity=1, optical_isomers=1)
        assert len(frequencies) == 3
        assert any(isinstance(mode, NonlinearRotor) for mode in conformer.modes)
        oscillator = [mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)][0]
        assert oscillator.frequencies.value_si.tolist() == pytest.approx(frequencies)
        with pytest.raises(LogError):
            log.load_conformer(symmetry=2, optical_isomers=1)  # no multiplicity given, and no job input next to it

    def test_spin_multiplicity_from_the_job_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copyfile(os.path.join(DATA, 'CH3_freq_uma.yml'), os.path.join(tmp, 'output.yml'))
            with open(os.path.join(tmp, 'input.yml'), 'w') as f:
                f.write('job_type: freq\ncharge: 0\nmultiplicity: 2\n')
            conformer, _ = ARCYAMLLog(os.path.join(tmp, 'output.yml')).load_conformer(symmetry=6, optical_isomers=1)
            assert conformer.spin_multiplicity == 2

    def test_force_constant_matrix_reproduces_the_frequencies(self):
        for name in ('H2O_freq_uma.yml', 'CH3_freq_uma.yml', 'H2O_freq_uma_ev_hessian.out'):
            log = ARCYAMLLog(os.path.join(DATA, name))
            force = log.load_force_constant_matrix()
            _, _, mass = log.load_geometry()
            stored = sorted(f for f in log.data['freqs'] if f > 0)
            # the unprojected Hessian has near-zero rotations too; the highest modes are the vibrations
            assert wavenumbers(force, mass)[-len(stored):].tolist() == pytest.approx(stored, abs=2.0), name

    def test_hessian_units_are_detected(self):
        assert ARCYAMLLog(os.path.join(DATA, 'H2O_freq_uma.yml'))._hessian_units(
            np.array(ARCYAMLLog(os.path.join(DATA, 'H2O_freq_uma.yml')).data['hessian'])) == 'hartree/angstrom^2'
        old = ARCYAMLLog(os.path.join(DATA, 'H2O_freq_uma_ev_hessian.out'))
        assert old._hessian_units(np.array(old.data['hessian'])) == 'ev/angstrom^2'

    def test_transition_state(self):
        log = ARCYAMLLog(os.path.join(DATA, 'TS_OHCH4_freq_uma.yml'))
        assert log.load_negative_frequency() == pytest.approx(-1058.2, abs=0.1)
        _, frequencies = log.load_conformer(symmetry=1, spin_multiplicity=2, optical_isomers=1)
        assert all(f > 0 for f in frequencies)
        assert len(frequencies) == 3 * log.get_number_of_atoms() - 7  # 3N-6 minus the reaction coordinate

    def test_zero_point_energy(self):
        log = ARCYAMLLog(os.path.join(DATA, 'H2O_freq_uma.yml'))
        expected = 0.5 * constants.h * constants.c * 100 * sum(log.data['freqs']) * constants.Na
        assert log.load_zero_point_energy() == pytest.approx(expected)
        assert log.load_zero_point_energy() / 1000 == pytest.approx(56.0, abs=0.5)  # kJ/mol

    def test_load_scan_energies(self):
        vlist, angle = ARCYAMLLog(os.path.join(DATA, 'C3H8_scan_uma.yml')).load_scan_energies()
        assert len(vlist) == len(angle) == 46
        assert min(vlist) == 0.0
        assert angle[0] == 0.0 and angle[-1] == pytest.approx(2 * math.pi)
        assert np.all(np.diff(angle) > 0)
        assert 5e3 < max(vlist) < 2e4  # a methyl rotor barrier, J/mol

    def test_a_yaml_scan_is_a_hindered_rotor_scan_log(self):
        """StatMechJob fits a hindered rotor to an ARC YAML scan, as it does to a Gaussian one"""
        from rmgpy.species import Species
        from rmgpy.statmech import HinderedRotor
        from arkane.statmech import StatMechJob
        freq, scan = os.path.join(DATA, 'C3H8_freq_uma.yml'), os.path.join(DATA, 'C3H8_scan_uma.yml')
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'C3H8.py')
            with open(path, 'w') as f:
                f.write(f"spin_multiplicity = 1\noptical_isomers = 1\nexternal_symmetry = 2\nenergy = -119.0\n"
                        f"geometry = Log('{freq}')\nfrequencies = Log('{freq}')\n"
                        f"rotors = [HinderedRotor(scanLog=Log('{scan}'), pivots=[1, 2], top=[1, 4, 5, 6], "
                        f"symmetry=3, fit='fourier')]\n")
            species = Species(label='C3H8', smiles='CCC')
            job = StatMechJob(species=species, path=path)
            job.applyAtomEnergyCorrections = job.applyBondEnergyCorrections = False
            job.load(pdep=False, plot=False)
        rotors = [mode for mode in species.conformer.modes if isinstance(mode, HinderedRotor)]
        assert len(rotors) == 1
        assert rotors[0].symmetry == 3
        # The methyl torsion leaves the harmonic modes, which keep 3N - 6 - 1 = 26 of propane's 27 frequencies.
        oscillator = [mode for mode in species.conformer.modes if isinstance(mode, HarmonicOscillator)][0]
        assert len(oscillator.frequencies.value_si) == 26

    def test_a_failed_job_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'output.yml')
            with open(path, 'w') as f:
                f.write("error: 'Optimization failed: SCF did not converge'\nopt_xyz: {symbols: [H], coords: [[0, 0, 0]]}\n")
            with pytest.raises(LogError):
                ARCYAMLLog(path)
