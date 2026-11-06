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
Arkane ESS adapter for the YAML output files that ARC writes for programs without a log of their own,
e.g., ARC's ASE adapter (machine-learned interatomic potentials such as UMA, TorchANI, or xTB through ASE).

The file is a YAML mapping. The keys read here are:

- ``opt_xyz`` or ``xyz``: the geometry, ``{'symbols': [...], 'coords': [[x, y, z], ...], 'isotopes': [...]}`` in Angstrom
- ``freqs``: the vibrational frequencies in cm^-1, imaginary ones as negative numbers
- ``hessian``: the Cartesian Hessian, (3N x 3N) or (N, N, 3, 3). ARC writes it in Hartree/Angstrom^2, but older
  ARC versions wrote eV/Angstrom^2 for some calculators, so the units are taken from ``hessian_units`` if given,
  and are otherwise the ones that reproduce ``freqs``.
- ``sp`` (kJ/mol) or ``energy`` (Hartree): the electronic energy; ``zpe`` (kJ/mol), optional
- ``energies`` (Hartree) and ``angles`` (degrees): a 1D rotor scan
- ``error``: set by ARC when the job failed
"""

import logging
import math
import os

import numpy as np
import yaml

import rmgpy.constants as constants
from rmgpy.statmech import Conformer, HarmonicOscillator, IdealGasTranslation, LinearRotor, NonlinearRotor

from arkane.common import check_conformer_energy, get_element_mass, get_principal_moments_of_inertia
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

logger = logging.getLogger(__name__)

# Keys that identify a YAML output written by ARC.
ARC_YAML_KEYS = ('freqs', 'frequencies', 'hessian', 'opt_xyz', 'xyz', 'sp', 'energies')

# Hessian units ARC has written, as factors to Hartree/Bohr^2.
HESSIAN_UNITS_TO_HARTREE_PER_BOHR2 = {
    'hartree/angstrom^2': (constants.a0 * 1e10) ** 2,
    'ev/angstrom^2': (constants.a0 * 1e10) ** 2 / (constants.E_h / constants.e),
    'hartree/bohr^2': 1.0,
}
INERTIA_ZERO_TOL = 1e-4  # amu*angstrom^2


def is_arc_yaml(path: str) -> bool:
    """
    Whether a file is a YAML output written by ARC, as recognized by its content rather than its extension
    (ARC copies these files to names such as ``freq.out``).

    Args:
        path (str): The path to the file.

    Returns:
        bool: Whether the file holds a YAML mapping with any of ARC's output keys.
    """
    try:
        content = load_arc_yaml(path)
    except Exception:
        return False
    return isinstance(content, dict) and any(key in content for key in ARC_YAML_KEYS)


def load_arc_yaml(path: str) -> dict:
    """
    Load an ARC YAML output. ARC writes coordinates as tuples, which the safe loader cannot construct.
    """
    with open(path, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


class ARCYAMLLog(ESSAdapter):
    """
    A class for parsing the YAML output files written by ARC (e.g., by its ASE adapter).
    """

    def __init__(self, path, check_for_errors=True, scratch_directory=None):
        self.data = load_arc_yaml(path)
        if not isinstance(self.data, dict):
            raise LogError(f'The file {path} is not an ARC YAML output.')
        super().__init__(path, check_for_errors=check_for_errors, scratch_directory=scratch_directory)

    def check_for_errors(self):
        """
        Raise an error if ARC recorded that the job failed.
        """
        if self.data.get('error'):
            raise LogError(f'The job of {self.path} failed: {self.data["error"]}')

    def _xyz(self) -> dict:
        xyz = self.data.get('opt_xyz') or self.data.get('xyz')
        if not xyz or not xyz.get('symbols'):
            raise LogError(f'No geometry found in {self.path}.')
        return xyz

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecule.
        """
        return len(self._xyz()['symbols'])

    def load_geometry(self):
        """
        Return the geometry as Cartesian coordinates (Angstrom), atomic numbers, and atomic masses (amu).
        """
        xyz = self._xyz()
        isotopes = xyz.get('isotopes') or [None] * len(xyz['symbols'])
        number, mass = list(), list()
        for symbol, isotope in zip(xyz['symbols'], isotopes):
            m, z = get_element_mass(symbol, isotope)
            number.append(z)
            mass.append(m)
        coord = np.array(xyz['coords'], dtype=np.float64).reshape(-1, 3)
        return coord, np.array(number, int), np.array(mass, np.float64)

    def _all_frequencies(self) -> list:
        freqs = self.data.get('freqs')
        if freqs is None:
            freqs = self.data.get('frequencies')
        return [float(f) for f in freqs] if freqs is not None else list()

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degrees of freedom.

        ARC's YAML output has no spin multiplicity; if none is given, it is taken from ARC's job input.
        The frequencies are the real ones; translations and rotations are not among them.

        Returns:
            tuple: The Conformer (E0 not set) and the unscaled real frequencies in cm^-1.
        """
        if optical_isomers is None or symmetry is None:
            _optical_isomers, _symmetry, _ = self.get_symmetry_properties()
            optical_isomers = _optical_isomers if optical_isomers is None else optical_isomers
            symmetry = _symmetry if symmetry is None else symmetry
        if not spin_multiplicity:
            spin_multiplicity = self._load_spin_multiplicity()
        coord, number, mass = self.load_geometry()
        modes = [IdealGasTranslation(mass=(float(sum(mass)), 'amu'))]
        if len(number) > 1:
            inertia = list(get_principal_moments_of_inertia(coord, numbers=number)[0])
            nonzero = [i for i in inertia if abs(i) >= INERTIA_ZERO_TOL]
            if len(nonzero) < len(inertia):
                modes.append(LinearRotor(inertia=(max(nonzero), 'amu*angstrom^2'), symmetry=symmetry))
            else:
                modes.append(NonlinearRotor(inertia=(nonzero, 'amu*angstrom^2'), symmetry=symmetry))
        unscaled_frequencies = [f for f in self._all_frequencies() if f > 0]
        if unscaled_frequencies:
            modes.append(HarmonicOscillator(frequencies=(unscaled_frequencies, 'cm^-1')))
        return Conformer(E0=(0.0, 'kJ/mol'), modes=modes, spin_multiplicity=spin_multiplicity,
                         optical_isomers=optical_isomers), unscaled_frequencies

    def _load_spin_multiplicity(self) -> int:
        """
        The spin multiplicity: ``multiplicity`` in the output if present, otherwise in the ``input.yml`` that ARC
        writes next to the output in the job directory.
        """
        if self.data.get('multiplicity'):
            return int(self.data['multiplicity'])
        input_path = os.path.join(os.path.dirname(self.path), 'input.yml')
        if os.path.isfile(input_path):
            multiplicity = (load_arc_yaml(input_path) or dict()).get('multiplicity')
            if multiplicity:
                return int(multiplicity)
        raise LogError(f'No spin multiplicity found for {self.path}: set spin_multiplicity in the species file.')

    def _hessian_units(self, hessian: np.ndarray) -> str:
        """
        The units of the stored Hessian: ``hessian_units`` if given, otherwise the units under which the Hessian
        reproduces the stored frequencies (Hartree/Angstrom^2 if there are none).
        """
        units = self.data.get('hessian_units')
        if units:
            units = str(units).lower().replace(' ', '')
            if units not in HESSIAN_UNITS_TO_HARTREE_PER_BOHR2:
                raise LogError(f'Unknown hessian_units {units} in {self.path}.')
            return units
        freqs = [f for f in self._all_frequencies() if f > 0]
        if not freqs:
            return 'hartree/angstrom^2'
        _, _, mass = self.load_geometry()
        best, best_error = 'hartree/angstrom^2', None
        for units, factor in HESSIAN_UNITS_TO_HARTREE_PER_BOHR2.items():
            computed = self._frequencies_from_hessian(hessian * factor, mass)
            if len(computed) < len(freqs):
                continue
            error = abs(computed[-1] - max(freqs)) / max(freqs)
            if best_error is None or error < best_error:
                best, best_error = units, error
        return best

    @staticmethod
    def _frequencies_from_hessian(hessian_hartree_bohr2: np.ndarray, mass: np.ndarray) -> np.ndarray:
        """
        The sorted positive vibrational wavenumbers (cm^-1) of a Hessian in Hartree/Bohr^2, without projection.
        """
        to_si = constants.E_h / constants.a0 ** 2  # Hartree/Bohr^2 to J/m^2
        weights = np.repeat(np.sqrt(mass * constants.amu), 3)
        eigenvalues = np.linalg.eigvalsh(hessian_hartree_bohr2 * to_si / np.outer(weights, weights))
        return np.sort(np.sqrt(eigenvalues[eigenvalues > 0]) / (2 * math.pi * constants.c * 100))

    def load_force_constant_matrix(self):
        """
        Return the Cartesian force constant matrix in J/m^2, or ``None`` if the file has no Hessian.
        """
        hessian = self.data.get('hessian')
        if hessian is None:
            return None
        n_coords = 3 * self.get_number_of_atoms()
        hessian = np.array(hessian, dtype=np.float64).reshape(n_coords, n_coords)
        hessian = hessian * HESSIAN_UNITS_TO_HARTREE_PER_BOHR2[self._hessian_units(hessian)]
        return hessian * constants.E_h / constants.a0 ** 2

    def load_energy(self, zpe_scale_factor=1.):
        """
        Return the electronic energy in J/mol, without the zero-point energy.
        """
        if self.data.get('sp') is not None:
            return float(self.data['sp']) * 1000.0
        if self.data.get('energy') is not None:
            return float(self.data['energy']) * constants.E_h * constants.Na
        raise LogError(f'No electronic energy found in {self.path}.')

    def load_zero_point_energy(self):
        """
        Return the zero-point energy in J/mol: ``zpe`` if stored, otherwise the harmonic ZPE of the real frequencies.
        """
        if self.data.get('zpe') is not None:
            return float(self.data['zpe']) * 1000.0
        freqs = [f for f in self._all_frequencies() if f > 0]
        if not freqs:
            raise LogError(f'No zero-point energy or frequencies found in {self.path}.')
        return 0.5 * constants.h * constants.c * 100.0 * sum(freqs) * constants.Na

    def load_scan_energies(self):
        """
        Return the energies of a 1D rotor scan in J/mol, relative to the lowest one, and the scan angles in radians
        from 0 to 2*pi, as the other adapters do.
        """
        energies, angles = self.data.get('energies'), self.data.get('angles')
        if not energies or not angles or len(energies) != len(angles):
            raise LogError(f'No rotor scan found in {self.path}.')
        vlist = np.array(energies, dtype=np.float64)
        check_conformer_energy(vlist, self.path)
        vlist = (vlist - np.min(vlist)) * constants.E_h * constants.Na
        angle = np.radians(np.array(angles, dtype=np.float64))
        angle = angle - angle[0]
        angle[1:] = np.where(angle[1:] <= 0, angle[1:] + 2 * math.pi, angle[1:])
        return vlist, angle

    def load_negative_frequency(self):
        """
        Return the imaginary frequency of a transition state in cm^-1 (a negative number).
        """
        negative = [f for f in self._all_frequencies() if f < 0]
        if not negative:
            raise LogError(f'No imaginary frequency found in {self.path}.')
        return min(negative)

    def load_scan_pivot_atoms(self):
        """
        ARC's YAML scan output does not record the pivots; ARC gives them in the species file.
        """
        raise NotImplementedError('ARC YAML scan outputs do not record the scan pivots; give them in the species file.')

    def load_scan_frozen_atoms(self):
        """
        ARC's YAML scan output does not record frozen atoms.
        """
        return dict()


register_ess_adapter('ARCYAMLLog', ARCYAMLLog)
