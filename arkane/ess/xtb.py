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
Arkane ESS adapter for xTB (extended Tight Binding) output files.
Supports geometry (V2000 and $coord/Turbomol formats), energy, frequencies,
ZPE, and conformer loading.
"""

import logging
import re

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import (
    Conformer,
    HarmonicOscillator,
    IdealGasTranslation,
    LinearRotor,
    NonlinearRotor,
)

from arkane.common import (
    get_element_mass,
    get_principal_moments_of_inertia,
    symbol_by_number,
)
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

logger = logging.getLogger(__name__)

BOHR_TO_ANGSTROM = constants.a0 * 1e10  # ~0.52918 Angstrom
INERTIA_ZERO_TOL = 1e-4  # amu*angstrom^2


def _normalize_symbol(raw: str) -> str:
    """Normalize an element symbol from xTB output (e.g. 'c' -> 'C', 'cl' -> 'Cl')."""
    if len(raw) == 1:
        return raw.upper()
    return raw[0].upper() + raw[1:].lower()


class XTBLog(ESSAdapter):
    """
    Arkane ESS adapter for xTB output files.

    Parses output from the xtb program (https://github.com/grimme-lab/xtb),
    including geometry, electronic energy, vibrational frequencies, and
    thermochemical data.

    Supported geometry formats:
      - V2000 (SDF/Molfile) block from ``--opt`` outputs
      - Turbomol ``$coord`` block (Bohr) from ``--opt`` outputs
    Frequency-only (``--hess``) outputs do not embed a geometry.
    """

    def check_for_errors(self):
        """Check for common xTB errors in the output file."""
        with open(self.path, 'r') as f:
            for line in f:
                if 'ERROR' in line and 'SETUP' not in line:
                    raise LogError(f'xTB error found in {self.path}: {line.strip()}')
                if 'abnormal termination' in line.lower():
                    raise LogError(f'xTB job terminated abnormally in {self.path}')

    def get_number_of_atoms(self) -> int:
        """
        Return the number of atoms from the xTB output.

        Tries the explicit 'number of atoms' line first, then falls
        back to counting atoms in the geometry block.
        """
        with open(self.path, 'r') as f:
            for line in f:
                if 'number of atoms' in line:
                    return int(line.split()[-1])
        try:
            coord, _, _ = self.load_geometry()
            return len(coord)
        except LogError:
            pass
        raise LogError(f'Could not determine the number of atoms in {self.path}')

    def load_geometry(self):
        """
        Load the molecular geometry from the xTB output file.

        Supports two formats found in xTB ``--opt`` outputs:
          - **V2000 (SDF/Molfile)**: Cartesian coordinates in Angstroms.
          - **Turbomol ``$coord``**: Cartesian coordinates in Bohr (converted to Angstroms).

        Frequency-only (``--hess``) outputs do not embed a geometry and will raise.

        Returns:
            Tuple of (coord, number, mass):
                coord: np.ndarray (n, 3) in Angstroms
                number: np.ndarray (n,) atomic numbers
                mass: np.ndarray (n,) atomic masses in amu
        """
        with open(self.path, 'r') as f:
            lines = f.readlines()
        coord, number, mass = self._parse_v2000(lines)
        if not coord:
            coord, number, mass = self._parse_turbomol_coord(lines)
        if not coord:
            raise LogError(f'Could not find geometry in {self.path}')
        return (np.array(coord, dtype=np.float64),
                np.array(number, dtype=int),
                np.array(mass, dtype=np.float64))

    @staticmethod
    def _parse_v2000(lines):
        """Parse V2000 (SDF/Molfile) geometry block. Returns (coord, number, mass) or empty lists."""
        coord, number, mass = [], [], []
        for i, line in enumerate(lines):
            if 'V2000' in line:
                n_atoms = int(line.split()[0])
                coord, number, mass = [], [], []
                for j in range(i + 1, i + 1 + n_atoms):
                    tokens = lines[j].split()
                    x, y, z = float(tokens[0]), float(tokens[1]), float(tokens[2])
                    symbol = tokens[3]
                    coord.append([x, y, z])
                    mass_i, num_i = get_element_mass(symbol)
                    number.append(num_i)
                    mass.append(mass_i)
        return coord, number, mass

    @staticmethod
    def _parse_turbomol_coord(lines):
        """Parse Turbomol $coord block (Bohr). Converts to Angstroms."""
        coord, number, mass = [], [], []
        in_coord = False
        for line in lines:
            stripped = line.strip()
            if stripped == '$coord':
                in_coord = True
                coord, number, mass = [], [], []
                continue
            if in_coord:
                if stripped.startswith('$'):
                    break
                tokens = stripped.split()
                if len(tokens) >= 4:
                    try:
                        x = float(tokens[0]) * BOHR_TO_ANGSTROM
                        y = float(tokens[1]) * BOHR_TO_ANGSTROM
                        z = float(tokens[2]) * BOHR_TO_ANGSTROM
                        symbol = _normalize_symbol(tokens[3])
                        coord.append([x, y, z])
                        mass_i, num_i = get_element_mass(symbol)
                        number.append(num_i)
                        mass.append(mass_i)
                    except (ValueError, KeyError):
                        break
        return coord, number, mass

    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the electronic energy in J/mol from the xTB output.

        Returns the last ``total energy`` value found. The zero-point energy
        is NOT included.
        """
        e_elect = None
        with open(self.path, 'r') as f:
            for line in f:
                if ':: total energy' in line:
                    match = re.search(r'(-?\d+\.\d+)\s+Eh', line)
                    if match:
                        e_elect = float(match.group(1))
                elif 'TOTAL ENERGY' in line and 'Eh' in line:
                    match = re.search(r'(-?\d+\.\d+)\s+Eh', line)
                    if match:
                        e_elect = float(match.group(1))
        if e_elect is None:
            raise LogError(f'Unable to find energy in xTB output file {self.path}')
        return e_elect * constants.E_h * constants.Na

    def load_zero_point_energy(self):
        """
        Load the zero-point energy in J/mol from the xTB output.
        Only available in frequency (``--hess``) calculations.
        """
        zpe = None
        with open(self.path, 'r') as f:
            for line in f:
                if ':: zero point energy' in line or 'zero-point vibrational energy' in line.lower():
                    match = re.search(r'(\d+\.\d+)\s+Eh', line)
                    if match:
                        zpe = float(match.group(1))
        if zpe is None:
            raise LogError(f'Unable to find zero-point energy in xTB output file {self.path}')
        return zpe * constants.E_h * constants.Na

    def _load_all_frequencies(self):
        """
        Load ALL vibrational frequencies from the last eigval block in the xTB output,
        including near-zero (translational/rotational) and negative (imaginary) modes.

        xTB prints frequencies twice (after Hessian and in Frequency Printout).
        Returns the last complete set.
        """
        all_blocks, current_block = [], []
        with open(self.path, 'r') as f:
            for line in f:
                if line.strip().startswith('eigval :'):
                    values = line.split(':')[1].split()
                    current_block.extend(float(v) for v in values)
                elif current_block and not line.strip().startswith('eigval'):
                    all_blocks.append(current_block)
                    current_block = []
        if current_block:
            all_blocks.append(current_block)
        return all_blocks[-1] if all_blocks else []

    def _load_frequencies(self):
        """
        Load positive (real) vibrational frequencies in cm^-1.
        Filters out translational/rotational modes (near-zero) and imaginary (negative).
        """
        return [f for f in self._load_all_frequencies() if f > 0.1]

    def _load_spin_multiplicity(self):
        """
        Determine spin multiplicity from the xTB output.

        xTB reports ``spin`` as S (not 2S+1), so multiplicity = 2*S + 1.
        Some xTB versions use ``--uhf N`` where N = number of unpaired electrons.
        Default to singlet.
        """
        with open(self.path, 'r') as f:
            for line in f:
                # Format: "  spin :  0.5"
                if 'spin' in line:
                    parts = line.split()
                    if parts and parts[0] == 'spin':
                        s = float(parts[-1])
                        return int(2 * s + 1)
                # Also check for unpaired electrons in setup block
                if 'unpaired' in line and ':' in line:
                    parts = line.split(':')
                    try:
                        n_unpaired = int(parts[-1].strip())
                        return n_unpaired + 1
                    except ValueError:
                        pass
        return 1

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degree of freedom data from an xTB output file.

        Requires geometry to be available in the log file (optimization outputs).
        For frequency-only outputs, geometry must be loaded from a separate file.

        Returns:
            Tuple of (Conformer, unscaled_frequencies).
        """
        if optical_isomers is None or symmetry is None:
            _optical_isomers, _symmetry, _ = self.get_symmetry_properties()
            if optical_isomers is None:
                optical_isomers = _optical_isomers
            if symmetry is None:
                symmetry = _symmetry

        if spin_multiplicity == 0:
            spin_multiplicity = self._load_spin_multiplicity()

        unscaled_frequencies = self._load_frequencies()
        modes = []

        # Translation
        coord, number, mass = self.load_geometry()
        modes.append(IdealGasTranslation(mass=(sum(mass), "amu")))

        # Rotation — use tolerance for linear molecule detection
        symbols = [symbol_by_number[i] for i in number]
        inertia = list(get_principal_moments_of_inertia(coord, numbers=number, symbols=symbols)[0])
        if inertia and not all(abs(i) < INERTIA_ZERO_TOL for i in inertia):
            nonzero = [i for i in inertia if abs(i) >= INERTIA_ZERO_TOL]
            if len(nonzero) < len(inertia):
                # Linear molecule: one or more moments are ~zero
                modes.append(LinearRotor(inertia=(max(nonzero), "amu*angstrom^2"), symmetry=symmetry))
            else:
                modes.append(NonlinearRotor(inertia=(nonzero, "amu*angstrom^2"), symmetry=symmetry))

        # Vibration
        if unscaled_frequencies:
            modes.append(HarmonicOscillator(frequencies=(unscaled_frequencies, "cm^-1")))

        return Conformer(E0=(0.0, "kJ/mol"),
                         modes=modes,
                         spin_multiplicity=spin_multiplicity,
                         optical_isomers=optical_isomers), unscaled_frequencies

    def load_negative_frequency(self):
        """Load the first imaginary (negative) frequency in cm^-1 for a transition state."""
        neg_freqs = [f for f in self._load_all_frequencies() if f < -0.1]
        if not neg_freqs:
            raise LogError(f'No imaginary frequencies found in {self.path}')
        return neg_freqs[0]

    def load_force_constant_matrix(self):
        """xTB writes a separate hessian file; not parsed from the main output."""
        return None

    def load_scan_energies(self):
        raise NotImplementedError('Rotor scans are not supported by the xTB adapter.')

    def load_scan_pivot_atoms(self):
        raise NotImplementedError('Rotor scans are not supported by the xTB adapter.')

    def load_scan_frozen_atoms(self):
        raise NotImplementedError('Rotor scans are not supported by the xTB adapter.')


register_ess_adapter('XTBLog', XTBLog)
