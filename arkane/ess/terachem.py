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
Arkane TeraChem module
Used to parse TeraChem output files
"""

import logging
import math
import os.path

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import HarmonicOscillator, Conformer

from arkane.common import check_conformer_energy, get_element_mass, symbol_by_number
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class TeraChemLog(ESSAdapter):
    """
    Represent a log file from TeraChem. The attribute `path` refers to the
    location on disk of the TeraChem log file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy arrays.
    TeraChemLog is an adapter for the abstract class ESSAdapter.
    """

    def check_for_errors(self):
        """
        Checks for common errors in a TeraChem log file.
        If any are found, this method will raise an error and crash.
        """
        with open(os.path.join(self.path), "r") as f:
            lines = f.readlines()
            error = None
            for line in reversed(lines):
                # check for common error messages
                if "incorrect method" in line.lower():
                    error = "incorrect method"
                    break
                elif "error: " in line.lower():
                    # e.g.: "ERROR: Closed shell calculations can't have spin multiplicity 0."
                    error = "multiplicity"
                    break
            if error:
                raise LogError(f"There was an error ({error}) with TeraChem output file {self.path} " f"due to line:\n{line}")

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in the TeraChem output file.
        Accepted output files: TeraChem's log file, xyz format file, TeraChem's output.geometry file.
        """
        n_atoms = 0
        with open(self.path, "r") as f:
            file_extension = os.path.splitext(self.path)[1]
            if file_extension == ".xyz":
                n_atoms = int(f.readline())
            else:
                line = f.readline()
                while line and n_atoms == 0:
                    if "Total atoms:" in line:
                        n_atoms = int(line.split()[-1])
                    elif "****** QM coordinates ******" in line or "Type         X              Y              Z            Mass" in line:
                        line = f.readline()
                        while line != "\n":
                            n_atoms += 1
                            line = f.readline()
                    line = f.readline()
        return n_atoms

    def load_force_constant_matrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates) from the
        TeraChem log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file, ``None`` is returned.
        """
        force = None
        n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3

        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Read force constant matrix
                if "*** Hessian Matrix (Hartree/Bohr^2) ***" in line:
                    force = np.zeros((n_rows, n_rows), float)
                    for i in range(int(math.ceil(n_rows / 6.0))):
                        # Matrix element rows
                        for j in range(n_rows):
                            line = f.readline()
                            while len(line.split()) not in [4, 7]:
                                # This is a header row
                                line = f.readline()
                            data = line.split()
                            for k in range(len(data) - 1):
                                force[j, i * 6 + k] = float(data[k + 1])
                    # Convert from atomic units (Hartree/Bohr^2) to SI (J/m^2)
                    force *= 4.35974417e-18 / 5.291772108e-11**2
                line = f.readline()

        return force

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        TeraChem log file. If multiple such geometries are identified, only the
        last is returned.
        """
        coords, numbers, masses = list(), list(), list()

        with open(self.path) as f:
            lines = f.readlines()

        num_of_atoms = None  # used to verify the result
        if os.path.splitext(self.path)[1] == ".xyz":
            skip_line = False
            for line in lines:
                if not skip_line and line.rstrip():
                    if len(line.split()) == 1 and line[0].isdigit():
                        num_of_atoms = int(line.rstrip())
                        skip_line = True  # the next line is just a comment, skip it
                        continue
                    splits = line.split()
                    coords.append([float(c) for c in splits[1:]])
                    mass, num = get_element_mass(splits[0])
                    masses.append(mass)
                    numbers.append(num)
                if skip_line:
                    skip_line = False
                    coords, numbers, masses = list(), list(), list()
        else:
            for i, line in enumerate(lines):
                if "Type         X              Y              Z            Mass" in line:
                    # this is an output.geometry file
                    j = i + 1
                    while lines[j].strip():
                        # example: '   C   0.6640965100   0.0039526500   0.0710079300  12.0000000000'
                        # or: ' C      0.512276     -0.516064      0.779232'
                        splits = lines[j].split()
                        coords.append([float(c) for c in splits[1:-1]])
                        masses.append(float(splits[-1]))
                        numbers.append(list(symbol_by_number.keys())[list(symbol_by_number.values()).index(splits[0])])
                        j += 1
                    break
                if "*** Reference Geometry ***" in line:
                    # this is an output.out file, e.g., from a freq run
                    j = i + 2
                    while lines[j].strip():
                        # example: ' C      0.512276     -0.516064      0.779232'
                        splits = lines[j].split()
                        coords.append([float(c) for c in splits[1:]])
                        mass, num = get_element_mass(splits[0])
                        masses.append(mass)
                        numbers.append(num)
                        j += 1
                    break

        coords = np.array(coords, float)
        numbers = np.array(numbers, np.int)
        masses = np.array(masses, float)
        if (
            len(coords) == 0
            or len(numbers) == 0
            or len(masses) == 0
            or ((len(coords) != num_of_atoms or len(numbers) != num_of_atoms or len(masses) != num_of_atoms) and num_of_atoms is not None)
        ):
            raise LogError(
                f"Unable to read atoms from TeraChem geometry output file {self.path}. "
                f"If this is a TeraChem optimization log file, try using either the "
                f"frequencies calculation log file (important if torsion modes exist) or "
                f'the "output.geometry" or a ".xyz" file instead.'
            )

        return coords, numbers, masses

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=""):
        """
        Load the molecular degree of freedom data from an output file created as the result of a
        TeraChem "Freq" calculation. As TeraChem's guess of the external symmetry number might not always correct,
        you can use the `symmetry` parameter to substitute your own value;
        if not provided, the value in the TeraChem output file will be adopted.
        """
        modes, unscaled_freqs = list(), list()
        converged = False
        if optical_isomers is None:
            _optical_isomers = self.get_symmetry_properties()[0]
            if optical_isomers is None:
                optical_isomers = _optical_isomers

        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Read spin multiplicity if not explicitly given
                if "Spin multiplicity" in line and spin_multiplicity == 0 and len(line.split()) == 3:
                    spin_multiplicity = int(float(line.split()[-1]))
                    logging.debug(f"Conformer {label} is assigned a spin multiplicity of {spin_multiplicity}")
                # Read vibrational modes
                elif "Mode      Eigenvalue(AU)      Frequency(cm-1)" in line:
                    line = f.readline()
                    while line != "\n":
                        # example:
                        # 'Mode  Eigenvalue(AU)  Frequency(cm-1)  Intensity(km/mol)   Vib.Temp(K)      ZPE(AU) ...'
                        # '  1     0.0331810528   170.5666870932      52.2294230772  245.3982965841   0.0003885795 ...'
                        if "i" not in line.split()[2]:
                            # only consider non-imaginary frequencies in this function
                            unscaled_freqs.append(float(line.split()[2]))
                        line = f.readline()
                if "Vibrational Frequencies/Thermochemical Analysis" in line:
                    converged = True
                line = f.readline()
            if not len(unscaled_freqs):
                raise LogError(f"Could not read frequencies from TeraChem log file {self.path}")
            if not converged:
                raise LogError(f"TeraChem job {self.path} did not converge.")
            modes.append(HarmonicOscillator(frequencies=(unscaled_freqs, "cm^-1")))

        return Conformer(E0=(0.0, "kJ/mol"), modes=modes, spin_multiplicity=spin_multiplicity, optical_isomers=optical_isomers), unscaled_freqs

    def load_energy(self, zpe_scale_factor=1.0):
        """
        Load the energy in J/mol from a TeraChem log file. Only the last energy
        in the file is returned, unless the log file represents a frequencies calculation,
        in which case the first energy is returned. The zero-point energy is *not* included
        in the returned value.
        """
        e_elect, return_first = None, False
        with open(self.path, "r") as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if "FREQUENCY ANALYSIS" in line:
                return_first = True
            if "Ground state energy (a.u.):" in line:
                e_elect = float(lines[i + 1].strip())
                if return_first:
                    break
            if "FINAL ENERGY:" in line:
                # example: 'FINAL ENERGY: -114.5008455547 a.u.'
                e_elect = float(line.split()[2])
                if return_first:
                    break
        if e_elect is None:
            raise LogError(f"Unable to find energy in TeraChem output file {self.path}.")
        return e_elect * constants.E_h * constants.Na

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a TeraChem log file.
        """
        zpe = None
        with open(self.path, "r") as f:
            for line in f:
                if "Vibrational zero-point energy (ZPE)" in line:
                    # example:
                    # 'Vibrational zero-point energy (ZPE) = 243113.467652369843563065 J/mol =     0.09259703 AU'
                    zpe = float(line.split("J/mol")[0].split()[-1])
                    logging.debug(f"ZPE is {zpe}")
        if zpe is not None:
            return zpe
        else:
            raise LogError(f"Unable to find zero-point energy in TeraChem output file {self.path}.")

    def load_scan_energies(self):
        """
        Extract the optimized energies in J/mol from a TeraChem torsional scan log file.
        """
        v_list = list()
        with open(self.path, "r") as f:
            lines = f.readlines()
            v_index, expected_num_of_points = 0, 0
            for line in lines:
                if "Scan Cycle" in line:
                    # example: '-=#=-     Scan Cycle     5/37        -=#=-'
                    v_index += 1
                    if not expected_num_of_points:
                        expected_num_of_points = int(line.split()[3].split("/")[1])
                if "Optimized Energy:" in line:
                    # example: '-=#=- Optimized Energy:    -155.0315243910 a.u.'
                    v = float(line.split()[3])
                    if len(v_list) == v_index - 1:
                        # append this point, it is in order
                        v_list.append(v)
                    elif len(v_list) < v_index - 1:
                        # seems like points in this scan are missing... add None's instead,
                        # later they'll be removed along with the corresponding angles
                        v_list.extend([None] * (v_index - 1 - len(v_list)))
                    else:
                        # we added more points that we should have, something is wrong with the log file or this method
                        raise LogError(f"Could not parse scan energies from {self.path}")
        logging.info("   Assuming {0} is the output from a TeraChem PES scan...".format(os.path.basename(self.path)))

        v_list = np.array(v_list, float)

        # check to see if the scanlog indicates that one of the reacting species may not be the lowest energy conformer
        check_conformer_energy(v_list, self.path)

        # Adjust energies to be relative to minimum energy conformer
        # Also convert units from Hartree/particle to J/mol
        v_list -= np.min(v_list)
        v_list *= constants.E_h * constants.Na
        angles = np.arange(0.0, 2 * math.pi + 0.00001, 2 * math.pi / (len(v_list) - 1), float)

        # remove None's:
        indices_to_pop = [v_list.index[entry] for entry in v_list if entry is None]
        for i in reversed(indices_to_pop):
            v_list.pop(i)
            angles.pop(i)
        if v_index != expected_num_of_points:
            raise LogError(f"Expected to find {expected_num_of_points} scan points in TeraChem scan log file " f"{self.path}, but found: {v_index}")

        return v_list, angles

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        frequencies = []
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Read vibrational modes
                if "Mode      Eigenvalue(AU)      Frequency(cm-1)" in line:
                    line = f.readline()
                    # example:
                    # 'Mode  Eigenvalue(AU)  Frequency(cm-1)  Intensity(km/mol)   Vib.Temp(K)      ZPE(AU) ...'
                    # '  1     0.0331810528   170.5666870932i     52.2294230772  245.3982965841   0.0003885795 ...'
                    while "i" in line:
                        frequencies.append(-1 * float(line.split()[2][:-1]))  # remove 'i'
                        line = f.readline()
                    break
                f.readline()
        if len(frequencies) == 1:
            return frequencies[0]
        elif len(frequencies) > 1:
            logging.info("More than one imaginary frequency in TeraChem output file {0}.".format(self.path))
            return frequencies[0]
        else:
            raise LogError(f"Unable to find imaginary frequency in TeraChem output file {self.path}.")

    def load_scan_pivot_atoms(self):
        """Not implemented for TeraChem"""
        raise NotImplementedError("The load_scan_pivot_atoms method is not implemented for TeraChem Logs")

    def load_scan_frozen_atoms(self):
        """Not implemented for TeraChem"""
        raise NotImplementedError("The load_scan_frozen_atoms method is not implemented for TeraChem Logs")


register_ess_adapter("TeraChemLog", TeraChemLog)
