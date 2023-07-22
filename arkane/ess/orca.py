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
Arkane Orca module
Used to parse Orca output files
"""
import os.path
import math
import numpy as np
import rmgpy.constants as constants
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.common import get_element_mass, get_principal_moments_of_inertia, symbol_by_number
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class OrcaLog(ESSAdapter):
    """
    Represent an output file from Orca. The attribute `path` refers to the
    location on disk of the Orca output file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays. OrcaLog is an adapter for the abstract class ESSAdapter.
    """

    def check_for_errors(self):
        """
        Checks for common errors in an Orca log file.
        If any are found, this method will raise an error and crash.
        """
        with open(os.path.join(self.path), "r") as f:
            lines = f.readlines()
            error = None
            for line in reversed(lines):
                # check for common error messages
                if "ORCA finished by error termination in SCF" in line:
                    error = "SCF"
                    break
                elif "ORCA finished by error termination in MDCI" in line:
                    error = "MDCI"
                    break
                elif "Error : multiplicity" in line:
                    error = f"The multiplicity and charge combination for species {species_label} are wrong."
                    break
                elif "ORCA TERMINATED NORMALLY" in line:
                    break
            if error:
                raise LogError(f"There was an error ({error}) with Orca output file {self.path} " f"due to line:\n{line}")

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the Orca output file.
        """
        natoms = 0

        with open(self.path, "r") as f:
            line = f.readline()
            while line != "" and natoms == 0:
                # Automatically determine the number of atoms
                if "CARTESIAN COORDINATES (ANGSTROEM)" in line and natoms == 0:
                    for i in range(2):
                        line = f.readline()

                    while "---------------------------------" not in line:
                        natoms += 1
                        line = f.readline()
                        if not line.strip():
                            f.close()
                            return natoms
                line = f.readline()

    def load_force_constant_matrix(self):
        """
        Return the force constant matrix from the Orca log file.
        Orca prints the Hessian matrix as a '.hess' file which must be located
        in the same folder as the log file for this method to parse it.
        The units of the returned force constants are J/m^2.
        If no force constant matrix can be found, ``None`` is returned.
        """
        hess_files = list()
        for _, _, files in os.walk(os.path.dirname(self.path)):
            for file_ in files:
                if file_.endswith(".hess"):
                    hess_files.append(file_)
            break
        if len(hess_files) == 1:
            hess_file = hess_files[0]
        else:
            expected_hess_name = f"{os.path.basename(self.path).split('.')[0]}.hess"
            for hess_file in hess_files:
                if hess_file == expected_hess_name:
                    break
            else:
                raise LogError(
                    f"Could not identify the intended .hess file in {os.path.dirname(self.path)}.\n"
                    f"Expected to find {expected_hess_name} in that folder."
                )

        force = None
        n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3
        with open(hess_file, "r") as f:
            line = f.readline()
            while line != "":
                if "$hessian" in line:
                    line = f.readline()
                    force = np.zeros((n_rows, n_rows), float)
                    for i in range(int(math.ceil(n_rows / 5.0))):
                        line = f.readline()
                        for j in range(n_rows):
                            data = f.readline().split()[1:]
                            for k in range(len(data)):
                                force[j, i * 5 + k] = float(data[k])
                    # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                    force *= 4.35974417e-18 / 5.291772108e-11**2
                line = f.readline()
        return force

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Orca log file. If multiple such geometries are identified, only the
        last is returned.
        """
        atoms, coords, numbers, mass = [], [], [], []

        with open(self.path) as f:
            log = f.readlines()

        # Now look for the geometry.
        # Will return the final geometry in the file under Standard Nuclear Orientation.
        geometry_flag = False
        for i in reversed(range(len(log))):
            line = log[i]
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                for line in log[(i + 2) :]:
                    if not line.strip():
                        break
                    if "---------------------------------" not in line:
                        data = line.split()
                        atoms.append(data[0])
                        coords.append([float(c) for c in data[1:]])
                        geometry_flag = True

                if geometry_flag:
                    break

        # Assign appropriate mass to each atom in the molecule
        for atom1 in atoms:
            mass1, num1 = get_element_mass(atom1)
            mass.append(mass1)
            numbers.append(num1)
        coord = np.array(coords, float)
        number = np.array(numbers, int)
        mass = np.array(mass, float)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise LogError(f"Unable to read atoms from orca geometry output file {self.path}.")

        return coord, number, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=""):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of an Orca "Freq" quantum chemistry calculation. As
        Orca's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value; if
        not provided, the value in the Orca log file will be adopted.
        """
        freq, mmass, rot, unscaled_frequencies = [], [], [], []
        e0 = 0.0

        if optical_isomers is None or symmetry is None:
            _optical_isomers, _symmetry, _ = self.get_symmetry_properties()
            if optical_isomers is None:
                optical_isomers = _optical_isomers
            if symmetry is None:
                symmetry = _symmetry

        with open(self.path) as f:
            log = f.readlines()
        for i, line in enumerate(log):
            if spin_multiplicity == 0 and " Multiplicity           Mult" in line:
                spin_multiplicity = int(float(line.split()[3]))

            if " Mode    freq" in line:
                frequencies = list()
                for line_ in log[(i + 2) :]:
                    if not line_.strip():
                        break
                    frequencies.extend([float(line_.split()[1])])
                else:
                    raise Exception(f"Frequencies not found in {self.path}")
                frequencies = [f for f in frequencies if f > 0.0]
                unscaled_frequencies = frequencies
                vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                freq.append(vibration)

        # Get moments of inertia from external rotational modes, given in atomic units.
        coord, number, mass = self.load_geometry()
        symbols = [symbol_by_number[i] for i in number]
        inertia = get_principal_moments_of_inertia(coord, numbers=number, symbols=symbols)
        inertia = list(inertia[0])
        if len(inertia):
            if any(i == 0.0 for i in inertia):
                inertia.remove(0.0)
                rot.append(LinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry))
            else:
                rot.append(NonlinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry))

        translation = IdealGasTranslation(mass=(sum(mass), "amu"))
        mmass.append(translation)

        # Take only the last modes found (in the event of multiple jobs).
        modes = [mmass[-1], rot[-1], freq[-1]]

        return (
            Conformer(E0=(e0 * 0.001, "kJ/mol"), modes=modes, spin_multiplicity=spin_multiplicity, optical_isomers=optical_isomers),
            unscaled_frequencies,
        )

    def load_energy(self, zpe_scale_factor=1.0):
        """
        Load the energy in J/ml from an Orca log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        e_elect = None
        with open(self.path, "r") as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:  # for all methods in Orca
                    e_elect = float(line.split()[-1])
        if e_elect is None:
            raise LogError("Unable to find energy in Orca output file.")
        return e_elect * constants.E_h * constants.Na

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a Orca output file.
        """
        zpe = None
        with open(self.path, "r") as f:
            for line in f:
                if "Zero point energy" in line:
                    zpe = float(line.split()[-4]) * constants.E_h * constants.Na
        if zpe is None:
            raise LogError("Unable to find zero-point energy in Orca output file.")
        return zpe

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency calculation in cm^-1.
        Since there can be many imaginary frequencies, only the first one is returned.
        """
        frequencies = list()
        with open(self.path) as f:
            log = f.readlines()
        for line in log:
            if "***imaginary mode***" in line:
                frequencies.append(float(line.split()[1]))

        if len(frequencies) == 1:
            return frequencies[0]
        elif len(frequencies) > 1:
            logging.info("More than one imaginary frequency in Orca output file {0}.".format(self.path))
            return frequencies[0]
        else:
            raise LogError(f"Unable to find imaginary frequency in Orca output file {self.path}")

    def get_T1_diagnostic(self):
        """
        Returns the T1 diagnostic from output log.
        If multiple occurrences exist, returns the last occurrence.
        T1 diagnostic only available in Coupled Cluster calculations.
        """
        with open(self.path) as f:
            log = f.readlines()
        for line in reversed(log):
            if "T1 diagnostic " in line:
                items = line.split()
                return float(items[-1])

    def load_scan_energies(self):
        """not implemented in Orca"""
        raise NotImplementedError("The load_scan_energies method is not implemented for Orca.")

    def load_scan_pivot_atoms(self):
        """Not implemented for Orca"""
        raise NotImplementedError("The load_scan_pivot_atoms method is not implemented for Orca.")

    def load_scan_frozen_atoms(self):
        """Not implemented for Orca"""
        raise NotImplementedError("The load_scan_frozen_atoms method is not implemented for Orca.")


register_ess_adapter("OrcaLog", OrcaLog)
