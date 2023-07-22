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
Arkane Molpro module
Used to parse Molpro output files
"""

import logging
import math
import os.path

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.common import get_element_mass
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class MolproLog(ESSAdapter):
    """
    Represents a Molpro log file. The attribute `path` refers to the
    location on disk of the Molpro log file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy arrays.
    MolproLog is an adapter for the abstract class ESSAdapter.
    """

    def check_for_errors(self):
        """
        Checks for common errors in a Molpro log file.
        If any are found, this method will raise an error and crash.
        """
        with open(os.path.join(self.path), "r") as f:
            lines = f.readlines()
            error = None
            for line in reversed(lines):
                if "molpro calculation terminated" in line.lower() or "variable memory released" in line.lower():
                    break
                # check for common error messages
                elif "No convergence" in line:
                    error = "Unconverged"
                elif "A further" in line and "Mwords of memory are needed" in line and "Increase memory to" in line:
                    # e.g.: `A further 246.03 Mwords of memory are needed for the triples to run.
                    # Increase memory to 996.31 Mwords.` (w/o the line break)
                    error = "Memory"
                elif "insufficient memory available - require" in line:
                    # e.g.: `insufficient memory available - require              228765625  have
                    #        62928590
                    #        the request was for real words`
                    # add_mem = (float(line.split()[-2]) - float(prev_line.split()[0])) / 1e6
                    error = "Memory"
                elif "Basis library exhausted" in line:
                    # e.g.:
                    # ` SETTING BASIS          =    6-311G**
                    #
                    #
                    #  Using spherical harmonics
                    #
                    #  LIBRARY EXHAUSTED
                    #   Searching for I  S 6-311G
                    #   Library contains the following bases:
                    #  ? Error
                    #  ? Basis library exhausted
                    #  ? The problem occurs in Binput`
                    basis_set = None
                    for line0 in reversed(line):
                        if "SETTING BASIS" in line0:
                            basis_set = line0.split()[-1]
                    error = f"Unrecognized basis set {basis_set}"
                if error:
                    raise LogError(f"There was an error ({error}) with Molpro output file {self.path} " f"due to line:\n{line}")

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the MolPro log file.
        """
        n_atoms = 0
        # Open Molpro log file for parsing
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "" and n_atoms == 0:
                # Automatically determine the number of atoms
                if "ATOMIC COORDINATES" in line and n_atoms == 0:
                    for i in range(4):
                        line = f.readline()
                    while "Bond lengths" not in line and "nuclear charge" not in line.lower():
                        n_atoms += 1
                        line = f.readline()
                line = f.readline()

        return n_atoms - 1

    def load_force_constant_matrix(self):
        """
        Print the force constant matrix by including the print, hessian command in the input file
        """

        fc = None
        n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3

        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Read force constant matrix
                if "Force Constants (Second Derivatives of the Energy) in [a.u.]" in line:
                    fc = np.zeros((n_rows, n_rows), float)
                    for i in range(int(math.ceil(n_rows / 5.0))):
                        # Header row
                        line = f.readline()
                        # Matrix element rows
                        for j in range(i * 5, n_rows):
                            data = f.readline().split()
                            for k in range(len(data) - 1):
                                fc[j, i * 5 + k] = float(data[k + 1].replace("D", "E"))
                                fc[i * 5 + k, j] = fc[j, i * 5 + k]
                    # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                    fc *= 4.35974417e-18 / 5.291772108e-11**2
                line = f.readline()

        return fc

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Molpro .out file. If multiple such geometries are identified, only the
        last is returned.
        """

        symbol, coord, mass, number = [], [], [], []

        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Automatically determine the number of atoms
                if "Current geometry" in line:
                    symbol, coord = [], []
                    while "ENERGY" not in line:
                        line = f.readline()
                    line = f.readline()
                    while line != "\n":
                        data = line.split()
                        symbol.append(str(data[0]))
                        coord.append([float(data[1]), float(data[2]), float(data[3])])
                        line = f.readline()
                    line = f.readline()
                line = f.readline()

        # If no optimized coordinates were found, uses the input geometry
        # (for example if reading the geometry from a frequency file)
        if not coord:
            with open(self.path, "r") as f:
                line = f.readline()
                while line != "":
                    if "atomic coordinates" in line.lower():
                        symbol, coord = [], []
                        for i in range(4):
                            line = f.readline()
                        while line != "\n":
                            data = line.split()
                            symbol.append(str(data[1]))
                            coord.append([float(data[3]), float(data[4]), float(data[5])])
                            line = f.readline()
                    line = f.readline()

        # Assign appropriate mass to each atom in the molecule
        for atom1 in symbol:
            mass1, num1 = get_element_mass(atom1)
            mass.append(mass1)
            number.append(num1)
        number = np.array(number, np.int)
        mass = np.array(mass, float)
        coord = np.array(coord, float)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise LogError("Unable to read atoms from Molpro geometry output file {0}".format(self.path))

        return coord, number, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=""):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of a MolPro "Freq" quantum chemistry calculation with the thermo printed.
        """

        modes = []
        unscaled_frequencies = []
        e0 = 0.0
        if optical_isomers is None or symmetry is None:
            _optical_isomers, _symmetry, _ = self.get_symmetry_properties()
            if optical_isomers is None:
                optical_isomers = _optical_isomers
            if symmetry is None:
                symmetry = _symmetry
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Read the spin multiplicity if not explicitly given
                if spin_multiplicity == 0 and "spin" in line:
                    splits = line.replace("=", " ").replace(",", " ").split(" ")
                    for i, s in enumerate(splits):
                        if "spin" in s:
                            spin_multiplicity = int(splits[i + 1]) + 1
                            logging.debug("Conformer {0} is assigned a spin multiplicity of {1}".format(label, spin_multiplicity))
                            break
                if spin_multiplicity == 0 and "SPIN SYMMETRY" in line:
                    spin_symmetry = line.split()[-1]
                    if spin_symmetry == "Singlet":
                        spin_multiplicity = 1
                    elif spin_symmetry == "Doublet":
                        spin_multiplicity = 2
                    elif spin_symmetry == "Triplet":
                        spin_multiplicity = 3
                    elif spin_symmetry == "Quartet":
                        spin_multiplicity = 4
                    elif spin_symmetry == "Quintet":
                        spin_multiplicity = 5
                    elif spin_symmetry == "Sextet":
                        spin_multiplicity = 6
                    if spin_multiplicity:
                        logging.debug("Conformer {0} is assigned a spin multiplicity of {1}".format(label, spin_multiplicity))
                        break

                # The data we want is in the Thermochemistry section of the output
                if "THERMODYNAMICAL" in line:
                    modes = []
                    line = f.readline()
                    while line != "":
                        # This marks the end of the thermochemistry section
                        if "*************************************************" in line:
                            break

                        # Read molecular mass for external translational modes
                        elif "Molecular Mass:" in line:
                            mass = float(line.split()[2])
                            translation = IdealGasTranslation(mass=(mass, "amu"))
                            modes.append(translation)

                        # Read moments of inertia for external rotational modes
                        elif "Rotational Constants" in line and line.split()[-1] == "[GHz]":
                            inertia = [float(d) for d in line.split()[-4:-1]]
                            for i in range(3):
                                inertia[i] = constants.h / (8 * constants.pi * constants.pi * inertia[i] * 1e9) * constants.Na * 1e23
                            rotation = NonlinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                            modes.append(rotation)

                        elif "Rotational Constant" in line and line.split()[3] == "[GHz]":
                            inertia = float(line.split()[2])
                            inertia = constants.h / (8 * constants.pi * constants.pi * inertia * 1e9) * constants.Na * 1e23
                            rotation = LinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                            modes.append(rotation)

                        # Read vibrational modes
                        elif "Vibrational Temperatures" in line:
                            frequencies = []
                            frequencies.extend([float(d) for d in line.split()[3:]])
                            line = f.readline()
                            while line.strip() != "":
                                frequencies.extend([float(d) for d in line.split()])
                                line = f.readline()
                            # Convert from K to cm^-1
                            if len(frequencies) > 0:
                                frequencies = [freq * 0.695039 for freq in frequencies]  # kB = 0.695039 cm^-1/K
                                unscaled_frequencies = frequencies
                                vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                                modes.append(vibration)

                        # Read the next line in the file
                        line = f.readline()

                # Read the next line in the file
                line = f.readline()

        return (
            Conformer(E0=(e0 * 0.001, "kJ/mol"), modes=modes, spin_multiplicity=spin_multiplicity, optical_isomers=optical_isomers),
            unscaled_frequencies,
        )

    def load_energy(self, zpe_scale_factor=1.0):
        """
        Return either the f12 or MRCI energy in J/mol from a Molpro Logfile.
        If the MRCI job outputted the MRCI+Davidson energy, the latter is returned.
        For CCSD(T)-f12, the function determines which energy (f12a or f12b) to use based on the basis set,
        which it will parse out of the Molpro file. For the vdz and vtz basis sets f12a is
        a better approximation, but for higher basis sets f12b is a better approximation.
        """
        e_elect = None
        with open(self.path, "r") as f:
            lines = f.readlines()
            # Determine whether the sp method is f12,
            # if so whether we should parse f12a or f12b according to the basis set.
            # Otherwise, check whether the sp method is MRCI.
            f12, f12a, f12b, mrci = False, False, False, False
            for line in lines:
                if "basis" in line.lower():
                    if "vtz" in line.lower() or "vdz" in line.lower():
                        f12a = True  # MRCI could also have a vdz/vtz basis, so don't break yet
                    elif any(high_basis in line.lower() for high_basis in ["vqz", "v5z", "v6z", "v7z", "v8z"]):
                        f12b = True  # MRCI could also have a v(4+)z basis, so don't break yet
                elif "ccsd" in line.lower() and "f12" in line.lower():
                    f12 = True
                elif "mrci" in line.lower():
                    mrci = True
                    f12a, f12b = False, False
                    break
                elif "point group" in line.lower():
                    # We should know the method by this point, so break if possible, but don't throw an error yet
                    if any([mrci, f12a, f12b]):
                        break
            else:
                raise LogError("Could not determine type of calculation. Currently, CCSD(T)-F12a, CCSD(T)-F12b, " "MRCI, MRCI+Davidson are supported")
            # Search for e_elect
            for line in lines:
                if f12 and f12a:
                    if "CCSD(T)-F12a" in line and "energy" in line:
                        e_elect = float(line.split()[-1])
                        break
                elif f12 and f12b:
                    if "CCSD(T)-F12b" in line and "energy" in line:
                        e_elect = float(line.split()[-1])
                        break
                elif mrci:
                    # First search for MRCI+Davidson energy
                    if "(Davidson, relaxed reference)" in line:
                        e_elect = float(line.split()[3])
                        logging.debug("Found MRCI+Davidson energy in molpro log file {0}, using this value".format(self.path))
                        break
                elif not f12:
                    if "Electronic Energy at 0" in line:
                        e_elect = float(line.split()[-2])
                        break
                    if "CCSD" in line and "energy=" in line:
                        e_elect = float(line.split()[-1])
                        break
            if e_elect is None and mrci:
                # No Davidson correction is given, search for MRCI energy
                read_e_elect = False
                for line in lines:
                    if read_e_elect:
                        e_elect = float(line.split()[0])
                        logging.debug("Found MRCI energy in molpro log file {0}, using this value" " (did NOT find MRCI+Davidson)".format(self.path))
                        break
                    if all(w in line for w in ("MRCI", "MULTI", "HF-SCF")):
                        read_e_elect = True
        logging.debug("Molpro energy found is {0} Hartree".format(e_elect))
        # multiply e_elect by correct constants
        if e_elect is not None:
            e_elect *= constants.E_h * constants.Na
            logging.debug("Molpro energy found is {0} J/mol".format(e_elect))
            return e_elect
        else:
            raise LogError("Unable to find energy in Molpro log file {0}.".format(self.path))

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a MolPro log file.
        """

        zpe = None

        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Do NOT read the ZPE from the "E(ZPE)=" line, as this is the scaled version!
                # We will read in the unscaled ZPE and later multiply the scaling factor
                # from the input file
                if "Electronic Energy at 0 [K]:" in line:
                    electronic_energy = float(line.split()[5])
                    line = f.readline()
                    ee_plus_zpe = float(line.split()[5])
                    zpe = (ee_plus_zpe - electronic_energy) * constants.E_h * constants.Na
                line = f.readline()

        if zpe is not None:
            return zpe
        else:
            raise LogError(
                "Unable to find zero-point energy in Molpro log file. Make sure that the "
                "keyword {frequencies, thermo, print,thermo} is included in the input file."
            )

    def load_negative_frequency(self):
        """
        Return the negative frequency from a transition state frequency calculation in cm^-1.
        """
        freqs = []
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Read vibrational frequencies
                if "Normal Modes of imaginary frequencies" in line:
                    for i in range(3):
                        line = f.readline()
                    freqs.append(line.split()[2])
                line = f.readline()

        if len(freqs) == 1:
            return -float(freqs[0])
        elif len(freqs) > 1:
            logging.info("More than one imaginary frequency in Molpro output file {0}.".format(self.path))
            return -float(freqs[0])
        else:
            raise LogError("Unable to find imaginary frequency in Molpro output file {0}".format(self.path))

    def load_scan_energies(self):
        """
        Rotor scans are not implemented in Molpro
        """
        raise NotImplementedError("Rotor scans not implemented in Molpro")

    def get_T1_diagnostic(self):
        """
        Returns the T1 diagnostic from output log.
        If multiple occurrences exist, returns the last occurrence
        """
        with open(self.path) as f:
            log = f.readlines()

        for line in reversed(log):
            if "T1 diagnostic:  " in line:
                items = line.split()
                return float(items[-1])
        raise LogError("Unable to find T1 diagnostic in energy file: {0}".format(self.path))

    def get_D1_diagnostic(self):
        """
        Returns the D1 diagnostic from output log.
        If multiple occurrences exist, returns the last occurrence
        """
        with open(self.path) as f:
            log = f.readlines()

        for line in reversed(log):
            if "D1 diagnostic:  " in line:
                items = line.split()
                return float(items[-1])
        raise LogError("Unable to find D1 diagnostic in energy file: {0}".format(self.path))

    def load_scan_pivot_atoms(self):
        """Not implemented for Molpro"""
        raise NotImplementedError("The load_scan_pivot_atoms method is not implemented for Molpro Logs")

    def load_scan_frozen_atoms(self):
        """Not implemented for Molpro"""
        raise NotImplementedError("The load_scan_frozen_atoms method is not implemented for Molpro Logs")


register_ess_adapter("MolproLog", MolproLog)
