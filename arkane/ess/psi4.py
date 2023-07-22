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
Arkane Psi4 module
Used to parse Psi4 output files
"""

import logging

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.common import get_element_mass, get_principal_moments_of_inertia, convert_imaginary_freq_to_negative_float
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class Psi4Log(ESSAdapter):
    """
    Represent an output file from Psi4. The attribute ``path`` refers to the
    location on disk of the Psi4 output file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays. Psi4Log is an adapter for the abstract class ESSAdapter.
    """

    def check_for_errors(self):
        """
        Checks for common errors in a Psi4 log file.
        If any are found, this method will raise an error and crash.
        """
        with open(self.path) as f:
            log = f.readlines()
        error = None
        terminated = False
        for line in reversed(log):
            if "Psi4 exiting successfully" in line:
                terminated = True
            elif "PSIO Error" in line:
                error = "I/O error"
            elif "Fatal Error" in line:
                error = "Fatal Error"
            elif "RuntimeError" in line:
                error = "runtime"
            if error is not None:
                raise LogError(f"There was an error ({error}) with the Psi4 output file {self.path} " f"due to line:\n{line}")
        if not terminated:
            raise LogError(f"Psi4 run in output file {self.path} did not successfully converged.")

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in the Psi4 output file.
        """
        n_atoms = 0
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "" and n_atoms == 0:
                if "Center              X                  Y                   Z               Mass " in line:
                    _ = f.readline()
                    line = f.readline()
                    while line != "\n":
                        n_atoms += 1
                        line = f.readline()
                    break
                line = f.readline()
        return n_atoms

    def load_force_constant_matrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates) from the
        Psi4 log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned along with a warning message directing the user
        how to obtain it.
        """
        force = None
        n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                if "Force constants in mass-weighted Cartesian coordinates" in line:
                    f.readline()
                    f_array = list()
                    while line != "\n":
                        line = f.readline()
                        # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                        f_array.extend([float(f) * 4.35974417e-18 / 5.291772108e-11**2 for f in line.replace("[", "").replace("]", "").split()])
                    force = np.array(f_array).reshape(n_rows, n_rows)
                line = f.readline()
        if force is None:
            logging.warning(
                f"Could not find a force constant matrix in the Psi4 log file {self.path}\n"
                f"To make sure Psi4 prints out the force constant matrix,"
                f'make sure to set the verbose print level in Psi4 ("set_print") to at least 3.'
            )
        return force

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Psi4 log file. If multiple such geometries are identified, only the
        last is returned.
        The final geometry in the file is returned as Standard Nuclear Orientation.
        """
        atoms, coord, number, mass = [], [], [], []

        with open(self.path) as f:
            log = f.readlines()

        geometry_flag = False
        for i in reversed(range(len(log))):
            line = log[i]
            if "Center              X                  Y                   Z               Mass" in line:
                for line in log[(i + 2) :]:
                    if line != "\n":
                        data = line.split()
                        atoms.append(data[0])
                        coord.append([float(c) for c in data[1:-1]])
                        geometry_flag = True
                    else:
                        break
                if geometry_flag:
                    break

        for atom in atoms:
            mass_, num_ = get_element_mass(atom)
            mass.append(mass_)
            number.append(num_)
        coord = np.array(coord, float)
        number = np.array(number, np.int)
        mass = np.array(mass, float)
        if any(len(param) == 0 for param in [number, coord, mass]):
            raise LogError(f"Unable to read atoms from Psi4 geometry output file {self.path}.")

        return coord, number, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=""):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of a Psi4 "Freq" quantum chemistry calculation. As
        Psi4's guess of the external symmetry number is not always correct,
        you can use the ``symmetry`` argument to substitute your own value; if
        not provided, the value in the Psi4 log file will be adopted.
        """
        freq, mmass, rot, unscaled_frequencies = [], [], [], []
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
                if spin_multiplicity == 0 and "Multiplicity =" in line:
                    spin_multiplicity = int(float(line.split()[2]))
                    logging.debug(f"Conformer {label} is assigned a spin multiplicity of {spin_multiplicity}")

                if "Harmonic Vibrational Analysis" in line:
                    frequencies = []
                    while "Thermochemistry Components" not in line:
                        if "Freq [cm^-1]" in line:
                            if len(line.split()) == 5:
                                frequencies.extend([float(convert_imaginary_freq_to_negative_float(d)) for d in line.split()[-3:]])
                            elif len(line.split()) == 4:
                                frequencies.extend([float(convert_imaginary_freq_to_negative_float(d)) for d in line.split()[-2:]])
                            elif len(line.split()) == 3:
                                frequencies.extend([float(convert_imaginary_freq_to_negative_float(d)) for d in line.split()[-1:]])
                        line = f.readline()

                    frequencies = [f for f in frequencies if f > 0.0]
                    unscaled_frequencies = frequencies
                    vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                    freq.append(vibration)
                line = f.readline()

        # Get moments of inertia from external rotational modes, given in atomic units.
        coord, number, mass = self.load_geometry()
        inertia = get_principal_moments_of_inertia(coord, numbers=number)
        inertia = list(inertia[0])
        if len(inertia):
            if any(i == 0.0 for i in inertia):
                inertia.remove(0.0)
                rot.append(LinearRotor(inertia=(inertia[0], "amu*angstrom^2"), symmetry=symmetry))
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
        Load the energy in J/mol from a Psi4 log file. Only the smallest energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value. Only DFT or HF values are supported at the moment.
        PSi4 computes the Hessian numerically, therefore many total energy values are reported.
        With large basis sets and higher levels of theory there is a high probability that
        the smallest values is the correct value.
        """
        a = list()
        with open(self.path, "r") as f:
            for line in f:
                if "Total Energy =" in line:
                    a.append(float(line.split()[3]) * constants.E_h * constants.Na)
        if not len(a):
            raise LogError(f"Unable to find energy in Psi4 output file {self.path}.")
        e_elect = min(a)
        return e_elect

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a Psi4 output file.
        """
        zpe = []
        with open(self.path, "r") as f:
            for line in f.readlines():
                if "Correction ZPE" in line:
                    zpe.append(float(line.split()[2]) * 4184)  # Convert kcal/mol to J/mol.
                    logging.debug(f"ZPE is {zpe}")
        if not len(zpe):
            raise LogError(f"Unable to find zero-point energy in Psi4 output file {self.path}.")
        return zpe[-1]

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency calculation in cm^-1.
        Since there can be many imaginary frequencies, only the first one is returned.
        """
        negative_frequencies, frequency = None, None
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                if "Harmonic Vibrational Analysis" in line:
                    frequencies = []
                    while "Thermochemistry Components" not in line:
                        if "Freq [cm^-1]" in line:
                            if len(line.split()) == 5:
                                frequencies.extend([float(convert_imaginary_freq_to_negative_float(d)) for d in line.split()[-3:]])
                            elif len(line.split()) == 4:
                                frequencies.extend([float(convert_imaginary_freq_to_negative_float(d)) for d in line.split()[-2:]])
                            elif len(line.split()) == 3:
                                frequencies.extend([float(convert_imaginary_freq_to_negative_float(d)) for d in line.split()[-1:]])
                        line = f.readline()

                    negative_frequencies = [f for f in frequencies if f < 0.0]
                line = f.readline()
        if negative_frequencies is None:
            raise LogError("Unable to find imaginary frequency in Psi4 output file {0}.".format(self.path))
        elif len(negative_frequencies) == 1:
            return negative_frequencies[0]
        else:
            logging.info("More than one imaginary frequency in Psi4 output file {0}.".format(self.path))
            return negative_frequencies[0]

    def load_scan_energies(self):
        """Not implemented for Psi4"""
        raise NotImplementedError("The load_scan_energies method is not implemented for Psi4.")

    def load_scan_pivot_atoms(self):
        """Not implemented for Psi4"""
        raise NotImplementedError("The load_scan_pivot_atoms method is not implemented for Psi4.")

    def load_scan_frozen_atoms(self):
        """Not implemented for Psi4"""
        raise NotImplementedError("The load_scan_frozen_atoms method is not implemented for Psi4.")

    def get_T1_diagnostic(self):
        """Not implemented for Psi4"""
        raise NotImplementedError("The get_T1_diagnostic method is not implemented for Psi4.")


register_ess_adapter("Psi4Log", Psi4Log)
