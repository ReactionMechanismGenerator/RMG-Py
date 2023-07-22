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
Arkane QChem module
Used to parse QChem output files
"""

import math
import logging
import os.path

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.common import check_conformer_energy, get_element_mass
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class QChemLog(ESSAdapter):
    """
    Represent an output file from QChem. The attribute `path` refers to the
    location on disk of the QChem output file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays. QChemLog is an adapter for the abstract class ESSAdapter.
    """

    def check_for_errors(self):
        """
        Checks for common errors in a QChem log file.
        If any are found, this method will raise an error and crash.
        """
        with open(os.path.join(self.path), "r") as f:
            lines = f.readlines()
            error, warning_line, warning, warn_message = None, None, None, None
            for line in reversed(lines):
                # check for common error messages
                if "SCF failed" in line:
                    error = "SCF failed"
                    break
                elif (
                    "error" in line and "DIIS" not in line and "gprntSymmMtrx" not in line and "Relative error" not in line and "zonesort" not in line
                ):
                    # these are **normal** lines that we should not capture:
                    # "SCF converges when DIIS error is below 1.0E-08", or
                    # "Cycle       Energy         DIIS Error" or
                    # "gprntSymmMtrx error report: "

                    # "Relative error" is captured as a warning later

                    # If running on NERSC, also want to avoid lines containing 'zonesort'
                    # Per the NERSC documenation:
                    # The error message means that running the zonesort kernel failed for some reason.
                    # The end result is that your code may have run less optimally.
                    # Other than that, the message is usually harmless.
                    # https://docs.nersc.gov/jobs/errors/
                    error = "SCF failed"
                    break
                elif "Invalid charge/multiplicity combination" in line:
                    error = "Invalid charge/multiplicity combination"
                    break
                elif "MAXIMUM OPTIMIZATION CYCLES REACHED" in line:
                    error = "Maximum optimization cycles reached."
                    break
                elif "Relative error" in line:
                    warn_message = """ 
                    Per the QChem version 5 documentation: https://manual.q-chem.com/pdf/qchem_manual_5.0.pdf
                    A warning message is printed whenever the relative error in the numerical electron count
                    reaches 0.01%, indicating that the numerical XC results may not be reliable. If the warning
                    appears on the first SCF cycle, it is probably not serious, because the initial-guess density
                    matrix is sometimes not idempotent.
                    """
                    warning = "Relative error"
                    warning_line = line
            if error:
                raise LogError(f"There was an error ({error}) with QChem output file {self.path} " f"due to line:\n{line}")
            if warning:
                logging.warning(f"{warning} with QChem output file {self.path} due to line:\n" f"{warning_line}\n" f"{warn_message}")

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the QChem output file.
        """
        n_atoms = 0

        with open(self.path, "r") as f:
            line = f.readline()
            while line != "" and n_atoms == 0:
                # Automatically determine the number of atoms
                if "Standard Nuclear Orientation" in line and n_atoms == 0:
                    for i in range(3):
                        line = f.readline()
                    while "----------------------------------------------------" not in line:
                        n_atoms += 1
                        line = f.readline()
                line = f.readline()

        return n_atoms

    def load_force_constant_matrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates) from the
        QChem log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned.
        """
        force = None

        n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3
        with open(self.path, "r") as f:
            line = f.readline()
            while line != "":
                # Read force constant matrix
                if "Final Hessian." in line or "Hessian of the SCF Energy" in line:
                    force = np.zeros((n_rows, n_rows), float)
                    for i in range(int(math.ceil(n_rows / 6.0))):
                        # Header row
                        line = f.readline()
                        # Matrix element rows
                        for j in range(n_rows):  # for j in range(i*6, Nrows):
                            data = f.readline().split()
                            for k in range(len(data) - 1):
                                force[j, i * 6 + k] = float(data[k + 1])
                                # F[i*5+k,j] = F[j,i*5+k]
                    # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                    force *= 4.35974417e-18 / 5.291772108e-11**2
                line = f.readline()

        return force

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        QChem log file. If multiple such geometries are identified, only the
        last is returned.
        """
        atom, coord, number, mass = [], [], [], []

        with open(self.path) as f:
            log = f.readlines()

        # Now look for the geometry.
        # Will return the final geometry in the file under Standard Nuclear Orientation.
        geometry_flag = False
        for i in reversed(range(len(log))):
            line = log[i]
            if "Standard Nuclear Orientation" in line:
                for line in log[(i + 3) :]:
                    if "------------" not in line:
                        data = line.split()
                        atom.append(data[1])
                        coord.append([float(c) for c in data[2:]])
                        geometry_flag = True
                    else:
                        break
                if geometry_flag:
                    break

        # Assign appropriate mass to each atom in the molecule
        for atom1 in atom:
            mass1, num1 = get_element_mass(atom1)
            mass.append(mass1)
            number.append(num1)
        coord = np.array(coord, float)
        number = np.array(number, int)
        mass = np.array(mass, float)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise LogError("Unable to read atoms from QChem geometry output file {0}.".format(self.path))

        return coord, number, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=""):
        """
        Load the molecular degree of freedom data from an output file created as the result of a
        QChem "Freq" calculation. As QChem's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value;
        if not provided, the value in the QChem output file will be adopted.
        """
        modes = []
        freq = []
        mmass = []
        rot = []
        inertia = []
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
                # Read spin multiplicity if not explicitly given
                if "$molecule" in line and spin_multiplicity == 0:
                    line = f.readline()
                    if len(line.split()) == 2:
                        spin_multiplicity = int(float(line.split()[1]))
                        logging.debug("Conformer {0} is assigned a spin multiplicity of {1}".format(label, spin_multiplicity))
                # The rest of the data we want is in the Thermochemistry section of the output
                elif "VIBRATIONAL ANALYSIS" in line:
                    modes = []
                    line = f.readline()
                    while line != "":
                        # This marks the end of the thermochemistry section
                        if "Thank you very much for using Q-Chem." in line:
                            break

                        # Read vibrational modes
                        elif "VIBRATIONAL FREQUENCIES (CM**-1)" in line:
                            frequencies = []
                            while "STANDARD THERMODYNAMIC QUANTITIES AT" not in line:
                                if " Frequency:" in line:
                                    if len(line.split()) == 4:
                                        frequencies.extend([float(d) for d in line.split()[-3:]])
                                    elif len(line.split()) == 3:
                                        frequencies.extend([float(d) for d in line.split()[-2:]])
                                    elif len(line.split()) == 2:
                                        frequencies.extend([float(d) for d in line.split()[-1:]])
                                line = f.readline()
                            line = f.readline()
                            # If there is an imaginary frequency, remove it
                            if frequencies[0] < 0.0:
                                frequencies = frequencies[1:]

                            unscaled_frequencies = frequencies
                            vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                            # modes.append(vibration)
                            freq.append(vibration)
                        # Read molecular mass for external translational modes
                        elif "Molecular Mass:" in line:
                            mass = float(line.split()[2])
                            translation = IdealGasTranslation(mass=(mass, "amu"))
                            # modes.append(translation)
                            mmass.append(translation)

                        # Read moments of inertia for external rotational modes, given in atomic units
                        elif "Eigenvalues --" in line:
                            inertia = [float(d) for d in line.split()[-3:]]

                        # Read the next line in the file
                        line = f.readline()

                # Read the next line in the file
                line = f.readline()

                if len(inertia):
                    if inertia[0] == 0.0:
                        # If the first eigenvalue is 0, the rotor is linear
                        inertia.remove(0.0)
                        logging.debug("inertia is {}".format(str(inertia)))
                        for i in range(2):
                            inertia[i] *= (constants.a0 / 1e-10) ** 2
                        inertia = np.sqrt(inertia[0] * inertia[1])
                        rotation = LinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                        rot.append(rotation)
                    else:
                        for i in range(3):
                            inertia[i] *= (constants.a0 / 1e-10) ** 2
                            rotation = NonlinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                            # modes.append(rotation)
                        rot.append(rotation)

                    inertia = []
        # Take only the last modes found (in the event of multiple jobs)
        modes = mmass[-1:] + rot[-1:] + freq[-1:]
        return (
            Conformer(E0=(e0 * 0.001, "kJ/mol"), modes=modes, spin_multiplicity=spin_multiplicity, optical_isomers=optical_isomers),
            unscaled_frequencies,
        )

    def load_energy(self, zpe_scale_factor=1.0):
        """
        Load the energy in J/mol from a QChem log file. Prioritize the energy from a converged
        geometry optimization. If the file does not contain an optimization job or if the optimization
        hit the maximum cycles, return the next equivalent source, such as from a frequency job.
        The zero-point energy is *not* included in the returned value.
        """
        e_elect = None
        with open(self.path, "r") as f:
            preferred_source = alternative_source = None
            for line in f:
                if "Final energy is" in line:
                    preferred_source = float(line.split()[-1]) * constants.E_h * constants.Na
                if "Total energy in the final basis set" in line:
                    alternative_source = float(line.split()[-1]) * constants.E_h * constants.Na
                e_elect = preferred_source or alternative_source
        if e_elect is None:
            raise LogError("Unable to find energy in QChem output file {0}.".format(self.path))
        return e_elect

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a QChem output file.
        """
        zpe = None
        with open(self.path, "r") as f:
            for line in f:
                if "Zero point vibrational energy" in line:
                    zpe = float(line.split()[4]) * 4184  # QChem's ZPE is in kcal/mol, convert to J/mol
                    logging.debug("ZPE is {}".format(str(zpe)))
        if zpe is not None:
            return zpe
        else:
            raise LogError("Unable to find zero-point energy in QChem output file {0}.".format(self.path))

    def load_scan_energies(self):
        """
        Extract the optimized energies in J/mol from a QChem log file, e.g. the
        result of a QChem "PES Scan" quantum chemistry calculation.
        """
        v_list = []
        angle = []
        read = False
        with open(self.path, "r") as f:
            for line in f:
                if "-----------------" in line:
                    read = False
                if read:
                    values = [float(item) for item in line.split()]
                    angle.append(values[0])
                    v_list.append(values[1])
                if "Summary of potential scan:" in line:
                    logging.info("found a successfully completed QChem Job")
                    read = True
        logging.info("   Assuming {0} is the output from a QChem PES scan...".format(os.path.basename(self.path)))

        v_list = np.array(v_list, float)
        # check to see if the scanlog indicates that one of your reacting species may not be the lowest energy conformer
        check_conformer_energy(v_list, self.path)

        # Adjust energies to be relative to minimum energy conformer
        # Also convert units from Hartree/particle to J/mol
        v_list -= np.min(v_list)
        v_list *= constants.E_h * constants.Na
        angle = np.arange(0.0, 2 * math.pi + 0.00001, 2 * math.pi / (len(v_list) - 1), float)
        return v_list, angle

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        read_freqs = False
        num_freq_blocks = 0
        with open(self.path, "r") as f:
            for line in f:
                # only read the first row from a frequency block
                if read_freqs and " Frequency:" in line:
                    freqs = np.array([float(freq) for freq in line.split()[1:4]])
                    num_freq_blocks += 1
                    read_freqs = False

                if "VIBRATIONAL ANALYSIS" in line:
                    read_freqs = True

        logging.info(f"Identified {num_freq_blocks} frequency block(s)...")
        logging.info("Only examining frequencies from the last block...")
        neg_idx = np.where(freqs < 0)[0]
        if len(neg_idx) == 1:
            return freqs[neg_idx[0]]
        elif len(neg_idx) > 1:
            logging.info("More than one imaginary frequency in QChem output file {0}.".format(self.path))
            return freqs[neg_idx[0]]
        else:
            raise LogError("Unable to find imaginary frequency in QChem output file {0}.".format(self.path))

    def load_scan_pivot_atoms(self):
        """Not implemented for QChem"""
        raise NotImplementedError("The load_scan_pivot_atoms method is not implemented for QChem Logs")

    def load_scan_frozen_atoms(self):
        """Not implemented for QChem"""
        raise NotImplementedError("The load_scan_frozen_atoms method is not implemented for QChem Logs")


register_ess_adapter("QChemLog", QChemLog)
