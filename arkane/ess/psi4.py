#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

from arkane.common import  get_element_mass, get_principal_moments_of_inertia, convert_i_to_neg
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class Psi4Log(ESSAdapter):
    """
    Represent an output file from Psi4. The attribute `path` refers to the
    location on disk of the Psi4 output file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays. Psi4Log is an adapter for the abstract class ESSAdapter.
    """

    def __init__(self, path):
        self.path = path

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the Psi4 output file.
        """
        n_atoms = 0

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '' and n_atoms == 0:
                # Automatically determine the number of atoms
                if 'Center              X                  Y                   Z               Mass ' in line:
                    _ = f.readline()
                    line = f.readline()
                    while line != '\n':
                        n_atoms += 1
                        line = f.readline()
                line = f.readline()

        return n_atoms

    def load_force_constant_matrix(self):
        """
        make sure print level(1..5) set to 3. set_print = 3
        Return the force constant matrix (in Cartesian coordinates) from the
        Psi4 log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned.
        """
        force = None
        f_array = []

        n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read force constant matrix
                if 'Force constants in Cartesian coordinates.' in line:
                    f.readline()
                    f_array = []  # we read the last one
                    while line != '\n':
                        line = f.readline()
                        # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                        f_array.extend([float(f) * 4.35974417e-18 / 5.291772108e-11 ** 2 for f in
                                        line.replace('[', '').replace(']', '').split()])
                    force = np.array(f_array).reshape(n_rows, n_rows)
                line = f.readline()

        return force

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Psi4 log file. If multiple such geometries are identified, only the
        last is returned.
        """
        atom, coord, number, mass = [], [], [], []

        with open(self.path) as f:
            log = f.readlines()

        # First check that the Psi4 job file (not necessarily a geometry optimization)
        # has successfully completed, if not an error is thrown
        completed_job = False
        for line in reversed(log):
            if 'Psi4 exiting successfully' in line:
                logging.debug('Found a successfully completed Psi4 Job')
                completed_job = True
                break

        if not completed_job:
            raise LogError('Could not find a successfully completed Psi4 job '
                           'in Psi4 output file {0}'.format(self.path))

        # Now look for the geometry.
        # Will return the final geometry in the file under Standard Nuclear Orientation.
        geometry_flag = False
        for i in reversed(range(len(log))):
            line = log[i]
            if 'Center              X                  Y                   Z               Mass' in line:
                for line in log[(i + 2):]:
                    if line != '\n':
                        data = line.split()
                        atom.append(data[0])
                        coord.append([float(c) for c in data[1:-1]])
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
        coord = np.array(coord, np.float64)
        number = np.array(number, np.int)
        mass = np.array(mass, np.float64)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise LogError('Unable to read atoms from Psi4 geometry output file {0}.'.format(self.path))

        return coord, number, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degree of freedom data from an output file created as the result of a
        Psi4 "Freq" calculation. As Psi4's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value;
        if not provided, the value in the Psi4 output file will be adopted.
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

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read spin multiplicity if not explicitly given
                if 'Charge       =' in line:
                    charge = int(float(line.split()[2]))
                    logging.debug(
                        'Conformer {0} is assigned a charge of {1}'.format(label, charge))
                if 'Multiplicity =' in line and spin_multiplicity == 0:
                    spin_multiplicity = int(float(line.split()[2]))
                    logging.debug(
                        'Conformer {0} is assigned a spin multiplicity of {1}'.format(label, spin_multiplicity))
                # The rest of the data we want is in the Thermochemistry section of the output
                if 'Harmonic Vibrational Analysis' in line:
                    modes = []
                    frequencies = []
                    self.neg_frequencies = None
                    while 'Thermochemistry Components' not in line:
                        # This marks the end of the thermochemistry section
                        if 'Thermochemistry Components' in line:
                            break

                        # Read vibrational modes
                        if 'Freq [cm^-1]' in line:

                            if len(line.split()) == 5:
                                frequencies.extend([float(convert_i_to_neg(d)) for d in line.split()[-3:]])
                            elif len(line.split()) == 4:
                                frequencies.extend([float(convert_i_to_neg(d)) for d in line.split()[-2:]])
                            elif len(line.split()) == 3:
                                frequencies.extend([float(convert_i_to_neg(d)) for d in line.split()[-1:]])
                        line = f.readline()
                    # If there is an imaginary frequency, remove it
                    if any([f < 0.0 for f in frequencies]):
                        neg_frequencies = [f for f in frequencies if f < 0.0]
                        frequencies = [f for f in frequencies if f > 0.0]
                        self.neg_frequencies = neg_frequencies

                    unscaled_frequencies = frequencies
                    vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                    # modes.append(vibration)
                    freq.append(vibration)
                line = f.readline()

        # get moments of inertia from external rotational modes, given in atomic units        line = f.readline()
        coord, number, mass = self.load_geometry()
        inertia = get_principal_moments_of_inertia(coord, numbers=number, symbols=None)
        inertia = list(inertia[0])
        if len(inertia):
            if inertia[0] == 0.0:
                # If the first eigenvalue is 0, the rotor is linear
                inertia.remove(0.0)
                logging.debug('inertia is {}'.format(str(inertia)))
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

        translation = IdealGasTranslation(mass=(sum(mass), "amu"))
        # modes.append(translation)
        mmass.append(translation)

        # Take only the last modes found (in the event of multiple jobs)
        modes = mmass[-1:] + rot[-1:] + freq[-1:]
        return Conformer(E0=(e0 * 0.001, "kJ/mol"), modes=modes, spin_multiplicity=spin_multiplicity,
                         optical_isomers=optical_isomers), unscaled_frequencies

    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the energy in J/mol from a Psi4 log file. Only the smallest energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value. only DFT or HF values are supported at the moment.
        """
        e_elect = None
        a = []
        with open(self.path, 'r') as f:
            for line in f:
                if 'Total Energy =' in line:
                    a.append(float(line.split()[3]) * constants.E_h * constants.Na)
        # PSi4 does numerical hessian, therefore many total energy values are reported.
        # with large basis sets and higher levels of theory there is ih probability that
        # the smallest values is the best value.
        e_elect = min(a)
        if e_elect is None:
            raise LogError('Unable to find energy in Psi4 output file {0}.'.format(self.path))
        return e_elect

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a Psi4 output file.
        """
        zpe = []
        with open(self.path, 'r') as f:
            for line in f:
                if 'Correction ZPE' in line:
                    zpe.append(float(line.split()[2]) * 4184)  # Psi4's ZPE is in kcal/mol, convert to J/mol
                    logging.debug('ZPE is {}'.format(str(zpe)))
        if len(zpe) > 0:
            return zpe[-1]
        else:
            raise LogError('Unable to find zero-point energy in Psi4 output file {0}.'.format(self.path))

    def load_scan_energies(self):
        """Not implemented for Psi4"""
        raise NotImplementedError('The load_scan_pivot_atoms method is not implemented for Psi4 Logs')

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        since there can be many imaginary frequencies, we collect the last occurrence
        """
        # Make sure the frequency is imaginary:
        self.load_conformer()
        if self.neg_frequencies is not None:
            if len(self.neg_frequencies) == 1:
                return self.neg_frequencies[0]
            else:
                raise LogError('There is more than one imaginary frequency in Psi4 output file {0}.'.format(self.path))
        else:
            raise LogError('Unable to find imaginary frequency in Psi4 output file {0}.'.format(self.path))

    def load_scan_pivot_atoms(self):
        """Not implemented for Psi4"""
        raise NotImplementedError('The load_scan_pivot_atoms method is not implemented for Psi4 Logs')

    def load_scan_frozen_atoms(self):
        """Not implemented for Psi4"""
        raise NotImplementedError('The load_scan_frozen_atoms method is not implemented for Psi4 Logs')

register_ess_adapter("Psi4Log", Psi4Log)
