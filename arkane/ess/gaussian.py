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
Arkane Gaussian module
Used to parse Gaussian output files
"""

import logging
import math
import os.path

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.common import check_conformer_energy, get_element_mass
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class GaussianLog(ESSAdapter):
    """
    Represent a log file from Gaussian. The attribute `path` refers to the
    location on disk of the Gaussian log file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays. GaussianLog is an adapter for the abstract class ESSAdapter.
    """

    def check_for_errors(self):
        """
        Checks for common errors in a Gaussian log file.
        If any are found, this method will raise an error and crash.
        """
        with open(self.path, 'r') as f:
            lines = f.readlines()[-100:]
            error = None
            terminated = False
            for line in reversed(lines):
                # check for common error messages
                if 'termination' in line:
                    terminated = True
                    if 'l9999.exe' in line or 'link 9999' in line:
                        error = 'Unconverged'
                    elif 'l101.exe' in line:
                        error = 'The blank line after the coordinate section is missing, ' \
                                'or charge/multiplicity was not specified correctly.'
                    elif 'l103.exe' in line:
                        error = 'Internal coordinate error'
                    elif 'l108.exe' in line:
                        error = 'There are two blank lines between z-matrix and ' \
                                'the variables, expected only one.'
                    elif 'l202.exe' in line:
                        error = 'During the optimization process, either the standard ' \
                                'orientation or the point group of the molecule has changed.'
                    elif 'l502.exe' in line:
                        error = 'Unconverged SCF.'
                    elif 'l716.exe' in line:
                        error = 'Angle in z-matrix outside the allowed range 0 < x < 180.'
                    elif 'l906.exe' in line:
                        error = 'The MP2 calculation has failed. It may be related to pseudopotential. ' \
                                'Basis sets (CEP-121G*) that are used with polarization functions, ' \
                                'where no polarization functions actually exist.'
                    elif 'l913.exe' in line:
                        error = 'Maximum optimization cycles reached.'
                    if error:
                        raise LogError(f'There was an error ({error}) with Gaussian output file {self.path} '
                                       f'due to line:\n{line}')
                    else:
                        # no need to continue parsing if terminated without errors
                        break
            if not terminated:
                raise LogError(f'Gaussian output file {self.path} did not terminate')

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the Gaussian log file.
        """
        n_atoms = 0

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '' and n_atoms == 0:
                # Automatically determine the number of atoms
                if 'Input orientation:' in line and n_atoms == 0:
                    for i in range(5):
                        line = f.readline()
                    while '---------------------------------------------------------------------' not in line:
                        n_atoms += 1
                        line = f.readline()
                line = f.readline()

        return n_atoms

    def load_force_constant_matrix(self):
        """
        Return the force constant matrix from the Gaussian log file. The job
        that generated the log file must have the option ``iop(7/33=1)`` in
        order for the proper force constant matrix (in Cartesian coordinates)
        to be printed in the log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned.
        """
        force = None

        n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read force constant matrix
                if 'Force constants in Cartesian coordinates:' in line:
                    force = np.zeros((n_rows, n_rows), np.float64)
                    for i in range(int(math.ceil(n_rows / 5.0))):
                        # Header row
                        line = f.readline()
                        # Matrix element rows
                        for j in range(i * 5, n_rows):
                            data = f.readline().split()
                            for k in range(len(data) - 1):
                                force[j, i * 5 + k] = float(data[k + 1].replace('D', 'E'))
                                force[i * 5 + k, j] = force[j, i * 5 + k]
                    # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                    force *= 4.35974417e-18 / 5.291772108e-11 ** 2
                line = f.readline()

        return force

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Gaussian log file. If multiple such geometries are identified, only the
        last is returned.
        """
        number, coord, mass = [], [], []

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Automatically determine the number of atoms
                if 'Input orientation:' in line:
                    number, coord = [], []
                    for i in range(5):
                        line = f.readline()
                    while '---------------------------------------------------------------------' not in line:
                        data = line.split()
                        number.append(int(data[1]))
                        coord.append([float(data[3]), float(data[4]), float(data[5])])
                        line = f.readline()
                line = f.readline()

        # Assign appropriate mass to each atom in the molecule
        mass = []
        for num in number:
            mass1, _ = get_element_mass(num)
            mass.append(mass1)
        coord = np.array(coord, np.float64)
        number = np.array(number, np.int)
        mass = np.array(mass, np.float64)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise LogError('Unable to read atoms from Gaussian geometry output file {0}. '
                           'Make sure the output file is not corrupt.\nNote: if your species has '
                           '50 or more atoms, you will need to add the `iop(2/9=2000)` keyword to your '
                           'input file so Gaussian will print the input orientation geometry.'.format(self.path))

        return coord, number, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of a Gaussian "Freq" quantum chemistry calculation. As
        Gaussian's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value; if
        not provided, the value in the Gaussian log file will be adopted. In a
        log file with multiple Thermochemistry sections, only the last one will
        be kept.
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
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':

                # Read the spin multiplicity if not explicitly given
                if spin_multiplicity == 0 and 'Multiplicity =' in line:
                    spin_multiplicity = int(line.split()[-1])
                    logging.debug('Conformer {0} is assigned a spin multiplicity of {1}'
                                  .format(label, spin_multiplicity))

                # The data we want is in the Thermochemistry section of the output
                if '- Thermochemistry -' in line:
                    modes = []
                    in_partition_functions = False
                    line = f.readline()
                    while line != '':

                        # This marks the end of the thermochemistry section
                        if '-------------------------------------------------------------------' in line:
                            break

                        # Read molecular mass for external translational modes
                        elif 'Molecular mass:' in line:
                            mass = float(line.split()[2])
                            translation = IdealGasTranslation(mass=(mass, "amu"))
                            modes.append(translation)

                        # Read moments of inertia for external rotational modes
                        elif 'Rotational constants (GHZ):' in line:
                            inertia = [float(d) for d in line.split()[-3:]]
                            for i in range(3):
                                inertia[i] = constants.h / (8 * constants.pi * constants.pi * inertia[i] * 1e9) \
                                             * constants.Na * 1e23
                            rotation = NonlinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                            modes.append(rotation)
                        elif 'Rotational constant (GHZ):' in line:
                            inertia = [float(line.split()[3])]
                            inertia[0] = constants.h / (8 * constants.pi * constants.pi * inertia[0] * 1e9) \
                                         * constants.Na * 1e23
                            rotation = LinearRotor(inertia=(inertia[0], "amu*angstrom^2"), symmetry=symmetry)
                            modes.append(rotation)

                        # Read vibrational modes
                        elif 'Vibrational temperatures:' in line:
                            frequencies = []
                            frequencies.extend([float(d) for d in line.split()[2:]])
                            line = f.readline()
                            frequencies.extend([float(d) for d in line.split()[1:]])
                            line = f.readline()
                            while line.strip() != '':
                                frequencies.extend([float(d) for d in line.split()])
                                line = f.readline()
                            # Convert from K to cm^-1
                            if len(frequencies) > 0:
                                frequencies = [freq * 0.695039 for freq in frequencies]  # kB = 0.695039 cm^-1/K
                                unscaled_frequencies = frequencies
                                vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                                modes.append(vibration)

                        # Read ground-state energy
                        elif 'Sum of electronic and zero-point Energies=' in line:
                            e0 = float(line.split()[6]) * 4.35974394e-18 * constants.Na

                        # Read spin multiplicity if above method was unsuccessful
                        elif 'Electronic' in line and in_partition_functions and spin_multiplicity == 0:
                            spin_multiplicity = int(float(line.split()[1].replace('D', 'E')))

                        elif 'Log10(Q)' in line:
                            in_partition_functions = True

                        # Read the next line in the file
                        line = f.readline()

                # Read the next line in the file
                line = f.readline()

        return Conformer(E0=(e0 * 0.001, "kJ/mol"), modes=modes, spin_multiplicity=spin_multiplicity,
                         optical_isomers=optical_isomers), unscaled_frequencies

    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the energy in J/mol from a Gaussian log file. The file is checked 
        for a complete basis set extrapolation; if found, that value is 
        returned. Only the last energy in the file is returned. The zero-point
        energy is *not* included in the returned value; it is removed from the
        CBS-QB3 value.
        """
        e_elect, e0_composite, scaled_zpe = None, None, None
        elect_energy_source = ''
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':

                if 'SCF Done:' in line:
                    e_elect = float(line.split()[4]) * constants.E_h * constants.Na
                    elect_energy_source = 'SCF'
                elif ' E2(' in line and ' E(' in line:
                    e_elect = float(line.split()[-1].replace('D', 'E')) * constants.E_h * constants.Na
                    elect_energy_source = 'doublehybrd or MP2'
                elif 'MP2 =' in line:
                    e_elect = float(line.split()[-1].replace('D', 'E')) * constants.E_h * constants.Na
                    elect_energy_source = 'MP2'
                elif 'E(CORR)=' in line:
                    e_elect = float(line.split()[3]) * constants.E_h * constants.Na
                    elect_energy_source = 'CCSD'
                elif 'CCSD(T)= ' in line:
                    e_elect = float(line.split()[1].replace('D', 'E')) * constants.E_h * constants.Na
                    elect_energy_source = 'CCSD(T)'
                elif 'CBS-QB3 (0 K)' in line:
                    e0_composite = float(line.split()[3]) * constants.E_h * constants.Na
                elif 'E(CBS-QB3)=' in line:
                    # CBS-QB3 calculation without opt and freq calculation
                    # Keyword in Gaussian CBS-QB3(SP), No zero-point or thermal energies are included.
                    e_elect = float(line.split()[1]) * constants.E_h * constants.Na
                elif 'CBS-4 (0 K)=' in line:
                    e0_composite = float(line.split()[3]) * constants.E_h * constants.Na
                elif 'G3(0 K)' in line:
                    e0_composite = float(line.split()[2]) * constants.E_h * constants.Na
                elif 'G3 Energy=' in line:
                    # G3 calculation without opt and freq calculation
                    # Keyword in Gaussian G3(SP), No zero-point or thermal energies are included.
                    e_elect = float(line.split()[2]) * constants.E_h * constants.Na
                elif 'G4(0 K)' in line:
                    e0_composite = float(line.split()[2]) * constants.E_h * constants.Na
                elif 'G4 Energy=' in line:
                    # G4 calculation without opt and freq calculation
                    # Keyword in Gaussian G4(SP), No zero-point or thermal energies are included.
                    e_elect = float(line.split()[2]) * constants.E_h * constants.Na
                elif 'G4MP2(0 K)' in line:
                    e0_composite = float(line.split()[2]) * constants.E_h * constants.Na
                elif 'G4MP2 Energy=' in line:
                    # G4MP2 calculation without opt and freq calculation
                    # Keyword in Gaussian G4MP2(SP), No zero-point or thermal energies are included.
                    e_elect = float(line.split()[2]) * constants.E_h * constants.Na

                # Read the ZPE from the "E(ZPE)=" line, as this is the scaled version.
                # Gaussian defines the following as
                # E (0 K) = Elec + E(ZPE),
                # The ZPE is the scaled ZPE given by E(ZPE) in the log file,
                # hence to get the correct Elec from E (0 K) we need to subtract the scaled ZPE

                elif 'E(ZPE)' in line:
                    scaled_zpe = float(line.split()[1]) * constants.E_h * constants.Na
                elif '\\ZeroPoint=' in line:
                    line = line.strip() + f.readline().strip()
                    start = line.find('\\ZeroPoint=') + 11
                    end = line.find('\\', start)
                    scaled_zpe = float(line[start:end]) * constants.E_h * constants.Na * zpe_scale_factor
                # Read the next line in the file
                line = f.readline()

        if e0_composite is not None:
            logging.debug("Using the composite energy from the gaussian output file")
            if scaled_zpe is None:
                raise LogError('Unable to find zero-point energy in Gaussian log file.')
            return e0_composite - scaled_zpe
        elif e_elect is not None:
            logging.debug("Using the {0} energy from the gaussian output file".format(elect_energy_source))
            return e_elect
        else:
            raise LogError('Unable to find energy in Gaussian log file.')

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a Gaussian log file.
        """
        zpe = None

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Do NOT read the ZPE from the "E(ZPE)=" line, as this is the scaled version!
                # We will read in the unscaled ZPE and later multiply the scaling factor
                # from the input file
                if 'Zero-point correction=' in line:
                    zpe = float(line.split()[2]) * constants.E_h * constants.Na
                elif '\\ZeroPoint=' in line:
                    line = line.strip() + f.readline().strip()
                    start = line.find('\\ZeroPoint=') + 11
                    end = line.find('\\', start)
                    zpe = float(line[start:end]) * constants.E_h * constants.Na
                line = f.readline()

        if zpe is not None:
            return zpe
        else:
            raise LogError('Unable to find zero-point energy in Gaussian log file.')

    def load_scan_energies(self):
        """
        Extract the optimized energies in J/mol from a log file, e.g. the 
        result of a Gaussian "Scan" quantum chemistry calculation.
        """
        opt_freq = False
        rigid_scan = False

        vlist = []  # The array of potentials at each scan angle
        non_optimized = []  # The array of indexes of non-optimized point

        internal_coord = 'D(' + ','.join(self.load_scan_pivot_atoms()) + ')'
        angle = []
        # Parse the Gaussian log file, extracting the energies of each
        # optimized conformer in the scan
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # If the job contains a "freq" then we want to ignore the last energy
                if ' Freq' in line and ' Geom=' in line:
                    opt_freq = True
                # if # scan is keyword instead of # opt, then this is a rigid scan job
                # and parsing the energies is done a little differently
                if '# scan' in line:
                    rigid_scan = True
                # The lines containing "SCF Done" give the energy at each
                # iteration (even the intermediate ones)
                if 'SCF Done:' in line:
                    energy = float(line.split()[4])
                    # rigid scans will only not optimize, so just append every time it finds an energy.
                    if rigid_scan:
                        vlist.append(energy)
                # We want to keep the values of energy that come most recently before
                # the line containing "Optimization completed", since it refers
                # to the optimized geometry
                if 'Optimization completed' in line:
                    vlist.append(energy)
                # In some cases, the optimization cannot converge within the given steps.
                # Then, the geometry is not optimized. we need to exclude these values.
                if 'Optimization stopped' in line:
                    non_optimized.append(len(vlist))
                    vlist.append(energy)
                # Read the optimized angle from optimized parameters
                if internal_coord in line and 'Scan' not in line:
                    # EXAMPLE:
                    # ! D9    D(1,2,3,15)            42.4441         -DE/DX =    0.0                 !
                    angle.append(float(line.strip().split()[3]))

                line = f.readline()

        # give warning in case this assumption is not true
        if rigid_scan:
            print(f'   Assuming {os.path.basename(self.path)} is the output from a rigid scan...')
            # For rigid scans, all of the angles are evenly spaced with a constant step size
            scan_res = math.pi / 180 * self._load_scan_angle()
            angle = np.arange(0.0, scan_res * (len(vlist) - 1) + 0.00001, scan_res, np.float64)
        else:
            angle = np.array(angle, np.float64)
            # Convert -180 ~ 180 degrees to 0 ~ 2pi rads
            angle = (angle - angle[0])
            angle[angle < 0] += 360.0
            # Adjust angle[-1] to make it close to 360 degrees
            angle[-1] = angle[-1] if angle[-1] > 2 * self._load_scan_angle() else angle[-1] + 360.0
            angle = angle * math.pi / 180

        vlist = np.array(vlist, np.float64)
        # check to see if the scanlog indicates that a one of your reacting species may not be
        # the lowest energy conformer
        check_conformer_energy(vlist, self.path)

        # Adjust energies to be relative to minimum energy conformer
        # Also convert units from Hartree/particle to J/mol
        vlist -= np.min(vlist)
        vlist *= constants.E_h * constants.Na

        if opt_freq:
            vlist = vlist[:-1]
            angle = angle[:-1]

        if non_optimized:
            logging.warning(f'Scan results for angles at {angle[non_optimized]} are discarded '
                            f'due to non-converged optimization.')
            vlist = np.delete(vlist, non_optimized)
            angle = np.delete(angle, non_optimized)

        return vlist, angle

    def _load_scan_specs(self, letter_spec, get_after_letter_spec=False):
        """
        This method reads the ouptput file for optional parameters
        sent to gaussian, and returns the list of optional parameters
        as a list of tuples.

        `letter_spec` is a character used to identify whether a specification
        defines pivot atoms ('S'), frozen atoms ('F') or other attributes.

        `get_after_letter_spec` is a boolean that, if True, will return the
        parameters after letter_spec is found. If not specified or False, it will
        return the preceeding letters, which are typically the atom numbers.

        More information about the syntax can be found http://gaussian.com/opt/
        """
        output = []
        reached_input_spec_section = False
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if reached_input_spec_section:
                    terms = line.split()
                    if len(terms) == 0:
                        # finished reading specs
                        break
                    if terms[0] == 'D':
                        action_index = 5  # dihedral angle with four terms
                    elif terms[0] == 'A':
                        action_index = 4  # valance angle with three terms
                    elif terms[0] == 'B':
                        action_index = 3  # bond length with 2 terms
                    elif terms[0] == 'L':
                        # Can be either L 1 2 3 B or L 1 2 3 -1 B
                        # It defines a linear bend which is helpful in calculating
                        # molecules with ~180 degree bond angles. As no other module
                        # now depends on this information, simply skipping this line.
                        line = f.readline()
                        continue
                    else:
                        raise LogError('This file has an option not supported by Arkane. '
                                       'Unable to read scan specs for line: {0}'.format(line))
                    if len(terms) > action_index:
                        # specified type explicitly
                        if terms[action_index] == letter_spec:
                            if get_after_letter_spec:
                                output.append(terms[action_index+1:])
                            else:
                                output.append(terms[1:action_index])
                    else:
                        # no specific specification, assume freezing
                        if letter_spec == 'F':
                            output.append(terms[1:action_index])
                if " The following ModRedundant input section has been read:" in line:
                    reached_input_spec_section = True
                line = f.readline()
        return output

    def load_scan_pivot_atoms(self):
        """
        Extract the atom numbers which the rotor scan pivots around
        Return a list of atom numbers starting with the first atom as 1
        """
        output = self._load_scan_specs('S')
        return output[0] if len(output) > 0 else []

    def load_scan_frozen_atoms(self):
        """
        Extract the atom numbers which were frozen during the scan
        Return a list of list of atom numbers starting with the first atom as 1
        Each element of the outer lists represents a frozen bond
        Inner lists with length 2 represent frozen bond lengths
        Inner lists with length 3 represent frozen bond angles
        Inner lists with length 4 represent frozen dihedral angles
        """
        return self._load_scan_specs('F')

    def _load_scan_angle(self):
        """
        Return the angle difference (degrees) for a gaussian scan.
        """
        output = self._load_scan_specs('S', get_after_letter_spec=True)
        return float(output[0][1])

    def _load_number_scans(self):
        """
        Return the number of scans for a gaussian scan specified in the input file
        """
        output = self._load_scan_specs('S', get_after_letter_spec=True)
        return int(output[0][0])

    def load_negative_frequency(self):
        """
        Return the negative frequency from a transition state frequency
        calculation in cm^-1.
        """
        frequency = None
        frequencies = []
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read vibrational frequencies
                if 'Frequencies --' in line:
                    frequencies.extend(line.split()[2:])
                line = f.readline()

        frequencies = [float(freq) for freq in frequencies]
        frequencies.sort()
        try:
            frequency = [freq for freq in frequencies if freq < 0][0]
        except IndexError:
            raise LogError(f'Unable to find imaginary frequency in Gaussian output file {self.path}')
        return frequency


register_ess_adapter("GaussianLog", GaussianLog)
