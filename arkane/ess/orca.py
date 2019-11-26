#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import numpy
import rmgpy.constants as constants

from arkane.common import get_element_mass
from arkane.logs.log import Log
from arkane.exceptions import LogError

################################################################################


class OrcaLog(Log):
    """
    Represent an output file from Orca. The attribute `path` refers to the
    location on disk of the Orca output file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays.
    """

    def __init__(self, path):
        super(OrcaLog, self).__init__(path)

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the Orca output file.
        """
        natoms = 0

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '' and natoms == 0:
                # Automatically determine the number of atoms
                if 'CARTESIAN COORDINATES (ANGSTROEM)' in line and natoms == 0:
                    for i in range(2):
                        line = f.readline()

                    while '---------------------------------' not in line:
                        natoms += 1
                        line = f.readline()
                        if not line.strip():
                            f.close()
                            return natoms
                line = f.readline()

    def load_force_constant_matrix(self):
        """not implemented in orca"""
        # Orca print the hessian to .hess file.  you need to provide .hess instead of .log

        raise LogError('The load_force_constant_matrix method is not implemented for Orca Logs')

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Orca log file. If multiple such geometries are identified, only the
        last is returned.
        """
        atoms, coords, numbers, mass = [], [], [], []

        with open(self.path) as f:
            log = f.readlines()

        # First check that the Orca job file (not necessarily a geometry optimization)
        # has successfully completed, if not an error is thrown
        completed_job = False
        for line in reversed(log):
            if 'ORCA TERMINATED NORMALLY' in line:
                logging.debug('Found a successfully completed Orca Job')
                completed_job = True
                break

        if not completed_job:
            raise LogError(
                'Could not find a successfully completed Orca job in Orca output file {0}'.format(self.path))

        # Now look for the geometry.
        # Will return the final geometry in the file under Standard Nuclear Orientation.
        geometry_flag = False
        for i in reversed(range(len(log))):
            line = log[i]
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                for line in log[(i + 2):]:
                    if not line.strip():
                        break
                    if '---------------------------------' not in line:
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
        coord = numpy.array(coords, numpy.float64)
        number = numpy.array(numbers, numpy.int)
        mass = numpy.array(mass, numpy.float64)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise LogError('Unable to read atoms from Orca geometry output file {0}'.format(self.path))

        return coords, numbers, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of a Orca "Freq" quantum chemistry calculation. As
        CAUTION: The rotational entropy is not quite correctly treated here
         because it includes a symmetry number that is not yet correctly
         implemented in ORCA. Orca does not provide Rotational temperatures. Only E(rot) and E(trans) energies
         Consider using another supported software (such as Gaussian or QChem) to calculate the frequencies and derive E0
        """

        raise NotImplementedError('Could not find a successfully load Orca scan:Parser is not inpemented')

    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the energy in J/ml from an Orca log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        e_elect = None
        with open(self.path, 'r') as f:
            for line in f:
                if 'FINAL SINGLE POINT ENERGY' in line:  # for all methods in Orca
                    e_elect = float(line.split()[-1])
        if e_elect is None:
            raise LogError('Unable to find energy in Orca output file.')
        return e_elect * constants.E_h * constants.Na

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a Orca output file.
        """
        zpe = None
        with open(self.path, 'r') as f:
            for line in f:
                if 'Zero point energy' in line:
                    zpe = float(line.split()[-4]) * constants.E_h * constants.Na
        if zpe is None:
            raise LogError('Unable to find zero-point energy in Orca output file.')
        return zpe

    def load_scan_energies(self):
        """not implemented in Orca"""

        raise NotImplementedError('Could not find a successfully load Orca scan:Parser is not implemented')

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        frequency = 0
        with open(self.path, 'r') as f:
            for line in f:
                # Read imaginary frequency
                if '***imaginary mode***' in line:
                    frequency = float((line.split()[1]))
                    break
        # Make sure the frequency is imaginary:
        if frequency < 0:
            return frequency
        else:
            raise LogError('Unable to find imaginary frequency in Orca output file {0}'.format(self.path))

    def load_scan_pivot_atoms(self):
        """Not implemented for Orca"""
        raise NotImplementedError('The load_scan_pivot_atoms method is not implemented for Orca Logs')

    def load_scan_frozen_atoms(self):
        """Not implemented for Orca"""
        raise NotImplementedError('The load_scan_frozen_atoms method is not implemented for Orca Logs')

    def get_D1_diagnostic(self):
        """Not implemented for Orca"""
        raise NotImplementedError('The get_D1_diagnostic method is not implemented for Orca Logs')

    def get_T1_diagnostic(self):
        """
        Returns the T1 diagnostic from output log.
        If multiple occurrences exist, returns the last occurrence
        """
        with open(self.path) as f:
            log = f.readlines()

        for line in reversed(log):
            if 'T1 diagnostic ' in line:
                items = line.split()
                return float(items[-1])
        raise LogError('Unable to find T1 diagnostic in energy file: {}'.format(self.path))
