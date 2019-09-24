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
A general class for parsing quantum mechanical log files
"""

import logging
import os.path
import shutil

from rmgpy.qm.qmdata import QMData
from rmgpy.qm.symmetry import PointGroupCalculator

################################################################################


class Log(object):
    """
    Represent a general log file.
    The attribute `path` refers to the location on disk of the log file of interest.
    """

    def __init__(self, path):
        self.path = path

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the MolPro log file.
        """
        raise NotImplementedError("get_number_of_atoms is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def load_force_constant_matrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates) from the
        QChem log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned.
        """
        raise NotImplementedError("load_force_constant_matrix is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        log file. If multiple such geometries are identified, only the
        last is returned.
        """
        raise NotImplementedError("load_geometry is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degree of freedom data from an output file created as the result of a
        QChem "Freq" calculation. As QChem's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value;
        if not provided, the value in the QChem output file will be adopted.
        """
        raise NotImplementedError("load_conformer is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the energy in J/mol from a QChem log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        raise NotImplementedError("load_energy is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a QChem output file.
        """
        raise NotImplementedError("load_zero_point_energy is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def load_scan_energies(self):
        """
        Extract the optimized energies in J/mol from a QChem log file, e.g. the
        result of a QChem "PES Scan" quantum chemistry calculation.
        """
        raise NotImplementedError("load_scan_energies is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def load_scan_pivot_atoms(self):
        """
        Extract the atom numbers which the rotor scan pivots around
        Return a list of atom numbers starting with the first atom as 1
        """
        raise NotImplementedError("load_scan_pivot_atoms is not implemented for the Log class")

    def load_scan_frozen_atoms(self):
        """
        Extract the atom numbers which were frozen during the scan
        Return a list of list of atom numbers starting with the first atom as 1
        Each element of the outer lists represents a frozen bond
        Inner lists with length 2 represent frozen bond lengths
        Inner lists with length 3 represent frozen bond angles
        Inner lists with length 4 represent frozen dihedral angles
        """
        raise NotImplementedError("load_scan_frozen_atoms is not implemented for the Log class")

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        raise NotImplementedError("load_negative_frequency is not implemented for the Log class. "
                                  "This method should be implemented by a subclass.")

    def get_symmetry_properties(self):
        """
        This method uses the symmetry package from RMG's QM module
        and returns a tuple where the first element is the number
        of optical isomers, the second element is the symmetry number,
        and the third element is the point group identified.
        """
        coordinates, atom_numbers, _ = self.load_geometry()
        unique_id = '0'  # Just some name that the SYMMETRY code gives to one of its jobs
        # Scratch directory that the SYMMETRY code writes its files in:
        scr_dir = os.path.join(os.path.abspath('.'), str('scratch'))
        if not os.path.exists(scr_dir):
            os.makedirs(scr_dir)
        try:
            qmdata = QMData(
                groundStateDegeneracy=1,  # Only needed to check if valid QMData
                numberOfAtoms=len(atom_numbers),
                atomicNumbers=atom_numbers,
                atomCoords=(coordinates, str('angstrom')),
                energy=(0.0, str('kcal/mol'))  # Only needed to avoid error
            )
            # Dynamically create custom class to store the settings needed for the point group calculation
            # Normally, it expects an rmgpy.qm.main.QMSettings object, but we don't need all of those settings
            settings = type(str(''), (),
                            dict(symmetryPath=str('symmetry'), scratchDirectory=scr_dir))()
            pgc = PointGroupCalculator(settings, unique_id, qmdata)
            pg = pgc.calculate()
            if pg is not None:
                optical_isomers = 2 if pg.chiral else 1
                symmetry = pg.symmetry_number
                logging.debug("Symmetry algorithm found {0} optical isomers and a symmetry number of {1}".format(
                    optical_isomers, symmetry))
            else:
                logging.error('Symmetry algorithm errored when computing point group\nfor log file located at{0}.\n'
                              'Manually provide values in Arkane input.'.format(self.path))
            return optical_isomers, symmetry, pg.point_group
        finally:
            shutil.rmtree(scr_dir)

    def get_D1_diagnostic(self):
        """
        This method returns the D1 diagnostic for certain quantum jobs
        """
        raise NotImplementedError("get_D1_diagnostic is not implemented for all Log subclasses.")

    def get_T1_diagnostic(self):
        """
        This method returns the T1 diagnostic for certain quantum jobs
        """
        raise NotImplementedError("get_T1_diagnostic is not implemented for all Log subclasses.")
