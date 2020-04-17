#!/usr/bin/env python3

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
from abc import ABC, abstractmethod

from rmgpy.qm.qmdata import QMData
from rmgpy.qm.symmetry import PointGroupCalculator

################################################################################


class Log(ABC):
    """
    Represent a general log file.
    The attribute `path` refers to the location on disk of the log file of interest.
    """

    @abstractmethod
    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_force_constant_matrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates).
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degree of freedom data from a frequency calculations.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the energy.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_scan_energies(self):
        """
        Extract the optimized energies from a potential energy scan.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_scan_pivot_atoms(self):
        """
        Extract the atom numbers which the rotor scan pivots around.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_scan_frozen_atoms(self):
        """
        Extract the atom numbers which were frozen during the scan.
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency calculation.
        Should be implemented by a subclass.
        """
        pass

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

    @abstractmethod
    def get_D1_diagnostic(self):
        """
        This method returns the D1 diagnostic for certain quantum jobs
        Should be implemented by a subclass.
        """
        pass

    @abstractmethod
    def get_T1_diagnostic(self):
        """
        This method returns the T1 diagnostic for certain quantum jobs
        Should be implemented by a subclass.
        """
        pass
