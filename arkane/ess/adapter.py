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
A module for the abstract ESSAdapter class
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional
import os
import shutil

from rmgpy.qm.qmdata import QMData
from rmgpy.qm.symmetry import PointGroupCalculator

if TYPE_CHECKING:
    from rmgpy.species import Species


class ESSAdapter(ABC):
    """
    An abstract ESS Adapter class
    """

    def __init__(self, path, check_for_errors=True):
        self.path = path
        if check_for_errors:
            self.check_for_errors()

    @abstractmethod
    def check_for_errors(self):
        """
        Checks the log file for common errors.
        Optionally runs when the class is initialized to catch errors before parsing relevant information.
        """
        pass

    @abstractmethod
    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration.
        """
        pass

    @abstractmethod
    def load_force_constant_matrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates). If multiple such matrices
        are identified, only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file, ``None`` is returned.
        """
        pass

    @abstractmethod
    def load_geometry(self):
        """
        Return the optimum geometry of the molecular configuration.
        If multiple such geometries are identified, only the last is returned.
        """
        pass

    @abstractmethod
    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, species=None, label=''):
        """
        Return the optimum geometry of the molecular configuration.
        """
        pass

    @abstractmethod
    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the energy in J/mol from a log file. Only the last energy in the file is returned.
        The zero-point energy is *not* included in the returned value.
        """
        pass

    @abstractmethod
    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a log file.
        """
        pass

    @abstractmethod
    def load_scan_energies(self):
        """
        Extract the optimized energies in J/mol from a torsional scan log file.
        """
        pass

    @abstractmethod
    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        pass

    @abstractmethod
    def load_scan_pivot_atoms(self):
        """
        Extract the atom numbers which the rotor scan pivots around.
        Return a list of atom numbers starting with the first atom as 1.
        """
        pass

    @abstractmethod
    def load_scan_frozen_atoms(self):
        """
        Extract the atom numbers which were frozen during the scan.
        Return a list of list of atom numbers starting with the first atom as 1.
        Each element of the outer lists represents a frozen bond.
        Inner lists with length 2 represent frozen bond lengths.
        Inner lists with length 3 represent frozen bond angles.
        Inner lists with length 4 represent frozen dihedral angles.
        """
        pass

    def get_software(self):
        """
        Return the name of the software. Should correspond to the class
        name without 'Log'.
        """
        return self.__class__.__name__.replace('Log', '')

    def get_D1_diagnostic(self):
        """
        Returns the D1 diagnostic from output log.
        If multiple occurrences exist, returns the last occurrence.
        Should be implemented by the relevant subclass.
        """
        raise NotImplementedError(f"get_D1_diagnostic failed for {self.path} "
                                  f"since the method is not implemented for all ESSAdapter subclasses.")

    def get_T1_diagnostic(self):
        """
        Returns the T1 diagnostic from output log.
        If multiple occurrences exist, returns the last occurrence.
        Should be implemented by the relevant subclass.
        """
        raise NotImplementedError(f"get_T1_diagnostic failed for {self.path} "
                                  f"since the method is not implemented for all ESSAdapter subclasses.")

    def get_symmetry_properties(self, species: Optional['Species'] = None):
        """
        This method uses the symmetry package from RMG's QM module
        and returns a tuple where the first element is the number
        of optical isomers, the second element is the symmetry number,
        and the third element is the point group identified.
        If ``species`` is given, then the ``Species.get_symmetry_number()`` module is called
        to determine the external symmetry.
        """
        symmetry = None
        if species is not None:
            spc = species.copy(deep=True)
            spc.generate_resonance_structures(keep_isomorphic=False, filter_structures=True, save_order=True)
            spc.get_symmetry_number(external=True)
            symmetry = spc.symmetry_number

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
                if symmetry is None:
                    symmetry = pg.symmetry_number
            else:
                raise ValueError(f'Symmetry algorithm errored when computing point group for log file located at\n'
                                 f'{self.path}\nManually provide the {"external symmetry and" if symmetry is None else ""} '
                                 f'optical isomer number.')
            return optical_isomers, symmetry, pg.point_group
        finally:
            shutil.rmtree(scr_dir)
