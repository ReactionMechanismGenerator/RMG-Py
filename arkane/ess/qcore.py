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
Arkane Qcore module
Used to parse Qcore JSON output files
"""

import logging
import json

import numpy as np

from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter
from arkane.common import (
    get_principal_moments_of_inertia,
)

from typing import Any, Dict, Tuple, List, Union

################################################################################

class QcoreJSON(ESSAdapter):

    def __init__(self, path) -> None:
        self.path = path

        with open(self.path,'r') as f:
            self.data = json.load(f)

    def get_number_of_atoms(self) -> int:
        return len(np.array(self.data["numbers"]))

    def load_force_constant_matrix(self) -> np.ndarray:
        """
        To be added at later date
        """
        # n_atoms = self.get_number_of_atoms()
        # return np.asarray(self.data["force_constants"]).reshape(n_atoms, n_atoms)
        return None

    def load_geometry(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (
            np.array(self.data["coords"]),
            np.array(self.data["numbers"]),
            np.array(self.data["masses"]),
        )

    def load_conformer(
        self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=""
    ) -> Any:
        if optical_isomers is None or symmetry is None:
            _optical_isomers, _symmetry, _ = self.get_symmetry_properties()
            if optical_isomers is None:
                optical_isomers = _optical_isomers
            if symmetry is None:
                symmetry = _symmetry

        if spin_multiplicity == 0:
            spin_multiplicity = int(self.data['multiplicity'])

        vibration = HarmonicOscillator(frequencies=(self.data["frequencies"], "cm^-1"))
        translation = IdealGasTranslation(mass=(sum(self.data["masses"]), "amu"))
        coord, number, mass = self.load_geometry()
        inertia = get_principal_moments_of_inertia(coord, numbers=number, symbols=None)
        inertia = list(inertia[0])
        if len(inertia):
            if len([f for f in inertia if f == 0.0]) == 1:
                # If the first eigenvalue is 0, the rotor is linear
                inertia.remove(0.0)
                logging.debug("inertia is {}".format(str(inertia)))
                rotation = LinearRotor(
                    inertia=(inertia[0], "amu*angstrom^2"), symmetry=symmetry
                )
            else:
                rotation = NonlinearRotor(
                    inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry
                )

        modes = [translation, rotation, vibration]

        return (
            Conformer(
                E0=(self.data["E0"],"kJ/mol"),
                modes=modes,
                spin_multiplicity=spin_multiplicity,
                optical_isomers=optical_isomers,
            ),
            self.data["frequencies"],
        )

    def load_energy(self, zpe_scale_factor=1.) -> float:
        return self.data["E0"] * 1000. # J/mol

    def load_zero_point_energy(self) -> float:
        return self.data["zero_point_energy"] * 1000. # J/mol

    def load_scan_energies(self):

        raise NotImplementedError(
            "The load_scan_pivot_atoms method is not implemented for qcore"
        )

    def load_negative_frequency(self):
        raise NotImplementedError(
            "The load_negative_frequency is not implemented for qcore"
        )

    def load_scan_pivot_atoms(self):
        raise NotImplementedError(
            "The load_scan_pivot_atoms method is not implemented for qcore"
        )

    def load_scan_frozen_atoms(self):
        raise NotImplementedError(
            "The load_scan_frozen_atoms method is not implemented for qcore"
        )

register_ess_adapter("QcoreJSON", QcoreJSON)
