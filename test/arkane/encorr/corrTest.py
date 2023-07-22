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
This script contains unit tests for the :mod:`arkane.encorr.corr` module.
"""

import unittest

import numpy as np

from arkane.encorr.corr import (
    assign_frequency_scale_factor,
    get_atom_correction,
    get_bac,
)
from arkane.exceptions import AtomEnergyCorrectionError
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory


class TestCorr(unittest.TestCase):
    """
    A class for testing the functions in corr.py.
    """

    @classmethod
    def setUpClass(cls):
        cls.freq_lot = LevelOfTheory(
            method="wb97X-D3", basis="def2-TZVP", software="Q-Chem"
        )
        cls.energy_lot = LevelOfTheory(
            method="CCSD(T)-F12", basis="cc-pVDZ-F12", software="MOLPRO"
        )
        cls.composite_lot = CompositeLevelOfTheory(
            freq=cls.freq_lot, energy=cls.energy_lot
        )
        cls.lot_nonexisting = LevelOfTheory("notamethod")

    def test_get_atom_correction(self):
        """
        Test that AECs can be assigned.
        It's possible these values are refit in the future so a loose tolerance
        is used to just test that the values can be queried.
        """
        atoms = {"H": 1}
        aec = get_atom_correction(level_of_theory=self.freq_lot, atoms=atoms)
        # test value is obtained by (atom_hf['H'] - atom_thermal['H']) * 4184 - H atom_energy * 4.35974394e-18 * rmgpy.constants.Na
        test_value = 1524327
        self.assertAlmostEqual(aec, test_value, places=None, delta=1000)

        with self.assertRaises(AtomEnergyCorrectionError):
            aec = get_atom_correction(level_of_theory=self.freq_lot, atoms={"X": 1})

    def test_get_bac(self):
        """
        Test that the BACs can be assigned.
        It's possible these values are refit in the future so a loose tolerance
        is used to just test that the values can be queried.
        """
        bonds = {"H-H": 1}
        # https://github.com/ReactionMechanismGenerator/RMG-database/blob/main/input/reference_sets/main/Dihydrogen.yml#L153
        CCCBDB_coords = np.array(
            [
                [0, 0, 0],
                [0, 0, 0.7414],
            ]
        )
        nums = (1, 1)

        # test Petersson BACs
        bac_type = "p"
        bac = get_bac(
            level_of_theory=self.freq_lot,
            bonds=bonds,
            coords=CCCBDB_coords,
            nums=nums,
            bac_type=bac_type,
        )
        # test value is obtained by BAC(self.freq_lot, bac_type=bac_type).get_correction(bonds=bonds, coords=CCCBDB_coords, nums=nums).value_si
        test_value = 700
        self.assertAlmostEqual(bac, test_value, places=None, delta=100)

        # test Melius BACs
        bac_type = "m"
        bac = get_bac(
            level_of_theory=self.freq_lot,
            bonds=bonds,
            coords=CCCBDB_coords,
            nums=nums,
            bac_type=bac_type,
        )
        # test value is obtained by BAC(self.freq_lot, bac_type=bac_type).get_correction(bonds=bonds, coords=CCCBDB_coords, nums=nums).value_si
        test_value = 949
        self.assertAlmostEqual(bac, test_value, places=None, delta=100)

    def test_assign_frequency_scale_factor(self):
        """
        Test that the frequency factor can be assigned.
        It's possible these values could change in the future so a large tolerance
        is used to just test that the values can be queried.
        """
        freq_scale_factor = assign_frequency_scale_factor(None)
        self.assertAlmostEqual(freq_scale_factor, 1, places=1)

        scaling_factor = assign_frequency_scale_factor(self.lot_nonexisting)
        self.assertEqual(scaling_factor, 1)

        freq_scale_factor = assign_frequency_scale_factor(self.freq_lot)
        self.assertAlmostEqual(freq_scale_factor, 0.984, places=1)

        freq_scale_factor = assign_frequency_scale_factor(self.energy_lot)
        self.assertAlmostEqual(freq_scale_factor, 0.997, places=1)

        freq_scale_factor = assign_frequency_scale_factor(self.composite_lot)
        self.assertAlmostEqual(freq_scale_factor, 0.984, places=1)
