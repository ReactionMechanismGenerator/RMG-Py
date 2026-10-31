#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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


import numpy as np

from rmgpy import constants

from arkane.encorr.corr import (
    assign_frequency_scale_factor,
    get_atom_correction,
    get_bac,
)
from arkane.exceptions import AtomEnergyCorrectionError
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory
import pytest


class TestCorr:
    """
    A class for testing the functions in corr.py.
    """

    @classmethod
    def setup_class(cls):
        cls.freq_lot = LevelOfTheory(method="wb97X-D3", basis="def2-TZVP", software="Q-Chem")
        cls.energy_lot = LevelOfTheory(method="CCSD(T)-F12", basis="cc-pVDZ-F12", software="MOLPRO")
        cls.composite_lot = CompositeLevelOfTheory(freq=cls.freq_lot, energy=cls.energy_lot)
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
        assert abs(aec - test_value) < 1000

        with pytest.raises(AtomEnergyCorrectionError):
            aec = get_atom_correction(level_of_theory=self.freq_lot, atoms={"X": 1})

    @pytest.mark.parametrize(
        "symbol, hf298, h298_h0",
        [
            ("He", 0.0, 6.197), ("Ne", 0.0, 6.197), ("K", 89.0, 6.197), ("Ca", 177.8, 6.197),
            ("Ti", 473.0, 7.539), ("Cu", 337.4, 6.197), ("Zn", 130.40, 6.197), ("Ge", 372.0, 7.398),
            ("Br", 111.87, 6.197), ("Kr", 0.0, 6.197), ("Rb", 80.9, 6.197), ("Ag", 284.9, 6.197),
            ("Cd", 111.80, 6.197), ("Sn", 301.2, 6.215), ("I", 106.76, 6.197), ("Xe", 0.0, 6.197),
            ("Cs", 76.5, 6.197), ("Hg", 61.38, 6.197), ("Pb", 195.2, 6.197),
        ],
    )
    def test_get_atom_correction_reproduces_codata_atom_hf(self, symbol, hf298, h298_h0):
        """
        Test that the AEC of a free gas-phase atom reproduces its CODATA enthalpy of formation at 298 K.
        CODATA values (kJ/mol) are ΔfH°(298 K) and H(298 K) - H(0) of the gas-phase atom, from
        Cox, Wagman, and Medvedev, CODATA Key Values for Thermodynamics (1989).
        """
        atom_energy = -1.0  # Hartree, arbitrary: it cancels against the AEC reference term
        aec = get_atom_correction(level_of_theory=self.freq_lot, atoms={symbol: 1}, atom_energies={symbol: atom_energy})
        hf298_calc = (atom_energy * 4.35974394e-18 * constants.Na + aec) / 1000 + h298_h0  # kJ/mol
        assert abs(hf298_calc - hf298) < 0.3

    def test_br_bac_compensates_aec(self):
        """
        Test that the total AEC + BAC for Br species at a fitted level of theory is unchanged by the
        correction of the Br atom_hf double subtraction, since the Br BACs were fitted against the old AEC.
        """
        lot = LevelOfTheory(method="b3lyp2023", basis="def2tzvp", software="gaussian")
        cases = [
            ({"H": 1, "Br": 1}, {"Br-H": 1}, np.array([[0, 0, 0], [0, 0, 1.4144]]), (35, 1),
             {"p": 6760049904.469657, "m": 6760046922.099125}),
            ({"Br": 2}, {"Br-Br": 1}, np.array([[0, 0, 0], [0, 0, 2.281]]), (35, 35),
             {"p": 13517043472.255623, "m": 13517038056.646053}),
        ]
        for atoms, bonds, coords, nums, expected in cases:
            aec = get_atom_correction(level_of_theory=lot, atoms=atoms)
            for bac_type, value in expected.items():
                bac = get_bac(level_of_theory=lot, bonds=bonds, coords=coords, nums=nums, bac_type=bac_type)
                assert abs(aec + bac - value) < 1

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
        assert abs(bac - test_value) < 100

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
        assert abs(bac - test_value) < 100

    def test_assign_frequency_scale_factor(self):
        """
        Test that the frequency factor can be assigned.
        It's possible these values could change in the future so a large tolerance
        is used to just test that the values can be queried.
        """
        freq_scale_factor = assign_frequency_scale_factor(None)
        assert round(abs(freq_scale_factor - 1), 1) == 0

        scaling_factor = assign_frequency_scale_factor(self.lot_nonexisting)
        assert scaling_factor == 1

        freq_scale_factor = assign_frequency_scale_factor(self.freq_lot)
        assert round(abs(freq_scale_factor - 0.984), 1) == 0

        freq_scale_factor = assign_frequency_scale_factor(self.energy_lot)
        assert round(abs(freq_scale_factor - 0.997), 1) == 0

        freq_scale_factor = assign_frequency_scale_factor(self.composite_lot)
        assert round(abs(freq_scale_factor - 0.984), 1) == 0
