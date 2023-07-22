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
This script contains unit tests of the :mod:`rmgpy.statmech.conformer` module.
"""


import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import (
    Conformer,
    HarmonicOscillator,
    HinderedRotor,
    IdealGasTranslation,
    LinearRotor,
    NonlinearRotor,
)


class TestConformer:
    """
    Contains unit tests of the :class:`Conformer` class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.ethylene = Conformer(
            E0=(0.0, "kJ/mol"),
            modes=[
                IdealGasTranslation(mass=(28.03, "amu")),
                NonlinearRotor(inertia=([3.41526, 16.6498, 20.065], "amu*angstrom^2"), symmetry=4),
                HarmonicOscillator(
                    frequencies=(
                        [
                            828.397,
                            970.652,
                            977.223,
                            1052.93,
                            1233.55,
                            1367.56,
                            1465.09,
                            1672.25,
                            3098.46,
                            3111.7,
                            3165.79,
                            3193.54,
                        ],
                        "cm^-1",
                    )
                ),
            ],
            spin_multiplicity=1,
            optical_isomers=1,
        )
        self.oxygen = Conformer(
            E0=(0.0, "kJ/mol"),
            modes=[
                IdealGasTranslation(mass=(31.99, "amu")),
                LinearRotor(inertia=(11.6056, "amu*angstrom^2"), symmetry=2),
                HarmonicOscillator(frequencies=([1621.54], "cm^-1")),
            ],
            spin_multiplicity=3,
            optical_isomers=1,
        )

        # The following data is for ethane at the CBS-QB3 level
        self.coordinates = np.array(
            [
                [0.0000, 0.0000, 0.0000],
                [-0.0000, -0.0000, 1.0936],
                [1.0430, -0.0000, -0.3288],
                [-0.4484, 0.9417, -0.3288],
                [-0.7609, -1.2051, -0.5580],
                [-0.7609, -1.2051, -1.6516],
                [-0.3125, -2.1468, -0.2292],
                [-1.8039, -1.2051, -0.2293],
            ]
        )
        self.number = np.array([6, 1, 1, 1, 6, 1, 1, 1])
        self.mass = np.array([12, 1.007825, 1.007825, 1.007825, 12, 1.007825, 1.007825, 1.007825])
        self.E0 = -93.5097
        self.conformer = Conformer(
            E0=(self.E0, "kJ/mol"),
            modes=[
                IdealGasTranslation(mass=(30.0469, "amu")),
                NonlinearRotor(inertia=([6.27071, 25.3832, 25.3833], "amu*angstrom^2"), symmetry=6),
                HarmonicOscillator(
                    frequencies=(
                        [
                            818.917,
                            819.479,
                            987.099,
                            1206.76,
                            1207.05,
                            1396,
                            1411.35,
                            1489.73,
                            1489.95,
                            1492.49,
                            1492.66,
                            2995.36,
                            2996.06,
                            3040.77,
                            3041,
                            3065.86,
                            3066.02,
                        ],
                        "cm^-1",
                    )
                ),
                HinderedRotor(
                    inertia=(1.56768, "amu*angstrom^2"),
                    symmetry=3,
                    barrier=(2.69401, "kcal/mol"),
                    quantum=False,
                    semiclassical=False,
                ),
            ],
            spin_multiplicity=1,
            optical_isomers=1,
            coordinates=(self.coordinates, "angstrom"),
            number=self.number,
            mass=(self.mass, "amu"),
        )

    def test_get_partition_function_ethylene(self):
        """
        Test the StatMech.get_partition_function() method for ethylene.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([4.05311e09, 4.19728e10, 2.82309e12, 7.51135e13, 1.16538e15])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.ethylene.get_partition_function(temperature)
            assert abs(q_exp - q_act) < 1e-4 * q_exp

    def test_get_heat_capacity_ethylene(self):
        """
        Test the StatMech.get_heat_capacity() method for ethylene.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([5.11186, 7.40447, 11.1659, 13.1221, 14.1617]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.ethylene.get_heat_capacity(temperature)
            assert round(abs(cv_exp - cv_act), 3) == 0

    def test_get_enthalpy_ethylene(self):
        """
        Test the StatMech.get_enthalpy() method for ethylene.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([4.23129, 5.04826, 7.27337, 8.93167, 10.1223]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.ethylene.get_enthalpy(temperature)
            assert abs(h_exp - h_act) < 1e-4 * h_exp

    def test_get_entropy_ethylene(self):
        """
        Test the StatMech.get_entropy() method for ethylene.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([26.3540, 29.5085, 35.9422, 40.8817, 44.8142]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.ethylene.get_entropy(temperature)
            assert round(abs(s_exp - s_act), 3) == 0

    def test_get_sum_of_states_ethylene(self):
        """
        Test the StatMech.get_sum_of_states() method for ethylene.
        """
        e_list = np.arange(0, 5000 * 11.96, 2 * 11.96)
        sum_states = self.ethylene.get_sum_of_states(e_list)
        dens_states = self.ethylene.get_density_of_states(e_list)
        for n in range(10, len(e_list)):
            assert 0.8 < np.sum(dens_states[0 : n + 1]) / sum_states[n] < 1.25, "{0} != {1}".format(np.sum(dens_states[0 : n + 1]), sum_states[n])

    def test_get_density_of_states_ethylene(self):
        """
        Test the StatMech.get_density_of_states() method for ethylene.
        """
        e_list = np.arange(0, 5000 * 11.96, 2 * 11.96)
        dens_states = self.ethylene.get_density_of_states(e_list)
        temperature = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / temperature))
        q_exp = self.ethylene.get_partition_function(temperature)
        assert abs(q_exp - q_act) < 1e-1 * q_exp

    def test_get_partition_function_oxygen(self):
        """
        Test the StatMech.get_partition_function() method for oxygen.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([1.55584e09, 9.38339e09, 1.16459e11, 5.51016e11, 1.72794e12])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.oxygen.get_partition_function(temperature)
            assert abs(q_exp - q_act) < 1e-4 * q_exp

    def test_get_heat_capacity_oxygen(self):
        """
        Test the StatMech.get_heat_capacity() method for oxygen.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([3.52538, 3.70877, 4.14751, 4.32063, 4.39392]) * constants.R
        for temperature, Cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.oxygen.get_heat_capacity(temperature)
            assert round(abs(Cv_exp - cv_act), 3) == 0

    def test_get_enthalpy_oxygen(self):
        """
        Test the StatMech.get_enthalpy() method for oxygen.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([3.50326, 3.54432, 3.75062, 3.91623, 4.02765]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.oxygen.get_enthalpy(temperature)
            assert abs(h_exp - h_act) < 1e-4 * h_exp

    def test_get_entropy_oxygen(self):
        """
        Test the StatMech.get_entropy() method for oxygen.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([24.6685, 26.5065, 29.2314, 30.9513, 32.2056]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.oxygen.get_entropy(temperature)
            assert round(abs(s_exp - s_act), 3) == 0

    def test_get_sum_of_states_oxygen(self):
        """
        Test the StatMech.get_sum_of_states() method for oxygen.
        """
        e_list = np.arange(0, 5000 * 11.96, 2 * 11.96)
        sum_states = self.oxygen.get_sum_of_states(e_list)
        dens_states = self.oxygen.get_density_of_states(e_list)
        for n in range(10, len(e_list)):
            assert 0.8 < np.sum(dens_states[0 : n + 1]) / sum_states[n] < 1.25, "{0} != {1}".format(np.sum(dens_states[0 : n + 1]), sum_states[n])

    def test_get_density_of_states_oxygen(self):
        """
        Test the StatMech.get_density_of_states() method for oxygen.
        """
        e_list = np.arange(0, 5000 * 11.96, 2 * 11.96)
        dens_states = self.oxygen.get_density_of_states(e_list)
        temperature = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / temperature))
        q_exp = self.oxygen.get_partition_function(temperature)
        assert abs(q_exp - q_act) < 1e-1 * q_exp

    def test_get_total_mass(self):
        """
        Test the Conformer.get_total_mass() method.
        """
        assert round(abs(self.conformer.get_total_mass() * constants.Na * 1000.0 - np.sum(self.mass)), 6) == 0

    def test_get_center_of_mass(self):
        """
        Test the Conformer.get_center_of_mass() method.
        """
        cm = self.conformer.get_center_of_mass()
        assert round(abs(cm[0] * 1e10 - -0.38045), 4) == 0
        assert round(abs(cm[1] * 1e10 - -0.60255), 4) == 0
        assert round(abs(cm[2] * 1e10 - -0.27900), 4) == 0

    def test_get_moment_of_inertia_tensor(self):
        """
        Test the Conformer.get_moment_of_inertia_tensor() method.
        """
        inertia = self.conformer.get_moment_of_inertia_tensor()
        assert round(abs(inertia[0, 0] * constants.Na * 1e23 - 20.65968), 4) == 0
        assert round(abs(inertia[0, 1] * constants.Na * 1e23 - -7.48115), 4) == 0
        assert round(abs(inertia[0, 2] * constants.Na * 1e23 - -3.46416), 4) == 0
        assert round(abs(inertia[1, 0] * constants.Na * 1e23 - -7.48115), 4) == 0
        assert round(abs(inertia[1, 1] * constants.Na * 1e23 - 13.53472), 4) == 0
        assert round(abs(inertia[1, 2] * constants.Na * 1e23 - -5.48630), 4) == 0
        assert round(abs(inertia[2, 0] * constants.Na * 1e23 - -3.46416), 4) == 0
        assert round(abs(inertia[2, 1] * constants.Na * 1e23 - -5.48630), 4) == 0
        assert round(abs(inertia[2, 2] * constants.Na * 1e23 - 22.84296), 4) == 0

    def test_get_principal_moments_of_inertia(self):
        """
        Test the Conformer.get_principal_moments_of_inertia() method.
        """
        inertia, axes = self.conformer.get_principal_moments_of_inertia()
        assert round(abs(inertia[0] * constants.Na * 1e23 - 6.27074), 4) == 0
        assert round(abs(inertia[1] * constants.Na * 1e23 - 25.38321), 3) == 0
        assert round(abs(inertia[2] * constants.Na * 1e23 - 25.38341), 3) == 0
        # For some reason the axes seem to jump around (positioning and signs change)
        # but the absolute values should be the same as we expect
        expected = sorted(
            [
                0.497140,
                0.610114,
                0.616938,
                0.787360,
                0.018454,
                0.616218,
                0.364578,
                0.792099,
                0.489554,
            ]
        )
        result = sorted(abs(axes).flat)
        for i, j in zip(expected, result):
            assert round(abs(i - j), 4) == 0
        return

    def test_get_internal_reduced_moment_of_inertia(self):
        """
        Test the Conformer.get_internal_reduced_moment_of_inertia() method.
        """
        inertia = self.conformer.get_internal_reduced_moment_of_inertia(pivots=[1, 5], top1=[1, 2, 3, 4])
        assert round(abs(inertia * constants.Na * 1e23 - 1.56768), 4) == 0

    def test_get_number_degrees_of_freedom(self):
        """
        Test the Conformer.get_number_degrees_of_freedom() method.
        """
        # this is for ethane:
        number_degrees_of_freedom = self.conformer.get_number_degrees_of_freedom()
        assert number_degrees_of_freedom == 24

        # this is for ethylene:
        # It doesn't check against 3 * n_atoms, because n_atoms is not declared.
        number_degrees_of_freedom = self.ethylene.get_number_degrees_of_freedom()
        assert number_degrees_of_freedom == 18

        # this is for CO
        # It doesn't check against 3 * n_atoms, because n_atoms is not declared.
        number_degrees_of_freedom = self.oxygen.get_number_degrees_of_freedom()
        assert number_degrees_of_freedom == 6
