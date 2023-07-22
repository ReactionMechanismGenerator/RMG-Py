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
This script contains unit tests of the :mod:`rmgpy.statmech.vibration` module.
"""


import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech.vibration import HarmonicOscillator


class TestHarmonicOscillator:
    """
    Contains unit tests of the HarmonicOscillator class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.frequencies = np.array([500, 1000, 2000])
        self.quantum = True
        self.mode = HarmonicOscillator(
            frequencies=(self.frequencies, "cm^-1"),
            quantum=self.quantum,
        )

    def test_get_partition_function_classical(self):
        """
        Test the HarmonicOscillator.get_partition_function() method for a set of
        classical oscillators.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([0.00906536, 0.04196925, 0.335754, 1.13316978, 2.68603])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.get_partition_function(temperature)
            assert abs(q_exp - q_act) < 1e-4 * q_exp

    def test_get_partition_function_quantum(self):
        """
        Test the HarmonicOscillator.get_partition_function() method for a set of
        quantum oscillators.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([1.10923, 1.39358, 2.70819, 4.98825, 8.459780])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.get_partition_function(temperature)
            assert abs(q_exp - q_act) < 1e-4 * q_exp

    def test_get_heat_capacity_classical(self):
        """
        Test the HarmonicOscillator.get_heat_capacity() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([3, 3, 3, 3, 3]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.get_heat_capacity(temperature)
            assert abs(cv_exp - cv_act) < 1e-4 * cv_exp

    def test_get_heat_capacity_quantum(self):
        """
        Test the HarmonicOscillator.get_heat_capacity() method using a set of
        quantum oscillators.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([0.832004, 1.47271, 2.32513, 2.65024, 2.79124]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.get_heat_capacity(temperature)
            assert abs(cv_exp - cv_act) < 1e-4 * cv_exp

    def test_get_enthalpy_classical(self):
        """
        Test the HarmonicOscillator.get_enthalpy() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([3, 3, 3, 3, 3]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.get_enthalpy(temperature)
            assert abs(h_exp - h_act) < 1e-4 * h_exp

    def test_get_enthalpy_quantum(self):
        """
        Test the HarmonicOscillator.get_enthalpy() method using a set of quantum
        oscillators.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([0.280395, 0.637310, 1.30209, 1.70542, 1.96142]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.get_enthalpy(temperature)
            assert abs(h_exp - h_act) < 1e-4 * h_exp

    def test_get_entropy_classical(self):
        """
        Test the HarmonicOscillator.get_entropy() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([-1.70329, -0.170818, 1.90862, 3.12502, 3.98807]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.get_entropy(temperature)
            assert abs(s_exp - s_act) < 1e-4 * abs(s_exp)

    def test_get_entropy_quantum(self):
        """
        Test the HarmonicOscillator.get_entropy() method using a set of quantum
        oscillators.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([0.384065, 0.969182, 2.29837, 3.31251, 4.09675]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.get_entropy(temperature)
            assert abs(s_exp - s_act) < 1e-4 * s_exp

    def test_get_sum_of_states_classical(self):
        """
        Test the HarmonicOscillator.get_sum_of_states() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        self.mode.frequencies = ([500, 1000], "cm^-1")
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        sum_states = self.mode.get_sum_of_states(e_list)
        dens_states = self.mode.get_density_of_states(e_list)
        for n in range(10, len(e_list)):
            assert 0.8 < np.sum(dens_states[0:n]) / sum_states[n] < 1.25, "{0} != {1}".format(np.sum(dens_states[0:n]), sum_states[n])

    def test_get_sum_of_states_quantum(self):
        """
        Test the HarmonicOscillator.get_sum_of_states() method using a set of
        quantum oscillators.
        """
        self.mode.quantum = True
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        sum_states = self.mode.get_sum_of_states(e_list)
        dens_states = self.mode.get_density_of_states(e_list)
        for n in range(1, len(e_list)):
            if sum_states[n - 1] == 0:
                assert np.sum(dens_states[0:n]) == 0
            else:
                assert 0.8 < np.sum(dens_states[0:n]) / sum_states[n - 1] < 1.25, "{0} != {1}".format(np.sum(dens_states[0:n]), sum_states[n])

    def test_get_density_of_states_classical(self):
        """
        Test the HarmonicOscillator.get_density_of_states() method using a set of
        classical oscillators.
        """
        self.mode.quantum = False
        factor = constants.h * constants.c * 100.0 * constants.Na  # cm^-1 to J/mol
        e_list = np.arange(0, 10000 * factor, 1 * factor)
        dens_states = self.mode.get_density_of_states(e_list)
        temperature = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / temperature))
        q_exp = self.mode.get_partition_function(temperature)
        assert abs(q_exp - q_act) < 1e-4 * q_exp

    def test_get_density_of_states_quantum(self):
        """
        Test the HarmonicOscillator.get_density_of_states() method using a set of
        quantum oscillators.
        """
        self.mode.quantum = True
        factor = constants.h * constants.c * 100.0 * constants.Na  # cm^-1 to J/mol
        e_list = np.arange(0, 10000 * factor, 1 * factor)
        dens_states = self.mode.get_density_of_states(e_list)
        for n in range(len(e_list)):
            if dens_states[n] != 0:
                # The peaks should occur near a multiple of 500 cm^-1
                energy = float(e_list[n]) / factor
                assert energy % 500 < 5 or energy % 500 > 495

    def test_repr(self):
        """
        Test that a HarmonicOscillator object can be reconstructed from its
        repr() output with no loss of information.
        """
        namespace = {}
        exec("mode = {0!r}".format(self.mode), globals(), namespace)
        assert "mode" in namespace
        mode = namespace["mode"]
        assert self.mode.frequencies.value.shape == mode.frequencies.value.shape
        for freq0, freq in zip(self.mode.frequencies.value, mode.frequencies.value):
            assert round(abs(freq0 - freq), 6) == 0
        assert self.mode.frequencies.units == mode.frequencies.units
        assert self.mode.quantum == mode.quantum

    def test_pickle(self):
        """
        Test that a HarmonicOscillator object can be pickled and unpickled
        with no loss of information.
        """
        import pickle

        mode = pickle.loads(pickle.dumps(self.mode, -1))
        assert self.mode.frequencies.value.shape == mode.frequencies.value.shape
        for freq0, freq in zip(self.mode.frequencies.value, mode.frequencies.value):
            assert round(abs(freq0 - freq), 6) == 0
        assert self.mode.frequencies.units == mode.frequencies.units
        assert self.mode.quantum == mode.quantum
