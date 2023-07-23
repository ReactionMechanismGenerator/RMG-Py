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
This script contains unit tests of the :mod:`rmgpy.statmech.translation` module.
"""


import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech.translation import IdealGasTranslation


class TestIdealGasTranslation:
    """
    Contains unit tests of the IdealGasTranslation class.
    """

    def setup_class(self):
        self.mass = 32.0
        self.quantum = False
        self.mode = IdealGasTranslation(
            mass=(self.mass, "amu"),
            quantum=self.quantum,
        )

    def test_get_partition_function_classical(self):
        """
        Test the IdealGasTranslation.get_partition_function() method for a
        classical translator.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([7.22597e06, 2.59130e07, 1.46586e08, 4.03944e08, 8.29217e08])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.get_partition_function(temperature)
            assert abs(q_exp - q_act) < 1e-4 * q_exp

    def test_get_heat_capacity_classical(self):
        """
        Test the IdealGasTranslation.get_heat_capacity() method using a classical
        translator.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([2.5, 2.5, 2.5, 2.5, 2.5]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.get_heat_capacity(temperature)
            assert abs(cv_exp - cv_act) < 1e-4 * cv_exp

    def test_get_enthalpy_classical(self):
        """
        Test the IdealGasTranslation.get_enthalpy() method using a classical
        translator.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([2.5, 2.5, 2.5, 2.5, 2.5]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.get_enthalpy(temperature)
            assert abs(h_exp - h_act) < 1e-4 * h_exp

    def test_get_entropy_classical(self):
        """
        Test the IdealGasTranslation.get_entropy() method using a classical
        translator.
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([18.2932, 19.5703, 21.3031, 22.3168, 23.0360]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.get_entropy(temperature)
            assert abs(s_exp - s_act) < 1e-4 * s_exp

    def test_get_sum_of_states_classical(self):
        """
        Test the IdealGasTranslation.get_sum_of_states() method using a classical
        translator.
        """
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        sum_states = self.mode.get_sum_of_states(e_list)
        dens_states = self.mode.get_density_of_states(e_list)
        for n in range(10, len(e_list)):
            assert 0.8 < np.sum(dens_states[0:n]) / sum_states[n - 1] < 1.25, "{0} != {1}".format(np.sum(dens_states[0:n]), sum_states[n])

    def test_get_density_of_states_classical(self):
        """
        Test the IdealGasTranslation.get_density_of_states() method using a
        classical translator.
        """
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        dens_states = self.mode.get_density_of_states(e_list)
        temperature = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / temperature))
        q_exp = self.mode.get_partition_function(temperature)
        assert abs(q_exp - q_act) < 1e-6 * q_exp

    def test_repr(self):
        """
        Test that an IdealGasTranslation object can be reconstructed from its
        repr() output with no loss of information.
        """
        namespace = {}
        exec("mode = {0!r}".format(self.mode), globals(), namespace)
        assert "mode" in namespace
        mode = namespace["mode"]
        assert round(abs(self.mode.mass.value - mode.mass.value), 6) == 0
        assert self.mode.mass.units == mode.mass.units
        assert self.mode.quantum == mode.quantum

    def test_pickle(self):
        """
        Test that a IdealGasTranslation object can be pickled and unpickled
        with no loss of information.
        """
        import pickle

        mode = pickle.loads(pickle.dumps(self.mode, -1))
        assert round(abs(self.mode.mass.value - mode.mass.value), 6) == 0
        assert self.mode.mass.units == mode.mass.units
        assert self.mode.quantum == mode.quantum
