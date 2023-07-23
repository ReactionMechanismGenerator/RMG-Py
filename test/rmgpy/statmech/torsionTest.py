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
This script contains unit tests of the :mod:`rmgpy.statmech.torsion` module.
"""


import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech.torsion import FreeRotor, HinderedRotor

import pytest


class TestHinderedRotor:
    """
    Contains unit tests of the HinderedRotor class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 1.56764
        self.symmetry = 3
        self.barrier = 11.373
        self.quantum = True
        self.mode = HinderedRotor(
            inertia=(self.inertia, "amu*angstrom^2"),
            symmetry=self.symmetry,
            barrier=(self.barrier, "kJ/mol"),
            fourier=(
                [
                    [4.58375, 0.841648, -5702.71, 6.02657, 4.7446],
                    [0.726951, -0.677255, 0.207032, 0.553307, -0.503303],
                ],
                "J/mol",
            ),
            quantum=self.quantum,
        )
        self.freemode = FreeRotor(
            inertia=(self.inertia, "amu*angstrom^2"),
            symmetry=self.symmetry,
        )

    def test_get_rotational_constant(self):
        """
        Test getting the HinderedRotor.rotationalConstant property.
        """
        b_exp = 10.7535
        b_act = self.mode.rotationalConstant.value_si
        assert round(abs(b_exp) - abs(b_act), 4) == 0
        b_act2 = self.freemode.rotationalConstant.value_si
        assert round(abs(b_exp) - abs(b_act2), 4) == 0

    def test_set_rotational_constant(self):
        """
        Test setting the HinderedRotor.rotationalConstant property.
        """
        rotational_constant = self.mode.rotationalConstant
        rotational_constant.value_si *= 2
        self.mode.rotationalConstant = rotational_constant
        self.freemode.rotationalConstant = rotational_constant
        i_exp = 0.5 * self.inertia
        i_act = self.mode.inertia.value_si * constants.Na * 1e23
        i_act2 = self.freemode.inertia.value_si * constants.Na * 1e23
        assert round(abs(i_exp) - abs(i_act), 4) == 0
        assert round(abs(i_exp) - abs(i_act2), 4) == 0

    def test_get_potential_cosine(self):
        """
        Test the HinderedRotor.get_potential() method for a cosine potential.
        """
        self.mode.fourier = None
        phi = np.arange(0.0, 2 * constants.pi + 0.0001, constants.pi / 24.0)
        potential = np.zeros_like(phi)
        for i in range(phi.shape[0]):
            potential[i] = self.mode.get_potential(phi[i])

    def test_get_potential_fourier(self):
        """
        Test the HinderedRotor.get_potential() method for a Fourier series
        potential.
        """
        phi = np.arange(0.0, 2 * constants.pi + 0.0001, constants.pi / 24.0)
        potential = np.zeros_like(phi)
        for i in range(phi.shape[0]):
            potential[i] = self.mode.get_potential(phi[i])

    def test_get_partition_function_free(self):
        """
        Test the FreeRotor.get_partition_function() method
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.sqrt(8 * np.pi**3 * constants.kB * t_list * self.freemode.inertia.value_si) / (self.symmetry * constants.h)
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.freemode.get_partition_function(temperature)
            assert abs(q_exp) - abs(q_act) < 1e-4 * q_exp

    def test_get_partition_function_classical_cosine(self):
        """
        Test the HinderedRotor.get_partition_function() method for a cosine
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([0.741953, 1.30465, 2.68553, 3.88146, 4.91235])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.get_partition_function(temperature)
            assert abs(q_exp) - abs(q_act) < 1e-4 * q_exp

    def test_get_partition_function_classical_fourier(self):
        """
        Test the HinderedRotor.get_partition_function() method for a Fourier
        series potential in the classical limit.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([0.745526, 1.30751, 2.68722, 3.88258, 4.91315])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.get_partition_function(temperature)
            assert abs(q_exp) - abs(q_act) < 1e-4 * q_exp

    def test_get_partition_function_quantum_cosine(self):
        """
        Test the HinderedRotor.get_partition_function() method for a cosine
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([1.39947, 1.94793, 3.30171, 4.45856, 5.45188])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.get_partition_function(temperature)
            assert abs(q_exp) - abs(q_act) < 1e-4 * q_exp

    def test_get_partition_function_quantum_fourier(self):
        """
        Test the HinderedRotor.get_partition_function() method for a Fourier
        series potential in the quantum limit.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([1.39364, 1.94182, 3.29509, 4.45205, 5.44563])
        for temperature, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.get_partition_function(temperature)
            assert abs(q_exp) - abs(q_act) < 5e-4 * q_exp

    def test_get_heat_capacity_free(self):
        """
        Test the FreeRotor.get_heat_capacity() method
        """
        cv_exp = constants.R / 2.0
        t_list = np.array([300, 500, 1000, 1500, 2000])
        for temperature in t_list:
            cv_act = self.freemode.get_heat_capacity(temperature)
            assert abs(cv_exp) - abs(cv_act) < 1e-4 * cv_exp

    def test_get_heat_capacity_classical_cosine(self):
        """
        Test the HinderedRotor.get_heat_capacity() method using a cosine
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1.01741, 0.951141, 0.681919, 0.589263, 0.552028]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.get_heat_capacity(temperature)
            assert abs(cv_exp) - abs(cv_act) < 1e-4 * cv_exp

    def test_get_heat_capacity_classical_fourier(self):
        """
        Test the HinderedRotor.get_heat_capacity() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1.17682, 1.01369, 0.698588, 0.596797, 0.556293]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.get_heat_capacity(temperature)
            assert abs(cv_exp) - abs(cv_act) < 1e-4 * cv_exp

    def test_get_heat_capacity_quantum_cosine(self):
        """
        Test the HinderedRotor.get_heat_capacity() method using a cosine
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1.01271, 0.945341, 0.684451, 0.591949, 0.554087]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.get_heat_capacity(temperature)
            assert abs(cv_exp) - abs(cv_act) < 1e-4 * cv_exp

    def test_get_heat_capacity_quantum_fourier(self):
        """
        Test the HinderedRotor.get_heat_capacity() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1.01263, 0.946618, 0.685345, 0.592427, 0.554374]) * constants.R
        for temperature, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.get_heat_capacity(temperature)
            assert abs(cv_exp) - abs(cv_act) < 1e-3 * cv_exp

    def test_get_enthalpy_free(self):
        """
        Test the FreeRotor.get_enthalpy() method
        """
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = constants.R * t_list / 2.0
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.freemode.get_enthalpy(temperature)
            assert abs(h_exp) - abs(h_act) < 1e-4 * h_exp

    def test_get_enthalpy_classical_cosine(self):
        """
        Test the HinderedRotor.get_enthalpy() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([1.09556, 1.09949, 0.962738, 0.854617, 0.784333]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.get_enthalpy(temperature)
            assert abs(h_exp) - abs(h_act) < 1e-4 * h_exp

    def test_get_enthalpy_classical_fourier(self):
        """
        Test the HinderedRotor.get_enthalpy() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([1.08882, 1.09584, 0.961543, 0.854054, 0.784009]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.get_enthalpy(temperature)
            assert abs(h_exp) - abs(h_act) < 1e-4 * h_exp

    def test_get_enthalpy_quantum_cosine(self):
        """
        Test the HinderedRotor.get_enthalpy() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([0.545814, 0.727200, 0.760918, 0.717496, 0.680767]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.get_enthalpy(temperature)
            assert abs(h_exp) - abs(h_act) < 1e-4 * h_exp

    def test_get_enthalpy_quantum_fourier(self):
        """
        Test the HinderedRotor.get_enthalpy() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([0.548251, 0.728974, 0.762396, 0.718702, 0.681764]) * constants.R * t_list
        for temperature, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.get_enthalpy(temperature)
            assert abs(h_exp) - abs(h_act) < 2e-3 * h_exp

    def test_get_entropy_free(self):
        t_list = np.array([300, 500, 1000, 1500, 2000])
        pf = np.array([self.freemode.get_partition_function(temperature) for temperature in t_list])
        s_exp_list = constants.R * (np.log(pf) + 0.5)
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.freemode.get_entropy(temperature)
            assert abs(s_exp) - abs(s_act) < 1e-4 * s_exp

    def test_get_entropy_classical_cosine(self):
        """
        Test the HinderedRotor.get_entropy() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([0.797089, 1.36543, 1.95062, 2.21083, 2.37608]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.get_entropy(temperature)
            assert abs(s_exp) - abs(s_act) < 1e-4 * s_exp

    def test_get_entropy_classical_fourier(self):
        """
        Test the HinderedRotor.get_entropy() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([0.795154, 1.36396, 1.95005, 2.21055, 2.37592]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.get_entropy(temperature)
            assert abs(s_exp) - abs(s_act) < 1e-4 * s_exp

    def test_get_entropy_quantum_cosine(self):
        """
        Test the HinderedRotor.get_entropy() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([0.881906, 1.39397, 1.95536, 2.21232, 2.37673]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.get_entropy(temperature)
            assert abs(s_exp) - abs(s_act) < 1e-4 * s_exp

    def test_get_entropy_quantum_fourier(self):
        """
        Test the HinderedRotor.get_entropy() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([0.880170, 1.39260, 1.95483, 2.21207, 2.37658]) * constants.R
        for temperature, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.get_entropy(temperature)
            assert abs(s_exp) - abs(s_act) < 1e-3 * s_exp

    def test_get_sum_of_states_classical_cosine(self):
        """
        Test the HinderedRotor.get_sum_of_states() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        sum_states = self.mode.get_sum_of_states(e_list)
        dens_states = self.mode.get_density_of_states(e_list)
        for n in range(10, len(e_list)):
            assert 0.8 < np.sum(dens_states[0:n]) / sum_states[n - 1] < 1.25, "{0} != {1}".format(np.sum(dens_states[0:n]), sum_states[n])

    def test_get_sum_of_states_classical_fourier(self):
        """
        Test the HinderedRotor.get_sum_of_states() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        try:
            sum_states = self.mode.get_sum_of_states(e_list)
            assert False, "HinderedRotor.get_sum_of_states() should raise a NotImplementedError"
        except NotImplementedError:
            pass

    def test_get_sum_of_states_quantum_cosine(self):
        """
        Test the HinderedRotor.get_sum_of_states() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        sum_states = self.mode.get_sum_of_states(e_list)
        dens_states = self.mode.get_density_of_states(e_list)
        for n in range(10, len(e_list)):
            assert 0.8 < np.sum(dens_states[0:n]) / sum_states[n - 1] < 1.25, "{0} != {1}".format(np.sum(dens_states[0:n]), sum_states[n])

    def test_get_sum_of_states_quantum_fourier(self):
        """
        Test the HinderedRotor.get_sum_of_states() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        sum_states = self.mode.get_sum_of_states(e_list)
        dens_states = self.mode.get_density_of_states(e_list)
        for n in range(10, len(e_list)):
            assert 0.8 < np.sum(dens_states[0:n]) / sum_states[n - 1] < 1.25, "{0} != {1}".format(np.sum(dens_states[0:n]), sum_states[n])

    def test_get_density_of_states_classical_cosine(self):
        """
        Test the HinderedRotor.get_density_of_states() method using a classical
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        dens_states = self.mode.get_density_of_states(e_list)
        temperature = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / temperature))
        q_exp = self.mode.get_partition_function(temperature)
        assert abs(q_exp) - abs(q_act) < 1e-2 * q_exp

    def test_get_density_of_states_classical_fourier(self):
        """
        Test the HinderedRotor.get_density_of_states() method using a Fourier
        series potential in the classical limit.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        try:
            dens_states = self.mode.get_density_of_states(e_list)
            assert False, "NotImplementedError not raised by HinderedRotor.get_density_of_states()"
        except NotImplementedError:
            pass

    def test_get_density_of_states_quantum_cosine(self):
        """
        Test the HinderedRotor.get_density_of_states() method using a classical
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        dens_states = self.mode.get_density_of_states(e_list)
        temperature = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / temperature))
        q_exp = self.mode.get_partition_function(temperature)
        assert abs(q_exp) - abs(q_act) < 1e-2 * q_exp

    def test_get_density_of_states_quantum_fourier(self):
        """
        Test the HinderedRotor.get_density_of_states() method using a Fourier
        series potential in the quantum limit.
        """
        self.mode.quantum = True
        e_list = np.arange(0, 10000 * 11.96, 1 * 11.96)
        dens_states = self.mode.get_density_of_states(e_list)
        temperature = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / temperature))
        q_exp = self.mode.get_partition_function(temperature)
        assert abs(q_exp) - abs(q_act) < 1e-2 * q_exp

    def test_repr(self):
        """
        Test that a HinderedRotor object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("mode = {0!r}".format(self.mode), globals(), namespace)
        assert "mode" in namespace
        mode = namespace["mode"]
        assert round(abs(self.mode.inertia.value) - abs(mode.inertia.value), 6) == 0
        assert self.mode.inertia.units == mode.inertia.units, 6
        assert self.mode.fourier.value.shape == mode.fourier.value.shape
        for A0, A in zip(self.mode.fourier.value.flat, mode.fourier.value.flat):
            assert round(abs(A0 / A - 1.0), 6) == 0
        assert self.mode.fourier.units == mode.fourier.units
        assert round(abs(self.mode.barrier.value) - abs(mode.barrier.value), 6) == 0
        assert self.mode.barrier.units == mode.barrier.units
        assert self.mode.symmetry == mode.symmetry
        assert self.mode.quantum == mode.quantum

    def test_pickle(self):
        """
        Test that a HinderedRotor object can be pickled and unpickled with no
        loss of information.
        """
        import pickle

        mode = pickle.loads(pickle.dumps(self.mode, -1))
        assert round(abs(self.mode.inertia.value) - abs(mode.inertia.value), 6) == 0
        assert self.mode.inertia.units == mode.inertia.units
        assert self.mode.fourier.value.shape == mode.fourier.value.shape
        for A0, A in zip(self.mode.fourier.value.flat, mode.fourier.value.flat):
            assert round(abs(A0 / A - 1.0), 6) == 0
        assert self.mode.fourier.units == mode.fourier.units
        assert round(abs(self.mode.barrier.value) - abs(mode.barrier.value), 6) == 0
        assert self.mode.barrier.units == mode.barrier.units
        assert self.mode.symmetry == mode.symmetry
        assert self.mode.quantum == mode.quantum
