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
This script contains unit tests of the :mod:`rmgpy.statmech.rotation` module.
"""

from __future__ import division
import unittest
import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor

################################################################################


class TestLinearRotor(unittest.TestCase):
    """
    Contains unit tests of the LinearRotor class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75
        self.symmetry = 2
        self.quantum = False
        self.mode = LinearRotor(
            inertia=(self.inertia, "amu*angstrom^2"),
            symmetry=self.symmetry,
            quantum=self.quantum,
        )

    def test_get_rotational_constant(self):
        """
        Test getting the LinearRotor.rotationalConstant property.
        """
        b_exp = 1.434692
        b_act = self.mode.rotationalConstant.value_si
        self.assertAlmostEqual(b_exp, b_act, 4)

    def test_set_rotational_constant(self):
        """
        Test setting the LinearRotor.rotationalConstant property.
        """
        rotational_constant = self.mode.rotationalConstant
        rotational_constant.value_si *= 2
        self.mode.rotationalConstant = rotational_constant
        i_exp = 0.5 * self.inertia
        i_act = self.mode.inertia.value_si * constants.Na * 1e23
        self.assertAlmostEqual(i_exp, i_act, 4)

    def test_get_level_energy(self):
        """
        Test the LinearRotor.getLevelEnergy() method.
        """
        rotational_constant = self.mode.rotationalConstant.value_si * constants.h * constants.c * 100.
        rotational_constant *= constants.Na
        for J in range(0, 100):
            e_exp = rotational_constant * J * (J + 1)
            e_act = self.mode.getLevelEnergy(J)
            if J == 0:
                self.assertEqual(e_act, 0)
            else:
                self.assertAlmostEqual(e_exp, e_act, delta=1e-4*e_exp)

    def test_get_level_degeneracy(self):
        """
        Test the LinearRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            g_exp = 2 * J + 1
            g_act = self.mode.getLevelDegeneracy(J)
            self.assertEqual(g_exp, g_act)

    def test_get_partition_function_classical(self):
        """
        Test the LinearRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([72.6691, 121.115, 242.230, 363.346, 484.461])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-4*q_exp)

    def test_get_partition_function_quantum(self):
        """
        Test the LinearRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([72.8360, 121.282, 242.391, 363.512, 484.627])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-4*q_exp)

    def test_get_heat_capacity_classical(self):
        """
        Test the LinearRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1, 1, 1, 1, 1]) * constants.R
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(cv_exp, cv_act, delta=1e-4*cv_exp)

    def test_get_heat_capacity_quantum(self):
        """
        Test the LinearRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1, 1, 1, 1, 1]) * constants.R
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(cv_exp, cv_act, delta=1e-4*cv_exp)

    def test_get_enthalpy_classical(self):
        """
        Test the LinearRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([1, 1, 1, 1, 1]) * constants.R * t_list
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(h_exp, h_act, delta=1e-4*h_exp)

    def test_get_enthalpy_quantum(self):
        """
        Test the LinearRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array(
            [0.997705, 0.998624, 0.999312, 0.999541, 0.999656]) * constants.R * t_list
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(h_exp, h_act, delta=1e-4*h_exp)

    def test_get_entropy_classical(self):
        """
        Test the LinearRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([5.28592, 5.79674, 6.48989, 6.89535, 7.18304]) * constants.R
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.getEntropy(T)
            self.assertAlmostEqual(s_exp, s_act, delta=1e-4*s_exp)

    def test_get_entropy_quantum(self):
        """
        Test the LinearRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([5.28592, 5.79674, 6.48989, 6.89535, 7.18304]) * constants.R
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.getEntropy(T)
            self.assertAlmostEqual(s_exp, s_act, delta=1e-4*s_exp)

    def test_get_sum_of_states_classical(self):
        """
        Test the LinearRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 2000 * 11.96, 1.0 * 11.96)
        dens_states = self.mode.getDensityOfStates(e_list)
        sum_states = self.mode.getSumOfStates(e_list)
        for n in range(1, len(e_list)):
            self.assertAlmostEqual(np.sum(dens_states[0:n]) / sum_states[n], 1.0, 3)

    def test_get_sum_of_states_quantum(self):
        """
        Test the LinearRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        e_list = np.arange(0, 4000. * 11.96, 2.0 * 11.96)
        dens_states = self.mode.getDensityOfStates(e_list)
        sum_states = self.mode.getSumOfStates(e_list)
        for n in range(1, len(e_list)):
            self.assertAlmostEqual(np.sum(dens_states[0: n + 1]) / sum_states[n], 1.0, 3)

    def test_get_dsensity_of_states_classical(self):
        """
        Test the LinearRotor.getDensityOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 400, 500])
        e_list = np.arange(0, 4000. * 11.96, 1.0 * 11.96)
        for T in t_list:
            dens_states = self.mode.getDensityOfStates(e_list)
            q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
            q_exp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-2 * q_exp)

    def test_get_dsensity_of_states_quantum(self):
        """
        Test the LinearRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 400, 500])
        e_list = np.arange(0, 4000. * 11.96, 2.0 * 11.96)
        for T in t_list:
            dens_states = self.mode.getDensityOfStates(e_list)
            q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
            q_exp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-2 * q_exp)

    def test_repr(self):
        """
        Test that a LinearRotor object can be reconstructed from its repr()
        output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

    def test_pickle(self):
        """
        Test that a LinearRotor object can be pickled and unpickled with no
        loss of information.
        """
        import pickle
        mode = pickle.loads(pickle.dumps(self.mode, -1))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

################################################################################


class TestNonlinearRotor(unittest.TestCase):
    """
    Contains unit tests of the NonlinearRotor class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = np.array([3.415, 16.65, 20.07])
        self.symmetry = 4
        self.quantum = False
        self.mode = NonlinearRotor(
            inertia=(self.inertia, "amu*angstrom^2"),
            symmetry=self.symmetry,
            quantum=self.quantum,
        )

    def test_get_rotational_constant(self):
        """
        Test getting the NonlinearRotor.rotationalConstant property.
        """
        b_exp = np.array([4.93635, 1.0125, 0.839942])
        b_act = self.mode.rotationalConstant.value_si
        for rotational_constant0, rotational_constant in zip(b_exp, b_act):
            self.assertAlmostEqual(rotational_constant0, rotational_constant, 4)

    def test_set_rotational_constant(self):
        """
        Test setting the NonlinearRotor.rotationalConstant property.
        """
        rotational_constant = self.mode.rotationalConstant
        rotational_constant.value_si *= 2
        self.mode.rotationalConstant = rotational_constant
        i_exp = 0.5 * self.inertia
        i_act = self.mode.inertia.value_si * constants.Na * 1e23
        for I0, I in zip(i_exp, i_act):
            self.assertAlmostEqual(I0, I, 4)

    def test_get_partition_function_classical(self):
        """
        Test the NonlinearRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([651.162, 1401.08, 3962.84, 7280.21, 11208.6])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-4*q_exp)

    def test_get_heat_capacity_classical(self):
        """
        Test the NonlinearRotor.getHeatCapacity() method using a classical
        rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(cv_exp, cv_act, delta=1e-4*cv_exp)

    def test_get_enthalpy_classical(self):
        """
        Test the NonlinearRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R * t_list
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(h_exp, h_act, delta=1e-4*h_exp)

    def test_get_entropy_classical(self):
        """
        Test the NonlinearRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([7.97876, 8.74500, 9.78472, 10.3929, 10.8244]) * constants.R
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.getEntropy(T)
            self.assertAlmostEqual(s_exp, s_act, delta=1e-4*s_exp)

    def test_get_sum_of_states_classical(self):
        """
        Test the NonlinearRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 1000 * 11.96, 1 * 11.96)
        sum_states = self.mode.getSumOfStates(e_list)
        dens_states = self.mode.getDensityOfStates(e_list)
        for n in range(10, len(e_list)):
            self.assertTrue(0.8 < np.sum(
                dens_states[0:n]) / sum_states[n] < 1.25, '{0} != {1}'.format(np.sum(dens_states[0:n]), sum_states[n]))

    def test_get_sensity_of_states_classical(self):
        """
        Test the NonlinearRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 1000 * 11.96, 1 * 11.96)
        dens_states = self.mode.getDensityOfStates(e_list)
        T = 100
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
        q_exp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(q_exp, q_act, delta=1e-2*q_exp)

    def test_repr(self):
        """
        Test that a NonlinearRotor object can be reconstructed from its
        repr() output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertEqual(self.mode.inertia.value.shape, mode.inertia.value.shape)
        for I0, I in zip(self.mode.inertia.value, mode.inertia.value):
            self.assertAlmostEqual(I0, I, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

    def test_pickle(self):
        """
        Test that a NonlinearRotor object can be pickled and unpickled with
        no loss of information.
        """
        import pickle
        mode = pickle.loads(pickle.dumps(self.mode, -1))
        self.assertEqual(self.mode.inertia.value.shape, mode.inertia.value.shape)
        for I0, I in zip(self.mode.inertia.value, mode.inertia.value):
            self.assertAlmostEqual(I0, I, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

################################################################################


class TestKRotor(unittest.TestCase):
    """
    Contains unit tests of the KRotor class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75
        self.symmetry = 2
        self.quantum = False
        self.mode = KRotor(
            inertia=(self.inertia, "amu*angstrom^2"),
            symmetry=self.symmetry,
            quantum=self.quantum,
        )

    def test_get_rotational_constant(self):
        """
        Test getting the KRotor.rotationalConstant property.
        """
        b_exp = 1.434692
        b_act = self.mode.rotationalConstant.value_si
        self.assertAlmostEqual(b_exp, b_act, 4)

    def test_set_rotational_constant(self):
        """
        Test setting the KRotor.rotationalConstant property.
        """
        rotational_constant = self.mode.rotationalConstant
        rotational_constant.value_si *= 2
        self.mode.rotationalConstant = rotational_constant
        i_exp = 0.5 * self.inertia
        i_act = self.mode.inertia.value_si * constants.Na * 1e23
        self.assertAlmostEqual(i_exp, i_act, 4)

    def test_get_level_energy(self):
        """
        Test the KRotor.getLevelEnergy() method.
        """
        rotational_constant = self.mode.rotationalConstant.value_si * constants.h * constants.c * 100.
        rotational_constant *= constants.Na
        for J in range(0, 100):
            e_exp = float(rotational_constant * J * J)
            e_act = float(self.mode.getLevelEnergy(J))
            if J == 0:
                self.assertEqual(e_act, 0)
            else:
                self.assertAlmostEqual(e_exp, e_act, delta=1e-4*e_exp)

    def test_get_level_degeneracy(self):
        """
        Test the KRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            g_exp = 1 if J == 0 else 2
            g_act = self.mode.getLevelDegeneracy(J)
            self.assertEqual(g_exp, g_act, '{0} != {1}'.format(g_act, g_exp))

    def test_get_partition_function_classical(self):
        """
        Test the KRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([10.6839, 13.7929, 19.5060, 23.8899, 27.5857])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-4*q_exp)

    def test_get_partition_function_quantum(self):
        """
        Test the KRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([10.6839, 13.7929, 19.5060, 23.8899, 27.5857])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-4*q_exp)

    def test_get_heat_capacity_classical(self):
        """
        Test the KRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(cv_exp, cv_act, delta=1e-4*cv_exp)

    def test_get_heat_capacity_quantum(self):
        """
        Test the KRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(cv_exp, cv_act, delta=1e-4*cv_exp)

    def test_get_enthalpy_classical(self):
        """
        Test the KRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R * t_list
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(h_exp, h_act, delta=1e-4*h_exp)

    def test_get_enthalpy_quantum(self):
        """
        Test the KRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R * t_list
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(h_exp, h_act, delta=1e-4*h_exp)

    def test_get_entropy_classical(self):
        """
        Test the KRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([2.86874, 3.12415, 3.47072, 3.67346, 3.81730]) * constants.R
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.getEntropy(T)
            self.assertAlmostEqual(s_exp, s_act, delta=1e-4*s_exp)

    def test_get_entropy_quantum(self):
        """
        Test the KRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([2.86874, 3.12415, 3.47072, 3.67346, 3.81730]) * constants.R
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.getEntropy(T)
            self.assertAlmostEqual(s_exp, s_act, delta=1e-4*s_exp)

    def test_get_sum_of_states_classical(self):
        """
        Test the KRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 1000 * 11.96, 1 * 11.96)
        sum_states = self.mode.getSumOfStates(e_list)
        dens_states = self.mode.getDensityOfStates(e_list)
        for n in range(10, len(e_list)):
            self.assertTrue(0.75 < np.sum(dens_states[0: n + 1]) / sum_states[n] < 1.3333, '{0} != {1}'.format(
                np.sum(dens_states[0: n + 1]), sum_states[n]))

    def test_get_sum_of_states_quantum(self):
        """
        Test the KRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        e_list = np.arange(0, 1000 * 11.96, 1 * 11.96)
        sum_states = self.mode.getSumOfStates(e_list)
        dens_states = self.mode.getDensityOfStates(e_list)
        for n in range(10, len(e_list)):
            self.assertTrue(0.8 < np.sum(dens_states[0: n + 1]) / sum_states[n] < 1.25, '{0} != {1}'.format(
                np.sum(dens_states[0: n + 1]), sum_states[n]))

    def test_get_density_of_states_classical(self):
        """
        Test the KRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 3000 * 11.96, 0.05 * 11.96)
        dens_states = self.mode.getDensityOfStates(e_list)
        T = 500
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
        q_exp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(q_exp, q_act, delta=1e-2*q_exp)

    def test_get_density_of_states_quantum(self):
        """
        Test the KRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        e_list = np.arange(0, 4000 * 11.96, 2 * 11.96)
        dens_states = self.mode.getDensityOfStates(e_list)
        T = 500
        q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
        q_exp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(q_exp, q_act, delta=1e-2*q_exp)

    def test_repr(self):
        """
        Test that a KRotor object can be reconstructed from its repr() output
        with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

    def test_pickle(self):
        """
        Test that a KRotor object can be pickled and unpickled with no loss
        of information.
        """
        import pickle
        mode = pickle.loads(pickle.dumps(self.mode, -1))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

################################################################################


class TestSphericalTopRotor(unittest.TestCase):
    """
    Contains unit tests of the SphericalTopRotor class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75
        self.symmetry = 2
        self.quantum = False
        self.mode = SphericalTopRotor(
            inertia=(self.inertia, "amu*angstrom^2"),
            symmetry=self.symmetry,
            quantum=self.quantum,
        )

    def test_get_rotational_constant(self):
        """
        Test getting the SphericalTopRotor.rotationalConstant property.
        """
        b_exp = 1.434692
        b_act = self.mode.rotationalConstant.value_si
        self.assertAlmostEqual(b_exp, b_act, 4)

    def test_set_rotational_constant(self):
        """
        Test setting the SphericalTopRotor.rotationalConstant property.
        """
        rotational_constant = self.mode.rotationalConstant
        rotational_constant.value_si *= 2
        self.mode.rotationalConstant = rotational_constant
        i_exp = 0.5 * self.inertia
        i_act = self.mode.inertia.value_si * constants.Na * 1e23
        self.assertAlmostEqual(i_exp, i_act, 4)

    def test_get_level_energy(self):
        """
        Test the SphericalTopRotor.getLevelEnergy() method.
        """
        rotational_constant = self.mode.rotationalConstant.value_si * constants.h * constants.c * 100.
        rotational_constant *= constants.Na
        for J in range(0, 100):
            e_exp = rotational_constant * J * (J + 1)
            e_act = self.mode.getLevelEnergy(J)
            if J == 0:
                self.assertEqual(e_act, 0)
            else:
                self.assertAlmostEqual(e_exp, e_act, delta=1e-4*e_exp)

    def test_get_level_degeneracy(self):
        """
        Test the SphericalTopRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            g_exp = (2 * J + 1)**2
            g_act = self.mode.getLevelDegeneracy(J)
            self.assertEqual(g_exp, g_act, '{0} != {1}'.format(g_act, g_exp))

    def test_get_partition_function_classical(self):
        """
        Test the SphericalTopRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([1552.74, 3340.97, 9449.69, 17360.2, 26727.8])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-4*q_exp)

    def test_get_partition_function_quantum(self):
        """
        Test the SphericalTopRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        q_exp_list = np.array([1555.42, 3344.42, 9454.57, 17366.2, 26734.7])
        for T, q_exp in zip(t_list, q_exp_list):
            q_act = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-4*q_exp)

    def test_get_heat_capacity_classical(self):
        """
        Test the SphericalTopRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(cv_exp, cv_act, delta=1e-4*cv_exp)

    def test_get_heat_capacity_quantum(self):
        """
        Test the SphericalTopRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        cv_exp_list = np.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, cv_exp in zip(t_list, cv_exp_list):
            cv_act = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(cv_exp, cv_act, delta=1e-4*cv_exp)

    def test_get_enthalpy_classical(self):
        """
        Test the SphericalTopRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R * t_list
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(h_exp, h_act, delta=1e-4*h_exp)

    def test_get_enthalpy_quantum(self):
        """
        Test the SphericalTopRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        h_exp_list = np.array([1.49828, 1.49897, 1.49948, 1.49966, 1.49974]) * constants.R * t_list
        for T, h_exp in zip(t_list, h_exp_list):
            h_act = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(h_exp, h_act, delta=1e-4*h_exp)

    def test_get_entropy_classical(self):
        """
        Test the SphericalTopRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([8.84778, 9.61402, 10.6537, 11.2619, 11.6935]) * constants.R
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.getEntropy(T)
            self.assertAlmostEqual(s_exp, s_act, delta=1e-4*s_exp)

    def test_get_entropy_quantum(self):
        """
        Test the SphericalTopRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 500, 1000, 1500, 2000])
        s_exp_list = np.array([8.84778, 9.61402, 10.6537, 11.2619, 11.6935]) * constants.R
        for T, s_exp in zip(t_list, s_exp_list):
            s_act = self.mode.getEntropy(T)
            self.assertAlmostEqual(s_exp, s_act, delta=1e-4*s_exp)

    def test_get_sum_of_states_classical(self):
        """
        Test the SphericalTopRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        e_list = np.arange(0, 2000 * 11.96, 1.0 * 11.96)
        dens_states = self.mode.getDensityOfStates(e_list)
        sum_states = self.mode.getSumOfStates(e_list)
        for n in range(20, len(e_list)):
            self.assertAlmostEqual(np.sum(dens_states[0: n + 1]) / sum_states[n], 1.0, 1)

    def test_get_sum_of_states_quantum(self):
        """
        Test the SphericalTopRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        e_list = np.arange(0, 2000 * 11.96, 1.0 * 11.96)
        dens_states = self.mode.getDensityOfStates(e_list)
        sum_states = self.mode.getSumOfStates(e_list)
        for n in range(1, len(e_list)):
            self.assertAlmostEqual(np.sum(dens_states[0: n + 1]) / sum_states[n], 1.0, 3)

    def test_get_density_of_states_classical(self):
        """
        Test the SphericalTopRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        t_list = np.array([300, 400, 500])
        e_list = np.arange(0, 2000 * 11.96, 1.0 * 11.96)
        for T in t_list:
            dens_states = self.mode.getDensityOfStates(e_list)
            q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
            q_exp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-2*q_exp)

    def test_get_density_of_states_quantum(self):
        """
        Test the SphericalTopRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        t_list = np.array([300, 400, 500])
        e_list = np.arange(0, 4000 * 11.96, 2.0 * 11.96)
        for T in t_list:
            dens_states = self.mode.getDensityOfStates(e_list)
            q_act = np.sum(dens_states * np.exp(-e_list / constants.R / T))
            q_exp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(q_exp, q_act, delta=1e-2*q_exp)

    def test_repr(self):
        """
        Test that a SphericalTopRotor object can be reconstructed from its
        repr() output with no loss of information.
        """
        mode = None
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)

    def test_pickle(self):
        """
        Test that a SphericalTopRotor object can be pickled and unpickled
        with no loss of information.
        """
        import pickle
        mode = pickle.loads(pickle.dumps(self.mode, -1))
        self.assertAlmostEqual(self.mode.inertia.value, mode.inertia.value, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
