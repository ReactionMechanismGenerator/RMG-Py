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
This script contains unit tests of the :mod:`rmgpy.kinetics.surface` module.
"""

import unittest

import numpy as np

from rmgpy.kinetics.surface import StickingCoefficient, SurfaceArrhenius, SurfaceChargeTransfer
from rmgpy.species import Species
from rmgpy.molecule import Molecule
import rmgpy.quantity as quantity
import rmgpy.constants as constants

################################################################################


class TestStickingCoefficient(unittest.TestCase):
    """
    Contains unit tests of the :class:`StickingCoefficient` class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.A = 4.3e-2
        self.n = -0.21
        self.Ea = 1.2
        self.T0 = 1.
        self.Tmin = 300.
        self.Tmax = 3000.
        s = Species().from_adjacency_list('1 X u0 p0 c0')
        s.label = 'X'
        self.coverage_dependence = {s: {'a': quantity.Dimensionless(0.0),
                                        'm': quantity.Dimensionless(-1.0),
                                        'E': quantity.Energy(0.0, 'J/mol')}}
        self.comment = 'O2 dissociative'
        self.stick = StickingCoefficient(
            A=self.A,
            n=self.n,
            Ea=(self.Ea, "kJ/mol"),
            T0=(self.T0, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
            coverage_dependence=self.coverage_dependence,
        )

    def test_A(self):
        """
        Test that the StickingCoefficient A property was properly set.
        """
        self.assertAlmostEqual(self.stick.A.value_si, self.A, delta=1e0)

    def test_n(self):
        """
        Test that the StickingCoefficient n property was properly set.
        """
        self.assertAlmostEqual(self.stick.n.value_si, self.n, 6)

    def test_Ea(self):
        """
        Test that the StickingCoefficient Ea property was properly set.
        """
        self.assertAlmostEqual(self.stick.Ea.value_si * 0.001, self.Ea, 6)

    def test_T0(self):
        """
        Test that the StickingCoefficient T0 property was properly set.
        """
        self.assertAlmostEqual(self.stick.T0.value_si, self.T0, 6)

    def test_Tmin(self):
        """
        Test that the StickingCoefficient Tmin property was properly set.
        """
        self.assertAlmostEqual(self.stick.Tmin.value_si, self.Tmin, 6)

    def test_Tmax(self):
        """
        Test that the StickingCoefficient Tmax property was properly set.
        """
        self.assertAlmostEqual(self.stick.Tmax.value_si, self.Tmax, 6)

    def test_comment(self):
        """
        Test that the StickingCoefficient comment property was properly set.
        """
        self.assertEqual(self.stick.comment, self.comment)

    def test_coverage_dependence(self):
        """
        Test that the coverage dependent parameters was properly set.
        """
        for key in self.stick.coverage_dependence.keys():
            match = False
            for key2 in self.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            self.assertTrue(match)
        for species, parameters in self.stick.coverage_dependence.items():
            match = False
            for species2 in self.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        self.assertEqual(value.value_si, self.coverage_dependence[species2][key].value_si)
            self.assertTrue(match)

    def test_is_temperature_valid(self):
        """
        Test the StickingCoefficient.is_temperature_valid() method.
        """
        T_data = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 4000])
        valid_data = np.array([False, True, True, True, True, True, True, True, True, False], np.bool)
        for T, valid in zip(T_data, valid_data):
            valid0 = self.stick.is_temperature_valid(T)
            self.assertEqual(valid0, valid)

    def test_pickle(self):
        """
        Test that an StickingCoefficient object can be pickled and unpickled with no loss
        of information.
        """
        import pickle
        stick = pickle.loads(pickle.dumps(self.stick, -1))
        self.assertAlmostEqual(self.stick.A.value, stick.A.value, delta=1e0)
        self.assertEqual(self.stick.A.units, stick.A.units)
        self.assertAlmostEqual(self.stick.n.value, stick.n.value, 4)
        self.assertAlmostEqual(self.stick.Ea.value, stick.Ea.value, 4)
        self.assertEqual(self.stick.Ea.units, stick.Ea.units)
        self.assertAlmostEqual(self.stick.T0.value, stick.T0.value, 4)
        self.assertEqual(self.stick.T0.units, stick.T0.units)
        self.assertAlmostEqual(self.stick.Tmin.value, stick.Tmin.value, 4)
        self.assertEqual(self.stick.Tmin.units, stick.Tmin.units)
        self.assertAlmostEqual(self.stick.Tmax.value, stick.Tmax.value, 4)
        self.assertEqual(self.stick.Tmax.units, stick.Tmax.units)
        self.assertEqual(self.stick.comment, stick.comment)
        for key in self.stick.coverage_dependence.keys():
            match = False
            for key2 in stick.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            self.assertTrue(match)
        for species, parameters in self.stick.coverage_dependence.items():
            match = False
            for species2 in stick.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        self.assertEqual(value.value_si, stick.coverage_dependence[species2][key].value_si)
            self.assertTrue(match)
        self.assertEqual(dir(self.stick), dir(stick))

    def test_repr(self):
        """
        Test that an StickingCoefficient object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec(f'stick = {self.stick!r}', globals(), namespace)
        self.assertIn('stick', namespace)
        stick = namespace['stick']
        self.assertAlmostEqual(self.stick.A.value, stick.A.value, delta=1e0)
        self.assertEqual(self.stick.A.units, stick.A.units)
        self.assertAlmostEqual(self.stick.n.value, stick.n.value, 4)
        self.assertAlmostEqual(self.stick.Ea.value, stick.Ea.value, 4)
        self.assertEqual(self.stick.Ea.units, stick.Ea.units)
        self.assertAlmostEqual(self.stick.T0.value, stick.T0.value, 4)
        self.assertEqual(self.stick.T0.units, stick.T0.units)
        self.assertAlmostEqual(self.stick.Tmin.value, stick.Tmin.value, 4)
        self.assertEqual(self.stick.Tmin.units, stick.Tmin.units)
        self.assertAlmostEqual(self.stick.Tmax.value, stick.Tmax.value, 4)
        self.assertEqual(self.stick.Tmax.units, stick.Tmax.units)
        self.assertEqual(self.stick.comment, stick.comment)
        self.assertEqual([m.label for m in self.stick.coverage_dependence.keys()], list(stick.coverage_dependence.keys()))
        for species, parameters in self.stick.coverage_dependence.items():
            for key, value in parameters.items():
                self.assertEqual(value.value_si, stick.coverage_dependence[species.label][key].value_si)
        self.assertEqual(dir(self.stick), dir(stick))

    def test_copy(self):
        """
        Test that an StickingCoefficient object can be copied with deepcopy
        with no loss of information.
        """
        import copy
        stick = copy.deepcopy(self.stick)
        self.assertAlmostEqual(self.stick.A.value, stick.A.value, delta=1e0)
        self.assertEqual(self.stick.A.units, stick.A.units)
        self.assertAlmostEqual(self.stick.n.value, stick.n.value, 4)
        self.assertAlmostEqual(self.stick.Ea.value, stick.Ea.value, 4)
        self.assertEqual(self.stick.Ea.units, stick.Ea.units)
        self.assertAlmostEqual(self.stick.T0.value, stick.T0.value, 4)
        self.assertEqual(self.stick.T0.units, stick.T0.units)
        self.assertAlmostEqual(self.stick.Tmin.value, stick.Tmin.value, 4)
        self.assertEqual(self.stick.Tmin.units, stick.Tmin.units)
        self.assertAlmostEqual(self.stick.Tmax.value, stick.Tmax.value, 4)
        self.assertEqual(self.stick.Tmax.units, stick.Tmax.units)
        self.assertEqual(self.stick.comment, stick.comment)
        for key in self.stick.coverage_dependence.keys():
            match = False
            for key2 in stick.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            self.assertTrue(match)
        for species, parameters in self.stick.coverage_dependence.items():
            match = False
            for species2 in stick.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        self.assertEqual(value.value_si, stick.coverage_dependence[species2][key].value_si)
            self.assertTrue(match)
        self.assertEqual(dir(self.stick), dir(stick))

    def test_is_identical_to(self):
        """
        Test that the StickingCoefficient.is_identical_to method works on itself
        """
        self.assertTrue(self.stick.is_identical_to(self.stick))

################################################################################


class TestSurfaceArrhenius(unittest.TestCase):
    """
    Contains unit tests of the :class:`SurfaceArrhenius` class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.A = 1.44e18
        self.n = -0.087
        self.Ea = 63.4
        self.T0 = 1.
        self.Tmin = 300.
        self.Tmax = 3000.
        s = Species().from_adjacency_list('1 X u0 p0 c0')
        s.label = 'X'
        self.coverage_dependence = {s: {'a': quantity.Dimensionless(0.0),
                                        'm': quantity.Dimensionless(-1.0),
                                        'E': quantity.Energy(0.0, 'J/mol')}}
        self.comment = 'CH3x + Hx <=> CH4 + x + x'
        self.surfarr = SurfaceArrhenius(
            A=(self.A, 'm^2/(mol*s)'),
            n=self.n,
            Ea=(self.Ea, "kJ/mol"),
            T0=(self.T0, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
            coverage_dependence=self.coverage_dependence,
        )

    def test_A(self):
        """
        Test that the SurfaceArrhenius A property was properly set.
        """
        self.assertAlmostEqual(self.surfarr.A.value_si, self.A, delta=1e0)

    def test_n(self):
        """
        Test that the SurfaceArrhenius n property was properly set.
        """
        self.assertAlmostEqual(self.surfarr.n.value_si, self.n, 6)

    def test_Ea(self):
        """
        Test that the SurfaceArrhenius Ea property was properly set.
        """
        self.assertAlmostEqual(self.surfarr.Ea.value_si * 0.001, self.Ea, 6)

    def test_T0(self):
        """
        Test that the SurfaceArrhenius T0 property was properly set.
        """
        self.assertAlmostEqual(self.surfarr.T0.value_si, self.T0, 6)

    def test_Tmin(self):
        """
        Test that the SurfaceArrhenius Tmin property was properly set.
        """
        self.assertAlmostEqual(self.surfarr.Tmin.value_si, self.Tmin, 6)

    def test_Tmax(self):
        """
        Test that the SurfaceArrhenius Tmax property was properly set.
        """
        self.assertAlmostEqual(self.surfarr.Tmax.value_si, self.Tmax, 6)

    def test_comment(self):
        """
        Test that the SurfaceArrhenius comment property was properly set.
        """
        self.assertEqual(self.surfarr.comment, self.comment)

    def test_coverage_dependence(self):
        """
        Test that the coverage dependent parameters was properly set.
        """
        for key in self.surfarr.coverage_dependence.keys():
            match = False
            for key2 in self.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            self.assertTrue(match)
        for species, parameters in self.surfarr.coverage_dependence.items():
            match = False
            for species2 in self.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        self.assertEqual(value.value_si, self.coverage_dependence[species2][key].value_si)
            self.assertTrue(match)

    def test_is_temperature_valid(self):
        """
        Test the SurfaceArrhenius.is_temperature_valid() method.
        """
        T_data = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 4000])
        valid_data = np.array([False, True, True, True, True, True, True, True, True, False], np.bool)
        for T, valid in zip(T_data, valid_data):
            valid0 = self.surfarr.is_temperature_valid(T)
            self.assertEqual(valid0, valid)

    def test_pickle(self):
        """
        Test that an SurfaceArrhenius object can be pickled and unpickled with no loss
        of information.
        """
        import pickle
        surfarr = pickle.loads(pickle.dumps(self.surfarr, -1))
        self.assertAlmostEqual(self.surfarr.A.value, surfarr.A.value, delta=1e0)
        self.assertEqual(self.surfarr.A.units, surfarr.A.units)
        self.assertAlmostEqual(self.surfarr.n.value, surfarr.n.value, 4)
        self.assertAlmostEqual(self.surfarr.Ea.value, surfarr.Ea.value, 4)
        self.assertEqual(self.surfarr.Ea.units, surfarr.Ea.units)
        self.assertAlmostEqual(self.surfarr.T0.value, surfarr.T0.value, 4)
        self.assertEqual(self.surfarr.T0.units, surfarr.T0.units)
        self.assertAlmostEqual(self.surfarr.Tmin.value, surfarr.Tmin.value, 4)
        self.assertEqual(self.surfarr.Tmin.units, surfarr.Tmin.units)
        self.assertAlmostEqual(self.surfarr.Tmax.value, surfarr.Tmax.value, 4)
        self.assertEqual(self.surfarr.Tmax.units, surfarr.Tmax.units)
        self.assertEqual(self.surfarr.comment, surfarr.comment)
        for key in self.surfarr.coverage_dependence.keys():
            match = False
            for key2 in surfarr.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            self.assertTrue(match)
        for species, parameters in self.surfarr.coverage_dependence.items():
            match = False
            for species2 in surfarr.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        self.assertEqual(value.value_si, surfarr.coverage_dependence[species2][key].value_si)
            self.assertTrue(match)
        self.assertEqual(dir(self.surfarr), dir(surfarr))

    def test_repr(self):
        """
        Test that an SurfaceArrhenius object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec('surfarr = {0!r}'.format(self.surfarr), globals(), namespace)
        self.assertIn('surfarr', namespace)
        surfarr = namespace['surfarr']
        self.assertAlmostEqual(self.surfarr.A.value, surfarr.A.value, delta=1e0)
        self.assertEqual(self.surfarr.A.units, surfarr.A.units)
        self.assertAlmostEqual(self.surfarr.n.value, surfarr.n.value, 4)
        self.assertAlmostEqual(self.surfarr.Ea.value, surfarr.Ea.value, 4)
        self.assertEqual(self.surfarr.Ea.units, surfarr.Ea.units)
        self.assertAlmostEqual(self.surfarr.T0.value, surfarr.T0.value, 4)
        self.assertEqual(self.surfarr.T0.units, surfarr.T0.units)
        self.assertAlmostEqual(self.surfarr.Tmin.value, surfarr.Tmin.value, 4)
        self.assertEqual(self.surfarr.Tmin.units, surfarr.Tmin.units)
        self.assertAlmostEqual(self.surfarr.Tmax.value, surfarr.Tmax.value, 4)
        self.assertEqual(self.surfarr.Tmax.units, surfarr.Tmax.units)
        self.assertEqual(self.surfarr.comment, surfarr.comment)
        self.assertEqual([m.label for m in self.surfarr.coverage_dependence.keys()], list(surfarr.coverage_dependence.keys()))
        for species, parameters in self.surfarr.coverage_dependence.items():
            for key, value in parameters.items():
                self.assertEqual(value.value_si, surfarr.coverage_dependence[species.label][key].value_si)
        self.assertEqual(dir(self.surfarr), dir(surfarr))

    def test_copy(self):
        """
        Test that an SurfaceArrhenius object can be copied with deepcopy
        with no loss of information.
        """
        import copy
        surfarr = copy.deepcopy(self.surfarr)
        self.assertAlmostEqual(self.surfarr.A.value, surfarr.A.value, delta=1e0)
        self.assertEqual(self.surfarr.A.units, surfarr.A.units)
        self.assertAlmostEqual(self.surfarr.n.value, surfarr.n.value, 4)
        self.assertAlmostEqual(self.surfarr.Ea.value, surfarr.Ea.value, 4)
        self.assertEqual(self.surfarr.Ea.units, surfarr.Ea.units)
        self.assertAlmostEqual(self.surfarr.T0.value, surfarr.T0.value, 4)
        self.assertEqual(self.surfarr.T0.units, surfarr.T0.units)
        self.assertAlmostEqual(self.surfarr.Tmin.value, surfarr.Tmin.value, 4)
        self.assertEqual(self.surfarr.Tmin.units, surfarr.Tmin.units)
        self.assertAlmostEqual(self.surfarr.Tmax.value, surfarr.Tmax.value, 4)
        self.assertEqual(self.surfarr.Tmax.units, surfarr.Tmax.units)
        self.assertEqual(self.surfarr.comment, surfarr.comment)
        for key in self.surfarr.coverage_dependence.keys():
            match = False
            for key2 in surfarr.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            self.assertTrue(match)
        for species, parameters in self.surfarr.coverage_dependence.items():
            match = False
            for species2 in surfarr.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        self.assertEqual(value.value_si, surfarr.coverage_dependence[species2][key].value_si)
            self.assertTrue(match)
        self.assertEqual(dir(self.surfarr), dir(surfarr))

    def test_is_identical_to(self):
        """
        Test that the SurfaceArrhenius.is_identical_to method works on itself
        """
        self.assertTrue(self.surfarr.is_identical_to(self.surfarr))

    def test_to_surface_charge_transfer(self):
        """
        Test that the SurfaceArrhenius.to_surface_charge_transfer method works
        """

        surface_charge_transfer = self.surfarr.to_surface_charge_transfer(2,-2)
        self.assertIsInstance(surface_charge_transfer,SurfaceChargeTransfer)
        surface_charge_transfer0 = SurfaceChargeTransfer(
            A = self.surfarr.A,
            n = self.surfarr.n,
            Ea = self.surfarr.Ea,
            T0 = self.surfarr.T0,
            Tmin = self.surfarr.Tmin,
            Tmax = self.surfarr.Tmax,
            electrons = -2,
            V0 = (2,'V')
        )
        self.assertTrue(surface_charge_transfer.is_identical_to(surface_charge_transfer0))


################################################################################


class TestSurfaceChargeTransfer((unittest.TestCase)):
    """
    Contains unit tests of the :class:`SurfaceChargeTransfer` class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.A = 1.44e18
        self.n = -0.087
        self.Ea = 63.4
        self.T0 = 1.
        self.Tmin = 300.
        self.Tmax = 3000.
        self.V0 = 1
        self.electrons = -1
        self.comment = 'CH3x + Hx <=> CH4 + x + x'
        self.surfchargerxn_reduction = SurfaceChargeTransfer(
            A=(self.A, 'm^2/(mol*s)'),
            n=self.n,
            electrons=self.electrons,
            V0=(self.V0, "V"),
            Ea=(self.Ea, "kJ/mol"),
            T0=(self.T0, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

        self.surfchargerxn_oxidation = SurfaceChargeTransfer(
            A=(self.A, 'm^2/(mol*s)'),
            n=self.n,
            electrons=1,
            V0=(self.V0, "V"),
            Ea=(self.Ea, "kJ/mol"),
            T0=(self.T0, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_A(self):
        """
        Test that the SurfaceChargeTransfer A property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.A.value_si, self.A, delta=1e0)

    def test_n(self):
        """
        Test that the SurfaceChargeTransfer n property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.n.value_si, self.n, 6)

    def test_ne(self):
        """
        Test that the SurfaceChargeTransfer electrons property was properly set.
        """
        self.assertEqual(self.surfchargerxn_reduction.electrons.value_si, -1.0)

    def test_alpha(self):
        """
        Test that the SurfaceChargeTransfer alpha property was properly set.
        """
        self.assertEqual(self.surfchargerxn_reduction.alpha.value_si, 0.5)

    def test_Ea(self):
        """
        Test that the SurfaceChargeTransfer Ea property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.Ea.value_si * 0.001, self.Ea, 6)

    def test_T0(self):
        """
        Test that the SurfaceChargeTransfer T0 property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.T0.value_si, self.T0, 6)

    def test_Tmin(self):
        """
        Test that the SurfaceChargeTransfer Tmin property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmin.value_si, self.Tmin, 6)

    def test_Tmax(self):
        """
        Test that the SurfaceChargeTransfer Tmax property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmax.value_si, self.Tmax, 6)

    def test_V0(self):
        """
        Test that the SurfaceChargeTransfer V0 property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.V0.value_si, self.V0, 1)

    def test_Tmax(self):
        """
        Test that the SurfaceChargeTransfer Tmax property was properly set.
        """
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmax.value_si, self.Tmax, 6)

    def test_comment(self):
        """
        Test that the SurfaceChargeTransfer comment property was properly set.
        """
        self.assertEqual(self.surfchargerxn_reduction.comment, self.comment)

    def test_is_temperature_valid(self):
        """
        Test the SurfaceChargeTransfer.is_temperature_valid() method.
        """
        T_data = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 4000])
        valid_data = np.array([False, True, True, True, True, True, True, True, True, False], np.bool)
        for T, valid in zip(T_data, valid_data):
            valid0 = self.surfchargerxn_reduction.is_temperature_valid(T)
            self.assertEqual(valid0, valid)

    def test_pickle(self):
        """
        Test that an SurfaceChargeTransfer object can be pickled and unpickled with no loss
        of information.
        """
        import pickle
        surfchargerxn_reduction = pickle.loads(pickle.dumps(self.surfchargerxn_reduction, -1))
        self.assertAlmostEqual(self.surfchargerxn_reduction.A.value, surfchargerxn_reduction.A.value, delta=1e0)
        self.assertEqual(self.surfchargerxn_reduction.A.units, surfchargerxn_reduction.A.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.n.value, surfchargerxn_reduction.n.value, 4)
        self.assertAlmostEqual(self.surfchargerxn_reduction.electrons.value, surfchargerxn_reduction.electrons.value, 4)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Ea.value, surfchargerxn_reduction.Ea.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Ea.units, surfchargerxn_reduction.Ea.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.T0.value, surfchargerxn_reduction.T0.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.T0.units, surfchargerxn_reduction.T0.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.V0.value, surfchargerxn_reduction.V0.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.V0.units, surfchargerxn_reduction.V0.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmin.value, surfchargerxn_reduction.Tmin.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Tmin.units, surfchargerxn_reduction.Tmin.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmax.value, surfchargerxn_reduction.Tmax.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Tmax.units, surfchargerxn_reduction.Tmax.units)
        self.assertEqual(self.surfchargerxn_reduction.comment, surfchargerxn_reduction.comment)
        self.assertEqual(dir(self.surfchargerxn_reduction), dir(surfchargerxn_reduction))

    def test_repr(self):
        """
        Test that an SurfaceChargeTransfer object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec('surfchargerxn_reduction = {0!r}'.format(self.surfchargerxn_reduction), globals(), namespace)
        self.assertIn('surfchargerxn_reduction', namespace)
        surfchargerxn_reduction = namespace['surfchargerxn_reduction']
        self.assertAlmostEqual(self.surfchargerxn_reduction.A.value, surfchargerxn_reduction.A.value, delta=1e0)
        self.assertEqual(self.surfchargerxn_reduction.A.units, surfchargerxn_reduction.A.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.n.value, surfchargerxn_reduction.n.value, 4)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Ea.value, surfchargerxn_reduction.Ea.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Ea.units, surfchargerxn_reduction.Ea.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.T0.value, surfchargerxn_reduction.T0.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.T0.units, surfchargerxn_reduction.T0.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmin.value, surfchargerxn_reduction.Tmin.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Tmin.units, surfchargerxn_reduction.Tmin.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmax.value, surfchargerxn_reduction.Tmax.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Tmax.units, surfchargerxn_reduction.Tmax.units)
        self.assertEqual(self.surfchargerxn_reduction.comment, surfchargerxn_reduction.comment)
        self.assertEqual(dir(self.surfchargerxn_reduction), dir(surfchargerxn_reduction))

    def test_copy(self):
        """
        Test that an SurfaceChargeTransfer object can be copied with deepcopy
        with no loss of information.
        """
        import copy
        surfchargerxn_reduction = copy.deepcopy(self.surfchargerxn_reduction)
        self.assertAlmostEqual(self.surfchargerxn_reduction.A.value, surfchargerxn_reduction.A.value, delta=1e0)
        self.assertEqual(self.surfchargerxn_reduction.A.units, surfchargerxn_reduction.A.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.n.value, surfchargerxn_reduction.n.value, 4)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Ea.value, surfchargerxn_reduction.Ea.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Ea.units, surfchargerxn_reduction.Ea.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.T0.value, surfchargerxn_reduction.T0.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.T0.units, surfchargerxn_reduction.T0.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmin.value, surfchargerxn_reduction.Tmin.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Tmin.units, surfchargerxn_reduction.Tmin.units)
        self.assertAlmostEqual(self.surfchargerxn_reduction.Tmax.value, surfchargerxn_reduction.Tmax.value, 4)
        self.assertEqual(self.surfchargerxn_reduction.Tmax.units, surfchargerxn_reduction.Tmax.units)
        self.assertEqual(self.surfchargerxn_reduction.comment, surfchargerxn_reduction.comment)
        self.assertEqual(dir(self.surfchargerxn_reduction), dir(surfchargerxn_reduction))

    def test_is_identical_to(self):
        """
        Test that the SurfaceChargeTransfer.is_identical_to method works on itself
        """
        self.assertTrue(self.surfchargerxn_reduction.is_identical_to(self.surfchargerxn_reduction))

    def test_to_surface_arrhenius(self):
        """
        Test that the SurfaceChargeTransfer.to_surface_arrhenius method works
        """
        surface_arr = self.surfchargerxn_reduction.to_surface_arrhenius()
        self.assertIsInstance(surface_arr,SurfaceArrhenius)
        surface_arrhenius0 = SurfaceArrhenius(
            A = self.surfchargerxn_reduction.A,
            n = self.surfchargerxn_reduction.n,
            Ea = self.surfchargerxn_reduction.Ea,
            T0 = self.surfchargerxn_reduction.T0,
            Tmin = self.surfchargerxn_reduction.Tmin,
            Tmax = self.surfchargerxn_reduction.Tmax,
        )

        self.assertTrue(surface_arr.is_identical_to(surface_arrhenius0))

    def test_get_activation_energy_from_potential(self):
        """
        Test that the SurfaceChargeTransfer.get_activation_energy_from_potential method works
        """
 
        electrons_ox = self.surfchargerxn_oxidation.electrons.value_si
        V0_ox = self.surfchargerxn_oxidation.V0.value_si
        Ea0_ox = self.surfchargerxn_oxidation.Ea.value_si
        alpha_ox = self.surfchargerxn_oxidation.alpha.value_si

        electrons_red = self.surfchargerxn_reduction.electrons.value_si
        V0_red = self.surfchargerxn_reduction.V0.value_si
        Ea0_red = self.surfchargerxn_reduction.Ea.value_si
        alpha_red = self.surfchargerxn_reduction.alpha.value_si

        Potentials = (V0_ox + 1, V0_ox, V0_ox - 1)        
        
        for V in Potentials:
            Ea = self.surfchargerxn_oxidation.get_activation_energy_from_potential(V, False)
            self.assertAlmostEqual(Ea0_ox - (alpha_ox * electrons_ox * constants.F * (V-V0_ox)), Ea, 6)
            Ea = self.surfchargerxn_oxidation.get_activation_energy_from_potential(V, True)
            self.assertTrue(Ea>=0)
            Ea = self.surfchargerxn_reduction.get_activation_energy_from_potential(V, False)
            self.assertAlmostEqual(Ea0_red - (alpha_red * electrons_red * constants.F * (V-V0_red)), Ea, 6)

    def test_get_rate_coefficient(self):
        """
        Test that the SurfaceChargeTransfer.to_surface_arrhenius method works
        """

        A_ox = self.surfchargerxn_oxidation.A.value_si
        A_red = self.surfchargerxn_reduction.A.value_si
        electrons_ox = self.surfchargerxn_oxidation.electrons.value_si
        electrons_red = self.surfchargerxn_reduction.electrons.value_si
        n_ox = self.surfchargerxn_oxidation.n.value_si
        n_red = self.surfchargerxn_reduction.n.value_si
        Ea_ox = self.surfchargerxn_oxidation.Ea.value_si
        Ea_red = self.surfchargerxn_reduction.Ea.value_si
        V0_ox = self.surfchargerxn_oxidation.V0.value_si
        V0_red = self.surfchargerxn_reduction.V0.value_si
        T0_ox = self.surfchargerxn_oxidation.T0.value_si
        T0_red = self.surfchargerxn_reduction.T0.value_si
        alpha_ox = self.surfchargerxn_oxidation.alpha.value
        alpha_red = self.surfchargerxn_reduction.alpha.value

        Potentials = (V0_ox + 1, V0_ox, V0_ox - 1) 
        for V in Potentials:
            for T in np.linspace(300,3000,10):
                Ea = Ea_ox - (alpha_ox * electrons_ox * constants.F * (V-V0_ox))
                k_oxidation = A_ox * (T / T0_ox) ** n_ox * np.exp(-Ea / (constants.R * T)) 
                Ea = Ea_red - (alpha_red * electrons_red * constants.F * (V-V0_red))
                k_reduction = A_red * (T / T0_red) ** n_red * np.exp(-Ea / (constants.R * T)) 
                self.assertAlmostEqual(k_oxidation,self.surfchargerxn_oxidation.get_rate_coefficient(T,V))
                self.assertAlmostEqual(k_reduction,self.surfchargerxn_reduction.get_rate_coefficient(T,V))

    def test_change_v0(self):
        
        V0 = self.surfchargerxn_oxidation.V0.value_si
        electrons = self.surfchargerxn_oxidation.electrons.value
        alpha = self.surfchargerxn_oxidation.alpha.value
        for V in (V0 + 1, V0, V0 - 1, V0):
            delta = V - self.surfchargerxn_oxidation.V0.value_si
            V_i = self.surfchargerxn_oxidation.V0.value_si
            Ea_i = self.surfchargerxn_oxidation.Ea.value_si
            self.surfchargerxn_oxidation.change_v0(V)
            self.assertEqual(self.surfchargerxn_oxidation.V0.value_si, V_i + delta)
            self.assertAlmostEqual(self.surfchargerxn_oxidation.Ea.value_si, Ea_i - (alpha *electrons * constants.F * delta), 6)

        V0 = self.surfchargerxn_reduction.V0.value_si
        electrons = self.surfchargerxn_reduction.electrons.value
        alpha = self.surfchargerxn_reduction.alpha.value
        for V in (V0 + 1, V0, V0 - 1, V0):
            delta = V - self.surfchargerxn_reduction.V0.value_si
            V_i = self.surfchargerxn_reduction.V0.value_si
            Ea_i = self.surfchargerxn_reduction.Ea.value_si
            self.surfchargerxn_reduction.change_v0(V)
            self.assertEqual(self.surfchargerxn_reduction.V0.value_si, V_i + delta)
            self.assertAlmostEqual(self.surfchargerxn_reduction.Ea.value_si, Ea_i - (alpha *electrons * constants.F * delta), 6)

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))