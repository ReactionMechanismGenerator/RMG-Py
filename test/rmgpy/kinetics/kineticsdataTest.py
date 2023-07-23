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
This script contains unit tests of the :mod:`rmgpy.kinetics.arrhenius` module.
"""


import numpy as np

import rmgpy.constants as constants
from rmgpy.kinetics.kineticsdata import KineticsData, PDepKineticsData


class TestKineticsData:
    """
    Contains unit tests of the :class:`KineticsData` class.
    """

    def setup_class(self):
        self.Tdata = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000], float)
        self.kdata = np.array(
            [
                4.73e-19,
                3.93e-17,
                6.51e-16,
                4.60e-15,
                2.03e-14,
                6.28e-14,
                1.58e-13,
                3.31e-13,
                3.72e-12,
                1.49e-11,
            ],
            float,
        )
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.comment = "H + CH4 <=> H2 + CH3 (RPMD)"
        self.kinetics = KineticsData(
            Tdata=(self.Tdata, "K"),
            kdata=(self.kdata, "cm^3/(molecule*s)"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_temperature_data(self):
        """
        Test that the KineticsData Tdata property was properly set.
        """
        assert self.kinetics.Tdata.value_si.shape == self.Tdata.shape
        for T, T0 in zip(self.kinetics.Tdata.value_si, self.Tdata):
            assert round(abs(T - T0), 4) == 0

    def test_kdata(self):
        """
        Test that the KineticsData kdata property was properly set.
        """
        assert self.kinetics.kdata.value_si.shape == self.kdata.shape
        for k, k0 in zip(self.kinetics.kdata.value_si, self.kdata):
            k0 *= constants.Na * 1e-6
            assert abs(k - k0) < 1e-6 * k0

    def test_temperature_min(self):
        """
        Test that the KineticsData Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the KineticsData Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the KineticsData comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the KineticsData.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the KineticsData.get_rate_coefficient() method.
        """
        Tlist = np.array([300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        kexplist = np.array(
            [
                2.84847e-01,
                2.36670e01,
                2.77019e03,
                3.78191e04,
                1.99333e05,
                5.24644e05,
                1.38086e06,
                2.95680e06,
                5.15086e06,
                8.97299e06,
            ]
        )
        for T, kexp in zip(Tlist, kexplist):
            kact = self.kinetics.get_rate_coefficient(T)
            assert abs(kexp - kact) < 1e-4 * kexp

    def test_pickle(self):
        """
        Test that a KineticsData object can be pickled and unpickled with no
        loss of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        assert self.kinetics.Tdata.value.shape == kinetics.Tdata.value.shape
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            assert round(abs(T - T0), 4) == 0
        assert self.kinetics.Tdata.units == kinetics.Tdata.units
        assert self.kinetics.kdata.value.shape == kinetics.kdata.value.shape
        for k, k0 in zip(self.kinetics.kdata.value, kinetics.kdata.value):
            assert round(abs(k - k0), 4) == 0
        assert self.kinetics.kdata.units == kinetics.kdata.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a KineticsData object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        assert self.kinetics.Tdata.value.shape == kinetics.Tdata.value.shape
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            assert round(abs(T - T0), 4) == 0
        assert self.kinetics.Tdata.units == kinetics.Tdata.units
        assert self.kinetics.kdata.value.shape == kinetics.kdata.value.shape
        for k, k0 in zip(self.kinetics.kdata.value, kinetics.kdata.value):
            assert round(abs(k - k0), 4) == 0
        assert self.kinetics.kdata.units == kinetics.kdata.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment


class TestPDepKineticsData:
    """
    Contains unit tests of the :class:`PDepKineticsData` class.
    """

    def setup_class(self):
        self.Tdata = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000], float)
        self.Pdata = np.array([1e-1, 1e0, 1e1], float)
        self.kdata = np.array(
            [
                [
                    4.73e-21,
                    3.93e-19,
                    6.51e-18,
                    4.60e-17,
                    2.03e-16,
                    6.28e-16,
                    1.58e-15,
                    3.31e-15,
                    3.72e-14,
                    1.49e-13,
                ],
                [
                    4.73e-20,
                    3.93e-18,
                    6.51e-17,
                    4.60e-16,
                    2.03e-15,
                    6.28e-15,
                    1.58e-14,
                    3.31e-14,
                    3.72e-13,
                    1.49e-12,
                ],
                [
                    4.73e-19,
                    3.93e-17,
                    6.51e-16,
                    4.60e-15,
                    2.03e-14,
                    6.28e-14,
                    1.58e-13,
                    3.31e-13,
                    3.72e-12,
                    1.49e-11,
                ],
            ],
            float,
        ).T
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.Pmin = 1e-1
        self.Pmax = 1e1
        self.comment = ""
        self.kinetics = PDepKineticsData(
            Tdata=(self.Tdata, "K"),
            Pdata=(self.Pdata, "bar"),
            kdata=(self.kdata, "cm^3/(molecule*s)"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )

    def test_temperature_data(self):
        """
        Test that the PDepKineticsData Tdata property was properly set.
        """
        assert self.kinetics.Tdata.value_si.shape == self.Tdata.shape
        for T, T0 in zip(self.kinetics.Tdata.value_si, self.Tdata):
            assert round(abs(T - T0), 4) == 0

    def test_pressure_data(self):
        """
        Test that the PDepKineticsData Pdata property was properly set.
        """
        assert self.kinetics.Pdata.value_si.shape == self.Pdata.shape
        for P, P0 in zip(self.kinetics.Pdata.value_si, self.Pdata):
            assert round(abs(P * 1e-5 - P0), 4) == 0

    def test_kdata(self):
        """
        Test that the PDepKineticsData kdata property was properly set.
        """
        assert self.kinetics.kdata.value_si.shape == self.kdata.shape
        for i in range(self.kdata.shape[0]):
            for j in range(self.kdata.shape[1]):
                k0 = self.kdata[i, j] * constants.Na * 1e-6
                k = self.kinetics.kdata.value_si[i, j]
                assert abs(k - k0) < 1e-6 * k0

    def test_temperature_min(self):
        """
        Test that the PDepKineticsData Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the PDepKineticsData Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_pressure_min(self):
        """
        Test that the PDepKineticsData Pmin property was properly set.
        """
        assert round(abs(self.kinetics.Pmin.value_si * 1e-5 - self.Pmin), 6) == 0

    def test_pressure_max(self):
        """
        Test that the PDepKineticsData Pmax property was properly set.
        """
        assert round(abs(self.kinetics.Pmax.value_si * 1e-5 - self.Pmax), 6) == 0

    def test_comment(self):
        """
        Test that the PDepKineticsData comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the PDepKineticsData.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.is_temperature_valid(T)
            assert valid0 == valid

    def test_is_pressure_valid(self):
        """
        Test the PDepKineticsData.is_pressure_valid() method.
        """
        Pdata = np.array([1e3, 1e4, 1e5, 1e6, 1e7])
        validdata = np.array([False, True, True, True, False], bool)
        for P, valid in zip(Pdata, validdata):
            valid0 = self.kinetics.is_pressure_valid(P)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the PDepKineticsData.get_rate_coefficient() method.
        """
        Tlist = np.array([300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        Plist = np.array([1e4, 1e5, 1e6])
        kexplist = np.array(
            [
                [
                    2.84847e-03,
                    2.36670e-01,
                    2.77019e01,
                    3.78191e02,
                    1.99333e03,
                    5.24644e03,
                    1.38086e04,
                    2.95680e04,
                    5.15086e04,
                    8.97299e04,
                ],
                [
                    2.84847e-02,
                    2.36670e00,
                    2.77019e02,
                    3.78191e03,
                    1.99333e04,
                    5.24644e04,
                    1.38086e05,
                    2.95680e05,
                    5.15086e05,
                    8.97299e05,
                ],
                [
                    2.84847e-01,
                    2.36670e01,
                    2.77019e03,
                    3.78191e04,
                    1.99333e05,
                    5.24644e05,
                    1.38086e06,
                    2.95680e06,
                    5.15086e06,
                    8.97299e06,
                ],
            ]
        ).T
        for i in range(Tlist.shape[0]):
            for j in range(Plist.shape[0]):
                kexp = kexplist[i, j]
                kact = self.kinetics.get_rate_coefficient(Tlist[i], Plist[j])
                assert abs(kexp - kact) < 1e-4 * kexp

    def test_pickle(self):
        """
        Test that a PDepKineticsData object can be pickled and unpickled with no
        loss of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        assert self.kinetics.Tdata.value.shape == kinetics.Tdata.value.shape
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            assert round(abs(T - T0), 4) == 0
        assert self.kinetics.Tdata.units == kinetics.Tdata.units
        assert self.kinetics.Pdata.value.shape == kinetics.Pdata.value.shape
        for P, P0 in zip(self.kinetics.Pdata.value, kinetics.Pdata.value):
            assert round(abs(P - P0), 4) == 0
        assert self.kinetics.Pdata.units == kinetics.Pdata.units
        assert self.kinetics.kdata.value.shape == kinetics.kdata.value.shape
        for i in range(self.kinetics.kdata.value.shape[0]):
            for j in range(self.kinetics.kdata.value.shape[1]):
                k0 = self.kinetics.kdata.value[i, j]
                k = kinetics.kdata.value[i, j]
                assert abs(k - k0) < 1e-6 * k0
        assert self.kinetics.kdata.units == kinetics.kdata.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a PDepKineticsData object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        assert self.kinetics.Tdata.value.shape == kinetics.Tdata.value.shape
        for T, T0 in zip(self.kinetics.Tdata.value, kinetics.Tdata.value):
            assert round(abs(T - T0), 4) == 0
        assert self.kinetics.Tdata.units == kinetics.Tdata.units
        assert self.kinetics.Pdata.value.shape == kinetics.Pdata.value.shape
        for P, P0 in zip(self.kinetics.Pdata.value, kinetics.Pdata.value):
            assert round(abs(P - P0), 4) == 0
        assert self.kinetics.Pdata.units == kinetics.Pdata.units
        assert self.kinetics.kdata.value.shape == kinetics.kdata.value.shape
        for i in range(self.kinetics.kdata.value.shape[0]):
            for j in range(self.kinetics.kdata.value.shape[1]):
                k0 = self.kinetics.kdata.value[i, j]
                k = kinetics.kdata.value[i, j]
                assert abs(k - k0) < 1e-6 * k0
        assert self.kinetics.kdata.units == kinetics.kdata.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment
