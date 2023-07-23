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


import numpy as np

from rmgpy.kinetics.surface import StickingCoefficient, SurfaceArrhenius
from rmgpy.species import Species
from rmgpy.molecule import Molecule
import rmgpy.quantity as quantity


class TestStickingCoefficient:
    """
    Contains unit tests of the :class:`StickingCoefficient` class.
    """

    def setup_class(self):
        self.A = 4.3e-2
        self.n = -0.21
        self.Ea = 1.2
        self.T0 = 1.0
        self.Tmin = 300.0
        self.Tmax = 3000.0
        s = Species().from_adjacency_list("1 X u0 p0 c0")
        s.label = "X"
        self.coverage_dependence = {
            s: {
                "a": quantity.Dimensionless(0.0),
                "m": quantity.Dimensionless(-1.0),
                "E": quantity.Energy(0.0, "J/mol"),
            }
        }
        self.comment = "O2 dissociative"
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
        assert abs(self.stick.A.value_si - self.A) < 1e0

    def test_n(self):
        """
        Test that the StickingCoefficient n property was properly set.
        """
        assert round(abs(self.stick.n.value_si - self.n), 6) == 0

    def test_Ea(self):
        """
        Test that the StickingCoefficient Ea property was properly set.
        """
        assert round(abs(self.stick.Ea.value_si * 0.001 - self.Ea), 6) == 0

    def test_T0(self):
        """
        Test that the StickingCoefficient T0 property was properly set.
        """
        assert round(abs(self.stick.T0.value_si - self.T0), 6) == 0

    def test_Tmin(self):
        """
        Test that the StickingCoefficient Tmin property was properly set.
        """
        assert round(abs(self.stick.Tmin.value_si - self.Tmin), 6) == 0

    def test_Tmax(self):
        """
        Test that the StickingCoefficient Tmax property was properly set.
        """
        assert round(abs(self.stick.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the StickingCoefficient comment property was properly set.
        """
        assert self.stick.comment == self.comment

    def test_coverage_dependence(self):
        """
        Test that the coverage dependent parameters was properly set.
        """
        for key in self.stick.coverage_dependence.keys():
            match = False
            for key2 in self.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            assert match
        for species, parameters in self.stick.coverage_dependence.items():
            match = False
            for species2 in self.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        assert value.value_si == self.coverage_dependence[species2][key].value_si
            assert match

    def test_is_temperature_valid(self):
        """
        Test the StickingCoefficient.is_temperature_valid() method.
        """
        T_data = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 4000])
        valid_data = np.array([False, True, True, True, True, True, True, True, True, False], bool)
        for T, valid in zip(T_data, valid_data):
            valid0 = self.stick.is_temperature_valid(T)
            assert valid0 == valid

    def test_pickle(self):
        """
        Test that an StickingCoefficient object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        stick = pickle.loads(pickle.dumps(self.stick, -1))
        assert abs(self.stick.A.value - stick.A.value) < 1e0
        assert self.stick.A.units == stick.A.units
        assert round(abs(self.stick.n.value - stick.n.value), 4) == 0
        assert round(abs(self.stick.Ea.value - stick.Ea.value), 4) == 0
        assert self.stick.Ea.units == stick.Ea.units
        assert round(abs(self.stick.T0.value - stick.T0.value), 4) == 0
        assert self.stick.T0.units == stick.T0.units
        assert round(abs(self.stick.Tmin.value - stick.Tmin.value), 4) == 0
        assert self.stick.Tmin.units == stick.Tmin.units
        assert round(abs(self.stick.Tmax.value - stick.Tmax.value), 4) == 0
        assert self.stick.Tmax.units == stick.Tmax.units
        assert self.stick.comment == stick.comment
        for key in self.stick.coverage_dependence.keys():
            match = False
            for key2 in stick.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            assert match
        for species, parameters in self.stick.coverage_dependence.items():
            match = False
            for species2 in stick.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        assert value.value_si == stick.coverage_dependence[species2][key].value_si
            assert match
        assert dir(self.stick) == dir(stick)

    def test_repr(self):
        """
        Test that an StickingCoefficient object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec(f"stick = {self.stick!r}", globals(), namespace)
        assert "stick" in namespace
        stick = namespace["stick"]
        assert abs(self.stick.A.value - stick.A.value) < 1e0
        assert self.stick.A.units == stick.A.units
        assert round(abs(self.stick.n.value - stick.n.value), 4) == 0
        assert round(abs(self.stick.Ea.value - stick.Ea.value), 4) == 0
        assert self.stick.Ea.units == stick.Ea.units
        assert round(abs(self.stick.T0.value - stick.T0.value), 4) == 0
        assert self.stick.T0.units == stick.T0.units
        assert round(abs(self.stick.Tmin.value - stick.Tmin.value), 4) == 0
        assert self.stick.Tmin.units == stick.Tmin.units
        assert round(abs(self.stick.Tmax.value - stick.Tmax.value), 4) == 0
        assert self.stick.Tmax.units == stick.Tmax.units
        assert self.stick.comment == stick.comment
        assert [m.label for m in self.stick.coverage_dependence.keys()] == list(stick.coverage_dependence.keys())
        for species, parameters in self.stick.coverage_dependence.items():
            for key, value in parameters.items():
                assert value.value_si == stick.coverage_dependence[species.label][key].value_si
        assert dir(self.stick) == dir(stick)

    def test_copy(self):
        """
        Test that an StickingCoefficient object can be copied with deepcopy
        with no loss of information.
        """
        import copy

        stick = copy.deepcopy(self.stick)
        assert abs(self.stick.A.value - stick.A.value) < 1e0
        assert self.stick.A.units == stick.A.units
        assert round(abs(self.stick.n.value - stick.n.value), 4) == 0
        assert round(abs(self.stick.Ea.value - stick.Ea.value), 4) == 0
        assert self.stick.Ea.units == stick.Ea.units
        assert round(abs(self.stick.T0.value - stick.T0.value), 4) == 0
        assert self.stick.T0.units == stick.T0.units
        assert round(abs(self.stick.Tmin.value - stick.Tmin.value), 4) == 0
        assert self.stick.Tmin.units == stick.Tmin.units
        assert round(abs(self.stick.Tmax.value - stick.Tmax.value), 4) == 0
        assert self.stick.Tmax.units == stick.Tmax.units
        assert self.stick.comment == stick.comment
        for key in self.stick.coverage_dependence.keys():
            match = False
            for key2 in stick.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            assert match
        for species, parameters in self.stick.coverage_dependence.items():
            match = False
            for species2 in stick.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        assert value.value_si == stick.coverage_dependence[species2][key].value_si
            assert match
        assert dir(self.stick) == dir(stick)

    def test_is_identical_to(self):
        """
        Test that the StickingCoefficient.is_identical_to method works on itself
        """
        assert self.stick.is_identical_to(self.stick)


class TestSurfaceArrhenius:
    """
    Contains unit tests of the :class:`SurfaceArrhenius` class.
    """

    def setup_class(self):
        self.A = 1.44e18
        self.n = -0.087
        self.Ea = 63.4
        self.T0 = 1.0
        self.Tmin = 300.0
        self.Tmax = 3000.0
        s = Species().from_adjacency_list("1 X u0 p0 c0")
        s.label = "X"
        self.coverage_dependence = {
            s: {
                "a": quantity.Dimensionless(0.0),
                "m": quantity.Dimensionless(-1.0),
                "E": quantity.Energy(0.0, "J/mol"),
            }
        }
        self.comment = "CH3x + Hx <=> CH4 + x + x"
        self.surfarr = SurfaceArrhenius(
            A=(self.A, "m^2/(mol*s)"),
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
        assert abs(self.surfarr.A.value_si - self.A) < 1e0

    def test_n(self):
        """
        Test that the SurfaceArrhenius n property was properly set.
        """
        assert round(abs(self.surfarr.n.value_si - self.n), 6) == 0

    def test_Ea(self):
        """
        Test that the SurfaceArrhenius Ea property was properly set.
        """
        assert round(abs(self.surfarr.Ea.value_si * 0.001 - self.Ea), 6) == 0

    def test_T0(self):
        """
        Test that the SurfaceArrhenius T0 property was properly set.
        """
        assert round(abs(self.surfarr.T0.value_si - self.T0), 6) == 0

    def test_Tmin(self):
        """
        Test that the SurfaceArrhenius Tmin property was properly set.
        """
        assert round(abs(self.surfarr.Tmin.value_si - self.Tmin), 6) == 0

    def test_Tmax(self):
        """
        Test that the SurfaceArrhenius Tmax property was properly set.
        """
        assert round(abs(self.surfarr.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the SurfaceArrhenius comment property was properly set.
        """
        assert self.surfarr.comment == self.comment

    def test_coverage_dependence(self):
        """
        Test that the coverage dependent parameters was properly set.
        """
        for key in self.surfarr.coverage_dependence.keys():
            match = False
            for key2 in self.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            assert match
        for species, parameters in self.surfarr.coverage_dependence.items():
            match = False
            for species2 in self.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        assert value.value_si == self.coverage_dependence[species2][key].value_si
            assert match

    def test_is_temperature_valid(self):
        """
        Test the SurfaceArrhenius.is_temperature_valid() method.
        """
        T_data = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 4000])
        valid_data = np.array([False, True, True, True, True, True, True, True, True, False], bool)
        for T, valid in zip(T_data, valid_data):
            valid0 = self.surfarr.is_temperature_valid(T)
            assert valid0 == valid

    def test_pickle(self):
        """
        Test that an SurfaceArrhenius object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        surfarr = pickle.loads(pickle.dumps(self.surfarr, -1))
        assert abs(self.surfarr.A.value - surfarr.A.value) < 1e0
        assert self.surfarr.A.units == surfarr.A.units
        assert round(abs(self.surfarr.n.value - surfarr.n.value), 4) == 0
        assert round(abs(self.surfarr.Ea.value - surfarr.Ea.value), 4) == 0
        assert self.surfarr.Ea.units == surfarr.Ea.units
        assert round(abs(self.surfarr.T0.value - surfarr.T0.value), 4) == 0
        assert self.surfarr.T0.units == surfarr.T0.units
        assert round(abs(self.surfarr.Tmin.value - surfarr.Tmin.value), 4) == 0
        assert self.surfarr.Tmin.units == surfarr.Tmin.units
        assert round(abs(self.surfarr.Tmax.value - surfarr.Tmax.value), 4) == 0
        assert self.surfarr.Tmax.units == surfarr.Tmax.units
        assert self.surfarr.comment == surfarr.comment
        for key in self.surfarr.coverage_dependence.keys():
            match = False
            for key2 in surfarr.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            assert match
        for species, parameters in self.surfarr.coverage_dependence.items():
            match = False
            for species2 in surfarr.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        assert value.value_si == surfarr.coverage_dependence[species2][key].value_si
            assert match
        assert dir(self.surfarr) == dir(surfarr)

    def test_repr(self):
        """
        Test that an SurfaceArrhenius object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("surfarr = {0!r}".format(self.surfarr), globals(), namespace)
        assert "surfarr" in namespace
        surfarr = namespace["surfarr"]
        assert abs(self.surfarr.A.value - surfarr.A.value) < 1e0
        assert self.surfarr.A.units == surfarr.A.units
        assert round(abs(self.surfarr.n.value - surfarr.n.value), 4) == 0
        assert round(abs(self.surfarr.Ea.value - surfarr.Ea.value), 4) == 0
        assert self.surfarr.Ea.units == surfarr.Ea.units
        assert round(abs(self.surfarr.T0.value - surfarr.T0.value), 4) == 0
        assert self.surfarr.T0.units == surfarr.T0.units
        assert round(abs(self.surfarr.Tmin.value - surfarr.Tmin.value), 4) == 0
        assert self.surfarr.Tmin.units == surfarr.Tmin.units
        assert round(abs(self.surfarr.Tmax.value - surfarr.Tmax.value), 4) == 0
        assert self.surfarr.Tmax.units == surfarr.Tmax.units
        assert self.surfarr.comment == surfarr.comment
        assert [m.label for m in self.surfarr.coverage_dependence.keys()] == list(surfarr.coverage_dependence.keys())
        for species, parameters in self.surfarr.coverage_dependence.items():
            for key, value in parameters.items():
                assert value.value_si == surfarr.coverage_dependence[species.label][key].value_si
        assert dir(self.surfarr) == dir(surfarr)

    def test_copy(self):
        """
        Test that an SurfaceArrhenius object can be copied with deepcopy
        with no loss of information.
        """
        import copy

        surfarr = copy.deepcopy(self.surfarr)
        assert abs(self.surfarr.A.value - surfarr.A.value) < 1e0
        assert self.surfarr.A.units == surfarr.A.units
        assert round(abs(self.surfarr.n.value - surfarr.n.value), 4) == 0
        assert round(abs(self.surfarr.Ea.value - surfarr.Ea.value), 4) == 0
        assert self.surfarr.Ea.units == surfarr.Ea.units
        assert round(abs(self.surfarr.T0.value - surfarr.T0.value), 4) == 0
        assert self.surfarr.T0.units == surfarr.T0.units
        assert round(abs(self.surfarr.Tmin.value - surfarr.Tmin.value), 4) == 0
        assert self.surfarr.Tmin.units == surfarr.Tmin.units
        assert round(abs(self.surfarr.Tmax.value - surfarr.Tmax.value), 4) == 0
        assert self.surfarr.Tmax.units == surfarr.Tmax.units
        assert self.surfarr.comment == surfarr.comment
        for key in self.surfarr.coverage_dependence.keys():
            match = False
            for key2 in surfarr.coverage_dependence.keys():
                if key.is_identical(key2):
                    match = True
            assert match
        for species, parameters in self.surfarr.coverage_dependence.items():
            match = False
            for species2 in surfarr.coverage_dependence.keys():
                if species.is_identical(species2):
                    match = True
                    for key, value in parameters.items():
                        assert value.value_si == surfarr.coverage_dependence[species2][key].value_si
            assert match
        assert dir(self.surfarr) == dir(surfarr)

    def test_is_identical_to(self):
        """
        Test that the SurfaceArrhenius.is_identical_to method works on itself
        """
        assert self.surfarr.is_identical_to(self.surfarr)
