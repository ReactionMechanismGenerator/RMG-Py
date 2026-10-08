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
This script contains unit tests of the :mod:`rmgpy.kinetics.model` module.
"""


import pytest
import quantities as pq

from rmgpy.kinetics.model import (
    get_reaction_order_from_rate_coefficient_units,
    get_rate_coefficient_units_from_reaction_order,
    KineticsModel,
)
from rmgpy.kinetics.uncertainties import RateUncertainty


class TestKineticsModel:
    """
    Contains unit tests of the KineticsModel class
    """

    def setup_class(self):
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.Pmin = 0.1
        self.Pmax = 100.0
        self.comment = "foo bar"
        self.uncertainty = RateUncertainty(mu=0.3, var=0.6, Tref=1000.0, N=1, correlation="ab")
        self.km = KineticsModel(
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            uncertainty=self.uncertainty,
            comment=self.comment,
        )

    def test_is_identical_to(self):
        """
        Test that the KineticsModel.is_identical_to method works on itself.

        This just checks the Temperature range
        """
        assert self.km.is_identical_to(self.km)

        import copy

        km = copy.deepcopy(self.km)
        assert self.km.is_identical_to(self.km)

        km.Tmax = (self.Tmax - 50, "K")  # discrepancy must be more than 1%!
        assert not self.km.is_identical_to(km)

    def test_repr(self):
        """
        Test that an KineticsModel object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("km = {0!r}".format(self.km), globals(), namespace)
        assert "km" in namespace
        km = namespace["km"]
        assert self.km.is_identical_to(km)
        assert dir(self.km) == dir(km)
        for att in "Tmax Tmin Pmax Pmin comment uncertainty".split():
            assert repr(getattr(self.km, att)) == repr(getattr(km, att))

    def test_pickle(self):
        """
        Test that an KineticsModel object can be pickled and unpickled
        with no loss of information.
        """
        import pickle

        km = pickle.loads(pickle.dumps(self.km, -1))
        assert self.km.is_identical_to(km)
        assert dir(self.km) == dir(km)
        for att in "Tmax Tmin Pmax Pmin comment uncertainty".split():
            assert repr(getattr(self.km, att)) == repr(getattr(km, att))


class TestOrder:
    """
    Contains unit tests of the functions for converting rate coefficient units
    to/from reaction orders.
    """

    def test_to_order_zeroth(self):
        """
        Test the conversion of zeroth-order rate coefficient units to an integer
        reaction order.
        """
        assert 0 == get_reaction_order_from_rate_coefficient_units("mol/(m^3*s)")
        assert 0 == get_reaction_order_from_rate_coefficient_units("mol/(cm^3*s)")
        assert 0 == get_reaction_order_from_rate_coefficient_units("molecule/(m^3*s)")
        assert 0 == get_reaction_order_from_rate_coefficient_units("molecule/(cm^3*s)")

    def test_to_order_first(self):
        """
        Test the conversion of first-order rate coefficient units to an integer
        reaction order.
        """
        assert 1 == get_reaction_order_from_rate_coefficient_units("s^-1")

    def test_to_order_second(self):
        """
        Test the conversion of second-order rate coefficient units to an integer
        reaction order.
        """
        assert 2 == get_reaction_order_from_rate_coefficient_units("m^3/(mol*s)")
        assert 2 == get_reaction_order_from_rate_coefficient_units("cm^3/(mol*s)")
        assert 2 == get_reaction_order_from_rate_coefficient_units("m^3/(molecule*s)")
        assert 2 == get_reaction_order_from_rate_coefficient_units("cm^3/(molecule*s)")

    def test_to_order_third(self):
        """
        Test the conversion of third-order rate coefficient units to an integer
        reaction order.
        """
        assert 3 == get_reaction_order_from_rate_coefficient_units("m^6/(mol^2*s)")
        assert 3 == get_reaction_order_from_rate_coefficient_units("cm^6/(mol^2*s)")
        assert 3 == get_reaction_order_from_rate_coefficient_units("m^6/(molecule^2*s)")
        assert 3 == get_reaction_order_from_rate_coefficient_units("cm^6/(molecule^2*s)")

    def test_to_order_fifth(self):
        """
        Test the conversion of fifth-order rate coefficient units to an integer
        reaction order.
        """
        assert 5 == get_reaction_order_from_rate_coefficient_units("m^12/(mol^4*s)")
        assert 5 == get_reaction_order_from_rate_coefficient_units("cm^12/(mol^4*s)")
        assert 5 == get_reaction_order_from_rate_coefficient_units("m^12/(molecule^4*s)")
        assert 5 == get_reaction_order_from_rate_coefficient_units("cm^12/(molecule^4*s)")

    def test_to_order_invalid(self):
        """
        Test that units which do not belong to a whole, non-negative reaction
        order raise a ValueError, and are not rounded to the nearest order.
        """
        for kunits in ["m^4/(mol*s)", "m^1.5/(mol^0.5*s)", "mol^2/(m^6*s)", "m^12/(mol^4*s^2)"]:
            with pytest.raises(ValueError):
                get_reaction_order_from_rate_coefficient_units(kunits)

    def test_to_units_zeroth(self):
        """
        Test the conversion of a reaction order of zero to rate coefficient
        units.
        """
        assert "mol/(m^3*s)" == get_rate_coefficient_units_from_reaction_order(0)

    def test_to_units_first(self):
        """
        Test the conversion of a reaction order of one to rate coefficient
        units.
        """
        assert "s^-1" == get_rate_coefficient_units_from_reaction_order(1)

    def test_to_units_second(self):
        """
        Test the conversion of a reaction order of two to rate coefficient
        units.
        """
        assert "m^3/(mol*s)" == get_rate_coefficient_units_from_reaction_order(2)

    def test_to_units_third(self):
        """
        Test the conversion of a reaction order of three to rate coefficient
        units.
        """
        assert "m^6/(mol^2*s)" == get_rate_coefficient_units_from_reaction_order(3)

    def test_to_units_fifth(self):
        """
        Test the conversion of a reaction order of five to rate coefficient
        units, as for the reverse of a reaction with five gas-phase products.
        """
        assert "m^12/(mol^4*s)" == get_rate_coefficient_units_from_reaction_order(5)

    def test_to_units_surface(self):
        """
        Test the conversion of gas and surface species counts that have no
        explicit case to rate coefficient units.
        """
        assert "m^12/(mol^4*s)" == get_rate_coefficient_units_from_reaction_order(4, 1)
        assert "m^14/(mol^5*s)" == get_rate_coefficient_units_from_reaction_order(4, 2)
        assert "m^6/(mol^3*s)" == get_rate_coefficient_units_from_reaction_order(0, 4)

    def test_to_units_negative(self):
        """
        Test that a negative number of gas or surface species raises a
        ValueError.
        """
        with pytest.raises(ValueError):
            get_rate_coefficient_units_from_reaction_order(-1)
        with pytest.raises(ValueError):
            get_rate_coefficient_units_from_reaction_order(2, -1)

    def test_to_units_any_order(self):
        """
        Test that the units for any number of gas and surface species are those
        of the rate, mol/(m^3*s) or mol/(m^2*s) once a surface species takes
        part, divided by the concentration of each species, mol/m^3 for a gas
        and mol/m^2 for a surface species. This covers the explicit cases and
        the general formula.
        """
        for n_gas in range(8):
            for n_surf in range(6):
                kunits = get_rate_coefficient_units_from_reaction_order(n_gas, n_surf)
                rate_m = -3 if n_surf == 0 else -2
                expected = {
                    pq.m: rate_m + 3 * n_gas + 2 * n_surf,
                    pq.mol: 1 - n_gas - n_surf,
                    pq.s: -1,
                }
                expected = {unit: power for unit, power in expected.items() if power != 0}
                assert expected == dict(pq.Quantity(1.0, kunits).simplified.dimensionality), kunits

    def test_round_trip(self):
        """
        Test that the units for a gas-phase reaction order convert back to the
        same order.
        """
        for order in range(10):
            kunits = get_rate_coefficient_units_from_reaction_order(order)
            assert order == get_reaction_order_from_rate_coefficient_units(kunits)
