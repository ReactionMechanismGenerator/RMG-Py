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
This script contains unit tests of the :mod:`rmgpy.kinetics.tunneling` module.
"""


import numpy as np

from rmgpy.kinetics.tunneling import Wigner, Eckart


class TestWigner:
    """
    Contains unit tests of the :class:`Wigner` class.
    """

    def setup_class(self):
        self.frequency = -2017.96
        self.tunneling = Wigner(
            frequency=(self.frequency, "cm^-1"),
        )

    def test_frequency(self):
        """
        Test that the Wigner frequency property was properly set.
        """
        assert round(abs(self.tunneling.frequency.value_si - self.frequency), 4) == 0

    def test_calculate_tunneling_factor(self):
        """
        Test the Wigner.calculate_tunneling_factor() method.
        """
        Tlist = np.array([300, 500, 1000, 1500, 2000])
        kexplist = np.array([4.90263, 2.40495, 1.35124, 1.15611, 1.08781])
        for T, kexp in zip(Tlist, kexplist):
            kact = self.tunneling.calculate_tunneling_factor(T)
            assert round(abs(kexp - kact), 4) == 0

    def test_pickle(self):
        """
        Test that a Wigner object can be successfully pickled and unpickled
        with no loss of information.
        """
        import pickle

        tunneling = pickle.loads(pickle.dumps(self.tunneling, -1))
        assert round(abs(self.tunneling.frequency.value - tunneling.frequency.value), 2) == 0
        assert self.tunneling.frequency.units == tunneling.frequency.units

    def test_repr(self):
        """
        Test that a Wigner object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        namespace = {}
        exec("tunneling = {0!r}".format(self.tunneling), globals(), namespace)
        assert "tunneling" in namespace
        tunneling = namespace["tunneling"]
        assert round(abs(self.tunneling.frequency.value - tunneling.frequency.value), 2) == 0
        assert self.tunneling.frequency.units == tunneling.frequency.units


class TestEckart:
    """
    Contains unit tests of the :class:`Eckart` class.
    """

    def setup_class(self):
        self.frequency = -2017.96
        self.E0_reac = -295.563
        self.E0_TS = -12.7411
        self.E0_prod = (-10.2664) + (-253.48)
        self.tunneling = Eckart(
            frequency=(self.frequency, "cm^-1"),
            E0_reac=(self.E0_reac, "kJ/mol"),
            E0_TS=(self.E0_TS, "kJ/mol"),
            E0_prod=(self.E0_prod, "kJ/mol"),
        )

    def test_frequency(self):
        """
        Test that the Eckart frequency property was properly set.
        """
        assert round(abs(self.tunneling.frequency.value_si - self.frequency), 4) == 0

    def test_e0_reac(self):
        """
        Test that the Eckart E0_reac property was properly set.
        """
        assert round(abs(self.tunneling.E0_reac.value_si * 0.001 - self.E0_reac), 4) == 0

    def test_e0_ts(self):
        """
        Test that the Eckart E0_TS property was properly set.
        """
        assert round(abs(self.tunneling.E0_TS.value_si * 0.001 - self.E0_TS), 4) == 0

    def test_e0_prod(self):
        """
        Test that the Eckart E0_prod property was properly set.
        """
        assert round(abs(self.tunneling.E0_prod.value_si * 0.001 - self.E0_prod), 4) == 0

    def test_calculate_tunneling_factor(self):
        """
        Test the Eckart.calculate_tunneling_factor() method.
        """
        Tlist = np.array([300, 500, 1000, 1500, 2000])
        kexplist = np.array([1623051.0, 7.69349, 1.46551, 1.18111, 1.09858])
        for T, kexp in zip(Tlist, kexplist):
            kact = self.tunneling.calculate_tunneling_factor(T)
            assert abs(kexp - kact) < 1e-3 * kexp

    def test_pickle(self):
        """
        Test that an Eckart object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        tunneling = pickle.loads(pickle.dumps(self.tunneling, -1))
        assert round(abs(self.tunneling.frequency.value - tunneling.frequency.value), 2) == 0
        assert self.tunneling.frequency.units == tunneling.frequency.units
        assert round(abs(self.tunneling.E0_reac.value - tunneling.E0_reac.value), 3) == 0
        assert self.tunneling.E0_reac.units == tunneling.E0_reac.units
        assert round(abs(self.tunneling.E0_TS.value - tunneling.E0_TS.value), 3) == 0
        assert self.tunneling.E0_TS.units == tunneling.E0_TS.units
        assert round(abs(self.tunneling.E0_prod.value - tunneling.E0_prod.value), 3) == 0
        assert self.tunneling.E0_prod.units == tunneling.E0_prod.units

    def test_repr(self):
        """
        Test that an Eckart object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        namespace = {}
        exec("tunneling = {0!r}".format(self.tunneling), globals(), namespace)
        assert "tunneling" in namespace
        tunneling = namespace["tunneling"]
        assert round(abs(self.tunneling.frequency.value - tunneling.frequency.value), 2) == 0
        assert self.tunneling.frequency.units == tunneling.frequency.units
        assert round(abs(self.tunneling.E0_reac.value - tunneling.E0_reac.value), 3) == 0
        assert self.tunneling.E0_reac.units == tunneling.E0_reac.units
        assert round(abs(self.tunneling.E0_TS.value - tunneling.E0_TS.value), 3) == 0
        assert self.tunneling.E0_TS.units == tunneling.E0_TS.units
        assert round(abs(self.tunneling.E0_prod.value - tunneling.E0_prod.value), 3) == 0
        assert self.tunneling.E0_prod.units == tunneling.E0_prod.units
