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
This module contains unit tests of the :mod:`rmgpy.pdep.collision` module.
"""


from rmgpy.pdep.collision import SingleExponentialDown


class TestSingleExponentialDown:
    """
    Contains unit tests of the SingleExponentialDown class.
    """

    def setup_class(self):
        self.alpha0 = 0.5
        self.T0 = 300.0
        self.n = 0.85
        self.singleExponentialDown = SingleExponentialDown(
            alpha0=(self.alpha0, "kJ/mol"),
            T0=(self.T0, "K"),
            n=self.n,
        )

    def test_alpha0(self):
        """
        Test the SingleExponentialDown.sigma attribute.
        """
        assert round(abs(self.singleExponentialDown.alpha0.value_si * 0.001 - self.alpha0), 4) == 0

    def test_temperature_0(self):
        """
        Test the SingleExponentialDown.T0 attribute.
        """
        assert round(abs(self.singleExponentialDown.T0.value_si - self.T0), 4) == 0

    def test_n(self):
        """
        Test the SingleExponentialDown.n attribute.
        """
        assert round(abs(self.singleExponentialDown.n - self.n), 4) == 0

    def test_get_alpha(self):
        """
        Test the SingleExponentialDown.get_alpha() method.
        """
        for T in [300, 400, 500, 600, 800, 1000, 1500, 2000]:
            dEdown0 = 1000.0 * self.alpha0 * (T / self.T0) ** self.n
            dEdown = self.singleExponentialDown.get_alpha(T)
            assert round(abs(dEdown0 - dEdown), 6) == 0

    def test_pickle(self):
        """
        Test that a SingleExponentialDown object can be successfully pickled
        and unpickled with no loss of information.
        """
        import pickle

        singleExponentialDown = pickle.loads(pickle.dumps(self.singleExponentialDown, -1))
        assert round(abs(self.singleExponentialDown.alpha0.value - singleExponentialDown.alpha0.value), 6) == 0
        assert self.singleExponentialDown.alpha0.units == singleExponentialDown.alpha0.units
        assert round(abs(self.singleExponentialDown.T0.value - singleExponentialDown.T0.value), 6) == 0
        assert self.singleExponentialDown.T0.units == singleExponentialDown.T0.units
        assert round(abs(self.singleExponentialDown.n - singleExponentialDown.n), 4) == 0

    def test_repr(self):
        """
        Test that a SingleExponentialDown object can be successfully
        reconstructed from its repr() with no loss of information.
        """
        namespace = {}
        exec(
            "singleExponentialDown = {0!r}".format(self.singleExponentialDown),
            globals(),
            namespace,
        )
        assert "singleExponentialDown" in namespace
        singleExponentialDown = namespace["singleExponentialDown"]
        assert round(abs(self.singleExponentialDown.alpha0.value - singleExponentialDown.alpha0.value), 6) == 0
        assert self.singleExponentialDown.alpha0.units == singleExponentialDown.alpha0.units
        assert round(abs(self.singleExponentialDown.T0.value - singleExponentialDown.T0.value), 6) == 0
        assert self.singleExponentialDown.T0.units == singleExponentialDown.T0.units
        assert round(abs(self.singleExponentialDown.n - singleExponentialDown.n), 4) == 0
