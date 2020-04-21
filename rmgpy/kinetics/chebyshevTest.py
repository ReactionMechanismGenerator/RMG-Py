#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This script contains unit tests of the :mod:`rmgpy.kinetics.chebyshev` module.
"""

import unittest

import numpy as np

from rmgpy.exceptions import KineticsError
from rmgpy.kinetics.chebyshev import Chebyshev


################################################################################

class TestChebyshev(unittest.TestCase):
    """
    Contains unit tests of the Chebyshev class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Tmin = 300.
        self.Tmax = 2000.
        self.Pmin = 0.01
        self.Pmax = 100.
        self.coeffs = np.array([
            [11.67723, 0.729281, -0.11984, 0.00882175],
            [-1.02669, 0.853639, -0.0323485, -0.027367],
            [-0.447011, 0.244144, 0.0559122, -0.0101723],
            [-0.128261, 0.0111596, 0.0281176, 0.00604353],
            [-0.0117034, -0.0235646, 0.00061009, 0.00401309],
            [0.0155433, -0.0136846, -0.00463048, -0.000261353],
        ])
        self.comment = """acetyl + O2 -> acetylperoxy"""
        self.chebyshev = Chebyshev(
            coeffs=self.coeffs,
            kunits="cm^3/(mol*s)",
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )

    def test_coeffs(self):
        """
        Test that the Chebyshev coeffs property was properly set.
        """
        self.assertEqual(self.chebyshev.coeffs.value.shape, self.coeffs.shape)
        for i in range(self.chebyshev.coeffs.value.shape[0]):
            for j in range(self.chebyshev.coeffs.value.shape[1]):
                C0 = float(self.coeffs[i, j])
                C = float(self.chebyshev.coeffs.value_si[i, j])
                if i == 0 and j == 0: C0 -= 6  # Unit conversion from cm^3/(mol*s) to m^3/(mol*s)
                self.assertAlmostEqual(C0, C, delta=1e-6 * C0)

    def test_temperature_min(self):
        """
        Test that the Chebyshev Tmin property was properly set.
        """
        self.assertAlmostEqual(self.chebyshev.Tmin.value_si, self.Tmin, 6)

    def test_temperature_max(self):
        """
        Test that the Chebyshev Tmax property was properly set.
        """
        self.assertAlmostEqual(self.chebyshev.Tmax.value_si, self.Tmax, 6)

    def test_pressure_min(self):
        """
        Test that the Chebyshev Pmin property was properly set.
        """
        self.assertAlmostEqual(self.chebyshev.Pmin.value_si * 1e-5, self.Pmin, 6)

    def test_pressure_max(self):
        """
        Test that the Chebyshev Pmax property was properly set.
        """
        self.assertAlmostEqual(self.chebyshev.Pmax.value_si * 1e-5, self.Pmax, 6)

    def test_comment(self):
        """
        Test that the Chebyshev comment property was properly set.
        """
        self.assertEqual(self.chebyshev.comment, self.comment)

    def test_is_pressure_dependent(self):
        """
        Test the Chebyshev.is_pressure_dependent() method.
        
        """
        self.assertTrue(self.chebyshev.is_pressure_dependent())

    def test_get_rate_coefficient(self):
        """
        Test the Chebyshev.get_rate_coefficient() method.
        """
        Tlist = np.array([300, 500, 1000, 1500])
        Plist = np.array([1e4, 1e5, 1e6])
        Kexp = np.array([
            [2.29100e+06, 2.58452e+06, 2.57204e+06],
            [1.10198e+06, 2.04037e+06, 2.57428e+06],
            [4.37919e+04, 2.36481e+05, 8.57727e+05],
            [5.20144e+03, 4.10123e+04, 2.50401e+05],
        ])
        for t in range(Tlist.shape[0]):
            for p in range(Plist.shape[0]):
                Kact = self.chebyshev.get_rate_coefficient(Tlist[t], Plist[p])
                self.assertAlmostEqual(Kact / Kexp[t, p], 1.0, 4, '{0} != {1} within 4 places'.format(Kexp[t, p], Kact))

    def test_fit_to_data(self):
        """
        Test the Chebyshev.fit_to_data() method.
        """
        Tdata = np.array(
            [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000])
        Pdata = np.array([3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e7])
        nT = len(Tdata)
        nP = len(Pdata)
        kdata = np.zeros((nT, nP))
        for t in range(nT):
            for p in range(nP):
                kdata[t, p] = self.chebyshev.get_rate_coefficient(Tdata[t], Pdata[p]) * 1e6
        chebyshev = Chebyshev().fit_to_data(Tdata, Pdata, kdata, kunits="cm^3/(mol*s)", degreeT=6, degreeP=4, Tmin=300,
                                          Tmax=2000, Pmin=0.1, Pmax=10.)
        for t in range(nT):
            for p in range(nP):
                kfit = chebyshev.get_rate_coefficient(Tdata[t], Pdata[p]) * 1e6
                self.assertAlmostEqual(kfit, kdata[t, p], delta=1e-4 * kdata[t, p])

    def test_fit_to_data2(self):
        """
        Test the Chebyshev.fit_to_data() method throws error without enough degrees of freedom.
        
        Here only 3 temperatures are given, but the polynomial desired has 6 parameters.
        """
        Tdata = np.array([300, 1200, 2000])
        Pdata = np.array([1e5, 3e5, 1e6, 3e7])
        nT = len(Tdata)
        nP = len(Pdata)
        kdata = np.zeros((nT, nP))
        for t in range(nT):
            for p in range(nP):
                kdata[t, p] = self.chebyshev.get_rate_coefficient(Tdata[t], Pdata[p])
        with self.assertRaises(KineticsError):
            Chebyshev().fit_to_data(Tdata, Pdata, kdata, kunits="cm^3/(mol*s)", degreeT=12, degreeP=8, Tmin=300,
                                  Tmax=2000, Pmin=0.1, Pmax=10.)

    def test_pickle(self):
        """
        Test that a Chebyshev object can be pickled and unpickled with no loss
        of information.
        """
        import pickle
        chebyshev = pickle.loads(pickle.dumps(self.chebyshev, -1))
        self.assertEqual(self.chebyshev.coeffs.value.shape[0], chebyshev.coeffs.value.shape[0])
        self.assertEqual(self.chebyshev.coeffs.value.shape[1], chebyshev.coeffs.value.shape[1])
        for i in range(self.chebyshev.coeffs.value.shape[0]):
            for j in range(self.chebyshev.coeffs.value.shape[1]):
                C0 = self.chebyshev.coeffs.value_si[i, j]
                C = chebyshev.coeffs.value_si[i, j]
                self.assertAlmostEqual(C0, C, delta=1e-4 * C0)
        self.assertAlmostEqual(self.chebyshev.Tmin.value, chebyshev.Tmin.value, 4)
        self.assertEqual(self.chebyshev.Tmin.units, chebyshev.Tmin.units)
        self.assertAlmostEqual(self.chebyshev.Tmax.value, chebyshev.Tmax.value, 4)
        self.assertEqual(self.chebyshev.Tmax.units, chebyshev.Tmax.units)
        self.assertAlmostEqual(self.chebyshev.Pmin.value, chebyshev.Pmin.value, 4)
        self.assertEqual(self.chebyshev.Pmin.units, chebyshev.Pmin.units)
        self.assertAlmostEqual(self.chebyshev.Pmax.value, chebyshev.Pmax.value, 4)
        self.assertEqual(self.chebyshev.Pmax.units, chebyshev.Pmax.units)
        self.assertEqual(self.chebyshev.comment, chebyshev.comment)

    def test_repr(self):
        """
        Test that a Chebyshev object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec('chebyshev = {0!r}'.format(self.chebyshev), globals(), namespace)
        self.assertIn('chebyshev', namespace)
        chebyshev = namespace['chebyshev']
        self.assertEqual(self.chebyshev.coeffs.value.shape[0], chebyshev.coeffs.value.shape[0])
        self.assertEqual(self.chebyshev.coeffs.value.shape[1], chebyshev.coeffs.value.shape[1])
        for i in range(self.chebyshev.coeffs.value.shape[0]):
            for j in range(self.chebyshev.coeffs.value.shape[1]):
                C0 = self.chebyshev.coeffs.value[i, j]
                C = chebyshev.coeffs.value[i, j]
                self.assertAlmostEqual(C0, C, delta=1e-4 * C0)
        self.assertAlmostEqual(self.chebyshev.Tmin.value, chebyshev.Tmin.value, 4)
        self.assertEqual(self.chebyshev.Tmin.units, chebyshev.Tmin.units)
        self.assertAlmostEqual(self.chebyshev.Tmax.value, chebyshev.Tmax.value, 4)
        self.assertEqual(self.chebyshev.Tmax.units, chebyshev.Tmax.units)
        self.assertAlmostEqual(self.chebyshev.Pmin.value, chebyshev.Pmin.value, 4)
        self.assertEqual(self.chebyshev.Pmin.units, chebyshev.Pmin.units)
        self.assertAlmostEqual(self.chebyshev.Pmax.value, chebyshev.Pmax.value, 4)
        self.assertEqual(self.chebyshev.Pmax.units, chebyshev.Pmax.units)
        self.assertEqual(self.chebyshev.comment, chebyshev.comment)

    def test_change_rate(self):
        """
        Test the Chebyshev.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.chebyshev.get_rate_coefficient(T, 1e5) for T in Tlist])
        self.chebyshev.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.chebyshev.get_rate_coefficient(T, 1e5)
            self.assertAlmostEqual(2 * kexp, kact, delta=1e-6 * kexp)

    def test_is_identical_to(self):
        """
        Test the Chebyshev.is_identical_to() method.
        """
        # Trivial case, compare to a KineticsModel
        from rmgpy.kinetics.model import KineticsModel
        self.assertFalse(self.chebyshev.is_identical_to(KineticsModel()))

        # Compare to identical Chebyshev
        new_chebyshev = Chebyshev(
            coeffs=self.coeffs,
            kunits="cm^3/(mol*s)",
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )
        self.assertTrue(self.chebyshev.is_identical_to(new_chebyshev))

        # Compare to Chebyshev with different Tmin/Tmax
        new_chebyshev = Chebyshev(
            coeffs=self.coeffs,
            kunits="cm^3/(mol*s)",
            Tmin=(200, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )
        self.assertFalse(self.chebyshev.is_identical_to(new_chebyshev))

        new_chebyshev = Chebyshev(
            coeffs=self.coeffs,
            kunits="cm^3/(mol*s)",
            Tmin=(self.Tmin, "K"),
            Tmax=(2500, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )
        self.assertFalse(self.chebyshev.is_identical_to(new_chebyshev))

        # Compare to Chebyshev with different degreeT/degreeP
        new_chebyshev = Chebyshev(
            coeffs=self.coeffs[0:-1, :],  # Remove one T dimension
            kunits="cm^3/(mol*s)",
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )
        self.assertFalse(self.chebyshev.is_identical_to(new_chebyshev))

        new_chebyshev = Chebyshev(
            coeffs=self.coeffs[:, 0:-1],  # Remove one P dimension
            kunits="cm^3/(mol*s)",
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )
        self.assertFalse(self.chebyshev.is_identical_to(new_chebyshev))

        # Compare to Chebyshev with different units
        new_chebyshev = Chebyshev(
            coeffs=self.coeffs,
            kunits="m^3/(mol*s)",
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )
        self.assertFalse(self.chebyshev.is_identical_to(new_chebyshev))

        # Compare to Chebyshev with slightly different coefficients
        new_chebyshev = Chebyshev(
            coeffs=np.copy(self.coeffs) * 0.01,
            kunits="cm^3/(mol*s)",
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )
        self.assertFalse(self.chebyshev.is_identical_to(new_chebyshev))


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
