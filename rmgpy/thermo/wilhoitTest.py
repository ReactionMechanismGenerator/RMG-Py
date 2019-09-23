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
This script contains unit tests of the :mod:`rmgpy.thermo.wilhoit` module.
"""

from __future__ import division

import os.path
import unittest

import numpy as np

import rmgpy.constants as constants
from rmgpy.quantity import ScalarQuantity
from rmgpy.thermo.wilhoit import Wilhoit

################################################################################


class TestWilhoit(unittest.TestCase):
    """
    Contains unit tests of the :class:`Wilhoit` class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Cp0 = 4.0
        self.CpInf = 21.5
        self.a0 = 0.0977518
        self.a1 = -16.3067
        self.a2 = 26.2524
        self.a3 = -12.6785
        self.B = 1068.68
        self.H0 = -94088.  # -782.292 kJ/mol / constants.R
        self.S0 = -118.46  # -984.932 J/mol*K / constants.R
        self.Tmin = 300.
        self.Tmax = 3000.
        self.comment = 'C2H6'
        self.wilhoit = Wilhoit(
            Cp0=(self.Cp0 * constants.R, "J/(mol*K)"),
            CpInf=(self.CpInf * constants.R, "J/(mol*K)"),
            a0=self.a0,
            a1=self.a1,
            a2=self.a2,
            a3=self.a3,
            B=(self.B, "K"),
            H0=(self.H0 * 0.001 * constants.R, "kJ/mol"),
            S0=(self.S0 * constants.R, "J/(mol*K)"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def tearDown(self):
        """
        Reset the database & liquid parameters for solution
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    def test_Cp0(self):
        """
        Test that the Wilhoit Cp0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.Cp0.value_si / constants.R, self.Cp0, 6)

    def test_CpInf(self):
        """
        Test that the Wilhoit CpInf property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.CpInf.value_si / constants.R, self.CpInf, 6)

    def test_a0(self):
        """
        Test that the Wilhoit a0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a0, self.a0, 6)

    def test_a1(self):
        """
        Test that the Wilhoit a1 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a1, self.a1, 6)

    def test_a2(self):
        """
        Test that the Wilhoit a2 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a2, self.a2, 6)

    def test_a3(self):
        """
        Test that the Wilhoit a3 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a3, self.a3, 6)

    def test_B(self):
        """
        Test that the Wilhoit B property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.B.value_si, self.B, 6)

    def test_H0(self):
        """
        Test that the Wilhoit H0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.H0.value_si / constants.R, self.H0, 6)

    def test_S0(self):
        """
        Test that the Wilhoit S0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.S0.value_si / constants.R, self.S0, 6)

    def test_Tmin(self):
        """
        Test that the Wilhoit Tmin property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.Tmin.value_si, self.Tmin, 6)

    def test_Tmax(self):
        """
        Test that the Wilhoit Tmax property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.Tmax.value_si, self.Tmax, 6)

    def test_E0(self):
        """
        Test that the Wilhoit E0 property is properly calculated from Enthalpy at 0.001 K
        """
        self.assertAlmostEqual(self.wilhoit.E0.value_si, self.wilhoit.getEnthalpy(0.001), 1)

    def test_comment(self):
        """
        Test that the Wilhoit comment property was properly set.
        """
        self.assertEqual(self.wilhoit.comment, self.comment)

    def test_isTemperatureValid(self):
        """
        Test the Wilhoit.isTemperatureValid() method.
        """
        Tdata = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
        valid_data = [False, True, True, True, True, True, True, True, True, True]
        for T, valid in zip(Tdata, valid_data):
            valid0 = self.wilhoit.isTemperatureValid(T)
            self.assertEqual(valid0, valid)

    def test_getHeatCapacity(self):
        """
        Test the Wilhoit.getHeatCapacity() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        cp_exp_list = np.array([5.12003, 7.80327, 10.5528, 12.8323, 14.6013,
                                15.9511, 16.9842, 17.7837, 18.4114, 18.9117]) * constants.R
        for T, cp_exp in zip(Tlist, cp_exp_list):
            cp_act = self.wilhoit.getHeatCapacity(T)
            self.assertAlmostEqual(cp_exp / cp_act, 1.0, 3, '{0} != {1} within 3 places'.format(cp_exp, cp_act))

    def test_getEnthalpy(self):
        """
        Test the Wilhoit.getEnthalpy() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        h_exp_list = np.array([-51.9303, -22.7609, -12.1050, -6.14444, -2.16433,
                               0.747500, 2.99646, 4.79698, 6.27618, 7.51564]) * constants.R * Tlist
        for T, h_exp in zip(Tlist, h_exp_list):
            h_act = self.wilhoit.getEnthalpy(T)
            self.assertAlmostEqual(h_exp / h_act, 1.0, 3, '{0} != {1}'.format(h_exp, h_act))

    def test_getEntropy(self):
        """
        Test the Wilhoit.getEntropy() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        s_exp_list = np.array([25.3095, 29.6445, 33.3398, 36.7006, 39.7629,
                               42.5499, 45.0898, 47.4122, 49.5445, 51.5112]) * constants.R
        for T, s_exp in zip(Tlist, s_exp_list):
            s_act = self.wilhoit.getEntropy(T)
            self.assertAlmostEqual(s_exp / s_act, 1.0, 4, '{0} != {1}'.format(s_exp, s_act))

    def test_getFreeEnergy(self):
        """
        Test the Wilhoit.getFreeEnergy() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        for T in Tlist:
            g_exp = self.wilhoit.getEnthalpy(T) - T * self.wilhoit.getEntropy(T)
            g_act = self.wilhoit.getFreeEnergy(T)
            self.assertAlmostEqual(g_exp / g_act, 1.0, 4, '{0} != {1}'.format(g_exp, g_act))

    def test_pickle(self):
        """
        Test that a Wilhoit object can be pickled and unpickled with no loss
        of information.
        """
        import pickle
        wilhoit = pickle.loads(pickle.dumps(self.wilhoit))
        self.assertAlmostEqual(self.wilhoit.Cp0.value, wilhoit.Cp0.value, 4)
        self.assertEqual(self.wilhoit.Cp0.units, wilhoit.Cp0.units)
        self.assertAlmostEqual(self.wilhoit.CpInf.value, wilhoit.CpInf.value, 3)
        self.assertEqual(self.wilhoit.CpInf.units, wilhoit.CpInf.units)
        self.assertAlmostEqual(self.wilhoit.a0, wilhoit.a0, 4)
        self.assertAlmostEqual(self.wilhoit.a1, wilhoit.a1, 4)
        self.assertAlmostEqual(self.wilhoit.a2, wilhoit.a2, 4)
        self.assertAlmostEqual(self.wilhoit.a3, wilhoit.a3, 4)
        self.assertAlmostEqual(self.wilhoit.B.value, wilhoit.B.value, 4)
        self.assertEqual(self.wilhoit.B.units, wilhoit.B.units)
        self.assertAlmostEqual(self.wilhoit.H0.value, wilhoit.H0.value, 4)
        self.assertEqual(self.wilhoit.H0.units, wilhoit.H0.units)
        self.assertAlmostEqual(self.wilhoit.S0.value, wilhoit.S0.value, 3)
        self.assertEqual(self.wilhoit.S0.units, wilhoit.S0.units)
        self.assertAlmostEqual(self.wilhoit.Tmin.value, wilhoit.Tmin.value, 4)
        self.assertEqual(self.wilhoit.Tmin.units, wilhoit.Tmin.units)
        self.assertAlmostEqual(self.wilhoit.Tmax.value, wilhoit.Tmax.value, 4)
        self.assertEqual(self.wilhoit.Tmax.units, wilhoit.Tmax.units)
        self.assertAlmostEqual(self.wilhoit.E0.value, wilhoit.E0.value, 4)
        self.assertEqual(self.wilhoit.E0.units, wilhoit.E0.units)
        self.assertEqual(self.wilhoit.comment, wilhoit.comment)

    def test_repr(self):
        """
        Test that a Wilhoit object can be reconstructed from its repr() output
        with no loss of information.
        """
        namespace = {}
        exec('wilhoit = {0!r}'.format(self.wilhoit), globals(), namespace)
        self.assertIn('wilhoit', namespace)
        wilhoit = namespace['wilhoit']
        self.assertAlmostEqual(self.wilhoit.Cp0.value, wilhoit.Cp0.value, 4)
        self.assertEqual(self.wilhoit.Cp0.units, wilhoit.Cp0.units)
        self.assertAlmostEqual(self.wilhoit.CpInf.value, wilhoit.CpInf.value, 3)
        self.assertEqual(self.wilhoit.CpInf.units, wilhoit.CpInf.units)
        self.assertAlmostEqual(self.wilhoit.a0, wilhoit.a0, 4)
        self.assertAlmostEqual(self.wilhoit.a1, wilhoit.a1, 4)
        self.assertAlmostEqual(self.wilhoit.a2, wilhoit.a2, 4)
        self.assertAlmostEqual(self.wilhoit.a3, wilhoit.a3, 4)
        self.assertAlmostEqual(self.wilhoit.B.value, wilhoit.B.value, 4)
        self.assertEqual(self.wilhoit.B.units, wilhoit.B.units)
        self.assertAlmostEqual(self.wilhoit.H0.value, wilhoit.H0.value, 4)
        self.assertEqual(self.wilhoit.H0.units, wilhoit.H0.units)
        self.assertAlmostEqual(self.wilhoit.S0.value, wilhoit.S0.value, 3)
        self.assertEqual(self.wilhoit.S0.units, wilhoit.S0.units)
        self.assertAlmostEqual(self.wilhoit.Tmin.value, wilhoit.Tmin.value, 4)
        self.assertEqual(self.wilhoit.Tmin.units, wilhoit.Tmin.units)
        self.assertAlmostEqual(self.wilhoit.Tmax.value, wilhoit.Tmax.value, 4)
        self.assertEqual(self.wilhoit.Tmax.units, wilhoit.Tmax.units)
        self.assertAlmostEqual(self.wilhoit.E0.value, wilhoit.E0.value, 1)
        self.assertEqual(self.wilhoit.E0.units, wilhoit.E0.units)
        self.assertEqual(self.wilhoit.comment, wilhoit.comment)

    def test_fitToData(self):
        """
        Test the Wilhoit.fitToData() method.
        """
        h298 = self.wilhoit.getEnthalpy(298)
        s298 = self.wilhoit.getEntropy(298)
        Tdata = np.array([300., 400., 500., 600., 800., 1000., 1500.])
        cp_data = np.zeros_like(Tdata)
        for i in range(Tdata.shape[0]):
            cp_data[i] = self.wilhoit.getHeatCapacity(Tdata[i])
        cp_0 = self.Cp0 * constants.R
        cp_inf = self.CpInf * constants.R

        # Fit the Wilhoit polynomial to the data
        wilhoit = Wilhoit().fitToData(Tdata, cp_data, cp_0, cp_inf, h298, s298)

        # Check that the fit reproduces the input data
        for T in Tdata:
            cp_exp = self.wilhoit.getHeatCapacity(T)
            cp_act = wilhoit.getHeatCapacity(T)
            self.assertAlmostEqual(cp_act, cp_exp, 4)
            h_exp = self.wilhoit.getEnthalpy(T)
            h_act = wilhoit.getEnthalpy(T)
            self.assertAlmostEqual(h_act, h_exp, 3)
            s_exp = self.wilhoit.getEntropy(T)
            s_act = wilhoit.getEntropy(T)
            self.assertAlmostEqual(s_act, s_exp, 4)

        # Check that the fit reproduces the input parameters 
        # Since we're fitting to data generated from a Wilhoit (and since the
        # fitting algorithm is linear least-squares), we should get the same
        # Wilhoit parameters (with a small allowance for fitting error)
        self.assertAlmostEqual(wilhoit.Cp0.value_si, self.wilhoit.Cp0.value_si, 6)
        self.assertAlmostEqual(wilhoit.CpInf.value_si, self.wilhoit.CpInf.value_si, 6)
        self.assertAlmostEqual(wilhoit.a0, self.wilhoit.a0, 2)
        self.assertAlmostEqual(wilhoit.a1, self.wilhoit.a1, 2)
        self.assertAlmostEqual(wilhoit.a2, self.wilhoit.a2, 2)
        self.assertAlmostEqual(wilhoit.a3, self.wilhoit.a3, 2)
        self.assertAlmostEqual(wilhoit.B.value_si, self.wilhoit.B.value_si, 2)
        self.assertAlmostEqual(wilhoit.H0.value_si, self.wilhoit.H0.value_si, 0)
        self.assertAlmostEqual(wilhoit.S0.value_si, self.wilhoit.S0.value_si, 2)

    def testToWilhoit(self):
        """
        Test if the entropy computed from other thermo implementations is close to what Wilhoit computes.
        """

        from rmgpy import settings
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.species import Species

        # Load databases
        database = RMGDatabase()
        database.load_thermo(os.path.join(settings['database.directory'], 'thermo'), thermoLibraries=['Narayanaswamy'])
        database.load_solvation(os.path.join(settings['database.directory'], 'solvation'))

        spc = Species().from_smiles('CC')
        spc.get_thermo_data()

        T = 1350.  # not 298K!

        # nasa to wilhoit
        nasa = spc.thermo
        s_nasa = nasa.getEntropy(T)

        nasa_to_wh = nasa.toWilhoit()
        s_nasa_to_wh = nasa_to_wh.getEntropy(T)

        self.assertAlmostEqual(s_nasa, s_nasa_to_wh, -1)
        self.assertEqual(nasa.comment, nasa_to_wh.comment)

        # wilhoit to nasa conversion done in nasaTest.py

        # thermo data to wilhoit:
        td = nasa.toThermoData()
        s_td = td.getEntropy(T)

        wilhoit = td.toWilhoit(B=1000.)
        s_wh = wilhoit.getEntropy(T)

        self.assertAlmostEqual(s_td, s_wh, -1)
        self.assertEqual(td.comment, wilhoit.comment)

        # wilhoit back to thermodata
        td = wilhoit.toThermoData()
        s_td = td.getEntropy(T)

        self.assertAlmostEqual(s_td, s_wh, -1)
        self.assertEqual(td.comment, wilhoit.comment)

    def testWilhoitAsDict(self):
        """
        Test that a Wilhoit object can be converted to a dictionary representation properly
        """
        wilhoit_dict = self.wilhoit.as_dict()
        self.assertEqual(wilhoit_dict, {'comment': 'C2H6',
                                        'B': {'units': 'K', 'class': 'ScalarQuantity', 'value': 1068.68},
                                        'Tmin': {'units': 'K', 'class': 'ScalarQuantity', 'value': 300.0},
                                        'H0': {'units': 'kJ/mol', 'class': 'ScalarQuantity', 'value': -782.292041536},
                                        'Tmax': {'units': 'K', 'class': 'ScalarQuantity', 'value': 3000.0},
                                        'S0': {'units': 'J/(mol*K)', 'class': 'ScalarQuantity', 'value': -984.93235312},
                                        'a1': -16.3067,
                                        'a0': 0.0977518,
                                        'a3': -12.6785,
                                        'a2': 26.2524,
                                        'Cp0': {'units': 'J/(mol*K)', 'class': 'ScalarQuantity', 'value': 33.257888},
                                        'CpInf': {'units': 'J/(mol*K)', 'class': 'ScalarQuantity',
                                                  'value': 178.76114800000002},
                                        'class': 'Wilhoit'}
                         )

    def testMakeWilhoit(self):
        """
        Test that a Wilhoit object can be created from a dictionary representation
        """
        wilhoit_dict = self.wilhoit.as_dict()
        new_wilhoit = Wilhoit.__new__(Wilhoit)
        class_dictionary = {'ScalarQuantity': ScalarQuantity,
                            'Wilhoit': Wilhoit}

        new_wilhoit.make_object(wilhoit_dict, class_dictionary)


################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
