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
This script contains unit tests of the :mod:`rmgpy.thermo.wilhoit` module.
"""

import os.path


import numpy as np

import rmgpy.constants as constants
from rmgpy.quantity import ScalarQuantity
from rmgpy.thermo.wilhoit import Wilhoit


class TestWilhoit:
    """
    Contains unit tests of the :class:`Wilhoit` class.
    """

    def setup_class(self):
        self.Cp0 = 4.0
        self.CpInf = 21.5
        self.a0 = 0.0977518
        self.a1 = -16.3067
        self.a2 = 26.2524
        self.a3 = -12.6785
        self.B = 1068.68
        self.H0 = -94088.0  # -782.292 kJ/mol / constants.R
        self.S0 = -118.46  # -984.932 J/mol*K / constants.R
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.comment = "C2H6"
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

    def teardown_class(self):
        """
        Reset the database & liquid parameters for solution
        """
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

    def test_cp0(self):
        """
        Test that the Wilhoit Cp0 property was properly set.
        """
        assert round(abs(self.wilhoit.Cp0.value_si / constants.R - self.Cp0), 6) == 0

    def test_cp_inf(self):
        """
        Test that the Wilhoit CpInf property was properly set.
        """
        assert round(abs(self.wilhoit.CpInf.value_si / constants.R - self.CpInf), 6) == 0

    def test_a0(self):
        """
        Test that the Wilhoit a0 property was properly set.
        """
        assert round(abs(self.wilhoit.a0 - self.a0), 6) == 0

    def test_a1(self):
        """
        Test that the Wilhoit a1 property was properly set.
        """
        assert round(abs(self.wilhoit.a1 - self.a1), 6) == 0

    def test_a2(self):
        """
        Test that the Wilhoit a2 property was properly set.
        """
        assert round(abs(self.wilhoit.a2 - self.a2), 6) == 0

    def test_a3(self):
        """
        Test that the Wilhoit a3 property was properly set.
        """
        assert round(abs(self.wilhoit.a3 - self.a3), 6) == 0

    def test_b(self):
        """
        Test that the Wilhoit B property was properly set.
        """
        assert round(abs(self.wilhoit.B.value_si - self.B), 6) == 0

    def test_h0(self):
        """
        Test that the Wilhoit H0 property was properly set.
        """
        assert round(abs(self.wilhoit.H0.value_si / constants.R - self.H0), 6) == 0

    def test_s0(self):
        """
        Test that the Wilhoit S0 property was properly set.
        """
        assert round(abs(self.wilhoit.S0.value_si / constants.R - self.S0), 6) == 0

    def test_temperature_min(self):
        """
        Test that the Wilhoit Tmin property was properly set.
        """
        assert round(abs(self.wilhoit.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the Wilhoit Tmax property was properly set.
        """
        assert round(abs(self.wilhoit.Tmax.value_si - self.Tmax), 6) == 0

    def test_e0(self):
        """
        Test that the Wilhoit E0 property is properly calculated from Enthalpy at 0.001 K
        """
        assert round(abs(self.wilhoit.E0.value_si - self.wilhoit.get_enthalpy(0.001)), 1) == 0

    def test_comment(self):
        """
        Test that the Wilhoit comment property was properly set.
        """
        assert self.wilhoit.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the Wilhoit.is_temperature_valid() method.
        """
        Tdata = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
        valid_data = [False, True, True, True, True, True, True, True, True, True]
        for T, valid in zip(Tdata, valid_data):
            valid0 = self.wilhoit.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_heat_capacity(self):
        """
        Test the Wilhoit.get_heat_capacity() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        cp_exp_list = (
            np.array(
                [
                    5.12003,
                    7.80327,
                    10.5528,
                    12.8323,
                    14.6013,
                    15.9511,
                    16.9842,
                    17.7837,
                    18.4114,
                    18.9117,
                ]
            )
            * constants.R
        )
        for T, cp_exp in zip(Tlist, cp_exp_list):
            cp_act = self.wilhoit.get_heat_capacity(T)
            assert round(abs(cp_exp / cp_act - 1.0), 3) == 0, "{0} != {1} within 3 places".format(cp_exp, cp_act)

    def test_get_enthalpy(self):
        """
        Test the Wilhoit.get_enthalpy() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        h_exp_list = (
            np.array(
                [
                    -51.9303,
                    -22.7609,
                    -12.1050,
                    -6.14444,
                    -2.16433,
                    0.747500,
                    2.99646,
                    4.79698,
                    6.27618,
                    7.51564,
                ]
            )
            * constants.R
            * Tlist
        )
        for T, h_exp in zip(Tlist, h_exp_list):
            h_act = self.wilhoit.get_enthalpy(T)
            assert round(abs(h_exp / h_act - 1.0), 3) == 0, "{0} != {1}".format(h_exp, h_act)

    def test_get_entropy(self):
        """
        Test the Wilhoit.get_entropy() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        s_exp_list = (
            np.array(
                [
                    25.3095,
                    29.6445,
                    33.3398,
                    36.7006,
                    39.7629,
                    42.5499,
                    45.0898,
                    47.4122,
                    49.5445,
                    51.5112,
                ]
            )
            * constants.R
        )
        for T, s_exp in zip(Tlist, s_exp_list):
            s_act = self.wilhoit.get_entropy(T)
            assert round(abs(s_exp / s_act - 1.0), 4) == 0, "{0} != {1}".format(s_exp, s_act)

    def test_get_free_energy(self):
        """
        Test the Wilhoit.get_free_energy() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        for T in Tlist:
            g_exp = self.wilhoit.get_enthalpy(T) - T * self.wilhoit.get_entropy(T)
            g_act = self.wilhoit.get_free_energy(T)
            assert round(abs(g_exp / g_act - 1.0), 4) == 0, "{0} != {1}".format(g_exp, g_act)

    def test_pickle(self):
        """
        Test that a Wilhoit object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        wilhoit = pickle.loads(pickle.dumps(self.wilhoit))
        assert round(abs(self.wilhoit.Cp0.value - wilhoit.Cp0.value), 4) == 0
        assert self.wilhoit.Cp0.units == wilhoit.Cp0.units
        assert round(abs(self.wilhoit.CpInf.value - wilhoit.CpInf.value), 3) == 0
        assert self.wilhoit.CpInf.units == wilhoit.CpInf.units
        assert round(abs(self.wilhoit.a0 - wilhoit.a0), 4) == 0
        assert round(abs(self.wilhoit.a1 - wilhoit.a1), 4) == 0
        assert round(abs(self.wilhoit.a2 - wilhoit.a2), 4) == 0
        assert round(abs(self.wilhoit.a3 - wilhoit.a3), 4) == 0
        assert round(abs(self.wilhoit.B.value - wilhoit.B.value), 4) == 0
        assert self.wilhoit.B.units == wilhoit.B.units
        assert round(abs(self.wilhoit.H0.value - wilhoit.H0.value), 4) == 0
        assert self.wilhoit.H0.units == wilhoit.H0.units
        assert round(abs(self.wilhoit.S0.value - wilhoit.S0.value), 3) == 0
        assert self.wilhoit.S0.units == wilhoit.S0.units
        assert round(abs(self.wilhoit.Tmin.value - wilhoit.Tmin.value), 4) == 0
        assert self.wilhoit.Tmin.units == wilhoit.Tmin.units
        assert round(abs(self.wilhoit.Tmax.value - wilhoit.Tmax.value), 4) == 0
        assert self.wilhoit.Tmax.units == wilhoit.Tmax.units
        assert round(abs(self.wilhoit.E0.value - wilhoit.E0.value), 4) == 0
        assert self.wilhoit.E0.units == wilhoit.E0.units
        assert self.wilhoit.comment == wilhoit.comment

    def test_repr(self):
        """
        Test that a Wilhoit object can be reconstructed from its repr() output
        with no loss of information.
        """
        namespace = {}
        exec("wilhoit = {0!r}".format(self.wilhoit), globals(), namespace)
        assert "wilhoit" in namespace
        wilhoit = namespace["wilhoit"]
        assert round(abs(self.wilhoit.Cp0.value - wilhoit.Cp0.value), 4) == 0
        assert self.wilhoit.Cp0.units == wilhoit.Cp0.units
        assert round(abs(self.wilhoit.CpInf.value - wilhoit.CpInf.value), 3) == 0
        assert self.wilhoit.CpInf.units == wilhoit.CpInf.units
        assert round(abs(self.wilhoit.a0 - wilhoit.a0), 4) == 0
        assert round(abs(self.wilhoit.a1 - wilhoit.a1), 4) == 0
        assert round(abs(self.wilhoit.a2 - wilhoit.a2), 4) == 0
        assert round(abs(self.wilhoit.a3 - wilhoit.a3), 4) == 0
        assert round(abs(self.wilhoit.B.value - wilhoit.B.value), 4) == 0
        assert self.wilhoit.B.units == wilhoit.B.units
        assert round(abs(self.wilhoit.H0.value - wilhoit.H0.value), 4) == 0
        assert self.wilhoit.H0.units == wilhoit.H0.units
        assert round(abs(self.wilhoit.S0.value - wilhoit.S0.value), 3) == 0
        assert self.wilhoit.S0.units == wilhoit.S0.units
        assert round(abs(self.wilhoit.Tmin.value - wilhoit.Tmin.value), 4) == 0
        assert self.wilhoit.Tmin.units == wilhoit.Tmin.units
        assert round(abs(self.wilhoit.Tmax.value - wilhoit.Tmax.value), 4) == 0
        assert self.wilhoit.Tmax.units == wilhoit.Tmax.units
        assert round(abs(self.wilhoit.E0.value - wilhoit.E0.value), 1) == 0
        assert self.wilhoit.E0.units == wilhoit.E0.units
        assert self.wilhoit.comment == wilhoit.comment

    def test_fit_to_data(self):
        """
        Test the Wilhoit.fit_to_data() method.
        """
        h298 = self.wilhoit.get_enthalpy(298)
        s298 = self.wilhoit.get_entropy(298)
        Tdata = np.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0])
        cp_data = np.zeros_like(Tdata)
        for i in range(Tdata.shape[0]):
            cp_data[i] = self.wilhoit.get_heat_capacity(Tdata[i])
        cp_0 = self.Cp0 * constants.R
        cp_inf = self.CpInf * constants.R

        # Fit the Wilhoit polynomial to the data
        wilhoit = Wilhoit().fit_to_data(Tdata, cp_data, cp_0, cp_inf, h298, s298)

        # Check that the fit reproduces the input data
        for T in Tdata:
            cp_exp = self.wilhoit.get_heat_capacity(T)
            cp_act = wilhoit.get_heat_capacity(T)
            assert round(abs(cp_act - cp_exp), 4) == 0
            h_exp = self.wilhoit.get_enthalpy(T)
            h_act = wilhoit.get_enthalpy(T)
            assert round(abs(h_act - h_exp), 3) == 0
            s_exp = self.wilhoit.get_entropy(T)
            s_act = wilhoit.get_entropy(T)
            assert round(abs(s_act - s_exp), 4) == 0

        # Check that the fit reproduces the input parameters
        # Since we're fitting to data generated from a Wilhoit (and since the
        # fitting algorithm is linear least-squares), we should get the same
        # Wilhoit parameters (with a small allowance for fitting error)
        assert round(abs(wilhoit.Cp0.value_si - self.wilhoit.Cp0.value_si), 6) == 0
        assert round(abs(wilhoit.CpInf.value_si - self.wilhoit.CpInf.value_si), 6) == 0
        assert round(abs(wilhoit.a0 - self.wilhoit.a0), 2) == 0
        assert round(abs(wilhoit.a1 - self.wilhoit.a1), 2) == 0
        assert round(abs(wilhoit.a2 - self.wilhoit.a2), 2) == 0
        assert round(abs(wilhoit.a3 - self.wilhoit.a3), 2) == 0
        assert round(abs(wilhoit.B.value_si - self.wilhoit.B.value_si), 2) == 0
        assert round(abs(wilhoit.H0.value_si - self.wilhoit.H0.value_si), 0) == 0
        assert round(abs(wilhoit.S0.value_si - self.wilhoit.S0.value_si), 2) == 0

    def test_to_wilhoit(self):
        """
        Test if the entropy computed from other thermo implementations is close to what Wilhoit computes.
        """

        from rmgpy import settings
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.species import Species

        # Load databases
        database = RMGDatabase()
        database.load_thermo(
            os.path.join(settings["database.directory"], "thermo"),
            thermo_libraries=["Narayanaswamy"],
        )
        database.load_solvation(os.path.join(settings["database.directory"], "solvation"))

        spc = Species().from_smiles("CC")
        spc.get_thermo_data()

        T = 1350.0  # not 298K!

        # nasa to wilhoit
        nasa = spc.thermo
        s_nasa = nasa.get_entropy(T)

        nasa_to_wh = nasa.to_wilhoit()
        s_nasa_to_wh = nasa_to_wh.get_entropy(T)

        assert round(abs(s_nasa - s_nasa_to_wh), -1) == 0
        assert nasa.comment == nasa_to_wh.comment

        # wilhoit to nasa conversion done in nasaTest.py

        # thermo data to wilhoit:
        td = nasa.to_thermo_data()
        s_td = td.get_entropy(T)

        wilhoit = td.to_wilhoit(B=1000.0)
        s_wh = wilhoit.get_entropy(T)

        assert round(abs(s_td - s_wh), -1) == 0
        assert td.comment == wilhoit.comment

        # wilhoit back to thermodata
        td = wilhoit.to_thermo_data()
        s_td = td.get_entropy(T)

        assert round(abs(s_td - s_wh), -1) == 0
        assert td.comment == wilhoit.comment

    def test_wilhoit_as_dict(self):
        """
        Test that a Wilhoit object can be converted to a dictionary representation properly
        """
        wilhoit_dict = self.wilhoit.as_dict()
        assert wilhoit_dict == {
            "comment": "C2H6",
            "B": {"units": "K", "class": "ScalarQuantity", "value": 1068.68},
            "Tmin": {"units": "K", "class": "ScalarQuantity", "value": 300.0},
            "H0": {
                "units": "kJ/mol",
                "class": "ScalarQuantity",
                "value": -782.292041536,
            },
            "Tmax": {"units": "K", "class": "ScalarQuantity", "value": 3000.0},
            "S0": {
                "units": "J/(mol*K)",
                "class": "ScalarQuantity",
                "value": -984.93235312,
            },
            "a1": -16.3067,
            "a0": 0.0977518,
            "a3": -12.6785,
            "a2": 26.2524,
            "Cp0": {
                "units": "J/(mol*K)",
                "class": "ScalarQuantity",
                "value": 33.257888,
            },
            "CpInf": {
                "units": "J/(mol*K)",
                "class": "ScalarQuantity",
                "value": 178.76114800000002,
            },
            "class": "Wilhoit",
        }

    def test_make_wilhoit(self):
        """
        Test that a Wilhoit object can be created from a dictionary representation
        """
        wilhoit_dict = self.wilhoit.as_dict()
        new_wilhoit = Wilhoit.__new__(Wilhoit)
        class_dictionary = {"ScalarQuantity": ScalarQuantity, "Wilhoit": Wilhoit}

        new_wilhoit.make_object(wilhoit_dict, class_dictionary)
