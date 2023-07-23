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
This script contains unit tests of the :mod:`rmgpy.thermo` object conversion
methods.
"""


import numpy as np

import rmgpy.constants as constants
from rmgpy.thermo import Wilhoit, NASA, NASAPolynomial, ThermoData


class TestConverter:
    """
    Contains unit tests of the thermodynamics model conversion functions.
    """

    def setup_class(self):
        self.wilhoit = Wilhoit(
            Cp0=(4.0 * constants.R, "J/(mol*K)"),
            CpInf=(21.5 * constants.R, "J/(mol*K)"),
            a0=0.0977518,
            a1=-16.3067,
            a2=26.2524,
            a3=-12.6785,
            B=(1068.68, "K"),
            H0=(-94.088 * constants.R, "kJ/mol"),
            S0=(-118.46 * constants.R, "J/(mol*K)"),
            Tmin=(10, "K"),
            Tmax=(3000, "K"),
            comment="C2H6",
        )
        self.nasa = NASA(
            polynomials=[
                NASAPolynomial(
                    coeffs=[
                        4.03055,
                        -0.00214171,
                        4.90611e-05,
                        -5.99027e-08,
                        2.38945e-11,
                        -11257.6,
                        3.5613,
                    ],
                    Tmin=(10, "K"),
                    Tmax=(650.73, "K"),
                ),
                NASAPolynomial(
                    coeffs=[
                        -0.307954,
                        0.0245269,
                        -1.2413e-05,
                        3.07724e-09,
                        -3.01467e-13,
                        -10693,
                        22.628,
                    ],
                    Tmin=(650.73, "K"),
                    Tmax=(3000, "K"),
                ),
            ],
            Tmin=(10, "K"),
            Tmax=(3000, "K"),
            E0=(-93.6077, "kJ/mol"),
            Cp0=(4.0 * constants.R, "J/(mol*K)"),
            CpInf=(21.5 * constants.R, "J/(mol*K)"),
            comment="C2H6",
        )
        self.thermodata = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=(
                np.array([6.38268, 7.80327, 9.22175, 10.5528, 12.8323, 14.6013, 17.40890]) * constants.R,
                "J/(mol*K)",
            ),
            H298=(-81.7, "kJ/mol"),
            S298=(27.5727 * constants.R, "J/(mol*K)"),
            Cp0=(4.0 * constants.R, "J/(mol*K)"),
            CpInf=(21.5 * constants.R, "J/(mol*K)"),
            Tmin=(10, "K"),
            Tmax=(3000, "K"),
            E0=(-93.6077, "kJ/mol"),
            comment="C2H6",
        )

    def test_convert_wilhoit_to_nasa(self):
        """
        Test the conversion of a Wilhoit model to a NASA model.
        """
        wilhoit = self.wilhoit
        nasa = wilhoit.to_nasa(Tmin=10, Tmax=3000, Tint=1000)
        Tlist = np.arange(10, 3000, 10)
        for T in Tlist:
            cp_wilhoit = wilhoit.get_heat_capacity(T)
            cp_nasa = nasa.get_heat_capacity(T)
            assert abs(cp_nasa - cp_wilhoit) < 1e0
            h_wilhoit = wilhoit.get_enthalpy(T)
            h_nasa = nasa.get_enthalpy(T)
            assert abs(h_nasa - h_wilhoit) < 2e1
            s_wilhoit = wilhoit.get_entropy(T)
            s_nasa = nasa.get_entropy(T)
            assert abs(s_nasa - s_wilhoit) < 1e0
        assert abs(wilhoit.E0.value_si - nasa.E0.value_si) < 1e1

    def test_convert_wilhoit_to_thermo_data(self):
        """
        Test the conversion of a Wilhoit model to a ThermoData model.
        """
        wilhoit = self.wilhoit
        thermodata = wilhoit.to_thermo_data()
        Tlist = np.array([300, 400, 500, 600, 800, 1000, 1500])
        for T in Tlist:
            cp_wilhoit = wilhoit.get_heat_capacity(T)
            cp_thermodata = thermodata.get_heat_capacity(T)
            assert round(abs(cp_thermodata - cp_wilhoit), 4) == 0
        T = 298
        h_wilhoit = wilhoit.get_enthalpy(T)
        h_thermodata = thermodata.get_enthalpy(T)
        assert round(abs(h_thermodata - h_wilhoit), 4) == 0
        s_wilhoit = wilhoit.get_entropy(T)
        s_thermodata = thermodata.get_entropy(T)
        assert round(abs(s_thermodata - s_wilhoit), 4) == 0
        assert abs(wilhoit.E0.value_si - thermodata.E0.value_si) < 1e1

    def test_convert_nasa_to_wilhoit(self):
        """
        Test the conversion of a NASA model to a Wilhoit model.
        """
        nasa = self.nasa
        wilhoit = nasa.to_wilhoit()
        Tlist = np.arange(10, 3000, 10)
        for T in Tlist:
            cp_wilhoit = wilhoit.get_heat_capacity(T)
            cp_nasa = nasa.get_heat_capacity(T)
            assert abs(cp_nasa - cp_wilhoit) < 1e0
            h_wilhoit = wilhoit.get_enthalpy(T)
            h_nasa = nasa.get_enthalpy(T)
            assert abs(h_nasa - h_wilhoit) < 2e1
            s_wilhoit = wilhoit.get_entropy(T)
            s_nasa = nasa.get_entropy(T)
            assert abs(s_nasa - s_wilhoit) < 1e0
        assert abs(nasa.E0.value_si - wilhoit.E0.value_si) < 2e1

    def test_convert_nasa_to_thermo_data(self):
        """
        Test the conversion of a NASA model to a ThermoData model.
        """
        nasa = self.nasa
        thermodata = nasa.to_thermo_data()
        Tlist = np.array([300, 400, 500, 600, 800, 1000, 1500])
        for T in Tlist:
            cp_thermodata = thermodata.get_heat_capacity(T)
            cp_nasa = nasa.get_heat_capacity(T)
            assert round(abs(cp_nasa - cp_thermodata), 4) == 0
        T = 298
        h_thermodata = thermodata.get_enthalpy(T)
        h_nasa = nasa.get_enthalpy(T)
        assert round(abs(h_nasa - h_thermodata), 4) == 0
        s_thermodata = thermodata.get_entropy(T)
        s_nasa = nasa.get_entropy(T)
        assert round(abs(s_nasa - s_thermodata), 4) == 0
        assert abs(nasa.E0.value_si - thermodata.E0.value_si) < 1e1

    def test_convert_thermo_data_to_wilhoit(self):
        """
        Test the conversion of a ThermoData model to a Wilhoit model.
        """
        thermodata = self.thermodata
        wilhoit = thermodata.to_wilhoit()
        Tlist = np.array([300, 400, 500, 600, 800, 1000, 1500])
        for T in Tlist:
            cp_wilhoit = wilhoit.get_heat_capacity(T)
            cp_thermodata = thermodata.get_heat_capacity(T)
            assert round(abs(cp_thermodata - cp_wilhoit), 3) == 0
        T = 298
        h_wilhoit = wilhoit.get_enthalpy(T)
        h_thermodata = thermodata.get_enthalpy(T)
        assert round(abs(h_thermodata - h_wilhoit), 3) == 0
        s_wilhoit = wilhoit.get_entropy(T)
        s_thermodata = thermodata.get_entropy(T)
        assert round(abs(s_thermodata - s_wilhoit), 3) == 0
        assert abs(thermodata.E0.value_si - wilhoit.E0.value_si) < 1e1

    def test_convert_thermo_data_to_nasa(self):
        """
        Test the conversion of a ThermoData model to a NASA model.
        """
        thermodata = self.thermodata
        nasa = thermodata.to_nasa(Tmin=10, Tmax=3000, Tint=1000)
        Tlist = np.array([300, 400, 500, 600, 800, 1000, 1500])
        for T in Tlist:
            cp_thermodata = thermodata.get_heat_capacity(T)
            cp_nasa = nasa.get_heat_capacity(T)
            assert abs(cp_nasa - cp_thermodata) < 1e0
        T = 298
        h_thermodata = thermodata.get_enthalpy(T)
        h_nasa = nasa.get_enthalpy(T)
        assert abs(h_nasa - h_thermodata) < 1e0
        s_thermodata = thermodata.get_entropy(T)
        s_nasa = nasa.get_entropy(T)
        assert abs(s_nasa - s_thermodata) < 1e0
        assert abs(thermodata.E0.value_si - nasa.E0.value_si) < 1e1

    def test_wilhoit_nasa_wilhoit(self):
        """
        Test round-trip conversion from Wilhoit to NASA and back
        """
        wilhoit1 = self.wilhoit
        nasa = wilhoit1.to_nasa(Tmin=10, Tmax=3000, Tint=1000)
        wilhoit2 = nasa.to_wilhoit()
        Tlist = np.arange(10, 3000, 10)
        for T in Tlist:
            cp_1 = wilhoit1.get_heat_capacity(T)
            cp_2 = wilhoit2.get_heat_capacity(T)
            assert abs(cp_1 - cp_2) < 1e0
            h_1 = wilhoit1.get_enthalpy(T)
            h_2 = wilhoit2.get_enthalpy(T)
            assert abs(h_1 - h_2) < 2e1
            s_1 = wilhoit1.get_entropy(T)
            s_2 = wilhoit2.get_entropy(T)
            assert abs(s_1 - s_2) < 1e0
        assert abs(wilhoit1.E0.value_si - wilhoit2.E0.value_si) < 1e1

    def test_wilhoit_thermo_data_wilhoit(self):
        """
        Test round-trip conversion from Wilhoit to ThermoData and back
        """
        wilhoit1 = self.wilhoit
        thermodata = wilhoit1.to_thermo_data()
        wilhoit2 = thermodata.to_wilhoit()
        Tlist = np.arange(10, 3000, 10)
        for T in Tlist:
            cp_1 = wilhoit1.get_heat_capacity(T)
            cp_2 = wilhoit2.get_heat_capacity(T)
            assert abs(cp_1 - cp_2) < 1e0
            h_1 = wilhoit1.get_enthalpy(T)
            h_2 = wilhoit2.get_enthalpy(T)
            assert abs(h_1 - h_2) < 2e1
            s_1 = wilhoit1.get_entropy(T)
            s_2 = wilhoit2.get_entropy(T)
            assert abs(s_1 - s_2) < 1e0
        assert abs(wilhoit1.E0.value_si - wilhoit2.E0.value_si) < 1e1
