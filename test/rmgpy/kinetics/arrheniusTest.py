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

import math


import numpy as np

import rmgpy.constants as constants
from rmgpy.kinetics.arrhenius import (
    Arrhenius,
    ArrheniusEP,
    ArrheniusBM,
    PDepArrhenius,
    MultiArrhenius,
    MultiPDepArrhenius,
)
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
import pytest


class TestArrhenius:
    """
    Contains unit tests of the :class:`Arrhenius` class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.A = 1.0e12
        self.n = 0.5
        self.Ea = 41.84
        self.T0 = 1.0
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.comment = "C2H6"
        self.arrhenius = Arrhenius(
            A=(self.A, "cm^3/(mol*s)"),
            n=self.n,
            Ea=(self.Ea, "kJ/mol"),
            T0=(self.T0, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_a_factor(self):
        """
        Test that the Arrhenius A property was properly set.
        """
        assert abs(self.arrhenius.A.value_si * 1e6 - self.A) < 1e0

    def test_n(self):
        """
        Test that the Arrhenius n property was properly set.
        """
        assert round(abs(self.arrhenius.n.value_si - self.n), 6) == 0

    def test_ea(self):
        """
        Test that the Arrhenius Ea property was properly set.
        """
        assert round(abs(self.arrhenius.Ea.value_si * 0.001 - self.Ea), 6) == 0

    def test_temperature0(self):
        """
        Test that the Arrhenius T0 property was properly set.
        """
        assert round(abs(self.arrhenius.T0.value_si - self.T0), 6) == 0

    def test_temperature_min(self):
        """
        Test that the Arrhenius Tmin property was properly set.
        """
        assert round(abs(self.arrhenius.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the Arrhenius Tmax property was properly set.
        """
        assert round(abs(self.arrhenius.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the Arrhenius comment property was properly set.
        """
        assert self.arrhenius.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the Arrhenius.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], np.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.arrhenius.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the Arrhenius.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        kexplist = np.array(
            [
                1.6721e-4,
                6.8770e1,
                5.5803e3,
                5.2448e4,
                2.0632e5,
                5.2285e5,
                1.0281e6,
                1.7225e6,
                2.5912e6,
                3.6123e6,
            ]
        )
        for T, kexp in zip(Tlist, kexplist):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(kexp - kact) < 1e-4 * kexp

    def test_change_t0(self):
        """
        Test the Arrhenius.change_t0() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tlist])
        self.arrhenius.change_t0(300)
        assert self.arrhenius.T0.value_si == 300
        for T, kexp in zip(Tlist, k0list):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(kexp - kact) < 1e-6 * kexp

    def test_fit_to_data(self):
        """
        Test the Arrhenius.fit_to_data() method.
        """
        Tdata = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        kdata = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tdata])
        arrhenius = Arrhenius().fit_to_data(Tdata, kdata, kunits="m^3/(mol*s)")
        assert float(self.arrhenius.T0.value_si) == 1
        for T, k in zip(Tdata, kdata):
            assert abs(k - arrhenius.get_rate_coefficient(T)) < 1e-6 * k
        assert abs(arrhenius.A.value_si - self.arrhenius.A.value_si) < 1e0
        assert round(abs(arrhenius.n.value_si - self.arrhenius.n.value_si), 1) == 0, 4
        assert round(abs(arrhenius.Ea.value_si - self.arrhenius.Ea.value_si), 2) == 0
        assert round(abs(arrhenius.T0.value_si - self.arrhenius.T0.value_si), 4) == 0

    def test_fit_to_negative_data(self):
        """
        Test the Arrhenius.fit_to_data() method on negative rates
        """
        Tdata = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        kdata = np.array([-1 * self.arrhenius.get_rate_coefficient(T) for T in Tdata])
        arrhenius = Arrhenius().fit_to_data(Tdata, kdata, kunits="m^3/(mol*s)")
        assert float(self.arrhenius.T0.value_si) == 1
        for T, k in zip(Tdata, kdata):
            assert abs(k - arrhenius.get_rate_coefficient(T)) < 1e-6 * abs(k)
        assert abs(arrhenius.A.value_si - -1 * self.arrhenius.A.value_si) < 1e0
        assert round(abs(arrhenius.n.value_si - self.arrhenius.n.value_si), 1) == 0, 4
        assert round(abs(arrhenius.Ea.value_si - self.arrhenius.Ea.value_si), 2) == 0
        assert round(abs(arrhenius.T0.value_si - self.arrhenius.T0.value_si), 4) == 0

    def test_pickle(self):
        """
        Test that an Arrhenius object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        arrhenius = pickle.loads(pickle.dumps(self.arrhenius, -1))
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.Ea.value - arrhenius.Ea.value), 4) == 0
        assert self.arrhenius.Ea.units == arrhenius.Ea.units
        assert round(abs(self.arrhenius.T0.value - arrhenius.T0.value), 4) == 0
        assert self.arrhenius.T0.units == arrhenius.T0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_repr(self):
        """
        Test that an Arrhenius object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("arrhenius = {0!r}".format(self.arrhenius), globals(), namespace)
        assert "arrhenius" in namespace
        arrhenius = namespace["arrhenius"]
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.Ea.value - arrhenius.Ea.value), 4) == 0
        assert self.arrhenius.Ea.units == arrhenius.Ea.units
        assert round(abs(self.arrhenius.T0.value - arrhenius.T0.value), 4) == 0
        assert self.arrhenius.T0.units == arrhenius.T0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_change_rate(self):
        """
        Test the Arrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tlist])
        self.arrhenius.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(2 * kexp - kact) < 1e-6 * kexp

    def test_to_cantera_kinetics(self):
        """
        Test that the Arrhenius cantera object can be set properly within
        a cantera Reaction object
        """
        ctArrhenius = self.arrhenius.to_cantera_kinetics()
        assert round(abs(ctArrhenius.pre_exponential_factor - 1e9), 6) == 0
        assert round(abs(ctArrhenius.temperature_exponent - 0.5), 7) == 0
        assert round(abs(ctArrhenius.activation_energy - 41.84e6), 7) == 0

    def test_to_arrhenius_ep(self):
        """
        Tests that the Arrhenius object can be converted to ArrheniusEP
        """
        arr_rate = self.arrhenius.get_rate_coefficient(500)
        arr_ep = self.arrhenius.to_arrhenius_ep()
        arr_ep_rate = arr_ep.get_rate_coefficient(500, 10)  # the second number should not matter
        assert round(abs(arr_rate - arr_ep_rate), 7) == 0

    def test_to_arrhenius_ep_with_alpha_and_hrxn(self):
        """
        Tests that the Arrhenius object can be converted to ArrheniusEP given parameters
        """
        hrxn = 5
        arr_rate = self.arrhenius.get_rate_coefficient(500)
        arr_ep = self.arrhenius.to_arrhenius_ep(alpha=1, dHrxn=hrxn)
        assert round(abs(1.0 - arr_ep.alpha.value_si), 7) == 0
        arr_ep_rate = arr_ep.get_rate_coefficient(500, hrxn)
        assert round(abs(arr_rate - arr_ep_rate), 7) == 0

    def test_to_arrhenius_ep_throws_error_with_just_alpha(self):
        with pytest.raises(Exception):
            self.arrhenius.to_arrhenius_ep(alpha=1)


class TestArrheniusEP:
    """
    Contains unit tests of the :class:`ArrheniusEP` class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.A = 1.0e12
        self.n = 0.5
        self.alpha = 0.5
        self.E0 = 41.84
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.comment = "C2H6"
        self.arrhenius = ArrheniusEP(
            A=(self.A, "cm^3/(mol*s)"),
            n=self.n,
            alpha=self.alpha,
            E0=(self.E0, "kJ/mol"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_a_factor(self):
        """
        Test that the ArrheniusEP A property was properly set.
        """
        assert abs(self.arrhenius.A.value_si * 1e6 - self.A) < 1e0

    def test_n(self):
        """
        Test that the ArrheniusEP n property was properly set.
        """
        assert round(abs(self.arrhenius.n.value_si - self.n), 6) == 0

    def test_alpha(self):
        """
        Test that the ArrheniusEP alpha property was properly set.
        """
        assert round(abs(self.arrhenius.alpha.value_si - self.alpha), 6) == 0

    def test_e0(self):
        """
        Test that the ArrheniusEP E0 property was properly set.
        """
        assert round(abs(self.arrhenius.E0.value_si * 0.001 - self.E0), 6) == 0

    def test_temperature_min(self):
        """
        Test that the ArrheniusEP Tmin property was properly set.
        """
        assert round(abs(self.arrhenius.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the ArrheniusEP Tmax property was properly set.
        """
        assert round(abs(self.arrhenius.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the ArrheniusEP comment property was properly set.
        """
        assert self.arrhenius.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the ArrheniusEP.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], np.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.arrhenius.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the ArrheniusEP.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        kexplist = np.array(
            [
                1.6721e-4,
                6.8770e1,
                5.5803e3,
                5.2448e4,
                2.0632e5,
                5.2285e5,
                1.0281e6,
                1.7225e6,
                2.5912e6,
                3.6123e6,
            ]
        )
        for T, kexp in zip(Tlist, kexplist):
            kact = self.arrhenius.get_rate_coefficient(
                T,
            )
            assert abs(kexp - kact) < 1e-4 * kexp

    def test_pickle(self):
        """
        Test that an ArrheniusEP object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        arrhenius = pickle.loads(pickle.dumps(self.arrhenius, -1))
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.alpha.value - arrhenius.alpha.value), 4) == 0
        assert round(abs(self.arrhenius.E0.value - arrhenius.E0.value), 4) == 0
        assert self.arrhenius.E0.units == arrhenius.E0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_repr(self):
        """
        Test that an ArrheniusEP object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("arrhenius = {0!r}".format(self.arrhenius), globals(), namespace)
        assert "arrhenius" in namespace
        arrhenius = namespace["arrhenius"]
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.alpha.value - arrhenius.alpha.value), 4) == 0
        assert round(abs(self.arrhenius.E0.value - arrhenius.E0.value), 4) == 0
        assert self.arrhenius.E0.units == arrhenius.E0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_change_rate(self):
        """
        Test the ArrheniusEP.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tlist])
        self.arrhenius.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(2 * kexp - kact) < 1e-6 * kexp


class TestArrheniusBM:
    """
    Contains unit tests of the :class:`ArrheniusBM` class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.A = 8.00037e12
        self.n = 0.391734
        self.w0 = 798000
        self.E0 = 115905
        self.Tmin = 300.0
        self.Tmax = 2000.0
        self.comment = "rxn001084"
        self.arrhenius_bm = ArrheniusBM(
            A=(self.A, "s^-1"),
            n=self.n,
            w0=(self.w0, "J/mol"),
            E0=(self.E0, "J/mol"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

        self.rsmi = "NC(=NC=O)O"
        self.psmi = "O=CNC(=O)N"
        self.arrhenius = Arrhenius(
            A=(8.00037e12, "s^-1"),
            n=0.391734,
            Ea=(94.5149, "kJ/mol"),
            T0=(1, "K"),
            Tmin=(300, "K"),
            Tmax=(2000, "K"),
            comment="""Fitted to 50 data points; dA = *|/ 1.18377, dn = +|- 0.0223855, dEa = +|- 0.115431 kJ/mol""",
        )

        self.r_thermo = NASA(
            polynomials=[
                NASAPolynomial(
                    coeffs=[
                        3.90453,
                        0.0068491,
                        0.000125755,
                        -2.92973e-07,
                        2.12971e-10,
                        -45444.2,
                        10.0669,
                    ],
                    Tmin=(10, "K"),
                    Tmax=(433.425, "K"),
                ),
                NASAPolynomial(
                    coeffs=[
                        2.09778,
                        0.0367646,
                        -2.36023e-05,
                        7.24527e-09,
                        -8.51275e-13,
                        -45412,
                        15.8381,
                    ],
                    Tmin=(433.425, "K"),
                    Tmax=(3000, "K"),
                ),
            ],
            Tmin=(10, "K"),
            Tmax=(3000, "K"),
            E0=(-377.851, "kJ/mol"),
            Cp0=(33.2579, "J/(mol*K)"),
            CpInf=(232.805, "J/(mol*K)"),
            comment="""Thermo library: Spiekermann_refining_elementary_reactions""",
        )
        self.p_thermo = NASA(
            polynomials=[
                NASAPolynomial(
                    coeffs=[
                        3.88423,
                        0.00825528,
                        0.000133399,
                        -3.31802e-07,
                        2.52823e-10,
                        -51045.1,
                        10.3937,
                    ],
                    Tmin=(10, "K"),
                    Tmax=(428.701, "K"),
                ),
                NASAPolynomial(
                    coeffs=[
                        2.89294,
                        0.0351772,
                        -2.26349e-05,
                        7.00331e-09,
                        -8.2982e-13,
                        -51122.5,
                        12.4424,
                    ],
                    Tmin=(428.701, "K"),
                    Tmax=(3000, "K"),
                ),
            ],
            Tmin=(10, "K"),
            Tmax=(3000, "K"),
            E0=(-424.419, "kJ/mol"),
            Cp0=(33.2579, "J/(mol*K)"),
            CpInf=(232.805, "J/(mol*K)"),
            comment="""Thermo library: Spiekermann_refining_elementary_reactions""",
        )

    def test_a_factor(self):
        """
        Test that the ArrheniusBM A property was properly set.
        """
        assert abs(self.arrhenius_bm.A.value_si - self.A) < 1e0

    def test_n(self):
        """
        Test that the ArrheniusBM n property was properly set.
        """
        assert round(abs(self.arrhenius_bm.n.value_si - self.n), 6) == 0

    def test_w0(self):
        """
        Test that the ArrheniusBM w0 property was properly set.
        """
        assert round(abs(self.arrhenius_bm.w0.value_si - self.w0), 6) == 0

    def test_e0(self):
        """
        Test that the ArrheniusBM E0 property was properly set.
        """
        assert round(abs(self.arrhenius_bm.E0.value_si - self.E0), 6) == 0

    def test_temperature_min(self):
        """
        Test that the ArrheniusBM Tmin property was properly set.
        """
        assert round(abs(self.arrhenius_bm.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the ArrheniusBM Tmax property was properly set.
        """
        assert round(abs(self.arrhenius_bm.Tmax.value_si - self.Tmax), 6) == 0

    def test_is_temperature_valid(self):
        """
        Test the ArrheniusBM.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], np.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.arrhenius_bm.is_temperature_valid(T)
            assert valid0 == valid

    def test_fit_to_data(self):
        """
        Test the ArrheniusBM.fit_to_data() method.
        """
        reactant = Molecule(smiles=self.rsmi)
        product = Molecule(smiles=self.psmi)
        reaction = Reaction(
            reactants=[
                Species(
                    molecule=[reactant],
                    thermo=self.r_thermo,
                )
            ],
            products=[Species(molecule=[product], thermo=self.p_thermo)],
            kinetics=self.arrhenius,
        )

        arrhenius_bm = ArrheniusBM().fit_to_reactions([reaction], w0=self.w0)
        assert abs(arrhenius_bm.A.value_si - self.arrhenius_bm.A.value_si) < 1.5e1
        assert round(abs(arrhenius_bm.n.value_si - self.arrhenius_bm.n.value_si), 1) == 0, 4
        assert round(abs(arrhenius_bm.E0.value_si - self.arrhenius_bm.E0.value_si), 1) == 0

    def test_get_activation_energy(self):
        """
        Test the ArrheniusBM.get_activation_energy() method.
        """
        Hrxn = -44000  # J/mol
        Ea = self.arrhenius_bm.get_activation_energy(Hrxn)
        assert abs(Ea - 95074) < 1e1


class TestPDepArrhenius:
    """
    Contains unit tests of the :class:`PDepArrhenius` class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.arrhenius0 = Arrhenius(
            A=(1.0e6, "s^-1"),
            n=1.0,
            Ea=(10.0, "kJ/mol"),
            T0=(300.0, "K"),
            Tmin=(300.0, "K"),
            Tmax=(2000.0, "K"),
            comment="""This data is completely made up""",
        )
        self.arrhenius1 = Arrhenius(
            A=(1.0e12, "s^-1"),
            n=1.0,
            Ea=(20.0, "kJ/mol"),
            T0=(300.0, "K"),
            Tmin=(300.0, "K"),
            Tmax=(2000.0, "K"),
            comment="""This data is completely made up""",
        )
        self.pressures = np.array([0.1, 10.0])
        self.arrhenius = [self.arrhenius0, self.arrhenius1]
        self.Tmin = 300.0
        self.Tmax = 2000.0
        self.Pmin = 0.1
        self.Pmax = 10.0
        self.comment = """This data is completely made up"""
        self.kinetics = PDepArrhenius(
            pressures=(self.pressures, "bar"),
            arrhenius=self.arrhenius,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )

    def test_pressures(self):
        """
        Test that the PDepArrhenius pressures property was properly set.
        """
        assert len(self.kinetics.pressures.value_si) == 2
        for i in range(2):
            assert round(abs(self.kinetics.pressures.value_si[i] * 1e-5 - self.pressures[i]), 4) == 0

    def test_arrhenius(self):
        """
        Test that the PDepArrhenius arrhenius property was properly set.
        """
        assert len(self.kinetics.arrhenius) == 2
        for i in range(2):
            assert abs(self.kinetics.arrhenius[i].A.value - self.arrhenius[i].A.value) < 1e0
            assert self.kinetics.arrhenius[i].A.units == self.arrhenius[i].A.units
            assert round(abs(self.kinetics.arrhenius[i].n.value - self.arrhenius[i].n.value), 4) == 0
            assert round(abs(self.kinetics.arrhenius[i].Ea.value - self.arrhenius[i].Ea.value), 4) == 0
            assert self.kinetics.arrhenius[i].Ea.units == self.arrhenius[i].Ea.units
            assert round(abs(self.kinetics.arrhenius[i].T0.value - self.arrhenius[i].T0.value), 4) == 0
            assert self.kinetics.arrhenius[i].T0.units == self.arrhenius[i].T0.units
            assert round(abs(self.kinetics.arrhenius[i].Tmin.value - self.arrhenius[i].Tmin.value), 4) == 0
            assert self.kinetics.arrhenius[i].Tmin.units == self.arrhenius[i].Tmin.units
            assert round(abs(self.kinetics.arrhenius[i].Tmax.value - self.arrhenius[i].Tmax.value), 4) == 0
            assert self.kinetics.arrhenius[i].Tmax.units == self.arrhenius[i].Tmax.units
            assert self.kinetics.arrhenius[i].comment == self.arrhenius[i].comment

    def test_temperature_min(self):
        """
        Test that the PDepArrhenius Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the PDepArrhenius Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_pressure_min(self):
        """
        Test that the PDepArrhenius Pmin property was properly set.
        """
        assert round(abs(self.kinetics.Pmin.value_si * 1e-5 - self.Pmin), 6) == 0

    def test_pressure_max(self):
        """
        Test that the PDepArrhenius Pmax property was properly set.
        """
        assert round(abs(self.kinetics.Pmax.value_si * 1e-5 - self.Pmax), 6) == 0

    def test_comment(self):
        """
        Test that the PDepArrhenius comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_pressure_dependent(self):
        """
        Test the PDepArrhenius.is_pressure_dependent() method.
        """
        assert self.kinetics.is_pressure_dependent()

    def test_get_rate_coefficient(self):
        """
        Test the PDepArrhenius.get_rate_coefficient() method.
        """
        P = 1e4
        for T in [
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
            1500,
        ]:
            k0 = self.kinetics.get_rate_coefficient(T, P)
            k1 = self.arrhenius0.get_rate_coefficient(T)
            assert abs(k0 - k1) < 1e-6 * k1
        P = 1e6
        for T in [
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
            1500,
        ]:
            k0 = self.kinetics.get_rate_coefficient(T, P)
            k1 = self.arrhenius1.get_rate_coefficient(T)
            assert abs(k0 - k1) < 1e-6 * k1
        P = 1e5
        for T in [
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
            1500,
        ]:
            k0 = self.kinetics.get_rate_coefficient(T, P)
            k1 = math.sqrt(self.arrhenius0.get_rate_coefficient(T) * self.arrhenius1.get_rate_coefficient(T))
            assert abs(k0 - k1) < 1e-6 * k1

    def test_fit_to_data(self):
        """
        Test the PDepArrhenius.fit_to_data() method.
        """
        Tdata = np.array(
            [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500],
            np.float,
        )
        Pdata = np.array([1e4, 3e4, 1e5, 3e5, 1e6], np.float)
        kdata = np.zeros([len(Tdata), len(Pdata)], np.float)
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                kdata[t, p] = self.kinetics.get_rate_coefficient(Tdata[t], Pdata[p])
        kinetics = PDepArrhenius().fit_to_data(Tdata, Pdata, kdata, kunits="s^-1")
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                assert abs(kinetics.get_rate_coefficient(Tdata[t], Pdata[p]) - kdata[t, p]) < 1e-6 * kdata[t, p]

    def test_pickle(self):
        """
        Test that a PDepArrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        Narrh = 2
        assert len(self.kinetics.pressures.value) == Narrh
        assert len(kinetics.pressures.value) == Narrh
        assert len(self.kinetics.arrhenius) == Narrh
        assert len(kinetics.arrhenius) == Narrh
        for i in range(Narrh):
            assert round(abs(self.kinetics.pressures.value[i] - kinetics.pressures.value[i]), 4) == 0
            assert abs(self.kinetics.arrhenius[i].A.value - kinetics.arrhenius[i].A.value) < 1e0
            assert self.kinetics.arrhenius[i].A.units == kinetics.arrhenius[i].A.units
            assert round(abs(self.kinetics.arrhenius[i].n.value - kinetics.arrhenius[i].n.value), 7) == 0
            assert round(abs(self.kinetics.arrhenius[i].T0.value - kinetics.arrhenius[i].T0.value), 4) == 0
            assert self.kinetics.arrhenius[i].T0.units == kinetics.arrhenius[i].T0.units
            assert round(abs(self.kinetics.arrhenius[i].Ea.value - kinetics.arrhenius[i].Ea.value), 4) == 0
            assert self.kinetics.arrhenius[i].Ea.units == kinetics.arrhenius[i].Ea.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert round(abs(self.kinetics.Pmin.value - kinetics.Pmin.value), 4) == 0
        assert self.kinetics.Pmin.units == kinetics.Pmin.units
        assert round(abs(self.kinetics.Pmax.value - kinetics.Pmax.value), 4) == 0
        assert self.kinetics.Pmax.units == kinetics.Pmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a PDepArrhenius object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        Narrh = 2
        assert len(self.kinetics.pressures.value) == Narrh
        assert len(kinetics.pressures.value) == Narrh
        assert len(self.kinetics.arrhenius) == Narrh
        assert len(kinetics.arrhenius) == Narrh
        for i in range(Narrh):
            assert round(abs(self.kinetics.pressures.value[i] - kinetics.pressures.value[i]), 4) == 0
            assert abs(self.kinetics.arrhenius[i].A.value - kinetics.arrhenius[i].A.value) < 1e0
            assert self.kinetics.arrhenius[i].A.units == kinetics.arrhenius[i].A.units
            assert round(abs(self.kinetics.arrhenius[i].n.value - kinetics.arrhenius[i].n.value), 7) == 0
            assert round(abs(self.kinetics.arrhenius[i].T0.value - kinetics.arrhenius[i].T0.value), 4) == 0
            assert self.kinetics.arrhenius[i].T0.units == kinetics.arrhenius[i].T0.units
            assert round(abs(self.kinetics.arrhenius[i].Ea.value - kinetics.arrhenius[i].Ea.value), 4) == 0
            assert self.kinetics.arrhenius[i].Ea.units == kinetics.arrhenius[i].Ea.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert round(abs(self.kinetics.Pmin.value - kinetics.Pmin.value), 4) == 0
        assert self.kinetics.Pmin.units == kinetics.Pmin.units
        assert round(abs(self.kinetics.Pmax.value - kinetics.Pmax.value), 4) == 0
        assert self.kinetics.Pmax.units == kinetics.Pmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_change_rate(self):
        """
        Test the PDepArrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.kinetics.get_rate_coefficient(T, 1e5) for T in Tlist])
        self.kinetics.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.kinetics.get_rate_coefficient(T, 1e5)
            assert abs(2 * kexp - kact) < 1e-6 * kexp


class TestMultiArrhenius:
    """
    Contains unit tests of the :class:`MultiArrhenius` class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.Tmin = 350.0
        self.Tmax = 1500.0
        self.comment = "Comment"
        self.arrhenius = [
            Arrhenius(
                A=(9.3e-14, "cm^3/(molecule*s)"),
                n=0.0,
                Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                comment=self.comment,
            ),
            Arrhenius(
                A=(1.4e-9, "cm^3/(molecule*s)"),
                n=0.0,
                Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                comment=self.comment,
            ),
        ]
        self.kinetics = MultiArrhenius(
            arrhenius=self.arrhenius,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )
        self.single_kinetics = MultiArrhenius(
            arrhenius=self.arrhenius[:1],
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_arrhenius(self):
        """
        Test that the MultiArrhenius A property was properly set.
        """
        assert self.kinetics.arrhenius == self.arrhenius

    def test_temperature_min(self):
        """
        Test that the MultiArrhenius Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the MultiArrhenius Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the MultiArrhenius comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the MultiArrhenius.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, False, False, False], np.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the MultiArrhenius.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        kexplist = np.array(
            [
                2.85400e-06,
                4.00384e-01,
                2.73563e01,
                8.50699e02,
                1.20181e04,
                7.56312e04,
                2.84724e05,
                7.71702e05,
                1.67743e06,
                3.12290e06,
            ]
        )
        for T, kexp in zip(Tlist, kexplist):
            kact = self.kinetics.get_rate_coefficient(T)
            assert abs(kexp - kact) < 1e-4 * kexp

    def test_pickle(self):
        """
        Test that a MultiArrhenius object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        for arrh0, arrh in zip(self.kinetics.arrhenius, kinetics.arrhenius):
            assert abs(arrh0.A.value - arrh.A.value) < 1e-18
            assert arrh0.A.units == arrh.A.units
            assert round(abs(arrh0.n.value - arrh.n.value), 4) == 0
            assert round(abs(arrh0.Ea.value - arrh.Ea.value), 4) == 0
            assert arrh0.Ea.units == arrh.Ea.units
            assert round(abs(arrh0.T0.value - arrh.T0.value), 4) == 0
            assert arrh0.T0.units == arrh.T0.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a MultiArrhenius object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        for arrh0, arrh in zip(self.kinetics.arrhenius, kinetics.arrhenius):
            assert abs(arrh0.A.value - arrh.A.value) < 1e-18
            assert arrh0.A.units == arrh.A.units
            assert round(abs(arrh0.n.value - arrh.n.value), 4) == 0
            assert round(abs(arrh0.Ea.value - arrh.Ea.value), 4) == 0
            assert arrh0.Ea.units == arrh.Ea.units
            assert round(abs(arrh0.T0.value - arrh.T0.value), 4) == 0
            assert arrh0.T0.units == arrh.T0.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_to_arrhenius(self):
        """
        Test that we can convert to an Arrhenius
        """
        answer = self.single_kinetics.arrhenius[0]
        fitted = self.single_kinetics.to_arrhenius()

        assert abs(fitted.A.value_si - answer.A.value_si) < 1e0
        assert round(abs(fitted.n.value_si - answer.n.value_si), 1) == 0, 4
        assert round(abs(fitted.Ea.value_si - answer.Ea.value_si), 2) == 0
        assert round(abs(fitted.T0.value_si - answer.T0.value_si), 4) == 0

    def test_to_arrhenius_temperature_range(self):
        """
        Test the to_arrhenius temperature range is set correctly.
        """
        answer = self.single_kinetics.arrhenius[0]
        fitted = self.single_kinetics.to_arrhenius(Tmin=800, Tmax=1200)
        assert round(abs(fitted.Tmin.value_si - 800.0), 7) == 0
        assert round(abs(fitted.Tmax.value_si - 1200.0), 7) == 0
        for T in [800, 1000, 1200]:
            assert round(abs(fitted.get_rate_coefficient(T) / answer.get_rate_coefficient(T) - 1.0), 7) == 0

    def test_to_arrhenius_multiple(self):
        """
        Test the to_arrhenius fitting multiple kinetics over a small range, see if we're within 5% at a few points
        """
        answer = self.kinetics
        fitted = self.kinetics.to_arrhenius(Tmin=800, Tmax=1200)
        assert round(abs(fitted.Tmin.value_si - 800.0), 7) == 0
        assert round(abs(fitted.Tmax.value_si - 1200.0), 7) == 0
        for T in [800, 1000, 1200]:
            assert abs(fitted.get_rate_coefficient(T) / answer.get_rate_coefficient(T) - 1.0) < 0.05

    def test_change_rate(self):
        """
        Test the MultiArrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.kinetics.get_rate_coefficient(T) for T in Tlist])
        self.kinetics.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.kinetics.get_rate_coefficient(T)
            assert abs(2 * kexp - kact) < 1e-6 * kexp


class TestMultiPDepArrhenius:
    """
    Contains unit tests of the :class:`MultiPDepArrhenius` class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.Tmin = 350.0
        self.Tmax = 1500.0
        self.Pmin = 1e-1
        self.Pmax = 1e1
        self.pressures = np.array([1e-1, 1e1])
        self.comment = "CH3 + C2H6 <=> CH4 + C2H5 (Baulch 2005)"
        self.arrhenius = [
            PDepArrhenius(
                pressures=(self.pressures, "bar"),
                arrhenius=[
                    Arrhenius(
                        A=(9.3e-16, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                    Arrhenius(
                        A=(9.3e-14, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                ],
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                Pmin=(self.Pmin, "bar"),
                Pmax=(self.Pmax, "bar"),
                comment=self.comment,
            ),
            PDepArrhenius(
                pressures=(self.pressures, "bar"),
                arrhenius=[
                    Arrhenius(
                        A=(1.4e-11, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                    Arrhenius(
                        A=(1.4e-9, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                ],
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                Pmin=(self.Pmin, "bar"),
                Pmax=(self.Pmax, "bar"),
                comment=self.comment,
            ),
        ]
        self.kinetics = MultiPDepArrhenius(
            arrhenius=self.arrhenius,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )

    def test_arrhenius(self):
        """
        Test that the MultiPDepArrhenius arrhenius property was properly set.
        """
        assert self.kinetics.arrhenius == self.arrhenius

    def test_temperature_min(self):
        """
        Test that the MultiPDepArrhenius Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the MultiPDepArrhenius Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_pressure_min(self):
        """
        Test that the MultiPDepArrhenius Pmin property was properly set.
        """
        assert round(abs(self.kinetics.Pmin.value_si * 1e-5 - self.Pmin), 6) == 0

    def test_pressure_max(self):
        """
        Test that the MultiPDepArrhenius Pmax property was properly set.
        """
        assert round(abs(self.kinetics.Pmax.value_si * 1e-5 - self.Pmax), 6) == 0

    def test_comment(self):
        """
        Test that the MultiPDepArrhenius comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the MultiPDepArrhenius.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, False, False, False], np.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.is_temperature_valid(T)
            assert valid0 == valid

    def test_is_pressure_valid(self):
        """
        Test the MultiPDepArrhenius.is_pressure_valid() method.
        """
        Pdata = np.array([1e3, 1e4, 1e5, 1e6, 1e7])
        validdata = np.array([False, True, True, True, False], np.bool)
        for P, valid in zip(Pdata, validdata):
            valid0 = self.kinetics.is_pressure_valid(P)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the MultiPDepArrhenius.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        Plist = np.array([1e4, 1e5, 1e6])
        kexplist = np.array(
            [
                [
                    2.85400e-08,
                    4.00384e-03,
                    2.73563e-01,
                    8.50699e00,
                    1.20181e02,
                    7.56312e02,
                    2.84724e03,
                    7.71702e03,
                    1.67743e04,
                    3.12290e04,
                ],
                [
                    2.85400e-07,
                    4.00384e-02,
                    2.73563e00,
                    8.50699e01,
                    1.20181e03,
                    7.56312e03,
                    2.84724e04,
                    7.71702e04,
                    1.67743e05,
                    3.12290e05,
                ],
                [
                    2.85400e-06,
                    4.00384e-01,
                    2.73563e01,
                    8.50699e02,
                    1.20181e04,
                    7.56312e04,
                    2.84724e05,
                    7.71702e05,
                    1.67743e06,
                    3.12290e06,
                ],
            ]
        ).T
        for i in range(Tlist.shape[0]):
            for j in range(Plist.shape[0]):
                kexp = kexplist[i, j]
                kact = self.kinetics.get_rate_coefficient(Tlist[i], Plist[j])
                assert abs(kexp - kact) < 1e-4 * kexp

    def test_get_rate_coefficient_diff_plist(self):
        """
        Test the MultiPDepArrhenius.get_rate_coefficient() when plists are different.
        """
        # modify the MultiPDepArrhenius object with an additional entry
        pressures = np.array([1e-1, 1e-1, 1e1])
        self.kinetics.arrhenius[0].pressures = (pressures, "bar")
        self.kinetics.arrhenius[0].arrhenius.insert(0, self.kinetics.arrhenius[0].arrhenius[0])

        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        Plist = np.array([1e4, 1e5, 1e6])
        kexplist = np.array(
            [
                [
                    2.85400e-08,
                    4.00384e-03,
                    2.73563e-01,
                    8.50699e00,
                    1.20181e02,
                    7.56312e02,
                    2.84724e03,
                    7.71702e03,
                    1.67743e04,
                    3.12290e04,
                ],
                [
                    2.85400e-07,
                    4.00384e-02,
                    2.73563e00,
                    8.50699e01,
                    1.20181e03,
                    7.56312e03,
                    2.84724e04,
                    7.71702e04,
                    1.67743e05,
                    3.12290e05,
                ],
                [
                    2.85400e-06,
                    4.00384e-01,
                    2.73563e01,
                    8.50699e02,
                    1.20181e04,
                    7.56312e04,
                    2.84724e05,
                    7.71702e05,
                    1.67743e06,
                    3.12290e06,
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
        Test that a MultiPDepArrhenius object can be pickled and unpickled with
        no loss of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a MultiPDepArrhenius object can be reconstructed from its
        repr() output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_change_rate(self):
        """
        Test the PDepMultiArrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.kinetics.get_rate_coefficient(T, 1e5) for T in Tlist])
        self.kinetics.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.kinetics.get_rate_coefficient(T, 1e5)
            assert abs(2 * kexp - kact) < 1e-6 * kexp

    def test_generate_reverse_rate_coefficient(self):
        """
        Test ability to reverse a reaction rate.

        This is a real example from an imported chemkin file.
        """
        from rmgpy.species import Species
        from rmgpy.molecule import Molecule
        from rmgpy.data.kinetics import LibraryReaction
        from rmgpy.thermo import NASA, NASAPolynomial

        test_reaction = LibraryReaction(
            reactants=[
                Species(
                    label="C2H3",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    3.12502,
                                    0.00235137,
                                    2.36803e-05,
                                    -3.35092e-08,
                                    1.39444e-11,
                                    34524.3,
                                    8.81538,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    4.37211,
                                    0.00746869,
                                    -2.64716e-06,
                                    4.22753e-10,
                                    -2.44958e-14,
                                    33805.2,
                                    0.428772,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(285.696, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(108.088, "J/mol/K"),
                        comment="""ATcT3E\nC2H3 <g> ATcT ver. 1.122, DHf298 = 296.91 ± 0.33 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="[CH]=C")],
                    molecular_weight=(27.0452, "amu"),
                ),
                Species(
                    label="CH2O",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    4.77187,
                                    -0.00976266,
                                    3.70122e-05,
                                    -3.76922e-08,
                                    1.31327e-11,
                                    -14379.8,
                                    0.696586,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    2.91333,
                                    0.0067004,
                                    -2.55521e-06,
                                    4.27795e-10,
                                    -2.44073e-14,
                                    -14462.2,
                                    7.43823,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(-119.527, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(83.1447, "J/mol/K"),
                        comment="""ATcT3E\nH2CO <g> ATcT ver. 1.122, DHf298 = -109.188 ± 0.099 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="C=O")],
                    molecular_weight=(30.026, "amu"),
                ),
            ],
            products=[
                Species(
                    label="C2H4",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    3.65151,
                                    -0.00535067,
                                    5.16486e-05,
                                    -6.36869e-08,
                                    2.50743e-11,
                                    5114.51,
                                    5.38561,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    4.14446,
                                    0.0102648,
                                    -3.61247e-06,
                                    5.74009e-10,
                                    -3.39296e-14,
                                    4190.59,
                                    -1.14778,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(42.06, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(133.032, "J/mol/K"),
                        comment="""ATcT3E\nC2H4 <g> ATcT ver. 1.122, DHf298 = 52.45 ± 0.13 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="C=C")],
                    molecular_weight=(28.0532, "amu"),
                ),
                Species(
                    label="HCO",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    3.97075,
                                    -0.00149122,
                                    9.54042e-06,
                                    -8.8272e-09,
                                    2.67645e-12,
                                    3842.03,
                                    4.4466,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    3.85781,
                                    0.00264114,
                                    -7.44177e-07,
                                    1.23313e-10,
                                    -8.88959e-15,
                                    3616.43,
                                    3.92451,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(32.0237, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(58.2013, "J/mol/K"),
                        comment="""HCO <g> ATcT ver. 1.122, DHf298 = 41.803 ± 0.099 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="[CH]=O")],
                    molecular_weight=(29.018, "amu"),
                ),
            ],
            kinetics=MultiPDepArrhenius(
                arrhenius=[
                    PDepArrhenius(
                        pressures=([0.001, 0.01, 0.1, 1, 10, 100, 1000], "atm"),
                        arrhenius=[
                            Arrhenius(
                                A=(1.1e07, "cm^3/(mol*s)"),
                                n=1.09,
                                Ea=(1807, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(2.5e07, "cm^3/(mol*s)"),
                                n=0.993,
                                Ea=(1995, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(2.5e08, "cm^3/(mol*s)"),
                                n=0.704,
                                Ea=(2596, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(1.4e10, "cm^3/(mol*s)"),
                                n=0.209,
                                Ea=(3934, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(3.5e13, "cm^3/(mol*s)"),
                                n=-0.726,
                                Ea=(6944, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(3.3e14, "cm^3/(mol*s)"),
                                n=-0.866,
                                Ea=(10966, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(17, "cm^3/(mol*s)"),
                                n=3.17,
                                Ea=(9400, "cal/mol"),
                                T0=(1, "K"),
                            ),
                        ],
                    ),
                    PDepArrhenius(
                        pressures=([0.001, 0.01, 0.1, 1, 10, 100, 1000], "atm"),
                        arrhenius=[
                            Arrhenius(
                                A=(-2.3e16, "cm^3/(mol*s)"),
                                n=-1.269,
                                Ea=(20617, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-5.2e16, "cm^3/(mol*s)"),
                                n=-1.366,
                                Ea=(20805, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-1.5e18, "cm^3/(mol*s)"),
                                n=-1.769,
                                Ea=(22524, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-8.5e19, "cm^3/(mol*s)"),
                                n=-2.264,
                                Ea=(23862, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-4.4e23, "cm^3/(mol*s)"),
                                n=-3.278,
                                Ea=(27795, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-4.2e24, "cm^3/(mol*s)"),
                                n=-3.418,
                                Ea=(31817, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-2.1e11, "cm^3/(mol*s)"),
                                n=0.618,
                                Ea=(30251, "cal/mol"),
                                T0=(1, "K"),
                            ),
                        ],
                    ),
                ]
            ),
            duplicate=True,
        )
        test_reaction.generate_reverse_rate_coefficient()
