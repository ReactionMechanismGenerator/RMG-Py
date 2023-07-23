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
This script contains unit tests of the :mod:`rmgpy.quantity` module.
"""


import numpy as np

import rmgpy.constants as constants
import rmgpy.quantity as quantity
import pytest


class TestAcceleration:
    """
    Contains unit tests of the Acceleration unit type object.
    """

    def test_mpers2(self):
        """
        Test the creation of an acceleration quantity with units of m/s^2.
        """
        q = quantity.Acceleration(1.0, "m/s^2")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m/s^2"

    def test_cmpers2(self):
        """
        Test the creation of an acceleration quantity with units of cm/s^2.
        """
        q = quantity.Acceleration(1.0, "cm/s^2")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 0.01) < 1e-8
        assert q.units == "cm/s^2"


class TestArea:
    """
    Contains unit tests of the Area unit type object.
    """

    def test_m2(self):
        """
        Test the creation of an area quantity with units of m^2.
        """
        q = quantity.Area(1.0, "m^2")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m^2"

    def test_cm2(self):
        """
        Test the creation of an area quantity with units of m^2.
        """
        q = quantity.Area(1.0, "cm^2")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-4) < 1e-10
        assert q.units == "cm^2"


class TestConcentration:
    """
    Contains unit tests of the Concentration unit type object.
    """

    def test_perm3(self):
        """
        Test the creation of an concentration quantity with units of m^-3.
        """
        try:
            quantity.Concentration(1.0, "m^-3")
            assert False, 'Allowed invalid unit type "m^-3".'
        except quantity.QuantityError:
            pass

    def test_molperm3(self):
        """
        Test the creation of an concentration quantity with units of mol/m^3.
        """
        q = quantity.Concentration(1.0, "mol/m^3")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "mol/m^3"

    def test_moleculesperm3(self):
        """
        Test the creation of an concentration quantity with units of molecules/m^3.
        """
        q = quantity.Concentration(1.0, "molecules/m^3")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * constants.Na - 1.0) < 1e-6
        assert q.units == "molecules/m^3"


class TestEnergy:
    """
    Contains unit tests of the Energy unit type object.
    """

    def test_joule(self):
        """
        Test the creation of an energy quantity with units of J.
        """
        try:
            quantity.Energy(1.0, "J")
            assert False, 'Allowed invalid unit type "J".'
        except quantity.QuantityError:
            pass

    def test_joule_per_mol(self):
        """
        Test the creation of an energy quantity with units of J/mol.
        """
        q = quantity.Energy(1.0, "J/mol")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "J/mol"

    def test_cal(self):
        """
        Test the creation of an energy quantity with units of cal.
        """
        try:
            quantity.Energy(1.0, "cal")
            assert False, 'Allowed invalid unit type "cal".'
        except quantity.QuantityError:
            pass

    def test_calpermol(self):
        """
        Test the creation of an energy quantity with units of cal/mol.
        """
        q = quantity.Energy(1.0, "cal/mol")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 4.184) < 1e-6
        assert q.units == "cal/mol"

    def test_kj(self):
        """
        Test the creation of an energy quantity with units of kJ.
        """
        try:
            quantity.Energy(1.0, "kJ")
            assert False, 'Allowed invalid unit type "kJ".'
        except quantity.QuantityError:
            pass

    def test_kj_per_mol(self):
        """
        Test the creation of an energy quantity with units of kJ/mol.
        """
        q = quantity.Energy(1.0, "kJ/mol")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1000.0) < 1e-6
        assert q.units == "kJ/mol"

    def test_kcal(self):
        """
        Test the creation of an energy quantity with units of kcal.
        """
        try:
            quantity.Energy(1.0, "kcal")
            assert False, 'Allowed invalid unit type "kcal".'
        except quantity.QuantityError:
            pass

    def test_kcal_per_mol(self):
        """
        Test the creation of an energy quantity with units of kcal/mol.
        """
        q = quantity.Energy(1.0, "kcal/mol")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 4184.0) < 1e-6
        assert q.units == "kcal/mol"

    def test_kelvin(self):
        """
        Test the creation of an energy quantity with units of K (not really an energy!).
        """
        q = quantity.Energy(10.0, "K")
        assert abs(q.value - 10 * 8.314472) < 1e-6
        assert q.units == "J/mol"


class TestDipoleMoment:
    """
    Contains unit tests of the DipoleMoment unit type object.
    """

    def test_coulomb_meter(self):
        """
        Test the creation of a dipole moment quantity with units of C*m.
        """
        q = quantity.DipoleMoment(1.0, "C*m")
        assert round(abs(q.value - 1.0), 6) == 0
        assert round(abs(q.value_si - 1.0), 6) == 0
        assert q.units == "C*m"

    def test_debye(self):
        """
        Test the creation of a dipole moment quantity with units of Debye.
        """
        q = quantity.DipoleMoment(1.0, "De")
        assert round(abs(q.value - 1.0), 6) == 0
        assert round(abs(q.value_si * constants.c * 1.0e21 - 1.0), 6) == 0
        assert q.units == "De"


class TestFlux:
    """
    Contains unit tests of the Flux unit type object.
    """

    def test_perm2pers(self):
        """
        Test the creation of a flux quantity with units of m^-2*s^-1.
        """
        try:
            quantity.Flux(1.0, "m^-2*s^-1")
            assert False, 'Allowed invalid unit type "m^-2*s^-1".'
        except quantity.QuantityError:
            pass

    def test_molperm3(self):
        """
        Test the creation of a flux quantity with units of mol/(m^2*s).
        """
        q = quantity.Flux(1.0, "mol/(m^2*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "mol/(m^2*s)"

    def test_moleculesperm3(self):
        """
        Test the creation of a flux quantity with units of molecules/(m^2*s).
        """
        q = quantity.Flux(1.0, "molecules/(m^2*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * constants.Na - 1.0) < 1e-6
        assert q.units == "molecules/(m^2*s)"


class TestForce:
    """
    Contains unit tests of the Force unit type object.
    """

    def test_newton(self):
        """
        Test the creation of an force quantity with units of N.
        """
        q = quantity.Force(1.0, "N")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "N"


class TestFrequency:
    """
    Contains unit tests of the Frequency unit type object. Note that, as a
    special case, frequencies can be read in several units, but are always
    stored internally as cm^-1.
    """

    def test_cm_1(self):
        """
        Test the creation of a frequency quantity with units of cm^-1.
        """
        q = quantity.Frequency(1.0, "cm^-1")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "cm^-1"

    def test_s_1(self):
        """
        Test the creation of a frequency quantity with units of s^-1.
        """
        q = quantity.Frequency(1.0, "s^-1")
        assert abs(q.value - 1.0 / (constants.c * 100.0)) < 1e-17
        assert abs(q.value_si - 1.0 / (constants.c * 100.0)) < 1e-17
        assert q.units == "cm^-1"

    def test_k(self):
        """
        Test the creation of a frequency quantity with units of K.
        """
        q = quantity.Frequency(1.0, "K")
        assert round(abs(q.value - constants.kB / (constants.h * constants.c * 100.0)), 6) == 0
        assert abs(q.value_si - constants.kB / (constants.h * constants.c * 100.0)) < 1e-6
        assert q.units == "cm^-1"

    def test_ev(self):
        """
        Test the creation of a frequency quantity with units of eV.
        """
        q = quantity.Frequency(1.0, "eV")
        assert round(abs(q.value - constants.e / (constants.h * constants.c * 100.0)), 2) == 0
        assert abs(q.value_si - constants.e / (constants.h * constants.c * 100.0)) < 1e-2
        assert q.units == "cm^-1"

    def test_hz(self):
        """
        Test the creation of a frequency quantity with units of Hz.
        """
        q = quantity.Frequency(1.0, "Hz")
        assert abs(q.value - 1.0 / (constants.c * 100.0)) < 1e-17
        assert abs(q.value_si - 1.0 / (constants.c * 100.0)) < 1e-17
        assert q.units == "cm^-1"

    def test_khz(self):
        """
        Test the creation of a frequency quantity with units of kHz.
        """
        q = quantity.Frequency(1.0, "kHz")
        assert abs(q.value - 1e3 / (constants.c * 100.0)) < 1e-14
        assert abs(q.value_si - 1e3 / (constants.c * 100.0)) < 1e-14
        assert q.units == "cm^-1"

    def test_mhz(self):
        """
        Test the creation of a frequency quantity with units of MHz.
        """
        q = quantity.Frequency(1.0, "MHz")
        assert abs(q.value - 1e6 / (constants.c * 100.0)) < 1e-11
        assert abs(q.value_si - 1e6 / (constants.c * 100.0)) < 1e-11
        assert q.units == "cm^-1"

    def test_ghz(self):
        """
        Test the creation of a frequency quantity with units of GHz.
        """
        q = quantity.Frequency(1.0, "GHz")
        assert abs(q.value - 1e9 / (constants.c * 100.0)) < 1e-08
        assert abs(q.value_si - 1e9 / (constants.c * 100.0)) < 1e-08
        assert q.units == "cm^-1"


class TestHeatCapacity:
    """
    Contains unit tests of the HeatCapacity unit type object.
    """

    def test_joule_per_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of J/K.
        """
        try:
            quantity.HeatCapacity(1.0, "J/K")
            assert False, 'Allowed invalid unit type "J/K".'
        except quantity.QuantityError:
            pass

    def test_joule_per_mol_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of J/(mol*K).
        """
        q = quantity.HeatCapacity(1.0, "J/(mol*K)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "J/(mol*K)"

    def test_cal_per_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of cal/K.
        """
        try:
            quantity.HeatCapacity(1.0, "cal/K")
            assert False, 'Allowed invalid unit type "cal/K".'
        except quantity.QuantityError:
            pass

    def test_cal_per_mol_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of cal/(mol*K).
        """
        q = quantity.HeatCapacity(1.0, "cal/(mol*K)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 4.184) < 1e-6
        assert q.units == "cal/(mol*K)"

    def test_kj_per_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of kJ/K.
        """
        try:
            quantity.HeatCapacity(1.0, "kJ/K")
            assert False, 'Allowed invalid unit type "kJ/K".'
        except quantity.QuantityError:
            pass

    def test_kj_per_mol_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of kJ/(mol*K).
        """
        q = quantity.HeatCapacity(1.0, "kJ/(mol*K)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1000.0) < 1e-6
        assert q.units == "kJ/(mol*K)"

    def test_kcal_per_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of kcal/K.
        """
        try:
            quantity.HeatCapacity(1.0, "kcal/K")
            assert False, 'Allowed invalid unit type "kcal/K".'
        except quantity.QuantityError:
            pass

    def test_kcal_per_mol_kelvin(self):
        """
        Test the creation of a heat capacity quantity with units of kcal/(mol*K).
        """
        q = quantity.HeatCapacity(1.0, "kcal/(mol*K)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 4184.0) < 1e-6
        assert q.units == "kcal/(mol*K)"


class TestInertia:
    """
    Contains unit tests of the Inertia unit type object.
    """

    def test_kg_m2(self):
        """
        Test the creation of a moment of inertia quantity with units of kg*m^2.
        """
        q = quantity.Inertia(1.0, "kg*m^2")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "kg*m^2"

    def test_amu_angstrom2(self):
        """
        Test the creation of a moment of inertia quantity with units of amu*angstrom^2.
        """
        q = quantity.Inertia(1.0, "amu*angstrom^2")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * constants.Na * 1e23 - 1.0) < 1e-6
        assert q.units == "amu*angstrom^2"


class TestLength:
    """
    Contains unit tests of the Length unit type object.
    """

    def test_m(self):
        """
        Test the creation of a length quantity with units of m.
        """
        q = quantity.Length(1.0, "m")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m"

    def test_km(self):
        """
        Test the creation of a length quantity with units of km.
        """
        q = quantity.Length(1.0, "km")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e3) < 1e-3
        assert q.units == "km"

    def test_cm(self):
        """
        Test the creation of a length quantity with units of cm.
        """
        q = quantity.Length(1.0, "cm")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-2) < 1e-8
        assert q.units == "cm"

    def test_mm(self):
        """
        Test the creation of a length quantity with units of mm.
        """
        q = quantity.Length(1.0, "mm")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-3) < 1e-9
        assert q.units == "mm"

    def test_um(self):
        """
        Test the creation of a length quantity with units of um.
        """
        q = quantity.Length(1.0, "um")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-6) < 1e-12
        assert q.units == "um"

    def test_nm(self):
        """
        Test the creation of a length quantity with units of nm.
        """
        q = quantity.Length(1.0, "nm")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-9) < 1e-15
        assert q.units == "nm"

    def test_pm(self):
        """
        Test the creation of a length quantity with units of pm.
        """
        q = quantity.Length(1.0, "pm")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-12) < 1e-18
        assert q.units == "pm"


class TestMass:
    """
    Contains unit tests of the Mass unit type object.

    Note that value_si is always kg (per molecule), not kg/mol.
    """

    def test_kg(self):
        """
        Test the creation of a mass quantity with units of kg.
        """
        q = quantity.Mass(1.0, "kg")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "kg"

    def test_gpermol(self):
        """
        Test the creation of a mass quantity with units of g/mol.
        Note that g/mol is automatically coerced to amu.
        """
        q = quantity.Mass(1.0, "g/mol")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - constants.amu) < 1e-32
        assert q.units == "amu"

    def test_kgpermol(self):
        """
        Test the creation of a mass quantity with units of kg/mol.
        Note that kg/mol is automatically coerced to amu.
        """
        q = quantity.Mass(1.0, "kg/mol")
        assert round(abs(q.value - 1000.0), 3) == 0
        assert abs(q.value_si - 1000.0 * constants.amu) < 1e-29
        assert q.units == "amu"

    def test_amu(self):
        """
        Test the creation of a mass quantity with units of amu.
        """
        q = quantity.Mass(1.0, "amu")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - constants.amu) < 1e-32
        assert q.units == "amu"


class TestMomentum:
    """
    Contains unit tests of the Momentum unit type object.
    """

    def test_kgmpers2(self):
        """
        Test the creation of a momentum quantity with units of kg*m/s^2.
        """
        q = quantity.Momentum(1.0, "kg*m/s^2")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "kg*m/s^2"


class TestPower:
    """
    Contains unit tests of the Power unit type object.
    """

    def test_w(self):
        """
        Test the creation of a power quantity with units of W.
        """
        q = quantity.Power(1.0, "W")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "W"


class TestPressure:
    """
    Contains unit tests of the Pressure unit type object.
    """

    def test_pa(self):
        """
        Test the creation of a pressure quantity with units of Pa.
        """
        q = quantity.Pressure(1.0, "Pa")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "Pa"

    def test_bar(self):
        """
        Test the creation of a pressure quantity with units of bar.
        """
        q = quantity.Pressure(1.0, "bar")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e5) < 1e-6
        assert q.units == "bar"

    def test_atm(self):
        """
        Test the creation of a pressure quantity with units of atm.
        """
        q = quantity.Pressure(1.0, "atm")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 101325.0) < 1e-6
        assert q.units == "atm"

    def test_torr(self):
        """
        Test the creation of a pressure quantity with units of torr.
        """
        q = quantity.Pressure(1.0, "torr")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 101325.0 / 760.0) < 1e-6
        assert q.units == "torr"

    def test_psi(self):
        """
        Test the creation of a pressure quantity with units of psi.
        """
        q = quantity.Pressure(1.0, "psi")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 101325.0 / 14.695949) < 1e-2
        assert q.units == "psi"


class TestRateCoefficient:
    """
    Contains unit tests of the RateCoefficient unit type object.
    """

    def test_s(self):
        """
        Test the creation of a rate coefficient quantity with units of s^-1.
        """
        q = quantity.RateCoefficient(1.0, "s^-1")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "s^-1"
        assert round(abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1.0), 1) == 0  # 1 /s  =  1 /s

    def test_m3permols(self):
        """
        Test the creation of a rate coefficient quantity with units of m^3/(mol*s).
        """
        q = quantity.RateCoefficient(1.0, "m^3/(mol*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m^3/(mol*s)"
        assert round(abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e6), 1) == 0  # 1 m3/mol/s  =  1e6  cm3/mol/s

    def test_m6permol2s(self):
        """
        Test the creation of a rate coefficient quantity with units of m^6/(mol^2*s).
        """
        q = quantity.RateCoefficient(1.0, "m^6/(mol^2*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m^6/(mol^2*s)"
        assert round(abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e12), 1) == 0  # 1 m6/mol2/s  =  1e12  cm6/mol2/s

    def test_m9permol3s(self):
        """
        Test the creation of a rate coefficient quantity with units of m^9/(mol^3*s).
        """
        q = quantity.RateCoefficient(1.0, "m^9/(mol^3*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m^9/(mol^3*s)"
        assert abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e18) < 1e3  # 1 m9/mol3/s  =  1e18  cm9/mol3/s

    def test_cm3permols(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^3/(mol*s).
        """
        q = quantity.RateCoefficient(1.0, "cm^3/(mol*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * 1e6 - 1.0) < 1e-6
        assert q.units == "cm^3/(mol*s)"
        assert round(abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e6), 1) == 0  # 1 m3/mol/s  =  1 cm3/mol/s

    def test_cm6permol2s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^6/(mol^2*s).
        """
        q = quantity.RateCoefficient(1.0, "cm^6/(mol^2*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * 1e6**2 - 1.0) < 1e-6
        assert q.units == "cm^6/(mol^2*s)"
        assert round(abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e12), 1) == 0  # 1 m6/mol2/s  =  1e12  cm6/mol2/s

    def test_cm9permol3s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^9/(mol^3*s).
        """
        q = quantity.RateCoefficient(1.0, "cm^9/(mol^3*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * 1e6**3 - 1.0) < 1e-6
        assert q.units == "cm^9/(mol^3*s)"
        assert abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e18) < 1e3  # 1 m9/mol3/s  =  1e18  cm9/mol3/s

    def test_cm3permolecules(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^3/(molecule*s).
        """
        q = quantity.RateCoefficient(1.0, "cm^3/(molecule*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * 1e6 / constants.Na - 1.0) < 1e-6
        assert q.units == "cm^3/(molecule*s)"
        assert abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e6) < 1e0  # 1 m3/mol/s  =  1e6 cm3/mol/s

    def test_cm6permolecule2s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^6/(molecule^2*s).
        """
        q = quantity.RateCoefficient(1.0, "cm^6/(molecule^2*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * (1e6 / constants.Na) ** 2 - 1.0) < 1e-6
        assert q.units == "cm^6/(molecule^2*s)"
        assert abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e12) < 1e0  # 1 m6/mol2/s  =  1e12 cm6/mol2/s

    def test_cm9permolecule3s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^9/(molecule^3*s).
        """
        q = quantity.RateCoefficient(1.0, "cm^9/(molecule^3*s)")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si * (1e6 / constants.Na) ** 3 - 1.0) < 1e-6
        assert q.units == "cm^9/(molecule^3*s)"
        assert abs(q.get_conversion_factor_from_si_to_cm_mol_s() - 1e18) < 1e3  # 1 m9/mole3/s  =  1e18 cm9/mol3/s


class TestTemperature:
    """
    Contains unit tests of the Temperature unit type object.
    """

    def test_k(self):
        """
        Test the creation of a temperature quantity with units of K.
        """
        q = quantity.Temperature(1.0, "K")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "K"

    def test_deg_c(self):
        """
        Test the creation of a temperature quantity with units of degrees C.
        """
        with pytest.raises(NotImplementedError):
            quantity.Temperature(1.0, "degC")

    def test_deg_f(self):
        """
        Test the creation of a temperature quantity with units of degrees F.
        """
        with pytest.raises(NotImplementedError):
            quantity.Temperature(1.0, "degF")

    def test_deg_r(self):
        """
        Test the creation of a temperature quantity with units of degrees R.
        """
        with pytest.raises(NotImplementedError):
            quantity.Temperature(1.0, "degR")


class TestTime:
    """
    Contains unit tests of the Time unit type object.
    """

    def test_s(self):
        """
        Test the creation of a time quantity with units of s.
        """
        q = quantity.Time(1.0, "s")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "s"

    def test_ms(self):
        """
        Test the creation of a time quantity with units of ms.
        """
        q = quantity.Time(1.0, "ms")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-3) < 1e-9
        assert q.units == "ms"

    def test_us(self):
        """
        Test the creation of a time quantity with units of us.
        """
        q = quantity.Time(1.0, "us")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-6) < 1e-12
        assert q.units == "us"

    def test_ns(self):
        """
        Test the creation of a time quantity with units of ns.
        """
        q = quantity.Time(1.0, "ns")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-9) < 1e-15
        assert q.units == "ns"

    def test_ps(self):
        """
        Test the creation of a time quantity with units of ps.
        """
        q = quantity.Time(1.0, "ps")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-12) < 1e-18
        assert q.units == "ps"

    def test_fs(self):
        """
        Test the creation of a time quantity with units of fs.
        """
        q = quantity.Time(1.0, "fs")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-15) < 1e-21
        assert q.units == "fs"

    def test_min(self):
        """
        Test the creation of a time quantity with units of min.
        """
        q = quantity.Time(1.0, "min")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 60.0) < 1e-6
        assert q.units == "min"

    def test_hr(self):
        """
        Test the creation of a time quantity with units of hr.
        """
        q = quantity.Time(1.0, "hr")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 3600.0) < 1e-6
        assert q.units == "hr"


class TestVelocity:
    """
    Contains unit tests of the Velocity unit type object.
    """

    def test_mpers(self):
        """
        Test the creation of an velocity quantity with units of m/s.
        """
        q = quantity.Velocity(1.0, "m/s")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m/s"

    def test_cmpers(self):
        """
        Test the creation of an velocity quantity with units of m/s.
        """
        q = quantity.Velocity(1.0, "cm/s")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 0.01) < 1e-8
        assert q.units == "cm/s"


class TestVolume:
    """
    Contains unit tests of the Volume unit type object.
    """

    def test_m3(self):
        """
        Test the creation of an volume quantity with units of m^3.
        """
        q = quantity.Volume(1.0, "m^3")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0) < 1e-6
        assert q.units == "m^3"

    def test_liters(self):
        """
        Test the creation of an volume quantity with units of L.
        """
        q = quantity.Volume(1.0, "L")
        assert round(abs(q.value - 1.0), 6) == 0
        assert abs(q.value_si - 1.0e-3) < 1e-9
        assert q.units == "L"


class TestQuantity:
    """
    Contains unit tests testing the value and uncertainty storage behavior for ScalarQuantity and ArrayQuantity objects
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.  This tests the creation of several both ScalarQuantity
        and ArrayQuantity objects
        """
        self.Cp = quantity.Quantity(
            [-6.51, -5.19333, -4.47333, -3.76, -3.44333, -2.94667, -2.47],
            "cal/(mol*K)",
            "+|-",
            [2.72057, 3.42407, 4.84068, 5.11681, 5.13207, 5.8757, 8.29108],
        )
        self.v = quantity.Quantity([5, 10, 12], "cm/s", "*|/", [1.2, 0.4, 1])
        self.H = quantity.Quantity(33.1097, "kcal/mol", "+|-", 24.8344)
        self.A = quantity.Quantity(7.25e13, "cm^3/(mol*s)", "*|/", 5)
        self.Cp_array = quantity.ArrayQuantity(
            [-6.51, -5.19333, -4.47333, -3.76, -3.44333, -2.94667, -2.47],
            "cal/(mol*K)",
            [2.72057, 3.42407, 4.84068, 5.11681, 5.13207, 5.8757, 8.29108],
            "+|-",
        )
        self.v_array = quantity.ArrayQuantity([5, 10, 12], "cm/s", [1.2, 0.4, 1], "*|/")
        self.H_scalar = quantity.ScalarQuantity(
            33.1097,
            "kcal/mol",
            24.8344,
            "+|-",
        )
        self.A_scalar = quantity.ScalarQuantity(7.25e13, "cm^3/(mol*s)", 5, "*|/")

    def test_scalar_conversion(self):
        """
        ScalarQuantity: test that the value and uncertainty get converted to the proper si value.
        """
        # Uncertainty of type +|- must be adjusted by units
        assert round(abs(self.H.value_si - self.H.value * 4184), 7) == 0
        assert round(abs(self.H.uncertainty_si - self.H.uncertainty * 4184), 7) == 0
        assert round(abs(self.H_scalar.value_si - self.H_scalar.value * 4184), 7) == 0
        assert round(abs(self.H_scalar.uncertainty_si - self.H_scalar.uncertainty * 4184), 7) == 0

        # Uncertainty of type *|/ does not need to be adjusted by units
        assert round(abs(self.A.value_si - self.A.value * 1e-6), 7) == 0
        assert round(abs(self.A.uncertainty_si - self.A.uncertainty), 7) == 0
        assert round(abs(self.A_scalar.value_si - self.A_scalar.value * 1e-6), 7) == 0
        assert round(abs(self.A_scalar.uncertainty_si - self.A_scalar.uncertainty), 7) == 0

    def test_array_conversion(self):
        """
        ArrayQuantity: test that the value and uncertainty get converted to the proper si value.
        """
        np.testing.assert_array_almost_equal(self.v.value_si, self.v.value * 1e-2)
        np.testing.assert_array_almost_equal(self.v.uncertainty_si, self.v.uncertainty)
        np.testing.assert_array_almost_equal(self.v_array.value_si, self.v.value * 1e-2)
        np.testing.assert_array_almost_equal(self.v_array.uncertainty_si, self.v.uncertainty)

        np.testing.assert_array_almost_equal(self.Cp.value_si, self.Cp.value * 4.184)
        np.testing.assert_array_almost_equal(self.Cp.uncertainty_si, self.Cp.uncertainty * 4.184)
        np.testing.assert_array_almost_equal(self.Cp_array.value_si, self.Cp.value * 4.184)
        np.testing.assert_array_almost_equal(self.Cp_array.uncertainty_si, self.Cp.uncertainty * 4.184)

    def test_scalar_repr(self):
        """
        Test that the ScalarQuantity objects can be recreated using their __repr__ function
        """
        # Test that the values can be reconstituted
        H = quantity.Quantity(eval(repr(self.H)))
        assert H.value_si == self.H.value_si
        assert H.uncertainty_si == self.H.uncertainty_si
        assert H.uncertainty_type == self.H.uncertainty_type
        assert H.units == self.H.units

        A = quantity.Quantity(eval(repr(self.A)))
        assert A.value_si == self.A.value_si
        assert A.uncertainty_si == self.A.uncertainty_si
        assert A.uncertainty_type == self.A.uncertainty_type
        assert A.units == self.A.units

        # Test that the __repr__ strings are the same
        assert repr(H) == repr(self.H)
        assert repr(self.H) == repr(self.H_scalar)
        assert repr(A) == repr(self.A)
        assert repr(self.A) == repr(self.A_scalar)

    def test_array_repr(self):
        """
        Test that the ArrayQuantity objects can be recreated using their __repr__ function
        """
        # Test that the values can be reconstituted
        Cp = quantity.Quantity(eval(repr(self.Cp)))
        np.testing.assert_array_almost_equal(Cp.value_si, self.Cp.value_si)
        np.testing.assert_array_almost_equal(Cp.uncertainty_si, self.Cp.uncertainty_si)
        assert Cp.uncertainty_type == self.Cp.uncertainty_type
        assert Cp.units == self.Cp.units

        v = quantity.Quantity(eval(repr(self.v)))
        np.testing.assert_array_almost_equal(v.value_si, self.v.value_si)
        np.testing.assert_array_almost_equal(v.uncertainty_si, self.v.uncertainty_si)
        assert v.uncertainty_type == self.v.uncertainty_type
        assert v.units == self.v.units

        # Test that the __repr__ strings are the same
        assert repr(Cp) == repr(self.Cp)
        assert repr(self.Cp) == repr(self.Cp_array)
        assert repr(v) == repr(self.v)
        assert repr(self.v) == repr(self.v_array)


class TestQuantityDictionaryConversion:
    """
    Test that Scalar and Array Quantity objects can be represented and reconstructed from dictionaries
    """

    def setup_class(self):
        """
        Initialize necessary variables for the TestQuantityDictionaryConversion unit test
        """
        self.class_dict = {
            "ScalarQuantity": quantity.ScalarQuantity,
            "ArrayQuantity": quantity.ArrayQuantity,
            "np_array": np.array,
        }

        self.empty_scalar = quantity.ScalarQuantity()
        self.minimal_scalar = quantity.ScalarQuantity(value=5)
        self.know_scalar = quantity.ScalarQuantity(value=2.4, units="kcal/mol")
        self.uncertain_scalar = quantity.ScalarQuantity(value=3, uncertainty=0.2)

        self.empty_array = quantity.ArrayQuantity()
        self.minimal_array = quantity.ArrayQuantity(value=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
        self.known_array = quantity.ArrayQuantity(
            value=np.array([[1.2, 2.4, 3.4], [4.8, 5.0, 6.0], [7.4, 8.6, 9]]),
            units="kcal/mol",
        )
        self.uncertain_array = quantity.ArrayQuantity(
            value=np.array([[1.2, 2.4, 3.4], [4.8, 5.0, 6.0], [7.4, 8.6, 9.0]]),
            uncertainty=np.array([[0.2, 0.4, 0.6], [0.6, 0.4, 0.2], [0.8, 0.2, 0.4]]),
        )

    def test_scalar_as_dict(self):
        """
        Test the `as_dict` method of ScalarQuantity objects
        """
        assert self.empty_scalar.as_dict() == {"class": "ScalarQuantity", "value": 0.0}
        assert self.minimal_scalar.as_dict() == {"class": "ScalarQuantity", "value": 5}
        assert self.know_scalar.as_dict() == {"class": "ScalarQuantity", "value": 2.4, "units": "kcal/mol"}
        assert self.uncertain_scalar.as_dict() == {
            "class": "ScalarQuantity",
            "value": 3,
            "uncertainty": 0.2,
            "uncertainty_type": "+|-",
        }

    def test_scalar_make_object(self):
        """
        Test the `make_object` method of ScalarQuantity objects
        """
        empty_scalar = quantity.ScalarQuantity()
        minimal_scalar = quantity.ScalarQuantity()
        known_scalar = quantity.ScalarQuantity()
        uncertain_scalar = quantity.ScalarQuantity()

        empty_scalar.make_object({}, self.class_dict)
        minimal_scalar.make_object({"value": 5}, self.class_dict)
        known_scalar.make_object({"value": 2.4, "units": "kcal/mol"}, self.class_dict)
        uncertain_scalar.make_object(
            {
                "class": "ScalarQuantity",
                "value": 3,
                "uncertainty": 0.2,
                "uncertainty_type": "+|-",
            },
            self.class_dict,
        )

        assert empty_scalar.as_dict() == self.empty_scalar.as_dict()
        assert minimal_scalar.as_dict() == self.minimal_scalar.as_dict()
        assert known_scalar.as_dict() == self.know_scalar.as_dict()
        assert uncertain_scalar.as_dict() == self.uncertain_scalar.as_dict()

    def test_array_as_dict(self):
        """
        Test the `as_dict` method of ArrayQuantity objects
        """
        assert self.empty_array.as_dict() == {"class": "ArrayQuantity", "value": {"class": "np_array", "object": [0.0]}}

        assert self.minimal_array.as_dict() == {
            "class": "ArrayQuantity",
            "value": {
                "class": "np_array",
                "object": [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            },
        }

        assert self.known_array.as_dict() == {
            "class": "ArrayQuantity",
            "value": {
                "class": "np_array",
                "object": [[1.2, 2.4, 3.4], [4.8, 5.0, 6.0], [7.4, 8.6, 9]],
            },
            "units": "kcal/mol",
        }

        assert self.uncertain_array.as_dict() == {
            "class": "ArrayQuantity",
            "value": {
                "class": "np_array",
                "object": [[1.2, 2.4, 3.4], [4.8, 5.0, 6.0], [7.4, 8.6, 9.0]],
            },
            "uncertainty": {
                "class": "np_array",
                "object": [[0.2, 0.4, 0.6], [0.6, 0.4, 0.2], [0.8, 0.2, 0.4]],
            },
            "uncertainty_type": "+|-",
        }

    def test_array_make_object(self):
        """
        Test the `make_object` method of ArrayQuantity objects
        """
        empty_array = quantity.ArrayQuantity()
        minimal_array = quantity.ArrayQuantity()
        known_array = quantity.ArrayQuantity()
        uncertain_array = quantity.ArrayQuantity()

        minimal_dict = {
            "class": "ArrayQuantity",
            "value": {"class": "np_array", "object": [[1, 2, 3], [4, 5, 6], [7, 8, 9]]},
        }

        known_dict = {
            "class": "ArrayQuantity",
            "value": {
                "class": "np_array",
                "object": [[1.2, 2.4, 3.4], [4.8, 5.0, 6.0], [7.4, 8.6, 9]],
            },
            "units": "kcal/mol",
        }

        uncertain_dict = {
            "class": "ArrayQuantity",
            "value": {
                "class": "np_array",
                "object": [[1.2, 2.4, 3.4], [4.8, 5.0, 6.0], [7.4, 8.6, 9.0]],
            },
            "uncertainty": {
                "class": "np_array",
                "object": [[0.2, 0.4, 0.6], [0.6, 0.4, 0.2], [0.8, 0.2, 0.4]],
            },
            "uncertainty_type": "+|-",
        }

        empty_array.make_object({}, self.class_dict)
        minimal_array.make_object(minimal_dict, self.class_dict)
        known_array.make_object(known_dict, self.class_dict)
        uncertain_array.make_object(uncertain_dict, self.class_dict)

        assert empty_array.as_dict() == self.empty_array.as_dict()
        assert minimal_array.as_dict() == self.minimal_array.as_dict()
        assert known_array.as_dict() == self.known_array.as_dict()
        assert uncertain_array.as_dict() == self.uncertain_array.as_dict()
