#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This script contains unit tests of the :mod:`rmgpy.quantity` module.
"""

import unittest
import math
import numpy

import rmgpy.constants as constants
import rmgpy.quantity as quantity

################################################################################

class TestAcceleration(unittest.TestCase):
    """
    Contains unit tests of the Acceleration unit type object.
    """
            
    def test_mpers2(self):
        """
        Test the creation of an acceleration quantity with units of m/s^2.
        """
        q = quantity.Acceleration(1.0,"m/s^2")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m/s^2")

    def test_cmpers2(self):
        """
        Test the creation of an acceleration quantity with units of cm/s^2.
        """
        q = quantity.Acceleration(1.0,"cm/s^2")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 0.01, delta=1e-8)
        self.assertEqual(q.units, "cm/s^2")

################################################################################

class TestArea(unittest.TestCase):
    """
    Contains unit tests of the Area unit type object.
    """
            
    def test_m2(self):
        """
        Test the creation of an area quantity with units of m^2.
        """
        q = quantity.Area(1.0,"m^2")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m^2")

    def test_cm2(self):
        """
        Test the creation of an area quantity with units of m^2.
        """
        q = quantity.Area(1.0,"cm^2")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-4, delta=1e-10)
        self.assertEqual(q.units, "cm^2")

################################################################################

class TestConcentration(unittest.TestCase):
    """
    Contains unit tests of the Concentration unit type object.
    """
            
    def test_perm3(self):
        """
        Test the creation of an concentration quantity with units of m^-3.
        """
        try:
            q = quantity.Concentration(1.0,"m^-3")
            self.fail('Allowed invalid unit type "m^-3".')
        except quantity.QuantityError:
            pass

    def test_molperm3(self):
        """
        Test the creation of an concentration quantity with units of mol/m^3.
        """
        q = quantity.Concentration(1.0,"mol/m^3")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "mol/m^3")

    def test_moleculesperm3(self):
        """
        Test the creation of an concentration quantity with units of molecules/m^3.
        """
        q = quantity.Concentration(1.0,"molecules/m^3")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*constants.Na, 1.0, delta=1e-6)
        self.assertEqual(q.units, "molecules/m^3")

################################################################################

class TestEnergy(unittest.TestCase):
    """
    Contains unit tests of the Energy unit type object.
    """
            
    def test_J(self):
        """
        Test the creation of an energy quantity with units of J.
        """
        try:
            q = quantity.Energy(1.0,"J")
            self.fail('Allowed invalid unit type "J".')
        except quantity.QuantityError:
            pass

    def test_Jpermol(self):
        """
        Test the creation of an energy quantity with units of J/mol.
        """
        q = quantity.Energy(1.0,"J/mol")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "J/mol")

    def test_cal(self):
        """
        Test the creation of an energy quantity with units of cal.
        """
        try:
            q = quantity.Energy(1.0,"cal")
            self.fail('Allowed invalid unit type "cal".')
        except quantity.QuantityError:
            pass

    def test_calpermol(self):
        """
        Test the creation of an energy quantity with units of cal/mol.
        """
        q = quantity.Energy(1.0,"cal/mol")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 4.184, delta=1e-6)
        self.assertEqual(q.units, "cal/mol")

    def test_kJ(self):
        """
        Test the creation of an energy quantity with units of kJ.
        """
        try:
            q = quantity.Energy(1.0,"kJ")
            self.fail('Allowed invalid unit type "kJ".')
        except quantity.QuantityError:
            pass

    def test_kJpermol(self):
        """
        Test the creation of an energy quantity with units of kJ/mol.
        """
        q = quantity.Energy(1.0,"kJ/mol")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1000., delta=1e-6)
        self.assertEqual(q.units, "kJ/mol")

    def test_kcal(self):
        """
        Test the creation of an energy quantity with units of kcal.
        """
        try:
            q = quantity.Energy(1.0,"kcal")
            self.fail('Allowed invalid unit type "kcal".')
        except quantity.QuantityError:
            pass

    def test_kcalpermol(self):
        """
        Test the creation of an energy quantity with units of kcal/mol.
        """
        q = quantity.Energy(1.0,"kcal/mol")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 4184., delta=1e-6)
        self.assertEqual(q.units, "kcal/mol")

################################################################################

class TestDipoleMoment(unittest.TestCase):
    """
    Contains unit tests of the DipoleMoment unit type object.
    """
            
    def test_Ctimesm(self):
        """
        Test the creation of a dipole moment quantity with units of C*m.
        """
        q = quantity.DipoleMoment(1.0,"C*m")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, 6)
        self.assertEqual(q.units, "C*m")

    def test_D(self):
        """
        Test the creation of a dipole moment quantity with units of J/mol.
        """
        q = quantity.DipoleMoment(1.0,"De")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*constants.c*1.0e21, 1.0, 6)
        self.assertEqual(q.units, "De")

################################################################################

class TestFlux(unittest.TestCase):
    """
    Contains unit tests of the Flux unit type object.
    """
            
    def test_perm2pers(self):
        """
        Test the creation of a flux quantity with units of m^-2*s^-1.
        """
        try:
            q = quantity.Flux(1.0,"m^-2*s^-1")
            self.fail('Allowed invalid unit type "m^-2*s^-1".')
        except quantity.QuantityError:
            pass

    def test_molperm3(self):
        """
        Test the creation of a flux quantity with units of mol/(m^2*s).
        """
        q = quantity.Flux(1.0,"mol/(m^2*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "mol/(m^2*s)")

    def test_moleculesperm3(self):
        """
        Test the creation of a flux quantity with units of molecules/(m^2*s).
        """
        q = quantity.Flux(1.0,"molecules/(m^2*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*constants.Na, 1.0, delta=1e-6)
        self.assertEqual(q.units, "molecules/(m^2*s)")

################################################################################

class TestForce(unittest.TestCase):
    """
    Contains unit tests of the Force unit type object.
    """
            
    def test_N(self):
        """
        Test the creation of an force quantity with units of N.
        """
        q = quantity.Force(1.0,"N")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "N")

################################################################################

class TestFrequency(unittest.TestCase):
    """
    Contains unit tests of the Frequency unit type object. Note that, as a
    special case, frequencies can be read in several units, but are always
    stored internally as cm^-1.
    """
            
    def test_cm_1(self):
        """
        Test the creation of a frequency quantity with units of cm^-1.
        """
        q = quantity.Frequency(1.0,"cm^-1")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "cm^-1")

    def test_s_1(self):
        """
        Test the creation of a frequency quantity with units of s^-1.
        """
        q = quantity.Frequency(1.0,"s^-1")
        self.assertAlmostEqual(q.value, 1./(constants.c*100.), delta=1e-17)
        self.assertAlmostEqual(q.value_si, 1./(constants.c*100.), delta=1e-17)
        self.assertEqual(q.units, "cm^-1")

    def test_K(self):
        """
        Test the creation of a frequency quantity with units of K.
        """
        q = quantity.Frequency(1.0,"K")
        self.assertAlmostEqual(q.value, constants.kB/(constants.h*constants.c*100.), 6)
        self.assertAlmostEqual(q.value_si, constants.kB/(constants.h*constants.c*100.), delta=1e-6)
        self.assertEqual(q.units, "cm^-1")

    def test_eV(self):
        """
        Test the creation of a frequency quantity with units of eV.
        """
        q = quantity.Frequency(1.0,"eV")
        self.assertAlmostEqual(q.value, constants.e/(constants.h*constants.c*100.), 2)
        self.assertAlmostEqual(q.value_si, constants.e/(constants.h*constants.c*100.), delta=1e-2)
        self.assertEqual(q.units, "cm^-1")

    def test_Hz(self):
        """
        Test the creation of a frequency quantity with units of Hz.
        """
        q = quantity.Frequency(1.0,"Hz")
        self.assertAlmostEqual(q.value, 1./(constants.c*100.), delta=1e-17)
        self.assertAlmostEqual(q.value_si, 1./(constants.c*100.), delta=1e-17)
        self.assertEqual(q.units, "cm^-1")

    def test_kHz(self):
        """
        Test the creation of a frequency quantity with units of kHz.
        """
        q = quantity.Frequency(1.0,"kHz")
        self.assertAlmostEqual(q.value, 1e3/(constants.c*100.), delta=1e-14)
        self.assertAlmostEqual(q.value_si, 1e3/(constants.c*100.), delta=1e-14)
        self.assertEqual(q.units, "cm^-1")

    def test_MHz(self):
        """
        Test the creation of a frequency quantity with units of MHz.
        """
        q = quantity.Frequency(1.0,"MHz")
        self.assertAlmostEqual(q.value, 1e6/(constants.c*100.), delta=1e-11)
        self.assertAlmostEqual(q.value_si, 1e6/(constants.c*100.), delta=1e-11)
        self.assertEqual(q.units, "cm^-1")

    def test_GHz(self):
        """
        Test the creation of a frequency quantity with units of GHz.
        """
        q = quantity.Frequency(1.0,"GHz")
        self.assertAlmostEqual(q.value, 1e9/(constants.c*100.), delta=1e-08)
        self.assertAlmostEqual(q.value_si, 1e9/(constants.c*100.), delta=1e-08)
        self.assertEqual(q.units, "cm^-1")

################################################################################

class TestHeatCapacity(unittest.TestCase):
    """
    Contains unit tests of the HeatCapacity unit type object.
    """
            
    def test_JperK(self):
        """
        Test the creation of a heat capacity quantity with units of J/K.
        """
        try:
            q = quantity.HeatCapacity(1.0,"J/K")
            self.fail('Allowed invalid unit type "J/K".')
        except quantity.QuantityError:
            pass

    def test_JpermolperK(self):
        """
        Test the creation of a heat capacity quantity with units of J/(mol*K).
        """
        q = quantity.HeatCapacity(1.0,"J/(mol*K)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "J/(mol*K)")

    def test_calperK(self):
        """
        Test the creation of a heat capacity quantity with units of cal/K.
        """
        try:
            q = quantity.HeatCapacity(1.0,"cal/K")
            self.fail('Allowed invalid unit type "cal/K".')
        except quantity.QuantityError:
            pass

    def test_calpermolperK(self):
        """
        Test the creation of a heat capacity quantity with units of cal/(mol*K).
        """
        q = quantity.HeatCapacity(1.0,"cal/(mol*K)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 4.184, delta=1e-6)
        self.assertEqual(q.units, "cal/(mol*K)")

    def test_kJperK(self):
        """
        Test the creation of a heat capacity quantity with units of kJ/K.
        """
        try:
            q = quantity.HeatCapacity(1.0,"kJ/K")
            self.fail('Allowed invalid unit type "kJ/K".')
        except quantity.QuantityError:
            pass

    def test_kJpermolperK(self):
        """
        Test the creation of a heat capacity quantity with units of kJ/(mol*K).
        """
        q = quantity.HeatCapacity(1.0,"kJ/(mol*K)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1000., delta=1e-6)
        self.assertEqual(q.units, "kJ/(mol*K)")

    def test_kcalperK(self):
        """
        Test the creation of a heat capacity quantity with units of kcal/K.
        """
        try:
            q = quantity.HeatCapacity(1.0,"kcal/K")
            self.fail('Allowed invalid unit type "kcal/K".')
        except quantity.QuantityError:
            pass

    def test_kcalpermolperK(self):
        """
        Test the creation of a heat capacity quantity with units of kcal/(mol*K).
        """
        q = quantity.HeatCapacity(1.0,"kcal/(mol*K)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 4184., delta=1e-6)
        self.assertEqual(q.units, "kcal/(mol*K)")

################################################################################

class TestInertia(unittest.TestCase):
    """
    Contains unit tests of the Inertia unit type object.
    """
            
    def test_kg_m2(self):
        """
        Test the creation of a moment of inertia quantity with units of kg*m^2.
        """
        q = quantity.Inertia(1.0,"kg*m^2")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "kg*m^2")

    def test_amu_angstrom2(self):
        """
        Test the creation of a moment of inertia quantity with units of amu*angstrom^2.
        """
        q = quantity.Inertia(1.0,"amu*angstrom^2")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*constants.Na*1e23, 1.0, delta=1e-6)
        self.assertEqual(q.units, "amu*angstrom^2")

################################################################################

class TestLength(unittest.TestCase):
    """
    Contains unit tests of the Length unit type object.
    """
            
    def test_m(self):
        """
        Test the creation of a length quantity with units of m.
        """
        q = quantity.Length(1.0,"m")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m")

    def test_km(self):
        """
        Test the creation of a length quantity with units of km.
        """
        q = quantity.Length(1.0,"km")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e3, delta=1e-3)
        self.assertEqual(q.units, "km")

    def test_cm(self):
        """
        Test the creation of a length quantity with units of cm.
        """
        q = quantity.Length(1.0,"cm")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-2, delta=1e-8)
        self.assertEqual(q.units, "cm")

    def test_mm(self):
        """
        Test the creation of a length quantity with units of mm.
        """
        q = quantity.Length(1.0,"mm")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-3, delta=1e-9)
        self.assertEqual(q.units, "mm")

    def test_um(self):
        """
        Test the creation of a length quantity with units of um.
        """
        q = quantity.Length(1.0,"um")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-6, delta=1e-12)
        self.assertEqual(q.units, "um")

    def test_nm(self):
        """
        Test the creation of a length quantity with units of nm.
        """
        q = quantity.Length(1.0,"nm")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-9, delta=1e-15)
        self.assertEqual(q.units, "nm")

    def test_pm(self):
        """
        Test the creation of a length quantity with units of pm.
        """
        q = quantity.Length(1.0,"pm")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-12, delta=1e-18)
        self.assertEqual(q.units, "pm")

################################################################################

class TestMass(unittest.TestCase):
    """
    Contains unit tests of the Mass unit type object.
    """
            
    def test_kg(self):
        """
        Test the creation of a mass quantity with units of kg.
        """
        q = quantity.Mass(1.0,"kg")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "kg")

    def test_gpermol(self):
        """
        Test the creation of a mass quantity with units of g/mol. Note that
        g/mol is automatically coerced to amu.
        """
        q = quantity.Mass(1.0,"g/mol")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, constants.amu, delta=1e-32)
        self.assertEqual(q.units, "amu")

    def test_kgpermol(self):
        """
        Test the creation of a mass quantity with units of kg/mol. Note that
        kg/mol is automatically coerced to amu.
        """
        q = quantity.Mass(1.0,"kg/mol")
        self.assertAlmostEqual(q.value, 1000.0, 3)
        self.assertAlmostEqual(q.value_si, 1000.*constants.amu, delta=1e-29)
        self.assertEqual(q.units, "amu")

    def test_amu(self):
        """
        Test the creation of a mass quantity with units of amu.
        """
        q = quantity.Mass(1.0,"amu")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, constants.amu, delta=1e-32)
        self.assertEqual(q.units, "amu")

################################################################################

class TestMomentum(unittest.TestCase):
    """
    Contains unit tests of the Momentum unit type object.
    """
            
    def test_kgmpers2(self):
        """
        Test the creation of a momentum quantity with units of kg*m/s^2.
        """
        q = quantity.Momentum(1.0,"kg*m/s^2")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "kg*m/s^2")

################################################################################

class TestPower(unittest.TestCase):
    """
    Contains unit tests of the Power unit type object.
    """
            
    def test_W(self):
        """
        Test the creation of a power quantity with units of W.
        """
        q = quantity.Power(1.0,"W")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "W")

################################################################################

class TestPressure(unittest.TestCase):
    """
    Contains unit tests of the Pressure unit type object.
    """
            
    def test_Pa(self):
        """
        Test the creation of a pressure quantity with units of Pa.
        """
        q = quantity.Pressure(1.0,"Pa")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "Pa")

    def test_bar(self):
        """
        Test the creation of a pressure quantity with units of bar.
        """
        q = quantity.Pressure(1.0,"bar")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e5, delta=1e-6)
        self.assertEqual(q.units, "bar")

    def test_atm(self):
        """
        Test the creation of a pressure quantity with units of atm.
        """
        q = quantity.Pressure(1.0,"atm")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 101325., delta=1e-6)
        self.assertEqual(q.units, "atm")

    def test_torr(self):
        """
        Test the creation of a pressure quantity with units of torr.
        """
        q = quantity.Pressure(1.0,"torr")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 101325./760., delta=1e-6)
        self.assertEqual(q.units, "torr")

    def test_psi(self):
        """
        Test the creation of a pressure quantity with units of psi.
        """
        q = quantity.Pressure(1.0,"psi")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 101325./14.695949, delta=1e-2)
        self.assertEqual(q.units, "psi")

################################################################################

class TestRateCoefficient(unittest.TestCase):
    """
    Contains unit tests of the RateCoefficient unit type object.
    """
            
    def test_s(self):
        """
        Test the creation of a rate coefficient quantity with units of s^-1.
        """
        q = quantity.RateCoefficient(1.0,"s^-1")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "s^-1")

    def test_m3permols(self):
        """
        Test the creation of a rate coefficient quantity with units of m^3/(mol*s).
        """
        q = quantity.RateCoefficient(1.0,"m^3/(mol*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m^3/(mol*s)")

    def test_m6permol2s(self):
        """
        Test the creation of a rate coefficient quantity with units of m^6/(mol^2*s).
        """
        q = quantity.RateCoefficient(1.0,"m^6/(mol^2*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m^6/(mol^2*s)")

    def test_m9permol3s(self):
        """
        Test the creation of a rate coefficient quantity with units of m^9/(mol^3*s).
        """
        q = quantity.RateCoefficient(1.0,"m^9/(mol^3*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m^9/(mol^3*s)")

    def test_cm3permols(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^3/(mol*s).
        """
        q = quantity.RateCoefficient(1.0,"cm^3/(mol*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*1e6, 1.0, delta=1e-6)
        self.assertEqual(q.units, "cm^3/(mol*s)")

    def test_cm6permol2s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^6/(mol^2*s).
        """
        q = quantity.RateCoefficient(1.0,"cm^6/(mol^2*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*(1e6)**2, 1.0, delta=1e-6)
        self.assertEqual(q.units, "cm^6/(mol^2*s)")

    def test_cm9permol3s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^9/(mol^3*s).
        """
        q = quantity.RateCoefficient(1.0,"cm^9/(mol^3*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*(1e6)**3, 1.0, delta=1e-6)
        self.assertEqual(q.units, "cm^9/(mol^3*s)")

    def test_cm3permolecules(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^3/(molecule*s).
        """
        q = quantity.RateCoefficient(1.0,"cm^3/(molecule*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*1e6/constants.Na, 1.0, delta=1e-6)
        self.assertEqual(q.units, "cm^3/(molecule*s)")

    def test_cm6permolecule2s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^6/(molecule^2*s).
        """
        q = quantity.RateCoefficient(1.0,"cm^6/(molecule^2*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*(1e6/constants.Na)**2, 1.0, delta=1e-6)
        self.assertEqual(q.units, "cm^6/(molecule^2*s)")

    def test_cm9permolecule3s(self):
        """
        Test the creation of a rate coefficient quantity with units of cm^9/(molecule^3*s).
        """
        q = quantity.RateCoefficient(1.0,"cm^9/(molecule^3*s)")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si*(1e6/constants.Na)**3, 1.0, delta=1e-6)
        self.assertEqual(q.units, "cm^9/(molecule^3*s)")

################################################################################

class TestTemperature(unittest.TestCase):
    """
    Contains unit tests of the Temperature unit type object.
    """
            
    def test_K(self):
        """
        Test the creation of a temperature quantity with units of K.
        """
        q = quantity.Temperature(1.0,"K")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "K")

    def test_degC(self):
        """
        Test the creation of a temperature quantity with units of degrees C.
        """
        q = quantity.Temperature(1.0,"degC")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "degC")

    def test_degF(self):
        """
        Test the creation of a temperature quantity with units of degrees F.
        """
        q = quantity.Temperature(1.0,"degF")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 5.0/9.0, delta=1e-6)
        self.assertEqual(q.units, "degF")

    def test_degR(self):
        """
        Test the creation of a temperature quantity with units of degrees R.
        """
        q = quantity.Temperature(1.0,"degR")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 5.0/9.0, delta=1e-6)
        self.assertEqual(q.units, "degR")

################################################################################

class TestTime(unittest.TestCase):
    """
    Contains unit tests of the Time unit type object.
    """
            
    def test_s(self):
        """
        Test the creation of a time quantity with units of s.
        """
        q = quantity.Time(1.0,"s")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "s")
        
    def test_ms(self):
        """
        Test the creation of a time quantity with units of ms.
        """
        q = quantity.Time(1.0,"ms")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-3, delta=1e-9)
        self.assertEqual(q.units, "ms")

    def test_us(self):
        """
        Test the creation of a time quantity with units of us.
        """
        q = quantity.Time(1.0,"us")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-6, delta=1e-12)
        self.assertEqual(q.units, "us")

    def test_ns(self):
        """
        Test the creation of a time quantity with units of ns.
        """
        q = quantity.Time(1.0,"ns")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-9, delta=1e-15)
        self.assertEqual(q.units, "ns")

    def test_ps(self):
        """
        Test the creation of a time quantity with units of ps.
        """
        q = quantity.Time(1.0,"ps")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-12, delta=1e-18)
        self.assertEqual(q.units, "ps")

    def test_fs(self):
        """
        Test the creation of a time quantity with units of fs.
        """
        q = quantity.Time(1.0,"fs")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-15, delta=1e-21)
        self.assertEqual(q.units, "fs")

    def test_min(self):
        """
        Test the creation of a time quantity with units of min.
        """
        q = quantity.Time(1.0,"min")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 60.0, delta=1e-6)
        self.assertEqual(q.units, "min")

    def test_hr(self):
        """
        Test the creation of a time quantity with units of hr.
        """
        q = quantity.Time(1.0,"hr")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 3600.0, delta=1e-6)
        self.assertEqual(q.units, "hr")

################################################################################

class TestVelocity(unittest.TestCase):
    """
    Contains unit tests of the Velocity unit type object.
    """
            
    def test_mpers(self):
        """
        Test the creation of an velocity quantity with units of m/s.
        """
        q = quantity.Velocity(1.0,"m/s")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m/s")

    def test_cmpers(self):
        """
        Test the creation of an velocity quantity with units of m/s.
        """
        q = quantity.Velocity(1.0,"cm/s")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 0.01, delta=1e-8)
        self.assertEqual(q.units, "cm/s")

################################################################################

class TestVolume(unittest.TestCase):
    """
    Contains unit tests of the Volume unit type object.
    """
            
    def test_m3(self):
        """
        Test the creation of an volume quantity with units of m^3.
        """
        q = quantity.Volume(1.0,"m^3")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0, delta=1e-6)
        self.assertEqual(q.units, "m^3")
        
    def test_L(self):
        """
        Test the creation of an volume quantity with units of L.
        """
        q = quantity.Volume(1.0,"L")
        self.assertAlmostEqual(q.value, 1.0, 6)
        self.assertAlmostEqual(q.value_si, 1.0e-3, delta=1e-9)
        self.assertEqual(q.units, "L")
