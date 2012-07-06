#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.thermo module.
"""

import unittest

from rmgpy.thermo import *

################################################################################

class TestThermoModel(unittest.TestCase):
    """
    Contains unit tests of the ThermoModel class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.thermo = ThermoModel(Tmin=(300,"K"), Tmax=(2000,"K"))
        
    def testTemperatureRange(self):
        """
        Test that the temperature range is set and handled appropriately.
        """
        self.assertEqual(self.thermo.Tmin.value, 300)
        self.assertEqual(self.thermo.Tmin.units, "K")
        self.assertEqual(self.thermo.Tmax.value, 2000)
        self.assertEqual(self.thermo.Tmax.units, "K")
        self.assertFalse(self.thermo.isTemperatureValid(200))
        self.assertTrue(self.thermo.isTemperatureValid(300))
        self.assertTrue(self.thermo.isTemperatureValid(400))
        self.assertTrue(self.thermo.isTemperatureValid(500))
        self.assertTrue(self.thermo.isTemperatureValid(600))
        self.assertTrue(self.thermo.isTemperatureValid(800))
        self.assertTrue(self.thermo.isTemperatureValid(1000))
        self.assertTrue(self.thermo.isTemperatureValid(1500))
        self.assertTrue(self.thermo.isTemperatureValid(2000))
        self.assertFalse(self.thermo.isTemperatureValid(2500))
        
    def testPickle(self):
        """
        Test that a ThermoModel object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(self.thermo))
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)
    
    def testOutput(self):
        """
        Test that we can reconstruct a ThermoModel object from its repr()
        output with no loss of information.
        """
        exec('thermo = {0!r}'.format(self.thermo))
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)
        
################################################################################

class TestThermoData(unittest.TestCase):
    """
    Contains unit tests of the ThermoData class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.thermo = ThermoData(
            Tdata = ([300.0,400.0,500.0,600.0,800.0,1000.0,1500.0],"K"),
            Cpdata = ([3.0,4.0,5.0,6.0,8.0,10.0,15.0],"J/(mol*K)"),
            H298 = (-2.0,"kJ/mol"),
            S298 = (50.0,"J/(mol*K)"),
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
    
    def testHeatCapacity(self):
        """
        Test the ThermoData.getHeatCapacity() method.
        """
        self.assertEqual(self.thermo.getHeatCapacity(300), 3)
        self.assertEqual(self.thermo.getHeatCapacity(350), 3.5)
        self.assertEqual(self.thermo.getHeatCapacity(400), 4)
        self.assertEqual(self.thermo.getHeatCapacity(450), 4.5)
        self.assertEqual(self.thermo.getHeatCapacity(500), 5)
        self.assertEqual(self.thermo.getHeatCapacity(550), 5.5)
        self.assertEqual(self.thermo.getHeatCapacity(600), 6)
        self.assertEqual(self.thermo.getHeatCapacity(700), 7)
        self.assertEqual(self.thermo.getHeatCapacity(800), 8)
        self.assertEqual(self.thermo.getHeatCapacity(900), 9)
        self.assertEqual(self.thermo.getHeatCapacity(1000), 10)
        self.assertEqual(self.thermo.getHeatCapacity(1100), 11)
        self.assertEqual(self.thermo.getHeatCapacity(1200), 12)
        self.assertEqual(self.thermo.getHeatCapacity(1300), 13)
        self.assertEqual(self.thermo.getHeatCapacity(1400), 14)
        self.assertEqual(self.thermo.getHeatCapacity(1500), 15)
        self.assertEqual(self.thermo.getHeatCapacity(1600), 15)
        self.assertEqual(self.thermo.getHeatCapacity(1700), 15)
        self.assertEqual(self.thermo.getHeatCapacity(1800), 15)
        self.assertEqual(self.thermo.getHeatCapacity(1900), 15)
        self.assertEqual(self.thermo.getHeatCapacity(2000), 15)
    
    def testEnthalpy(self):
        """
        Test the ThermoData.getEnthalpy() method.
        """
        self.assertEqual(self.thermo.getEnthalpy(300), -1994.0)
        self.assertEqual(self.thermo.getEnthalpy(400), -1644.0)
        self.assertEqual(self.thermo.getEnthalpy(500), -1194.0)
        self.assertEqual(self.thermo.getEnthalpy(600), -644.0)
        self.assertEqual(self.thermo.getEnthalpy(800), 756.0)
        self.assertEqual(self.thermo.getEnthalpy(1000), 2556.0)
        self.assertEqual(self.thermo.getEnthalpy(1500), 8806.0)
        self.assertEqual(self.thermo.getEnthalpy(2000), 16306.0)
        
    def testEntropy(self):
        """
        Test the ThermoData.getEntropy() method.
        """
        self.assertAlmostEqual(self.thermo.getEntropy(300), 50.02, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(400), 51.02, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(500), 52.02, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(600), 53.02, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(800), 55.02, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(1000), 57.02, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(1500), 62.02, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(2000), 66.3352, 3)

    def testFreeEnergy(self):
        """
        Test the ThermoData.getFreeEnergy() method.
        """
        for T in [300,400,500,600,800,1000,1500,2000]:
            self.assertAlmostEqual(self.thermo.getFreeEnergy(T), self.thermo.getEnthalpy(T) - T * self.thermo.getEntropy(T), 9)
    
    def testPickle(self):
        """
        Test that a ThermoData object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(self.thermo))
        self.assertEqual(len(self.thermo.Tdata.values), len(thermo.Tdata.values))
        for T0, T in zip(self.thermo.Tdata.values, thermo.Tdata.values):
            self.assertEqual(T0, T)
        self.assertEqual(self.thermo.Tdata.units, thermo.Tdata.units)
        self.assertEqual(len(self.thermo.Cpdata.values), len(thermo.Cpdata.values))
        for Cp0, Cp in zip(self.thermo.Cpdata.values, thermo.Cpdata.values):
            self.assertEqual(Cp0, Cp)
        self.assertEqual(self.thermo.Cpdata.units, thermo.Cpdata.units)
        self.assertEqual(self.thermo.H298.value, thermo.H298.value)
        self.assertEqual(self.thermo.H298.units, thermo.H298.units)
        self.assertEqual(self.thermo.S298.value, thermo.S298.value)
        self.assertEqual(self.thermo.S298.units, thermo.S298.units)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

    def testOutput(self):
        """
        Test that we can reconstruct a ThermoData object from its repr()
        output with no loss of information.
        """
        exec('thermo = {0!r}'.format(self.thermo))
        self.assertEqual(len(self.thermo.Tdata.values), len(thermo.Tdata.values))
        for T0, T in zip(self.thermo.Tdata.values, thermo.Tdata.values):
            self.assertEqual(T0, T)
        self.assertEqual(self.thermo.Tdata.units, thermo.Tdata.units)
        self.assertEqual(len(self.thermo.Cpdata.values), len(thermo.Cpdata.values))
        for Cp0, Cp in zip(self.thermo.Cpdata.values, thermo.Cpdata.values):
            self.assertEqual(Cp0, Cp)
        self.assertEqual(self.thermo.Cpdata.units, thermo.Cpdata.units)
        self.assertEqual(self.thermo.H298.value, thermo.H298.value)
        self.assertEqual(self.thermo.H298.units, thermo.H298.units)
        self.assertEqual(self.thermo.S298.value, thermo.S298.value)
        self.assertEqual(self.thermo.S298.units, thermo.S298.units)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

################################################################################

class TestWilhoit(unittest.TestCase):
    """
    Contains unit tests of the ThermoData class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.thermo = Wilhoit(
            cp0 = (4.0*8.314472,"J/(mol*K)"), 
            cpInf = (21.0*8.314472,"J/(mol*K)"), 
            a0 = -3.95, 
            a1 = 9.26, 
            a2 = -15.6, 
            a3 = 8.55, 
            B = (500.0,"K"), 
            H0 = (-6.151e+04,"J/mol"), 
            S0 = (-790.2,"J/(mol*K)"),
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
    
    def testHeatCapacity(self):
        """
        Test the Wilhoit.getHeatCapacity() method.
        """
        self.assertAlmostEqual(self.thermo.getHeatCapacity(200), 64.398, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(400), 94.765, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(600), 116.464, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(800), 131.392, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1000), 141.658, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1200), 148.830, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1400), 153.948, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1600), 157.683, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1800), 160.469, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(2000), 162.589, 2)
    
    def testEnthalpy(self):
        """
        Test the Wilhoit.getEnthalpy() method.
        """
        self.assertAlmostEqual(self.thermo.getEnthalpy(200) / 1000., -166.312, 2)
        self.assertAlmostEqual(self.thermo.getEnthalpy(400) / 1000., -150.244, 2)
        self.assertAlmostEqual(self.thermo.getEnthalpy(600) / 1000., -128.990, 2)
        self.assertAlmostEqual(self.thermo.getEnthalpy(800) / 1000., -104.110, 2)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1000) / 1000., -76.7429, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1200) / 1000., -47.6526, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1400) / 1000., -17.3471, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1600) / 1000., 13.8348, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1800) / 1000., 45.6630, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(2000) / 1000., 77.9781, 3)
        
    def testEntropy(self):
        """
        Test the Wilhoit.getEntropy() method.
        """
        self.assertAlmostEqual(self.thermo.getEntropy(200), 287.421, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(400), 341.892, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(600), 384.685, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(800), 420.369, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(1000), 450.861, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(1200), 477.360, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(1400), 500.708, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(1600), 521.521, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(1800), 540.262, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(2000), 557.284, 2)

    def testFreeEnergy(self):
        """
        Test the Wilhoit.getFreeEnergy() method.
        """
        for T in [300,400,500,600,800,1000,1500,2000]:
            self.assertAlmostEqual(self.thermo.getFreeEnergy(T), self.thermo.getEnthalpy(T) - T * self.thermo.getEntropy(T), 9)
    
    def testPickle(self):
        """
        Test that a ThermoData object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(self.thermo))
        self.assertAlmostEqual(self.thermo.cp0.value, thermo.cp0.value, 4)
        self.assertEqual(self.thermo.cp0.units, thermo.cp0.units)
        self.assertAlmostEqual(self.thermo.cpInf.value, thermo.cpInf.value, 3)
        self.assertEqual(self.thermo.cpInf.units, thermo.cpInf.units)
        self.assertAlmostEqual(self.thermo.a0, thermo.a0, 4)
        self.assertAlmostEqual(self.thermo.a1, thermo.a1, 4)
        self.assertAlmostEqual(self.thermo.a2, thermo.a2, 4)
        self.assertAlmostEqual(self.thermo.a3, thermo.a3, 4)
        self.assertAlmostEqual(self.thermo.H0.value, thermo.H0.value, 4)
        self.assertEqual(self.thermo.H0.units, thermo.H0.units)
        self.assertAlmostEqual(self.thermo.S0.value, thermo.S0.value, 4)
        self.assertEqual(self.thermo.S0.units, thermo.S0.units)
        self.assertAlmostEqual(self.thermo.B.value, thermo.B.value, 4)
        self.assertEqual(self.thermo.B.units, thermo.B.units)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

    def testOutput(self):
        """
        Test that we can reconstruct a ThermoData object from its repr()
        output with no loss of information.
        """
        exec('thermo = {0!r}'.format(self.thermo))
        self.assertAlmostEqual(self.thermo.cp0.value, thermo.cp0.value, 4)
        self.assertEqual(self.thermo.cp0.units, thermo.cp0.units)
        self.assertAlmostEqual(self.thermo.cpInf.value, thermo.cpInf.value, 3)
        self.assertEqual(self.thermo.cpInf.units, thermo.cpInf.units)
        self.assertAlmostEqual(self.thermo.a0, thermo.a0, 4)
        self.assertAlmostEqual(self.thermo.a1, thermo.a1, 4)
        self.assertAlmostEqual(self.thermo.a2, thermo.a2, 4)
        self.assertAlmostEqual(self.thermo.a3, thermo.a3, 4)
        self.assertAlmostEqual(self.thermo.H0.value, thermo.H0.value, 4)
        self.assertEqual(self.thermo.H0.units, thermo.H0.units)
        self.assertAlmostEqual(self.thermo.S0.value, thermo.S0.value, 4)
        self.assertEqual(self.thermo.S0.units, thermo.S0.units)
        self.assertAlmostEqual(self.thermo.B.value, thermo.B.value, 4)
        self.assertEqual(self.thermo.B.units, thermo.B.units)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

################################################################################

class TestNASA(unittest.TestCase):
    """
    Contains unit tests of the NASA class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.thermo = NASA(
            coeffs = [1.0e0,1.0e-3,1.0e-6,1.0e-9,1.0e-12,1.0e3,1.0e0],
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
    
    def testHeatCapacity(self):
        """
        Test the NASA.getHeatCapacity() method.
        """
        self.assertAlmostEqual(self.thermo.getHeatCapacity(200), 10.3898, 3)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(400), 13.7156, 3)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(600), 19.1698, 3)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(800), 27.9499, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1000), 41.5724, 3)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1200), 61.8730, 3)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1400), 91.0069, 3)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1600), 131.448, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(1800), 185.991, 2)
        self.assertAlmostEqual(self.thermo.getHeatCapacity(2000), 257.749, 2)
    
    def testEnthalpy(self):
        """
        Test the NASA.getEnthalpy() method.
        """
        self.assertAlmostEqual(self.thermo.getEnthalpy(200) / 1000., 10.1697, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(400) / 1000., 12.5530, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(600) / 1000., 15.7971, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(800) / 1000., 20.4420, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1000) / 1000., 27.2992, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1200) / 1000., 37.5154, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1400) / 1000., 52.6365, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1600) / 1000., 74.6713, 3)
        self.assertAlmostEqual(self.thermo.getEnthalpy(1800) / 1000., 106.155, 2)
        self.assertAlmostEqual(self.thermo.getEnthalpy(2000) / 1000., 150.215, 2)
        
    def testEntropy(self):
        """
        Test the NASA.getEntropy() method.
        """
        self.assertAlmostEqual(self.thermo.getEntropy(200), 54.2219, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(400), 62.3519, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(600), 68.8549, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(800), 75.4761, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(1000), 83.0706, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(1200), 92.3279, 3)
        self.assertAlmostEqual(self.thermo.getEntropy(1400), 103.925, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(1600), 118.577, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(1800), 137.055, 2)
        self.assertAlmostEqual(self.thermo.getEntropy(2000), 160.200, 2)

    def testFreeEnergy(self):
        """
        Test the NASA.getFreeEnergy() method.
        """
        for T in [300,400,500,600,800,1000,1500,2000]:
            self.assertAlmostEqual(self.thermo.getFreeEnergy(T), self.thermo.getEnthalpy(T) - T * self.thermo.getEntropy(T), 9)
    
    def testPickle(self):
        """
        Test that a NASA object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(self.thermo))
        self.assertAlmostEqual(self.thermo.cm2, thermo.cm2, 4)
        self.assertAlmostEqual(self.thermo.cm1, thermo.cm1, 4)
        self.assertAlmostEqual(self.thermo.c0 / thermo.c0, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c1 / thermo.c1, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c2 / thermo.c2, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c3 / thermo.c3, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c4 / thermo.c4, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c5 / thermo.c5, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c6 / thermo.c6, 1.0, 4)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

    def testOutput(self):
        """
        Test that we can reconstruct a NASA object from its repr()
        output with no loss of information.
        """
        exec('thermo = {0!r}'.format(self.thermo))
        self.assertAlmostEqual(self.thermo.cm2, thermo.cm2, 4)
        self.assertAlmostEqual(self.thermo.cm1, thermo.cm1, 4)
        self.assertAlmostEqual(self.thermo.c0 / thermo.c0, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c1 / thermo.c1, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c2 / thermo.c2, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c3 / thermo.c3, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c4 / thermo.c4, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c5 / thermo.c5, 1.0, 4)
        self.assertAlmostEqual(self.thermo.c6 / thermo.c6, 1.0, 4)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

################################################################################

class TestMultiNASA(unittest.TestCase):
    """
    Contains unit tests of the MultiNASA class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        nasa0 = NASA(coeffs=[11.0,12.0,13.0,14.0,15.0,16.0,17.0], Tmin=(300.0,"K"), Tmax=(1000.0,"K"), comment='This data is completely made up and unphysical')
        nasa1 = NASA(coeffs=[21.0,22.0,23.0,24.0,25.0,26.0,27.0], Tmin=(1000.0,"K"), Tmax=(6000.0,"K"), comment='This data is also completely made up and unphysical')
        self.thermo = MultiNASA(
            polynomials=[nasa0, nasa1],
            Tmin = (300.0,"K"), 
            Tmax = (6000.0,"K"), 
            comment = """This data is completely made up and unphysical""",
        )
    
    def testPickle(self):
        """
        Test that a MultiNASA object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(self.thermo))
        self.assertEqual(len(self.thermo.polynomials), len(thermo.polynomials))
        for poly0, poly in zip(self.thermo.polynomials, thermo.polynomials):
            self.assertEqual(poly0.cm2, poly.cm2)
            self.assertEqual(poly0.cm1, poly.cm1)
            self.assertEqual(poly0.c0, poly.c0)
            self.assertEqual(poly0.c1, poly.c1)
            self.assertEqual(poly0.c2, poly.c2)
            self.assertEqual(poly0.c3, poly.c3)
            self.assertEqual(poly0.c4, poly.c4)
            self.assertEqual(poly0.c5, poly.c5)
            self.assertEqual(poly0.c6, poly.c6)
            self.assertEqual(poly0.Tmin.value, poly.Tmin.value)
            self.assertEqual(poly0.Tmin.units, poly.Tmin.units)
            self.assertEqual(poly0.Tmax.value, poly.Tmax.value)
            self.assertEqual(poly0.Tmax.units, poly.Tmax.units)
            self.assertEqual(poly0.comment, poly.comment)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

    def testOutput(self):
        """
        Test that we can reconstruct a MultiNASA object from its repr()
        output with no loss of information.
        """
        exec('thermo = {0!r}'.format(self.thermo))
        self.assertEqual(len(self.thermo.polynomials), len(thermo.polynomials))
        for poly0, poly in zip(self.thermo.polynomials, thermo.polynomials):
            self.assertEqual(poly0.cm2, poly.cm2)
            self.assertEqual(poly0.cm1, poly.cm1)
            self.assertEqual(poly0.c0, poly.c0)
            self.assertEqual(poly0.c1, poly.c1)
            self.assertEqual(poly0.c2, poly.c2)
            self.assertEqual(poly0.c3, poly.c3)
            self.assertEqual(poly0.c4, poly.c4)
            self.assertEqual(poly0.c5, poly.c5)
            self.assertEqual(poly0.c6, poly.c6)
            self.assertEqual(poly0.Tmin.value, poly.Tmin.value)
            self.assertEqual(poly0.Tmin.units, poly.Tmin.units)
            self.assertEqual(poly0.Tmax.value, poly.Tmax.value)
            self.assertEqual(poly0.Tmax.units, poly.Tmax.units)
            self.assertEqual(poly0.comment, poly.comment)
        self.assertEqual(self.thermo.Tmin.value, thermo.Tmin.value)
        self.assertEqual(self.thermo.Tmin.units, thermo.Tmin.units)
        self.assertEqual(self.thermo.Tmax.value, thermo.Tmax.value)
        self.assertEqual(self.thermo.Tmax.units, thermo.Tmax.units)
        self.assertEqual(self.thermo.comment, thermo.comment)

################################################################################

class TestConversion(unittest.TestCase):
    """
    Contains unit tests involving conversion from one thermodynamics model to
    another.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.thermoData = ThermoData(
            Tdata=([303.231,385.553,467.875,550.197,632.519,714.841,797.163,879.485,961.807,1044.13,1126.45,1208.77,1291.09,1373.42,1455.74,1538.06,1620.38,1702.7,1785.03,1867.35],"K"), 
            Cpdata=([50.4468,57.4605,64.3653,70.7093,76.4,81.4756,86.0017,90.0394,93.6412,96.853,99.7159,102.268,104.543,106.574,108.388,110.011,111.466,112.772,113.947,115.007],"J/(mol*K)"), 
            H298=(12.3533,"kJ/mol"), 
            S298=(261.187,"J/(mol*K)"),
        )
        self.wilhoit = Wilhoit(
            cp0 = (4.0*8.314472,"J/(mol*K)"), 
            cpInf = (21.0*8.314472,"J/(mol*K)"), 
            a0 = -3.95, 
            a1 = 9.26, 
            a2 = -15.6, 
            a3 = 8.55, 
            B = (500.0,"K"), 
            H0 = (-6.151e+04,"J/mol"), 
            S0 = (-790.2,"J/(mol*K)"),
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        nasa0 = NASA(coeffs=[0.93355381E+00, 0.26424579E-01, 0.61059727E-05,-0.21977499E-07, 0.95149253E-11,-0.13958520E+05, 0.19201691E+02], Tmin=(298.0,"K"), Tmax=(1000.0,"K"), comment='Low-temperature polynomial for C3H8 from GRI-Mech 3.0')
        nasa1 = NASA(coeffs=[0.75341368E+01, 0.18872239E-01,-0.62718491E-05, 0.91475649E-09,-0.47838069E-13,-0.16467516E+05,-0.17892349E+02], Tmin=(1000.0,"K"), Tmax=(5000.0,"K"), comment='High-temperature polynomial for C3H8 from GRI-Mech 3.0')
        self.multiNASA = MultiNASA(
            polynomials=[nasa0, nasa1],
            Tmin = (298.0,"K"), 
            Tmax = (5000.0,"K"), 
            comment = """C3H8 from GRI-Mech 3.0""",
        )
    
    def testThermoDataToWilhoit(self):
        """
        Test that a ThermoData object can be successfully converted to a
        Wilhoit object.
        """
        thermo0 = self.thermoData
        thermo = convertThermoModel(thermo0, Wilhoit, linear=False, nFreq=11, nRotors=1)
        self.assertNotEqual(thermo, None)
        for T in [300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(thermo0.getHeatCapacity(T) / thermo.getHeatCapacity(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEnthalpy(T) / thermo.getEnthalpy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEntropy(T) / thermo.getEntropy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getFreeEnergy(T) / thermo.getFreeEnergy(T), 1.0, 0)
        
    def testWilhoitToThermoData(self):
        """
        Test that a Wilhoit object can be successfully converted to a
        ThermoData object.
        """
        thermo0 = self.wilhoit
        thermo = convertThermoModel(thermo0, ThermoData, Tdata=[300,400,500,600,800,1000,1500])
        self.assertNotEqual(thermo, None)
        for T in [300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(thermo0.getHeatCapacity(T) / thermo.getHeatCapacity(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEnthalpy(T) / thermo.getEnthalpy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEntropy(T) / thermo.getEntropy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getFreeEnergy(T) / thermo.getFreeEnergy(T), 1.0, 0)
    
    def testThermoDataToMultiNASA(self):
        """
        Test that a ThermoData object can be successfully converted to a
        MultiNASA object.
        """
        thermo0 = self.thermoData
        thermo = convertThermoModel(thermo0, MultiNASA, linear=False, nFreq=11, nRotors=1, Tmin=300, Tmax=6000, Tint=1000)
        self.assertNotEqual(thermo, None)
        for T in [300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(thermo0.getHeatCapacity(T) / thermo.getHeatCapacity(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEnthalpy(T) / thermo.getEnthalpy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEntropy(T) / thermo.getEntropy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getFreeEnergy(T) / thermo.getFreeEnergy(T), 1.0, 0)
    
    def testMultiNASAToThermoData(self):
        """
        Test that a MultiNASA object can be successfully converted to a
        ThermoData object.
        """
        thermo0 = self.multiNASA
        thermo = convertThermoModel(thermo0, ThermoData, Tdata=[300,400,500,600,800,1000,1500])
        self.assertNotEqual(thermo, None)
        for T in [300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(thermo0.getHeatCapacity(T) / thermo.getHeatCapacity(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEnthalpy(T) / thermo.getEnthalpy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEntropy(T) / thermo.getEntropy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getFreeEnergy(T) / thermo.getFreeEnergy(T), 1.0, -1)
    
    def testWilhoitToMultiNASA(self):
        """
        Test that a Wilhoit object can be successfully converted to a
        MultiNASA object.
        """
        thermo0 = self.wilhoit
        thermo = convertThermoModel(thermo0, MultiNASA, Tmin=300, Tmax=6000, Tint=1000)
        self.assertNotEqual(thermo, None)
        for T in [300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(thermo0.getHeatCapacity(T) / thermo.getHeatCapacity(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEnthalpy(T) / thermo.getEnthalpy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEntropy(T) / thermo.getEntropy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getFreeEnergy(T) / thermo.getFreeEnergy(T), 1.0, 0)
    
    def testMultiNASAToWilhoit(self):
        """
        Test that a MultiNASA object can be successfully converted to a
        Wilhoit object.
        """
        thermo0 = self.multiNASA
        thermo = convertThermoModel(thermo0, Wilhoit, linear=False, nFreq=25, nRotors=2)
        self.assertNotEqual(thermo, None)
        for T in [300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(thermo0.getHeatCapacity(T) / thermo.getHeatCapacity(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEnthalpy(T) / thermo.getEnthalpy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getEntropy(T) / thermo.getEntropy(T), 1.0, 0)
            self.assertAlmostEqual(thermo0.getFreeEnergy(T) / thermo.getFreeEnergy(T), 1.0, 0)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
