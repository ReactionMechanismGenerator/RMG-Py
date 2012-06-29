#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.kinetics module.
"""

import numpy
import unittest

from rmgpy.kinetics import *

################################################################################

class TestKineticsModel(unittest.TestCase):
    """
    Contains unit tests of the KineticsModel class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.kinetics = KineticsModel(Tmin=(300,"K"), Tmax=(2000,"K"), Pmin=(0.01,"bar"), Pmax=(100,"bar"))
        
    def testTemperatureRange(self):
        """
        Test that the temperature range is set and handled appropriately.
        """
        self.assertEqual(self.kinetics.Tmin.value, 300)
        self.assertEqual(self.kinetics.Tmin.units, "K")
        self.assertEqual(self.kinetics.Tmax.value, 2000)
        self.assertEqual(self.kinetics.Tmax.units, "K")
        self.assertFalse(self.kinetics.isTemperatureValid(200))
        self.assertTrue(self.kinetics.isTemperatureValid(300))
        self.assertTrue(self.kinetics.isTemperatureValid(400))
        self.assertTrue(self.kinetics.isTemperatureValid(500))
        self.assertTrue(self.kinetics.isTemperatureValid(600))
        self.assertTrue(self.kinetics.isTemperatureValid(800))
        self.assertTrue(self.kinetics.isTemperatureValid(1000))
        self.assertTrue(self.kinetics.isTemperatureValid(1500))
        self.assertTrue(self.kinetics.isTemperatureValid(2000))
        self.assertFalse(self.kinetics.isTemperatureValid(2500))
        
    def testPressureRange(self):
        """
        Test that the pressure range is set and handled appropriately.
        """
        self.assertEqual(self.kinetics.Pmin.value, 1e3)
        self.assertEqual(self.kinetics.Pmin.units, "bar")
        self.assertEqual(self.kinetics.Pmax.value, 1e7)
        self.assertEqual(self.kinetics.Pmax.units, "bar")
        self.assertFalse(self.kinetics.isPressureValid(1e2))
        self.assertTrue(self.kinetics.isPressureValid(1e3))
        self.assertTrue(self.kinetics.isPressureValid(1e4))
        self.assertTrue(self.kinetics.isPressureValid(1e5))
        self.assertTrue(self.kinetics.isPressureValid(1e6))
        self.assertTrue(self.kinetics.isPressureValid(1e7))
        self.assertFalse(self.kinetics.isPressureValid(1e8))
        
    def testPickle(self):
        """
        Test that a KineticsModel object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
    
#    def testOutput(self):
#        """
#        Test that we can reconstruct a KineticsModel object from its repr()
#        output with no loss of information.
#        """
#        # KineticsModel does not have repr since it must be from derived class.
#        exec('kinetics = {0!r}'.format(self.kinetics))
#        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
#        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
#        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
#        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
#        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
#        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
#        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
#        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
#        self.assertEqual(self.kinetics.comment, kinetics.comment)
        
################################################################################

class TestKineticsData(unittest.TestCase):
    """
    Contains unit tests of the KineticsData class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        Tlist = [300.0,400.0,500.0,600.0,800.0,1000.0,1500.0]
        self.kinetics = KineticsData(
            Tdata = (Tlist,"K"),
            kdata = ([10**(1000/T) for T in Tlist],"m^3/(mol*s)"),
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
    
    def testIsPressureDependent(self):
        """
        Test the KineticsData.isPressureDependent() method.
        
        """
        self.assertFalse(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the KineticsData.getRateCoefficient() method.
        """
        for T in xrange(300,1500,50):
            self.assertAlmostEqual(self.kinetics.getRateCoefficient(T), 10**(1000.0/T),4)
    
    def testPickle(self):
        """
        Test that a KineticsData object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        self.assertEqual(len(self.kinetics.Tdata.values), len(kinetics.Tdata.values))
        for T0, T in zip(self.kinetics.Tdata.values, kinetics.Tdata.values):
            self.assertEqual(T0, T)
        self.assertEqual(self.kinetics.Tdata.units, kinetics.Tdata.units)
        self.assertEqual(len(self.kinetics.kdata.values), len(kinetics.kdata.values))
        for k0, k in zip(self.kinetics.kdata.values, kinetics.kdata.values):
            self.assertEqual(k0, k)
        self.assertEqual(self.kinetics.kdata.units, kinetics.kdata.units)
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testOutput(self):
        """
        Test that we can reconstruct a KineticsData object from its repr()
        output with no loss of information.
        """
        exec('kinetics = {0!r}'.format(self.kinetics))
        self.assertEqual(len(self.kinetics.Tdata.values), len(kinetics.Tdata.values))
        for T0, T in zip(self.kinetics.Tdata.values, kinetics.Tdata.values):
            self.assertEqual(T0, T)
        self.assertEqual(self.kinetics.Tdata.units, kinetics.Tdata.units)
        self.assertEqual(len(self.kinetics.kdata.values), len(kinetics.kdata.values))
        for k0, k in zip(self.kinetics.kdata.values, kinetics.kdata.values):
            self.assertAlmostEqual(k0, k, 2)
        self.assertEqual(self.kinetics.kdata.units, kinetics.kdata.units)
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
    
################################################################################

class TestArrhenius(unittest.TestCase):
    """
    Contains unit tests of the Arrhenius class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.kinetics = Arrhenius(
            A = (1.0e6,"cm^3/(mol*s)"),
            n = 1,
            Ea = (20,"kJ/mol"),
            T0 = (300,"K"),
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
    
    def testIsPressureDependent(self):
        """
        Test the Arrhenius.isPressureDependent() method.
        
        """
        self.assertFalse(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the Arrhenius.getRateCoefficient() method.
        """
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(300) / 3.2943e-4, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(500) / 1.3568e-2, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1000) / 3.0075e-1, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1500) / 1.0058e0, 1, 3)
    
    def testChangeT0(self):
        """
        Test the Arrhenius.changeT0() method.
        """
        Tlist = [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]
        k0list = [self.kinetics.getRateCoefficient(T) for T in Tlist]
        self.kinetics.changeT0(1)
        self.assertEqual(self.kinetics.T0.value, 1)
        self.assertEqual(self.kinetics.T0.units, "K")
        for T, k0 in zip(Tlist, k0list):
            self.assertAlmostEqual(k0 / self.kinetics.getRateCoefficient(T), 1, 6)
        
    def testFitToData(self):
        """
        Test the Arrhenius.fitToData() method.
        """
        Tdata = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500], numpy.float)
        kdata = numpy.array([self.kinetics.getRateCoefficient(T) for T in Tdata], numpy.float)
        kinetics = Arrhenius().fitToData(Tdata, kdata, kunits="m^3/(mol*s)", T0=300)
        self.assertEqual(self.kinetics.T0.value, 300)
        self.assertEqual(self.kinetics.T0.units, "K")
        for T, k in zip(Tdata, kdata):
            self.assertAlmostEqual(k / kinetics.getRateCoefficient(T), 1, 6)
        self.assertAlmostEqual(kinetics.A.value / self.kinetics.A.value, 1, 4)
        self.assertAlmostEqual(kinetics.n.value / self.kinetics.n.value, 1, 4)
        self.assertAlmostEqual(kinetics.Ea.value / self.kinetics.Ea.value, 1, 4)
        self.assertAlmostEqual(kinetics.T0.value / self.kinetics.T0.value, 1, 4)
        
    def testPickle(self):
        """
        Test that an Arrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        self.assertEqual(self.kinetics.A.value, kinetics.A.value)
        self.assertEqual(self.kinetics.A.units, kinetics.A.units)
        self.assertEqual(self.kinetics.n.value, kinetics.n.value)
        self.assertEqual(self.kinetics.n.units, kinetics.n.units)
        self.assertEqual(self.kinetics.T0.value, kinetics.T0.value)
        self.assertEqual(self.kinetics.T0.units, kinetics.T0.units)
        self.assertEqual(self.kinetics.Ea.value, kinetics.Ea.value)
        self.assertEqual(self.kinetics.Ea.units, kinetics.Ea.units)
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testOutput(self):
        """
        Test that we can reconstruct an Arrhenius object from its repr()
        output with no loss of information.
        """
        exec('kinetics = {0!r}'.format(self.kinetics))
        self.assertEqual(self.kinetics.A.value, kinetics.A.value)
        self.assertEqual(self.kinetics.A.units, kinetics.A.units)
        self.assertEqual(self.kinetics.n.value, kinetics.n.value)
        self.assertEqual(self.kinetics.n.units, kinetics.n.units)
        self.assertEqual(self.kinetics.T0.value, kinetics.T0.value)
        self.assertEqual(self.kinetics.T0.units, kinetics.T0.units)
        self.assertEqual(self.kinetics.Ea.value, kinetics.Ea.value)
        self.assertEqual(self.kinetics.Ea.units, kinetics.Ea.units)
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
    
################################################################################

class TestArrheniusEP(unittest.TestCase):
    """
    Contains unit tests of the ArrheniusEP class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.kinetics = ArrheniusEP(
            A = (1.0e6,"cm^3/(mol*s)"),
            n = 1,
            alpha = 0.5,
            E0 = (20,"kJ/mol"),
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        
    def testIsPressureDependent(self):
        """
        Test the ArrheniusEP.isPressureDependent() method.
        
        """
        self.assertFalse(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the ArrheniusEP.getRateCoefficient() method.
        """
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(300, 0) / 9.8829e-2, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(500, 0) / 4.0703e0, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1000, 0) / 9.0225e1, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1500, 0) / 3.0175e2, 1, 3)
    
    def testGetActivationEnergy(self):
        """
        Test the ArrheniusEP.getActivationEnergy() method.
        """
        for dHrxn in [-100000, 0, 100000]:
            Ea0 = self.kinetics.getActivationEnergy(dHrxn)
            Ea = self.kinetics.alpha.value * dHrxn + self.kinetics.E0.value
            self.assertAlmostEqual(Ea0 / Ea, 1, 6)
    
    def TestToArrhenius(self):
        """
        Test the ArrheniusEP.toArrhenius() method.
        """
        for dHrxn in [-100000, 0, 100000]:
            kinetics = self.kinetics.toArrhenius(dHrxn)
            Ea0 = self.kinetics.getActivationEnergy(dHrxn)
            Ea = self.kinetics.alpha.value * dHrxn + self.kinetics.E0.value
            self.assertAlmostEqual(Ea0 / Ea, 1, 6)
            self.assertAlmostEqual(Ea0 / kinetics.Ea, 1, 6)
            self.assertAlmostEqual(self.kinetics.A.value / kinetics.A.value, 1, 6)
            self.assertAlmostEqual(self.kinetics.n.value / kinetics.n.value, 1, 6)
            self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
            self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
            self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
            self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
            self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
            self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
            self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testPickle(self):
        """
        Test that an ArrheniusEP object can be successfully pickled and
        unpickled with no loss of information.
        """

        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        
        self.assertAlmostEqual(self.kinetics.A.value, kinetics.A.value, 4)
        self.assertAlmostEqual(self.kinetics.n.value, kinetics.n.value, 4)
        self.assertAlmostEqual(self.kinetics.alpha.value, kinetics.alpha.value, 4)
        self.assertAlmostEqual(self.kinetics.E0.value, kinetics.E0.value, 4)

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testOutput(self):
        """
        Test that an ArrheniusEP object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        exec('kinetics = {0!r}'.format(self.kinetics))
        
        self.assertAlmostEqual(self.kinetics.A.value, kinetics.A.value, 4)
        self.assertAlmostEqual(self.kinetics.n.value, kinetics.n.value, 4)
        self.assertAlmostEqual(self.kinetics.alpha.value, kinetics.alpha.value, 4)
        self.assertAlmostEqual(self.kinetics.E0.value, kinetics.E0.value, 4)

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

################################################################################

class TestMultiKinetics(unittest.TestCase):
    """
    Contains unit tests of the MultiKinetics class: very simple test with two Arrhenius rates.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        arrh0 = Arrhenius(
            A = (1.0e6,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (10,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        arrh1 = Arrhenius(
            A = (1.0e12,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (20,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        self.kinetics = MultiKinetics(
            kineticsList = [arrh0, arrh1],
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
    
    def testIsPressureDependent(self):
        """
        Test the MultiKinetics.isPressureDependent() method.
        
        """
        self.assertFalse(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the MultiArrhenius.getRateCoefficient() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            k0 = self.kinetics.getRateCoefficient(T)
            k1 = sum([arrh.getRateCoefficient(T) for arrh in self.kinetics.kineticsList])
            self.assertAlmostEqual(k0 / k1, 1, 6)
    
    def testPickle(self):
        """
        Test that a MultiKinetics object can be successfully pickled and
        unpickled with no loss of information.
        """

        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        
        Narrh = 2
        self.assertEqual(len(self.kinetics.kineticsList), Narrh)
        self.assertEqual(len(kinetics.kineticsList), Narrh)
        for i in range(Narrh):
            self.assertEqual(self.kinetics.kineticsList[i].A.value, kinetics.kineticsList[i].A.value)
            self.assertEqual(self.kinetics.kineticsList[i].n.value, kinetics.kineticsList[i].n.value)
            self.assertEqual(self.kinetics.kineticsList[i].T0.value, kinetics.kineticsList[i].T0.value)
            self.assertEqual(self.kinetics.kineticsList[i].Ea.value, kinetics.kineticsList[i].Ea.value)

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
    
    def testOutput(self):
        """
        Test that a MultiKinetics object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        exec('kinetics = {0!r}'.format(self.kinetics))
        
        Narrh = 2
        self.assertEqual(len(self.kinetics.kineticsList), Narrh)
        self.assertEqual(len(kinetics.kineticsList), Narrh)
        for i in range(Narrh):
            self.assertEqual(self.kinetics.kineticsList[i].A.value, kinetics.kineticsList[i].A.value)
            self.assertEqual(self.kinetics.kineticsList[i].n.value, kinetics.kineticsList[i].n.value)
            self.assertEqual(self.kinetics.kineticsList[i].T0.value, kinetics.kineticsList[i].T0.value)
            self.assertEqual(self.kinetics.kineticsList[i].Ea.value, kinetics.kineticsList[i].Ea.value)

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

################################################################################

class TestPDepArrhenius(unittest.TestCase):
    """
    Contains unit tests of the PDepArrhenius class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        arrh0 = Arrhenius(
            A = (1.0e6,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (10,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        arrh1 = Arrhenius(
            A = (1.0e12,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (20,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        self.kinetics = PDepArrhenius(
            pressures = ([0.1, 10.0],"bar"),
            arrhenius = [arrh0, arrh1],
            highPlimit = arrh1,
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            Pmin = (0.1,"bar"), 
            Pmax = (10.0,"bar"), 
            comment = """This data is completely made up""",
        )

    def testIsPressureDependent(self):
        """
        Test the PDepArrhenius.isPressureDependent() method.
        
        """
        self.assertTrue(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the PDepArrhenius.getRateCoefficient() method.
        """
        P = 1e4
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            k0 = self.kinetics.getRateCoefficient(T, P)
            k1 = self.kinetics.arrhenius[0].getRateCoefficient(T)
            self.assertAlmostEqual(k0 / k1, 1, 6)
        P = 1e6
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            k0 = self.kinetics.getRateCoefficient(T, P)
            k1 = self.kinetics.arrhenius[1].getRateCoefficient(T)
            self.assertAlmostEqual(k0 / k1, 1, 6)
        P = 1e5
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            k0 = self.kinetics.getRateCoefficient(T, P)
            k1 = math.sqrt(self.kinetics.arrhenius[0].getRateCoefficient(T) * self.kinetics.arrhenius[1].getRateCoefficient(T))
            self.assertAlmostEqual(k0 / k1, 1, 6)
        
    def testFitToData(self):
        """
        Test the PDepArrhenius.fitToData() method.
        """
        Tdata = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500], numpy.float)
        Pdata = numpy.array([1e4,3e4,1e5,3e5,1e6], numpy.float)
        kdata = numpy.zeros([len(Tdata),len(Pdata)], numpy.float)
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                kdata[t,p] = self.kinetics.getRateCoefficient(Tdata[t], Pdata[p])
        kinetics = PDepArrhenius().fitToData(Tdata, Pdata, kdata, kunits='m^3/(mol*s)')
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                self.assertAlmostEqual(kinetics.getRateCoefficient(Tdata[t], Pdata[p]) / kdata[t,p], 1, 4)
        
    def testPickle(self):
        """
        Test that a PDepArrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """

        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        
        Narrh = 2
        self.assertEqual(len(self.kinetics.pressures.values), Narrh)
        self.assertEqual(len(kinetics.pressures.values), Narrh)
        self.assertEqual(len(self.kinetics.arrhenius), Narrh)
        self.assertEqual(len(kinetics.arrhenius), Narrh)
        for i in range(Narrh):
            self.assertEqual(self.kinetics.pressures.values[i], kinetics.pressures.values[i])
            self.assertEqual(self.kinetics.arrhenius[i].A.value, kinetics.arrhenius[i].A.value)
            self.assertEqual(self.kinetics.arrhenius[i].n.value, kinetics.arrhenius[i].n.value)
            self.assertEqual(self.kinetics.arrhenius[i].T0.value, kinetics.arrhenius[i].T0.value)
            self.assertEqual(self.kinetics.arrhenius[i].Ea.value, kinetics.arrhenius[i].Ea.value)

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
        self.assertEqual(self.kinetics.highPlimit.A.value,  kinetics.highPlimit.A.value)
        self.assertEqual(self.kinetics.highPlimit.n.value,  kinetics.highPlimit.n.value)
        self.assertEqual(self.kinetics.highPlimit.T0.value, kinetics.highPlimit.T0.value)
        self.assertEqual(self.kinetics.highPlimit.Ea.value, kinetics.highPlimit.Ea.value)

    def testOutput(self):
        """
        Test that a PDepArrhenius object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        exec('kinetics = {0!r}'.format(self.kinetics))
        
        Narrh = 2
        self.assertEqual(len(self.kinetics.pressures.values), Narrh)
        self.assertEqual(len(kinetics.pressures.values), Narrh)
        self.assertEqual(len(self.kinetics.arrhenius), Narrh)
        self.assertEqual(len(kinetics.arrhenius), Narrh)
        for i in range(Narrh):
            self.assertEqual(self.kinetics.pressures.values[i], kinetics.pressures.values[i])
            self.assertEqual(self.kinetics.arrhenius[i].A.value, kinetics.arrhenius[i].A.value)
            self.assertEqual(self.kinetics.arrhenius[i].n.value, kinetics.arrhenius[i].n.value)
            self.assertEqual(self.kinetics.arrhenius[i].T0.value, kinetics.arrhenius[i].T0.value)
            self.assertEqual(self.kinetics.arrhenius[i].Ea.value, kinetics.arrhenius[i].Ea.value)

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
        self.assertEqual(self.kinetics.highPlimit.A.value,  kinetics.highPlimit.A.value)
        self.assertEqual(self.kinetics.highPlimit.n.value,  kinetics.highPlimit.n.value)
        self.assertEqual(self.kinetics.highPlimit.T0.value, kinetics.highPlimit.T0.value)
        self.assertEqual(self.kinetics.highPlimit.Ea.value, kinetics.highPlimit.Ea.value)

################################################################################

class TestChebyshev(unittest.TestCase):
    """
    Contains unit tests of the Chebyshev class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        coeffs = numpy.array([[1.0,2.0,3.0], [4.0,5.0,6.0], [7.0,8.0,9.0], [10.0,11.0,12.0]], numpy.float64)
        self.kinetics = Chebyshev(
            coeffs = coeffs,
            kunits = 'm^3/(mol*s)',
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            Pmin = (0.1,"bar"), 
            Pmax = (10.0,"bar"), 
            comment = """This data is completely made up""",
        )
        
    def testIsPressureDependent(self):
        """
        Test the Chebyshev.isPressureDependent() method.
        
        """
        self.assertTrue(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the Chebyshev.getRateCoefficient() method.
        """
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(300, 1e4) / 1e-6, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(300, 1e5) / 1e0, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(300, 1e6) / 1e-18, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(500, 1e4) / 4.9370e-5, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(500, 1e5) / 5.6558e-1, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(500, 1e6) / 1.2034e-13, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1000, 1e4) / 3.1734e-6, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1000, 1e5) / 5.5742e-2, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1000, 1e6) / 3.1958e-17, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1500, 1e4) / 9.4797e11, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1500, 1e5) / 8.31e-6, 1, 3)
        self.assertAlmostEqual(self.kinetics.getRateCoefficient(1500, 1e6) / 8.519e35, 1, 3)
        
    def testFitToData(self):
        """
        Test the Chebyshev.fitToData() method.
        """
        Tdata = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500], numpy.float)
        Pdata = numpy.array([1e4,3e4,1e5,3e5,1e6], numpy.float)
        kdata = numpy.zeros([len(Tdata),len(Pdata)], numpy.float)
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                kdata[t,p] = self.kinetics.getRateCoefficient(Tdata[t], Pdata[p])
        kinetics = Chebyshev().fitToData(Tdata, Pdata, kdata, kunits=self.kinetics.kunits, degreeT=6, degreeP=4, Tmin=300, Tmax=1500, Pmin=1e4, Pmax=1e6)
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                self.assertAlmostEqual(kinetics.getRateCoefficient(Tdata[t], Pdata[p]) / kdata[t,p], 1, 4)
        
    def testPickle(self):
        """
        Test that a Chebyshev object can be successfully pickled and
        unpickled with no loss of information.
        """

        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        
        degreeT = 4; degreeP = 3
        self.assertEqual(self.kinetics.coeffs.shape[0], degreeT)
        self.assertEqual(self.kinetics.coeffs.shape[1], degreeP)
        self.assertEqual(self.kinetics.kunits, kinetics.kunits)
        self.assertEqual(self.kinetics.degreeT, degreeT)
        self.assertEqual(self.kinetics.degreeP, degreeP)
        self.assertEqual(kinetics.coeffs.shape[0], degreeT)
        self.assertEqual(kinetics.coeffs.shape[1], degreeP)
        self.assertEqual(kinetics.degreeT, degreeT)
        self.assertEqual(kinetics.degreeP, degreeP)
        for i in range(degreeT):
            for j in range(degreeP):
                self.assertEqual(self.kinetics.coeffs[i,j], kinetics.coeffs[i,j])
        
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testOutput(self):
        """
        Test that a Chebyshev object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        exec('kinetics = {0!r}'.format(self.kinetics))
        
        degreeT = 4; degreeP = 3
        self.assertEqual(self.kinetics.coeffs.shape[0], degreeT)
        self.assertEqual(self.kinetics.coeffs.shape[1], degreeP)
        self.assertEqual(self.kinetics.kunits, kinetics.kunits)
        self.assertEqual(self.kinetics.degreeT, degreeT)
        self.assertEqual(self.kinetics.degreeP, degreeP)
        self.assertEqual(kinetics.coeffs.shape[0], degreeT)
        self.assertEqual(kinetics.coeffs.shape[1], degreeP)
        self.assertEqual(kinetics.degreeT, degreeT)
        self.assertEqual(kinetics.degreeP, degreeP)
        for i in range(degreeT):
            for j in range(degreeP):
                self.assertEqual(self.kinetics.coeffs[i,j], kinetics.coeffs[i,j])

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
        
################################################################################

class TestThirdBody(unittest.TestCase):
    """
    Contains unit tests of the ThirdBody class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        arrhHigh = Arrhenius(
            A = (1.0e6,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (10,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        efficiencies = {'N#N': 0.5, '[Ar]': 1.5}
        self.kinetics = ThirdBody(
            arrheniusHigh = arrhHigh, 
            efficiencies = efficiencies, 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            Pmin = (0.1,"bar"), 
            Pmax = (10.0,"bar"), 
            comment = """This data is completely made up""",
        )
        
    def testIsPressureDependent(self):
        """
        Test the ThirdBody.isPressureDependent() method.
        
        """
        self.assertTrue(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the ThirdBody.getRateCoefficient() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            for P in [1e4,1e5,1e6]:
                k0 = self.kinetics.getRateCoefficient(T, P)
                k = self.kinetics.arrheniusHigh.getRateCoefficient(T) * (P / 8.314472 / T)
                self.assertAlmostEqual(k / k0, 1, 6)
        
    def testGetColliderEfficiency(self):
        """
        Test the ThirdBody.getColliderEfficiency() method.
        """
        colliders = {}
        for collider in self.kinetics.efficiencies:
            colliders[collider] = 1.0 / len(self.kinetics.efficiencies)
            efficiency = self.kinetics.getColliderEfficiency(collider)
            self.assertEqual(efficiency, self.kinetics.efficiencies[collider])
        
        efficiency = self.kinetics.getColliderEfficiency(colliders)
        efficiency0 = sum(self.kinetics.efficiencies.values()) / len(self.kinetics.efficiencies)
        self.assertEqual(efficiency, efficiency0)
            
    def testPickle(self):
        """
        Test that a ThirdBody object can be successfully pickled and
        unpickled with no loss of information.
        """

        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        
        self.assertEqual(self.kinetics.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(self.kinetics.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(self.kinetics.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(self.kinetics.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in self.kinetics.efficiencies.iteritems():
            for collider0 in kinetics.efficiencies:
                if collider.toSMILES() == collider0.toSMILES():
                    self.assertEqual(efficiency, kinetics.efficiencies[collider0])
                    break
            else:
                self.fail('Unable to find collider "{0}" in reconstructed object.'.format(collider))

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testOutput(self):
        """
        Test that a ThirdBody object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        exec('kinetics = {0!r}'.format(self.kinetics))
        
        self.assertEqual(self.kinetics.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(self.kinetics.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(self.kinetics.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(self.kinetics.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in self.kinetics.efficiencies.iteritems():
            for collider0 in kinetics.efficiencies:
                if collider.toSMILES() == collider0.toSMILES():
                    self.assertEqual(efficiency, kinetics.efficiencies[collider0])
                    break
            else:
                self.fail('Unable to find collider "{0}" in reconstructed object.'.format(collider))
        
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

################################################################################

class TestLindemann(unittest.TestCase):
    """
    Contains unit tests of the Lindemann class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        arrhHigh = Arrhenius(
            A = (1.0e6,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (10,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        arrhLow = Arrhenius(
            A = (1.0e12,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (20,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        efficiencies = {'N#N': 0.5, '[Ar]': 1.5}
        self.kinetics = Lindemann(
            arrheniusLow = arrhLow, 
            arrheniusHigh = arrhHigh, 
            efficiencies = efficiencies, 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            Pmin = (0.1,"bar"), 
            Pmax = (10.0,"bar"), 
            comment = """This data is completely made up""",
        )
        
    def testIsPressureDependent(self):
        """
        Test the Lindemann.isPressureDependent() method.
        
        """
        self.assertTrue(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the Lindemann.getRateCoefficient() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            for P in [1e4,1e5,1e6]:
                k0 = self.kinetics.getRateCoefficient(T, P)
                khigh = self.kinetics.arrheniusHigh.getRateCoefficient(T)
                klow = self.kinetics.arrheniusLow.getRateCoefficient(T)
                C = P / 8.314472 / T
                Pr = klow * C / khigh
                k = khigh * (Pr / (1 + Pr))
                self.assertAlmostEqual(k / k0, 1, 6)
        
    def testPickle(self):
        """
        Test that a Lindemann object can be successfully pickled and
        unpickled with no loss of information.
        """

        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        
        self.assertEqual(self.kinetics.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(self.kinetics.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(self.kinetics.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(self.kinetics.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(self.kinetics.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(self.kinetics.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(self.kinetics.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(self.kinetics.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in self.kinetics.efficiencies.iteritems():
            for collider0 in kinetics.efficiencies:
                if collider.toSMILES() == collider0.toSMILES():
                    self.assertEqual(efficiency, kinetics.efficiencies[collider0])
                    break
            else:
                self.fail('Unable to find collider "{0}" in reconstructed object.'.format(collider))
        
        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testOutput(self):
        """
        Test that a Lindemann object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        exec('kinetics = {0!r}'.format(self.kinetics))
        
        self.assertEqual(self.kinetics.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(self.kinetics.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(self.kinetics.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(self.kinetics.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(self.kinetics.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(self.kinetics.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(self.kinetics.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(self.kinetics.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in self.kinetics.efficiencies.iteritems():
            for collider0 in kinetics.efficiencies:
                if collider.toSMILES() == collider0.toSMILES():
                    self.assertEqual(efficiency, kinetics.efficiencies[collider0])
                    break
            else:
                self.fail('Unable to find collider "{0}" in reconstructed object.'.format(collider))

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

################################################################################

class TestTroe(unittest.TestCase):
    """
    Contains unit tests of the Troe class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        arrhHigh = Arrhenius(
            A = (1.0e6,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (10,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        arrhLow = Arrhenius(
            A = (1.0e12,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (20,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        efficiencies = {'N#N': 0.5, '[Ar]': 1.5}
        self.kinetics = Troe(
            arrheniusLow = arrhLow, 
            arrheniusHigh = arrhHigh, 
            efficiencies = efficiencies, 
            alpha = 0.6, 
            T3 = (1000.0,"K"), 
            T1 = (500.0,"K"), 
            T2 = (300.0,"K"),
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            Pmin = (0.1,"bar"), 
            Pmax = (10.0,"bar"), 
            comment = """This data is completely made up""",
        )
        
    def testIsPressureDependent(self):
        """
        Test the Troe.isPressureDependent() method.
        
        """
        self.assertTrue(self.kinetics.isPressureDependent())
    
    def testRateCoefficient(self):
        """
        Test the Troe.getRateCoefficient() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            for P in [1e4,1e5,1e6]:
                k0 = self.kinetics.getRateCoefficient(T, P)
                khigh = self.kinetics.arrheniusHigh.getRateCoefficient(T)
                klow = self.kinetics.arrheniusLow.getRateCoefficient(T)
                C = P / 8.314472 / T
                Pr = klow * C / khigh
                Fcent = (1 - self.kinetics.alpha.value) * math.exp(-T / self.kinetics.T3.value) + self.kinetics.alpha.value * math.exp(-T / self.kinetics.T1.value)
                if self.kinetics.T2 is not None: Fcent += math.exp(-self.kinetics.T2.value / T)
                d = 0.14
                n = 0.75 - 1.27 * math.log10(Fcent)
                c = -0.4 - 0.67 * math.log10(Fcent)
                F = 10.0**(math.log10(Fcent)/(1 + ((math.log10(Pr) + c)/(n - d * (math.log10(Pr))))**2))
                k = khigh * (Pr / (1 + Pr)) * F
                self.assertAlmostEqual(k / k0, 1, 6)
        
    def testPickle(self):
        """
        Test that a Troe object can be successfully pickled and
        unpickled with no loss of information.
        """

        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))

        self.assertEqual(self.kinetics.alpha.value, kinetics.alpha.value)
        self.assertEqual(self.kinetics.T3.value, kinetics.T3.value)
        self.assertEqual(self.kinetics.T1.value, kinetics.T1.value)
        self.assertEqual(self.kinetics.T2.value, kinetics.T2.value)
        self.assertEqual(self.kinetics.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(self.kinetics.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(self.kinetics.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(self.kinetics.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(self.kinetics.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(self.kinetics.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(self.kinetics.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(self.kinetics.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in self.kinetics.efficiencies.iteritems():
            for collider0 in kinetics.efficiencies:
                if collider.toSMILES() == collider0.toSMILES():
                    self.assertEqual(efficiency, kinetics.efficiencies[collider0])
                    break
            else:
                self.fail('Unable to find collider "{0}" in reconstructed object.'.format(collider))

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

    def testOutput(self):
        """
        Test that a Troe object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        exec('kinetics = {0!r}'.format(self.kinetics))
        
        self.assertEqual(self.kinetics.alpha.value, kinetics.alpha.value)
        self.assertEqual(self.kinetics.T3.value, kinetics.T3.value)
        self.assertEqual(self.kinetics.T1.value, kinetics.T1.value)
        self.assertEqual(self.kinetics.T2.value, kinetics.T2.value)
        self.assertEqual(self.kinetics.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(self.kinetics.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(self.kinetics.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(self.kinetics.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(self.kinetics.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(self.kinetics.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(self.kinetics.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(self.kinetics.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in self.kinetics.efficiencies.iteritems():
            for collider0 in kinetics.efficiencies:
                if collider.toSMILES() == collider0.toSMILES():
                    self.assertEqual(efficiency, kinetics.efficiencies[collider0])
                    break
            else:
                self.fail('Unable to find collider "{0}" in reconstructed object.'.format(collider))

        self.assertEqual(self.kinetics.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(self.kinetics.Tmin.units, kinetics.Tmin.units)
        self.assertEqual(self.kinetics.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(self.kinetics.Tmax.units, kinetics.Tmax.units)
        self.assertEqual(self.kinetics.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(self.kinetics.Pmin.units, kinetics.Pmin.units)
        self.assertEqual(self.kinetics.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(self.kinetics.Pmax.units, kinetics.Pmax.units)
        self.assertEqual(self.kinetics.comment, kinetics.comment)

################################################################################

class TestIdentity(unittest.TestCase):
    """
    Tests the isIdenticalTo function for the different ypes of kinetics:
    KineticsModel, KineticsData, Arrhenius, ArrheniusEP, PDepArrhenius
    Chebyshev, ThirdBody, Lindemann, Troe
    """

    def testIdentity(self):
        arrhenius1 = Arrhenius(A=(1.21e+13,"cm^3/(mol*s)"), n=0, Ea=(0,"kcal/mol"), T0=(1,"K"))
        arrhenius1b = Arrhenius(A=(1.21e+13,"cm^3/(mol*s)"), n=0, Ea=(0,"kcal/mol"), T0=(1,"K"))
        arrhenius2 = Arrhenius(A=(1.21e+13,"cm^3/(mol*s)"), n=0, Ea=(0,"kcal/mol"), T0=(1,"K"), Tmin=(353,"K"),Tmax =(400,"K"))
        arrhenius3= Arrhenius(A=(1.75e+10,"cm^3/(mol*s)"), n=0, Ea=(-3.28,"kcal/mol"), T0=(1,"K"))
        chebyshev1=Chebyshev(coeffs=[[12.71,0.052197,-0.0026878,5.7448e-06], [-0.43104,0.099138,-0.0048957,-1.404e-05], [-0.281,0.084727,-0.0036625,-7.4128e-05], [-0.19107,0.064542,-0.0021801,-0.0001291]], Tmin=(300,"K"), Tmax=(2500,"K"), Pmin=(4.93462,"atm"), Pmax=(59.2154,"atm"), comment="""NetReaction from PDepNetwork #2 (DMP) High-P Limit: R_Recombination exact:  [ C_rad/H/NonDeC , CO_rad/NonDe ] For the above reaction, deltaHrxn(T=298K) = -87.3 kcal/mol""")
        chebyshev1b=Chebyshev(coeffs=[[12.71,0.052197,-0.0026878,5.7448e-06], [-0.43104,0.099138,-0.0048957,-1.404e-05], [-0.281,0.084727,-0.0036625,-7.4128e-05], [-0.19107,0.064542,-0.0021801,-0.0001291]], Tmin=(300,"K"), Tmax=(2500,"K"), Pmin=(4.93462,"atm"), Pmax=(59.2154,"atm"), comment="""NetReaction from PDepNetwork #2 (DMP) High-P Limit: R_Recombination exact:  [ C_rad/H/NonDeC , CO_rad/NonDe ]""")
        chebyshev2=Chebyshev(coeffs=[[7.9506,0.53959,-1.586e-06,-3.1831e-07], [-0.60789,4.1618e-06,1.2168e-06,2.4422e-07], [-0.217,7.4673e-07,2.1833e-07,4.382e-08], [-0.075306,2.4397e-08,7.1345e-09,1.4325e-09]], Tmin=(300,"K"), Tmax=(2500,"K"), Pmin=(4.93462,"atm"), Pmax=(59.2154,"atm"), comment="""NetReaction from PDepNetwork #5 (HO2) High-P Limit: Oa_R_Recombination estimate: (Average:) [ O_pri_rad , Oa ] For the above reaction, deltaHrxn(T=298K) = -66.0 kcal/mol""")
        thirdbody1=ThirdBody(arrheniusHigh=Arrhenius(A=(1.2e+17,"cm^6/(mol^2*s)"),n=-1,Ea=(0,"kcal/mol"),T0=(1,"K"),),efficiencies={
            "C":2,"C(=O)=O":3.6,"CC":3,"O":15.4,"[Ar]":0.83,"[C]=O":1.75,"[H][H]":2.4,},comment="""ReactionLibrary:GRI-Mech3.0""")
        thirdbody1b=ThirdBody(arrheniusHigh=Arrhenius(A=(1.2e+17,"cm^6/(mol^2*s)"),n=-1,Ea=(0,"kcal/mol"),T0=(1,"K"),),efficiencies={
            "C":2,"C(=O)=O":3.6,"CC":3,"O":15.4,"[Ar]":0.83,"[C]=O":1.75,"[H][H]":2.4,},comment="""ReactionLibrary""")
        thirdbody2=ThirdBody(arrheniusHigh=Arrhenius(A=(1.2e+17,"cm^6/(mol^2*s)"),n=-.5,Ea=(0,"kcal/mol"),T0=(1,"K"),),efficiencies={
                "C":2,"C(=O)=O":3.6,"CC":3,"O":15.4,"[Ar]":0.83,"[C]=O":1.75,"[H][H]":2.4,},comment="""ReactionLibrary:GRI-Mech3.0""")
        thirdbody3=ThirdBody(arrheniusHigh=Arrhenius(A=(1.2e+17,"cm^6/(mol^2*s)"),n=-1,Ea=(0,"kcal/mol"),T0=(1,"K"),),efficiencies={
            "C":2,"CC=O":3.6,"CC":3,"O":15.4,"[Ar]":0.83,"[C]=O":1.75,"[H][H]":2.4,},comment="""ReactionLibrary:GRI-Mech3.0""")
        troe1 = Troe(arrheniusHigh=Arrhenius(A=(7.4e+13,"cm^3/(mol*s)"),n=-0.37,Ea=(0,"kcal/mol"),T0=(1,"K"), ), arrheniusLow=Arrhenius(A=(2.3e+18,"cm^6/(mol^2*s)"),n=-0.9,Ea=(-1.7,"kcal/mol"),T0=(1,"K"), ), alpha=0.7346, T3=(94,"K"), T1=(1756,"K"), T2=(5182,"K"), efficiencies={
            "C": 2,"C(=O)=O": 2,"CC": 3,"O": 6,"[Ar]": 0.7,"[C]=O": 1.5,"[H][H]": 2, }, comment="""ReactionLibrary: GRI-Mech3.0   """)
        troe1b = Troe(arrheniusHigh=Arrhenius(A=(7.4e+13,"cm^3/(mol*s)"),n=-0.37,Ea=(0,"kcal/mol"),T0=(1,"K"), ), arrheniusLow=Arrhenius(A=(2.3e+18,"cm^6/(mol^2*s)"),n=-0.9,Ea=(-1.7,"kcal/mol"),T0=(1,"K"), ), alpha=0.7346, T3=(94,"K"), T1=(1756,"K"), T2=(5182,"K"), efficiencies={
            "C": 2,"C(=O)=O": 2,"CC": 3,"O": 6,"[Ar]": 0.7,"[C]=O": 1.5,"[H][H]": 2, }, comment="""Reaction""")
        troe2 = Troe(arrheniusHigh=Arrhenius(A=(7.4e+13,"cm^3/(mol*s)"),n=-0.57,Ea=(0,"kcal/mol"),T0=(1,"K"), ), arrheniusLow=Arrhenius(A=(2.3e+18,"cm^6/(mol^2*s)"),n=-0.9,Ea=(-1.7,"kcal/mol"),T0=(1,"K"), ), alpha=0.7346, T3=(94,"K"), T1=(1756,"K"), T2=(5182,"K"), efficiencies={
            "C": 2,"C(=O)=O": 2,"CC": 3,"O": 6,"[Ar]": 0.7,"[C]=O": 1.5,"[H][H]": 2, }, comment="""Reaction""")
        troe3 = Troe(arrheniusHigh=Arrhenius(A=(7.4e+13,"cm^3/(mol*s)"),n=-0.37,Ea=(0,"kcal/mol"),T0=(1,"K"), ), arrheniusLow=Arrhenius(A=(2.3e+18,"cm^6/(mol^2*s)"),n=-0.9,Ea=(-1.7,"kcal/mol"),T0=(1,"K"), ), alpha=0.7346, T3=(94,"K"), T1=(1756,"K"), T2=(5182,"K"), efficiencies={
            "C": 2,"CC=O": 2,"CC": 3,"O": 6,"[Ar]": 0.7,"[C]=O": 1.5,"[H][H]": 2, }, comment="""Reaction""")
        plog1 = PDepArrhenius(pressures=[3000, 10000],arrhenius=[arrhenius1,arrhenius3],highPlimit=arrhenius1)
        plog1b = PDepArrhenius(pressures=[3000, 10000],arrhenius=[arrhenius1b,arrhenius3],highPlimit=arrhenius1b)
        plog2 = PDepArrhenius(pressures=[3000, 10000],arrhenius=[arrhenius1,arrhenius3],highPlimit=arrhenius1, Tmin=(300,"K"),Tmax=(1500,"K"))
        plog3 = PDepArrhenius(pressures=[3000, 10000, 20000],arrhenius=[arrhenius1,arrhenius3, arrhenius2],highPlimit=arrhenius1)
        multikinetics1 = MultiKinetics(kineticsList=[arrhenius1,troe1])
        multikinetics1b = MultiKinetics(kineticsList=[arrhenius1b,troe1b])
        multikinetics2 = MultiKinetics(kineticsList=[thirdbody1,troe2])  

        self.assertTrue(arrhenius1.isIdenticalTo(arrhenius1b))
        self.assertFalse(arrhenius1.isIdenticalTo(arrhenius2))
        self.assertFalse(arrhenius1.isIdenticalTo(arrhenius3))
        self.assertFalse(arrhenius2.isIdenticalTo(arrhenius3))
        self.assertTrue(chebyshev1.isIdenticalTo(chebyshev1b))
        self.assertFalse(chebyshev1.isIdenticalTo(chebyshev2))
        self.assertTrue(thirdbody1.isIdenticalTo(thirdbody1b))
        self.assertFalse(thirdbody1.isIdenticalTo(thirdbody2))
        self.assertFalse(thirdbody1.isIdenticalTo(thirdbody3))
        self.assertFalse(thirdbody2.isIdenticalTo(thirdbody3))
        self.assertTrue(troe1.isIdenticalTo(troe1b))
        self.assertFalse(troe1.isIdenticalTo(troe2))
        self.assertFalse(troe1.isIdenticalTo(troe3))
        self.assertFalse(troe2.isIdenticalTo(troe3))
        self.assertTrue(plog1.isIdenticalTo(plog1b))
        self.assertFalse(plog1.isIdenticalTo(plog2))
        self.assertFalse(plog1.isIdenticalTo(plog3))
        self.assertTrue(multikinetics1.isIdenticalTo(multikinetics1b))
        self.assertFalse(multikinetics1.isIdenticalTo(multikinetics2))

        self.assertFalse(arrhenius1.isIdenticalTo(chebyshev1))
        self.assertFalse(arrhenius1.isIdenticalTo(thirdbody1))
        self.assertFalse(arrhenius1.isIdenticalTo(troe1))
        self.assertFalse(arrhenius1.isIdenticalTo(plog1))
        self.assertFalse(arrhenius1.isIdenticalTo(multikinetics1))
        self.assertFalse(chebyshev1.isIdenticalTo(thirdbody1))
        self.assertFalse(chebyshev1.isIdenticalTo(troe1))
        self.assertFalse(thirdbody1.isIdenticalTo(troe1))
        self.assertFalse(troe1.isIdenticalTo(thirdbody1))

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
