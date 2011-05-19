#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.statmech module.
"""

import math
import numpy
import unittest

from rmgpy.statmech import *
from rmgpy.quantity import constants

################################################################################

class TestTranslation(unittest.TestCase):
    """
    Contains unit tests of the Translation class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.mode = Translation(mass=(28.03,"u"))

    def testPartitionFunction(self):
        """
        Test the Translation.getPartitionFunction() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Q = self.mode.getPartitionFunction(T)
            
    def testHeatCapacity(self):
        """
        Test the Translation.getHeatCapacity() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(self.mode.getHeatCapacity(T), 1.5 * constants.R, 6)

    def testEnthalpy(self):
        """
        Test the Translation.getEnthalpy() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(self.mode.getEnthalpy(T), 1.5 * constants.R * T, 6)

    def testEntropy(self):
        """
        Test the Translation.getEntropy() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Q = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(self.mode.getEntropy(T), (math.log(Q) + 2.5) * constants.R, 6)

    def testDensityOfStates(self):
        """
        Test the Translation.getDensityOfStates() method.
        """
        Elist = numpy.arange(0, 100000, 2000, numpy.float)
        densStates = self.mode.getDensityOfStates(Elist)
    
    def testOutput(self):
        """
        Test that a Translation object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = %r' % self.mode)
        self.assertAlmostEqual(self.mode.mass.value, mode.mass.value, 6)
        self.assertEqual(self.mode.mass.units, mode.mass.units)
        
    def testPickle(self):
        """
        Test that a Translation object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertAlmostEqual(self.mode.mass.value, mode.mass.value, 6)
        self.assertEqual(self.mode.mass.units, mode.mass.units)
        
################################################################################

class TestRigidRotor(unittest.TestCase):
    """
    Contains unit tests of the RigidRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.mode = RigidRotor(linear=False, inertia=([1.0,10.0,100.0],"amu*angstrom^2"), symmetry=3)

    def testPartitionFunction(self):
        """
        Test the RigidRotor.getPartitionFunction() method.
        """
        pass
            
    def testHeatCapacity(self):
        """
        Test the RigidRotor.getHeatCapacity() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(self.mode.getHeatCapacity(T), 1.5 * constants.R, 6)

    def testEnthalpy(self):
        """
        Test the RigidRotor.getEnthalpy() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            self.assertAlmostEqual(self.mode.getEnthalpy(T), 1.5 * constants.R * T, 6)

    def testEntropy(self):
        """
        Test the RigidRotor.getEntropy() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Q = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(self.mode.getEntropy(T), (math.log(Q) + 1.5) * constants.R, 6)

    def testDensityOfStates(self):
        """
        Test the RigidRotor.getDensityOfStates() method.
        """
        Elist = numpy.arange(0, 100000, 2000, numpy.float)
        densStates = self.mode.getDensityOfStates(Elist)
    
    def testOutput(self):
        """
        Test that a RigidRotor object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = %r' % self.mode)
        self.assertEqual(self.mode.linear, mode.linear)
        self.assertEqual(len(self.mode.inertia.values), len(mode.inertia.values))
        for I1, I2 in zip(self.mode.inertia.values, mode.inertia.values):
            self.assertAlmostEqual(I1 * 6.022e46, I2 * 6.022e46, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        
    def testPickle(self):
        """
        Test that a RigidRotor object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertEqual(self.mode.linear, mode.linear)
        self.assertEqual(len(self.mode.inertia.values), len(mode.inertia.values))
        for I1, I2 in zip(self.mode.inertia.values, mode.inertia.values):
            self.assertAlmostEqual(I1 * 6.022e46, I2 * 6.022e46, 6)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        
################################################################################

class TestHarmonicOscillator(unittest.TestCase):
    """
    Contains unit tests of the HarmonicOscillator class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.mode = HarmonicOscillator(frequencies=[834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3, 3121.6, 3136.7, 3192.5, 3221.0])
        
    def testPartitionFunction(self):
        """
        Test the HarmonicOscillator.getPartitionFunction() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Q = self.mode.getPartitionFunction(T)
            
    def testHeatCapacity(self):
        """
        Test the HarmonicOscillator.getHeatCapacity() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Cp = self.mode.getHeatCapacity(T)

    def testEnthalpy(self):
        """
        Test the HarmonicOscillator.getEnthalpy() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            H = self.mode.getEnthalpy(T)

    def testEntropy(self):
        """
        Test the HarmonicOscillator.getEntropy() method.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            S = self.mode.getEntropy(T)

    def testDensityOfStates(self):
        """
        Test the HarmonicOscillator.getDensityOfStates() method.
        """
        Elist = numpy.arange(0, 100000, 2000, numpy.float)
        densStates = self.mode.getDensityOfStates(Elist)
    
    def testOutput(self):
        """
        Test that a HarmonicOscillator object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = %r' % self.mode)
        self.assertEqual(len(self.mode.frequencies.values), len(mode.frequencies.values))
        for freq1, freq2 in zip(self.mode.frequencies.values, mode.frequencies.values):
            self.assertAlmostEqual(freq1, freq2, 4)
        self.assertEqual(self.mode.frequencies.units, mode.frequencies.units)
        
    def testPickle(self):
        """
        Test that a HarmonicOscillator object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertEqual(len(self.mode.frequencies.values), len(mode.frequencies.values))
        for freq1, freq2 in zip(self.mode.frequencies.values, mode.frequencies.values):
            self.assertAlmostEqual(freq1, freq2, 4)
        self.assertEqual(self.mode.frequencies.units, mode.frequencies.units)
        
################################################################################

class TestHinderedRotor(unittest.TestCase):
    """
    Contains unit tests of the HinderedRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        fourier = numpy.array([ [-4.683e-01, 8.767e-05], [-2.827e+00, 1.048e-03], [ 1.751e-01,-9.278e-05], [-1.355e-02, 1.916e-06], [-1.128e-01, 1.025e-04] ], numpy.float64) * 4184
        self.mode = HinderedRotor(inertia=(7.38359,"amu*angstrom^2"), barrier=(3.20429,"kcal/mol"), symmetry=1, fourier=fourier)

    def testPartitionFunction(self):
        """
        Test the HinderedRotor.getPartitionFunction() method, both with a 
        Fourier series potential and a cosine potential.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Q = self.mode.getPartitionFunction(T)
        self.mode.fourier = None; self.mode.energies = None
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Q = self.mode.getPartitionFunction(T)
            
    def testHeatCapacity(self):
        """
        Test the HinderedRotor.getHeatCapacity() method, both with a Fourier 
        series potential and a cosine potential.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Cp = self.mode.getHeatCapacity(T)
        self.mode.fourier = None; self.mode.energies = None
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            Cp = self.mode.getHeatCapacity(T)
            
    def testEnthalpy(self):
        """
        Test the HinderedRotor.getEnthalpy() method, both with a Fourier series
        potential and a cosine potential.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            H = self.mode.getEnthalpy(T)
        self.mode.fourier = None; self.mode.energies = None
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            H = self.mode.getEnthalpy(T)
        
    def testEntropy(self):
        """
        Test the HinderedRotor.getEntropy() method, both with a Fourier series
        potential and a cosine potential.
        """
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            S = self.mode.getEntropy(T)
        self.mode.fourier = None; self.mode.energies = None
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            S = self.mode.getEntropy(T)
        
    def testDensityOfStates(self):
        """
        Test the HinderedRotor.getDensityOfStates() method.
        """
        Elist = numpy.arange(0, 100000, 2000, numpy.float)
        densStates = self.mode.getDensityOfStates(Elist)
        self.mode.fourier = None; self.mode.energies = None
        densStates = self.mode.getDensityOfStates(Elist)
        
    def testFrequency(self):
        """
        Test the HinderedRotor.getFrequency() method.
        """
        freq = self.mode.getFrequency()
        self.assertAlmostEqual(freq, 50.9164, 3)
    
    def testOutput(self):
        """
        Test that a HinderedRotor object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = %r' % self.mode)
        self.assertAlmostEqual(self.mode.barrier.value / 1000., mode.barrier.value / 1000., 4)
        self.assertEqual(self.mode.barrier.units, mode.barrier.units)
        self.assertAlmostEqual(self.mode.inertia.value * 6.022e46, mode.inertia.value * 6.022e46, 4)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.fourier.values.shape[0], mode.fourier.values.shape[0])
        self.assertEqual(self.mode.fourier.values.shape[1], mode.fourier.values.shape[1])
        for i in range(self.mode.fourier.values.shape[0]):
            for j in range(self.mode.fourier.values.shape[1]):
                self.assertAlmostEqual(self.mode.fourier.values[i,j] / 1000., mode.fourier.values[i,j] / 1000., 4)
        
    def testPickle(self):
        """
        Test that a HinderedRotor object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertAlmostEqual(self.mode.barrier.value / 1000., mode.barrier.value / 1000., 4)
        self.assertEqual(self.mode.barrier.units, mode.barrier.units)
        self.assertAlmostEqual(self.mode.inertia.value * 6.022e46, mode.inertia.value * 6.022e46, 4)
        self.assertEqual(self.mode.inertia.units, mode.inertia.units)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.fourier.values.shape[0], mode.fourier.values.shape[0])
        self.assertEqual(self.mode.fourier.values.shape[1], mode.fourier.values.shape[1])
        for i in range(self.mode.fourier.values.shape[0]):
            for j in range(self.mode.fourier.values.shape[1]):
                self.assertAlmostEqual(self.mode.fourier.values[i,j] / 1000., mode.fourier.values[i,j] / 1000., 4)
    
    def testModerateBarrierRotor(self):
        """
        Compare the Fourier series and cosine potentials for a hindered rotor
        with a moderate barrier.
        """
        fourier = numpy.array([ [-4.683e-01,-2.827e+00, 1.751e-01,-1.355e-02,-1.128e-01], [ 8.767e-05, 1.048e-03,-9.278e-05, 1.916e-06, 1.025e-04] ], numpy.float64) * 4184
        hr1 = HinderedRotor(inertia=(7.38359,"amu*angstrom^2"), barrier=(2139.3*11.96,"J/mol"), symmetry=2)
        hr2 = HinderedRotor(inertia=(7.38359,"amu*angstrom^2"), barrier=(3.20429*4184,"J/mol"), symmetry=1, fourier=fourier)
        ho = HarmonicOscillator(frequencies=[hr1.getFrequency()])

        # Check that it matches the harmonic oscillator model at low T
        Tlist = numpy.arange(10, 41.0, 1.0, numpy.float64)
        Q1 = hr1.getPartitionFunctions(Tlist)
        Q2 = hr2.getPartitionFunctions(Tlist)
        Q0 = ho.getPartitionFunctions(Tlist)
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Q1[i] / Q0[i], 1.0, 2)
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Q2[i] / Q0[i], 1.0, 2)

    def testLowBarrierRotor(self):
        """
        Compare the Fourier series and cosine potentials for a hindered rotor
        with a low barrier.
        """

        fourier = numpy.array([ [ 1.377e-02,-3.481e-03,-2.511e-01, 6.786e-04,-1.191e-02], [-2.226e-05, 1.859e-05, 2.025e-04,-3.212e-05, 2.027e-05] ], numpy.float64) * 4184
        hr1 = HinderedRotor(inertia=(1.60779,"amu*angstrom^2"), barrier=(176.4*11.96,"J/mol"), symmetry=3)
        hr2 = HinderedRotor(inertia=(1.60779,"amu*angstrom^2"), barrier=(0.233317*4184,"J/mol"), symmetry=3, fourier=fourier)
        
        # Check that the potentials between the two rotors are approximately consistent
        phi = numpy.arange(0, 2*math.pi, math.pi/48.0, numpy.float64)
        V1 = hr1.getPotential(phi)
        V2 = hr2.getPotential(phi)
        
        Vmax = hr1.barrier.value
        for i in range(len(phi)):
            self.assertTrue(abs(V2[i] - V1[i]) / Vmax < 0.1)

        # Check that it matches the harmonic oscillator model at low T
        Tlist = numpy.arange(100.0, 2001.0, 10.0, numpy.float64)
        Q1 = hr1.getPartitionFunctions(Tlist)
        Q2 = hr2.getPartitionFunctions(Tlist)
        C1 = hr1.getHeatCapacities(Tlist)
        C2 = hr2.getHeatCapacities(Tlist)
        H1 = hr1.getEnthalpies(Tlist)
        H2 = hr2.getEnthalpies(Tlist)
        S1 = hr1.getEntropies(Tlist)
        S2 = hr2.getEntropies(Tlist)
        for i in range(len(Tlist)):
            self.assertTrue(abs(Q2[i] - Q1[i]) < 0.1 * Q1[i])
            self.assertTrue(abs(C2[i] - C1[i]) < 0.1 * C1[i])
            self.assertTrue(abs(H2[i] - H1[i]) < 0.1 * H1[i])
            self.assertTrue(abs(S2[i] - S1[i]) < 0.1 * S1[i])

        #import pylab
        #pylab.plot(Tlist, Q1, '-r', Tlist, Q2, '-b')
        #pylab.plot(Tlist, C1, '-r', Tlist, C2, '-b')
        #pylab.plot(Tlist, H1, '-r', Tlist, H2, '-b')
        #pylab.plot(Tlist, S1, '-r', Tlist, S2, '-b')
        #pylab.show()

    
################################################################################

class TestStatmechModel(unittest.TestCase):
    """
    Contains unit tests of the StatmechModel class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.mode = Translation(mass=(28.03,"u"))

    def testOutput(self):
        """
        Test that a Translation object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = %r' % self.mode)
        self.assertAlmostEqual(self.mode.mass.value, mode.mass.value, 6)
        self.assertEqual(self.mode.mass.units, mode.mass.units)
        
    def testPickle(self):
        """
        Test that a Translation object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertAlmostEqual(self.mode.mass.value, mode.mass.value, 6)
        self.assertEqual(self.mode.mass.units, mode.mass.units)
        
################################################################################

class StatesModelTest(unittest.TestCase):
    """
    Contains unit tests for the StatesModel class.
    """
    
    def testModesForEthylene(self):
        """
        Uses data for ethylene (C2H4) to test the various modes. The data comes
        from a CBS-QB3 calculation using Gaussian03.
        """

        T = 298.15

        trans = Translation(mass=0.02803)
        rot = RigidRotor(linear=False, inertia=[5.6952e-47, 2.7758e-46, 3.3454e-46], symmetry=1)
        vib = HarmonicOscillator(frequencies=[834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3, 3121.6, 3136.7, 3192.5, 3221.0])

        self.assertAlmostEqual(trans.getPartitionFunction(T) / 5.83338e6, 1.0, 3)
        self.assertAlmostEqual(rot.getPartitionFunction(T) / 2.59622e3, 1.0, 3)
        self.assertAlmostEqual(vib.getPartitionFunction(T) / 1.0481e0, 1.0, 3)

        self.assertAlmostEqual(trans.getHeatCapacity(T) / 4.184 / 2.981, 1.0, 3)
        self.assertAlmostEqual(rot.getHeatCapacity(T) / 4.184 / 2.981, 1.0, 3)
        self.assertAlmostEqual(vib.getHeatCapacity(T) / 4.184 / 2.133, 1.0, 3)

        self.assertAlmostEqual(trans.getEnthalpy(T) / 8.314472 / T / 1.5, 1.0, 3)
        self.assertAlmostEqual(rot.getEnthalpy(T) / 8.314472 / T / 1.5, 1.0, 3)
        self.assertAlmostEqual(vib.getEnthalpy(T) / 8.314472 / T / 0.221258, 1.0, 3)

        self.assertAlmostEqual(trans.getEntropy(T) / 4.184 / 35.927, 1.0, 2)
        self.assertAlmostEqual(rot.getEntropy(T) / 4.184 / 18.604, 1.0, 3)
        self.assertAlmostEqual(vib.getEntropy(T) / 4.184 / 0.533, 1.0, 3)

        states = StatesModel(modes=[rot, vib], spinMultiplicity=1)

        dE = 10.0
        Elist = numpy.arange(0, 100001, dE, numpy.float64)
        rho = states.getDensityOfStates(Elist)
        self.assertAlmostEqual(numpy.sum(rho * numpy.exp(-Elist / 8.314472 / 298.15) * dE) / states.getPartitionFunction(T), 1.0, 2)

    def testModesForOxygen(self):
        """
        Uses data for oxygen (O2) to test the various modes. The data comes
        from a CBS-QB3 calculation using Gaussian03.
        """

        T = 298.15

        trans = Translation(mass=0.03199)
        rot = RigidRotor(linear=True, inertia=[1.9271e-46], symmetry=2)
        vib = HarmonicOscillator(frequencies=[1637.9])

        self.assertAlmostEqual(trans.getPartitionFunction(T) / 7.11169e6, 1.0, 3)
        self.assertAlmostEqual(rot.getPartitionFunction(T) / 7.13316e1, 1.0, 3)
        self.assertAlmostEqual(vib.getPartitionFunction(T) / 1.000037e0, 1.0, 3)

        self.assertAlmostEqual(trans.getHeatCapacity(T) / 4.184 / 2.981, 1.0, 3)
        self.assertAlmostEqual(rot.getHeatCapacity(T) / 4.184 / 1.987, 1.0, 3)
        self.assertAlmostEqual(vib.getHeatCapacity(T) / 4.184 / 0.046, 1.0, 2)

        self.assertAlmostEqual(trans.getEnthalpy(T) / 8.314472 / T / 1.5, 1.0, 3)
        self.assertAlmostEqual(rot.getEnthalpy(T) / 8.314472 / T / 1.0, 1.0, 3)
        self.assertAlmostEqual(vib.getEnthalpy(T) / 8.314472 / T / 0.0029199, 1.0, 3)

        self.assertAlmostEqual(trans.getEntropy(T) / 4.184 / 36.321, 1.0, 2)
        self.assertAlmostEqual(rot.getEntropy(T) / 4.184 / 10.467, 1.0, 3)
        self.assertAlmostEqual(vib.getEntropy(T) / 4.184 / 0.00654, 1.0, 2)

        states = StatesModel(modes=[rot, vib], spinMultiplicity=3)

        dE = 10.0
        Elist = numpy.arange(0, 100001, dE, numpy.float64)
        rho = states.getDensityOfStates(Elist)
        self.assertAlmostEqual(numpy.sum(rho * numpy.exp(-Elist / 8.314472 / 298.15) * dE) / states.getPartitionFunction(T), 1.0, 2)

    def testDensityOfStatesILTTranslation(self):
        """
        Test that the density of states as obtained via inverse Laplace
        transform of the partition function is equivalent to that obtained
        directly (via convolution) for translation.
        """
        trans = Translation(mass=0.02803)
        Elist = numpy.arange(0.0, 200000.0, 500.0, numpy.float64)
        states = StatesModel(modes=[trans])
        densStates0 = states.getDensityOfStates(Elist)
        densStates1 = states.getDensityOfStatesILT(Elist)
        for i in range(10, len(Elist)):
            self.assertTrue(0.8 < densStates1[i] / densStates0[i] < 1.25)
    
    def testDensityOfStatesILTRotation(self):
        """
        Test that the density of states as obtained via inverse Laplace
        transform of the partition function is equivalent to that obtained
        directly (via convolution) for external rotation.
        """
        rot = RigidRotor(linear=False, inertia=[5.6952e-47, 2.7758e-46, 3.3454e-46], symmetry=1)
        Elist = numpy.arange(0.0, 200000.0, 500.0, numpy.float64)
        states = StatesModel(modes=[rot])
        densStates0 = states.getDensityOfStates(Elist)
        densStates1 = states.getDensityOfStatesILT(Elist)
        for i in range(10, len(Elist)):
            self.assertTrue(0.8 < densStates1[i] / densStates0[i] < 1.25)
    
    def testDensityOfStatesILTVibration(self):
        """
        Test that the density of states as obtained via inverse Laplace
        transform of the partition function is equivalent to that obtained
        directly (via convolution) for vibration. The rotational modes are
        also included to smooth the density of states as obtained directly
        via convolution.
        """
        rot = RigidRotor(linear=False, inertia=[5.6952e-47, 2.7758e-46, 3.3454e-46], symmetry=1)
        vib = HarmonicOscillator(frequencies=[834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3, 3121.6, 3136.7, 3192.5, 3221.0])
        Elist = numpy.arange(0.0, 200000.0, 500.0, numpy.float64)
        states = StatesModel(modes=[rot, vib])
        densStates0 = states.getDensityOfStates(Elist)
        densStates1 = states.getDensityOfStatesILT(Elist)
        for i in range(25, len(Elist)):
            self.assertTrue(0.8 < densStates1[i] / densStates0[i] < 1.25)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
