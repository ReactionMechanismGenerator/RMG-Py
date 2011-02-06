#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

import math
import numpy
import logging
import pylab

from rmgpy.chem.molecule import Molecule
from rmgpy.chem.states import *
from rmgpy.chem.thermo import ThermoGAModel
import rmgpy.chem.constants as constants

from rmgpy.data.statesfit import fitStatesToHeatCapacity

showPlots = False

################################################################################

def makeHinderedRotor(freq, barr, symm=1):
    symmetry = symm
    barrier = barr * constants.h * constants.c * 100. * constants.Na
    frequency = freq * constants.c * 100.
    inertia = symmetry * symmetry * (barrier / constants.Na) / (8 * math.pi * math.pi * frequency * frequency)
    return HinderedRotor(inertia=inertia, barrier=barrier, symmetry=symmetry)

################################################################################

class StatesfitCheck(unittest.TestCase):
    """
    Unit tests of the rmgpy.statesfit package. These tests are designed to see
    how well a known set of vibrations and/or hindered rotors are recreated by
    the fitting algorithm.
    """

    def testOneVibration(self):
        """
        Unit tests involving fitting of one vibrational frequency and zero
        hindered rotors.
        """
        frequencies = [500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
        Tlist = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)

        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        for i, freq in enumerate(frequencies):
            ho = HarmonicOscillator(frequencies=[freq])
            Cv_data = ho.getHeatCapacities(Tlist) / constants.R

            modes = fitStatesToHeatCapacity(Tlist, Cv_data, Nvib=1, Nrot=0)
            
            self.assertTrue(len(modes) == 1)
            self.assertTrue(isinstance(modes[0], HarmonicOscillator))
            self.assertTrue(len(modes[0].frequencies) == 1)
            self.assertAlmostEqual(modes[0].frequencies[0] / freq, 1.0, 4, 'Fitted vibrational frequency of %g cm^-1 does not match expected value of %g cm^-1.' % (modes[0].frequencies[0], freq))

            Cv_model = numpy.zeros_like(Tlist)
            for mode in modes:
                Cv_model += mode.getHeatCapacities(Tlist) / constants.R

            pylab.plot(Tlist, Cv_data, 'o%s' % colors[i], Tlist, Cv_model, '-%s' % colors[i])

        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Heat capacity / R')
        if showPlots: pylab.show()

    def testOneRotor(self):
        """
        Unit tests involving fitting of zero vibrational frequencies and one
        hindered rotor.
        """
        rotors = [(300.0, 100.0), (300.0, 1000.0), (300.0, 3000.0), (1000.0, 100.0), (1000.0, 500.0), (1000.0, 3000.0)]
        Tlist = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)

        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        for i, rotor in enumerate(rotors):
            freq, barr = rotor

            hr = makeHinderedRotor(freq, barr)
            Cv_data = hr.getHeatCapacities(Tlist) / constants.R

            modes = fitStatesToHeatCapacity(Tlist, Cv_data, Nvib=0, Nrot=1)
            
            self.assertTrue(len(modes) == 1)
            self.assertTrue(isinstance(modes[0], HinderedRotor))
            freq0 = modes[0].getFrequency()
            barr0 = modes[0].barrier / (constants.h * constants.c * 100.0 * constants.Na)
            self.assertAlmostEqual(freq0 / freq, 1.0, 4, 'Fitted rotor frequency of %g cm^-1 does not match expected value of %g cm^-1.' % (freq0, freq))
            self.assertAlmostEqual(barr0 / barr, 1.0, 4, 'Fitted rotor barrier of %g cm^-1 does not match expected value of %g cm^-1.' % (barr0, barr))

            Cv_model = numpy.zeros_like(Tlist)
            for mode in modes:
                Cv_model += mode.getHeatCapacities(Tlist) / constants.R

            pylab.plot(Tlist, Cv_data, 'o%s' % colors[i], Tlist, Cv_model, '-%s' % colors[i])

        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Heat capacity / R')
        if showPlots: pylab.show()

    def testTwoVibrations(self):
        """
        Unit tests involving fitting of two vibrational frequencies and zero
        hindered rotors.
        """
        frequencies = [[450.0, 550.0], [450.0, 1550.0], [450.0, 2550.0], [1450.0, 1550.0], [1450.0, 2550.0], [2450.0, 2550.0]]
        Tlist = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)

        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        for i, freq in enumerate(frequencies):
            ho = HarmonicOscillator(frequencies=freq)
            Cv_data = ho.getHeatCapacities(Tlist) / constants.R

            modes = fitStatesToHeatCapacity(Tlist, Cv_data, Nvib=2, Nrot=0)
            
            self.assertTrue(len(modes) == 1)
            self.assertTrue(isinstance(modes[0], HarmonicOscillator))
            self.assertTrue(len(modes[0].frequencies) == 2)
            self.assertAlmostEqual(modes[0].frequencies[0] / freq[0], 1.0, 4, 'Fitted vibrational frequency of %g cm^-1 does not match expected value of %g cm^-1.' % (modes[0].frequencies[0], freq[0]))
            self.assertAlmostEqual(modes[0].frequencies[1] / freq[1], 1.0, 4, 'Fitted vibrational frequency of %g cm^-1 does not match expected value of %g cm^-1.' % (modes[0].frequencies[1], freq[1]))

            Cv_model = numpy.zeros_like(Tlist)
            for mode in modes:
                Cv_model += mode.getHeatCapacities(Tlist) / constants.R

            pylab.plot(Tlist, Cv_data, 'o%s' % colors[i], Tlist, Cv_model, '-%s' % colors[i])

        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Heat capacity / R')
        if showPlots: pylab.show()

    def testOneFrequencyAndOneRotor(self):
        """
        Unit tests involving fitting of one vibrational frequency and zero
        hindered rotors.
        """
        frequencies = [1500.0, 1500.0, 1500.0, 1500.0]
        rotors = [(300.0, 100.0), (300.0, 1000.0), (1000.0, 100.0), (1000.0, 500.0)]
        Tlist = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)

        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        for i in range(len(frequencies)):
            vibFreq = frequencies[i]
            freq, barr = rotors[i]

            ho = HarmonicOscillator(frequencies=[vibFreq])
            hr = makeHinderedRotor(freq, barr)
            Cv_data = ho.getHeatCapacities(Tlist) / constants.R + hr.getHeatCapacities(Tlist) / constants.R

            modes = fitStatesToHeatCapacity(Tlist, Cv_data, Nvib=1, Nrot=1)
            
            self.assertTrue(len(modes) == 2)
            self.assertTrue(isinstance(modes[0], HarmonicOscillator))
            self.assertTrue(len(modes[0].frequencies) == 1)
            self.assertTrue(isinstance(modes[1], HinderedRotor))
            freq0 = modes[1].getFrequency()
            barr0 = modes[1].barrier / (constants.h * constants.c * 100.0 * constants.Na)
            self.assertAlmostEqual(modes[0].frequencies[0] / vibFreq, 1.0, 3, 'Fitted vibrational frequency of %g cm^-1 does not match expected value of %g cm^-1.' % (modes[0].frequencies[0], vibFreq))
            self.assertAlmostEqual(freq0 / freq, 1.0, 4, 'Fitted rotor frequency of %g cm^-1 does not match expected value of %g cm^-1.' % (freq0, freq))
            self.assertAlmostEqual(barr0 / barr, 1.0, 4, 'Fitted rotor barrier of %g cm^-1 does not match expected value of %g cm^-1.' % (barr0, barr))

            Cv_model = numpy.zeros_like(Tlist)
            for mode in modes:
                Cv_model += mode.getHeatCapacities(Tlist) / constants.R

            pylab.plot(Tlist, Cv_data, 'o%s' % colors[i], Tlist, Cv_model, '-%s' % colors[i])

        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Heat capacity / R')
        if showPlots: pylab.show()

    def testManyVibrations(self):
        """
        Unit tests involving fitting of many vibrational frequency and zero
        hindered rotors.
        """
        frequencies = numpy.arange(300.0, 3001.0, 300.0, numpy.float64)
        Tlist = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)

        ho = HarmonicOscillator(frequencies=list(frequencies))
        Cv_data = ho.getHeatCapacities(Tlist) / constants.R

        modes = fitStatesToHeatCapacity(Tlist, Cv_data, Nvib=len(frequencies), Nrot=0)
        
        self.assertTrue(len(modes) == 1)
        self.assertTrue(isinstance(modes[0], HarmonicOscillator))
        self.assertTrue(len(modes[0].frequencies) == len(frequencies))
        for i in range(len(frequencies)):
            self.assertAlmostEqual(modes[0].frequencies[i] / frequencies[i], 1.0, 2, 'Fitted vibrational frequency of %g cm^-1 does not match expected value of %g cm^-1.' % (modes[0].frequencies[i], frequencies[i]))
        
        T_model = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)
        Cv_model = numpy.zeros_like(T_model)
        for mode in modes:
            Cv_model += mode.getHeatCapacities(T_model) / constants.R

        pylab.plot(Tlist, Cv_data, 'ok', T_model, Cv_model, '-k')

        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Heat capacity / R')
        if showPlots: pylab.show()

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
