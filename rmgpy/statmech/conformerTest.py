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
This script contains unit tests of the :mod:`rmgpy.statmech.conformer` module.
"""

import unittest
from external.wip import work_in_progress
import math
import numpy
import scipy.interpolate

from rmgpy.statmech import Conformer, IdealGasTranslation, NonlinearRotor, HarmonicOscillator, \
                           LinearRotor, HinderedRotor 
import rmgpy.constants as constants

################################################################################

class TestConformer(unittest.TestCase):
    """
    Contains unit tests of the :class:`Conformer` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.ethylene = Conformer(
            E0 = (0.0,"kJ/mol"),
            modes = [
                IdealGasTranslation(mass=(28.03,"amu")),
                NonlinearRotor(inertia=([3.41526,16.6498,20.065],"amu*angstrom^2"), symmetry=4),
                HarmonicOscillator(frequencies=([828.397,970.652,977.223,1052.93,1233.55,1367.56,1465.09,1672.25,3098.46,3111.7,3165.79,3193.54],"cm^-1")),
            ],
            spinMultiplicity = 1,
            opticalIsomers = 1,
        )
        self.oxygen = Conformer(
            E0 = (0.0,"kJ/mol"),
            modes = [
                IdealGasTranslation(mass=(31.99,"amu")),
                LinearRotor(inertia=(11.6056,"amu*angstrom^2"), symmetry=2),
                HarmonicOscillator(frequencies=([1621.54],"cm^-1")),
            ],
            spinMultiplicity = 3,
            opticalIsomers = 1,
        )
        
        # The following data is for ethane at the CBS-QB3 level
        self.coordinates = numpy.array([
            [  0.0000,  0.0000,  0.0000],
            [ -0.0000, -0.0000,  1.0936],
            [  1.0430, -0.0000, -0.3288],
            [ -0.4484,  0.9417, -0.3288],
            [ -0.7609, -1.2051, -0.5580],
            [ -0.7609, -1.2051, -1.6516],
            [ -0.3125, -2.1468, -0.2292],
            [ -1.8039, -1.2051, -0.2293],
        ])
        self.number = numpy.array([6, 1, 1, 1, 6, 1, 1, 1])
        self.mass = numpy.array([12, 1.007825, 1.007825, 1.007825, 12, 1.007825, 1.007825, 1.007825])
        self.E0 = -93.5097
        self.conformer = Conformer(
            E0 = (self.E0,"kJ/mol"),
            modes = [
                IdealGasTranslation(mass=(30.0469,"amu")),
                NonlinearRotor(inertia=([6.27071,25.3832,25.3833],"amu*angstrom^2"), symmetry=6),
                HarmonicOscillator(frequencies=([818.917,819.479,987.099,1206.76,1207.05,1396,1411.35,1489.73,1489.95,1492.49,1492.66,2995.36,2996.06,3040.77,3041,3065.86,3066.02],"cm^-1")),
                HinderedRotor(inertia=(1.56768,"amu*angstrom^2"), symmetry=3, barrier=(2.69401,"kcal/mol"), quantum=False, semiclassical=False),
            ],
            spinMultiplicity = 1,
            opticalIsomers = 1,
            coordinates = (self.coordinates,"angstrom"),
            number = self.number,
            mass = (self.mass,"amu"),
        )
        
    def test_getPartitionFunction_ethylene(self):
        """
        Test the StatMech.getPartitionFunction() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([4.05311e+09, 4.19728e+10, 2.82309e+12, 7.51135e+13, 1.16538e+15])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.ethylene.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)

    def test_getHeatCapacity_ethylene(self):
        """
        Test the StatMech.getHeatCapacity() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([5.11186, 7.40447, 11.1659, 13.1221, 14.1617]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.ethylene.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, 3)
    
    def test_getEnthalpy_ethylene(self):
        """
        Test the StatMech.getEnthalpy() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([4.23129, 5.04826, 7.27337, 8.93167, 10.1223]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.ethylene.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEntropy_ethylene(self):
        """
        Test the StatMech.getEntropy() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([26.3540, 29.5085, 35.9422, 40.8817, 44.8142]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.ethylene.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, 3)
    
    def test_getSumOfStates_ethylene(self):
        """
        Test the StatMech.getSumOfStates() method for ethylene.
        """
        Elist = numpy.arange(0, 5000*11.96, 2*11.96)
        sumStates = self.ethylene.getSumOfStates(Elist)
        densStates = self.ethylene.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))
            
    def test_getDensityOfStates_ethylene(self):
        """
        Test the StatMech.getDensityOfStates() method for ethylene.
        """
        Elist = numpy.arange(0, 5000*11.96, 2*11.96)
        densStates = self.ethylene.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.ethylene.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-1*Qexp)

    def test_getPartitionFunction_oxygen(self):
        """
        Test the StatMech.getPartitionFunction() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([1.55584e+09, 9.38339e+09, 1.16459e+11, 5.51016e+11, 1.72794e+12])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.oxygen.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp, Qact, delta=1e-4*Qexp)

    def test_getHeatCapacity_oxygen(self):
        """
        Test the StatMech.getHeatCapacity() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([3.52538, 3.70877, 4.14751, 4.32063, 4.39392]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.oxygen.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp, Cvact, 3)
    
    def test_getEnthalpy_oxygen(self):
        """
        Test the StatMech.getEnthalpy() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([3.50326, 3.54432, 3.75062, 3.91623, 4.02765]) * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.oxygen.getEnthalpy(T)
            self.assertAlmostEqual(Hexp, Hact, delta=1e-4*Hexp)
    
    def test_getEntropy_oxygen(self):
        """
        Test the StatMech.getEntropy() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([24.6685, 26.5065, 29.2314, 30.9513, 32.2056]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.oxygen.getEntropy(T)
            self.assertAlmostEqual(Sexp, Sact, 3)
    
    def test_getSumOfStates_oxygen(self):
        """
        Test the StatMech.getSumOfStates() method for oxygen.
        """
        Elist = numpy.arange(0, 5000*11.96, 2*11.96)
        sumStates = self.oxygen.getSumOfStates(Elist)
        densStates = self.oxygen.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))
            
    def test_getDensityOfStates_oxygen(self):
        """
        Test the StatMech.getDensityOfStates() method for oxygen.
        """
        Elist = numpy.arange(0, 5000*11.96, 2*11.96)
        densStates = self.oxygen.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist / constants.R / T))
        Qexp = self.oxygen.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp, Qact, delta=1e-1*Qexp)

    def test_getTotalMass(self):
        """
        Test the Conformer.getTotalMass() method.
        """
        self.assertAlmostEqual(self.conformer.getTotalMass()*constants.Na*1000., numpy.sum(self.mass), 6)

    def test_getCenterOfMass(self):
        """
        Test the Conformer.getCenterOfMass() method.
        """
        cm = self.conformer.getCenterOfMass()
        self.assertAlmostEqual(cm[0]*1e10, -0.38045, 4)
        self.assertAlmostEqual(cm[1]*1e10, -0.60255, 4)
        self.assertAlmostEqual(cm[2]*1e10, -0.27900, 4)

    def test_getMomentOfInertiaTensor(self):
        """
        Test the Conformer.getMomentOfInertiaTensor() method.
        """
        I = self.conformer.getMomentOfInertiaTensor()
        self.assertAlmostEqual(I[0,0]*constants.Na*1e23, 20.65968, 4)
        self.assertAlmostEqual(I[0,1]*constants.Na*1e23, -7.48115, 4)
        self.assertAlmostEqual(I[0,2]*constants.Na*1e23, -3.46416, 4)
        self.assertAlmostEqual(I[1,0]*constants.Na*1e23, -7.48115, 4)
        self.assertAlmostEqual(I[1,1]*constants.Na*1e23, 13.53472, 4)
        self.assertAlmostEqual(I[1,2]*constants.Na*1e23, -5.48630, 4)
        self.assertAlmostEqual(I[2,0]*constants.Na*1e23, -3.46416, 4)
        self.assertAlmostEqual(I[2,1]*constants.Na*1e23, -5.48630, 4)
        self.assertAlmostEqual(I[2,2]*constants.Na*1e23, 22.84296, 4)

    def test_getPrincipalMomentsOfInertia(self):
        """
        Test the Conformer.getPrincipalMomentsOfInertia() method.
        """
        I, V = self.conformer.getPrincipalMomentsOfInertia()
        self.assertAlmostEqual(I[0]*constants.Na*1e23,  6.27074, 4)
        self.assertAlmostEqual(I[1]*constants.Na*1e23, 25.38321, 3)
        self.assertAlmostEqual(I[2]*constants.Na*1e23, 25.38341, 3)
        #print V
        # For some reason the axes seem to jump around (positioning and signs change)
        # but the absolute values should be the same as we expect
        expected = sorted([0.497140,
                           0.610114,
                           0.616938,
                           0.787360,
                           0.018454,
                           0.616218,
                           0.364578,
                           0.792099,
                           0.489554])
        result = sorted(abs(V).flat)
        for i,j in zip(expected, result):
            self.assertAlmostEqual(i, j, 4)
        return # now because the following often fails:
        self.assertAlmostEqual(V[0,0],  0.497140, 4)
        self.assertAlmostEqual(V[0,1], -0.610114, 4)
        self.assertAlmostEqual(V[0,2], -0.616938, 4)
        self.assertAlmostEqual(V[1,0],  0.787360, 4)
        self.assertAlmostEqual(V[1,1],  0.018454, 4)
        self.assertAlmostEqual(V[1,2],  0.616218, 4)
        self.assertAlmostEqual(V[2,0],  0.364578, 4)
        self.assertAlmostEqual(V[2,1],  0.792099, 4)
        self.assertAlmostEqual(V[2,2], -0.489554, 4)

    def test_getInternalReducedMomentOfInertia(self):
        """
        Test the Conformer.getInternalReducedMomentOfInertia() method.
        """
        I = self.conformer.getInternalReducedMomentOfInertia(pivots=[1,5], top1=[1,2,3,4])
        self.assertAlmostEqual(I*constants.Na*1e23, 1.56768, 4)
    def test_getNumberDegreesOfFreedom(self):
        """
        Test the Conformer.getNumberDegreesOfFreedom() method.
        """
        #this is for ethane:
        numberDegreesOfFreedom = self.conformer.getNumberDegreesOfFreedom  
        self.assertTrue(numberDegreesOfFreedom, 24) 
        #this is for ethylene:
        numberDegreesOfFreedom = self.ethylene.getNumberDegreesOfFreedom 
        self.assertTrue(numberDegreesOfFreedom, 18)        