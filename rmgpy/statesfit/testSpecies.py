#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

import math
import numpy
import logging

from rmgpy.chem.molecule import Molecule
from rmgpy.chem.states import *
from rmgpy.chem.thermo import ThermoGAModel
import rmgpy.chem.constants as constants

from rmgpy.statesfit.fit import fitSpectralDataToHeatCapacity

################################################################################

class GroupFrequencyCheck(unittest.TestCase):

    Tlist = numpy.arange(300.0, 2001.0, 50.0, numpy.float64)

    def testAcetyl(self):
        """
        Tests the acetyl radical.
        """

        molecule = Molecule(SMILES='C[C]=O')
        linear = False
        numVibrations = 11
        numRotors = 1
        thermoData = ThermoGAModel(
            Tdata=numpy.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], numpy.float64),
            Cpdata=numpy.array([12.06*4.184, 14.11*4.184, 16.08*4.184, 17.82*4.184, 20.67*4.184, 22.81*4.184, 26.16*4.184], numpy.float64),
            H298=-1.95*4184,
            S298=63.74*4.184,
        )

        # Subtract out contributions to heat capacity from the group frequencies
        Tlist = GroupFrequencyCheck.Tlist
        Cv = numpy.array([thermoData.getHeatCapacity(T) / constants.R for T in Tlist], numpy.float64)
        ho = HarmonicOscillator(frequencies=[2750.0, 2800.0, 2850.0, 1350.0, 1500.0, 750.0, 1050.0, 1375.0, 1855.0, 455.0])
        Cv -= ho.getHeatCapacities(Tlist) / constants.R
        # Subtract out translational modes
        Cv -= 1.5
        # Subtract out external rotational modes
        Cv -= (1.5 if not linear else 1.0)
        # Subtract out PV term (Cp -> Cv)
        Cv -= 1.0
        # Check that all Cv values are still positive (should we do this?)
        for C in Cv:
            self.assertTrue(C > 0)

        # Fit remaining frequencies and hindered rotors to the heat capacity data
        vib, hind = fitSpectralDataToHeatCapacity(molecule, Tlist, Cv, numVibrations - len(ho.frequencies), numRotors)
        print vib, hind

        model = StatesModel()
        model.modes.append(HarmonicOscillator(frequencies=[freq for freq, degen in vib for d in range(degen)]))
        for freq, barr, degen in hind:
            inertia = (barr*constants.c*100.0*constants.h) / (2 * (freq*constants.c*100.0)**2)
            barrier = barr*constants.c*100.0*constants.h*constants.Na
            for d in range(degen):
                model.modes.append(HinderedRotor(inertia=inertia, barrier=barrier, symmetry=1))

        Cp_fit = model.getHeatCapacities(Tlist) / constants.R - 1.0
        Cp_data = Cv
        rms = numpy.sqrt(numpy.sum((Cp_fit - Cp_data) * (Cp_fit - Cp_data))) / len(Tlist)
        print Cp_data
        print Cp_fit
        print rms
        self.assertTrue(rms < 0.05, "RMS error per point of heat capacity = %g" % rms)

        # vib = [2761.6297]
        # hind = [(448.9323, 1941.5089, 1)]

    def testOxygen(self):
        """
        Tests the oxygen molecule.
        """

        molecule = Molecule(SMILES='[O][O]')
        linear = True
        numVibrations = 1
        numRotors = 0

        ho = HarmonicOscillator(frequencies=[1482.2752])
        Tlist = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)
        Cv = numpy.array([ho.getHeatCapacity(T) / constants.R for T in Tlist], numpy.float64)

        # Fit remaining frequencies and hindered rotors to the heat capacity data
        vib, hind = fitSpectralDataToHeatCapacity(molecule, Tlist, Cv, numVibrations, numRotors)
        print vib, hind
        self.assertTrue(len(vib) == 1)
        self.assertTrue(len(hind) == 0)
        self.assertAlmostEqual(vib[0][0] / ho.frequencies[0], 1.0, 3)

        model = StatesModel()
        model.modes.append(HarmonicOscillator(frequencies=[freq for freq, degen in vib for d in range(degen)]))
        for freq, barr, degen in hind:
            inertia = (barr*constants.c*100.0*constants.h) / (2 * (freq*constants.c*100.0)**2)
            barrier = barr*constants.c*100.0*constants.h*constants.Na
            for d in range(degen):
                model.modes.append(HinderedRotor(inertia=inertia, barrier=barrier, symmetry=1))

        Cp_fit = model.getHeatCapacities(Tlist) / constants.R - 1.0
        Cp_data = Cv
        rms = numpy.sqrt(numpy.sum((Cp_fit - Cp_data) * (Cp_fit - Cp_data))) / len(Tlist)
        print Cp_data
        print Cp_fit
        print rms
        self.assertTrue(rms < 0.005, "RMS error per point of heat capacity = %g" % rms)

    def testAcetylperoxy(self):
        """
        Tests the acetylperoxy radical.
        """

        molecule = Molecule(SMILES='CC(=O)O[O]')
        linear = False
        numVibrations = 16
        numRotors = 2
        thermoData = ThermoGAModel(
            Tdata=numpy.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], numpy.float64),
            Cpdata=numpy.array([19.91*4.184, 23.07*4.184, 25.85*4.184, 28.15*4.184, 31.56*4.184, 33.90*4.184, 37.28*4.184], numpy.float64),
            H298=-36.67*4184,
            S298=76.21*4.184
        )

        # Subtract out contributions to heat capacity from the group frequencies
        Tlist = GroupFrequencyCheck.Tlist
        Cv = numpy.array([thermoData.getHeatCapacity(T) / constants.R for T in Tlist], numpy.float64)
        ho = HarmonicOscillator(frequencies=[2750.0, 2800.0, 2850.0, 1350.0, 1500.0, 750.0, 1050.0, 1375.0, 492.5, 1135.0])
        Cv -= ho.getHeatCapacities(Tlist) / constants.R
        # Subtract out translational modes
        Cv -= 1.5
        # Subtract out external rotational modes
        Cv -= (1.5 if not linear else 1.0)
        # Subtract out PV term (Cp -> Cv)
        Cv -= 1.0
        # Check that all Cv values are still positive (should we do this?)
        for C in Cv:
            self.assertTrue(C > 0)

        # Fit remaining frequencies and hindered rotors to the heat capacity data
        vib, hind = fitSpectralDataToHeatCapacity(molecule, Tlist, Cv, numVibrations - len(ho.frequencies), numRotors)
        print vib, hind

        model = StatesModel()
        model.modes.append(HarmonicOscillator(frequencies=[freq for freq, degen in vib for d in range(degen)]))
        for freq, barr, degen in hind:
            inertia = (barr*constants.c*100.0*constants.h) / (2 * (freq*constants.c*100.0)**2)
            barrier = barr*constants.c*100.0*constants.h*constants.Na
            for d in range(degen):
                model.modes.append(HinderedRotor(inertia=inertia, barrier=barrier, symmetry=1))

        Cp_fit = model.getHeatCapacities(Tlist) / constants.R - 1.0
        Cp_data = Cv
        rms = numpy.sqrt(numpy.sum((Cp_fit - Cp_data) * (Cp_fit - Cp_data))) / len(Tlist)
        print Cp_data
        print Cp_fit
        print rms
        self.assertTrue(rms < 0.05, "RMS error per point of heat capacity = %g" % rms)

        # vib = [(586.74145,4), (1444.6377, 2)]
        # hind = [(40.0, 150.0, 1), (1120.2859, 2173.5399, 1)]

    def testHydroperoxylvinoxy(self):
        """
        Tests the hydroperoxylvinoxy radical.
        """

        molecule = Molecule(SMILES='[CH2]C(=O)OO')
        linear = False
        numVibrations = 16
        numRotors = 2
        thermoData = ThermoGAModel(
            Tdata=numpy.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], numpy.float64),
            Cpdata=numpy.array([21.04*4.184, 25.44*4.184, 28.98*4.184, 31.67*4.184, 35.20*4.184, 37.28*4.184, 39.60*4.184], numpy.float64),
            H298=-34.67*4184,
            S298=73.73*4.184,
        )

        # Subtract out contributions to heat capacity from the group frequencies
        Tlist = GroupFrequencyCheck.Tlist
        Cv = numpy.array([thermoData.getHeatCapacity(T) / constants.R for T in Tlist], numpy.float64)
        ho = HarmonicOscillator(frequencies=[3000.0, 3100.0, 440.0, 815.0, 1455.0, 3615.0, 1310.0, 387.5, 850.0])
        Cv -= ho.getHeatCapacities(Tlist) / constants.R
        # Subtract out translational modes
        Cv -= 1.5
        # Subtract out external rotational modes
        Cv -= (1.5 if not linear else 1.0)
        # Subtract out PV term (Cp -> Cv)
        Cv -= 1.0
        # Check that all Cv values are still positive (should we do this?)
        for C in Cv:
            self.assertTrue(C > 0)

        # Fit remaining frequencies and hindered rotors to the heat capacity data
        vib, hind = fitSpectralDataToHeatCapacity(molecule, Tlist, Cv, numVibrations - len(ho.frequencies), numRotors)
        print vib, hind

        model = StatesModel()
        model.modes.append(HarmonicOscillator(frequencies=[freq for freq, degen in vib for d in range(degen)]))
        for freq, barr, degen in hind:
            inertia = (barr*constants.c*100.0*constants.h) / (2 * (freq*constants.c*100.0)**2)
            barrier = barr*constants.c*100.0*constants.h*constants.Na
            for d in range(degen):
                model.modes.append(HinderedRotor(inertia=inertia, barrier=barrier, symmetry=1))

        Cp_fit = model.getHeatCapacities(Tlist) / constants.R - 1.0
        Cp_data = Cv
        rms = numpy.sqrt(numpy.sum((Cp_fit - Cp_data) * (Cp_fit - Cp_data))) / len(Tlist)
        print Cp_data
        print Cp_fit
        print rms
        self.assertTrue(rms < 0.05, "RMS error per point of heat capacity = %g" % rms)

        # vib = [(660.71709,3), (1167.9948, 3)]
        # hind = [(100.0, 779.05087, 1), (300.0, 2311.3714, 2)]

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
