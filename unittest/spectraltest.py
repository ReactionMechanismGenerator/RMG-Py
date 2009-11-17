#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

import unittest

import sys
sys.path.append('../source')

import logging

import rmg.spectral as spectral
import rmg.structure as structure
import rmg.thermo as thermo

################################################################################

def loadFrequencyDatabase(databasePath):
	"""
	Create and load the frequencies databases.
	"""
	import os.path
	databasePath = os.path.join(databasePath, 'frequencies')

	# Create and load thermo databases
	spectral.frequencyDatabase = spectral.FrequencyDatabase()
	logging.debug('\tFrequencies database')
	spectral.frequencyDatabase.load(
		dictstr=os.path.join(databasePath, 'Dictionary.txt'),
		treestr=os.path.join(databasePath, 'Tree.txt'),
		libstr=os.path.join(databasePath, 'Library.txt'))

################################################################################

class GroupFrequencyCheck(unittest.TestCase):

	def testAcetyl(self):
		"""
		Tests the acetyl radical.
		"""

		struct = structure.Structure(SMILES='C[C]=O')
		thermoData = thermo.ThermoGAData(H298=-3.08*4184, S298=64.27*4.184,
			Cp=[12.28*4.184, 14.34*4.184, 16.30*4.184, 18.05*4.184, 20.92*4.184, 23.08*4.184, 26.39*4.184])

		spectralData = spectral.generateSpectralData(struct, thermoData)
		
		frequencies0 = [mode.frequency for mode in spectralData.modes if isinstance(mode, spectral.HarmonicOscillator)]
		frequencies0.sort()
		frequencies = [2750.0, 2800.0, 2850.0, 1350.0, 1500.0, 750.0, 1050.0, 1375.0, 1855.0, 455.0, 2761.6297]
		frequencies.sort()

		self.assertTrue(len(frequencies0) == len(frequencies))
		for freq1, freq2 in zip(frequencies0, frequencies):
			self.assertAlmostEqual(freq1 / freq2, 1.0, 1, 'Harmonic oscillator frequencies %s and %s do not match.' % (freq1, freq2))

		frequencies0 = [mode.frequency for mode in spectralData.modes if isinstance(mode, spectral.HinderedRotor)]
		frequencies0.sort()
		frequencies = [448.9323]
		frequencies.sort()
		barriers0 = [mode.barrier for mode in spectralData.modes if isinstance(mode, spectral.HinderedRotor)]
		barriers0.sort()
		barriers = [1941.5089]
		barriers.sort()

		self.assertTrue(len(frequencies0) == len(frequencies))
		self.assertTrue(len(barriers0) == len(barriers))
		for freq1, freq2 in zip(frequencies0, frequencies):
			self.assertAlmostEqual(freq1 / freq2, 1.0, 1, 'Hindered rotor frequencies %s and %s do not match.' % (freq1, freq2))
		for barr1, barr2 in zip(barriers0, barriers):
			self.assertAlmostEqual(barr1 / barr2, 1.0, 1, 'Hindered rotor barriers %s and %s do not match.' % (barr1, barr2))

	def testAcetylperoxy(self):
		"""
		Tests the acetylperoxy radical.
		"""

		struct = structure.Structure(SMILES='CC(=O)O[O]')
		thermoData = thermo.ThermoGAData(H298=-38.57*4184, S298=75.56*4.184,
			Cp=[19.48*4.184, 23.22*4.184, 26.33*4.184, 28.81*4.184, 32.39*4.184, 34.76*4.184, 37.99*4.184])

		spectralData = spectral.generateSpectralData(struct, thermoData)
		
		frequencies0 = [mode.frequency for mode in spectralData.modes if isinstance(mode, spectral.HarmonicOscillator)]
		frequencies0.sort()
		frequencies = [2750.0, 2800.0, 2850.0, 1350.0, 1500.0, 750.0, 1050.0, 1375.0, 492.5, 1135.0, 586.74145, 586.74145, 586.74145, 586.74145, 1444.6377, 1444.6377]
		frequencies.sort()

		self.assertTrue(len(frequencies0) == len(frequencies))
		for freq1, freq2 in zip(frequencies0, frequencies):
			self.assertAlmostEqual(freq1 / freq2, 1.0, 1, 'Harmonic oscillator frequencies %s and %s do not match.' % (freq1, freq2))

		frequencies0 = [mode.frequency for mode in spectralData.modes if isinstance(mode, spectral.HinderedRotor)]
		frequencies0.sort()
		frequencies = [40.0, 150.0]
		frequencies.sort()
		barriers0 = [mode.barrier for mode in spectralData.modes if isinstance(mode, spectral.HinderedRotor)]
		barriers0.sort()
		barriers = [1120.2859, 2173.5399]
		barriers.sort()

		self.assertTrue(len(frequencies0) == len(frequencies))
		self.assertTrue(len(barriers0) == len(barriers))
		for freq1, freq2 in zip(frequencies0, frequencies):
			self.assertAlmostEqual(freq1 / freq2, 1.0, 1, 'Hindered rotor frequencies %s and %s do not match.' % (freq1, freq2))
		for barr1, barr2 in zip(barriers0, barriers):
			self.assertAlmostEqual(barr1 / barr2, 1.0, 1, 'Hindered rotor barriers %s and %s do not match.' % (barr1, barr2))

	def testHydroperoxylvinoxy(self):
		"""
		Tests the hydroperoxylvinoxy radical.
		"""

		struct = structure.Structure(SMILES='[CH2]C(=O)OO')
		thermoData = thermo.ThermoGAData(H298=-32.95*4184, S298=79.25*4.184,
			Cp=[21.79*4.184, 25.53*4.184, 28.35*4.184, 30.47*4.184, 33.34*4.184, 35.12*4.184, 37.57*4.184])

		spectralData = spectral.generateSpectralData(struct, thermoData)
		
		frequencies0 = [mode.frequency for mode in spectralData.modes if isinstance(mode, spectral.HarmonicOscillator)]
		frequencies0.sort()
		frequencies = [3000.0, 3100.0, 440.0, 815.0, 1455.0, 3615.0, 1310.0, 387.5, 850.0, 660.71709, 660.71709, 660.71709, 1167.9948, 1167.9948, 1167.9948]
		frequencies.sort()

		self.assertTrue(len(frequencies0) == len(frequencies))
		for freq1, freq2 in zip(frequencies0, frequencies):
			self.assertAlmostEqual(freq1 / freq2, 1.0, 1, 'Harmonic oscillator frequencies %s and %s do not match.' % (freq1, freq2))

		frequencies0 = [mode.frequency for mode in spectralData.modes if isinstance(mode, spectral.HinderedRotor)]
		frequencies0.sort()
		frequencies = [100.0, 300.0, 300.0]
		frequencies.sort()
		barriers0 = [mode.barrier for mode in spectralData.modes if isinstance(mode, spectral.HinderedRotor)]
		barriers0.sort()
		barriers = [779.05087, 2311.3714, 2311.3714]
		barriers.sort()

		self.assertTrue(len(frequencies0) == len(frequencies))
		self.assertTrue(len(barriers0) == len(barriers))
		for freq1, freq2 in zip(frequencies0, frequencies):
			self.assertAlmostEqual(freq1 / freq2, 1.0, 1, 'Hindered rotor frequencies %s and %s do not match.' % (freq1, freq2))
		for barr1, barr2 in zip(barriers0, barriers):
			self.assertAlmostEqual(barr1 / barr2, 1.0, 1, 'Hindered rotor barriers %s and %s do not match.' % (barr1, barr2))

	def testOxygen(self):
		"""
		Tests the oxygen molecule.
		"""

		struct = structure.Structure(SMILES='O=O')
		thermoData = thermo.ThermoGAData(H298=0.0*4184, S298=49.00*4.184,
			Cp=[7.00*4.184, 7.22*4.184, 7.44*4.184, 7.65*4.184, 8.07*4.184, 8.35*4.184, 8.72*4.184])

		spectralData = spectral.generateSpectralData(struct, thermoData)

		frequencies0 = [mode.frequency for mode in spectralData.modes if isinstance(mode, spectral.HarmonicOscillator)]
		frequencies = [1482.2752]
		
		self.assertTrue(len(frequencies0) == len(frequencies))
		for freq1, freq2 in zip(frequencies0, frequencies):
			self.assertAlmostEqual(freq1 / freq2, 1.0, 1, 'Harmonic oscillator frequencies %s and %s do not match.' % (freq1, freq2))

		
################################################################################

if __name__ == '__main__':

	# Load databases
	databasePath = '../data/RMG_database'
	loadFrequencyDatabase(databasePath)

	# Run unit tests
	unittest.main()
	