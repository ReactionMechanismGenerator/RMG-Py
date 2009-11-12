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
import rmg.main as main
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
		thermoData = thermo.ThermoGAData(H298=-2.28*4184, S298=64.25*4.184,
			Cp=[12.39*4.184, 14.28*4.184, 16.26*4.184, 18.05*4.184, 21.06*4.184, 23.24*4.184, 25.14*4.184])

		frequencies0 = spectral.generateSpectralData(struct, thermoData)
		frequencies0.sort()
		frequencies = [2750.0, 2800.0, 2850.0, 1350.0, 1500.0, 750.0, 1050.0, 1375.0, 1855.0, 455.0]
		frequencies.sort()
		
		self.assertTrue(frequencies0 == frequencies)

#		import numpy
#		modes = [spectral.HarmonicOscillator(frequency=freq) for freq in frequencies0]
#		Tlist = numpy.arange(300.0, 10000.0, 100.0)
#		Cp = numpy.zeros(len(Tlist), numpy.float64)
#		for mode in modes:
#			Cp += mode.heatCapacity(Tlist)
#		import pylab
#		pylab.plot(Tlist, Cp)
#		pylab.show()

	def testAcetylperoxy(self):
		"""
		Tests the acetylperoxy radical.
		"""

		struct = structure.Structure(SMILES='CC(=O)O[O]')
		thermoData = thermo.ThermoGAData(H298=-38.47*4184, S298=75.43*4.184,
			Cp=[18.98*4.184, 22.44*4.184, 25.21*4.184, 27.46*4.184, 30.86*4.184, 33.26*4.184, 36.78*4.184])

		frequencies0 = spectral.generateSpectralData(struct, thermoData)
		frequencies0.sort()
		frequencies = [2750.0, 2800.0, 2850.0, 1350.0, 1500.0, 750.0, 1050.0, 1375.0, 492.5, 1135.0]
		frequencies.sort()

		self.assertTrue(frequencies0 == frequencies)

	def testHydroperoxylvinoxy(self):
		"""
		Tests the hydroperoxylvinoxy radical.
		"""

		struct = structure.Structure(SMILES='[CH2]C(=O)OO')
		thermoData = thermo.ThermoGAData(H298=-35.70*4184, S298=78.62*4.184,
			Cp=[20.26*4.184, 23.92*4.184, 26.85*4.184, 29.15*4.184, 32.42*4.184, 34.49*4.184, 37.20*4.184])

		frequencies0 = spectral.generateSpectralData(struct, thermoData)
		frequencies0.sort()
		frequencies = [3000.0, 3100.0, 440.0, 815.0, 1455.0, 3615.0, 1310.0, 387.5, 850.0]
		frequencies.sort()

		self.assertTrue(frequencies0 == frequencies)

	def testOxygen(self):
		"""
		Tests the oxygen molecule.
		"""

		struct = structure.Structure(SMILES='O=O')
		thermoData = thermo.ThermoGAData(H298=0.0*4184, S298=49.00*4.184,
			Cp=[7.00*4.184, 7.22*4.184, 7.44*4.184, 7.65*4.184, 8.07*4.184, 8.35*4.184, 8.72*4.184])

		frequencies0 = spectral.generateSpectralData(struct, thermoData)
		frequencies0.sort()
		frequencies = []
		frequencies.sort()

		self.assertTrue(frequencies0 == frequencies)
		
################################################################################

if __name__ == '__main__':

	# Load databases
	databasePath = '../data/RMG_database'
	loadFrequencyDatabase(databasePath)

	# Run unit tests
	unittest.main()
	