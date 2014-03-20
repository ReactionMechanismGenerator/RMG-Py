#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import os
import numpy as np

from rmgpy.qm.main import QMCalculator
from rmgpy.molecule import Molecule
from rmgpy.qm.gaussian import GaussianMolPM3
import qmdata

class TestGaussianMolPM3(unittest.TestCase):
	"""
	Contains unit tests for the GaussianMolPM3 class.
	"""
	gaussEnv = os.getenv('GAUSS_EXEDIR') or os.getenv('g09root') or os.getenv('g03root') or ""
	if os.path.exists(os.path.join(gaussEnv , 'g09')):
	    executablePath = os.path.join(gaussEnv , 'g09')
	elif os.path.exists(os.path.join(gaussEnv , 'g03')):
	    executablePath = os.path.join(gaussEnv , 'g03')
	else:
	    executablePath = os.path.join(gaussEnv , '(g03 or g09)')
	
	@unittest.skipIf(os.path.exists(executablePath)==False, "Gaussian not found. Try resetting your environment variables if you want to use it.")
	
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""
		qm = QMCalculator()
		qm.settings.software = 'gaussian'
		qm.settings.fileStore = os.path.join(os.getenv('RMGpy'), 'QMfiles')
		qm.settings.scratchDirectory = None
		qm.settings.onlyCyclics = False
		qm.settings.maxRadicalNumber = 0

		if not os.path.exists(qm.settings.fileStore):
			os.mkdir(qm.settings.fileStore)

		mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		self.qmmol1 = GaussianMolPM3(mol1, qm.settings)

	def testGenerateQMData(self):
		"""
		Test that generateQMData() works correctly.
		"""
		try:
			# Remove the output file so we can test the gaussian script
			os.remove(self.qmmol1.outputFilePath)
		except OSError:
			pass

		result = self.qmmol1.generateQMData()
		self.assertTrue(os.path.exists(self.qmmol1.inputFilePath))
		self.assertTrue(os.path.exists(self.qmmol1.outputFilePath))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		self.assertAlmostEqual(result.energy.value_si, 169708.01906637018, 1)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)

	def testGenerateThermoData(self):
		"""
		Test that generateThermoData() works correctly. The testGenerateQMData should
		have tested the running of Gaussian to generate a `.thermo` file, so this will now test
		the loading reading of an existing `.thermo` file.
		"""
		self.qmmol1.generateThermoData()
		self.assertIsNotNone(self.qmmol1.thermo)
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place



################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )