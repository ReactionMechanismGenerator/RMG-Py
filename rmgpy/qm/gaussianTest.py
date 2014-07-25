#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import os
import numpy as np

from rmgpy import getPath
from rmgpy.qm.main import QMCalculator
from rmgpy.molecule import Molecule
from rmgpy.qm.gaussian import GaussianMolPM3, GaussianMolPM6


gaussEnv = os.getenv('GAUSS_EXEDIR') or os.getenv('g09root') or os.getenv('g03root') or ""
# GAUSS_EXEDIR may be a list like "path1:path2:path3"
for possibleDir in gaussEnv.split(':'):
	if os.path.exists(os.path.join(possibleDir , 'g09')):
		executablePath = os.path.join(possibleDir , 'g09')
		break
	elif os.path.exists(os.path.join(possibleDir , 'g03')):
		executablePath = os.path.join(possibleDir , 'g03')
		break
else:
	executablePath = os.path.join(gaussEnv , '(g03 or g09)')

qm = QMCalculator()
qm.settings.software = 'gaussian'
RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
qm.settings.fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles')
qm.settings.scratchDirectory = None
qm.settings.onlyCyclics = False
qm.settings.maxRadicalNumber = 0
	
mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')

class TestGaussianMolPM3(unittest.TestCase):
	"""
	Contains unit tests for the Geometry class.
	"""

	@unittest.skipIf(os.path.exists(executablePath)==False, "Gaussian not found. Try resetting your environment variables if you want to use it.")
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""

		if not os.path.exists(qm.settings.fileStore):
			os.makedirs(qm.settings.fileStore)

		self.qmmol1 = GaussianMolPM3(mol1, qm.settings)

	def testGenerateThermoData(self):
		"""
		Test that generateThermoData() works correctly.
		"""
		try:
			fileList = os.listdir(self.qmmol1.settings.fileStore)
			for fileName in fileList:
				os.remove(os.path.join(self.qmmol1.settings.fileStore, fileName))
		except OSError:
			pass

		self.qmmol1.generateThermoData()
		result = self.qmmol1.qmData

		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM GaussianMolPM3 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)

		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place

	def testLoadThermoData(self):
		"""
		Test that generateThermoData() can load thermo from a previous run.

		Check that it loaded, and the values are the same as above.
		"""

		self.qmmol1.generateThermoData()
		result = self.qmmol1.qmData

		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM GaussianMolPM3 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		self.assertAlmostEqual(result.energy.value_si, 169708.01906637018, 1)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)

		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place
		
class TestGaussianMolPM6(unittest.TestCase):
	"""
	Contains unit tests for the Geometry class.
	"""

	@unittest.skipIf(os.path.exists(executablePath)==False, "Gaussian not found. Try resetting your environment variables if you want to use it.")
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""

		if not os.path.exists(qm.settings.fileStore):
			os.makedirs(qm.settings.fileStore)

		self.qmmol1 = GaussianMolPM6(mol1, qm.settings)

	def testGenerateThermoData(self):
		"""
		Test that generateThermoData() works correctly.
		"""
		try:
			fileList = os.listdir(self.qmmol1.settings.fileStore)
			for fileName in fileList:
				os.remove(os.path.join(self.qmmol1.settings.fileStore, fileName))
		except OSError:
			pass

		self.qmmol1.generateThermoData()
		result = self.qmmol1.qmData

		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM GaussianMolPM6 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)

		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place

	def testLoadThermoData(self):
		"""
		Test that generateThermoData() can load thermo from a previous run.

		Check that it loaded, and the values are the same as above.
		"""

		self.qmmol1.generateThermoData()
		result = self.qmmol1.qmData

		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM GaussianMolPM6 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		self.assertAlmostEqual(result.energy.value_si, 169708.01906637018, 1)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)

		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place



################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )