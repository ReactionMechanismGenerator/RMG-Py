#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import distutils.spawn
import itertools
import logging
import numpy as np
import os

from rmgpy import getPath
from rmgpy.qm.main import QMCalculator
from rmgpy.molecule import Molecule
from rmgpy.qm.gaussian import GaussianMolPM3, GaussianMolPM6


executablesToTry = ('g09', 'g03')

for exe in executablesToTry:
    try:
        executablePath = distutils.spawn.find_executable(exe)
    except:
        executablePath = None
    if executablePath is not None:
        break
else:  # didn't break
    logging.debug("Did not find Gaussian on path, checking if it exists in a declared GAUSS_EXEDIR, g09root or g03root...")
    gaussEnv = os.getenv('GAUSS_EXEDIR') or os.getenv('g09root') or os.getenv('g03root') or ""
    possibleDirs = gaussEnv.split(':')# GAUSS_EXEDIR may be a list like "path1:path2:path3"
    for exe, possibleDir in itertools.product(executablesToTry, possibleDirs):
        executablePath = os.path.join(possibleDir, exe)
        if os.path.exists(executablePath):
            break
    else:  # didn't break
        executablePath = os.path.join(gaussEnv , '(Gaussian 2003 or 2009)')
	
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
		RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
		
		qm = QMCalculator(software = 'gaussian',
						  method = 'pm3',
						  fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles'),
						  scratchDirectory = os.path.join(RMGpy_path, 'testing', 'qm', 'QMscratch'),
						  )

		
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
			self.assertAlmostEqual(result.molecularMass.value, 128.0626, 3)

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
		if result.molecularMass.units=='amu':
			self.assertAlmostEqual(result.molecularMass.value, 128.0626, 3)
		
class TestGaussianMolPM6(unittest.TestCase):
	"""
	Contains unit tests for the Geometry class.
	"""

	@unittest.skipIf(os.path.exists(executablePath)==False, "Gaussian not found. Try resetting your environment variables if you want to use it.")
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""
		RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
		
		qm = QMCalculator(software = 'gaussian',
						  method = 'pm6',
						  fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles'),
						  scratchDirectory = os.path.join(RMGpy_path, 'testing', 'qm', 'QMscratch'),
						  )

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
			self.assertAlmostEqual(result.molecularMass.value, 128.0626, 3)

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
		if result.molecularMass.units=='amu':
			self.assertAlmostEqual(result.molecularMass.value, 128.0626, 3)


################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )