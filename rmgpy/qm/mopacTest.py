#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import os
import numpy as np

from rmgpy import getPath
from rmgpy.qm.main import QMCalculator
from rmgpy.molecule import Molecule
from rmgpy.qm.mopac import MopacMolPM3, MopacMolPM6, MopacMolPM7


mopacEnv = os.getenv('MOPAC_DIR', default="/opt/mopac")
if os.path.exists(os.path.join(mopacEnv , 'MOPAC2012.exe')):
	executablePath = os.path.join(mopacEnv , 'MOPAC2012.exe')
elif os.path.exists(os.path.join(mopacEnv , 'MOPAC2009.exe')):
	executablePath = os.path.join(mopacEnv , 'MOPAC2009.exe')
else:
	executablePath = os.path.join(mopacEnv , '(MOPAC 2009 or 2012)')
	
qm = QMCalculator()
qm.settings.software = 'mopac'
RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
qm.settings.fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles')
qm.settings.scratchDirectory = None
qm.settings.onlyCyclics = False
qm.settings.maxRadicalNumber = 0

mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')

class TestMopacMolPM3(unittest.TestCase):
	"""
	Contains unit tests for the Geometry class.
	"""
	
	@unittest.skipIf(os.path.exists(executablePath)==False, "MOPAC not found. Try resetting your environment variables if you want to use it.")
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""
		
		if not os.path.exists(qm.settings.fileStore):
			os.makedirs(qm.settings.fileStore)
		
		self.qmmol1 = MopacMolPM3(mol1, qm.settings)
	
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
		
		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM3 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)
		
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 0) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place
	
	def testLoadThermoData(self):
		"""
		Test that generateThermoData() can load thermo from a previous run.
		
		Check that it loaded, and the values are the same as above.
		"""
		
		self.qmmol1.generateThermoData()
		result = self.qmmol1.qmData
		
		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM3 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)
		
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 0) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place
			
class TestMopacMolPM6(unittest.TestCase):
	"""
	Contains unit tests for the Geometry class.
	"""

	@unittest.skipIf(os.path.exists(executablePath)==False, "MOPAC not found. Try resetting your environment variables if you want to use it.")
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""

		if not os.path.exists(qm.settings.fileStore):
			os.makedirs(qm.settings.fileStore)

		self.qmmol1 = MopacMolPM6(mol1, qm.settings)

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
		
		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM6 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)
		
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 167704.4270, 0) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 338.0999241, 1) # to 1 decimal place
	
	def testLoadThermoData(self):
		"""
		Test that generateThermoData() can load thermo from a previous run.
		
		Check that it loaded, and the values are the same as above.
		"""
		
		self.qmmol1.generateThermoData()
		result = self.qmmol1.qmData
		
		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM6 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)
		
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 167704.0681, 0) # to 0 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 338.0999241, 1) # to 1 decimal place

class TestMopacMolPM7(unittest.TestCase):
	"""
	Contains unit tests for the Geometry class.
	"""
	
	@unittest.skipIf(os.path.exists(executablePath)==False, "MOPAC not found. Try resetting your environment variables if you want to use it.")
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""

		if not os.path.exists(qm.settings.fileStore):
			os.makedirs(qm.settings.fileStore)

		mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		self.qmmol1 = MopacMolPM7(mol1, qm.settings)

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
		
		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM7 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)
		
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 166168.9863, 0) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 336.3330406, 1) # to 1 decimal place
	
	def testLoadThermoData(self):
		"""
		Test that generateThermoData() can load thermo from a previous run.
		
		Check that it loaded, and the values are the same as above.
		"""
		
		self.qmmol1.generateThermoData()
		result = self.qmmol1.qmData
		
		self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM7 calculation'))
		self.assertEqual(result.numberOfAtoms, 18)
		self.assertIsInstance(result.atomicNumbers, np.ndarray)
		if result.molecularMass.units=='amu':
			self.assertEqual(result.molecularMass.value, 128.173)
		
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 166168.8571, 0) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 336.3330406, 1) # to 1 decimal place	

		
################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )