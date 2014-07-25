#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import os

from rmgpy import getPath
from rmgpy.qm.main import QMSettings, QMCalculator
from rmgpy.molecule import Molecule
import qmdata

class TestQMSettings(unittest.TestCase):
	"""
	Contains unit tests for the QMSettings class.
	"""
	
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""
		RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
		
		self.settings1 = QMSettings()
		self.settings1.software = 'mopac'
		self.settings1.method = 'pm3'
		self.settings1.fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles')
		self.settings1.scratchDirectory = None
		self.settings1.onlyCyclics = False
		self.settings1.maxRadicalNumber = 0
		
		self.settings2 = QMSettings()

	def testCheckAllSet(self):
		"""
		Test that checkAllSet() works correctly.
		"""
		try:
			self.settings1.checkAllSet()
		except AssertionError:
			self.fail("checkAllSet() raised unexpected AssertionError.")
		
		with self.assertRaises(AssertionError):
			self.settings2.checkAllSet()

class TestQMCalculator(unittest.TestCase):
	"""
	Contains unit tests for the QMSettings class.
	"""
	
	mopacEnv = os.getenv('MOPAC_DIR', default="/opt/mopac")
	if os.path.exists(os.path.join(mopacEnv , 'MOPAC2012.exe')):
		mopExecutablePath = os.path.join(mopacEnv , 'MOPAC2012.exe')
	elif os.path.exists(os.path.join(mopacEnv , 'MOPAC2009.exe')):
		mopExecutablePath = os.path.join(mopacEnv , 'MOPAC2009.exe')
	else:
		mopExecutablePath = os.path.join(mopacEnv , '(MOPAC 2009 or 2012)')
	
	gaussEnv = os.getenv('GAUSS_EXEDIR') or os.getenv('g09root') or os.getenv('g03root') or ""
	# GAUSS_EXEDIR may be a list like "path1:path2:path3"
	for possibleDir in gaussEnv.split(':'):
		if os.path.exists(os.path.join(possibleDir , 'g09')):
			gaussExecutablePath = os.path.join(possibleDir , 'g09')
			break
		elif os.path.exists(os.path.join(possibleDir , 'g03')):
			gaussExecutablePath = os.path.join(possibleDir , 'g03')
			break
	else:
		gaussExecutablePath = os.path.join(gaussEnv , '(g03 or g09)')
	
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""
		RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
		
		fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles')
		if not os.path.exists(fileStore):
			os.makedirs(fileStore)
		
		self.mop1 = QMCalculator(fileStore=fileStore)
		self.mop1.settings.software = 'mopac'
		self.mop1.settings.method = 'pm3'
		self.mop1.settings.onlyCyclics = False
		self.mop1.settings.maxRadicalNumber = 0
				
		self.mop2 = QMCalculator()
		self.mop2.settings.software = 'mopac'
		self.mop2.settings.method = 'pm6'
		self.mop2.settings.onlyCyclics = False
		self.mop2.settings.maxRadicalNumber = 0

		
		self.mop3 = QMCalculator(fileStore=fileStore)
		self.mop3.settings.software = 'mopac'
		self.mop3.settings.method = 'pm7'
		self.mop3.settings.onlyCyclics = False
		self.mop3.settings.maxRadicalNumber = 0
		
		self.mop4 = QMCalculator(fileStore=fileStore)
		self.mop4.settings.software = 'mopac'
		self.mop4.settings.method = 'pm8'
		self.mop4.settings.onlyCyclics = False
		self.mop4.settings.maxRadicalNumber = 0
		
		self.gauss1 = QMCalculator()
		self.gauss1.settings.software = 'gaussian'
		self.gauss1.settings.method = 'pm3'
		self.gauss1.settings.onlyCyclics = False
		self.gauss1.settings.maxRadicalNumber = 0	
		
		self.gauss2 = QMCalculator(fileStore=fileStore)
		self.gauss2.settings.software = 'gaussian'
		self.gauss2.settings.method = 'pm6'
		self.gauss2.settings.onlyCyclics = False
		self.gauss2.settings.maxRadicalNumber = 0
		
		self.gauss3 = QMCalculator(fileStore=fileStore)
		self.gauss3.settings.software = 'gaussian'
		self.gauss3.settings.method = 'pm7'
		self.gauss3.settings.onlyCyclics = False
		self.gauss3.settings.maxRadicalNumber = 0
		
		self.molpro1 = QMCalculator(fileStore=fileStore)
		self.molpro1.settings.software = 'molpro'
		self.molpro1.settings.method = 'mp2'
		self.molpro1.settings.onlyCyclics = False
		self.molpro1.settings.maxRadicalNumber = 0
		
		self.qmmol1 = QMCalculator(fileStore=fileStore)
		self.qmmol1.RMG_bin_path = os.path.join(RMGpy_path, 'testing', 'qm', 'bin')
		
		self.qmmol2 = QMCalculator(fileStore=fileStore)
		self.qmmol2.RMG_bin_path = os.path.join(RMGpy_path, 'testing', 'hexadiene', 'input.py')

	def testSetDefaultOutputDirectory(self):
		"""
		Test that setDefaultOutputDirectory() works correctly.
		"""
		
		self.assertIsNotNone(self.mop1.settings.fileStore)
		self.assertIsNotNone(self.mop3.settings.fileStore)
		self.assertIsNotNone(self.gauss2.settings.fileStore)
		
		self.assertIsNone(self.mop2.settings.fileStore)
		self.assertIsNone(self.gauss1.settings.fileStore)
		
		self.assertIsNone(self.mop1.settings.scratchDirectory)
		self.assertIsNone(self.mop2.settings.scratchDirectory)
		self.assertIsNone(self.mop3.settings.scratchDirectory)
		self.assertIsNone(self.gauss1.settings.scratchDirectory)
		self.assertIsNone(self.gauss2.settings.scratchDirectory)
		
		# Now set the default directories for those not set
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..')
		self.mop1.setDefaultOutputDirectory(outputDirectory)
		self.mop2.setDefaultOutputDirectory(outputDirectory)
		self.mop3.setDefaultOutputDirectory(outputDirectory)
		self.gauss1.setDefaultOutputDirectory(outputDirectory)
		self.gauss2.setDefaultOutputDirectory(outputDirectory)
		
		self.assertIsNotNone(self.mop1.settings.fileStore)
		self.assertIsNotNone(self.mop2.settings.fileStore)
		self.assertIsNotNone(self.mop3.settings.fileStore)
		self.assertIsNotNone(self.gauss1.settings.fileStore)
		self.assertIsNotNone(self.gauss2.settings.fileStore)
		self.assertIsNotNone(self.mop1.settings.scratchDirectory)
		self.assertIsNotNone(self.mop2.settings.scratchDirectory)
		self.assertIsNotNone(self.mop3.settings.scratchDirectory)
		self.assertIsNotNone(self.gauss1.settings.scratchDirectory)
		self.assertIsNotNone(self.gauss2.settings.scratchDirectory)

	def testCheckPaths(self):
		"""
		Test that checkPaths() works correctly.
		"""
		
		with self.assertRaises(Exception):
			self.qmmol1.checkPaths()
			self.qmmol2.checkPaths()
	
	def testInitialize(self):
		"""
		Test that initialize() works correctly.
		"""
		
		# Now set the default directories for those not set
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..')
		self.mop1.setDefaultOutputDirectory(outputDirectory)
		self.mop2.setDefaultOutputDirectory(outputDirectory)
		self.mop3.setDefaultOutputDirectory(outputDirectory)
		self.gauss1.setDefaultOutputDirectory(outputDirectory)
		self.gauss2.setDefaultOutputDirectory(outputDirectory)
		
		try:
			self.mop1.initialize()
			self.mop2.initialize()
			self.mop3.initialize()
			self.gauss1.initialize()
			self.gauss2.initialize()
		except AssertionError:
			self.fail("checkAllSet() raised unexpected AssertionError.")
		except Exception:
			self.fail("initialize() raised Exception. Output file paths not correctly set.")
	
	def testGetThermoData(self):
		"""
		Test that getThermoData() fails when expected.
		"""
		outputDirectory = os.path.join(self.mop4.settings.fileStore, '..')
		self.mop4.setDefaultOutputDirectory(outputDirectory)
		self.gauss3.setDefaultOutputDirectory(outputDirectory)
		self.molpro1.setDefaultOutputDirectory(outputDirectory)
		
		mol = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		
		with self.assertRaises(Exception):
			self.mop4.getThermoData(mol)
			self.gauss3.getThermoData(mol)
			self.molpro1.getThermoData(mol)
		
	@unittest.skipIf(os.path.exists(mopExecutablePath)==False, "If MOPAC installed, try checking your environment variables.")
	def testGetThermoDataMopac(self):
		"""
		Test that getThermoData() works correctly.
		"""
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..')
		self.mop1.setDefaultOutputDirectory(outputDirectory)
		self.mop2.setDefaultOutputDirectory(outputDirectory)
		self.mop3.setDefaultOutputDirectory(outputDirectory)
		
		mol = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		
		fileList = os.listdir(self.mop1.settings.fileStore)
		for fileName in fileList:
			os.remove(os.path.join(self.mop1.settings.fileStore, fileName))
		thermo1 = self.mop1.getThermoData(mol)
		
		fileList = os.listdir(self.mop2.settings.fileStore)
		for fileName in fileList:
			os.remove(os.path.join(self.mop2.settings.fileStore, fileName))
		thermo2 = self.mop2.getThermoData(mol)
		
		fileList = os.listdir(self.mop3.settings.fileStore)
		for fileName in fileList:
			os.remove(os.path.join(self.mop3.settings.fileStore, fileName))
		thermo3 = self.mop3.getThermoData(mol)
			
		self.assertTrue(thermo1.comment.startswith('QM MopacMolPM3'))
		self.assertTrue(thermo2.comment.startswith('QM MopacMolPM6'))
		self.assertTrue(thermo3.comment.startswith('QM MopacMolPM7'))
		
		self.assertAlmostEqual(thermo1.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo1.S298.value_si, 334.5007584, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo2.H298.value_si, 167704.4270, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo2.S298.value_si, 338.0999241, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo3.H298.value_si, 166168.8571, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo3.S298.value_si, 336.3330406, 1) # to 1 decimal place
		
	@unittest.skipIf(os.path.exists(gaussExecutablePath)==False, "If GAUSSIAN installed, try checking your environment variables.")
	def testGetThermoDataGaussian(self):
		"""
		Test that getThermoData() works correctly.
		"""
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..')
		self.gauss1.setDefaultOutputDirectory(outputDirectory)
		self.gauss2.setDefaultOutputDirectory(outputDirectory)
		
		mol = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		
		thermo1 = self.gauss1.getThermoData(mol)
		fileList = os.listdir(self.gauss1.settings.fileStore)
		for fileName in fileList:
			os.remove(os.path.join(self.gauss1.settings.fileStore, fileName))
		
		thermo2 = self.gauss2.getThermoData(mol)
		fileList = os.listdir(self.gauss2.settings.fileStore)
		for fileName in fileList:
			os.remove(os.path.join(self.gauss2.settings.fileStore, fileName))
		
		self.assertTrue(thermo1.comment.startswith('QM GaussianMolPM3'))
		self.assertTrue(thermo1.comment.startswith('QM GaussianMolPM6'))
		
		self.assertAlmostEqual(thermo1.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo1.S298.value_si, 334.5007584, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo2.H298.value_si, 167704.4270, 1) # to 1 decimal place
		self.assertAlmostEqual(thermo2.S298.value_si, 338.0999241, 1) # to 1 decimal place

################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )