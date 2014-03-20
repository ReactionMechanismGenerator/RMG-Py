#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import os
import numpy as np

from rmgpy.qm.main import QMCalculator
from rmgpy.molecule import Molecule
from rmgpy.qm.mopac import MopacMolPM3
import qmdata

class TestMopacMolPM3(unittest.TestCase):
	"""
	Contains unit tests for the Geometry class.
	"""
	mopacEnv = os.getenv('MOPAC_DIR', default="/opt/mopac")
	if os.path.exists(os.path.join(mopacEnv , 'MOPAC2012.exe')):
		executablePath = os.path.join(mopacEnv , 'MOPAC2012.exe')
	elif os.path.exists(os.path.join(mopacEnv , 'MOPAC2009.exe')):
		executablePath = os.path.join(mopacEnv , 'MOPAC2009.exe')
	else:
		executablePath = os.path.join(mopacEnv , '(MOPAC 2009 or 2012)')
		
	@unittest.skipIf(os.path.exists(executablePath)==False, "MOPAC not found. Try resetting your environment variables if you want to use it.")
	
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""
		qm = QMCalculator()
		qm.settings.software = 'mopac'
		qm.settings.fileStore = os.path.join(os.getenv('RMGpy'), 'QMfiles')
		qm.settings.scratchDirectory = None
		qm.settings.onlyCyclics = False
		qm.settings.maxRadicalNumber = 0
		
		if not os.path.exists(qm.settings.fileStore):
			os.mkdir(qm.settings.fileStore)
		
		mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		self.qmmol1 = MopacMolPM3(mol1, qm.settings)
		
	def testGenerateQMData(self):
		"""
		Test that generateQMData() works correctly.
		"""
		try:
			# Remove the output file so we can test the mopac script
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
		have tested the running of MOPAC to generate a `.thermo` file, so this will now test
		the loading reading of an existing `.thermo` file.
		"""
		self.qmmol1.generateThermoData()
		self.assertIsNotNone(self.qmmol1.thermo)
		self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 1) # to 1 decimal place
		self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1) # to 1 decimal place
		

		
################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )