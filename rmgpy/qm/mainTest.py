#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import unittest
import subprocess
import os
import shutil

from rmgpy import getPath
from rmgpy.qm.main import QMSettings, QMCalculator
from rmgpy.molecule import Molecule

from rmgpy.qm.gaussian import Gaussian
from rmgpy.qm.mopac import Mopac

class TestQMSettings(unittest.TestCase):
    """
    Contains unit tests for the QMSettings class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
        
        self.settings1 = QMSettings(software = 'mopac',
                                   method = 'pm3',
                                   fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles'),
                                   scratchDirectory = None,
                                   onlyCyclics = False,
                                   maxRadicalNumber = 0,
                                   )
        
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
	
	mopExecutablePath = Mopac.executablePath
	if not os.path.exists(mopExecutablePath):
		NO_MOPAC = NO_LICENCE = True
	else:
		NO_MOPAC = False
		process = subprocess.Popen(mopExecutablePath,
	                               stdin=subprocess.PIPE,
	                               stdout=subprocess.PIPE,
	                               stderr=subprocess.PIPE)
		stdut, stderr = process.communicate("\n")
		NO_LICENCE = 'To install the MOPAC license' in stderr

	gaussExecutablePath = Gaussian.executablePath
	NO_GAUSSIAN = not os.path.exists(gaussExecutablePath)

	
	def setUp(self):
		"""
		A function run before each unit test in this class.
		"""
		RMGpy_path = os.path.normpath(os.path.join(getPath(),'..'))
		
		fileStore = os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles')
		
		self.mop1 = QMCalculator(software = 'mopac',
								method = 'pm3',
								fileStore = fileStore
								)
				
		self.mop2 = QMCalculator(software = 'mopac',
								method = 'pm6',
								)
		
		self.mop3 = QMCalculator(software = 'mopac',
								method = 'pm7',
								fileStore = fileStore
								)
		
		self.mop4 = QMCalculator(software = 'mopac',
								method = 'pm8',
								fileStore = fileStore
								)
		
		self.gauss1 = QMCalculator(software = 'gaussian',
								  method = 'pm3',
								  )	
		
		self.gauss2 = QMCalculator(software = 'gaussian',
								  method = 'pm6',
								  fileStore = fileStore
								  )
		
		self.gauss3 = QMCalculator(software = 'gaussian',
								  method = 'pm7',
								  fileStore = fileStore
								  )
		
		self.molpro1 = QMCalculator(software = 'molpro',
								   method = 'mp2',
								   fileStore = fileStore
								   )
		
		self.qmmol1 = QMCalculator(fileStore=fileStore)
		
		self.qmmol2 = QMCalculator(fileStore=fileStore)

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
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..','..')
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
	
	def testInitialize(self):
		"""
		Test that initialize() works correctly.
		"""
		
		# Now set the default directories for those not set
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..', '..')
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
			self.fail("initialize() raised unexpected AssertionError.")
		except Exception:
			self.fail("initialize() raised Exception. Output file paths not correctly set.")
	
	def testGetThermoData(self):
		"""
		Test that getThermoData() fails when expected.
		"""
		outputDirectory = os.path.join(self.mop4.settings.fileStore, '..', '..')
		self.mop4.setDefaultOutputDirectory(outputDirectory)
		self.gauss3.setDefaultOutputDirectory(outputDirectory)
		self.molpro1.setDefaultOutputDirectory(outputDirectory)
		
		mol = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		
		with self.assertRaises(Exception):
			self.mop4.getThermoData(mol)
			self.gauss3.getThermoData(mol)
			self.molpro1.getThermoData(mol)
		
	@unittest.skipIf(NO_MOPAC, "MOPAC not found. Try resetting your environment variables if you want to use it.")
	@unittest.skipIf(NO_LICENCE, "MOPAC license not installed. Run mopac for instructions")
	def testGetThermoDataMopac(self):
		"""
		Test that Mocpac getThermoData() works correctly.
		"""
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..', '..')
		self.mop1.setDefaultOutputDirectory(outputDirectory)
		self.mop2.setDefaultOutputDirectory(outputDirectory)
		self.mop3.setDefaultOutputDirectory(outputDirectory)
		
		mol = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		
		for directory in (self.mop1.settings.fileStore, self.mop1.settings.scratchDirectory):
			shutil.rmtree(directory, ignore_errors=True)
		
		for directory in (self.mop2.settings.fileStore, self.mop2.settings.scratchDirectory):
			shutil.rmtree(directory, ignore_errors=True)
		
		for directory in (self.mop3.settings.fileStore, self.mop3.settings.scratchDirectory):
			shutil.rmtree(directory, ignore_errors=True)
			
		thermo1 = self.mop1.getThermoData(mol)
		thermo2 = self.mop2.getThermoData(mol)
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
		
	@unittest.skipIf(NO_GAUSSIAN, "Gaussian not found. Try resetting your environment variables if you want to use it.")
	def testGetThermoDataGaussian(self):
		"""
		Test that Gaussian getThermoData() works correctly.
		"""
		outputDirectory = os.path.join(self.mop1.settings.fileStore, '..', '..')
		self.gauss1.setDefaultOutputDirectory(outputDirectory)
		self.gauss2.setDefaultOutputDirectory(outputDirectory)
		
		mol = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
		
		for directory in (self.gauss1.settings.fileStore, self.gauss1.settings.scratchDirectory):
			shutil.rmtree(directory, ignore_errors=True)
		
		for directory in (self.gauss1.settings.fileStore, self.gauss2.settings.scratchDirectory):
			shutil.rmtree(directory, ignore_errors=True)
					
		thermo1 = self.gauss1.getThermoData(mol)
		thermo2 = self.gauss2.getThermoData(mol)
		
		self.assertTrue(thermo1.comment.startswith('QM GaussianMolPM3'))
		self.assertTrue(thermo2.comment.startswith('QM GaussianMolPM6'))
		
		self.assertAlmostEqual(thermo1.H298.value_si, 169908.3376, 0) # to 1 decimal place
		self.assertAlmostEqual(thermo1.S298.value_si, 335.5438748, 0) # to 1 decimal place
		self.assertAlmostEqual(thermo2.H298.value_si, 169326.2504, 0) # to 1 decimal place
		self.assertAlmostEqual(thermo2.S298.value_si, 338.2696063, 0) # to 1 decimal place

################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
