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

import numpy as np
import os
import shutil
import subprocess

from rmgpy import getPath
from rmgpy.qm.main import QMCalculator
from rmgpy.molecule import Molecule
from rmgpy.qm.mopac import Mopac, MopacMolPM3, MopacMolPM6, MopacMolPM7
from rmgpy.exceptions import DependencyError


NO_MOPAC = NO_LICENCE = False
try:
    Mopac().testReady()
except DependencyError, e:
    if "Couldn't find MOPAC executable" in e.message:
        NO_MOPAC = NO_LICENCE = True
    elif 'To install the MOPAC license' in e.message or 'MOPAC_LICENSE' in e.message:
        NO_LICENCE = True
    else:
        raise

mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')

class TestMopacMolPM3(unittest.TestCase):
    """
    Contains unit tests for the Geometry class.
    """

    @unittest.skipIf(NO_MOPAC, "MOPAC not found. Try resetting your environment variables if you want to use it.")
    @unittest.skipIf(NO_LICENCE, "MOPAC license not installed. Run mopac for instructions")
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        RMGpy_path = os.path.normpath(os.path.join(getPath(), '..'))

        qm = QMCalculator(software='mopac',
                          method='pm3',
                          fileStore=os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles'),
                          scratchDirectory=os.path.join(RMGpy_path, 'testing', 'qm', 'QMscratch'),
                          )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        self.qmmol1 = MopacMolPM3(mol1, qm.settings)

    def testGenerateThermoData(self):
        """
        Test that generateThermoData() works correctly for MOPAC PM3
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (self.qmmol1.settings.fileStore, self.qmmol1.settings.scratchDirectory):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData
        self.assertTrue(self.qmmol1.verifyOutputFile())
        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM3 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == 'amu':
            self.assertAlmostEqual(result.molecularMass.value, 128.173, 2)

        self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 0)  # to 1 decimal place
        self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1)  # to 1 decimal place

    def testLoadThermoData(self):
        """
        Test that generateThermoData() can load thermo from the previous MOPAC PM3 run.
        
        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData

        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM3 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == 'amu':
            self.assertAlmostEqual(result.molecularMass.value, 128.173, 2)

        self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 169708.0608, 0)  # to 1 decimal place
        self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 334.5007584, 1)  # to 1 decimal place

class TestMopacMolPM6(unittest.TestCase):
    """
    Contains unit tests for the Geometry class.
    """

    @unittest.skipIf(NO_MOPAC, "MOPAC not found. Try resetting your environment variables if you want to use it.")
    @unittest.skipIf(NO_LICENCE, "MOPAC license not installed. Run mopac for instructions")
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        RMGpy_path = os.path.normpath(os.path.join(getPath(), '..'))

        qm = QMCalculator(software='mopac',
                          method='pm6',
                          fileStore=os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles'),
                          scratchDirectory=os.path.join(RMGpy_path, 'testing', 'qm', 'QMscratch'),
                          )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        self.qmmol1 = MopacMolPM6(mol1, qm.settings)

    def testGenerateThermoData(self):
        """
        Test that generateThermoData() works correctly for MOPAC PM6
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (self.qmmol1.settings.fileStore, self.qmmol1.settings.scratchDirectory):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData
        self.assertTrue(self.qmmol1.verifyOutputFile())
        
        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM6 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == 'amu':
            self.assertAlmostEqual(result.molecularMass.value, 128.173, 2)

        self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 167704.4270, 0)  # to 1 decimal place
        self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 338.0999241, 1)  # to 1 decimal place

    def testLoadThermoData(self):
        """
        Test that generateThermoData() can load thermo from the previous MOPAC PM6 run.
        
        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData

        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM6 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == 'amu':
            self.assertEqual(result.molecularMass.value, 128.173)

        self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 167704.0681, 0)  # to 0 decimal place
        self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 338.0999241, 1)  # to 1 decimal place

class TestMopacMolPM7(unittest.TestCase):
    """
    Contains unit tests for the Geometry class.
    """

    @unittest.skipIf(NO_MOPAC, "MOPAC not found. Try resetting your environment variables if you want to use it.")
    @unittest.skipIf(NO_LICENCE, "MOPAC license not installed. Run mopac for instructions")
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        RMGpy_path = os.path.normpath(os.path.join(getPath(), '..'))

        qm = QMCalculator(software='mopac',
                          method='pm7',
                          fileStore=os.path.join(RMGpy_path, 'testing', 'qm', 'QMfiles'),
                          scratchDirectory=os.path.join(RMGpy_path, 'testing', 'qm', 'QMscratch'),
                          )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
        self.qmmol1 = MopacMolPM7(mol1, qm.settings)

    def testGenerateThermoData(self):
        """
        Test that generateThermoData() works correctly for MOPAC PM7
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (self.qmmol1.settings.fileStore, self.qmmol1.settings.scratchDirectory):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData
        self.assertTrue(self.qmmol1.verifyOutputFile())
        
        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM7 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == 'amu':
            self.assertAlmostEqual(result.molecularMass.value, 128.173, 2)

        self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 166168.9863, 0)  # to 1 decimal place
        self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 336.3330406, 1)  # to 1 decimal place
        
    def testLoadThermoData(self):
        """
        Test that generateThermoData() can load thermo from the previous MOPAC PM7 run.
        
        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData

        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM MopacMolPM7 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == 'amu':
            self.assertAlmostEqual(result.molecularMass.value, 128.173, 2)

        self.assertAlmostEqual(self.qmmol1.thermo.H298.value_si, 166168.8571, 0)  # to 1 decimal place
        self.assertAlmostEqual(self.qmmol1.thermo.S298.value_si, 336.3330406, 1)  # to 1 decimal place
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
