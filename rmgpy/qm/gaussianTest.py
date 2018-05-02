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

import distutils.spawn
import itertools
import logging
import numpy as np
import os
import shutil

from rmgpy import getPath
from rmgpy.qm.main import QMCalculator
from rmgpy.molecule import Molecule
from rmgpy.qm.gaussian import Gaussian, GaussianMolPM3, GaussianMolPM6


executablePath = Gaussian.executablePath
NO_GAUSSIAN = not os.path.exists(executablePath)
    
mol1 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')

class TestGaussianMolPM3(unittest.TestCase):
    """
    Contains unit tests for the Geometry class.
    """

    @unittest.skipIf(NO_GAUSSIAN, "Gaussian not found. Try resetting your environment variables if you want to use it.")
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
        Test that generateThermoData() works correctly on gaussian PM3.
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (self.qmmol1.settings.fileStore, self.qmmol1.settings.scratchDirectory):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData

        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM GaussianMolPM3 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units=='amu':
            self.assertAlmostEqual(result.molecularMass.value, 128.0626, 3)

    def testLoadThermoData(self):
        """
        Test that generateThermoData() can load thermo from the previous gaussian PM3 run.

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

    @unittest.skipIf(NO_GAUSSIAN, "Gaussian not found. Try resetting your environment variables if you want to use it.")
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

    @unittest.skipIf('g03' in executablePath, "This test was shown not to work on g03.")
    def testGenerateThermoData(self):
        """
        Test that generateThermoData() works correctly for gaussian PM6.
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (self.qmmol1.settings.fileStore, self.qmmol1.settings.scratchDirectory):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generateThermoData()
        result = self.qmmol1.qmData

        self.assertTrue(self.qmmol1.thermo.comment.startswith('QM GaussianMolPM6 calculation'))
        self.assertEqual(result.numberOfAtoms, 18)
        self.assertIsInstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units=='amu':
            self.assertAlmostEqual(result.molecularMass.value, 128.0626, 3)

    @unittest.skipIf('g03' in executablePath, "This test was shown not to work on g03.")
    def testLoadThermoData(self):
        """
        Test that generateThermoData() can load thermo from the previous gaussian PM6 run.

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
