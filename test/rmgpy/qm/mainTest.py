#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import os
import shutil
import subprocess
import unittest

from rmgpy import get_path
from rmgpy.molecule import Molecule
from rmgpy.qm.gaussian import Gaussian
from rmgpy.qm.main import QMSettings, QMCalculator
from rmgpy.qm.mopac import Mopac
from rmgpy.species import Species


MOPAC_CLOSE_ENOUGH_PERCENT = 0.005  # 0.5%


class TestQMSettings(unittest.TestCase):
    """
    Contains unit tests for the QMSettings class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        rmg_path = os.path.normpath(os.path.join(get_path(), ".."))

        self.settings1 = QMSettings(
            software="mopac",
            method="pm3",
            fileStore=os.path.join(rmg_path, "testing", "qm", "QMfiles"),
            scratchDirectory=None,
            onlyCyclics=False,
            maxRadicalNumber=0,
        )

        self.settings2 = QMSettings()

    def test_check_all_set(self):
        """
        Test that check_all_set() works correctly.
        """
        try:
            self.settings1.check_all_set()
        except AssertionError:
            self.fail("check_all_set() raised unexpected AssertionError.")

        with self.assertRaises(AssertionError):
            self.settings2.check_all_set()


class TestQMCalculator(unittest.TestCase):
    """
    Contains unit tests for the QMSettings class.
    """

    mopexecutable_path = Mopac.executable_path

    gaussexecutable_path = Gaussian.executable_path
    NO_GAUSSIAN = not os.path.exists(gaussexecutable_path)

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        rmg_path = os.path.normpath(os.path.join(get_path(), ".."))
        self.fileStore = os.path.join(rmg_path, "testing", "qm", "QMfiles")

        self.mop1 = QMCalculator(
            software="mopac", method="pm3", fileStore=self.fileStore
        )

        self.mop2 = QMCalculator(
            software="mopac",
            method="pm6",
        )

        self.mop3 = QMCalculator(
            software="mopac", method="pm7", fileStore=self.fileStore
        )

        self.mop4 = QMCalculator(
            software="mopac", method="pm8", fileStore=self.fileStore
        )

        self.gauss1 = QMCalculator(
            software="gaussian",
            method="pm3",
        )

        self.gauss2 = QMCalculator(
            software="gaussian", method="pm6", fileStore=self.fileStore
        )

        self.gauss3 = QMCalculator(
            software="gaussian", method="pm7", fileStore=self.fileStore
        )

        self.molpro1 = QMCalculator(
            software="molpro", method="mp2", fileStore=self.fileStore
        )

        self.qmmol1 = QMCalculator(fileStore=self.fileStore)

        self.qmmol2 = QMCalculator(fileStore=self.fileStore)

    def test_set_default_output_directory(self):
        """
        Test that set_default_output_directory() works correctly.
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
        outputDirectory = os.path.join(self.mop1.settings.fileStore, "..", "..")
        self.mop1.set_default_output_directory(outputDirectory)
        self.mop2.set_default_output_directory(outputDirectory)
        self.mop3.set_default_output_directory(outputDirectory)
        self.gauss1.set_default_output_directory(outputDirectory)
        self.gauss2.set_default_output_directory(outputDirectory)

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

    def test_initialize(self):
        """
        Test that initialize() works correctly.
        """

        # Now set the default directories for those not set
        outputDirectory = os.path.join(self.mop1.settings.fileStore, "..", "..")
        self.mop1.set_default_output_directory(outputDirectory)
        self.mop2.set_default_output_directory(outputDirectory)
        self.mop3.set_default_output_directory(outputDirectory)
        self.gauss1.set_default_output_directory(outputDirectory)
        self.gauss2.set_default_output_directory(outputDirectory)

        try:
            self.mop1.initialize()
            self.mop2.initialize()
            self.mop3.initialize()
            self.gauss1.initialize()
            self.gauss2.initialize()
        except AssertionError:
            self.fail("initialize() raised unexpected AssertionError.")
        except Exception:
            self.fail(
                "initialize() raised Exception. Output file paths not correctly set."
            )

    def test_get_thermo_data(self):
        """
        Test that get_thermo_data() fails when expected.
        """
        output_directory = os.path.join(self.mop4.settings.fileStore, "..", "..")
        self.mop4.set_default_output_directory(output_directory)
        self.gauss3.set_default_output_directory(output_directory)
        self.molpro1.set_default_output_directory(output_directory)

        mol = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")

        with self.assertRaises(Exception):
            self.mop4.get_thermo_data(mol)
            self.gauss3.get_thermo_data(mol)
            self.molpro1.get_thermo_data(mol)

    def test_get_thermo_data_mopac(self):
        """
        Test that Mocpac get_thermo_data() works correctly.
        """
        output_directory = os.path.join(self.mop1.settings.fileStore, "..", "..")
        self.mop1.set_default_output_directory(output_directory)
        self.mop2.set_default_output_directory(output_directory)
        self.mop3.set_default_output_directory(output_directory)

        mol = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")

        for directory in (
            self.mop1.settings.fileStore,
            self.mop1.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        for directory in (
            self.mop2.settings.fileStore,
            self.mop2.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        for directory in (
            self.mop3.settings.fileStore,
            self.mop3.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        thermo1 = self.mop1.get_thermo_data(mol)
        thermo2 = self.mop2.get_thermo_data(mol)
        thermo3 = self.mop3.get_thermo_data(mol)

        self.assertTrue(thermo1.comment.startswith("QM MopacMolPM3"))
        self.assertTrue(thermo2.comment.startswith("QM MopacMolPM6"))
        self.assertTrue(thermo3.comment.startswith("QM MopacMolPM7"))

        self.assertAlmostEqual(
            thermo1.H298.value_si, 169708, delta=169708 * MOPAC_CLOSE_ENOUGH_PERCENT
        )
        self.assertAlmostEqual(
            thermo1.S298.value_si, 334.500, delta=334.500 * MOPAC_CLOSE_ENOUGH_PERCENT
        )
        self.assertAlmostEqual(
            thermo2.H298.value_si, 167704, delta=167704 * MOPAC_CLOSE_ENOUGH_PERCENT
        )
        self.assertAlmostEqual(
            thermo2.S298.value_si, 338.099, delta=338.099 * MOPAC_CLOSE_ENOUGH_PERCENT
        )
        self.assertAlmostEqual(
            thermo3.H298.value_si, 166168, delta=166168 * MOPAC_CLOSE_ENOUGH_PERCENT
        )
        self.assertAlmostEqual(
            thermo3.S298.value_si, 336.333, delta=336.333 * MOPAC_CLOSE_ENOUGH_PERCENT
        )

    @unittest.skipIf(
        NO_GAUSSIAN,
        "Gaussian not found. Try resetting your environment variables if you want to use it.",
    )
    def test_get_thermo_data_gaussian(self):
        """
        Test that Gaussian get_thermo_data() works correctly.
        """
        output_directory = os.path.join(self.mop1.settings.fileStore, "..", "..")
        self.gauss1.set_default_output_directory(output_directory)
        self.gauss2.set_default_output_directory(output_directory)

        mol = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")

        for directory in (
            self.gauss1.settings.fileStore,
            self.gauss1.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        for directory in (
            self.gauss1.settings.fileStore,
            self.gauss2.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        thermo1 = self.gauss1.get_thermo_data(mol)
        thermo2 = self.gauss2.get_thermo_data(mol)

        self.assertTrue(thermo1.comment.startswith("QM GaussianMolPM3"))
        self.assertTrue(thermo2.comment.startswith("QM GaussianMolPM6"))

        self.assertAlmostEqual(
            thermo1.H298.value_si, 169908.3376, 0
        )  # to 1 decimal place
        self.assertAlmostEqual(
            thermo1.S298.value_si, 335.5438748, 0
        )  # to 1 decimal place
        self.assertAlmostEqual(
            thermo2.H298.value_si, 169326.2504, 0
        )  # to 1 decimal place
        self.assertAlmostEqual(
            thermo2.S298.value_si, 338.2696063, 0
        )  # to 1 decimal place

    def test_run_jobs(self):
        """Test that run_jobs() works properly."""
        qm = QMCalculator(
            software="mopac",
            method="pm3",
            fileStore=self.fileStore,
            onlyCyclics=True,
            maxRadicalNumber=0,
        )
        output_directory = os.path.join(qm.settings.fileStore, "..", "..")
        qm.set_default_output_directory(output_directory)

        spc1 = Species().from_smiles("c1ccccc1")
        spc2 = Species().from_smiles("CC1C=CC=CC=1")
        spc_list = [spc1, spc2]

        qm.run_jobs(spc_list, procnum=1)


################################################################################
