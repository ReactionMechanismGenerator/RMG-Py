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


import numpy as np

from rmgpy import get_path
from rmgpy.molecule.molecule import Molecule
from rmgpy.qm.gaussian import Gaussian, GaussianMolPM3, GaussianMolPM6
from rmgpy.qm.main import QMCalculator

import pytest

EXECUTABLE_PATH = Gaussian.executable_path
NO_GAUSSIAN = not os.path.exists(EXECUTABLE_PATH)

mol1 = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")


@pytest.mark.skipif(
    NO_GAUSSIAN,
    reason="Gaussian not found. Try resetting your environment variables if you want to use it.",
)
class TestGaussianMolPM3:
    """
    Contains unit tests for the Geometry class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        rmgpy_path = os.path.normpath(os.path.join(get_path(), ".."))

        qm = QMCalculator(
            software="gaussian",
            method="pm3",
            fileStore=os.path.join(rmgpy_path, "testing", "qm", "QMfiles"),
            scratchDirectory=os.path.join(rmgpy_path, "testing", "qm", "QMscratch"),
        )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        self.qmmol1 = GaussianMolPM3(mol1, qm.settings)

    def test_generate_thermo_data(self):
        """
        Test that generate_thermo_data() works correctly on gaussian PM3.
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (
            self.qmmol1.settings.fileStore,
            self.qmmol1.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data

        assert self.qmmol1.thermo.comment.startswith("QM GaussianMolPM3 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.0626), 3) == 0

    def test_load_thermo_data(self):
        """
        Test that generate_thermo_data() can load thermo from the previous gaussian PM3 run.

        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data

        assert self.qmmol1.thermo.comment.startswith("QM GaussianMolPM3 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.0626), 3) == 0


@pytest.mark.skipif(
    NO_GAUSSIAN,
    reason="Gaussian not found. Try resetting your environment variables if you want to use it.",
)
class TestGaussianMolPM6:
    """
    Contains unit tests for the Geometry class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        rmgpy_path = os.path.normpath(os.path.join(get_path(), ".."))

        qm = QMCalculator(
            software="gaussian",
            method="pm6",
            fileStore=os.path.join(rmgpy_path, "testing", "qm", "QMfiles"),
            scratchDirectory=os.path.join(rmgpy_path, "testing", "qm", "QMscratch"),
        )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        self.qmmol1 = GaussianMolPM6(mol1, qm.settings)

    @pytest.mark.skipif("g03" in EXECUTABLE_PATH, reason="This test was shown not to work on g03.")
    def test_generate_thermo_data(self):
        """
        Test that generate_thermo_data() works correctly for gaussian PM6.
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (
            self.qmmol1.settings.fileStore,
            self.qmmol1.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data

        assert self.qmmol1.thermo.comment.startswith("QM GaussianMolPM6 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.0626), 3) == 0

    @pytest.mark.skipif("g03" in EXECUTABLE_PATH, reason="This test was shown not to work on g03.")
    def test_load_thermo_data(self):
        """
        Test that generate_thermo_data() can load thermo from the previous gaussian PM6 run.

        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data

        assert self.qmmol1.thermo.comment.startswith("QM GaussianMolPM6 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.0626), 3) == 0
