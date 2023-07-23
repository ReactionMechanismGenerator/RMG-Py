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
from rmgpy.exceptions import DependencyError
from rmgpy.molecule.molecule import Molecule
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.mopac import Mopac, MopacMolPM3, MopacMolPM6, MopacMolPM7

mol1 = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")
MOPAC_CLOSE_ENOUGH_PERCENT = 0.005  # 0.5%


class TestMopacMolPM3:
    """
    Contains unit tests for the Geometry class.
    """

    def setup_class(self):
        rmgpy_path = os.path.normpath(os.path.join(get_path(), ".."))

        qm = QMCalculator(
            software="mopac",
            method="pm3",
            fileStore=os.path.join(rmgpy_path, "testing", "qm", "QMfiles"),
            scratchDirectory=os.path.join(rmgpy_path, "testing", "qm", "QMscratch"),
        )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        self.qmmol1 = MopacMolPM3(mol1, qm.settings)

    def test_generate_thermo_data(self):
        """
        Test that generate_thermo_data() works correctly for MOPAC PM3
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (
            self.qmmol1.settings.fileStore,
            self.qmmol1.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data
        assert self.qmmol1.verify_output_file()
        assert self.qmmol1.thermo.comment.startswith("QM MopacMolPM3 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.173), 2) == 0

        assert abs(self.qmmol1.thermo.H298.value_si - 169708) < 169708 * MOPAC_CLOSE_ENOUGH_PERCENT
        assert abs(self.qmmol1.thermo.S298.value_si - 334.500) < 334.500 * MOPAC_CLOSE_ENOUGH_PERCENT

    def test_load_thermo_data(self):
        """
        Test that generate_thermo_data() can load thermo from the previous MOPAC PM3 run.

        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data

        assert self.qmmol1.thermo.comment.startswith("QM MopacMolPM3 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.173), 2) == 0

        assert abs(self.qmmol1.thermo.H298.value_si - 169708) < 169708 * MOPAC_CLOSE_ENOUGH_PERCENT
        assert abs(self.qmmol1.thermo.S298.value_si - 334.500) < 334.500 * MOPAC_CLOSE_ENOUGH_PERCENT


class TestMopacMolPM6:
    """
    Contains unit tests for the Geometry class.
    """

    def setup_class(self):
        rmgpy_path = os.path.normpath(os.path.join(get_path(), ".."))

        qm = QMCalculator(
            software="mopac",
            method="pm6",
            fileStore=os.path.join(rmgpy_path, "testing", "qm", "QMfiles"),
            scratchDirectory=os.path.join(rmgpy_path, "testing", "qm", "QMscratch"),
        )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        self.qmmol1 = MopacMolPM6(mol1, qm.settings)

    def test_generate_thermo_data(self):
        """
        Test that generate_thermo_data() works correctly for MOPAC PM6
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (
            self.qmmol1.settings.fileStore,
            self.qmmol1.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data
        assert self.qmmol1.verify_output_file()

        assert self.qmmol1.thermo.comment.startswith("QM MopacMolPM6 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.173), 2) == 0

        assert abs(self.qmmol1.thermo.H298.value_si - 167704) < 167704 * MOPAC_CLOSE_ENOUGH_PERCENT
        assert abs(self.qmmol1.thermo.S298.value_si - 338.099) < 338.099 * MOPAC_CLOSE_ENOUGH_PERCENT

    def test_load_thermo_data(self):
        """
        Test that generate_thermo_data() can load thermo from the previous MOPAC PM6 run.

        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data

        assert self.qmmol1.thermo.comment.startswith("QM MopacMolPM6 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert result.molecularMass.value == 128.173

        assert abs(self.qmmol1.thermo.H298.value_si - 167704) < 167704 * MOPAC_CLOSE_ENOUGH_PERCENT  # to 0 decimal place
        assert abs(self.qmmol1.thermo.S298.value_si - 338.099) < 338.099 * MOPAC_CLOSE_ENOUGH_PERCENT


class TestMopacMolPM7:
    """
    Contains unit tests for the Geometry class.
    """

    def setup_class(self):
        rmgpy_path = os.path.normpath(os.path.join(get_path(), ".."))

        qm = QMCalculator(
            software="mopac",
            method="pm7",
            fileStore=os.path.join(rmgpy_path, "testing", "qm", "QMfiles"),
            scratchDirectory=os.path.join(rmgpy_path, "testing", "qm", "QMscratch"),
        )

        if not os.path.exists(qm.settings.fileStore):
            os.makedirs(qm.settings.fileStore)

        mol1 = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")
        self.qmmol1 = MopacMolPM7(mol1, qm.settings)

    def test_generate_thermo_data(self):
        """
        Test that generate_thermo_data() works correctly for MOPAC PM7
        """
        # First ensure any old data are removed, or else they'll be reused!
        for directory in (
            self.qmmol1.settings.fileStore,
            self.qmmol1.settings.scratchDirectory,
        ):
            shutil.rmtree(directory, ignore_errors=True)

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data
        assert self.qmmol1.verify_output_file()

        assert self.qmmol1.thermo.comment.startswith("QM MopacMolPM7 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.173), 2) == 0

        assert abs(self.qmmol1.thermo.H298.value_si - 166169) < 166169 * MOPAC_CLOSE_ENOUGH_PERCENT
        assert abs(self.qmmol1.thermo.S298.value_si - 336.333) < 336.333 * MOPAC_CLOSE_ENOUGH_PERCENT

    def test_load_thermo_data(self):
        """
        Test that generate_thermo_data() can load thermo from the previous MOPAC PM7 run.

        Check that it loaded, and the values are the same as above.
        """

        self.qmmol1.generate_thermo_data()
        result = self.qmmol1.qm_data

        assert self.qmmol1.thermo.comment.startswith("QM MopacMolPM7 calculation")
        assert result.numberOfAtoms == 18
        assert isinstance(result.atomicNumbers, np.ndarray)
        if result.molecularMass.units == "amu":
            assert round(abs(result.molecularMass.value - 128.173), 2) == 0

        assert abs(self.qmmol1.thermo.H298.value_si - 166169) < 166169 * MOPAC_CLOSE_ENOUGH_PERCENT
        assert abs(self.qmmol1.thermo.S298.value_si - 336.333) < 336.333 * MOPAC_CLOSE_ENOUGH_PERCENT
