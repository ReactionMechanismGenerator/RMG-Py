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

"""
This script contains unit tests for the :mod:`arkane.encorr.ae` module.
"""

import importlib
import os
import tempfile


from arkane.encorr.ae import AE, SPECIES_LABELS
from arkane.modelchem import LevelOfTheory
import pytest


class TestAE:
    """
    A class for testing that the AEJob class functions properly.
    """

    @classmethod
    def setup_class(cls):
        cls.species_energies = {lbl: i + 1 for i, lbl in enumerate(SPECIES_LABELS)}
        cls.ae = AE(cls.species_energies)

    def test_load_refdata(self):
        """
        Test that the species for fitting can be loaded.
        """
        self.ae._load_refdata()
        assert self.ae.ref_data is not None
        for spc in self.ae.ref_data.values():
            spc_ref_data = spc.reference_data[self.ae.ref_data_src]
            assert spc_ref_data.atomization_energy is not None
            assert spc_ref_data.zpe is not None

        # Test that new instance already has data loaded
        ae = AE(self.species_energies)
        assert ae.ref_data is not None

    def test_fit(self):
        """
        Test that atom energies can be fitted.
        """
        assert self.ae.atom_energies is None
        assert self.ae.confidence_intervals is None

        self.ae.fit()
        assert self.ae.atom_energies is not None
        assert self.ae.confidence_intervals is not None

    def test_write_to_database(self):
        """
        Test that results can be written to the database.
        """
        # Check that error is raised when no energies are available
        self.ae.atom_energies = None
        with pytest.raises(ValueError) as e:
            self.ae.write_to_database("test")
        assert "No atom energies" in str(e.exconly())

        # Check that error is raised if energies already exist
        self.ae.atom_energies = {"H": 1.0, "C": 2.0}
        tmp_datafile_fd, tmp_datafile_path = tempfile.mkstemp(suffix=".py")

        lot = LevelOfTheory(method="wb97m-v", basis="def2-tzvpd", software="Q-Chem")
        with pytest.raises(ValueError) as e:
            self.ae.write_to_database(lot, alternate_path=tmp_datafile_path)
        assert "overwrite" in str(e.exconly())

        # Dynamically set data file as module
        spec = importlib.util.spec_from_file_location(os.path.basename(tmp_datafile_path), tmp_datafile_path)
        module = importlib.util.module_from_spec(spec)

        # Check that existing energies can be overwritten
        self.ae.write_to_database(lot, overwrite=True, alternate_path=tmp_datafile_path)
        spec.loader.exec_module(module)  # Load data as module
        assert self.ae.atom_energies == module.atom_energies[repr(lot)]

        # Check that new energies can be written
        lot = LevelOfTheory("test")
        self.ae.write_to_database(lot, alternate_path=tmp_datafile_path)
        spec.loader.exec_module(module)  # Reload data module
        assert self.ae.atom_energies == module.atom_energies[repr(lot)]

        os.close(tmp_datafile_fd)
        os.remove(tmp_datafile_path)
