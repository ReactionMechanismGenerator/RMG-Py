#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
This script contains unit tests for the :mod:`arkane.encorr.bac` module.
"""

import importlib
import os
import tempfile
import unittest
from collections import Counter

import numpy as np
from openbabel import pybel

from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.quantity import ScalarQuantity

from arkane.encorr.bac import BAC
from arkane.encorr.data import BACDataset, BOND_SYMBOLS, _pybel_to_rmg
from arkane.encorr.reference import ReferenceDatabase
from arkane.exceptions import BondAdditivityCorrectionError
from arkane.modelchem import LevelOfTheory


class TestBAC(unittest.TestCase):
    """
    A class for testing that the BAC class functions properly.
    """

    @classmethod
    def setUpClass(cls):
        cls.lot_get = LevelOfTheory(method='CCSD(T)-F12', basis='cc-pVTZ-F12', software='Molpro')
        cls.lot_fit = LevelOfTheory(method='wB97M-V', basis='def2-TZVPD', software='Q-Chem')
        cls.lot_nonexisting = LevelOfTheory('notamethod')

        cls.bac = BAC(cls.lot_get)

        cls.tmp_melius_params = {
            'atom_corr': {'H': 1.0, 'C': 2.0, 'N': 3.0, 'O': 4.0, 'S': 5.0, 'F': 6.0, 'Cl': 7.0, 'Br': 8.0},
            'bond_corr_length': {'H': 1.0, 'C': 2.0, 'N': 3.0, 'O': 4.0, 'S': 5.0, 'F': 6.0, 'Cl': 7.0, 'Br': 8.0},
            'bond_corr_neighbor': {'H': 1.0, 'C': 2.0, 'N': 3.0, 'O': 4.0, 'S': 5.0, 'F': 6.0, 'Cl': 7.0, 'Br': 8.0},
            'mol_corr': 1.0
        }
        cls.tmp_petersson_params = {'C-H': 1.0, 'C-C': 2.0, 'C=C': 3.0, 'C-O': 4.0}

        # Set molecule, bonds, nums, and coords for testing Petersson and Melius BACs
        cls.multiplicity = 1
        smi = 'C=C(OSC=S)C#CC1C(=O)N=CNC1SSC(O)C#N'

        mol = Molecule(smiles=smi, multiplicity=cls.multiplicity)
        cls.bonds = Counter(
            f'{b.atom1.element.symbol}{BOND_SYMBOLS[b.order]}{b.atom2.element.symbol}' for b in mol.get_all_edges()
        )

        pybel_mol = pybel.readstring('smi', smi)
        pybel_mol.addh()
        pybel_mol.make3D()
        mol_3d = _pybel_to_rmg(pybel_mol)
        cls.nums = [atom.number for atom in mol_3d.atoms]
        cls.coords = np.array([atom.coords for atom in mol_3d.atoms])

    def test_loading_parameters(self):
        """
        Test that BAC parameters for levels of theory are loaded
        correctly and that errors are raised otherwise.
        """
        self.bac.level_of_theory = self.lot_get
        self.bac.bac_type = 'p'
        self.assertIsInstance(self.bac.bacs, dict)

        self.bac.level_of_theory = self.lot_nonexisting
        self.bac.bac_type = 'm'
        self.assertIsNone(self.bac.bacs)

        with self.assertRaises(BondAdditivityCorrectionError):
            self.bac.bac_type = ''

    def test_load_database(self):
        """
        Test that reference database can be loaded.
        """
        key = self.bac.load_database(names='main')
        expected_key = (os.path.join(settings['database.directory'], 'reference_sets', 'main'),)
        self.assertEqual(key, expected_key)
        self.assertIsInstance(self.bac.ref_databases[key], ReferenceDatabase)

        # Test that other instance already has loaded database
        bac = BAC(self.lot_fit)
        self.assertIsInstance(bac.ref_databases[key], ReferenceDatabase)

    def test_get_correction(self):
        """
        Test that BAC corrections can be obtained.
        """
        self.bac.level_of_theory = self.lot_get
        self.bac.bac_type = 'p'
        corr = self.bac.get_correction(bonds=self.bonds)
        self.assertIsInstance(corr, ScalarQuantity)

        # Can use actual Melius parameters once they're available in database
        self.bac.bac_type = 'm'
        with self.assertRaises(BondAdditivityCorrectionError):
            # No multiplicity specified
            self.bac._get_melius_correction(coords=self.coords, nums=self.nums, params=self.tmp_melius_params)
        corr1 = self.bac._get_melius_correction(coords=self.coords, nums=self.nums, multiplicity=self.multiplicity,
                                                params=self.tmp_melius_params)
        self.assertIsInstance(corr1, ScalarQuantity)

        self.bac.level_of_theory = self.lot_nonexisting
        with self.assertRaises(BondAdditivityCorrectionError):
            self.bac.get_correction()

    def _clear_bac_data(self):
        self.bac.bacs = None
        self.bac.species = self.bac.ref_data = self.bac.calc_data = self.bac.bac_data = None

    def _check_bac_data(self):
        self.assertIsInstance(self.bac.bacs, dict)
        self.assertIsInstance(self.bac.dataset, BACDataset)
        self.assertIsNotNone(self.bac.database_key)
        self.assertLess(self.bac.dataset.calculate_stats(for_bac_data=True).rmse,
                        self.bac.dataset.calculate_stats().rmse)

    def test_fit_petersson(self):
        """
        Test that Petersson BAC parameters can be derived.
        """
        self.bac.level_of_theory = self.lot_fit
        self.bac.bac_type = 'p'
        self._clear_bac_data()
        self.bac.fit()

        self._check_bac_data()
        self.assertIn('C-H', self.bac.bacs)

        self.bac.level_of_theory = self.lot_nonexisting
        with self.assertRaises(BondAdditivityCorrectionError):
            self.bac.fit()

    def test_fit_melius(self):
        """
        Test that Melius BAC parameters can be derived.
        """
        self.bac.level_of_theory = self.lot_fit
        self.bac.bac_type = 'm'
        self._clear_bac_data()

        # With molecular correction, no global opt
        self.bac.fit(fit_mol_corr=True, global_opt=False, lsq_max_nfev=50)
        self._check_bac_data()
        self.assertSetEqual(set(self.bac.bacs.keys()),
                            {'atom_corr', 'bond_corr_length', 'bond_corr_neighbor', 'mol_corr'})
        self.assertNotAlmostEqual(self.bac.bacs['mol_corr'], 0.0)

        # Without molecular correction, with global opt
        self._clear_bac_data()
        self.bac.fit(fit_mol_corr=False, global_opt=True, global_opt_iter=1, lsq_max_nfev=50)
        self._check_bac_data()
        self.assertAlmostEqual(self.bac.bacs['mol_corr'], 0.0)

    def test_test(self):
        """
        Test that enthalpies of formation can be evaluated.
        """
        with self.assertRaises(BondAdditivityCorrectionError) as e:
            self.bac.test(species=[], db_names=[])
        self.assertIn('several data sources', str(e.exception))

        with self.assertRaises(BondAdditivityCorrectionError) as e:
            self.bac.test(species=[])
        self.assertIn('No data', str(e.exception))

        self.bac.level_of_theory = self.lot_fit
        self.bac.bac_type = 'm'
        self.bac.bacs = self.tmp_melius_params

        # Get a few species to test on
        key = self.bac.load_database(names='main')
        species = self.bac.ref_databases[key].extract_level_of_theory(self.bac.level_of_theory,
                                                                      as_error_canceling_species=False)[:10]

        dataset = self.bac.test(species=species)
        self.assertIsInstance(dataset, BACDataset)
        self.assertIsNotNone(dataset.bac_data)

    def test_write_to_database(self):
        """
        Test that BAC parameters can be written to a file.
        """
        # Check that error is raised when no BACs are available
        self.bac.bacs = None
        with self.assertRaises(BondAdditivityCorrectionError) as e:
            self.bac.write_to_database()
        self.assertIn('No BACs', str(e.exception))

        self.bac.level_of_theory = self.lot_get
        self.bac.bac_type = 'p'
        self.bac.bacs = self.tmp_petersson_params

        tmp_datafile_fd, tmp_datafile_path = tempfile.mkstemp(suffix='.py')

        # Check that error is raised if BACs already exist and overwrite is False
        with self.assertRaises(IOError) as e:
            self.bac.write_to_database(alternate_path=tmp_datafile_path)
        self.assertIn('overwrite', str(e.exception))

        # Dynamically set data file as module
        spec = importlib.util.spec_from_file_location(os.path.basename(tmp_datafile_path), tmp_datafile_path)
        module = importlib.util.module_from_spec(spec)

        # Check that existing Petersson BACs can be overwritten
        self.bac.write_to_database(overwrite=True, alternate_path=tmp_datafile_path)
        spec.loader.exec_module(module)  # Load data as module
        self.assertEqual(self.bac.bacs, module.pbac[repr(self.bac.level_of_theory)])

        # Check that new Petersson BACs can be written
        self.bac.level_of_theory = self.lot_nonexisting
        self.bac.bacs = self.tmp_petersson_params
        self.bac.write_to_database(alternate_path=tmp_datafile_path)
        spec.loader.exec_module(module)  # Reload data module
        self.assertEqual(self.bac.bacs, module.pbac[repr(self.bac.level_of_theory)])

        # Check that new Melius BACs can be written
        self.bac.bac_type = 'm'
        self.bac.bacs = self.tmp_melius_params
        self.bac.write_to_database(alternate_path=tmp_datafile_path)
        spec.loader.exec_module(module)
        self.assertEqual(self.bac.bacs, module.mbac[repr(self.bac.level_of_theory)])

        os.close(tmp_datafile_fd)
        os.remove(tmp_datafile_path)

    def test_save_correlation_mat(self):
        """
        Test that visual of correlation matrix can be generated.
        """
        self.bac.correlation = None
        with self.assertRaises(BondAdditivityCorrectionError):
            self.bac.save_correlation_mat('')

        self.bac.bacs = self.tmp_melius_params
        self.bac.correlation = np.random.uniform(size=(24, 24))

        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_corr_path = os.path.join(tmpdirname, 'corr.pdf')
            self.bac.save_correlation_mat(tmp_corr_path)
            self.assertTrue(os.path.exists(tmp_corr_path))


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
