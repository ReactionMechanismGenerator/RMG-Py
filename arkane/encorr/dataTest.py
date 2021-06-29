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
This script contains unit tests for the :mod:`arkane.encorr.data` module.
"""

import unittest
from collections import Counter

import numpy as np
from openbabel import pybel

from rmgpy.molecule import Molecule as RMGMolecule

import arkane.encorr.data as data
from arkane.encorr.data import (Molecule, Stats, BACDatapoint, DatasetProperty, BACDataset,
                                extract_dataset, geo_to_mol, _pybel_to_rmg)
from arkane.encorr.reference import ReferenceDatabase
from arkane.exceptions import BondAdditivityCorrectionError
from arkane.modelchem import LOT, LevelOfTheory

DATABASE = ReferenceDatabase()
DATABASE.load()
LEVEL_OF_THEORY = LevelOfTheory(method='wb97m-v', basis='def2-tzvpd', software='qchem')


class TestDataLoading(unittest.TestCase):
    """
    A class for testing that the quantum correction data is loaded
    correctly from the RMG database.
    """

    def test_contains_data(self):
        """
        Test that the necessary dictionaries are available.
        """
        self.assertTrue(hasattr(data, 'atom_hf'))
        self.assertTrue(hasattr(data, 'atom_thermal'))
        self.assertTrue(hasattr(data, 'SOC'))
        self.assertTrue(hasattr(data, 'atom_energies'))
        self.assertTrue(hasattr(data, 'pbac'))
        self.assertTrue(hasattr(data, 'mbac'))
        self.assertTrue(hasattr(data, 'freq_dict'))

    def test_level_of_theory(self):
        """
        Test that level of theory objects were created.
        """
        for lot in data.atom_energies.keys():
            self.assertIsInstance(lot, LOT)
        for lot in data.pbac.keys():
            self.assertIsInstance(lot, LOT)
        for lot in data.mbac.keys():
            self.assertIsInstance(lot, LOT)
        for lot in data.freq_dict.keys():
            self.assertIsInstance(lot, LOT)


class TestMolecule(unittest.TestCase):
    """
    A class for testing that the Molecule wrapper class functions
    properly.
    """

    def test_molecule(self):
        """
        Test that a Molecule contains the `id` attribute.
        """
        rmg_mol = RMGMolecule(smiles='C')
        mol = Molecule(smiles='C')
        self.assertIsInstance(mol, RMGMolecule)
        self.assertTrue(hasattr(mol, 'id'))
        self.assertFalse(hasattr(rmg_mol, 'id'))


class TestStats(unittest.TestCase):
    """
    A class for testing that the Stats class functions properly.
    """

    def test_stats(self):
        """
        Test that a Stats instance contains the correct attributes.
        """
        stats = Stats(1.0, 2.0)
        self.assertTrue(hasattr(stats, 'rmse'))
        self.assertTrue(hasattr(stats, 'mae'))


class TestBACDatapoint(unittest.TestCase):
    """
    A class for testing that the BACDatapoint class functions properly.
    """

    @classmethod
    def setUpClass(cls):
        cls.spc = list(DATABASE.reference_sets.values())[0][0]

    def setUp(self):
        self.datapoint = BACDatapoint(self.spc, level_of_theory=LEVEL_OF_THEORY)

    def test_assert_level_of_theory(self):
        """
        Test that decorator correctly determines when a level of theory
        is not defined.
        """
        self.datapoint.level_of_theory = None
        with self.assertRaises(BondAdditivityCorrectionError):
            _ = self.datapoint.calc_data

    def test_weight(self):
        """
        Test that weight is initialized to 1.
        """
        self.assertEqual(self.datapoint.weight, 1)

    def test_mol(self):
        """
        Test that BACDatapoint can be converted to a Molecule.
        """
        with self.assertRaises(ValueError):
            _ = self.datapoint.mol

        # From adjacency list
        mol_adj = self.datapoint.to_mol(from_geo=False)
        self.assertIsInstance(mol_adj, Molecule)
        self.assertIs(mol_adj, self.datapoint.mol)
        mol_adj2 = self.datapoint.to_mol(from_geo=False)
        self.assertIs(mol_adj, mol_adj2)  # Check that cached molecule is used

        # From geometry
        mol_geo = self.datapoint.to_mol(from_geo=True)
        self.assertIsNot(mol_geo, mol_adj)  # Check that cached molecule is NOT used
        coords_spc = np.vstack(tuple(a.coords for a in mol_geo.atoms))
        coords_dp = self.spc.calculated_data[LEVEL_OF_THEORY].xyz_dict['coords']
        self.assertIsNone(np.testing.assert_allclose(coords_dp, coords_spc))
        self.assertIsInstance(mol_geo, Molecule)
        self.assertIs(mol_geo, self.datapoint.mol)
        mol_geo2 = self.datapoint.to_mol(from_geo=True)
        self.assertIs(mol_geo, mol_geo2)  # Check that cached molecule is used

    def test_bonds(self):
        """
        Test that bonds can be obtained.
        """
        bonds = self.datapoint.bonds
        self.assertIsInstance(bonds, Counter)
        bonds2 = self.datapoint.bonds
        self.assertIs(bonds, bonds2)  # Check that cached bonds are used

    def test_ref_data(self):
        """
        Test that reference data can be obtained.
        """
        ref_data = self.datapoint.ref_data
        self.assertIsInstance(ref_data, float)

    def test_calc_data(self):
        """
        Test that calculated data can be obtained.
        """
        calc_data = self.datapoint.calc_data
        self.assertIsInstance(calc_data, float)

    def test_bac_data(self):
        """
        Test that `bac_data` can be used.
        """
        with self.assertRaises(ValueError):
            _ = self.datapoint.bac_data

        self.datapoint.bac_data = 1.0
        self.assertIsInstance(self.datapoint.bac_data, float)

    def test_substructs(self):
        """
        Test that BACDatapoint can be decomposed into substructures.
        """
        substructs = self.datapoint.substructs
        self.assertIsInstance(substructs, Counter)

        # Check that exactly one of 'neutral', 'cation', or 'anion' is set
        # and same for 'singlet', 'doublet', 'triplet+'.
        self.assertEqual(sum(substructs[k] for k in ('neutral', 'cation', 'anion')), 1)  # Can only be one of these
        self.assertEqual(sum(substructs[k] for k in ('singlet', 'doublet', 'triplet+')), 1)

        substructs2 = self.datapoint.substructs
        self.assertIs(substructs, substructs2)  # Check that cached substructures are used


class TestDatasetProperty(unittest.TestCase):
    """
    A class for testing that the DatasetProperty descriptor functions
    properly.
    """

    @staticmethod
    def make_set(asarray=False, settable=False):
        class Point:
            def __init__(self, val):
                self.val = val
                self.val2 = None

        class Set:
            def __init__(self):
                self.data = [Point(i) for i in range(10)]

            def __len__(self):
                return len(self.data)

            val = DatasetProperty('val', asarray=asarray, settable=settable)
            val2 = DatasetProperty('val2', asarray=asarray, settable=settable)

        return Set()

    def test_init(self):
        """
        Test that the descriptor is initialized properly.
        """
        s = self.make_set()
        dset_prop = type(s).__dict__['val']
        self.assertIsInstance(dset_prop, DatasetProperty)
        self.assertEqual(dset_prop.pub_attr, 'val')
        self.assertEqual(dset_prop.priv_attr, '_val')
        self.assertFalse(dset_prop.asarray)
        self.assertFalse(dset_prop.settable)

    def test_get(self):
        """
        Test that the descriptor can be used to get attributes of the
        data contained in the set.
        """
        s = self.make_set(asarray=False)
        self.assertListEqual(s.val, list(range(10)))
        self.assertListEqual(s._val, list(range(10)))  # Check that list was cached
        self.assertIsNone(s.val2)

        s = self.make_set(asarray=True)
        self.assertIsInstance(s.val, np.ndarray)

    def test_set(self):
        """
        Test that the descriptor can be used to set attributes of the
        data contained in the set.
        """
        s = self.make_set(settable=False)
        with self.assertRaises(AttributeError):
            s.val = list(range(9, -1, -1))

        s = self.make_set(asarray=True, settable=True)
        with self.assertRaises(ValueError):  # Try setting with wrong length data
            s.val = list(range(9))
        s.val = list(range(9, -1, -1))
        self.assertIsInstance(s.val, np.ndarray)
        with self.assertRaises(AttributeError):  # Check that cache is not available
            _ = s._val
        self.assertTrue(all(d.val == v for d, v in zip(s.data, list(range(9, -1, -1)))))


class TestBACDataset(unittest.TestCase):
    """
    A class for testing that the BACDataset class functions properly.
    """

    @classmethod
    def setUpClass(cls):
        cls.species = list(DATABASE.reference_sets.values())[0][:5]

    def setUp(self):
        self.dataset = BACDataset([BACDatapoint(spc, level_of_theory=LEVEL_OF_THEORY) for spc in self.species])

    def test_append(self):
        """
        Test that a datapoint can be appended.
        """
        self.dataset.append(BACDatapoint(self.species[0]))
        self.assertEqual(len(self.dataset), len(self.species) + 1)

    def test_sort(self):
        """
        Test that the dataset can be sorted.
        """
        self.dataset.sort(key=lambda d: d.spc.smiles)  # Sort by SMILES
        smiles = sorted(spc.smiles for spc in self.species)
        smiles_from_dset = [d.spc.smiles for d in self.dataset]
        self.assertListEqual(smiles, smiles_from_dset)

    def test_attrs(self):
        """
        Test that DatasetProperty attributes behave properly.
        """
        self.assertIsInstance(self.dataset.bonds, list)
        self.assertEqual(len(self.dataset.bonds), len(self.species))
        self.assertIsInstance(self.dataset.ref_data, np.ndarray)
        self.assertEqual(len(self.dataset.ref_data), len(self.species))
        self.assertIsInstance(self.dataset.calc_data, np.ndarray)
        self.assertEqual(len(self.dataset.calc_data), len(self.species))
        self.assertIsInstance(self.dataset.substructs, list)
        self.assertEqual(len(self.dataset.substructs), len(self.species))

        with self.assertRaises(ValueError):
            _ = self.dataset.bac_data
        self.dataset.bac_data = list(range(len(self.species)))
        self.assertIsInstance(self.dataset.bac_data, np.ndarray)
        self.assertEqual(len(self.dataset.bac_data), len(self.species))

        self.assertIsInstance(self.dataset.weight, np.ndarray)
        self.assertTrue(all(w1 == w2 for w1, w2 in zip(self.dataset.weight, self.dataset.weights)))
        self.assertEqual(len(self.dataset.weights), len(self.species))

    def test_get_mols(self):
        """
        Test that molecules can be retrieved.
        """
        mols = self.dataset.get_mols()
        self.assertIsInstance(mols, list)
        self.assertEqual(len(mols), len(self.species))

    def test_calculate_stats(self):
        """
        Test that RMSE and MAE are calculated correctly.
        """
        stats_calc = self.dataset.calculate_stats()
        self.assertLessEqual(stats_calc.mae, stats_calc.rmse)

        with self.assertRaises(ValueError):
            _ = self.dataset.calculate_stats(for_bac_data=True)

        self.dataset.bac_data = list(range(len(self.species)))
        stats_bac = self.dataset.calculate_stats(for_bac_data=True)
        self.assertLessEqual(stats_bac.mae, stats_bac.rmse)
        self.assertNotEqual(stats_calc, stats_bac)

    def test_compute_weights(self):
        """
        Test that weights can be computed.
        """
        with self.assertRaises(NotImplementedError):
            self.dataset.compute_weights(weight_type='')

        self.dataset.compute_weights()
        self.assertTrue(all(0 <= w < 1 for w in self.dataset.weights))


class TestFuncs(unittest.TestCase):
    """
    A class for testing that the functions in the data module work.
    """

    def test_extract_dataset(self):
        """
        Test that a reference dataset can be extracted.
        """
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY)
        self.assertIsInstance(dataset, BACDataset)

        # Test excluding elements
        elements = 'N'
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, exclude_elements=elements)
        for d in dataset:
            self.assertFalse(elements in d.spc.formula)
        elements = ['N', 'O']
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, exclude_elements=elements)
        for d in dataset:
            self.assertFalse(any(e in d.spc.formula for e in elements))

        # Test specifying charge
        charges = 'negative'
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            self.assertTrue(d.spc.charge < 0)
        charges = 1
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            self.assertTrue(d.spc.charge == 1)
        charges = ['positive', 'neutral']
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            self.assertTrue(d.spc.charge >= 0)
        charges = [-1, 'positive']
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            self.assertTrue(d.spc.charge > 0 or d.spc.charge == -1)

        # Test specifying multiplicity
        multiplicities = 2
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, multiplicity=multiplicities)
        for d in dataset:
            self.assertTrue(d.spc.multiplicity == 2)
        multiplicities = [2, 3]
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, multiplicity=multiplicities)
        for d in dataset:
            self.assertTrue(d.spc.multiplicity in {2, 3})

    def test_geo_to_mol(self):
        """
        Test that a geometry can be converted to an RMG molecule.
        """
        # Hydrogen
        nums = (1, 1)
        coords = np.array([[0, 0, 0], [0, 0, 0.74]])
        mol = geo_to_mol(coords, nums=nums)
        self.assertIsInstance(mol, Molecule)
        self.assertTrue(mol.to_single_bonds().is_isomorphic(Molecule(smiles='[H][H]').to_single_bonds()))

        # Methane
        symbols = ('H', 'C', 'H', 'H', 'H')
        coords = np.array([
            [0.5288, 0.1610, 0.9359],
            [0.0000, 0.0000, 0.0000],
            [0.2051, 0.8240, -0.6786],
            [0.3345, -0.9314, -0.4496],
            [-1.0685, -0.0537, 0.1921]
        ])
        mol = geo_to_mol(coords, symbols=symbols)
        self.assertIsInstance(mol, Molecule)
        self.assertTrue(mol.to_single_bonds().is_isomorphic(Molecule(smiles='C').to_single_bonds()))

    def test_pybel_to_rmg(self):
        """
        Test that a Pybel molecule can be converted to an RMG molecule.
        """
        pybel_mol = pybel.readstring('smi', 'O=CCN')
        pybel_mol.addh()
        pybel_mol.make3D()
        mol = _pybel_to_rmg(pybel_mol)

        for atom, pybel_atom in zip(mol.atoms, pybel_mol.atoms):
            self.assertEqual(atom.number, pybel_atom.atomicnum)
            self.assertIsNone(np.testing.assert_allclose(atom.coords, pybel_atom.coords))
            self.assertEqual(Counter(a.number for a in atom.bonds),
                             Counter(a.GetAtomicNum() for a in pybel.ob.OBAtomAtomIter(pybel_atom.OBAtom)))
