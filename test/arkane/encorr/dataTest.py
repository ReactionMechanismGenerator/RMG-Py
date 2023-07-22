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
This script contains unit tests for the :mod:`arkane.encorr.data` module.
"""


from collections import Counter

import numpy as np
from openbabel import pybel

from rmgpy.molecule import Molecule as RMGMolecule

import arkane.encorr.data as data
from arkane.encorr.data import (
    Molecule,
    Stats,
    BACDatapoint,
    DatasetProperty,
    BACDataset,
    extract_dataset,
    geo_to_mol,
    _pybel_to_rmg,
)
from arkane.encorr.reference import ReferenceDatabase
from arkane.exceptions import BondAdditivityCorrectionError
from arkane.modelchem import LOT, LevelOfTheory
import pytest

DATABASE = ReferenceDatabase()
DATABASE.load()
LEVEL_OF_THEORY = LevelOfTheory(method="wb97m-v", basis="def2-tzvpd", software="qchem")


class TestDataLoading:
    """
    A class for testing that the quantum correction data is loaded
    correctly from the RMG database.
    """

    def test_contains_data(self):
        """
        Test that the necessary dictionaries are available.
        """
        assert hasattr(data, "atom_hf")
        assert hasattr(data, "atom_thermal")
        assert hasattr(data, "SOC")
        assert hasattr(data, "atom_energies")
        assert hasattr(data, "pbac")
        assert hasattr(data, "mbac")
        assert hasattr(data, "freq_dict")

    def test_level_of_theory(self):
        """
        Test that level of theory objects were created.
        """
        for lot in data.atom_energies.keys():
            assert isinstance(lot, LOT)
        for lot in data.pbac.keys():
            assert isinstance(lot, LOT)
        for lot in data.mbac.keys():
            assert isinstance(lot, LOT)
        for lot in data.freq_dict.keys():
            assert isinstance(lot, LOT)


class TestMolecule:
    """
    A class for testing that the Molecule wrapper class functions
    properly.
    """

    def test_molecule(self):
        """
        Test that a Molecule contains the `id` attribute.
        """
        rmg_mol = RMGMolecule(smiles="C")
        mol = Molecule(smiles="C")
        assert isinstance(mol, RMGMolecule)
        assert hasattr(mol, "id")
        assert not hasattr(rmg_mol, "id")


class TestStats:
    """
    A class for testing that the Stats class functions properly.
    """

    def test_stats(self):
        """
        Test that a Stats instance contains the correct attributes.
        """
        stats = Stats(1.0, 2.0)
        assert hasattr(stats, "rmse")
        assert hasattr(stats, "mae")


class TestBACDatapoint:
    """
    A class for testing that the BACDatapoint class functions properly.
    """

    @classmethod
    def setup_class(cls):
        cls.spc = list(DATABASE.reference_sets.values())[0][0]

    def setup_class(self):
        self.datapoint = BACDatapoint(self.spc, level_of_theory=LEVEL_OF_THEORY)

    def test_assert_level_of_theory(self):
        """
        Test that decorator correctly determines when a level of theory
        is not defined.
        """
        self.datapoint.level_of_theory = None
        with pytest.raises(BondAdditivityCorrectionError):
            _ = self.datapoint.calc_data

    def test_weight(self):
        """
        Test that weight is initialized to 1.
        """
        assert self.datapoint.weight == 1

    def test_mol(self):
        """
        Test that BACDatapoint can be converted to a Molecule.
        """
        with pytest.raises(ValueError):
            _ = self.datapoint.mol

        # From adjacency list
        mol_adj = self.datapoint.to_mol(from_geo=False)
        assert isinstance(mol_adj, Molecule)
        assert mol_adj is self.datapoint.mol
        mol_adj2 = self.datapoint.to_mol(from_geo=False)
        assert mol_adj is mol_adj2  # Check that cached molecule is used

        # From geometry
        mol_geo = self.datapoint.to_mol(from_geo=True)
        assert mol_geo is not mol_adj  # Check that cached molecule is NOT used
        coords_spc = np.vstack(tuple(a.coords for a in mol_geo.atoms))
        coords_dp = self.spc.calculated_data[LEVEL_OF_THEORY].xyz_dict["coords"]
        assert np.testing.assert_allclose(coords_dp, coords_spc) is None
        assert isinstance(mol_geo, Molecule)
        assert mol_geo is self.datapoint.mol
        mol_geo2 = self.datapoint.to_mol(from_geo=True)
        assert mol_geo is mol_geo2  # Check that cached molecule is used

    def test_bonds(self):
        """
        Test that bonds can be obtained.
        """
        bonds = self.datapoint.bonds
        assert isinstance(bonds, Counter)
        bonds2 = self.datapoint.bonds
        assert bonds is bonds2  # Check that cached bonds are used

    def test_ref_data(self):
        """
        Test that reference data can be obtained.
        """
        ref_data = self.datapoint.ref_data
        assert isinstance(ref_data, float)

    def test_calc_data(self):
        """
        Test that calculated data can be obtained.
        """
        calc_data = self.datapoint.calc_data
        assert isinstance(calc_data, float)

    def test_bac_data(self):
        """
        Test that `bac_data` can be used.
        """
        with pytest.raises(ValueError):
            _ = self.datapoint.bac_data

        self.datapoint.bac_data = 1.0
        assert isinstance(self.datapoint.bac_data, float)

    def test_substructs(self):
        """
        Test that BACDatapoint can be decomposed into substructures.
        """
        substructs = self.datapoint.substructs
        assert isinstance(substructs, Counter)

        # Check that exactly one of 'neutral', 'cation', or 'anion' is set
        # and same for 'singlet', 'doublet', 'triplet+'.
        assert sum(substructs[k] for k in ("neutral", "cation", "anion")) == 1  # Can only be one of these
        assert sum(substructs[k] for k in ("singlet", "doublet", "triplet+")) == 1

        substructs2 = self.datapoint.substructs
        assert substructs is substructs2  # Check that cached substructures are used


class TestDatasetProperty:
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

            val = DatasetProperty("val", asarray=asarray, settable=settable)
            val2 = DatasetProperty("val2", asarray=asarray, settable=settable)

        return Set()

    def test_init(self):
        """
        Test that the descriptor is initialized properly.
        """
        s = self.make_set()
        dset_prop = type(s).__dict__["val"]
        assert isinstance(dset_prop, DatasetProperty)
        assert dset_prop.pub_attr == "val"
        assert dset_prop.priv_attr == "_val"
        assert not dset_prop.asarray
        assert not dset_prop.settable

    def test_get(self):
        """
        Test that the descriptor can be used to get attributes of the
        data contained in the set.
        """
        s = self.make_set(asarray=False)
        assert s.val == list(range(10))
        assert s._val == list(range(10))  # Check that list was cached
        assert s.val2 is None

        s = self.make_set(asarray=True)
        assert isinstance(s.val, np.ndarray)

    def test_set(self):
        """
        Test that the descriptor can be used to set attributes of the
        data contained in the set.
        """
        s = self.make_set(settable=False)
        with pytest.raises(AttributeError):
            s.val = list(range(9, -1, -1))

        s = self.make_set(asarray=True, settable=True)
        with pytest.raises(ValueError):  # Try setting with wrong length data
            s.val = list(range(9))
        s.val = list(range(9, -1, -1))
        assert isinstance(s.val, np.ndarray)
        with pytest.raises(AttributeError):  # Check that cache is not available
            _ = s._val
        assert all(d.val == v for d, v in zip(s.data, list(range(9, -1, -1))))


class TestBACDataset:
    """
    A class for testing that the BACDataset class functions properly.
    """

    @classmethod
    def setup_class(cls):
        cls.species = list(DATABASE.reference_sets.values())[0][:5]

    def setup_class(self):
        self.dataset = BACDataset([BACDatapoint(spc, level_of_theory=LEVEL_OF_THEORY) for spc in self.species])

    def test_append(self):
        """
        Test that a datapoint can be appended.
        """
        self.dataset.append(BACDatapoint(self.species[0]))
        assert len(self.dataset) == len(self.species) + 1

    def test_sort(self):
        """
        Test that the dataset can be sorted.
        """
        self.dataset.sort(key=lambda d: d.spc.smiles)  # Sort by SMILES
        smiles = sorted(spc.smiles for spc in self.species)
        smiles_from_dset = [d.spc.smiles for d in self.dataset]
        assert smiles == smiles_from_dset

    def test_attrs(self):
        """
        Test that DatasetProperty attributes behave properly.
        """
        assert isinstance(self.dataset.bonds, list)
        assert len(self.dataset.bonds) == len(self.species)
        assert isinstance(self.dataset.ref_data, np.ndarray)
        assert len(self.dataset.ref_data) == len(self.species)
        assert isinstance(self.dataset.calc_data, np.ndarray)
        assert len(self.dataset.calc_data) == len(self.species)
        assert isinstance(self.dataset.substructs, list)
        assert len(self.dataset.substructs) == len(self.species)

        with pytest.raises(ValueError):
            _ = self.dataset.bac_data
        self.dataset.bac_data = list(range(len(self.species)))
        assert isinstance(self.dataset.bac_data, np.ndarray)
        assert len(self.dataset.bac_data) == len(self.species)

        assert isinstance(self.dataset.weight, np.ndarray)
        assert all(w1 == w2 for w1, w2 in zip(self.dataset.weight, self.dataset.weights))
        assert len(self.dataset.weights) == len(self.species)

    def test_get_mols(self):
        """
        Test that molecules can be retrieved.
        """
        mols = self.dataset.get_mols()
        assert isinstance(mols, list)
        assert len(mols) == len(self.species)

    def test_calculate_stats(self):
        """
        Test that RMSE and MAE are calculated correctly.
        """
        stats_calc = self.dataset.calculate_stats()
        assert stats_calc.mae <= stats_calc.rmse

        with pytest.raises(ValueError):
            _ = self.dataset.calculate_stats(for_bac_data=True)

        self.dataset.bac_data = list(range(len(self.species)))
        stats_bac = self.dataset.calculate_stats(for_bac_data=True)
        assert stats_bac.mae <= stats_bac.rmse
        assert stats_calc != stats_bac

    def test_compute_weights(self):
        """
        Test that weights can be computed.
        """
        with pytest.raises(NotImplementedError):
            self.dataset.compute_weights(weight_type="")

        self.dataset.compute_weights()
        assert all(0 <= w < 1 for w in self.dataset.weights)


class TestFuncs:
    """
    A class for testing that the functions in the data module work.
    """

    def test_extract_dataset(self):
        """
        Test that a reference dataset can be extracted.
        """
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY)
        assert isinstance(dataset, BACDataset)

        # Test only retrieving specific indices
        idxs = 211
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, idxs=idxs)
        assert len(dataset) == 1
        assert dataset[0].spc.index == idxs
        idxs = [211, 362]
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, idxs=idxs)
        assert len(dataset) == 2
        for d in dataset:
            assert d.spc.index in {211, 362}

        # Test excluding indices
        idxs = 211
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, exclude_idxs=idxs)
        for d in dataset:
            assert d.spc.index != idxs
        idxs = [211, 362]
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, exclude_idxs=idxs)
        for d in dataset:
            assert d.spc.index not in {211, 362}

        # Test excluding elements
        elements = "N"
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, exclude_elements=elements)
        for d in dataset:
            assert not (elements in d.spc.formula)
        elements = ["N", "O"]
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, exclude_elements=elements)
        for d in dataset:
            assert not any(e in d.spc.formula for e in elements)

        # Test specifying charge
        charges = "negative"
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            assert d.spc.charge < 0
        charges = 1
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            assert d.spc.charge == 1
        charges = ["positive", "neutral"]
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            assert d.spc.charge >= 0
        charges = [-1, "positive"]
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, charge=charges)
        for d in dataset:
            assert d.spc.charge > 0 or d.spc.charge == -1

        # Test specifying multiplicity
        multiplicities = 2
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, multiplicity=multiplicities)
        for d in dataset:
            assert d.spc.multiplicity == 2
        multiplicities = [2, 3]
        dataset = extract_dataset(DATABASE, LEVEL_OF_THEORY, multiplicity=multiplicities)
        for d in dataset:
            assert d.spc.multiplicity in {2, 3}

    def test_geo_to_mol(self):
        """
        Test that a geometry can be converted to an RMG molecule.
        """
        # Hydrogen
        nums = (1, 1)
        coords = np.array([[0, 0, 0], [0, 0, 0.74]])
        mol = geo_to_mol(coords, nums=nums)
        assert isinstance(mol, Molecule)
        assert mol.to_single_bonds().is_isomorphic(Molecule(smiles="[H][H]").to_single_bonds())

        # Methane
        symbols = ("H", "C", "H", "H", "H")
        coords = np.array(
            [
                [0.5288, 0.1610, 0.9359],
                [0.0000, 0.0000, 0.0000],
                [0.2051, 0.8240, -0.6786],
                [0.3345, -0.9314, -0.4496],
                [-1.0685, -0.0537, 0.1921],
            ]
        )
        mol = geo_to_mol(coords, symbols=symbols)
        assert isinstance(mol, Molecule)
        assert mol.to_single_bonds().is_isomorphic(Molecule(smiles="C").to_single_bonds())

    def test_pybel_to_rmg(self):
        """
        Test that a Pybel molecule can be converted to an RMG molecule.
        """
        pybel_mol = pybel.readstring("smi", "O=CCN")
        pybel_mol.addh()
        pybel_mol.make3D()
        mol = _pybel_to_rmg(pybel_mol)

        for atom, pybel_atom in zip(mol.atoms, pybel_mol.atoms):
            assert atom.number == pybel_atom.atomicnum
            assert np.testing.assert_allclose(atom.coords, pybel_atom.coords) is None
            assert Counter(a.number for a in atom.bonds) == Counter(a.GetAtomicNum() for a in pybel.ob.OBAtomAtomIter(pybel_atom.OBAtom))
