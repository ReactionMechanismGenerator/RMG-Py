#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This file loads in the atomic energies, frequency scale factors, and BAC parameters for several model chemistries. The
data is stored in the RMG-database in input/quantum_corrections/data.py

This module also provides classes to store data for BAC fitting and evaluating.
"""

import functools
import importlib.util
import os
from collections import Counter
from dataclasses import dataclass
from typing import Callable, Iterable, List, Sequence, Set, Union

import numpy as np
from openbabel import pybel

from rmgpy import settings
from rmgpy.molecule import Atom, Bond, get_element
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.rmgobject import recursive_make_object

from arkane.encorr.decomp import get_substructs
from arkane.encorr.reference import ReferenceSpecies, ReferenceDatabase
from arkane.exceptions import BondAdditivityCorrectionError
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory

# ######## Database loading ##########
quantum_corrections_path = os.path.join(settings['database.directory'], 'quantum_corrections', 'data.py')
spec = importlib.util.spec_from_file_location("quantum_calculations", quantum_corrections_path)
data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data)

# Convert module to dictionary and create level of theory objects from their string representations
data_dict = recursive_make_object(
    {k: v for k, v in vars(data).items() if not k.startswith('__')},
    {'LevelOfTheory': LevelOfTheory, 'CompositeLevelOfTheory': CompositeLevelOfTheory}
)

# Assign the data here so that it can be imported
for k, v in data_dict.items():
    vars()[k] = v
######################################


# ########## Data classes ############
BOND_SYMBOLS = {1: '-', 2: '=', 3: '#'}


class Molecule(RMGMolecule):
    """Wrapper for RMG Molecule to add ID attribute"""
    def __init__(self, *args, mol_id=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.id = mol_id


@dataclass
class Stats:
    """Small class to store BAC fitting statistics"""
    rmse: Union[float, np.ndarray]
    mae: Union[float, np.ndarray]


class BACDatapoint:
    """A BACDatapoint contains a single species"""

    class _Decorators:
        @staticmethod
        def assert_level_of_theory(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                if args[0].level_of_theory is None:  # args[0] is the instance
                    raise BondAdditivityCorrectionError('Level of theory is not defined')
                return func(*args, **kwargs)
            return wrapper

    def __init__(self, spc: ReferenceSpecies, level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory] = None):
        self.spc = spc
        self.level_of_theory = level_of_theory

        self._mol = None
        self._mol_type = None
        self._bonds = None
        self._ref_data = None
        self._calc_data = None
        self._bac_data = None
        self._substructs = None
        self.weight = 1

    @property
    def mol(self) -> Molecule:
        if self._mol is None:
            raise ValueError('Call `BACDatapoint.to_mol` first')
        return self._mol

    def to_mol(self, from_geo: bool = False) -> Molecule:
        """
        Convert to RMG molecule. If `from_geo` is True, a single-bonded
        molecule is perceived from 3D coordinates with Open Babel.
        """
        if self._mol is None:
            if from_geo:
                self._mol_from_geo()
            else:
                self._mol_from_adjlist()
        else:  # Use cached molecule if possible
            if from_geo and self._mol_type != 'geo':
                self._mol_from_geo()
            elif not from_geo and self._mol_type != 'adj':
                self._mol_from_adjlist()

        return self._mol

    @_Decorators.assert_level_of_theory
    def _mol_from_geo(self):
        xyz = self.spc.calculated_data[self.level_of_theory].xyz_dict
        self._mol = geo_to_mol(xyz['coords'], symbols=xyz['symbols'])
        self._mol_type = 'geo'

    def _mol_from_adjlist(self):
        self._mol = Molecule().from_adjacency_list(self.spc.adjacency_list,
                                                   raise_atomtype_exception=False,
                                                   raise_charge_exception=False)
        self._mol_type = 'adj'

    @property
    def bonds(self) -> Counter:
        """Get bond counts"""
        if self._bonds is None:
            mol = self.to_mol(from_geo=False)
            bonds = Counter()
            for bond in mol.get_all_edges():
                symbols = [bond.atom1.element.symbol, bond.atom2.element.symbol]
                symbols.sort()
                symbol = symbols[0] + BOND_SYMBOLS[bond.order] + symbols[1]
                bonds[symbol] += 1
            self._bonds = bonds
        return self._bonds

    @property
    def ref_data(self) -> float:
        """Get reference enthalpy in kcal/mol"""
        if self._ref_data is None:
            self._ref_data = self.spc.get_reference_enthalpy().h298.value_si / 4184
        return self._ref_data

    @property
    @_Decorators.assert_level_of_theory
    def calc_data(self) -> float:
        """Get calculated enthalpy in kcal/mol"""
        if self._calc_data is None:
            self._calc_data = self.spc.calculated_data[self.level_of_theory].thermo_data.H298.value_si / 4184
        return self._calc_data

    @property
    def bac_data(self) -> float:
        if self._bac_data is None:
            raise ValueError('No BAC data available')
        return self._bac_data

    @bac_data.setter
    def bac_data(self, val: float):
        self._bac_data = val

    @property
    def substructs(self) -> Counter:
        """Decompose into substructures"""
        if self._substructs is None:
            self._substructs = get_substructs(self.spc.smiles)

            # Add charge and multiplicity "substructures"
            if self.spc.charge == 0:
                self._substructs['neutral'] = 1
            elif self.spc.charge > 0:
                self._substructs['cation'] = 1
            elif self.spc.charge < 0:
                self._substructs['anion'] = 1

            if self.spc.multiplicity == 1:
                self._substructs['singlet'] = 1
            elif self.spc.multiplicity == 2:
                self._substructs['doublet'] = 1
            elif self.spc.multiplicity >= 3:
                self._substructs['triplet+'] = 1

        return self._substructs


class DatasetProperty:
    """
    Descriptor to simplify BACDataset properties.

    This descriptor is essentially a specialized version of Python
    properties. An instance of this descriptor can be defined as a
    class attribute in a class containing a `data` attribute, where
    each item in the `data` sequence contains an attribute with its
    name corresponding to the first argument passed to `__init__`.

    The descriptor is used by accessing the attribute with the <attr>
    name, which implicitly retrieves the value cached in the <_attr>
    attribute. If no cached value exists, the <attr> attributes from
    all items in `data` are obtained as a list or as a Numpy array (if
    `asarray` is `True`), cached in <_attr> (only if <attr> in each
    `data` item cannot be changed, i.e., `settable` should be `False`),
    and returned. If all <attr>'s in `data` are `None`, `None` is
    returned (instead of a sequence of `None`s).

    The descriptor can also be used to set the <attr> attributes in the
    `data` items. This requires that `settable` is `True`, otherwise an
    `AttributeError` is raised. It also requires that the length of the
    sequence used for setting the <attr> attributes is the same length
    as `data`.
    """

    def __init__(self, attr, asarray=False, settable=False):
        self.pub_attr = attr  # Name of attribute defined in BACDataset and available in BACDatapoint
        self.priv_attr = '_' + attr  # Private attribute that is set in BACDataset and used to retrieve cached values
        self.asarray = asarray  # Whether the data should be get/set as Numpy arrays
        self.settable = settable  # Whether the BACDatapoint attributes can be set

    def __get__(self, obj, objtype=None):
        if hasattr(obj, self.priv_attr):  # Return cached value if available
            return getattr(obj, self.priv_attr)
        val = [getattr(d, self.pub_attr) for d in obj.data]  # Retrieve the attributes from the items in data
        if all(v is None for v in val):  # Instead of returning list of None's, just return None
            val = None
        elif self.asarray:
            val = np.array(val)
        if not self.settable:  # Only cache sequence if it cannot be changed in BACDatapoint
            setattr(obj, self.priv_attr, val)
        return val

    def __set__(self, obj, val):
        if not self.settable:
            raise AttributeError(f'Cannot set {self.pub_attr}')
        if len(val) != len(obj):  # Requires that __len__ is defined in obj
            raise ValueError('Number of data do not match number of datapoints')
        for d, v in zip(obj.data, val):
            setattr(d, self.pub_attr, v)


class BACDataset:
    """A BACDataset contains a list of BACDatapoints"""
    def __init__(self, data: List[BACDatapoint]):
        self.data = data

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, item) -> Union[BACDatapoint, List[BACDatapoint]]:
        return self.data[item]

    def append(self, item: BACDatapoint):
        self.data.append(item)

    def sort(self, key: Callable, reverse: bool = False):
        self.data.sort(key=key, reverse=reverse)

    bonds = DatasetProperty('bonds')
    ref_data = DatasetProperty('ref_data', asarray=True)
    calc_data = DatasetProperty('calc_data', asarray=True)
    bac_data = DatasetProperty('bac_data', asarray=True, settable=True)
    substructs = DatasetProperty('substructs')
    weight = DatasetProperty('weight', asarray=True, settable=True)
    weights = weight  # Alias for weight

    def get_mols(self, from_geo: bool = False) -> List[Molecule]:
        return [d.to_mol(from_geo=from_geo) for d in self.data]

    def calculate_stats(self, for_bac_data: bool = False) -> Stats:
        """
        Calculate RMSE and MAE with respect to `calc_data` or `bac_data`.

        Args:
            for_bac_data: Use the fitted BAC data to calculate RMSE/MAE if set to True.
                          Otherwise, calculate RMSE/MAE with respect to the calculated (quantum) data.

        Returns:
            Stats dataclass containing RMSE and MAE.
        """
        other_data = self.bac_data if for_bac_data else self.calc_data
        diff = self.ref_data - other_data
        rmse = np.sqrt(np.dot(diff, diff) / len(self.ref_data))
        mae = np.sum(np.abs(diff)) / len(self.ref_data)
        return Stats(rmse, mae)

    def compute_weights(self, weight_type: str = 'substructs'):
        """
        Set weights for each datapoint such that molecular diversity is
        maximized. I.e., effectively balance dataset by having higher
        weights for molecules with underrepresented substructures.

        Args:
            weight_type: Currently only supports 'substructs'.
        """
        if weight_type == 'substructs':
            # Counts of substructures across all molecules
            all_substructs = Counter()
            for s in self.substructs:
                all_substructs += s

            # Compute weight for each molecule as average across substructure frequencies
            self.weights = [
                sum(1 / all_substructs[s] for s in substructs.elements())  # Sum of frequencies
                / sum(substructs.values())  # Divide by number of substructures in molecule
                for substructs in self.substructs  # For each molecule
            ]
        else:
            raise NotImplementedError(f'{weight_type} weight type is unavailable')


def extract_dataset(ref_database: ReferenceDatabase,
                    level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory],
                    exclude_elements: Union[Sequence[str], Set[str], str] = None,
                    charge: Union[Sequence[Union[str, int]], Set[Union[str, int]], str, int] = 'all',
                    multiplicity: Union[Sequence[int], Set[int], int, str] = 'all') -> BACDataset:
    """
    Extract species for a given model chemistry from a reference
    database and convert to a BACDataset.

    Args:
         ref_database: Reference database.
         level_of_theory: Level of theory.
         exclude_elements: Sequence of element symbols to exclude.
         charge: Allowable charges. Possible values are 'all'; a combination of 'neutral, 'positive', and 'negative';
                 or a sequence of integers.
         multiplicity: Allowable multiplicites. Possible values are 'all' or positive integers.

    Returns:
        BACDataset containing species with data available at given level of theory.
    """
    species = ref_database.extract_level_of_theory(level_of_theory, as_error_canceling_species=False)

    if exclude_elements is not None:
        elements = {exclude_elements} if isinstance(exclude_elements, str) else set(exclude_elements)
        species = [spc for spc in species if not any(e in spc.formula for e in elements)]
    if charge != 'all':
        charges = {charge} if isinstance(charge, (str, int)) else set(charge)
        species = [spc for spc in species
                   if spc.charge == 0 and 'neutral' in charges
                   or spc.charge > 0 and 'positive' in charges
                   or spc.charge < 0 and 'negative' in charges
                   or spc.charge in charges]
    if multiplicity != 'all':
        multiplicities = {multiplicity} if isinstance(multiplicity, int) else set(multiplicity)
        species = [spc for spc in species if spc.multiplicity in multiplicities]

    return BACDataset([BACDatapoint(spc, level_of_theory=level_of_theory) for spc in species])


def geo_to_mol(coords: np.ndarray, symbols: Iterable[str] = None, nums: Iterable[int] = None) -> Molecule:
    """
    Convert molecular geometry specified by atomic coordinates and
    atomic symbols/numbers to RMG molecule.

    Use Open Babel because it's better at recognizing long bonds.
    """
    if nums is None and symbols is None:
        raise ValueError('Must specify nums or symbols')

    symbols = [get_element(int(n)).symbol for n in nums] if symbols is None else list(symbols)
    xyz = f'{len(symbols)}\n\n'
    xyz += '\n'.join(f'{s}  {c[0]: .10f}  {c[1]: .10f}  {c[2]: .10f}' for s, c in zip(symbols, coords))
    mol = pybel.readstring('xyz', xyz)
    return _pybel_to_rmg(mol)


def _pybel_to_rmg(pybel_mol: pybel.Molecule) -> Molecule:
    """
    Convert Pybel molecule to RMG molecule but ignore charge,
    multiplicity, and bond orders.
    """
    mol = Molecule()
    for pybel_atom in pybel_mol:
        element = get_element(pybel_atom.atomicnum)
        atom = Atom(element=element, coords=np.array(pybel_atom.coords))
        mol.vertices.append(atom)
    for obbond in pybel.ob.OBMolBondIter(pybel_mol.OBMol):
        begin_idx = obbond.GetBeginAtomIdx() - 1  # Open Babel indexes atoms starting at 1
        end_idx = obbond.GetEndAtomIdx() - 1
        bond = Bond(mol.vertices[begin_idx], mol.vertices[end_idx])
        mol.add_bond(bond)
    return mol
