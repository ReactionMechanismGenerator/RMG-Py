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
This module provides a class for deriving and applying two types of
bond additivity corrections (BACs).

The first type, Petersson-type BACs, are described in:
Petersson et al., J. Chem. Phys. 1998, 109, 10570-10579

The second type, Melius-type BACs, are described in:
Anantharaman and Melius, J. Phys. Chem. A 2005, 109, 1734-1747
"""

import logging
import re
from typing import Dict, Iterable, Union

import numpy as np
import pybel

from rmgpy.molecule import Atom, Bond, Molecule, get_element

import arkane.encorr.data as data
from arkane.exceptions import BondAdditivityCorrectionError


class BAC:
    """
    A class for deriving and applying bond additivity corrections.
    """

    atom_spins = {
        'H': 0.5, 'C': 1.0, 'N': 1.5, 'O': 1.0, 'F': 0.5,
        'Si': 1.0, 'P': 1.5, 'S': 1.0, 'Cl': 0.5, 'Br': 0.5, 'I': 0.5
    }
    exp_coeff = 3.0  # Melius-type parameter (Angstrom^-1)

    def __init__(self, model_chemistry: str, bac_type: str = 'p'):
        """
        Initialize a BAC instance.

        There are two implemented BAC types:
            Petersson-type: Petersson et al., J. Chem. Phys. 1998, 109, 10570-10579
            Melius-type: Anantharaman and Melius, J. Phys. Chem. A 2005, 109, 1734-1747

        Args:
            model_chemistry: Model chemistry to get preexisting BACs or data from reference database.
            bac_type: Type of BACs to get/fit ('p' for Petersson and 'm' for Melius).
        """
        self._model_chemistry = self._bac_type = None  # Set these first to avoid errors in setters
        self.model_chemistry = model_chemistry
        self.bac_type = bac_type

    @property
    def bac_type(self) -> str:
        return self._bac_type

    @bac_type.setter
    def bac_type(self, val: str):
        """Check validity and update BACs every time the BAC type is changed."""
        if val not in {'m', 'p'}:
            raise BondAdditivityCorrectionError(f'Invalid BAC type: {val}')
        self._bac_type = val
        self._update_bacs()

    @property
    def model_chemistry(self) -> str:
        return self._model_chemistry

    @model_chemistry.setter
    def model_chemistry(self, val: str):
        """Update BACs every time the model chemistry is changed."""
        self._model_chemistry = val
        self._update_bacs()

    def _update_bacs(self):
        self.bacs = None
        try:
            if self.bac_type == 'm':
                self.bacs = data.mbac[self.model_chemistry]
            elif self.bac_type == 'p':
                self.bacs = data.pbac[self.model_chemistry]
        except KeyError:
            pass

    def get_correction(self,
                       bonds: Dict[str, int] = None,
                       coords: np.ndarray = None,
                       nums: Iterable[int] = None,
                       multiplicity: int = None) -> float:
        """
        Returns the bond additivity correction in J/mol.

        There are two bond additivity corrections currently supported.
        Peterson-type corrections can be specified by setting
        `self.bac_type` to 'p'. This will use the `bonds` variable,
        which is a dictionary associating bond types with the number of
        that bond in the molecule.

        The Melius-type BAC is specified with 'm' and utilizes the atom
        coordinates in `coords` and the structure's multiplicity.

        Args:
            bonds: A dictionary of bond types (e.g., 'C=O') with their associated counts.
            coords: A Numpy array of Cartesian molecular coordinates.
            nums: A sequence of atomic numbers.
            multiplicity: The spin multiplicity of the molecule.

        Returns:
            The bond correction to the electronic energy in J/mol.
        """
        if self.bacs is None:
            bac_type_str = 'Melius' if self.bac_type == 'm' else 'Petersson'
            raise BondAdditivityCorrectionError(
                f'Missing {bac_type_str}-type BAC parameters for model chemistry {self.model_chemistry}'
            )

        if self.bac_type == 'm':
            return self._get_melius_correction(coords, nums, multiplicity=multiplicity)
        elif self.bac_type == 'p':
            return self._get_petersson_correction(bonds)

    def _get_petersson_correction(self, bonds: Dict[str, int]) -> float:
        """
        Given the model_chemistry and a dictionary of bonds, return the
        total BAC.

        Args:
            bonds: Dictionary of bonds with the following format:
                bonds = {
                    'C-H': C-H_bond_count,
                    'C-C': C-C_bond_count,
                    'C=C': C=C_bond_count,
                    ...
                }

        Returns:
            Petersson-type bond additivity correction in J/mol.
        """

        # Sum up corrections for all bonds
        bac = 0.0
        for symbol, count in bonds.items():
            if symbol in self.bacs:
                bac += count * self.bacs[symbol]
            else:
                symbol_flipped = ''.join(re.findall('[a-zA-Z]+|[^a-zA-Z]+', symbol)[::-1])  # Check reversed symbol
                if symbol_flipped in self.bacs:
                    bac += count * self.bacs[symbol_flipped]
                else:
                    logging.warning(f'Bond correction not applied for unknown bond type {symbol}.')

        return bac * 4184.0  # Convert kcal/mol to J/mol

    def _get_melius_correction(self,
                               coords: np.ndarray = None,
                               nums: Iterable[int] = None,
                               multiplicity: int = None,
                               params: Dict[str, Union[float, Dict[str, float]]] = None) -> float:
        """
        Given the model chemistry, molecular coordinates, atomic numbers,
        and dictionaries of BAC parameters, return the total BAC.

        Notes:
            A molecular correction term other than 0 destroys the size
            consistency of the quantum chemistry method. This correction
            also requires the multiplicity of the molecule.

            The negative of the total correction described in
            Anantharaman and Melius (JPCA 2005) is returned so that it
            can be added to the energy.

        Args:
            coords: Numpy array of Cartesian atomic coordinates.
            nums: Sequence of atomic numbers.
            multiplicity: Multiplicity of the molecule.
            params: Optionally provide parameters other than those stored in self.

        Returns:
            Melius-type bond additivity correction in J/mol.
        """
        if params is None:
            params = self.bacs
        atom_corr = params['atom_corr']
        bond_corr_length = params['bond_corr_length']
        bond_corr_neighbor = params['bond_corr_neighbor']
        mol_corr = params.get('mol_corr', 0.0)

        # Get single-bonded RMG molecule
        mol = _geo_to_mol(nums, coords)

        # Molecular correction
        if mol_corr != 0 and multiplicity is None:
            raise BondAdditivityCorrectionError(f'Missing multiplicity for {mol}')
        spin = 0.5 * (multiplicity - 1)
        bac_mol = mol_corr * (spin - sum(self.atom_spins[atom.element.symbol] for atom in mol.atoms))

        # Atomic correction
        bac_atom = sum(atom_corr[atom.element.symbol] for atom in mol.atoms)

        # Bond correction
        bac_bond = 0.0
        for bond in mol.get_all_edges():
            atom1 = bond.atom1
            atom2 = bond.atom2
            symbol1 = atom1.element.symbol
            symbol2 = atom2.element.symbol

            # Bond length correction
            length_corr = (bond_corr_length[symbol1] * bond_corr_length[symbol2]) ** 0.5
            length = np.linalg.norm(atom1.coords - atom2.coords)
            bac_bond += length_corr * np.exp(-self.exp_coeff * length)

            # Neighbor correction
            for other_atom, other_bond in mol.get_bonds(atom1).items():  # Atoms adjacent to atom1
                if other_bond is not bond:
                    other_symbol = other_atom.element.symbol
                    bac_bond += bond_corr_neighbor[symbol1] + bond_corr_neighbor[other_symbol]
            for other_atom, other_bond in mol.get_bonds(atom2).items():  # Atoms adjacent to atom2
                if other_bond is not bond:
                    other_symbol = other_atom.element.symbol
                    bac_bond += bond_corr_neighbor[symbol2] + bond_corr_neighbor[other_symbol]

        # Note the minus sign
        return -(bac_mol + bac_atom + bac_bond) * 4184.0  # Convert kcal/mol to J/mol


def _geo_to_mol(nums: Iterable[int], coords: np.ndarray) -> Molecule:
    """
    Convert molecular geometry specified by atomic coordinates and
    atomic numbers to RMG molecule.

    Use Open Babel for most cases because it's better at recognizing
    long bonds. Use RMG for hydrogen because Open Babel can't do it for
    mysterious reasons.
    """
    if list(nums) == [1, 1]:
        mol = Molecule()
        mol.from_xyz(nums, coords)
    else:
        symbols = [get_element(int(n)).symbol for n in nums]
        xyz = f'{len(symbols)}\n\n'
        xyz += '\n'.join(f'{s}  {c[0]: .10f}  {c[1]: .10f}  {c[2]: .10f}' for s, c in zip(symbols, coords))
        mol = pybel.readstring('xyz', xyz)
        mol = _pybel_to_rmg(mol)
    return mol


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
