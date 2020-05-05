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
This module provides methods to decompose a molecule into
substructures/motifs. The presence of these motifs in the molecules in
the reference data can be used to weight the BAC fitting and thus fit
to maximally diversified data.
"""

from collections import Counter
from typing import List

from rdkit import Chem


def get_substructs(smi: str) -> Counter:
    """
    Convert SMILES to counts of substructures.

    Args:
        smi: SMILES.

    Returns:
        Counter containing number of times each substructure SMILES occurs.
    """
    mol = _smi_to_mol(smi)
    cliques = substruct_decomp(mol)
    return Counter(_get_substruct(mol, c) for c in cliques)


def substruct_decomp(mol: Chem.Mol) -> List[List[int]]:
    """
    Decompose molecule into substructures. The substructures are
    defined as atoms, bonds, two-bond and three-bond substructures that
    are not in rings, rings, and bridged ring systems that share at
    least three atoms.

    Note:
        Adapted from junction tree decomposition described in
        https://arxiv.org/abs/1802.04364.

    Args:
        mol: RDKit molecule.

    Returns:
        List of indices contained in each substructure.
    """
    n_atoms = mol.GetNumAtoms()
    if n_atoms == 1:
        return [[0]]

    cliques = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        if not bond.IsInRing():
            cliques.append([a1, a2])

    ssr = [list(x) for x in Chem.GetSymmSSSR(mol)]
    cliques.extend(ssr)

    nei_list = [[] for _ in range(n_atoms)]
    for i, c in enumerate(cliques):
        for atom in c:
            nei_list[atom].append(i)

    # Merge Rings with intersection > 2 atoms
    for i, c in enumerate(cliques):
        if len(c) > 2:
            for atom in c:
                for j in nei_list[atom]:
                    if i >= j or len(cliques[j]) <= 2:
                        continue
                    inter = set(c) & set(cliques[j])
                    if len(inter) > 2:
                        cliques[i].extend(cliques[j])
                        cliques[i] = list(set(cliques[i]))
                        cliques[j] = []

    cliques = [c for c in cliques if len(c) > 0]
    nei_list = [[] for _ in range(n_atoms)]
    for i, c in enumerate(cliques):
        for atom in c:
            nei_list[atom].append(i)

    # Add singleton cliques
    for atom in range(n_atoms):
        cnei = nei_list[atom]
        if len(cnei) > 1:
            bonds = [c for c in cnei if len(cliques[c]) == 2]
            if len(bonds) > 1 or (len(bonds) <= 1 and len(cnei) > 2):
                cliques.append([atom])

    return cliques


def _get_substruct(mol: Chem.Mol, atoms: List[int]) -> str:
    """Convert a list of atom indices to a substructure."""
    if mol.GetNumAtoms() == 1:
        smiles = _mol_to_smi(mol)
    else:
        # For single-atom cliques, we want the substructure to contain its neighbors
        if len(atoms) == 1:
            atoms = atoms[:]
            atoms.extend([nei.GetIdx() for nei in mol.GetAtomWithIdx(atoms[0]).GetNeighbors()])
        smiles = Chem.MolFragmentToSmiles(mol, atoms, kekuleSmiles=True)
    return _mol_to_smi(_copy_mol(Chem.MolFromSmiles(smiles, sanitize=False)))


def _mol_to_smi(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, kekuleSmiles=True)


def _smi_to_mol(smi: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol)
    return mol


def _copy_mol(mol: Chem.Mol) -> Chem.Mol:
    new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = _copy_atom(atom)
        new_mol.AddAtom(new_atom)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)
    return new_mol.GetMol()


def _copy_atom(atom: Chem.Atom) -> Chem.Atom:
    new_atom = Chem.Atom(atom.GetSymbol())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons())
    return new_atom
