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
This module provides methods for converting molecules between RMG, RDKit, and OpenBabel.
"""

import logging
import sys

import cython
# Assume that rdkit is installed
from rdkit import Chem
# Test if openbabel is installed
try:
    from openbabel import openbabel
except ImportError:
    openbabel = None

import rmgpy.molecule.element as elements
import rmgpy.molecule.molecule as mm
from rmgpy.exceptions import DependencyError


def to_rdkit_mol(mol, remove_h=True, return_mapping=False, sanitize=True):
    """
    Convert a molecular structure to a RDKit rdmol object. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    Perceives aromaticity and, unless remove_h==False, removes Hydrogen atoms.

    If return_mapping==True then it also returns a dictionary mapping the
    atoms to RDKit's atom indices.
    """

    # Sort the atoms before converting to ensure output is consistent
    # between different runs
    mol.sort_atoms()
    atoms = mol.vertices
    rd_atom_indices = {}  # dictionary of RDKit atom indices
    rdkitmol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    for index, atom in enumerate(mol.vertices):
        if atom.element.symbol == 'X':
            rd_atom = Chem.rdchem.Atom('Pt')  # not sure how to do this with linear scaling when this might not be Pt
        else:
            rd_atom = Chem.rdchem.Atom(atom.element.symbol)
        if atom.element.isotope != -1:
            rd_atom.SetIsotope(atom.element.isotope)
        rd_atom.SetNumRadicalElectrons(atom.radical_electrons)
        rd_atom.SetFormalCharge(atom.charge)
        if atom.element.symbol == 'C' and atom.lone_pairs == 1 and mol.multiplicity == 1: rd_atom.SetNumRadicalElectrons(
            2)
        rdkitmol.AddAtom(rd_atom)
        if remove_h and atom.symbol == 'H':
            pass
        else:
            rd_atom_indices[atom] = index

    rd_bonds = Chem.rdchem.BondType
    orders = {'S': rd_bonds.SINGLE, 'D': rd_bonds.DOUBLE, 'T': rd_bonds.TRIPLE, 'B': rd_bonds.AROMATIC,
              'Q': rd_bonds.QUADRUPLE}
    # Add the bonds
    for atom1 in mol.vertices:
        for atom2, bond in atom1.edges.items():
            if bond.is_hydrogen_bond():
                continue
            index1 = atoms.index(atom1)
            index2 = atoms.index(atom2)
            if index1 < index2:
                order_string = bond.get_order_str()
                order = orders[order_string]
                rdkitmol.AddBond(index1, index2, order)

    # Make editable mol into a mol and rectify the molecule
    rdkitmol = rdkitmol.GetMol()
    if sanitize:
        Chem.SanitizeMol(rdkitmol)
    if remove_h:
        rdkitmol = Chem.RemoveHs(rdkitmol, sanitize=sanitize)
    if return_mapping:
        return rdkitmol, rd_atom_indices
    return rdkitmol


def from_rdkit_mol(mol, rdkitmol, raise_atomtype_exception=True):
    """
    Convert a RDKit Mol object `rdkitmol` to a molecular structure. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    This Kekulizes everything, removing all aromatic atom types.
    """
    cython.declare(i=cython.int,
                   radical_electrons=cython.int,
                   charge=cython.int,
                   lone_pairs=cython.int,
                   number=cython.int,
                   order=cython.float,
                   atom=mm.Atom,
                   atom1=mm.Atom,
                   atom2=mm.Atom,
                   bond=mm.Bond)

    mol.vertices = []

    # Add hydrogen atoms to complete molecule if needed
    rdkitmol.UpdatePropertyCache(strict=False)
    rdkitmol = Chem.AddHs(rdkitmol)
    Chem.rdmolops.Kekulize(rdkitmol, clearAromaticFlags=True)

    # iterate through atoms in rdkitmol
    for i in range(rdkitmol.GetNumAtoms()):
        rdkitatom = rdkitmol.GetAtomWithIdx(i)

        # Use atomic number as key for element
        number = rdkitatom.GetAtomicNum()
        isotope = rdkitatom.GetIsotope()
        element = elements.get_element(number, isotope or -1)

        # Process charge
        charge = rdkitatom.GetFormalCharge()
        radical_electrons = rdkitatom.GetNumRadicalElectrons()

        atom = mm.Atom(element, radical_electrons, charge, '', 0)
        mol.vertices.append(atom)

        # Add bonds by iterating again through atoms
        for j in range(0, i):
            rdkitbond = rdkitmol.GetBondBetweenAtoms(i, j)
            if rdkitbond is not None:
                order = 0

                # Process bond type
                rdbondtype = rdkitbond.GetBondType()
                if rdbondtype.name == 'SINGLE':
                    order = 1
                elif rdbondtype.name == 'DOUBLE':
                    order = 2
                elif rdbondtype.name == 'TRIPLE':
                    order = 3
                elif rdbondtype.name == 'QUADRUPLE':
                    order = 4
                elif rdbondtype.name == 'AROMATIC':
                    order = 1.5

                bond = mm.Bond(mol.vertices[i], mol.vertices[j], order)
                mol.add_bond(bond)

    # We need to update lone pairs first because the charge was set by RDKit
    mol.update_lone_pairs()
    # Set atom types and connectivity values
    mol.update(raise_atomtype_exception=raise_atomtype_exception)

    # Assume this is always true
    # There are cases where 2 radical_electrons is a singlet, but
    # the triplet is often more stable,
    mol.multiplicity = mol.get_radical_count() + 1
    # mol.update_atomtypes()

    return mol


def debug_rdkit_mol(rdmol, level=logging.INFO):
    """
    Takes an rdkit molecule object and logs some debugging information
    equivalent to calling rdmol.Debug() but uses our logging framework.
    Default logging level is INFO but can be controlled with the `level` parameter.
    Also returns the message as a string, should you want it for something.
    """
    import tempfile
    import os
    my_temp_file = tempfile.NamedTemporaryFile()
    try:
        old_stdout_file_descriptor = os.dup(sys.stdout.fileno())
    except:
        message = "Can't access the sys.stdout file descriptor, so can't capture RDKit debug info"
        print(message)
        rdmol.Debug()
        return message
    os.dup2(my_temp_file.fileno(), sys.stdout.fileno())
    rdmol.Debug()
    os.dup2(old_stdout_file_descriptor, sys.stdout.fileno())
    my_temp_file.file.seek(0)
    message = my_temp_file.file.read().decode('utf-8')
    message = "RDKit Molecule debugging information:\n" + message
    logging.log(level, message)
    return message


def to_ob_mol(mol, return_mapping=False):
    """
    Convert a molecular structure to an OpenBabel OBMol object. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    """
    if openbabel is None:
        raise DependencyError('OpenBabel is not installed. Please install or use RDKit.')

    # Sort the atoms to ensure consistent output
    mol.sort_atoms()
    atoms = mol.vertices

    ob_atom_ids = {}  # dictionary of OB atom IDs
    obmol = openbabel.OBMol()
    for atom in atoms:
        a = obmol.NewAtom()
        if atom.element.symbol == 'X':
            a.SetAtomicNum(78)  # not sure how to do this with linear scaling when this might not be Pt
        else:
            a.SetAtomicNum(atom.number)
        if atom.element.isotope != -1:
            a.SetIsotope(atom.element.isotope)
        a.SetFormalCharge(atom.charge)
        ob_atom_ids[atom] = a.GetId()
    orders = {1: 1, 2: 2, 3: 3, 4: 4, 1.5: 5}
    for atom1 in mol.vertices:
        for atom2, bond in atom1.edges.items():
            if bond.is_hydrogen_bond():
                continue
            index1 = atoms.index(atom1)
            index2 = atoms.index(atom2)
            if index1 < index2:
                order = orders[bond.order]
                obmol.AddBond(index1 + 1, index2 + 1, order)

    obmol.AssignSpinMultiplicity(True)

    if return_mapping:
        return obmol, ob_atom_ids

    return obmol


def from_ob_mol(mol, obmol, raise_atomtype_exception=True):
    """
    Convert a OpenBabel Mol object `obmol` to a molecular structure. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.

    It estimates radical placement based on undervalence of atoms,
    and assumes overall spin multiplicity is radical count + 1
    """
    # Below are the declared variables for cythonizing the module
    # cython.declare(i=cython.int)
    # cython.declare(radical_electrons=cython.int, charge=cython.int, lone_pairs=cython.int)
    # cython.declare(atom=mm.Atom, atom1=mm.Atom, atom2=mm.Atom, bond=mm.Bond)
    if openbabel is None:
        raise DependencyError('OpenBabel is not installed. Please install or use RDKit.')

    mol.vertices = []

    # Add hydrogen atoms to complete molecule if needed
    obmol.AddHydrogens()
    # TODO Chem.rdmolops.Kekulize(obmol, clearAromaticFlags=True)

    # iterate through atoms in obmol
    for obatom in openbabel.OBMolAtomIter(obmol):
        # Use atomic number as key for element
        number = obatom.GetAtomicNum()
        isotope = obatom.GetIsotope()
        element = elements.get_element(number, isotope or -1)
        # Process charge
        charge = obatom.GetFormalCharge()
        # Calculate the radical electrons due to undervalence,
        # ignoring whatever may be set on obatom.GetSpinMultiplicity()
        valence = obatom.GetTotalValence()
        radical_electrons = openbabel.GetTypicalValence(number, valence, charge) - valence

        atom = mm.Atom(element, radical_electrons, charge, '', 0)
        mol.vertices.append(atom)

    # iterate through bonds in obmol
    for obbond in openbabel.OBMolBondIter(obmol):
        # Process bond type
        oborder = obbond.GetBondOrder()
        if oborder not in [1, 2, 3, 4] and obbond.IsAromatic():
            oborder = 1.5

        bond = mm.Bond(mol.vertices[obbond.GetBeginAtomIdx() - 1],
                       mol.vertices[obbond.GetEndAtomIdx() - 1],
                       oborder)  # python array indices start at 0
        mol.add_bond(bond)

    # Set atom types and connectivity values
    mol.update_connectivity_values()
    mol.update_atomtypes(log_species=True, raise_exception=raise_atomtype_exception)
    mol.update_multiplicity()
    mol.identify_ring_membership()

    # Assume this is always true
    # There are cases where 2 radical_electrons is a singlet, but
    # the triplet is often more stable,
    mol.multiplicity = mol.get_radical_count() + 1

    return mol
