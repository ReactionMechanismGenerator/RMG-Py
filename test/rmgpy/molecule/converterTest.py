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
This module contains unit test for the converter module.
"""


from rmgpy.molecule.converter import (
    debug_rdkit_mol,
    to_rdkit_mol,
    from_rdkit_mol,
    to_ob_mol,
    from_ob_mol,
)
from rmgpy.molecule.molecule import Molecule


class RDKitTest:
    def test_debugger(self):
        """Test the debug_rdkit_mol(rdmol) function doesn't crash

        We can't really test it in the unit testing framework, because
        that already captures and redirects standard output, and that
        conflicts with the function, but this checks it doesn't crash.
        """
        import rdkit.Chem
        import logging

        rdmol = rdkit.Chem.MolFromSmiles("CCC")
        message = debug_rdkit_mol(rdmol, level=logging.INFO)
        assert message is not None

    def test_lone_pair_retention(self):
        """Test that we don't lose any lone pairs on round trip RDKit conversion."""
        mol = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""
        )
        rdmol = to_rdkit_mol(mol)

        mol2 = from_rdkit_mol(Molecule(), rdmol)
        assert mol.is_isomorphic(mol2)

    def test_atom_mapping_1(self):
        """Test that to_rdkit_mol returns correct indices and atom mappings."""
        bond_order_dict = {"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3, "AROMATIC": 1.5}
        mol = Molecule().from_smiles("C1CCC=C1C=O")
        rdkitmol, rd_atom_indices = to_rdkit_mol(mol, remove_h=False, return_mapping=True)
        for atom in mol.atoms:
            # Check that all atoms are found in mapping
            assert atom in rd_atom_indices
            # Check that all bonds are in rdkitmol with correct mapping and order
            for connected_atom, bond in atom.bonds.items():
                bond_type = str(rdkitmol.GetBondBetweenAtoms(rd_atom_indices[atom], rd_atom_indices[connected_atom]).GetBondType())
                rdkit_bond_order = bond_order_dict[bond_type]
                assert bond.order == rdkit_bond_order

        # Test for remove_h = True
        rdkitmol2, rd_atom_indices2 = to_rdkit_mol(mol, remove_h=True, return_mapping=True)
        for atom in mol.atoms:
            # Check that all non-hydrogen atoms are found in mapping
            if atom.symbol != "H":
                assert atom in rd_atom_indices2
                # Check that all bonds connected to non-hydrogen have the correct mapping and order
                for connected_atom, bond in atom.bonds.items():
                    if connected_atom.symbol != "H":
                        bond_type = str(rdkitmol2.GetBondBetweenAtoms(rd_atom_indices2[atom], rd_atom_indices2[connected_atom]).GetBondType())
                        rdkit_bond_order = bond_order_dict[bond_type]
                        assert bond.order == rdkit_bond_order

    def test_atom_mapping_2(self):
        """Test that to_rdkit_mol returns correct indices and atom mappings when hydrogens are removed."""
        adjlist = """
1 H u0 p0 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 O u0 p2 c0 {2,S} {6,S}
6 H u0 p0 c0 {5,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        rdkitmol, rd_atom_indices = to_rdkit_mol(mol, remove_h=True, return_mapping=True)

        heavy_atoms = [at for at in mol.atoms if at.number != 1]
        for at1 in heavy_atoms:
            for at2 in heavy_atoms:
                if mol.has_bond(at1, at2):
                    try:
                        rdkitmol.GetBondBetweenAtoms(rd_atom_indices[at1], rd_atom_indices[at2])
                    except RuntimeError:
                        assert False, "RDKit failed in finding the bond in the original atom!"

    def test_atom_mapping_3(self):
        """Test that to_rdkit_mol with save_order=True retains the atom order and create the correct RDKit Molecule"""
        adjlist = """1  H u0 p0 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,T}
3  N u0 p1 c0 {2,T}
"""
        mol = Molecule().from_adjacency_list(adjlist)
        rdkitmol, _ = to_rdkit_mol(mol, remove_h=False, return_mapping=True, save_order=True)

        assert [atom.number for atom in mol.atoms] == [1, 6, 7]
        assert [rdkitmol.GetAtomWithIdx(idx).GetAtomicNum() for idx in range(3)] == [1, 6, 7]


class ConverterTest:
    def setup_class(self):
        """Function run before each test in this class."""
        self.test_mols = [
            Molecule().from_smiles("C"),
            Molecule().from_smiles("O"),
            Molecule().from_smiles("N"),
            Molecule().from_smiles("S"),
            Molecule().from_smiles("[CH2]C"),
            Molecule().from_smiles("[CH]C"),
            Molecule().from_smiles("C=CC=C"),
            Molecule().from_smiles("C#C[CH2]"),
            Molecule().from_smiles("c1ccccc1"),
            Molecule().from_smiles("[13CH3]C"),
            Molecule().from_smiles("O=CCO").generate_h_bonded_structures()[0],
        ]
        self.test_Hbond_free_mol = Molecule().from_smiles("O=CCO")

    def test_rdkit_round_trip(self):
        """Test conversion to and from RDKitMol"""
        for mol in self.test_mols:
            rdkit_mol = to_rdkit_mol(mol)
            new_mol = from_rdkit_mol(Molecule(), rdkit_mol)
            assert mol.is_isomorphic(new_mol) or self.test_Hbond_free_mol.is_isomorphic(new_mol)
            assert mol.get_element_count() == new_mol.get_element_count()

    def test_ob_round_trip(self):
        """Test conversion to and from OBMol"""
        for mol in self.test_mols:
            ob_mol = to_ob_mol(mol)
            new_mol = from_ob_mol(Molecule(), ob_mol)
            assert mol.is_isomorphic(new_mol) or self.test_Hbond_free_mol.is_isomorphic(new_mol)
            assert mol.get_element_count() == new_mol.get_element_count()
