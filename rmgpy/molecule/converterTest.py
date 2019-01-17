#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import unittest

from rmgpy.exceptions import AtomTypeError
from rmgpy.molecule.converter import debugRDKitMol, toRDKitMol, fromRDKitMol, toOBMol, fromOBMol
from rmgpy.molecule.molecule import Molecule


class RDKitTest(unittest.TestCase):

    def testDebugger(self):
        """Test the debugRDKitMol(rdmol) function doesn't crash

        We can't really test it in the unit testing framework, because
        that already captures and redirects standard output, and that
        conflicts with the function, but this checks it doesn't crash.
        """
        import rdkit.Chem
        import logging
        rdmol = rdkit.Chem.MolFromSmiles('CCC')
        message = debugRDKitMol(rdmol, level=logging.INFO)
        self.assertIsNotNone(message)

    def test_lone_pair_retention(self):
        """Test that we don't lose any lone pairs on round trip RDKit conversion."""
        mol = Molecule().fromAdjacencyList(
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""")
        rdmol = toRDKitMol(mol)

        try:
            mol2 = fromRDKitMol(Molecule(), rdmol)
        except AtomTypeError as e:
            self.fail('Could not convert from RDKitMol: ' + e.message)
        else:
            self.assertTrue(mol.isIsomorphic(mol2))

    def test_atom_mapping_1(self):
        """Test that toRDKitMol returns correct indices and atom mappings."""
        bondOrderDict = {'SINGLE': 1, 'DOUBLE': 2, 'TRIPLE': 3, 'AROMATIC': 1.5}
        mol = Molecule().fromSMILES('C1CCC=C1C=O')
        rdkitmol, rdAtomIndices = toRDKitMol(mol, removeHs=False, returnMapping=True)
        for atom in mol.atoms:
            # Check that all atoms are found in mapping
            self.assertTrue(atom in rdAtomIndices)
            # Check that all bonds are in rdkitmol with correct mapping and order
            for connectedAtom, bond in atom.bonds.iteritems():
                bondType = str(
                    rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom], rdAtomIndices[connectedAtom]).GetBondType())
                rdkitBondOrder = bondOrderDict[bondType]
                self.assertEqual(bond.order, rdkitBondOrder)

        # Test for removeHs = True
        rdkitmol2, rdAtomIndices2 = toRDKitMol(mol, removeHs=True, returnMapping=True)
        for atom in mol.atoms:
            # Check that all non-hydrogen atoms are found in mapping
            if atom.symbol != 'H':
                self.assertTrue(atom in rdAtomIndices2)
                # Check that all bonds connected to non-hydrogen have the correct mapping and order
                for connectedAtom, bond in atom.bonds.iteritems():
                    if connectedAtom.symbol != 'H':
                        bondType = str(rdkitmol2.GetBondBetweenAtoms(rdAtomIndices2[atom],
                                                                     rdAtomIndices2[connectedAtom]).GetBondType())
                        rdkitBondOrder = bondOrderDict[bondType]
                        self.assertEqual(bond.order, rdkitBondOrder)

    def test_atom_mapping_2(self):
        """Test that toRDKitMol returns correct indices and atom mappings when hydrogens are removed."""
        adjlist = """
1 H u0 p0 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 O u0 p2 c0 {2,S} {6,S}
6 H u0 p0 c0 {5,S}
        """

        mol = Molecule().fromAdjacencyList(adjlist)
        rdkitmol, rdAtomIndices = toRDKitMol(mol, removeHs=True, returnMapping=True)

        heavy_atoms = [at for at in mol.atoms if at.number != 1]
        for at1 in heavy_atoms:
            for at2 in heavy_atoms:
                if mol.hasBond(at1, at2):
                    try:
                        rdkitmol.GetBondBetweenAtoms(rdAtomIndices[at1], rdAtomIndices[at2])
                    except RuntimeError:
                        self.fail("RDKit failed in finding the bond in the original atom!")

class ConverterTest(unittest.TestCase):

    def setUp(self):
        """Function run before each test in this class."""
        self.test_mols = [
            Molecule().fromSMILES('C'),
            Molecule().fromSMILES('O'),
            Molecule().fromSMILES('N'),
            Molecule().fromSMILES('S'),
            Molecule().fromSMILES('[CH2]C'),
            Molecule().fromSMILES('[CH]C'),
            Molecule().fromSMILES('C=CC=C'),
            Molecule().fromSMILES('C#C[CH2]'),
            Molecule().fromSMILES('c1ccccc1'),
            Molecule().fromSMILES('[13CH3]C')
        ]

    def test_rdkit_round_trip(self):
        """Test conversion to and from RDKitMol"""
        for mol in self.test_mols:
            rdkit_mol = toRDKitMol(mol)
            new_mol = fromRDKitMol(Molecule(), rdkit_mol)

            self.assertTrue(mol.isIsomorphic(new_mol))
            self.assertEqual(mol.get_element_count(), new_mol.get_element_count())

    def test_ob_round_trip(self):
        """Test conversion to and from OBMol"""
        for mol in self.test_mols:
            ob_mol = toOBMol(mol)
            new_mol = fromOBMol(Molecule(), ob_mol)

            self.assertTrue(mol.isIsomorphic(new_mol))
            self.assertEqual(mol.get_element_count(), new_mol.get_element_count())
