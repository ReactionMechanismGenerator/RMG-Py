import unittest

import rmgpy.molecule.parser as parser
from rmgpy.molecule.molecule import Molecule


class parserTest(unittest.TestCase):

    def test_fromAugmentedInChI(self):
        aug_inchi = 'InChI=1S/CH4/h1H4'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(not mol.InChI == '')

        aug_inchi = 'InChI=1/CH4/h1H4'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(not mol.InChI == '')
 
    def test_OptionalMultLayer(self):
        """
        Test that multiplicity layer be optional in cases where the layer is not needed.
        """

        aug_inchi = 'InChI=1S/CH3/h1H3'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 2)
        
        aug_inchi = 'InChI=1S/CH3/h1H3'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 2)
        
        aug_inchi = 'InChI=1S/CH3/h1H3/u0'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 2)
        
        aug_inchi = 'InChI=1S/CH2/h1H2'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 3)
        
        aug_inchi = 'InChI=1S/CH2/h1H2'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 1)
        
        aug_inchi = 'InChI=1S/CH2/h1H2/u0'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 3)
        
        aug_inchi = 'InChI=1S/CH2/h1H2/u0'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 1)

    def test_toRDKitMol(self):
        """
        Test that toRDKitMol returns correct indices and atom mappings.
        """
        bondOrderDict = {'SINGLE':'S','DOUBLE':'D','TRIPLE':'T','AROMATIC':'B'}
        mol = parser.fromSMILES(Molecule(), 'C1CCC=C1C=O')
        rdkitmol, rdAtomIndices = mol.toRDKitMol(removeHs=False, returnMapping=True, sanitize=True)
        for atom in mol.atoms:
            # Check that all atoms are found in mapping
            self.assertTrue(atom in rdAtomIndices)
            # Check that all bonds are in rdkitmol with correct mapping and order
            for connectedAtom, bond in atom.bonds.iteritems():
                bondType = str(rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom],rdAtomIndices[connectedAtom]).GetBondType())
                rdkitBondOrder = bondOrderDict[bondType]
                self.assertEqual(bond.order, rdkitBondOrder)
        
        # Test for removeHs = True        
        rdkitmol2, rdAtomIndices2 = mol.toRDKitMol(removeHs=True, returnMapping=True, sanitize=True)
        
        for atom in mol.atoms:
            # Check that all non-hydrogen atoms are found in mapping
            if atom.symbol != 'H':
                self.assertTrue(atom in rdAtomIndices)
                # Check that all bonds connected to non-hydrogen have the correct mapping and order
                for connectedAtom, bond in atom.bonds.iteritems():
                    if connectedAtom.symbol != 'H':
                        bondType = str(rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom],rdAtomIndices[connectedAtom]).GetBondType())
                        rdkitBondOrder = bondOrderDict[bondType]
                        self.assertEqual(bond.order, rdkitBondOrder)   