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
        
        aug_inchi = 'InChI=1S/CH3/h1H3/mult2'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 2)
        
        aug_inchi = 'InChI=1S/CH3/h1H3/mult2/u0'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 2)
        
        aug_inchi = 'InChI=1S/CH2/h1H2/mult3'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 3)
        
        aug_inchi = 'InChI=1S/CH2/h1H2/mult1'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 1)
        
        aug_inchi = 'InChI=1S/CH2/h1H2/mult3/u0'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 3)
        
        aug_inchi = 'InChI=1S/CH2/h1H2/mult1/u0'
        mol = parser.fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(mol.multiplicity == 1)
