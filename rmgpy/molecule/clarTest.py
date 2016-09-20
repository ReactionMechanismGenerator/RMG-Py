#!/usr/bin/env python
# encoding: utf-8

import unittest
from external.wip import work_in_progress

import clar
from .molecule import Molecule

################################################################################

class TestClar(unittest.TestCase):
    """
    Contains unit tests for Clar structure methods.
    """

    def setUp(self):
        pass

    def testClarTransformation(self):
        """
        Basic test that aromatic ring is generated.
        """
        mol = Molecule().fromSMILES('c1ccccc1')
        sssr = mol.getSmallestSetOfSmallestRings()
        clar.clarTransformation(mol, sssr[0])

        self.assertTrue(mol.isAromatic())

    def testPhenanthrene(self):
        """
        Test phenanthrene, which is a basic case that should work.
        """
        mol = Molecule().fromSMILES('C1=CC=C2C(C=CC3=CC=CC=C32)=C1')
        newmol = clar.generateClarStructures(mol)

        self.assertEqual(len(newmol), 1)
        self.assertTrue(newmol[0].isAromatic())

    def testPhenalene(self):
        """
        Test phenalene, which currently does not have feasible starting point.
        """
        mol = Molecule().fromSMILES('C1=CC2=CC=CC3CC=CC(=C1)C=32')
        newmol = clar.generateClarStructures(mol)

        self.assertEqual(len(newmol), 2)
        self.assertTrue(newmol[0].isAromatic())
        self.assertTrue(newmol[1].isAromatic())

    def testCorannulene(self):
        """
        Test corannulene, which does not give integer results after initial optimization.
        """
        mol = Molecule().fromSMILES('C1=CC2=CC=C3C=CC4=C5C6=C(C2=C35)C1=CC=C6C=C4')
        newmol = clar.generateClarStructures(mol)

        self.assertEqual(len(newmol), 5)
        self.assertTrue(newmol[0].isAromatic())
        self.assertTrue(newmol[1].isAromatic())
        self.assertTrue(newmol[2].isAromatic())
        self.assertTrue(newmol[3].isAromatic())
        self.assertTrue(newmol[4].isAromatic())

