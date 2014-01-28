#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.molecule.atomtype module.
"""

import unittest
import os
import os.path

import rmgpy.molecule
from rmgpy.molecule import  Molecule
from rmgpy.molecule.draw import *

################################################################################

class TestMoleculeDrawer(unittest.TestCase):
    """
    Contains unit tests of the MoleculeDrawer class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.drawer = MoleculeDrawer()
        self.molecule = Molecule(SMILES='CC(=O)CC')
        
    def testDrawPNG(self):
        """
        Test we can create PNG files from molecules.
        """
        path = 'test_molecule.png'
        from cairo import ImageSurface
        if os.path.exists(path):
            os.unlink(path)
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='png', path=path)
        self.assertTrue(os.path.exists(path), "File doesn't exist")
        os.unlink(path)
        self.assertIsInstance(surface, ImageSurface)

    def testDrawPDF(self):
        """
        Test we can create PDF files from molecules.
        """
        from cairo import PDFSurface
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='pdf')
        self.assertIsInstance(surface, PDFSurface)
        self.assertGreater(width, height)

    def testDrawPolycycle(self):
        """
        Test we can draw a polycyclic molecule
        """
        from cairo import PDFSurface
        polycycle = Molecule(SMILES="C123CC4CC1COCC2CCC34")
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='pdf')
        self.assertIsInstance(surface, PDFSurface)
        self.assertGreater(width, height)


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
