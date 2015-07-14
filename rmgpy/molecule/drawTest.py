#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.molecule.atomtype module.
"""

import unittest
import os
import os.path

from rmgpy.molecule import  Molecule
from rmgpy.molecule.draw import MoleculeDrawer

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
        try:
            from cairocffi import ImageSurface
        except ImportError:
            from cairo import ImageSurface
        path = 'test_molecule.png'
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
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        path = 'test_molecule.pdf'
        if os.path.exists(path):
            os.unlink(path)
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='pdf', path=path)
        self.assertIsInstance(surface, PDFSurface)
        self.assertGreater(width, height)
        os.unlink(path)

    def testDrawPolycycle(self):
        """
        Test we can draw a polycyclic molecule
        """
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        path = 'test_molecule.pdf'
        if os.path.exists(path):
            os.unlink(path)
        polycycle = Molecule(SMILES="C123CC4CC1COCC2CCC34")
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='pdf', path=path)
        self.assertIsInstance(surface, PDFSurface)
        self.assertGreater(width, height)
        os.unlink(path)

    def testDrawPDFwithoutFile(self):
        """
        Test we can create PDF surface without a temporary file (newer versions of PyCairo?)
        """
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='pdf')
        self.assertIsInstance(surface, PDFSurface)
        self.assertGreater(width, height)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
