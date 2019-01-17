#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains unit tests of the rmgpy.molecule.atomtype module.
"""

import unittest
import os
import os.path

from rmgpy.molecule import  Molecule
from rmgpy.molecule.draw import MoleculeDrawer
from rmgpy.species import Species
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
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='png', target=path)
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
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='pdf', target=path)
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
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(self.molecule, format='pdf', target=path)
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

    def testDrawNonStandardBonds(self):
        
        spec = Species().fromSMILES('[CH2]C=C[CH2]')
        hybrid = spec.getResonanceHybrid()
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        surface, cr, (xoff, yoff, width, height) = self.drawer.draw(hybrid, format='pdf')
        self.assertIsInstance(surface, PDFSurface)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
