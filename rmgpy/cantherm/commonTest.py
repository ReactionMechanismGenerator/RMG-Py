#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This script contains unit tests of the :mod:`rmgpy.quantity` module.
"""
import unittest
import numpy
import rmgpy.constants as constants
################################################################################

class CommonTest(unittest.TestCase):
    """
    Contains unit tests of the Cantherm common functions.
    """
    
    def test_checkConformerEnergy(self):
        """
        test the checkConformerEnergy function with an list of energies.
        """    
        Vlist = [-272.2779012225, -272.2774933703, -272.2768397635, -272.2778432059, -272.278645477, -272.2789602654, -272.2788749196, -272.278496709, -272.2779350675, -272.2777008843, -272.2777167286, -272.2780937643, -272.2784838846, -272.2788050464, -272.2787865352, -272.2785091607, -272.2779977452, -272.2777957743, -272.2779134906, -272.2781827547, -272.278443339, -272.2788244214, -272.2787748749]
        Vlist = numpy.array(Vlist, numpy.float64)
        Vdiff = (Vlist[0] - numpy.min(Vlist))*constants.E_h*constants.Na/1000
        self.assertAlmostEqual(Vdiff / 2.7805169838282797, 1, 5)

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
