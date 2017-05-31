#!/usr/bin/env python
# encoding: utf-8 -*-

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

import numpy
import unittest
import os

from rmgpy.cantherm.molepro import MoleProLog
import rmgpy.constants as constants

################################################################################

class MoleProTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.io.gaussian module, used for reading
    and writing Molepro files.
    """
    
    def testLoadDzFromMoleProLog_F12(self):
        """
        Uses a Molepro log file for ethylene_dz (C2H4) to test that F12a
        energy can be properly read.
        """
        
        log=MoleProLog(os.path.join(os.path.dirname(__file__),'data','ethylene_f12_dz.out'))
        E0=log.loadEnergy()
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -78.474353559604, 5)
    
    def testLoadQzFromMoleProLog_F12(self):
        """
        Uses a Molepro log file for ethylene_qz (C2H4) to test that F12b
        energy can be properly read.
        """
        
        log=MoleProLog(os.path.join(os.path.dirname(__file__),'data','ethylene_f12_qz.out'))
        E0=log.loadEnergy()
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -78.472682755635, 5)

    def testLoadRadFromMoleProLog_F12(self):
        """
        Uses a Molepro log file for OH (C2H4) to test that radical
        energy can be properly read.
        """
        
        log=MoleProLog(os.path.join(os.path.dirname(__file__),'data','OH_f12.out'))
        E0=log.loadEnergy()
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -75.663696424380, 5)
