#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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

import os
import unittest 

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import *

###################################################


class TestSpecies(unittest.TestCase):
    """
    Contains unit tests of the Species class.
    """

    def setUp(self):
        """
        A method that is run before each unit test in this class.
        """
        # set-up RMG object
        self.rmg = RMG()

        # load kinetic database and forbidden structures
        self.rmg.database = RMGDatabase()
        path = os.path.join(settings['database.directory'])

        # forbidden structure loading
        self.rmg.database.loadThermo(os.path.join(path, 'thermo'))
        

    def testGetThermoData(self):
        """
        Test that getThermoData method of Species works.
        """
        spc = Species().fromSMILES('CCC')

        self.assertFalse(spc.thermo)
        spc.getThermoData(self.rmg.database)
        self.assertTrue(spc.thermo)
        thermo = spc.thermo
        spc.getThermoData(self.rmg.database)

        self.assertEquals(id(thermo), id(spc.thermo))
        
        spc.thermo = None
        spc.getThermoData(self.rmg.database)
        self.assertNotEquals(id(thermo), id(spc.thermo))

    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

if __name__ == '__main__':
    unittest.main()
