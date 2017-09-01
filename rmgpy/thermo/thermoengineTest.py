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

"""
This module contains unit tests of the rmgpy.parallel module.
"""

import os
import sys
import unittest
import random
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import RMG
from rmgpy.scoop_framework.framework import TestScoopCommon

from rmgpy.species import Species
from rmgpy.thermo.thermoengine import submit

try:
    from scoop import futures, _control, shared
except ImportError, e:
    import logging as logging
    logging.debug("Could not properly import SCOOP.")


def load():
    tearDown()
    rmg = RMG()#for solvent
    database = RMGDatabase()
    database.loadThermo(os.path.join(settings['database.directory'], 'thermo'))
    database.loadTransport(os.path.join(settings['database.directory'], 'transport'))
    database.loadSolvation(os.path.join(settings['database.directory'], 'solvation'))

def tearDown():
    """
    Reset the loaded database
    """
    import rmgpy.data.rmg
    rmgpy.data.rmg.database = None

def funcSubmit():
    """
    Test that we can submit a number of species.
    """
    load()

    spcs = [
            Species().fromSMILES('C'),\
            Species().fromSMILES('CC'), \
            Species().fromSMILES('CCC')
            ]
    
    for spc in spcs:
        submit(spc)

    return True

def funcGet():
    """
    Test if we can retrieve thermo of species even before we have submitted them explicitly.
    """
    load()

    spcs = [
            Species().fromSMILES('C'),
            Species().fromSMILES('CC'), \
            Species().fromSMILES('CCC')
            ]
    
    output = []
    for spc in spcs:
        data = spc.getThermoData()
        output.append((spc, data))

    for spc, data in output:
        if not data:
            return False

    return True

def funcSubmitGet():
    """
    Test if we can retrieve thermo of species after submitting some of them.
    """
    load()

    spcs = [
            Species().fromSMILES('C'),\
            Species().fromSMILES('CC'), \
            Species().fromSMILES('CCC')
            ]
    
    for spc in spcs:
        submit(spc)

    absent = Species().fromSMILES('[CH3]')
    data = absent.getThermoData()
    if not data: return False

    present = Species().fromSMILES('CC')
    data = present.getThermoData()
    if not data: return False

    random.shuffle(spcs)
    for spc in spcs:
        data = spc.getThermoData()
        if not data: return False        

    return True

@work_in_progress
class AsyncThermoTest(TestScoopCommon):

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "test currently only runs on linux")
    def testSubmit(self):
        """
        Test that we can submit a request to generate
        thermo/transport for a number of species.
        """
        result = futures._startup(funcSubmit)
        self.assertEquals(result, True)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "test currently only runs on linux")
    def testGet(self):
        """
        Test that we can get the data of a number of species.
        """
        result = futures._startup(funcGet)
        self.assertEquals(result, True)

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main()
