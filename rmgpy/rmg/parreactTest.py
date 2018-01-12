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
import logging
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.data.kinetics import TemplateReaction
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import Species

from rmgpy.scoop_framework.framework import TestScoopCommon

from rmgpy.rmg.react import *

try:
    from scoop import futures, _control, shared
except ImportError, e:
    import logging as logging
    logging.debug("Could not properly import SCOOP.")

TESTFAMILY = 'H_Abstraction'

def tearDown():
    """
    Reset the loaded database
    """
    import rmgpy.data.rmg
    rmgpy.data.rmg.database = None

def load():
    """
    A method that is run before each unit test in this class.
    """
    tearDown()
    # set-up RMG object
    rmg = RMG()

    # load kinetic database and forbidden structures
    rmg.database = RMGDatabase()
    path = os.path.join(settings['database.directory'])

    # forbidden structure loading
    rmg.database.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
    # kinetics family loading
    rmg.database.loadKinetics(os.path.join(path, 'kinetics'),
                                   kineticsFamilies=[TESTFAMILY],
                                   reactionLibraries=[]
                                   )

def generate():
    """
    Test that reaction generation from the available families works.
    """
    load()
    spcA = Species().fromSMILES('[OH]')
    spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
    spcTuples = [(spcA, spc) for spc in spcs]

    reactionList = list(react(*spcTuples))

    if not reactionList: return False

    for rxn in reactionList:
        if not isinstance(rxn, TemplateReaction): return False

    return True

@work_in_progress
class ParallelReactTest(TestScoopCommon):

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()
    
    @unittest.skipUnless(sys.platform.startswith("linux"),
                          "test currently only runs on linux")
    def test(self):
        """
        Test that we can generate reactions in parallel.
        """

        result = futures._startup(generate)
        self.assertEquals(result, True)        

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main()
