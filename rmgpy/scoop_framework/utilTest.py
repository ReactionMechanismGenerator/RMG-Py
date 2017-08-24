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
from external.wip import work_in_progress

from rmgpy.scoop_framework.framework import TestScoopCommon

try:
    from scoop import futures, _control, shared
except ImportError, e:
    import logging as logging
    logging.debug("Could not properly import SCOOP.")

from .util import *

def boom():
    return 0 / 0

class WorkerWrapperTest(unittest.TestCase):

    def test_WorkerWrapper(self):
        """
        Test that the WorkerWrapper correctly redirects the error
        message.
        """
        f = WorkerWrapper(boom)    
        with self.assertRaises(ZeroDivisionError):
            f()


def funcBroadcast():
    """
    Broadcast the data with the given key, 
    and retrieve it again by querying the key again.
    """
    data = 'foo'
    key = 'bar'

    broadcast(data, key)

    try:
        assert data == shared.getConst(key)
    except AssertionError:
        return False
    
    return True

def funcRetrieve():
    """
    Broadcast the data with the given key, 
    retrieve it again by querying the key again.
    """
    
    data = 'foo'
    key = 'bar'


    broadcast(data, key)

    try:
        assert data == get(key)
    except AssertionError:
        return False
    
    return True    

@work_in_progress
class BroadcastTest(TestScoopCommon):

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "test currently only runs on linux")
    def test_generic(self):
        """
        Test that we can broadcast a simple string.
        """

        result = futures._startup(funcBroadcast)
        self.assertEquals(result, True)

@work_in_progress
class GetTest(TestScoopCommon):

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    def test_generic(self):
        """
        Test that we can retrieve a simple shared string.
        """

        result = futures._startup(funcRetrieve)
        self.assertEquals(result, True)        

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main()
