#!/usr/bin/env python
# encoding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.parallel module.
"""

import os
import sys
import unittest

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


def funcShared():
    """
    Broadcast the data with the given key, 
    and retrieve it again by querying the key again.
    """
    data = 'foo'
    key = 'bar'

    broadcast(data, key)

    result = True
    try:
        assert data == shared.getConst(key)
    except AssertionError:
        return False
    
    return True

class BroadcastTest(TestScoopCommon):

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    def test_generic(self):
        """
        Test that we can broadcast a simple string.
        """

        result = futures._startup(funcShared)
        self.assertEquals(result, True)

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main()

