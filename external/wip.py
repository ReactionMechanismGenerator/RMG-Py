#!/usr/bin/env python
# encoding: utf-8
#
# Decorator to mark a unit test as a "work_in_progress"
# From http://www.natpryce.com/articles/000788.html
# Copyright 2011 Nat Pryce. Posted 2011-05-30
from functools import wraps
from nose.plugins.attrib import attr
from nose.plugins.skip import SkipTest

def fail(message):
    raise AssertionError(message)

def work_in_progress(f):
    @wraps(f)
    def run_test(*args, **kwargs):
        try:
            f(*args, **kwargs)
        except Exception as e:
            raise SkipTest("WIP test failed: " + str(e))
        fail("test passed but marked as work in progress")

    return attr('work_in_progress')(run_test)