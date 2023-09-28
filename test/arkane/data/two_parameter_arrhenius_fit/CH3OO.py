#!/usr/bin/env python
# -*- coding: utf-8 -*-

energy = Log('ch3oo_sp.log')

geometry = Log('ch3oo_freq.log')

# the `check_for_errors` keyword argument is not needed for this file. It is added to
# 1) demonstrate the syntax for instructing Arkane to ignore any error messages in the file
# 2) test that this syntax works properly during a functional test
# this argument is optional, and it is True by default
frequencies = Log('ch3oo_freq.log', check_for_errors=False)
