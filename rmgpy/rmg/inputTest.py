#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import unittest 

from rmgpy.rmg.main import RMG
from rmgpy.rmg import input as inp

###################################################

def setUpModule(self):
    """
    A method that is run before the class.
    """
    # set-up RMG object and get global rmg object in input.py file
    # so methods can be tested
    global rmg
    rmg = RMG()
    inp.setGlobalRMG(rmg)

def tearDownModule(self):
    # remove RMG object
    global rmg
    rmg = None

class TestInputDatabase(unittest.TestCase):
    """
    Contains unit tests rmgpy.rmg.input.database
    """

    def tearDown(self):
        # remove the reactionLibraries value
        global rmg
        rmg.reactionLibraries = None

    def testImportingDatabaseReactionLibrariesFromString(self):
        """
        Test that we can import Reaction Libraries using the non-tuple form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=['test'])
        self.assertIsInstance(rmg.reactionLibraries[0], tuple)
        self.assertFalse(rmg.reactionLibraries[0][1])
        
    def testImportingDatabaseReactionLibrariesFromFalseTuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple False form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[('test',False)])
        self.assertIsInstance(rmg.reactionLibraries[0], tuple)
        self.assertFalse(rmg.reactionLibraries[0][1])
        
    def testImportingDatabaseReactionLibrariesFromTrueTuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple True form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[('test',True)])
        self.assertIsInstance(rmg.reactionLibraries[0], tuple)
        self.assertTrue(rmg.reactionLibraries[0][1])

class TestInputMLEstimator(unittest.TestCase):
    """
    Contains unit tests rmgpy.rmg.input.mlEstimator
    """
    def tearDown(self):
        # remove the reactionLibraries value
        global rmg
        rmg.ml_estimator = None

    def testMLEstimator(self):
        """
        Test that we can input.
        """
        from rmgpy.ml.estimator import MLEstimator
        global rmg
        # add database properties to RMG
        inp.mlEstimator(thermo=True)
        self.assertIsInstance(rmg.ml_estimator, MLEstimator)
        self.assertIsInstance(rmg.ml_settings, dict)

class TestInputThemoCentralDatabase(unittest.TestCase):
    """
    Contains unit tests rmgpy.rmg.input.thermoCentralDatabase
    """
    def tearDown(self):
        # remove the reactionLibraries value
        global rmg
        rmg.thermoCentralDatabase = None
        
    def testThemoCentralDatabase(self):
        """
        Test that we can input.
        """
        global rmg
        # add database properties to RMG
        inp.thermoCentralDatabase(
                host='some_host',
                port=0,
                username='some_usr',
                password='some_pw',
                application='some_app'
        )
        self.assertEqual(rmg.thermoCentralDatabase.host, 'some_host')
        self.assertEqual(rmg.thermoCentralDatabase.port, 0)
        self.assertEqual(rmg.thermoCentralDatabase.username, 'some_usr')
        self.assertEqual(rmg.thermoCentralDatabase.password, 'some_pw')
        self.assertEqual(rmg.thermoCentralDatabase.application, 'some_app')
        self.assertEqual(rmg.thermoCentralDatabase.client, None)


if __name__ == '__main__':
    unittest.main()
