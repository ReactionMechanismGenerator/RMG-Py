#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

import rmgpy.rmg.input as inp
from rmgpy.rmg.main import RMG


###################################################

def setUpModule(self):
    """
    A method that is run before the class.
    """
    # set-up RMG object and get global rmg object in input.py file
    # so methods can be tested
    global rmg
    rmg = RMG()
    inp.set_global_rmg(rmg)


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
        rmg.reaction_libraries = None

    def test_importing_database_reaction_libraries_from_string(self):
        """
        Test that we can import Reaction Libraries using the non-tuple form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=['test'])
        self.assertIsInstance(rmg.reaction_libraries[0], tuple)
        self.assertFalse(rmg.reaction_libraries[0][1])

    def test_importing_database_reaction_libraries_from_false_tuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple False form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[('test', False)])
        self.assertIsInstance(rmg.reaction_libraries[0], tuple)
        self.assertFalse(rmg.reaction_libraries[0][1])

    def test_importing_database_reaction_libraries_from_true_tuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple True form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[('test', True)])
        self.assertIsInstance(rmg.reaction_libraries[0], tuple)
        self.assertTrue(rmg.reaction_libraries[0][1])


class TestInputMLEstimator(unittest.TestCase):
    """
    Contains unit tests rmgpy.rmg.input.mlEstimator
    """

    def tearDown(self):
        # remove the reactionLibraries value
        global rmg
        rmg.ml_estimator = None

    def test_ml_estimator(self):
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
    Contains unit tests rmgpy.rmg.input.thermo_central_database
    """

    def tearDown(self):
        # remove the reactionLibraries value
        global rmg
        rmg.thermo_central_database = None

    def test_themo_central_database(self):
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
        self.assertEqual(rmg.thermo_central_database.host, 'some_host')
        self.assertEqual(rmg.thermo_central_database.port, 0)
        self.assertEqual(rmg.thermo_central_database.username, 'some_usr')
        self.assertEqual(rmg.thermo_central_database.password, 'some_pw')
        self.assertEqual(rmg.thermo_central_database.application, 'some_app')
        self.assertEqual(rmg.thermo_central_database.client, None)


if __name__ == '__main__':
    unittest.main()
