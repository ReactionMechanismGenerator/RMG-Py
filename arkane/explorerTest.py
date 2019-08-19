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

"""
This module contains unit tests of the :mod:`arkane.explorer` module.
"""

import os
import unittest

from nose.plugins.attrib import attr

from arkane import Arkane
from arkane.explorer import ExplorerJob

################################################################################


@attr('functional')
class TestExplorerJob(unittest.TestCase):
    """
    Contains tests for ExplorerJob class execute method
    """

    @classmethod
    def setUpClass(cls):
        """A method that is run before each unit test in this class"""
        arkane = Arkane()

        cls.jobList = arkane.loadInputFile(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'methoxy_explore.py'))
        for job in cls.jobList:
            if not isinstance(job, ExplorerJob):
                job.execute(outputFile=None, plot=None)
            else:
                thermo_library, kinetics_library, species_list = arkane.getLibraries()
                job.execute(outputFile=None, plot=None, speciesList=species_list, thermoLibrary=thermo_library,
                            kineticsLibrary=kinetics_library)

        cls.thermoLibrary = thermo_library
        cls.kineticsLibrary = kinetics_library
        cls.explorerjob = cls.jobList[-1]
        cls.pdepjob = cls.jobList[-2]

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg
        rmgpy.data.rmg.database.kinetics = None

    def test_reactions(self):
        """
        test that the right number of reactions are in output network
        """
        self.assertEqual(len(self.explorerjob.networks[0].pathReactions), 6)

    def test_isomers(self):
        """
        test that the right number of isomers are in the output network
        """
        self.assertEqual(len(self.explorerjob.networks[0].isomers), 2)

    def test_job_rxns(self):
        """
        test that in this case all the reactions in the job
        ended up in the final network
        """
        for rxn in self.explorerjob.jobRxns:
            self.assertIn(rxn, self.explorerjob.networks[0].pathReactions)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
