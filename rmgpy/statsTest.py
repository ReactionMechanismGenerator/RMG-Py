#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This script contains unit tests of the :mod:`rmgpy.stats` module.
"""

import os
import os.path
import shutil
import unittest

from rmgpy.rmg.main import RMG, CoreEdgeReactionModel
from rmgpy.stats import ExecutionStatsWriter


################################################################################

class TestExecutionStatsWriter(unittest.TestCase):
    """
    Contains unit tests of the ExecutionStatsWriter.
    """

    def setUp(self):
        """
        Set up an RMG object
        """

        folder = os.path.join(os.getcwd(), 'rmgpy/output')
        if not os.path.isdir(folder):
            os.mkdir(folder)

        self.rmg = RMG(output_directory=folder)
        self.rmg.reaction_model = CoreEdgeReactionModel()

        self.rmg.save_everything()

    def test_save(self):
        """
        Tests if the statistics output file can be found.
        """

        folder = self.rmg.output_directory

        writer = ExecutionStatsWriter(folder)
        writer.update(self.rmg)

        statsfile = os.path.join(folder, 'statistics.xls')

        self.assertTrue(os.path.isfile(statsfile))

    def tearDown(self):
        shutil.rmtree(self.rmg.output_directory)
