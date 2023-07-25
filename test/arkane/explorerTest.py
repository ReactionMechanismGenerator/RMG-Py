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
This module contains unit tests of the :mod:`arkane.explorer` module.
"""

import os


import pytest

from arkane import Arkane
from arkane.explorer import ExplorerJob
import rmgpy.data.rmg

ADMONITION = (
    "This functional test has been known to fail for apparently no reason if run alongside "
    "other tests (when run alone it passes, but with others it fails). This is likely due "
    "to a global state issue that was implicitly handled differently by the old nose testing "
    "framework. \nThe best solution for this problem is to remove the abuse of the global "
    "state in Arkane and RMG-Py rather than try and fix this single functional test."
)


@pytest.mark.skip(reason=ADMONITION)
@pytest.mark.functional
class TestExplorerJob:
    """
    Contains tests for ExplorerJob class execute method
    """

    def test_explorer(self):
        """A method that is run before each unit test in this class"""
        arkane = Arkane()
        job_list = arkane.load_input_file(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "arkane", "data", "methoxy_explore.py")
        )
        for job in job_list:
            if not isinstance(job, ExplorerJob):
                job.execute(output_file=None, plot=None)
            else:
                thermo_library, kinetics_library, species_list = arkane.get_libraries()
                job.execute(
                    output_file=None,
                    plot=None,
                    thermo_library=thermo_library,
                    kinetics_library=kinetics_library,
                )
                explorer_job = job_list[-1]
                assert len(explorer_job.networks[0].path_reactions) in [6, 7]
                assert len(explorer_job.networks[0].isomers) == 2
                for rxn in explorer_job.job_rxns:
                    assert rxn in explorer_job.networks[0].path_reactions

        rmgpy.data.rmg.database.kinetics = None
