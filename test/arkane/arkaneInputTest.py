#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains unit tests of the :mod:`arkane.input` module.
"""

import os

from rmgpy.species import Species, TransitionState

from arkane import input
from arkane.input import job_list
from arkane.modelchem import LevelOfTheory
from arkane.statmech import StatMechJob
import pytest


ADMONITION = (
    "This unit test fails in the new pytest framework despite other similar tests passing. "
    "This is likely due to a global state issue that was implicitly handled differently by the old nose testing "
    "framework. \nThe best solution for this problem is to remove the abuse of the global "
    "state in Arkane and RMG-Py rather than try and fix this single test."
)


@pytest.mark.skip(reason=ADMONITION)
class TestArkaneInput:
    """
    Contains unit tests for loading and processing Arkane input files.
    """

    def setup_class(cls):
        """Preparation for all unit tests in this class."""
        cls.directory = os.path.join(os.path.dirname(os.path.dirname(__file__)), "..", "examples", "arkane")
        cls.level_of_theory = LevelOfTheory("cbs-qb3")
        cls.frequencyScaleFactor = 0.99
        cls.useHinderedRotors = False
        cls.useBondCorrections = True

    def test_species(self):
        """Test loading of species input file."""
        spec = input.species("C2H4", os.path.join(self.directory, "species", "C2H4", "ethene.py"))
        assert isinstance(spec, Species)
        assert len(spec.molecule) == 0
        # statmech job test
        job = job_list[-1]
        assert isinstance(job, StatMechJob)
        job.level_of_theory = self.level_of_theory
        job.frequencyScaleFactor = self.frequencyScaleFactor
        job.includeHinderedRotors = self.useHinderedRotors
        job.applyBondEnergyCorrections = self.useBondCorrections
        job.load()
        assert isinstance(job.species.props["element_counts"], dict)
        assert job.species.props["element_counts"]["C"] == 2
        assert job.species.props["element_counts"]["H"] == 4
        # thermo job tests
        input.thermo("C2H4", "NASA")
        job = job_list[-1]
        filepath = os.path.join(self.directory, "reactions", "H+C2H4=C2H5")
        job.execute(output_directory=filepath)
        assert os.path.isfile(os.path.join(filepath, "output.py"))
        assert os.path.isfile(os.path.join(filepath, "chem.inp"))
        os.remove(os.path.join(filepath, "output.py"))
        os.remove(os.path.join(filepath, "chem.inp"))

    def test_transition_state(self):
        """Test loading of transition state input file."""
        ts = input.transitionState("TS", os.path.join(self.directory, "reactions", "H+C2H4=C2H5", "TS.py"))
        assert isinstance(ts, TransitionState)
        # stat mech job tests
        job = job_list[-1]
        assert isinstance(job, StatMechJob)
        job.level_of_theory = self.level_of_theory
        job.frequencyScaleFactor = self.frequencyScaleFactor
        job.includeHinderedRotors = self.useHinderedRotors
        job.applyBondEnergyCorrections = self.useBondCorrections
        job.load()
