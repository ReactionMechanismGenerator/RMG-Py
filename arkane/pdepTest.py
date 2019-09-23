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
This module contains unit tests of the :mod:`arkane.pdep` module.
"""

import logging
import os
import shutil
import unittest

from nose.plugins.attrib import attr

from rmgpy import settings
from rmgpy.chemkin import readReactionsBlock
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.species import Species

from arkane.main import Arkane

################################################################################


@attr('functional')
class ArkaneTest(unittest.TestCase):
    """
    Contains unit tests for the sensitivity module in Arkane
    """

    @classmethod
    def setUp(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.directory = os.path.join(settings['test_data.directory'], 'arkane', 'tst1', '')
        cls.input_file = os.path.join(cls.directory, 'pdep_sa.py')

        # clean working folder from all previous test output
        dirs = [d for d in os.listdir(cls.directory) if not os.path.isfile(os.path.join(cls.directory, d))]
        for d in dirs:
            shutil.rmtree(os.path.join(settings['test_data.directory'], 'arkane', 'tst1', d, ''))
        files = [f for f in os.listdir(cls.directory) if os.path.isfile(os.path.join(cls.directory, f))]
        for f in files:
            if 'pdep_sa' not in f:
                os.remove(os.path.join(settings['test_data.directory'], 'arkane', 'tst1', f))

    def testPDepJob(self):
        """
        A general test for a PDep job in Arkane
        """
        self.tst1 = Arkane()
        self.tst1.inputFile = self.input_file
        self.tst1.outputDirectory = self.directory
        self.tst1.verbose = logging.WARN
        self.tst1.plot = False
        self.tst1.jobList = []
        self.tst1.jobList = self.tst1.loadInputFile(self.tst1.inputFile)
        self.tst1.execute()

        job = self.tst1.jobList[0]
        self.assertEquals(job.Tmin.value_si, 300.0)
        self.assertEquals(job.minimumGrainCount, 100)
        self.assertFalse(job.rmgmode)
        self.assertTrue(job.activeJRotor)
        self.assertEquals(job.network.pathReactions[0].label, 'acetylperoxy <=> hydroperoxylvinoxy')
        self.assertAlmostEquals(job.network.pathReactions[0].transition_state.tunneling.E0_TS.value_si, -24267.2)
        self.assertAlmostEquals(job.network.pathReactions[0].transition_state.tunneling.frequency.value_si, -1679.04)
        self.assertEquals(len(job.network.netReactions[0].reactants[0].conformer.modes), 6)
        # self.assertEquals(self.tst1.frequencyScaleFactor, 0.947)

        # test that a network pdf was generated
        files = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        self.assertTrue(any(f == 'network.pdf' for f in files))

        # Test the generated network reaction
        dictionary = {'hydroperoxylvinoxy': Species().fromSMILES('[CH2]C(=O)OO'),
                      'acetylperoxy': Species().fromSMILES('CC(=O)O[O]')}
        with open(os.path.join(self.directory, 'chem.inp'), 'r') as chem:
            reaction_list = readReactionsBlock(chem, dictionary)
        rxn = reaction_list[0]
        self.assertIsInstance(rxn.kinetics, Chebyshev)
        self.assertAlmostEquals(rxn.kinetics.getRateCoefficient(1000.0, 1.0), 88.88253229631246)

        files = [f for f in os.listdir(os.path.join(self.directory, 'sensitivity', ''))
                 if os.path.isfile(os.path.join(self.directory, 'sensitivity', f))]
        self.assertTrue(any('hydroperoxylvinoxy.pdf' in f for f in files))

        with open(os.path.join(self.directory, 'sensitivity', 'network1.txt'), 'r') as f:
            lines = f.readlines()
            for line in lines:
                if '1000.0' in line:
                    break
        sa_coeff = line.split()[-2]
        self.assertEquals(float(sa_coeff), -8.23e-6)

    @classmethod
    def tearDown(cls):
        """A function that is run ONCE after all unit tests in this class."""
        cls.directory = os.path.join(settings['test_data.directory'], 'arkane', 'tst1', '')
        cls.input_file = os.path.join(cls.directory, 'pdep_sa.py')

        # clean working folder from all previous test output
        dirs = [d for d in os.listdir(cls.directory) if not os.path.isfile(os.path.join(cls.directory, d))]
        for d in dirs:
            shutil.rmtree(os.path.join(settings['test_data.directory'], 'arkane', 'tst1', d, ''))
        files = [f for f in os.listdir(cls.directory) if os.path.isfile(os.path.join(cls.directory, f))]
        for f in files:
            if 'pdep_sa' not in f:
                os.remove(os.path.join(settings['test_data.directory'], 'arkane', 'tst1', f))
