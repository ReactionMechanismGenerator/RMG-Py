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
This module contains unit tests of the :mod:`arkane.pdep` module.
"""

import logging
import os
import shutil
import unittest

from nose.plugins.attrib import attr

from rmgpy import settings
from rmgpy.chemkin import read_reactions_block
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

    def test_pdep_job(self):
        """
        A general test for a PDep job in Arkane
        """
        self.tst1 = Arkane()
        self.tst1.input_file = self.input_file
        self.tst1.output_directory = self.directory
        self.tst1.verbose = logging.WARN
        self.tst1.plot = False
        self.tst1.job_list = []
        self.tst1.job_list = self.tst1.load_input_file(self.tst1.input_file)
        self.tst1.execute()

        job = self.tst1.job_list[0]
        self.assertEqual(job.Tmin.value_si, 300.0)
        self.assertEqual(job.minimum_grain_count, 100)
        self.assertFalse(job.rmgmode)
        self.assertTrue(job.active_j_rotor)
        self.assertEqual(job.network.path_reactions[0].label, 'acetylperoxy <=> hydroperoxylvinoxy')
        self.assertAlmostEqual(job.network.path_reactions[0].transition_state.tunneling.E0_TS.value_si, -24267.2)
        self.assertAlmostEqual(job.network.path_reactions[0].transition_state.tunneling.frequency.value_si, -1679.04)
        self.assertEqual(len(job.network.net_reactions[0].reactants[0].conformer.modes), 6)
        # self.assertEqual(self.tst1.frequencyScaleFactor, 0.947)

        # test that a network pdf was generated
        files = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        self.assertTrue(any(f == 'network.pdf' for f in files))

        # Test the generated network reaction
        dictionary = {'hydroperoxylvinoxy': Species().from_smiles('[CH2]C(=O)OO'),
                      'acetylperoxy': Species().from_smiles('CC(=O)O[O]')}
        with open(os.path.join(self.directory, 'chem.inp'), 'r') as chem:
            reaction_list = read_reactions_block(chem, dictionary)
        rxn = reaction_list[0]
        self.assertIsInstance(rxn.kinetics, Chebyshev)
        # Accept a delta of 0.2, which could result from numerical discrepancies
        # See RMG-Py #1682 on GitHub for discussion
        self.assertAlmostEqual(rxn.kinetics.get_rate_coefficient(1000.0, 1.0), 88.88253229631246, delta=0.2)

        files = [f for f in os.listdir(os.path.join(self.directory, 'sensitivity', ''))
                 if os.path.isfile(os.path.join(self.directory, 'sensitivity', f))]
        self.assertTrue(any('hydroperoxylvinoxy.pdf' in f for f in files))

        with open(os.path.join(self.directory, 'sensitivity', 'network1.txt'), 'r') as f:
            lines = f.readlines()
            for line in lines:
                if '1000.0' in line:
                    break
        sa_coeff = line.split()[-2]
        self.assertAlmostEqual(float(sa_coeff), -7.02e-07, delta=0.02e-6)

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

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
