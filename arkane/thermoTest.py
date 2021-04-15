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
This module contains unit tests of the :mod:`arkane.thermo` module.
"""

import os
import unittest

from rmgpy.species import Species

from arkane.ess.gaussian import GaussianLog
from arkane.thermo import ThermoJob

################################################################################


class TestThermo(unittest.TestCase):
    """
    Contains unit tests of the ThermoJob class.
    """

    @classmethod
    def setUp(cls):
        """A method that is run before each unit test in this class"""
        spc = Species().from_smiles('CCO')
        log = GaussianLog(os.path.join(os.path.dirname(__file__), 'data', 'gaussian', 'ethylene.log'))
        spc.conformer = log.load_conformer()[0]
        coords, numbers, masses = log.load_geometry()
        spc.conformer.coordinates = coords, 'angstroms'
        spc.conformer.number = numbers
        spc.conformer.mass = masses, 'amu'
        cls.thermo_job = ThermoJob(species=spc, thermo_class='NASA')

    def test_element_count_from_conformer(self):
        """Test Getting an element count dictionary from the species.conformer attribute"""
        element_count = self.thermo_job.element_count_from_conformer()
        self.assertEqual(element_count, {'H': 4, 'C': 2})

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
