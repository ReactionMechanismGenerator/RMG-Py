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
This module contains unit tests of the :mod:`arkane.kinetics` module.
"""

import unittest

from rmgpy.reaction import Reaction
from rmgpy.species import TransitionState

from arkane.kinetics import KineticsJob

################################################################################


class KineticsTest(unittest.TestCase):
    """
    Contains unit tests for the Arkane Kinetics module
    """

    def test_give_tlist_for_kineticsjob(self):
        """
        Ensures that the proper temperature ranges are set when Tlist is specified
        """
        rxn = Reaction(transition_state=TransitionState())
        t_list = [50.7, 100, 300, 800, 1255]
        kjob = KineticsJob(rxn, Tlist=(t_list, 'K'))
        self.assertEqual(min(t_list), kjob.Tmin.value_si)
        self.assertEqual(max(t_list), kjob.Tmax.value_si)
        self.assertEqual(len(t_list), kjob.Tcount)

    def test_give_temperature_range_for_kineticsjob(self):
        """
        Ensures that Tlist is set when a range of temperatures is specified
        """
        rxn = Reaction(transition_state=TransitionState())
        kjob = KineticsJob(rxn, Tmin=(50, 'K'), Tmax=(4000, 'K'), Tcount=5)
        self.assertEqual(5, len(kjob.Tlist.value_si))
        self.assertEqual(50, min(kjob.Tlist.value_si))
        self.assertAlmostEqual(4000, max(kjob.Tlist.value_si))
        inverse_tlist = 1 / kjob.Tlist.value_si
        self.assertAlmostEqual(inverse_tlist[1] - inverse_tlist[0],
                               inverse_tlist[-1] - inverse_tlist[-2],
                               msg='The points for fitting do not appear 1/T spaced. '
                                   'Obtained values of {0} and {1}'.format(inverse_tlist[1] - inverse_tlist[0],
                                                                           inverse_tlist[-1] - inverse_tlist[-2]))

    def test_get_tlist_for_kineticsjob(self):
        """
        Ensures that Tlist is set when no range is specified
        """
        rxn = Reaction(transition_state=TransitionState())
        kjob = KineticsJob(rxn)
        self.assertAlmostEqual(298, kjob.Tmin.value_si)
        self.assertAlmostEqual(2500, kjob.Tmax.value_si)
        self.assertEqual(50, kjob.Tcount)
        self.assertEqual(50, len(kjob.Tlist.value_si))
        self.assertAlmostEqual(298, min(kjob.Tlist.value_si))
        self.assertAlmostEqual(2500, max(kjob.Tlist.value_si))
        inverse_tlist = 1 / kjob.Tlist.value_si
        self.assertAlmostEqual(inverse_tlist[1] - inverse_tlist[0],
                               inverse_tlist[-1] - inverse_tlist[-2],
                               msg='The points for fitting do not appear 1/T spaced. '
                                   'Obtained values of {0} and {1}'.format(inverse_tlist[1] - inverse_tlist[0],
                                                                           inverse_tlist[-1] - inverse_tlist[-2]))

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
