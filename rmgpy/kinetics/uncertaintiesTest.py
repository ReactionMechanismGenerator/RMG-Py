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

"""
This script contains unit tests of the :mod:`rmgpy.kinetics.arrhenius` module.
"""

import unittest
import numpy as np
from rmgpy.kinetics.uncertainties import RateUncertainty
from rmgpy.constants import R

class TestUncertainties(unittest.TestCase):
    """
    Contains unit tests for the RateUncertainty class
    """
    def testFactorConstructed(self):
        """
        Test RateUncertainty constructed from factor
        """
        unc = RateUncertainty(Tref=1000.0,f=3.0)
        fp = unc.getEnergyUncertaintyFactor(1000.0)
        self.assertAlmostEqual(fp,3.0)

    def testEnergyConstructed(self):
        """
        Test RateUncertainty constructed from energy
        """
        unc = RateUncertainty(Tref=1000.0,dE=1.0)
        self.assertAlmostEqual(unc.f,np.exp(1.0/(R*1000.0)))
        fp = unc.getEnergyUncertaintyFactor(700.0)
        self.assertAlmostEqual(fp,unc.f*np.exp(1.0/(R*700.0)-1.0/(R*1000.0)))
