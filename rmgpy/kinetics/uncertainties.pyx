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
This module contains information related to kinetic uncertainties
"""

from __future__ import division

import numpy as np
cimport numpy as np

from rmgpy.quantity import Quantity


cdef class RateUncertainty(object):
    """
    Class for storing kinetic uncertainty information describes identically a normal distribution 
    on Log(k) and a lognormal distribution on k by mu and var at a single temperature Tref also 
    includes potentially useful uncertainty treatment variables
    N = number of samples used to generate the distribution 
    correlation = label identifying source of estimate for correlated error treatment
    Note that correlated errors are expected to be associated only with mu (the bias of the distribution)
    """

    def __init__(self, mu, var, Tref, N=None, correlation=None):
        self.Tref = Tref
        self.correlation = correlation
        self.mu = mu
        self.var = var
        self.N = N

    def __repr__(self):
        s = "RateUncertainty(mu={mu}, var={var}, Tref={Tref},".format(mu=self.mu, var=self.var, Tref=self.Tref)
        if self.N is not None:
            s += " N={0!r},".format(self.N)
        if self.correlation is not None:
            s += " correlation={0!r},".format(self.correlation)
        s += ")"
        return s

    cpdef double get_expected_log_uncertainty(self):
        """
        The expected uncertainty in Log(k) at Tref
        """
        return np.sqrt(self.var * 2.0 / np.pi) + abs(self.mu)

rank_accuracy_map = {1: (0.0, 'kcal/mol'),
                     2: (0.5, 'kcal/mol'),
                     3: (1.0, 'kcal/mol'),
                     4: (1.5, 'kcal/mol'),
                     5: (2.5, 'kcal/mol'),
                     6: (3.5, 'kcal/mol'),
                     7: (4.0, 'kcal/mol'),
                     8: (5.0, 'kcal/mol'),
                     9: (14.0, 'kcal/mol'),
                     10: (14.0, 'kcal/mol'),
                     None: (14.0, 'kcal/mol'),
                     0: (14.0, 'kcal/mol'),
                     '': (14.0, 'kcal/mol'),
                     11: (14.0, 'kcal/mol'),
                     }
rank_accuracy_map = {i: Quantity(rank_accuracy_map[i]) for i in rank_accuracy_map.keys()}
