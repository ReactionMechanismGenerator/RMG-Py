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
from libc.math cimport exp, log, sqrt, log10
from rmgpy.constants import R
from rmgpy.quantity import Quantity

"""
This module contains information related to kinetic uncertainties
"""

cdef class RateUncertainty(object):
    """
    Class for storing kinetic uncertainty information
    this provides uncertainty estimates from two different models:
    Constant uncertainty factor (self.f)
    and
    Energy uncertainty factor (self.getEnergyUncertaintyFactor(T))
    dE is the assumed 2 sigma error in the activation energy of the rate
    or alternatively f is the assumed 2 sigma log-error in the prefactor of the rate
    """

    def __init__(self,Tref,f=None,dE=None):
        self.Tref = Tref
        if dE:
            self.dE = dE
            self.f = exp(dE/(R*Tref))
        elif f:
            self.f = f
            self.dE = log(f)*R*Tref
        else:
            raise ValueError("Need Tref and either f or dE as input for RateUncertainty")

    property f:
      def __get__(self):
          return self._f
      def __set__(self, value):
          self._f = value

    property dE:
      def __get__(self):
          return self._dE
      def __set__(self, value):
          self._dE = value

    property Tref:
      def __get__(self):
          return self._Tref
      def __set__(self, value):
          self._Tref = value

    cpdef double getEnergyUncertaintyFactor(self,double T):
        return exp(self._dE/(R*T))

    cpdef double getConstantUncertaintyFactor(self):
        return self._f



rank_accuracy_map ={1:(0.0,'kcal/mol'),
                  2:(0.5,'kcal/mol'),
                  3:(1.0,'kcal/mol'),
                  4:(1.5,'kcal/mol'),
                  5:(2.5,'kcal/mol'),
                  6:(3.5,'kcal/mol'),
                  7:(4.0,'kcal/mol'),
                  8:(5.0,'kcal/mol'),
                  9:(14.0,'kcal/mol'),
                  10:(14.0,'kcal/mol'),
                  None:(14.0,'kcal/mol'),
                  0:(14.0,'kcal/mol'),
                  '':(14.0,'kcal/mol'),
                  11:(14.0,'kcal/mol'),
                  }
rank_accuracy_map = {i:Quantity(rank_accuracy_map[i]) for i in rank_accuracy_map.keys()}
