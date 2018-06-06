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

import os
import numpy as np
import logging
import shutil
from copy import deepcopy

from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.data.rmg import getDB

class ExplorerJob(object):
    def __init__(self, source, pdepjob, explore_tol, energy_tol=np.inf, flux_tol=0.0):
        self.source = source
        self.explore_tol = explore_tol
        self.energy_tol = energy_tol
        self.flux_tol = flux_tol
        
        self.pdepjob = pdepjob
        
        if not hasattr(self.pdepjob,'outputFile'):
            self.pdepjob.outputFile = None

    def copy(self):
        """
        Return a copy of the explorer job.
        """
        return ExplorerJob(
            source=deepcopy(self.source),
            pdepjob=self.pdepjob,
            explore_tol=self.explore_tol,
            energy_tol=self.energy_tol,
            flux_tol=self.flux_tol
        )
        

        
        
        
        