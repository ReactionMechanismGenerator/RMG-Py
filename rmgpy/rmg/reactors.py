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
Contains classes for building RMS simulations
Note since rms currently can't be imported inside nosetests
This most of this can't be tested using nosetests
"""
import numpy as np
import sys
import logging
import itertools
try:
    from pyrms import rms
    from diffeqpy import de
except:
    pass

from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.thermo.thermodata import ThermoData
from rmgpy.kinetics.arrhenius import Arrhenius, ArrheniusEP, ArrheniusBM, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius
from rmgpy.kinetics.kineticsdata import KineticsData
from rmgpy.kinetics.falloff import Troe, ThirdBody, Lindemann
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.data.solvation import SolventData
from rmgpy.kinetics.surface import StickingCoefficient
from rmgpy.solver.termination import TerminationTime, TerminationConversion, TerminationRateRatio
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.depository import DepositoryReaction
