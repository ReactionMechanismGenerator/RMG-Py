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

from rmgpy.quantity import Quantity
"""
Contains classes for termination criteria
"""
class TerminationTime:
    """
    Represent a time at which the simulation should be terminated. This class
    has one attribute: the termination `time` in seconds.
    """

    def __init__(self, time=(0.0, 's')):
        self.time = Quantity(time)


################################################################################

class TerminationConversion:
    """
    Represent a conversion at which the simulation should be terminated. This
    class has two attributes: the `species` to monitor and the fractional
    `conversion` at which to terminate.
    """

    def __init__(self, spec=None, conv=0.0):
        self.species = spec
        self.conversion = conv


class TerminationRateRatio:
    """
    Represent a fraction of the maximum characteristic rate of the simulation
    at which the simulation should be terminated.  This class has one attribute
    the ratio between the current and maximum characteristic rates at which
    to terminate
    """

    def __init__(self, ratio=0.01):
        self.ratio = ratio

class TerminationPolymerConversion:
    """
    Terminate when the DEFECT-ADJUSTED condensed polymer mass summed over ALL
    solver polymer pools (configured, spawned, homolysis-daughter,
    deprop-daughter and side-group feature pools, plus explicit-DP oligomer
    carriers) has dropped by the given fraction of its value at simulation
    initialization (r86):

        M_pool(t)   = max(0, mu1_p(t)*monomer_mw_g_mol_p
                             - mu0_p(t)*chain_mass_defect_g_mol_p)
        M_poly(t)   = sum over pools of M_pool(t)
        X_polymer(t) = 1 - M_poly(t)/M_poly(0)

    One attribute: the fractional `conversion` at which to terminate,
    validated strictly inside (0, 1) at this single chokepoint (the deck
    reader and the solver both construct through here).
    """

    def __init__(self, conv=0.0):
        conv = float(conv)
        if not 0.0 < conv < 1.0:
            raise ValueError(
                "terminationPolymerConversion must satisfy 0 < f < 1, got "
                f"{conv!r}.")
        self.conversion = conv
