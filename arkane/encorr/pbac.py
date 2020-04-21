#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This module provides methods for applying Petersson-type bond additivity
corrections (P-BAC) as described in:
Petersson et al., J. Chem. Phys. 1998, 109, 10570-10579
"""

import logging
import re

import arkane.encorr.data as data
from arkane.exceptions import BondAdditivityCorrectionError

################################################################################


def get_bac(model_chemistry, bonds):
    """
    Given the model_chemistry and a dictionary of bonds, return the
    total BAC (should be ADDED to energy).

    The dictionary of bonds should have the following form:

    bonds = {
        'C-H': bac1,
        'C-C': bac2,
        'C=C': bac3,
        ...
    }
    """

    # Get BAC parameters
    try:
        params = data.pbac[model_chemistry]
    except KeyError:
        raise BondAdditivityCorrectionError(
            'Missing Petersson-type BAC parameters for model chemistry {}'.format(model_chemistry)
        )

    # Sum corrections
    bac = 0.0
    for symbol, count in bonds.items():
        if symbol in params:
            bac += count * params[symbol]
        else:
            symbol_flipped = ''.join(re.findall('[a-zA-Z]+|[^a-zA-Z]+', symbol)[::-1])  # Check reversed symbol
            if symbol_flipped in params:
                bac += count * params[symbol_flipped]
            else:
                logging.warning('Ignored unknown bond type {}.'.format(symbol))

    return bac * 4184.0  # Convert kcal/mol to J/mol
