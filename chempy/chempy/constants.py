#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains a number of physical constants to be made available
throughout ChemPy. ChemPy uses SI units throughout; accordingly, all of the
constants in this module are stored in combinations of meters, seconds,
kilograms, moles, etc.

The constants available are listed below. All values were taken from
`NIST <http://physics.nist.gov/cuu/Constants/index.html>`_

"""

import math
import cython
import numpy

################################################################################

#: The Avogadro constant
Na = 6.02214179e23

#: The Boltzmann constant
kB = 1.3806504e-23

#: The gas law constant
R = 8.314472

#: The Planck constant
h = 6.62606896e-34

#: The speed of light in a vacuum
c = 299792458

#: pi
pi = float(math.pi)

################################################################################

def processQuantity(quantity):
    """
    Convert a numeric quantity or set of quantities in the input file to a
    consistent set of units (SI). The parameter `quantity` is usually a 2-tuple
    or 2-list, with the first element a single number or a list or tuple of
    numbers, and the second element is the units associated with the number(s)
    in the first element. If `quantity` is a number or a list or tuple of
    numbers, then they are assumed to already be in SI units.

    .. note::

        The ``quantities`` package is used to convert your numeric parameters
        into SI units, and therefore inherits all of the idiosyncracies from
        that package. In particular, the ``quantities`` package does *not*
        follow the SI convention that all units after the circumflex are in the
        denominator. For example, ``J/mol*K`` would be interpreted as
        ``(J/mol)*K`` rather than ``J/(mol*K)``. Thus we recommend using
        parentheses where necessary to make your intentions explicit.

    """
    import quantities

    # Do nothing if the parameter is invalid
    if quantity is None: return None, ''
    # If the parameter is a number or a numpy array, then immediately return it
    # (so we avoid the slow calls to quantities)
    if isinstance(quantity, float) or isinstance(quantity, int) or isinstance(quantity, numpy.ndarray):
        return quantity, ''

    if (isinstance(quantity, tuple) or isinstance(quantity, list)) and len(quantity) == 2:
        value, units = quantity
    else:
        value = quantity; units = ''

    # Get the output (SI) units corresponding to the input units
    factor = quantities.Quantity(1.0, units).simplified
    newUnits = str(factor.units).split()[1]
    factor = float(factor)

    if isinstance(value, tuple) or isinstance(value, list):
        return numpy.array([v * factor for v in value], numpy.float64), newUnits
    elif isinstance(value, numpy.ndarray):
        return value, newUnits
    else:
        return float(value) * factor, newUnits

