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
