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

The constants available are listed below.

"""

import quantities as pq
import math
import cython

################################################################################

# RMG uses SI units throughout
pq.set_default_units('si')

value = pq.constants.Avogadro_constant.simplified
#: The Avogadro constant
Na = float(value)

value = pq.constants.Boltzmann_constant.simplified
#: The Boltzmann constant
kB = float(value)

value = pq.constants.R.simplified
#: The gas law constant
R = float(value)

value = pq.constants.Planck_constant.simplified
#: The Planck constant
h = float(value)

value = pq.constants.natural_unit_of_velocity.simplified
#: The speed of light in a vacuum
c = float(value)

#: pi
pi = float(math.pi)
