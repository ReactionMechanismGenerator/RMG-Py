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

The constants available are:

=============== ================================================================
Variable        Description
=============== ================================================================
`Na`            Avogadro constant
`kB`            Boltzmann constant
`R`             Gas law constant
`h`             Planck constant
`c`             Speed of light in vacuum
=============== ================================================================

"""

import quantities as pq
import math
import cython

################################################################################

# RMG uses SI units throughout
pq.set_default_units('si')

if __name__ == "__main__":
	print 'Constants available:'
	print ''

#: The Avogadro constant
value = pq.constants.Avogadro_constant.simplified
if __name__ == "__main__":
	print 'Avogadro constant: Na = %s' % (value)
Na = float(value)

#: The Boltzmann constant
value = pq.constants.Boltzmann_constant.simplified
if __name__ == "__main__":
	print 'Boltzmann constant: kB = %s' % (value)
kB = float(value)

#: The gas law constant
value = pq.constants.R.simplified
if __name__ == "__main__":
	print 'Gas law constant: R = %s' % (value)
R = float(value)

#: The Planck constant
value = pq.constants.Planck_constant.simplified
if __name__ == "__main__":
	print 'Planck constant: h = %s' % (value)
h = float(value)

#: The speed of light in a vacuum
value = pq.constants.natural_unit_of_velocity.simplified
if __name__ == "__main__":
	print 'Speed of light in a vacuum: c = %s' % (value)
c = float(value)

#: pi
pi = float(math.pi)

################################################################################

if __name__ == "__main__":
	print ''
