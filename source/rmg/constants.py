#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""Contains a number of constants to be used throughout RMG."""

import quantities as pq

################################################################################

# Global variables: important directories
outputDir = '.'
scratchDir = '.'
libraryDir = '.'

# Global variables: options
drawMolecules = False
generatePlots = False

################################################################################

# RMG uses SI units throughout
pq.set_default_units('si')

if __name__ == "__main__":
	print 'Constants available:'
	print ''

# putting a comment with a colon just before the thing is defined
# makes it show up in the documentation with autodoc. Clever hey?
#: The Avogadro constant
Na = pq.constants.Avogadro_constant.simplified
if __name__ == "__main__":
	print 'Avogadro constant: Na = %s' % (Na)
Na = float(Na)

#: The Boltzmann constant
kB = pq.constants.Boltzmann_constant.simplified
if __name__ == "__main__":
	print 'Boltzmann constant: kB = %s' % (kB)
kB = float(kB)

#: The gas law constant
R = pq.constants.R.simplified
if __name__ == "__main__":
	print 'Gas law constant: R = %s' % (R)
R = float(R)
"""The molar gas constant, in J/mol/s"""
# you can also put a docstring straight after the definition 
# and autodoc will find it

#: The Planck constant
h = pq.constants.Planck_constant.simplified
if __name__ == "__main__":
	print 'Planck constant: h = %s' % (h)
h = float(h)

#: The speed of light in a vacuum
c = pq.constants.natural_unit_of_velocity.simplified
if __name__ == "__main__":
	print 'Speed of light in a vacuum: c = %s' % (c)
c = float(c)

if __name__ == "__main__":
	print ''
