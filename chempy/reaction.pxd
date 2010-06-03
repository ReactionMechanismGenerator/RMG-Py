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

from species cimport Species
from kinetics cimport KineticsModel

################################################################################

cdef class Reaction:
	
	cdef public int index
	cdef public list reactants
	cdef public list products
	cdef public KineticsModel kinetics

	cpdef double getEnthalpyOfReaction(self, double T)

	cpdef double getEntropyOfReaction(self, double T)

	cpdef double getFreeEnergyOfReaction(self, double T)

	cpdef double getEquilibriumConstant(self, double T, str type=?)

	cpdef double getRateConstant(self, double T, double P)

	cpdef int getStoichiometricCoefficient(self, Species spec)

	cpdef double getRate(self, double T, double P, dict conc)
