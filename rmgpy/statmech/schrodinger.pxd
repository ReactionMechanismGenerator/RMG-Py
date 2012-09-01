################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

cimport numpy

################################################################################

cpdef double getPartitionFunction(double T, energy, degeneracy=?, int n0=?, int nmax=?, double tol=?) except -1

cpdef double getHeatCapacity(double T, energy, degeneracy=?, int n0=?, int nmax=?, double tol=?) except -100000000

cpdef double getEnthalpy(double T, energy, degeneracy=?, int n0=?, int nmax=?, double tol=?) except 100000000

cpdef double getEntropy(double T, energy, degeneracy=?, int n0=?, int nmax=?, double tol=?) except -100000000

cpdef numpy.ndarray getSumOfStates(numpy.ndarray Elist, energy, degeneracy=?, int n0=?, numpy.ndarray sumStates0=?)

cpdef numpy.ndarray getDensityOfStates(numpy.ndarray Elist, energy, degeneracy=?, int n0=?, numpy.ndarray densStates0=?)
