# -*- coding: utf-8 -*-

################################################################################
#
#   PyDAS - A Python wrapper for several differential algebraic system solvers
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

import numpy
cimport numpy

cdef class DASSL:
	
	cdef public int maxOrder
	cdef public object tstop
	cdef public double initialStep
	cdef public double maximumStep
	cdef public object bandwidths
	cdef public bint nonnegative
	
	cdef public double t
	cdef public numpy.ndarray y
	cdef public numpy.ndarray dydt
	
	cdef numpy.ndarray info
	cdef numpy.ndarray atol
	cdef numpy.ndarray rtol
	cdef numpy.ndarray rwork
	cdef numpy.ndarray iwork
	cdef numpy.ndarray rpar
	cdef numpy.ndarray ipar
	cdef int idid
	
	cpdef initialize(self, double t0, y0, dydt0=?, atol=?, rtol=?)
	
	cpdef advance(self, double tout)
		
	cpdef step(self, double tout)
		
	cdef solve(self, double tout)
	