# -*- coding: utf-8 -*-

################################################################################
#
#   PyDAS - A Python wrapper for several differential algebraic system solvers
#
#   Copyright (c) 2010-2014 by Joshua W. Allen (jwallen@mit.edu) 
#   and Connie W. Gao (connieg@mit.edu)
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


import numpy as np
cimport numpy as np

cdef class DASPK:

    cdef public int maxOrder
    cdef public object tstop
    cdef public double initialStep
    cdef public double maximumStep
    cdef public object bandwidths
    cdef public bint nonnegative
    cdef public bint sensitivity
    cdef public int sensmethod
    
    cdef public double t
    cdef public np.ndarray y
    cdef public np.ndarray dydt
    cdef public np.ndarray senpar
    
    cdef np.ndarray info
    cdef np.ndarray atol
    cdef np.ndarray rtol
    cdef np.ndarray rwork
    cdef np.ndarray iwork
    cdef np.ndarray rpar
    cdef np.ndarray ipar
    cdef int idid
    
    cpdef initialize(self, double t0, np.ndarray y0, np.ndarray dydt0=?, np.ndarray senpar=?, atol=?, rtol=?)
    
    cpdef advance(self, double tout)
    
    cpdef step(self, double tout)
    
    cdef solve(self, double tout)