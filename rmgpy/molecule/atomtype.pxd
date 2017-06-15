################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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

cdef class AtomType:

    cdef public str label
    cdef public list generic
    cdef public list specific

    cdef public list incrementBond
    cdef public list decrementBond
    cdef public list formBond
    cdef public list breakBond
    cdef public list incrementRadical
    cdef public list decrementRadical
    cdef public list incrementLonePair
    cdef public list decrementLonePair

    cdef public list single
    cdef public list allDouble
    cdef public list rDouble
    cdef public list oDouble
    cdef public list sDouble
    cdef public list triple
    cdef public list benzene
    cdef public list lonePairs
    cdef public list charge

    cpdef bint isSpecificCaseOf(self, AtomType other)

    cpdef bint equivalent(self, AtomType other)

    cpdef list getFeatures(self)

cpdef list getFeatures(atom, dict bonds)

cpdef AtomType getAtomType(atom, dict bonds)
