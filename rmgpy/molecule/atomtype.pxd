###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

cdef class AtomType:

    cdef public str label
    cdef public list generic
    cdef public list specific

    cdef public list increment_bond
    cdef public list decrement_bond
    cdef public list form_bond
    cdef public list break_bond
    cdef public list increment_radical
    cdef public list decrement_radical
    cdef public list increment_lone_pair
    cdef public list decrement_lone_pair
    cdef public list increment_charge
    cdef public list decrement_charge

    cdef public list single
    cdef public list all_double
    cdef public list r_double
    cdef public list o_double
    cdef public list s_double
    cdef public list triple
    cdef public list quadruple
    cdef public list benzene
    cdef public list lone_pairs
    cdef public list charge

    cpdef bint is_specific_case_of(self, AtomType other)

    cpdef bint equivalent(self, AtomType other)

    cpdef list get_features(self)

cpdef list get_features(atom, dict bonds)

cpdef AtomType get_atomtype(atom, dict bonds)
