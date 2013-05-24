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

from rmgpy.statmech.mode cimport Mode
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Torsion(Mode):
    pass

################################################################################

cdef class HinderedRotor(Torsion):

    cdef public ScalarQuantity _inertia
    cdef public ArrayQuantity _fourier
    cdef public ScalarQuantity _barrier
    cdef public int symmetry
    cdef public bint semiclassical
    cdef public double frequency
    cdef public numpy.ndarray energies
    
    cpdef double getFrequency(self) except -1

    cpdef double getLevelEnergy(self, int J) except -1
    
    cpdef int getLevelDegeneracy(self, int J) except -1
    
    cpdef numpy.ndarray solveSchrodingerEquation(self, int Nbasis=?)
    
    cpdef numpy.ndarray getHamiltonian(self, int Nbasis)
    
    cdef double getRotationalConstantEnergy(self)
    
    cpdef double getPotential(self, double phi) except -100000000

    cpdef double getPartitionFunction(self, double T) except -1
        
    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=?)
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=?)

    cpdef fitFourierPotentialToData(self, numpy.ndarray angle, numpy.ndarray V)
    
    cpdef fitCosinePotentialToData(self, numpy.ndarray angle, numpy.ndarray V)
