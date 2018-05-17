###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

cimport rmgpy.constants as constants
from rmgpy.species cimport Species, TransitionState
from rmgpy.molecule.molecule cimport Atom, Molecule
from rmgpy.molecule.element cimport Element
from rmgpy.kinetics.model cimport KineticsModel
from rmgpy.kinetics.arrhenius cimport Arrhenius

cimport numpy

################################################################################

cdef class Reaction:
    
    cdef public int index
    cdef public str label
    cdef public list reactants
    cdef public list products
    cdef public Species specificCollider
    cdef public bint reversible
    cdef public TransitionState transitionState
    cdef public KineticsModel kinetics
    cdef public Arrhenius network_kinetics
    cdef public bint duplicate
    cdef public float _degeneracy
    cdef public list pairs
    cdef public bint allow_pdep_route
    cdef public bint elementary_high_p
    cdef public str comment
    cdef public dict k_effective_cache
    
    cpdef bint isIsomerization(self)

    cpdef bint isDissociation(self)

    cpdef bint isAssociation(self)
    
    cpdef bint isUnimolecular(self)

    cpdef bint hasTemplate(self, list reactants, list products)
    
    cpdef bint matchesSpecies(self, list reactants, list products=?)

    cpdef bint isIsomorphic(self, Reaction other, bint eitherDirection=?, bint checkIdentical=?, bint checkOnlyLabel=?, bint checkTemplateRxnProducts=?)

    cpdef double getEnthalpyOfReaction(self, double T)

    cpdef double getEntropyOfReaction(self, double T)

    cpdef double getFreeEnergyOfReaction(self, double T)

    cpdef double getEquilibriumConstant(self, double T, str type=?)

    cpdef numpy.ndarray getEnthalpiesOfReaction(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEntropiesOfReaction(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getFreeEnergiesOfReaction(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEquilibriumConstants(self, numpy.ndarray Tlist, str type=?)

    cpdef int getStoichiometricCoefficient(self, Species spec)

    cpdef double getRateCoefficient(self, double T, double P=?)

    cpdef fixBarrierHeight(self, bint forcePositive=?)

    cpdef reverseThisArrheniusRate(self, Arrhenius kForward, str reverseUnits)

    cpdef generateReverseRateCoefficient(self, bint network_kinetics=?)

    cpdef numpy.ndarray calculateTSTRateCoefficients(self, numpy.ndarray Tlist)

    cpdef double calculateTSTRateCoefficient(self, double T) except -2

    cpdef bint canTST(self) except -2

    cpdef calculateMicrocanonicalRateCoefficient(self, numpy.ndarray Elist, numpy.ndarray Jlist, numpy.ndarray reacDensStates, numpy.ndarray prodDensStates=?, double T=?)

    cpdef bint isBalanced(self)
    
    cpdef generatePairs(self)
    
    cpdef copy(self)

    cpdef ensure_species(self, bint reactant_resonance=?, bint product_resonance=?)

cpdef bint _isomorphicSpeciesList(list list1, list list2, bint checkIdentical=?, bint checkOnlyLabel=?)