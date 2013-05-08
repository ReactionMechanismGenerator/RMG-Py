# cython: embedsignature=True, cdivision=True

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

"""
This module contains functions for computing the microcanonical rate 
coefficient :math:`k(E)` using RRKM theory (the microcanonical analogue of
TST) and the inverse Laplace transform method. The former is more accurate,
but requires detailed information about the transition state and reactants.
"""

import logging
import numpy
cimport numpy
cimport cython
from libc.math cimport abs, exp, sqrt, cosh

cimport rmgpy.constants as constants
from rmgpy.kinetics.arrhenius cimport Arrhenius
from rmgpy.statmech.schrodinger import convolve

################################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
def calculateMicrocanonicalRateCoefficient(reaction,                                    
    numpy.ndarray[numpy.float64_t,ndim=1] Elist, 
    numpy.ndarray[numpy.int_t,ndim=1] Jlist, 
    numpy.ndarray[numpy.float64_t,ndim=2] reacDensStates, 
    numpy.ndarray[numpy.float64_t,ndim=2] prodDensStates=None,
    double T=0.0):
    """
    Calculate the microcanonical rate coefficient :math:`k(E)` for the reaction
    `reaction` at the energies `Elist` in J/mol. `reacDensStates` and 
    `prodDensStates` are the densities of states of the reactant and product
    configurations for this reaction. If the reaction is irreversible, only the
    reactant density of states is required; if the reaction is reversible, then
    both are required. This function will try to use the best method that it
    can based on the input data available:
    
    * If detailed information has been provided for the transition state (i.e.
      the molecular degrees of freedom), then RRKM theory will be used.
    
    * If the above is not possible but high-pressure limit kinetics
      :math:`k_\\infty(T)` have been provided, then the inverse Laplace 
      transform method will be used.

    The density of states for the product `prodDensStates` and the temperature
    of interest `T` in K can also be provided. For isomerization and association
    reactions `prodDensStates` is required; for dissociation reactions it is
    optional. The temperature is used if provided in the detailed balance
    expression to determine the reverse kinetics, and in certain cases in the
    inverse Laplace transform method.
    """        
    cdef int Ngrains, NJ, r, s
    cdef numpy.ndarray[numpy.float64_t,ndim=2] kf, kr
    cdef double C0inv
    cdef list modes
    cdef bint reactantStatesKnown, productStatesKnown, forward
    
    Ngrains = Elist.shape[0]
    NJ = Jlist.shape[0]
    kf = numpy.zeros((Ngrains,NJ))
    kr = numpy.zeros_like(kf)
    activeJRotor = Jlist is None
    activeKRotor = False
    
    reactantStatesKnown = reacDensStates is not None and reacDensStates.any()
    productStatesKnown = prodDensStates is not None and prodDensStates.any()
    forward = True
    
    C0inv = constants.R * T / 1e5
    
    if reaction.canTST():
        
        modes = reaction.transitionState.conformer.getActiveModes(activeJRotor=activeJRotor, activeKRotor=activeKRotor)
        
        # We've been provided with molecular degree of freedom data for the
        # transition state, so let's use the more accurate RRKM theory
        logging.debug('Calculating microcanonical rate coefficient using RRKM theory for {0}...'.format(reaction))
        if reactantStatesKnown and (reaction.isIsomerization() or reaction.isDissociation()):
            kf = applyRRKMTheory(reaction.transitionState, Elist, Jlist, reacDensStates)
            kf *= C0inv**(len(reaction.reactants) - 1)
            forward = True
        elif productStatesKnown and reaction.isAssociation():
            kr = applyRRKMTheory(reaction.transitionState, Elist, Jlist, prodDensStates)
            kr *= C0inv**(len(reaction.products) - 1)        
            forward = False
        else:
            raise Exception('Unable to compute k(E) values via RRKM theory for path reaction "{0}".'.format(reaction))           
        
    elif reaction.kinetics is not None:
        # We've been provided with high-pressure-limit rate coefficient data,
        # so let's use the less accurate inverse Laplace transform method
        logging.debug('Calculating microcanonical rate coefficient using ILT method for {0}...'.format(reaction))
        if reactantStatesKnown:
            kinetics = reaction.kinetics
            kf = applyInverseLaplaceTransformMethod(reaction.transitionState, kinetics, Elist, Jlist, reacDensStates, T)
            forward = True
        elif productStatesKnown:
            kinetics = reaction.generateReverseRateCoefficient()
            kr = applyInverseLaplaceTransformMethod(reaction.transitionState, kinetics, Elist, Jlist, prodDensStates, T)
            forward = False
        else:
            raise Exception('Unable to compute k(E) values via ILT method for path reaction "{0}".'.format(reaction))
    
    else:
        raise Exception('Unable to compute k(E) values for path reaction "{0}".'.format(reaction))

    # If the reaction is endothermic and barrierless, it is possible that the
    # forward k(E) will have a nonzero value at an energy where the product
    # density of states is zero (but the reactant density of states is not),
    # which violates detailed balance
    # To fix, we set the forward k(E) to zero wherever this is true
    # (This is correct within the accuracy of discretizing the energy grains)
    if reactantStatesKnown and productStatesKnown:
        for r in range(Ngrains):
            for s in range(NJ):
                if reacDensStates[r,s] != 0 and prodDensStates[r,s] != 0:
                    break
                kf[r,s] = 0
                kr[r,s] = 0
    
    # Get the reverse microcanonical rate coefficient if possible
    if forward and productStatesKnown:
        # We computed the forward rate coefficient above
        # Thus we need to compute the reverse rate coefficient here
        kr = numpy.zeros_like(kf)        
        for s in range(NJ):
            for r in range(Ngrains):
                if prodDensStates[r,s] != 0:
                    kr[r,s] = kf[r,s] * reacDensStates[r,s] / prodDensStates[r,s]
        kr *= C0inv**(len(reaction.products) - len(reaction.reactants))
    elif not forward and reactantStatesKnown:
        # We computed the reverse rate coefficient above
        # Thus we need to compute the forward rate coefficient here
        kf = numpy.zeros_like(kr)        
        for s in range(NJ):
            for r in range(Ngrains):
                if reacDensStates[r,s] != 0:
                    kf[r,s] = kr[r,s] * prodDensStates[r,s] / reacDensStates[r,s]
        kf *= C0inv**(len(reaction.reactants) - len(reaction.products))
     
    return kf, kr

@cython.boundscheck(False)
@cython.wraparound(False)
def applyRRKMTheory(transitionState,
                    numpy.ndarray[numpy.float64_t,ndim=1] Elist, 
                    numpy.ndarray[numpy.int_t,ndim=1] Jlist, 
                    numpy.ndarray[numpy.float64_t,ndim=2] densStates):
    """
    Calculate the microcanonical rate coefficient for a reaction using RRKM
    theory, where `transitionState` is the transition state of the reaction,
    `Elist` is the array of energies in J/mol at which to evaluate the
    microcanonial rate, and `densStates` is the density of states of the
    reactant.
    """
    cdef numpy.ndarray[numpy.float64_t,ndim=2] k0, k, sumStates
    cdef int Ngrains, NJ
    cdef bint activeJRotor
    cdef double dE, E0_TS
    cdef int r, s

    from rmgpy.pdep import Configuration
    
    Ngrains = Elist.shape[0]
    NJ = Jlist.shape[0]
    activeJRotor = (NJ == 1)
    k0 = numpy.zeros((Ngrains,NJ))
    k = numpy.zeros_like(k0)
    dE = Elist[1] - Elist[0]
    E0_TS = transitionState.conformer.E0.value_si
    
    conf = Configuration(transitionState)
    conf.calculateDensityOfStates(Elist - Elist[0], activeJRotor=activeJRotor)
    
    # Compute tunneling function
    kappa = transitionState.calculateTunnelingFunction(Elist)
    
    # Convolve with transition state density of states to get new transition
    # state sum of states that includes tunneling
    conf.sumStates = convolve(conf.densStates, kappa)
    conf.Elist += Elist[0] - E0_TS
    
    for r in range(Ngrains):
        if conf.sumStates[r] > 0:
            E0 = conf.Elist[r]
            break
    conf.Elist -= E0
    
    sumStates = conf.mapSumOfStates(Elist - E0, Jlist)
    
    # Generate k(E) using RRKM formula (with tunneling)
    dE /= constants.h * constants.Na        # J/mol -> s^-1
    for s in range(NJ):
        for r in range(Ngrains):
            if sumStates[r,s] > 0 and densStates[r,s] > 0:
                k[r,s] = sumStates[r,s] / densStates[r,s] * dE
            
    return k

@cython.boundscheck(False)
@cython.wraparound(False)
def applyInverseLaplaceTransformMethod(transitionState,
                                       Arrhenius kinetics,
                                       numpy.ndarray[numpy.float64_t,ndim=1] Elist, 
                                       numpy.ndarray[numpy.int_t,ndim=1] Jlist, 
                                       numpy.ndarray[numpy.float64_t,ndim=2] densStates, 
                                       double T=0.0):
    """
    Calculate the microcanonical rate coefficient for a reaction using the
    inverse Laplace transform method, where `kinetics` is the high pressure 
    limit rate coefficient, `E0` is the ground-state energy of the transition
    state, `Elist` is the array of energies in kJ/mol at which to evaluate the
    microcanonial rate, and `densStates` is the density of states of the
    reactant. The temperature `T` in K is not required, and is only used when
    the temperature exponent of the Arrhenius expression is negative (for which
    the inverse transform is undefined).
    """
    cdef numpy.ndarray[numpy.float64_t,ndim=2] k
    cdef numpy.ndarray[numpy.float64_t,ndim=1] phi0, phi
    cdef int Ngrains, NJ
    cdef bint activeJRotor
    cdef double dE, R, A, n, Ea, m0, rem, E0, num, E, n_crit
    cdef int r, s, m

    Ngrains = Elist.shape[0]
    NJ = Jlist.shape[0]
    k = numpy.zeros((Ngrains,NJ))
    dE = Elist[1] - Elist[0]
    R = constants.R
    E0 = transitionState.conformer.E0.value_si
    
    n = kinetics._n.value_si
    A = kinetics._A.value_si / (kinetics._T0.value_si**n)
    Ea = kinetics._Ea.value_si

    if isinstance(kinetics, Arrhenius) and (T != 0.0 or (Ea >= 0 and n >= 0)):
        
        # The inverse Laplace transform is not defined for Ea < 0 or n < 0
        # In these cases we move the offending portion into the preexponential
        # at the temperature of interest
        # This is an approximation, but it's not worth a more robust procedure
        
        # Including the T^n piece explicitly also has numerical difficulties
        # for small positive n; it turns out that using this approximation is
        # actually more accurate than trying to handle the T^n piece "properly"
        # For now the implementation is to use this approximation for all n
        # below some critical value, which is purposely placed a bit above zero
        n_crit = 0.25
        
        if Ea < 0:
            A *= exp(-Ea / R / T)
            Ea = 0.0
        if n < n_crit:
            A *= T**n
            n = 0.0

        if n < n_crit:
            # Determine the microcanonical rate directly
            m0, rem = divmod(Ea, dE)
            m = int(m0)
            if rem == 0:
                for s in range(NJ):
                    for r in range(m, Ngrains):
                        if Elist[r] > E0 and densStates[r,s] != 0:
                            k[r,s] = A * densStates[r-m,s] / densStates[r,s]
            else:
                for s in range(NJ):
                    for r in range(m+1, Ngrains):
                        if Elist[r] > E0 and densStates[r,s] != 0 and abs(densStates[r-m,s]) > 1e-12 and abs(densStates[r-m-1,s]) > 1e-12:
                            num = densStates[r-m,s] * (densStates[r-m-1,s] / densStates[r-m,s]) ** (-rem / (Elist[r-m-1] - Elist[r-m]))
                            k[r,s] = A * num / densStates[r,s]
                    
        elif n >= n_crit:
            import scipy.special
            # Evaluate the inverse Laplace transform of the T**n exp(-Ea/RT) piece, which only
            # exists for n >= 0
            phi0 = numpy.zeros(Ngrains, numpy.float64)
            for r in range(Ngrains):
                E = Elist[r] - Elist[0] - Ea
                if E > 0:
                    phi0[r] = (E/R)**(n-1.0)
            phi0 = phi0 * (dE / R) / scipy.special.gamma(n)
            # Evaluate the convolution
            for s in range(NJ):
                phi = convolve(phi0, densStates[:,s])
                # Apply to determine the microcanonical rate
                for r in range(Ngrains):
                    if densStates[r,s] != 0:
                        k[r,s] = A * phi[r] / densStates[r,s]
                            
    else:
        raise Exception('Unable to use inverse Laplace transform method for non-Arrhenius kinetics or for n < 0.')
    
    return k
