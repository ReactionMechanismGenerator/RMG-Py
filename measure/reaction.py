#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   MEASURE - Master Equation Automatic Solver for Unimolecular REactions
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

"""
Contains functions for working with reaction objects, particularly for
calculating microcanonical rate coefficients using various methods.
"""

import numpy
import logging

import chempy.constants as constants
import chempy.reaction
from chempy.kinetics import *
from chempy.states import convolve

################################################################################

class ReactionError(Exception): 
    """
    An exception raised when working with reactions causes exceptional behavior
    for any reason. Pass a string describing the cause of the exceptional 
    behavior.
    """
    pass

################################################################################

def calculateMicrocanonicalRateCoefficient(reaction, Elist, reacDensStates, prodDensStates=None, T=None):
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
      
    """
    
    kf = numpy.zeros_like(Elist)
    kr = numpy.zeros_like(Elist)
    
    if reaction.transitionState.states is not None:
        # We've been provided with molecular degree of freedom data for the
        # transition state, so let's use the more accurate RRKM theory
        logging.debug('Using RRKM theory for reaction "%s"' % reaction)
        kf = applyRRKMTheory(reaction.transitionState, Elist, reacDensStates)
    elif reaction.kinetics is not None:
        # We've been provided with high-pressure-limit rate coefficient data,
        # so let's use the less accurate inverse Laplace transform method
        logging.debug('Using ILT method for reaction "%s"' % reaction)
        kf = applyInverseLaplaceTransformMethod(reaction.kinetics, reaction.transitionState.E0, Elist, reacDensStates, T)
    
    # If the reaction is reversible, calculate the reverse microcanonical rate
    # using detailed balance
    if reaction.reversible:
        for r in range(len(Elist)):
            if prodDensStates[r] > 0: 
                kr[r] = kf[r] * reacDensStates[r] / prodDensStates[r]
    
    return kf, kr

################################################################################

def applyRRKMTheory(transitionState, Elist, densStates):
    """
    Calculate the microcanonical rate coefficient for a reaction using RRKM
    theory, where `transitionState` is the transition state of the reaction,
    `Elist` is the array of energies in J/mol at which to evaluate the
    microcanonial rate, and `densStates` is the density of states of the
    reactant.
    """
    
    k = numpy.zeros_like((Elist))
    Ngrains = len(Elist)
    
    # Calculate sum of states of transition state
    sumStates0 = transitionState.states.getSumOfStates(Elist)
    # Shift to common zero of energy
    r0 = int(round(transitionState.E0 / dE))
    sumStates[n+Nisom,r0:] = sumStates0[:-r0+len(densStates0)]
    
    # Generate k(E) using RRKM formula
    for r in range(len(Elist)):
        if densStates[r] > 0:
            k[r] = sumStates[r] / constants.h / densStates[r]
    
    return k

################################################################################

def applyInverseLaplaceTransformMethod(kinetics, E0, Elist, densStates, T=None):
    """
    Calculate the microcanonical rate coefficient for a reaction using the
    inverse Laplace transform method, where `kinetics` is the high pressure 
    limit rate coefficient, `E0` is the ground-state energy of the transition
    state, `Elist` is the array of energies in J/mol at which to evaluate the
    microcanonial rate, and `densStates` is the density of states of the
    reactant.
    """
    
    k = numpy.zeros_like((Elist))
    
    if isinstance(kinetics, ArrheniusModel) and (T is not None or (kinetics.Ea >= 0 and kinetics.n >= 0)):
        A = kinetics.A
        n = kinetics.n
        Ea = kinetics.Ea
        dE = Elist[1] - Elist[0]

        # The inverse Laplace transform is not defined for Ea < 0 or n < 0
        # In these cases we move the offending portion into the preexponential
        # at the temperature of interest
        # This is an approximation, but it's not worth a more robust procedure
        if Ea < 0:
            A *= math.exp(-Ea / constants.R / T)
            Ea = 0.0
        if n < 0:
            A *= T**n
            n = 0.0

        if n == 0:
            # Determine the microcanonical rate directly
            s = int(math.floor(Ea / dE))
            for r in range(len(Elist)):
                if Elist[r] > E0 and densStates[r] != 0:
                    k[r] = A * densStates[r - s] / densStates[r]
                    
        elif n > 0.0:
            import scipy.special
            # Evaluate the inverse Laplace transform of the T**n piece, which only
            # exists for n >= 0
            phi = numpy.zeros(len(Elist), numpy.float64)
            for i, E in enumerate(Elist):
                if E == 0.0:
                    phi[i] = 0.0
                else:
                    phi[i] = E**(n-1) / (constants.R**n * scipy.special.gamma(n))
            # Evaluate the convolution
            phi = convolve(phi, densStates, Elist)
            # Apply to determine the microcanonical rate
            s = int(math.floor(Ea / dE))
            for r in range(len(Elist)):
                if Elist[r] > E0 and densStates[r] != 0:
                    k[r] = A * phi[r - s] / densStates[r]

    else:
        raise ReactionError('Unable to use inverse Laplace transform method for non-Arrhenius kinetics or for n < 0.')
    
    return k
    