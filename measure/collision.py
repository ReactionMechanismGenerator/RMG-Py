#!/usr/bin/python
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
Contains classes that represent the collision models available in MEASURE.
Each collision model provides a collisional energy transfer probability function
that returns the value of :math:`P(E, E^\prime)` for that model.
"""

import math
import numpy

import chempy.constants as constants

################################################################################

class CollisionError(Exception): 
    """
    An exception raised when working with collision models causes exceptional 
    behavior for any reason. Pass a string describing the cause of the  
    exceptional behavior.
    """
    pass

################################################################################

def calculateCollisionFrequency(species, T, P, bathGas):
    """
    Calculate the collision frequency for a given `species` with a bath gas
    `bathGas` at a given temperature `T` in K and pressure `P` in Pa. The
    Lennard-Jones collision model is used, which generally is a slight 
    underestimate, but reasonable enough.
    """
    
    gasConc = P / constants.kB / T
    mu = 1 / (1/species.molecularWeight + 1/bathGas.molecularWeight) / 6.022e23
    sigma = 0.5 * (species.lennardJones.sigma + bathGas.lennardJones.sigma)
    epsilon = math.sqrt(species.lennardJones.epsilon * bathGas.lennardJones.epsilon)

    # Evaluate configuration integral
    Tred = constants.kB * T / epsilon
    omega22 = 1.16145 * Tred**(-0.14874) + 0.52487 * math.exp(-0.77320 * Tred) + 2.16178 * math.exp(-2.43787 * Tred)

    # Evaluate collision frequency
    return omega22 * math.sqrt(8 * constants.kB * T / math.pi / mu) * math.pi * sigma**2 * gasConc

################################################################################

def calculateCollisionEfficiency(species, T, Elist, densStates, collisionModel, E0, Ereac):
    """
    Calculate an efficiency factor for collisions, particularly useful for the
    modified strong collision method. The collisions involve the given 
    `species` with density of states `densStates` in mol/J corresponding to
    energies `Elist` in J/mol, ground-state energy `E0` in J/mol, and first 
    reactive energy `Ereac` in J/mol. The collisions occur at temperature `T` 
    in K and are described by the collision model `collisionModel`. The
    algorithm here is implemented as described by Chang, Bozzelli, and Dean.
    """

    if not isinstance(collisionModel, SingleExponentialDownModel):
        raise CollisionError('Modified strong collision method requires the single exponential down collision model.')
    alpha = collisionModel.alpha
    
    # Ensure that the barrier height is sufficiently above the ground state
    # Otherwise invalid efficiencies are observed
    if Ereac - E0 < 100000:
        Ereac = E0 + 100000

    Ngrains = len(Elist)
    dE = Elist[1] - Elist[0]
    FeNum = 0; FeDen = 0
    Delta1 = 0; Delta2 = 0; DeltaN = 0; Delta = 1

    for r in range(Ngrains):
        value = densStates[r] * math.exp(-Elist[r] / constants.R / T)
        if Elist[r] > Ereac:
            FeNum += value * dE
            if FeDen == 0:
                FeDen = value * constants.R * T
    if FeDen == 0: return 1.0
    Fe = FeNum / FeDen

    # Chang, Bozzelli, and Dean recommend "freezing out" Fe at values greater
    # than 1e6 to avoid issues of roundoff error
    # They claim that the collision efficiency isn't too temperature-dependent
    # in this regime, so it's an okay approximation to use
    if Fe > 1e6: Fe = 1e6
    
    for r in range(Ngrains):
        value = densStates[r] * math.exp(-Elist[r] / constants.R / T)
        # Delta
        if Elist[r] < Ereac:
            Delta1 += value * dE
            Delta2 += value * dE * math.exp(-(Ereac - Elist[r]) / (Fe * constants.R * T))
        DeltaN += value * dE

    Delta1 /= DeltaN
    Delta2 /= DeltaN

    Delta = Delta1 - (Fe * constants.R * T) / (alpha + Fe * constants.R * T) * Delta2

    beta = (alpha / (alpha + Fe * constants.R * T))**2 / Delta

    if beta > 1:
        logging.warning('Collision efficiency %s calculated at %s K is greater than unity, so it will be set to unity..' % (beta, T))
    if beta < 0:
        raise CollisionError('Invalid collision efficiency %s calculated at %s K.' % (beta, T))
    
    return beta

################################################################################

class CollisionModel:
    """
    A base class for collision models.
    """
    pass

################################################################################

class SingleExponentialDownModel(CollisionModel):
    """
    A single exponential down collision model, based around the collisional 
    energy transfer probability function
    
    .. math:: P(E, E^\prime) = C(E^\prime) \exp \left( - \frac{E^\prime - E}{\alpha} \right) \hspace{40pt} E < E^\prime
    
    where the parameter :math:`\alpha = \left< \Delta E_\mathrm{d} \right>`
    represents the average energy transferred in a deactivating collision.
    
    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `alpha`         ``double``      The average energy transferred in a deactivating collision in J/mol
    =============== =============== ============================================
    
    """

    def __init__(self, alpha=0.0):
        self.alpha = alpha

    def generateCollisionMatrix(self, Elist, T, densStates):
        """
        Generate and return the collisional transfer probability matrix 
        :math:`P(E, E^\prime)` for this model for a given
        set of energies `Elist` in J/mol, temperature `T` in K, and isomer 
        density of states `densStates`.
        """
        Ngrains = len(Elist)
        P = numpy.zeros((Ngrains,Ngrains), numpy.float64)
        
        start = -1
        for i in range(Ngrains):
            if densStates[i] > 0 and start == -1:
                start = i
                break
        
        # Determine unnormalized entries in collisional transfer probability matrix
        for r in range(start, Ngrains):
            for s in range(start, r+1):
                # Er >= Es
                P[s,r] = math.exp(-(Elist[r] - Elist[s]) / self.alpha)
            for s in range(r+1, Ngrains):
                # Er < Es
                P[s,r] = math.exp(-(Elist[s] - Elist[r]) / self.alpha) * densStates[s] / densStates[r] * math.exp(-(Elist[s] - Elist[r]) / (constants.R * T))
        
        # Normalize using detailed balance
        # This method is much more robust, and corresponds to:
        #    [ 1 1 1 1 ...]
        #    [ 1 2 2 2 ...]
        #    [ 1 2 3 3 ...]
        #    [ 1 2 3 4 ...]
        for r in range(start, Ngrains):
            C = (1 - numpy.sum(P[start:r,r])) / sum(P[r:Ngrains,r])
            # Check for normalization consistency (i.e. all numbers are positive)
            if C < 0: raise ChemPyError('Encountered negative normalization coefficient while normalizing collisional transfer probabilities matrix.')
            P[r,r+1:Ngrains] *= C
            P[r:Ngrains,r] *= C
            P[r,r] -= 1
        # This method is described by Pilling and Holbrook, and corresponds to:
        #    [ ... 4 3 2 1 ]
        #    [ ... 3 3 2 1 ]
        #    [ ... 2 2 2 1 ]
        #    [ ... 1 1 1 1 ]
        #for r in range(Ngrains, start, -1):
            #C = (1 - numpy.sum(M[r:Ngrains,r])) / sum(M[0:r,r])
            ## Check for normalization consistency (i.e. all numbers are positive)
            #if C < 0: raise ChemPyError('Encountered negative normalization coefficient while normalizing collisional transfer probabilities matrix.')
            #P[r,0:r-1] *= C
            #P[0:r,r] *= C
            #P[r,r] -= 1

        return P

################################################################################
