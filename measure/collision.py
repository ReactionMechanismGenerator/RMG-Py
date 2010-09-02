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
Contains classes that represent the collision models available in MEASURE,
and methods for calculating various collision parameters. Each collision model
provides a method :meth:`generateCollisionMatrix()` that generates the collision
matrix :math:`\\matrix{M}_\\mathrm{coll} / \\omega = \\matrix{P} - \\matrix{I}`
corresponding to the collisional energy transfer probability function
:math:`P(E, E^\\prime)` for that model. The available collision models are:

* :class:`SingleExponentialDownModel`

"""

import math
import numpy
import logging

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
    Calculate the Lennard-Jones collision frequency for a given `species` with
    a dictionary of bath gases and their mole fractions `bathGas` at a given
    temperature `T` in K and pressure `P` in Pa. The Lennard-Jones model is
    generally a slight underestimate, but reasonable enough. If the bath gas
    is a mixture, arithmetic means are used to compute its effective
    Lennard-Jones :math:`\\sigma` parameter and molecular weight, while a
    geometric mean is used to calculate its effective Lennard-Jones
    :math:`\\epsilon` parameter.
    """
    
    bathGasSigma = 0.0; bathGasEpsilon = 1.0; bathGasMW = 0.0
    for key, value in bathGas.iteritems():
        try:
            bathGasSigma += key.lennardJones.sigma * value
            bathGasEpsilon *= key.lennardJones.epsilon ** value
        except AttributeError:
            raise CollisionError('No Lennard-Jones parameters specified for component "%s" in bath gas.' % bathGas)
        bathGasMW += key.molecularWeight * value

    gasConc = P / constants.kB / T
    mu = 1.0 / (1.0/species.molecularWeight + 1.0/bathGasMW) / 6.022e23
    try:
        sigma = 0.5 * (species.lennardJones.sigma + bathGasSigma)
        epsilon = math.sqrt(species.lennardJones.epsilon * bathGasEpsilon)
    except AttributeError:
        raise CollisionError('No Lennard-Jones parameters specified for species "%s".' % species)

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
    in K and are described by the collision model `collisionModel`, which
    currently must be a :class:`SingleExponentialDownModel` object. The
    algorithm here is implemented as described by Chang, Bozzelli, and Dean
    [Chang2000]_.

    .. [Chang2000] A. Y. Chang, J. W. Bozzelli, and A. M. Dean.
       *Z. Phys. Chem.* **214**, p. 1533-1568 (2000).
       `doi: 10.1524/zpch.2000.214.11.1533 <http://dx.doi.org/10.1524/zpch.2000.214.11.1533>`_

    """

    if not isinstance(collisionModel, SingleExponentialDownModel):
        raise CollisionError('Calculation of collision efficients requires the single exponential down collision model.')
    alpha = collisionModel.getAlpha(T)
    
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
        logging.warning('Collision efficiency %s calculated at %s K is greater than unity, so it will be set to unity.' % (beta, T))
    if beta < 0:
        raise CollisionError('Invalid collision efficiency %s calculated at %s K.' % (beta, T))
    
    return beta

################################################################################

class CollisionModel:
    """
    A base class for collision models. To create a custom collision model,
    derive from this class and implement the :meth:`generateCollisionMatrix()`
    method, which returns the collision matrix for the collision model you are
    implementing.

    .. note:: As with all collision models, you can only specify either the
        deactivating direction or the activating direction of the collisional
        transfer probabilities function :math:`P(E, E^\\prime)`, as the other
        is constrained by detailed balance.
       
    """
    pass

################################################################################

class SingleExponentialDownModel(CollisionModel):
    r"""
    A single exponential down collision model, based around the collisional 
    energy transfer probability function
    
    .. math:: P(E, E^\prime) = C(E^\prime) \exp \left( - \frac{E^\prime - E}{\alpha} \right) \hspace{40pt} E < E^\prime
    
    where the parameter :math:`\alpha = \left< \Delta E_\mathrm{d} \right>`
    represents the average energy transferred in a deactivating collision. This
    is the most commonly-used collision model, simply because it only has one
    parameter to determine. The parameter :math:`\alpha` is specified using the
    equation

    .. math:: \alpha = \alpha_0 \left( \frac{T}{T_0} \right)^n

    where :math:`\alpha_0` is the value of :math:`\alpha` at temperature
    :math:`T_0` in K. Set the exponent :math:`n` to zero to obtain a
    temperature-independent value for :math:`\alpha`.
    
    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `alpha`         ``double``      The average energy transferred in a deactivating collision in J/mol
    =============== =============== ============================================
    
    """

    def __init__(self, alpha0=0.0, T0=1.0, n=0.0):
        self.alpha0 = alpha0
        self.T0 = T0
        self.n = n

    def getAlpha(self, T):
        """
        Return the value of the :math:`\\alpha` parameter at temperature `T` in
        K. The :math:`\\alpha` parameter represents the average energy
        transferred in a deactivating collision
        :math:`\\left< \\Delta E_\\mathrm{d} \\right>`, and has units of J/mol.
        """
        return self.alpha0 * (T / self.T0) ** self.n

    def generateCollisionMatrix(self, Elist, T, densStates):
        """
        Generate and return the collision matrix
        :math:`\\matrix{M}_\\mathrm{coll} / \\omega = \\matrix{P} - \\matrix{I}`
        corresponding to this collision model for a given set of energies
        `Elist` in J/mol, temperature `T` in K, and isomer density of states
        `densStates`.
        """
        Ngrains = len(Elist)
        P = numpy.zeros((Ngrains,Ngrains), numpy.float64)
        
        start = -1
        for i in range(Ngrains):
            if densStates[i] > 0 and start == -1:
                start = i
                break

        # Determine value of parameters at this temperature
        alpha = self.getAlpha(T)

        # Determine unnormalized entries in collisional transfer probability matrix
        for r in range(start, Ngrains):
            P[0:r+1,r] = numpy.exp(-(Elist[r] - Elist[0:r+1]) / alpha)
            P[r+1:,r] = numpy.exp(-(Elist[r+1:] - Elist[r]) / alpha) * densStates[r+1:] / densStates[r] * numpy.exp(-(Elist[r+1:] - Elist[r]) / (constants.R * T))
        
        # Normalize using detailed balance
        # This method is much more robust, and corresponds to:
        #    [ 1 1 1 1 ...]
        #    [ 1 2 2 2 ...]
        #    [ 1 2 3 3 ...]
        #    [ 1 2 3 4 ...]
        for r in range(start, Ngrains):
            C = (1 - numpy.sum(P[start:r,r])) / numpy.sum(P[r:Ngrains,r])
            # Check for normalization consistency (i.e. all numbers are positive)
            if C < 0: raise CollisionError('Encountered negative normalization coefficient while normalizing collisional transfer probabilities matrix.')
            P[r,r+1:Ngrains] *= C
            P[r:Ngrains,r] *= C
            P[r,r] -= 1
        # This method is described by Pilling and Holbrook, and corresponds to:
        #    [ ... 4 3 2 1 ]
        #    [ ... 3 3 2 1 ]
        #    [ ... 2 2 2 1 ]
        #    [ ... 1 1 1 1 ]
        #for r in range(Ngrains, start, -1):
            #C = (1 - numpy.sum(M[r:Ngrains,r])) / numpy.sum(M[0:r,r])
            ## Check for normalization consistency (i.e. all numbers are positive)
            #if C < 0: raise CollisionError('Encountered negative normalization coefficient while normalizing collisional transfer probabilities matrix.')
            #P[r,0:r-1] *= C
            #P[0:r,r] *= C
            #P[r,r] -= 1

        return P

################################################################################
