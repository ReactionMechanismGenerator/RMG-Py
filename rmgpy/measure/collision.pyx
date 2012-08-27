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

cdef extern from "math.h":
    cdef double exp(double x)
    cdef double sqrt(double x)

import math
import numpy
cimport numpy
import logging
import cython

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity, ScalarQuantity, ArrayQuantity
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

class CollisionError(Exception): 
    """
    An exception raised when working with collision models causes exceptional 
    behavior for any reason. Pass a string describing the cause of the  
    exceptional behavior.
    """
    pass

################################################################################

cpdef double calculateCollisionFrequency(species, double T, double P, bathGas):
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
   
    cdef double gasConc, sigma, epsilon, Tred, omega22, value
    cdef double bathGasSigma, bathGasEpsilon, bathGasMW
    cdef double kB = constants.kB, pi = math.pi

    bathGasSigma = 0.0; bathGasEpsilon = 1.0; bathGasMW = 0.0
    for key, value in bathGas.iteritems():
        try:
            bathGasSigma += key.lennardJones.sigma.value_si * value
            bathGasEpsilon *= key.lennardJones.epsilon.value_si ** value
        except AttributeError:
            raise CollisionError('No Lennard-Jones parameters specified for component "{0}" in bath gas.'.format(bathGas))
        bathGasMW += key.molecularWeight.value_si * value
        
    gasConc = P / constants.kB / T
    mu = 1.0 / (1.0/species.molecularWeight.value_si + 1.0/bathGasMW) / 6.022e23
    try:
        sigma = 0.5 * (species.lennardJones.sigma.value_si + bathGasSigma)
        epsilon = math.sqrt(species.lennardJones.epsilon.value_si * bathGasEpsilon)
    except AttributeError:
        raise CollisionError('No Lennard-Jones parameters specified for species "{0}".'.format(species))

    # Evaluate configuration integral
    Tred = kB * T / epsilon
    omega22 = 1.16145 * Tred**(-0.14874) + 0.52487 * exp(-0.77320 * Tred) + 2.16178 * exp(-2.43787 * Tred)

    # Evaluate collision frequency
    return omega22 * sqrt(8 * kB * T / pi / mu) * pi * sigma*sigma * gasConc

################################################################################

def calculateCollisionEfficiency(double T,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=1] densStates,
    double dEdown, double E0, double Ereac):
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

    cdef double dE, FeNum, FeDen, Delta1, Delta2, DeltaN, Delta, value, beta
    cdef double R = constants.R
    cdef int Ngrains, r

    # Ensure that the barrier height is sufficiently above the ground state
    # Otherwise invalid efficiencies are observed
    if Ereac - E0 < 100000:
        Ereac = E0 + 100000

    Ngrains = len(Elist)
    dE = Elist[1] - Elist[0]
    FeNum = 0; FeDen = 0
    Delta1 = 0; Delta2 = 0; DeltaN = 0; Delta = 1

    for r in range(Ngrains):
        value = densStates[r] * exp(-Elist[r] / R / T)
        if Elist[r] > Ereac:
            FeNum += value * dE
            if FeDen == 0:
                FeDen = value * R * T
    if FeDen == 0: return 1.0
    Fe = FeNum / FeDen

    # Chang, Bozzelli, and Dean recommend "freezing out" Fe at values greater
    # than 1e6 to avoid issues of roundoff error
    # They claim that the collision efficiency isn't too temperature-dependent
    # in this regime, so it's an okay approximation to use
    if Fe > 1e6: Fe = 1e6
    
    for r in range(Ngrains):
        value = densStates[r] * exp(-Elist[r] / R / T)
        # Delta
        if Elist[r] < Ereac:
            Delta1 += value * dE
            Delta2 += value * dE * exp(-(Ereac - Elist[r]) / (Fe * R * T))
        DeltaN += value * dE

    Delta1 /= DeltaN
    Delta2 /= DeltaN

    Delta = Delta1 - (Fe * R * T) / (dEdown + Fe * R * T) * Delta2

    beta = (dEdown / (dEdown + Fe * R * T))**2 / Delta

    if beta > 1:
        logging.warning('Collision efficiency {0:.3f} calculated at {1:g} K is greater than unity, so it will be set to unity.'.format(beta, T))
        beta = 1
    if beta < 0:
        raise CollisionError('Invalid collision efficiency {0:.3f} calculated at {1:g} K.'.format(beta, T))
    
    return beta

################################################################################

cdef class CollisionModel:
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

cdef class SingleExponentialDown(CollisionModel):
    r"""
    Refactoring of collision and reaction modules to full Cython syntax.
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

    cdef public ScalarQuantity alpha0, T0, n

    def __init__(self, alpha0=None, T0=None, n=None):
        if alpha0 is not None:
            self.alpha0 = Quantity(alpha0)
            if self.alpha0.units == 'cm^-1':
                self.alpha0.value_si *= 11.96
                self.alpha0.units = 'J/mol'
        else:
            self.alpha0 = None
        if T0 is not None:
            self.T0 = Quantity(T0)
        else:
            self.T0 = None
        if n is not None:
            self.n = Quantity(n)
        else:
            self.n = None

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        SingleExponentialDownModel object.
        """
        return 'SingleExponentialDown(alpha0={0!r}, T0={1!r}, n={2!r})'.format(self.alpha0, self.T0, self.n)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (SingleExponentialDown, (self.alpha0, self.T0, self.n))

    cpdef double getAlpha(self, double T):
        """
        Return the value of the :math:`\\alpha` parameter at temperature `T` in
        K. The :math:`\\alpha` parameter represents the average energy
        transferred in a deactivating collision
        :math:`\\left< \\Delta E_\\mathrm{d} \\right>`, and has units of J/mol.
        """
        if self.alpha0 is None:
            return 0.0
        elif self.T0 is None or self.n is None:
            return self.alpha0.value_si
        else:
            return self.alpha0.value_si * (T / self.T0.value_si) ** self.n.value_si

    def generateCollisionMatrix(self,
        numpy.ndarray[numpy.float64_t,ndim=1] Elist,
        double T,
        numpy.ndarray[numpy.float64_t,ndim=1] densStates,
        double dEdown):
        """
        Generate and return the collision matrix
        :math:`\\matrix{M}_\\mathrm{coll} / \\omega = \\matrix{P} - \\matrix{I}`
        corresponding to this collision model for a given set of energies
        `Elist` in J/mol, temperature `T` in K, and isomer density of states
        `densStates`. You must also supply the average energy transferred in
        a deactivating collision `dEdown` in J/mol.
        """

        cdef double alpha = 1.0/dEdown, beta = 1.0 / (constants.R * T)
        cdef double C, left, right
        cdef int Ngrains, start, i, r
        cdef numpy.ndarray[numpy.float64_t,ndim=2] P

        Ngrains = len(Elist)
        P = numpy.zeros((Ngrains,Ngrains), numpy.float64)

        start = -1
        for i in range(Ngrains):
            if densStates[i] > 0 and start == -1:
                start = i
                break

        # Determine unnormalized entries in collisional transfer probability matrix
        for r in range(start, Ngrains):
            for s in range(start,r+1):
                P[s,r] = exp(-(Elist[r] - Elist[s]) * alpha)
            for s in range(r+1,Ngrains):
                P[s,r] = exp(-(Elist[s] - Elist[r]) * alpha) * densStates[s] / densStates[r] * exp(-(Elist[s] - Elist[r]) * beta)
        
        # Normalize using detailed balance
        # This method is much more robust, and corresponds to:
        #    [ 1 1 1 1 ...]
        #    [ 1 2 2 2 ...]
        #    [ 1 2 3 3 ...]
        #    [ 1 2 3 4 ...]
        for r in range(start, Ngrains):
            left = 0.0; right = 0.0
            for s in range(start, r): left += P[s,r]
            for s in range(r, Ngrains): right += P[s,r]
            C = (1 - left) / right
            # Check for normalization consistency (i.e. all numbers are positive)
            if C < 0: raise CollisionError('Encountered negative normalization coefficient while normalizing collisional transfer probabilities matrix.')
            for s in range(r+1,Ngrains):
                P[r,s] *= C
                P[s,r] *= C
            P[r,r] = P[r,r] * C - 1
        # This method is described by Pilling and Holbrook, and corresponds to:
        #    [ ... 4 3 2 1 ]
        #    [ ... 3 3 2 1 ]
        #    [ ... 2 2 2 1 ]
        #    [ ... 1 1 1 1 ]
        #for r in range(Ngrains, start, -1):
            #left = 0.0; right = 0.0
            #for s in range(start, r): left += P[s,r]
            #for s in range(r, Ngrains): right += P[s,r]
            #C = (1 - right) / left
            ## Check for normalization consistency (i.e. all numbers are positive)
            #if C < 0: raise CollisionError('Encountered negative normalization coefficient while normalizing collisional transfer probabilities matrix.')
            #for s in range(r-1):
                #P[r,s] *= C
                #P[s,r] *= C
            #P[r,r] = P[r,r] * C - 1

        return P

################################################################################
