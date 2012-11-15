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
Contains functions for working with reaction rate coefficients. The
:meth:`calculateMicrocanonicalRateCoefficients()` function is used to calculate
the forward and reverse microcanonical rate coefficients :math:`k(E)` for each
path reaction. The :meth:`fitInterpolationModel()` function is used to fit an
interpolation model to a phenomenological rate coefficient :meth:`k(T,P)` for
each net reaction.
"""

cdef extern from "math.h":
    cdef double exp(double x)
    cdef double floor(double x)
    cdef double log(double x)
    cdef double sqrt(double x)

import numpy
cimport numpy
import logging

cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity
import rmgpy.reaction
from rmgpy.kinetics import Arrhenius, PDepArrhenius, Chebyshev

def convolve(numpy.ndarray[numpy.float64_t,ndim=1] rho1,
    numpy.ndarray[numpy.float64_t,ndim=1] rho2,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist):
    """
    Convolutes two density of states arrays `rho1` and `rho2` with corresponding
    energies `Elist` together using the equation

    .. math:: \\rho(E) = \\int_0^E \\rho_1(x) \\rho_2(E-x) \\, dx

    The units of the parameters do not matter so long as they are consistent.
    """

    cdef numpy.ndarray[numpy.float64_t,ndim=1] rho
    cdef double dE
    cdef int nE, i, j

    rho = numpy.zeros_like(Elist)
    dE = Elist[1] - Elist[0]
    nE = len(Elist)
    for i in range(nE):
        for j in range(i+1):
            rho[i] += rho2[i-j] * rho1[i] * dE

    return rho

################################################################################

class ReactionError(Exception): 
    """
    An exception raised when working with reactions causes exceptional behavior
    for any reason. Pass a string describing the cause of the exceptional 
    behavior.
    """
    pass

################################################################################

def calculateMicrocanonicalRateCoefficient(reaction, 
    numpy.ndarray[numpy.float64_t,ndim=1] Elist, 
    numpy.ndarray[numpy.float64_t,ndim=1] reacDensStates, 
    numpy.ndarray[numpy.float64_t,ndim=1] prodDensStates=None, 
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

    cdef numpy.ndarray[numpy.float64_t,ndim=1] kf, kr
    cdef double Keq, reacQ, prodQ, R = constants.R
    cdef int r
    cdef bint reactantStatesKnown, productStatesKnown

    kf = numpy.zeros_like(Elist)
    kr = numpy.zeros_like(Elist)

    reactantStatesKnown = reacDensStates is not None and reacDensStates.any()
    productStatesKnown = prodDensStates is not None and prodDensStates.any()

    if reaction.transitionState.states is not None:
        # We've been provided with molecular degree of freedom data for the
        # transition state, so let's use the more accurate RRKM theory
        if reactantStatesKnown and (reaction.isIsomerization() or reaction.isDissociation()):
            kf = applyRRKMTheory(reaction.transitionState, Elist, reacDensStates)
        elif productStatesKnown and reaction.isAssociation():
            kr = applyRRKMTheory(reaction.transitionState, Elist, prodDensStates)
        else:
            raise ReactionError('Unable to compute k(E) values via RRKM theory for path reaction "{0}".'.format(reaction))
    
    elif reaction.kinetics is not None:
        # We've been provided with high-pressure-limit rate coefficient data,
        # so let's use the less accurate inverse Laplace transform method
        if reactantStatesKnown:
            kf = applyInverseLaplaceTransformMethod(reaction.kinetics, reaction.transitionState.E0.value_si, Elist, reacDensStates, T)
        elif productStatesKnown:
            kinetics = reaction.generateReverseRateCoefficient()
            kr = applyInverseLaplaceTransformMethod(kinetics, reaction.transitionState.E0.value_si, Elist, prodDensStates, T)
        else:
            raise ReactionError('Unable to compute k(E) values via ILT method for path reaction "{0}".'.format(reaction))
    
    else:
        raise ReactionError('Unable to compute k(E) values for path reaction "{0}".'.format(reaction))
    
    # If the reaction is endothermic and barrierless, it is possible that the
    # forward k(E) will have a nonzero value at an energy where the product
    # density of states is zero (but the reactant density of states is not),
    # which violates detailed balance
    # To fix, we set the forward k(E) to zero wherever this is true
    # (This is correct within the accuracy of discretizing the energy grains)
    if kf.any() and productStatesKnown:
        for r in range(len(Elist)):
            if reacDensStates[r] != 0 and prodDensStates[r] != 0:
                break
            kf[r] = 0
    if kr.any() and reactantStatesKnown:
        for r in range(len(Elist)):
            if reacDensStates[r] != 0 and prodDensStates[r] != 0:
                break
            kr[r] = 0
    
    # Get the reverse microcanonical rate coefficient if possible
    if (all([reactant.thermo is not None for reactant in reaction.reactants]) and 
        all([product.thermo is not None for product in reaction.products])):

        # Determine the parameters we will need to compute the reverse k(E)
        Keq = reaction.getEquilibriumConstant(T, 'Kc')
        if reactantStatesKnown:
            reacEqDist = reacDensStates * numpy.exp(-Elist / R / T)
            reacQ = numpy.sum(reacEqDist)
        if productStatesKnown:
            prodEqDist = prodDensStates * numpy.exp(-Elist / R / T)
            prodQ = numpy.sum(prodEqDist)
            
        if kf.any():
            # We computed the forward rate coefficient above
            # Thus we need to compute the reverse rate coefficient here
            if reaction.isIsomerization() and productStatesKnown:
                for r in range(len(Elist)):
                    if prodDensStates[r] > 0: break
                kr[r:] = kf[r:] * (reacDensStates[r:] / prodDensStates[r:]) * (prodQ / reacQ) / Keq
            elif reaction.isDissociation():
                kr = kf * (reacEqDist / reacQ) / Keq
            elif reaction.isAssociation():
                if productStatesKnown:
                    for r in range(len(Elist)):
                        if prodDensStates[r] > 0: break
                    kr[r:] = kf[r:] * (reacDensStates[r:] / prodDensStates[r:]) * (prodQ / reacQ) / Keq
                kf = kf * reacEqDist / reacQ

        elif kr.any():
            # We computed the reverse rate coefficient above
            # Thus we need to compute the forward rate coefficient here
            if reaction.isIsomerization() and reactantStatesKnown:
                for r in range(len(Elist)):
                    if reacDensStates[r] > 0: break
                kf[r:] = kr[r:] * (prodDensStates[r:] / reacDensStates[r:]) * (reacQ / prodQ) * Keq
            elif reaction.isAssociation():
                kf = kr * (prodEqDist / prodQ) * Keq
            elif reaction.isDissociation():
                if reactantStatesKnown:
                    for r in range(len(Elist)):
                        if reacDensStates[r] > 0: break
                    kf[r:] = kr[r:] * (prodDensStates[r:] / reacDensStates[r:]) * (reacQ / prodQ) * Keq
                kr = kr * prodEqDist / prodQ

    return kf, kr

################################################################################

def applyRRKMTheory(transitionState, 
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=1] densStates):
    """
    Calculate the microcanonical rate coefficient for a reaction using RRKM
    theory, where `transitionState` is the transition state of the reaction,
    `Elist` is the array of energies in J/mol at which to evaluate the
    microcanonial rate, and `densStates` is the density of states of the
    reactant.
    """
    
    k = numpy.zeros_like((Elist))
    Ngrains = len(Elist)
    dE = Elist[1] - Elist[0]
    
    # Calculate sum of states of transition state
    sumStates0 = transitionState.states.getSumOfStates(Elist - Elist[0])
    # Shift to common zero of energy
    r0 = int(round((transitionState.E0.value_si - Elist[0]) / dE))
    sumStates = numpy.zeros_like(Elist)
    sumStates[r0:] = sumStates0[:-r0+len(sumStates0)]
    
    # Generate k(E) using RRKM formula
    for r in range(len(Elist)):
        if densStates[r] > 0:
            k[r] = sumStates[r] / constants.h / densStates[r] / constants.Na
    
    return k

################################################################################

def applyInverseLaplaceTransformMethod(kinetics, double E0,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=1] densStates,
    double T=0.0):
    """
    Calculate the microcanonical rate coefficient for a reaction using the
    inverse Laplace transform method, where `kinetics` is the high pressure 
    limit rate coefficient, `E0` is the ground-state energy of the transition
    state, `Elist` is the array of energies in J/mol at which to evaluate the
    microcanonial rate, and `densStates` is the density of states of the
    reactant. The temperature `T` in K is not required, and is only used when
    the temperature exponent of the Arrhenius expression is negative (for which
    the inverse transform is undefined).
    """

    cdef double A, n, T0, Ea, dE, s0, rem, R = constants.R
    cdef numpy.ndarray[numpy.float64_t,ndim=1] k, phi
    cdef int i, r, s, Ngrains = len(Elist)

    k = numpy.zeros_like((Elist))
    
    if isinstance(kinetics, Arrhenius) and (T != 0.0 or (kinetics.Ea >= 0 and kinetics.n >= 0)):
        A = kinetics.A.value_si / (kinetics.T0.value_si**kinetics.n.value_si)
        n = kinetics.n.value_si
        Ea = kinetics.Ea.value_si
        dE = Elist[1] - Elist[0]

        # The inverse Laplace transform is not defined for Ea < 0 or n < 0
        # In these cases we move the offending portion into the preexponential
        # at the temperature of interest
        # This is an approximation, but it's not worth a more robust procedure
        if Ea < 0:
            A *= exp(-Ea / R / T)
            Ea = 0.0
        if n < 0:
            A *= T**n
            n = 0.0

        if n < 0.001:
            # Determine the microcanonical rate directly
            s0, rem = divmod(Ea, dE)
            s = int(s0)
            if rem == 0:
                for r in range(s, Ngrains):
                    if Elist[r] > E0 and densStates[r] != 0:
                        k[r] = A * densStates[r-s] / densStates[r]
            else:
                for r in range(s+1, Ngrains):
                    if Elist[r] > E0 and densStates[r] != 0 and densStates[r-s] != 0:
                        num = densStates[r-s] * (densStates[r-s-1] / densStates[r-s]) ** (-rem / (Elist[r-s-1] - Elist[r-s]))
                        k[r] = A * num / densStates[r]
                    
        elif n >= 0.001:
            import scipy.special
            # Evaluate the inverse Laplace transform of the T**n piece, which only
            # exists for n >= 0
            phi = numpy.zeros(Ngrains, numpy.float64)
            for i in range(1, Ngrains):
                phi[i] = (Elist[i] - Elist[0])**(n-1) / (R**n * scipy.special.gamma(n))
            # Evaluate the convolution
            phi = convolve(phi, densStates, Elist)
            # Apply to determine the microcanonical rate
            s0, rem = divmod(Ea, dE)
            s = int(s0)
            if rem == 0:
                for r in range(s, Ngrains):
                    if Elist[r] > E0 and densStates[r] != 0:
                        k[r] = A * phi[r-s] / densStates[r]
            else:
                for r in range(s+1, Ngrains):
                    if Elist[r] > E0 and densStates[r] != 0 and phi[r-s] != 0:
                        num = phi[r-s] * (phi[r-s-1] / phi[r-s]) ** (-rem / (Elist[r-s-1] - Elist[r-s]))
                        k[r] = A * num / densStates[r]

    else:
        raise ReactionError('Unable to use inverse Laplace transform method for non-Arrhenius kinetics or for n < 0.')
    
    return k

################################################################################

def fitInterpolationModel(reaction, Tlist, Plist, K, model, Tmin, Tmax, Pmin, Pmax, errorCheck=False):
    """
    For a set of phenomenological rate coefficients `K` computed at a grid of
    temperatures `Tlist` in K and pressures `Plist` in Pa, fit a :math:`k(T,P)`
    interpolation `model`, a tuple where the first item is a string describing
    the type of model - either ``'chebyshev'`` or ``'pdeparrhenius'`` - and the
    remaining elements contain parameters for that model. For Chebyshev
    polynomials, the parameters are the number of terms to use in each of the
    temperature and pressure dimensions. For pressure-dependent Arrhenius models
    there are no additional parameters. `Tmin`, `Tmax`, `Pmin`, and `Pmax`
    specify the temperature and pressure ranges in K and Pa, respectively,
    over which the interpolation model is valid. If `errorCheck` is ``True``,
    a check will be performed to ensure that the interpolation model does not
    deviate too much from the data; as this is not necessarily a fast process,
    it is optional.
    """

    # Determine units for k(T,P)
    if len(reaction.reactants) == 1:
        kunits = 's^-1'
    elif len(reaction.reactants) == 2:
        kunits = 'm^3/(mol*s)'
    elif len(reaction.reactants) == 3:
        kunits = 'm^6/(mol^2*s)'
    else:
        kunits = ''

    # Set/update the net reaction kinetics using interpolation model
    if model[0].lower() == 'chebyshev':
        modelType, degreeT, degreeP = model
        chebyshev = Chebyshev()
        chebyshev.fitToData(Tlist, Plist, K, kunits, degreeT, degreeP, Tmin, Tmax, Pmin, Pmax)
        kinetics = chebyshev
    elif model[0].lower() == 'pdeparrhenius':
        pDepArrhenius = PDepArrhenius()
        pDepArrhenius.fitToData(Tlist, Plist, K, kunits=kunits, T0=298.0)
        kinetics = pDepArrhenius
    else:
        return None

    # Set temperature and pressure ranges explicitly (as they may be different
    # from min(Tlist), max(Tlist), min(Plist), max(Plist))
    kinetics.Tmin = Quantity(Tmin,"K")
    kinetics.Tmax = Quantity(Tmax,"K")
    kinetics.Pmin = Quantity(Pmin/1e5,"bar")
    kinetics.Pmax = Quantity(Pmax/1e5,"bar")

    # Compute log RMS error for fit
    if errorCheck:
        logRMS = 0.0
        # Check that fit is within an order of magnitude at all points
        for t, T in enumerate(Tlist):
            for p, P in enumerate(Plist):
                logkmodel = log(kinetics.getRateCoefficient(T, P))
                logkdata = log(K[t,p])
                logRMS += (logkmodel - logkdata) * (logkmodel - logkdata)
        logRMS = sqrt(logRMS / len(Tlist) / len(Plist))
        if logRMS > 0.5:
            logging.warning('RMS error for k(T,P) fit = {0:g} for reaction {1}.'.format(logRMS, reaction))
    
    return kinetics
