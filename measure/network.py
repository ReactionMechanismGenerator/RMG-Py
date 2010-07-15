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
Contains classes that define an internal representation of a unimolecular
reaction network.
"""

import math
import numpy
import cython
import logging

import chempy.constants as constants
import chempy.states as states

from reaction import *
from collision import *

################################################################################

class NetworkError(Exception):
    """
    An exception raised while manipulating unimolecular reaction networks for
    any reason. Pass a string describing the cause of the exceptional behavior.
    """
    pass

################################################################################

class Network:
    """
    A representation of a unimolecular reaction network. The attributes are:

    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `isomers`           ``list``                A list of the unimolecular isomers in the network
    `reactants`         ``list``                A list of the bimolecular reactant channels in the network
    `products`          ``list``                A list of the bimolecular product channels in the network
    `pathReactions`     ``list``                A list of reaction objects that connect adjacent isomers (the high-pressure-limit)
    `bathGas`           :class:`Species`        The bath gas
    `collisionModel`    :class:`CollisionModel` The collision model to use
    `netReactions`      ``list``                A list of reaction objects that connect any pair of isomers
    =================== ======================= ================================

    """

    def __init__(self, isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None):
        self.isomers = isomers or []
        self.reactants = reactants or []
        self.products = products or []
        self.pathReactions = pathReactions or []
        self.bathGas = bathGas
        self.netReactions = []
    
    def getEnergyGrains(self, Emin, Emax, dE=0.0, Ngrains=0):
        """
        Return an array of energy grains that have a minimum of `Emin`, a
        maximum of `Emax`, and either a spacing of `dE` or have number of
        grains `nGrains`. The first three parameters are in J/mol, as is the
        returned array of energy grains.
        """
        useGrainSize = False

        if Ngrains <= 0 and dE <= 0.0:
            # Neither grain size nor number of grains specified, so raise exception
            raise NetworkError('You must specify a positive value for either dE or Ngrains.')
        elif Ngrains <= 0 and dE > 0.0:
            # Only grain size was specified, so we must use it
            useGrainSize = True
        elif Ngrains > 0 and dE <= 0.0:
            # Only number of grains was specified, so we must use it
            useGrainSize = False
        else:
            # Both were specified, so we choose the tighter constraint
            # (i.e. the one that will give more grains, and so better accuracy)
            dE0 = (Emax - Emin) / (Ngrains - 1)
            useGrainSize = (dE0 > dE)
        
        # Generate the array of energies
        if useGrainSize:
            return numpy.arange(Emin, Emax + dE, dE, numpy.float64)
        else:
            return numpy.linspace(Emin, Emax, Ngrains, numpy.float64)
        
    def autoGenerateEnergyGrains(self, Tmax, grainSize=0.0, Ngrains=0):
        """
        Select a suitable list of energies to use for subsequent calculations.
        The procedure is:

        1. Calculate the equilibrium distribution of the highest-energy isomer
           at the largest temperature of interest (to get the broadest
           distribution)

        2. Calculate the energy at which the tail of the distribution is some
           fraction of the maximum

        3. Add the difference between the ground-state energy of the isomer and
           the highest ground-state energy in the system (either isomer or
           transition state)

        You must specify either the desired grain spacing `grainSize` in J/mol 
        or the desired number of grains `Ngrains`, as well as a temperature 
        `Tmax` in K to use for the equilibrium calculation, which should be the
        highest temperature of interest. You can specify both `grainSize` and 
        `Ngrains`, in which case the one that gives the more accurate result 
        will be used (i.e. they represent a maximum grain size and a minimum
        number of grains). An array containing the energy grains in J/mol is 
        returned.
        """
        
        if grainSize == 0.0 and Ngrains == 0:
            raise NetworkError('Must provide either grainSize or Ngrains parameter to Network.determineEnergyGrains().')

        # For the purposes of finding the maximum energy we will use 251 grains
        nE = 251; dE = 0.0

        # The minimum energy is the lowest isomer energy on the PES
        Emin = 1.0e25
        for species in self.isomers:
            if species.E0 < Emin: Emin = species.E0
        Emin = math.floor(Emin) # Round to nearest whole number

        # Determine the isomer with the maximum ground-state energy
        isomer = None
        for species in self.isomers:
            if isomer is None: isomer = species
            elif species.E0 > isomer.E0: isomer = species
        Emax0 = isomer.E0

        # (Try to) purposely overestimate Emax using arbitrary multiplier
        # to (hopefully) avoid multiple density of states calculations
        mult = 50
        done = False
        maxIter = 5
        iterCount = 0
        while not done and iterCount < maxIter:

            iterCount += 1
            
            Emax = math.ceil(Emax0 + mult * constants.R * Tmax)

            Elist = self.getEnergyGrains(0.0, Emax-Emin, dE, nE)
            densStates = isomer.states.getDensityOfStates(Elist)
            eqDist = densStates * numpy.exp(-Elist / constants.R / Tmax)
            eqDist /= numpy.sum(eqDist)
            
            # Find maximum of distribution
            maxIndex = eqDist.argmax()
            
            # If tail of distribution is much lower than the maximum, then we've found bounds for Emax
            tol = 1e-4
            if eqDist[-1] / eqDist[maxIndex] < tol:
                r = nE - 1
                while r > maxIndex and not done:
                    if eqDist[r] / eqDist[maxIndex] > tol: done = True
                    else: r -= 1
                Emax = Elist[r] + Emin
                # A final check to ensure we've captured almost all of the equilibrium distribution
                if abs(1.0 - numpy.sum(eqDist[0:r]) / numpy.sum(eqDist)) > tol:
                    done = False
                    mult += 50
            else:
                mult += 50

        # Add difference between isomer ground-state energy and highest
        # transition state or reactant channel energy
        Emax0_iso = Emin
        for species in self.reactants:
            E = sum([spec.E0 for spec in species])
            if Emax0_iso < E: Emax0_iso = E
        Emax0_rxn = Emin
        for rxn in self.pathReactions:
            if rxn.transitionState is not None:
                E = rxn.transitionState.E0
                if Emax0_rxn < E: Emax0_rxn = E
        Emax += max([isomer.E0, Emax0_iso, Emax0_rxn]) - isomer.E0

        # Round Emax up to nearest integer
        Emax = math.ceil(Emax)

        # Return the chosen energy grains
        return self.getEnergyGrains(Emin, Emax, grainSize, Ngrains)

    def calculateDensitiesOfStates(self, Elist, E0):
        """
        Calculate and return an array containing the density of states for each
        isomer and reactant channel in the network. `Elist` represents the
        array of energies in J/mol at which to compute each density of states.
        The ground-state energies `E0` in J/mol are used to shift each density
        of states for each configuration to the same zero of energy. The 
        returned density of states is in units of mol/J.
        """
        
        Ngrains = len(Elist)
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        densStates = numpy.zeros((Nisom+Nreac, Ngrains), numpy.float64)
        dE = Elist[1] - Elist[0]
        
        logging.info('Calculating densities of states...')
        
        # Densities of states for isomers
        for i in range(Nisom):
            logging.info('Calculating density of states for isomer "%s"' % self.isomers[i])
            densStates0 = self.isomers[i].states.getDensityOfStates(Elist)
            # Shift to common zero of energy
            r0 = int(round(E0[i] / dE))
            densStates[i,r0:] = densStates0[:-r0+len(densStates0)]
        
        # Densities of states for reactant channels
        for n in range(Nreac):
            if self.reactants[n][0].states is not None and self.reactants[n][1].states is not None:
                logging.debug('Calculating density of states for reactant channel "%s"' % (' + '.join([str(spec) for spec in self.reactants[n]])))
                densStates0 = self.reactants[n][0].states.getDensityOfStates(Elist)
                densStates1 = self.reactants[n][1].states.getDensityOfStates(Elist)
                densStates0 = states.convolve(densStates0, densStates1, Elist)
                # Shift to common zero of energy
                r0 = int(round(E0[n+Nisom] / dE))
                densStates[n+Nisom,r0:] = densStates0[:-r0+len(densStates0)]
            else:
                logging.debug('NOT calculating density of states for reactant channel "%s"' % (' + '.join([str(spec) for spec in self.reactants[n]])))
        logging.debug('')
        
        return densStates

    def calculateMicrocanonicalRates(self, Elist, densStates):
        """
        Calculate and return arrays containing the microcanonical rate 
        coefficients :math:`k(E)` for the isomerization, dissociation, and
        association path reactions in the network. `Elist` represents the
        array of energies in J/mol at which to compute each density of states,
        while `densStates` represents the density of states of each isomer and
        reactant channel in mol/J.
        """
        
        Ngrains = len(Elist)
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        Kij = numpy.zeros([Nisom,Nisom,Ngrains], numpy.float64)
        Gnj = numpy.zeros([Nreac+Nprod,Nisom,Ngrains], numpy.float64)
        Fim = numpy.zeros([Nisom,Nreac,Ngrains], numpy.float64)
        
        logging.info('Calculating microcanonical rate coefficients k(E)...')
        
        for rxn in self.pathReactions:
            if rxn.reactants[0] in self.isomers and rxn.products[0] in self.isomers:
                # Isomerization
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.isomers.index(rxn.products[0])
                Kij[prod,reac,:], Kij[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], densStates[prod,:])
            elif rxn.reactants[0] in self.isomers and rxn.products in self.reactants:
                # Dissociation (reversible)
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.reactants.index(rxn.products)
                Gnj[prod,reac,:], Fim[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], densStates[prod+Nisom,:])
            elif rxn.reactants[0] in self.isomers and rxn.products in self.products:
                # Dissociation (irreversible)
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.products.index(rxn.products) + Nreac
                Gnj[prod,reac,:], dummy = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], None)
            elif rxn.reactants in self.reactants and rxn.products[0] in self.isomers:
                # Association
                reac = self.reactants.index(rxn.reactants)
                prod = self.isomers.index(rxn.products[0])
                Fim[prod,reac,:], Gnj[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac+Nisom,:], densStates[prod,:])
            else:
                raise NetworkError('Unexpected type of path reaction "%s"' % rxn)
        logging.debug('')
        
        return Kij, Gnj, Fim
        
    def calculateRateCoefficients(self, Tlist, Plist, Elist, method):
        """
        Calculate the phenomenological rate coefficients :math:`k(T,P)` for the
        network at the given temperatures `Tlist` in K and pressures `Plist` in
        Pa. The `method` string is used to indicate the method to use, and
        should be one of "modified strong collision", "reservoir state", or
        "chemically-significant eigenvalues".
        """

        # Determine the values of some counters
        Ngrains = len(Elist)
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        # Get ground-state energies of all isomers and each reactant channel
        # that has the necessary parameters
        # An exception will be raised if a unimolecular isomer is missing
        # this information
        E0 = numpy.zeros((Nisom+Nreac), numpy.float64)
        for i in range(Nisom):
            E0[i] = self.isomers[i].E0
        for n in range(Nreac):
            E0[n+Nisom] = sum([spec.E0 for spec in self.reactants[n]])
        
        # Get first reactive grain for each isomer
        Ereac = numpy.ones(Nisom, numpy.float64) * 1e20
        for i in range(Nisom):
            for rxn in self.pathReactions:
                if rxn.reactants[0] == self.isomers[i] or rxn.products[0] == self.isomers[i]:
                    if rxn.transitionState.E0 < Ereac[i]: 
                        Ereac[i] = rxn.transitionState.E0
        
        # Shift energy grains such that lowest is zero
        for rxn in self.pathReactions:
            rxn.transitionState.E0 -= Elist[0]
        E0 -= Elist[0]
        Ereac -= Elist[0]
        Elist -= Elist[0]

        # Calculate density of states for each isomer and each reactant channel
        # that has the necessary parameters
        densStates = self.calculateDensitiesOfStates(Elist, E0)
        
        # Calculate microcanonical rate coefficients for each path reaction
        # If degree of freedom data is provided for the transition state, then RRKM theory is used
        # If high-pressure limit Arrhenius data is provided, then the inverse Laplace transform method is used
        # Otherwise an exception is raised
        Kij, Gnj, Fim = self.calculateMicrocanonicalRates(Elist, densStates)

        K = numpy.zeros((len(Tlist),len(Plist),Nisom+Nreac+Nprod,Nisom+Nreac+Nprod), numpy.float64)
        
        for t, T in enumerate(Tlist):
            
            # Rescale densities of states such that, when they are integrated 
            # using the Boltzmann factor as a weighting factor, the result is unity
            for i in range(Nisom+Nreac):
                Q = numpy.sum(densStates[i,:] * numpy.exp(-Elist / constants.R / T))
                densStates[i,:] /= Q
        
            for p, P in enumerate(Plist):
                
                logging.info('Calculating k(T,P) values at %g K, %g bar...' % (T, P/1e5))
                
                # Calculate collision frequencies
                collFreq = numpy.zeros(Nisom, numpy.float64)
                for i in range(Nisom):
                    collFreq[i] = calculateCollisionFrequency(self.isomers[i], T, P, self.bathGas)
                
                # Apply method
                if method.lower() == 'modified strong collision':
                    # Modify collision frequencies using efficiency factor
                    for i in range(Nisom): 
                        collFreq[i] *= calculateCollisionEfficiency(self.isomers[i], T, Elist, densStates[i,:], self.collisionModel, E0[i], Ereac[i])
                    # Apply modified strong collision method
                    import msc
                    K[t,p,:,:], p0 = msc.applyModifiedStrongCollisionMethod(T, P, Elist, densStates, collFreq, Kij, Fim, Gnj, Ereac, Nisom, Nreac, Nprod)
                elif method.lower() == 'reservoir state':
                    # The full collision matrix for each isomer
                    Mcoll = numpy.zeros((Nisom,Ngrains,Ngrains), numpy.float64)
                    for i in range(Nisom):
                        Mcoll[i,:,:] = collFreq[i] * self.collisionModel.generateCollisionMatrix(Elist, T, densStates[i,:])
                    # Apply reservoir state method
                    import rs
                    K[t,p,:,:], p0 = rs.applyReservoirStateMethod(T, P, Elist, densStates, Mcoll, Kij, Fim, Gnj, Ereac, Nisom, Nreac, Nprod)
                elif method.lower() == 'chemically-significant eigenvalues':
                    # The full collision matrix for each isomer
                    Mcoll = numpy.zeros((Nisom,Ngrains,Ngrains), numpy.float64)
                    for i in range(Nisom):
                        Mcoll[i,:,:] = collFreq[i] * self.collisionModel.generateCollisionMatrix(Elist, T, densStates[i,:])
                    # Equilibrium ratios of each isomer and reactant channel (for symmetrization)
                    eqRatios = numpy.zeros(Nisom+Nreac, numpy.float64)
                    conc = 1e5 / constants.R / T
                    for i in range(Nisom):
                        G = self.isomers[i].getFreeEnergy(T)
                        eqRatios[i] = math.exp(-G / constants.R / T)
                    for n in range(Nreac):
                        G = sum([spec.getFreeEnergy(T) for spec in self.reactants[n]])
                        eqRatios[n+Nisom] = math.exp(-G / constants.R / T) * conc ** (len(self.reactants[n]) - 1)
                    eqRatios /= numpy.sum(eqRatios)
                    # Apply chemically-significant eigenvalues method
                    import cse
                    K[t,p,:,:], p0 = cse.applyChemicallySignificantEigenvaluesMethod(T, P, Elist, densStates, Mcoll, Kij, Fim, Gnj, eqRatios, Nisom, Nreac, Nprod)
                else:
                    raise NetworkError('Unknown method "%s".' % method)
                
                logging.debug('')
        
        return K
