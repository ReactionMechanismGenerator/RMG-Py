#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
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

import chempy.constants as constants
import chempy.states as states

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

    

################################################################################
