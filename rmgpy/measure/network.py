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
Contains the :class:`Network` class, which defines an internal representation
of a unimolecular reaction network. This provides a convienent means of
keeping track of things while we are performing the master equation calculation.
"""

import math
import numpy
import cython
import logging
import os.path

import rmgpy.constants as constants
from rmgpy.reaction import Reaction
import rmgpy.statmech as states
from rmgpy.reaction import Reaction

from reaction import *
from collision import *
import settings

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
    `index`             ``int``                 A unique integer identifier for the network
    `isomers`           ``list``                A list of the unimolecular isomers in the network
    `reactants`         ``list``                A list of the bimolecular reactant channels in the network
    `products`          ``list``                A list of the bimolecular product channels in the network
    `pathReactions`     ``list``                A list of reaction objects that connect adjacent isomers (the high-pressure-limit)
    `bathGas`           ``dict``                A dictionary of the bath gas species (keys) and their mole fractions (values)
    `collisionModel`    :class:`CollisionModel` The collision model to use
    `netReactions`      ``list``                A list of reaction objects that connect any pair of isomers
    `valid`             ``bool``                ``True`` if the rate coefficients for the network have been computed, ``False`` if not
    `title`             ``str``                 A short title for the network
    `description`       ``str``                 A verbose description of the network
    ------------------- ----------------------- --------------------------------
    `Tmin`              ``float``               The minimum temperature of interest in K
    `Tmax`              ``float``               The maximum temperature of interest in K
    `Pmin`              ``float``               The minimum pressure of interest in Pa
    `Pmax`              ``float``               The maximum pressure of interest in Pa
    `grainSize`         ``float``               The maximum grain size in J/mol
    `grainCount`        ``int``                 The minimum number of grains
    `E0`                :class:`numpy.ndarray`  The ground-state energies in J/mol for each configuration
    `Ereac`             :class:`numpy.ndarray`  The first reactive energy in J/mol for each isomer
    `Elist_full`        :class:`numpy.ndarray`  The energies used to compute the "full" densities of states
    `densStates_full`   :class:`numpy.ndarray`  The "full" densities of states for each configuration, referenced to the ground-state energy of that configuration
    `T`                 ``float``               The current temperature of interest in K
    `P`                 ``float``               The current pressure of interest in Pa
    `Elist`             :class:`numpy.ndarray`  The energy grains used at the current conditions
    `densStates`        :class:`numpy.ndarray`  The density of states for each configuration at the current conditions
    `Kij`               :class:`numpy.ndarray`  The microcanonical isomerization rate coefficients at the current conditions
    `Gnj`               :class:`numpy.ndarray`  The microcanonical dissociation rate coefficients at the current conditions
    `Fim`               :class:`numpy.ndarray`  The microcanonical association rate coefficients times equilibrium distributions at the current conditions
    `eqRatios`          :class:`numpy.ndarray`  The ratio of each configuration expected at equilibrium at the current conditions
    `dEdown`            :class:`numpy.ndarray`  The average energy transferred in a deactivating collision for each isomer at the current conditions
    `collFreq`          :class:`numpy.ndarray`  The collision frequency for each isomer at the current conditions
    `Mcoll`             :class:`numpy.ndarray`  The full collsion matrix for each isomer at the current conditions
    `K`                 :class:`numpy.ndarray`  The computed phenomenological rate coefficients at the current conditions
    `p0`                :class:`numpy.ndarray`  The computed reduced-space manifold at the current conditions
    =================== ======================= ================================

    """

    def __init__(self, index=-1, isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None, collisionModel=None, title='', description=''):
        self.index = index
        self.title = title
        self.description = description
        self.isomers = isomers or []
        self.reactants = reactants or []
        self.products = products or []
        self.pathReactions = pathReactions or []
        self.bathGas = bathGas or {}
        self.collisionModel = collisionModel
        self.netReactions = []
        self.valid = False
        self.errorString = ''
        
        self.clear()

    def __str__(self):
        if self.index != -1:
            return 'Network #{0:d}'.format(self.index)
        else:
            return 'Network "{0}"'.format(self.title)
    
    def clear(self):
        """
        Clear the state variables used in the pressure dependence calculation.
        """
        self.Tmin = 0
        self.Tmax = 0
        self.Pmin = 0
        self.Pmax = 0
        self.grainSize = 0.0
        self.grainCount = 0
        
        self.E0 = None
        self.Ereac = None
        self.Elist_full = None
        self.densStates_full = None
        
        self.T = 0
        self.P = 0
        
        self.Elist = None
        self.densStates = None
        self.Kij = None
        self.Gnj = None
        self.Fim = None
        self.eqRatios = None
        self.dEdown = None
        
        self.collFreq = None
        self.Mcoll = None
        
        self.K = None
        self.p0 = None
    
    def invalidate(self):
        """
        Mark the network as in need of a new calculation to determine the
        pressure-dependent rate coefficients
        """
        self.valid = False

    def merge(self, other):
        """
        Merge another :class:`Network` object `other` into this object.
        """
        self.isomers.extend(other.isomers)
        self.pathReactions.extend(other.pathReactions)
        self.explored.extend(other.explored)
        self.invalidate()

    def containsSpecies(self, species):
        """
        Return :data:`True` if the `species` is a *unimolecular* isomer in the
        network, and :data:`False` if not.
        """
        for rxn in self.pathReactions:
            if len(rxn.reactants) == 1:
                if rxn.reactants[0] == species and species in self.explored: return True
            if len(rxn.products) == 1:
                if rxn.products[0] == species and species in self.explored: return True
        return False

    def getAllSpecies(self):
        """
        Return a list of all unique species in the network, including all
        isomers, reactant and product channels, and bath gas species.
        """
        speciesList = []
        for isomer in self.isomers:
            if isomer not in speciesList: speciesList.append(isomer)
        for reactants in self.reactants:
            for spec in reactants:
                if spec not in speciesList: speciesList.append(spec)
        for products in self.products:
            for spec in products:
                if spec not in speciesList: speciesList.append(spec)
        for spec in self.bathGas:
            if spec not in speciesList: speciesList.append(spec)
        return speciesList
    
    def deleteSpecies(self, species):
        """
        Completely remove the given `species` from the network, including any
        path and net reactions it participates in. 
        """
        # Remove the species from the list of isomers
        if species in self.isomers:
            self.isomers.remove(species)
        
        # Remove any reactant channels containing the species
        reactantsToRemove = []
        for reactants in self.reactants:
            if species in reactants:
                reactantsToRemove.append(reactants)
        for reactants in reactantsToRemove:
            self.reactants.remove(reactants)
            
        # Remove any product channels containing the species
        productsToRemove = []
        for products in self.products:
            if species in products:
                productsToRemove.append(products)
        for products in productsToRemove:
            self.products.remove(products)
        
        # Remove any path reactions containing the species
        pathReactionsToRemove = []
        for reaction in self.pathReactions:
            if species in reaction.reactants or species in reaction.products:
                pathReactionsToRemove.append(reaction)
        for reaction in pathReactionsToRemove:
            self.pathReactions.remove(reaction)
        
        # Remove any net reactions containing the species
        netReactionsToRemove = []
        for reaction in self.netReactions:
            if species in reaction.reactants or species in reaction.products:
                netReactionsToRemove.append(reaction)
        for reaction in netReactionsToRemove:
            self.netReactions.remove(reaction)
    
    def deletePathReaction(self, reaction):
        """
        Remove the given path `reaction` from the network. Does not remove the
        species involved in the reaction.
        """
        self.pathReactions.remove(reaction)
    
    def printSummary(self, level=logging.DEBUG):
        """
        Print a formatted list of information about the current network. Each
        molecular configuration - unimolecular isomers, bimolecular reactant
        channels, and bimolecular product channels - is given along with its
        energy on the potential energy surface. The path reactions connecting
        adjacent molecular configurations are also given, along with their
        energies on the potential energy surface. The `level` parameter controls
        the level of logging to which the summary is written, and is DEBUG by
        default.
        """

        logging.log(level, '')
        logging.log(level, '========================================================================')
        logging.log(level, 'Network Information (Network #{0})'.format(self.index))
        logging.log(level, '-------------------')
        logging.log(level, 'Isomers:')
        for isomer in self.isomers:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(str(isomer), isomer.E0.value_si / 1000.0))
        logging.log(level, 'Reactant channels:')
        for reactants in self.reactants:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(' + '.join([str(spec) for spec in reactants]), sum([spec.E0.value_si for spec in reactants]) / 1000.0))
        logging.log(level, 'Product channels:')
        for products in self.products:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(' + '.join([str(spec) for spec in products]), sum([spec.E0.value_si for spec in products]) / 1000.0))
        logging.log(level, 'Path reactions:')
        for rxn in self.pathReactions:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(rxn, rxn.transitionState.E0.value_si / 1000.0))
        logging.log(level, '========================================================================')
        logging.log(level, '')

    def __getEnergyGrains(self, Emin, Emax, grainSize=0.0, grainCount=0):
        """
        Return an array of energy grains that have a minimum of `Emin`, a
        maximum of `Emax`, and either a spacing of `grainSize` or have number of
        grains `grainCount`. The first three parameters are in J/mol, as is the
        returned array of energy grains.
        """
        
        # Now determine the grain size and number of grains to use
        if grainCount <= 0 and grainSize <= 0.0:
            # Neither grain size nor number of grains specified, so raise exception
            raise NetworkError('You must specify a positive value for either dE or Ngrains.')
        elif grainCount <= 0 and grainSize > 0.0:
            # Only grain size was specified, so we must use it
            useGrainSize = True
        elif grainCount > 0 and grainSize <= 0.0:
            # Only number of grains was specified, so we must use it
            useGrainSize = False
        else:
            # Both were specified, so we choose the tighter constraint
            # (i.e. the one that will give more grains, and so better accuracy)
            grainSize0 = (Emax - Emin) / (grainCount - 1)
            useGrainSize = (grainSize0 > grainSize)

        # Generate the array of energies
        if useGrainSize:
            Elist = numpy.arange(Emin, Emax + grainSize, grainSize, numpy.float64)
        else:
            Elist = numpy.linspace(Emin, Emax, grainCount, numpy.float64)

        return Elist

    def selectEnergyGrains(self, T, grainSize=0.0, grainCount=0):
        """
        Select a suitable list of energies to use for subsequent calculations.
        This is done by finding the minimum and maximum energies on the 
        potential energy surface, then adding a multiple of 
        :math:`k_\\mathrm{B} T` onto the maximum energy.

        You must specify either the desired grain spacing `grainSize` in J/mol
        or the desired number of grains `Ngrains`, as well as a temperature
        `T` in K to use for the equilibrium calculation. You can specify both
        `grainSize` and `grainCount`, in which case the one that gives the more
        accurate result will be used (i.e. they represent a maximum grain size
        and a minimum number of grains). An array containing the energy grains
        in J/mol is returned.
        """

        if grainSize == 0.0 and grainCount == 0:
            raise NetworkError('Must provide either grainSize or Ngrains parameter to Network.determineEnergyGrains().')

        # The minimum energy is the lowest isomer or reactant or product energy on the PES
        Emin = 1.0e25
        for isomer in self.isomers:
            if isomer.E0.value_si < Emin: Emin = isomer.E0.value_si
        for reactants in self.reactants:
            E0 = sum([reactant.E0.value_si for reactant in reactants])
            if E0 < Emin: Emin = E0
        for products in self.products:
            E0 = sum([product.E0.value_si for product in products])
            if E0 < Emin: Emin = E0
        Emin = math.floor(Emin) # Round to nearest whole number

        # Use the highest energy on the PES as the initial guess for Emax0
        Emax = -1.0e25
        for isomer in self.isomers:
            if isomer.E0.value_si > Emax: Emax = isomer.E0.value_si
        for reactants in self.reactants:
            E0 = sum([reactant.E0.value_si for reactant in reactants])
            if E0 > Emax: Emax = E0
        for products in self.products:
            E0 = sum([product.E0.value_si for product in products])
            if E0 > Emax: Emax = E0
        for rxn in self.pathReactions:
            if rxn.transitionState is not None:
                E0 = rxn.transitionState.E0.value_si
                if E0 > Emax: Emax = E0
        
        # Choose the actual Emax as many kB * T above the maximum energy on the PES
        # You should check that this is high enough so that the Boltzmann distributions have trailed off to negligible values
        Emax += 40 * constants.R * T

        return self.__getEnergyGrains(Emin, Emax, grainSize, grainCount)

    def calculateDensitiesOfStates(self):
        """
        Calculate the densities of states of each configuration that has states
        data. The densities of states are computed such that they can be
        applied to each temperature in the range of interest by interpolation.
        """
        
        Tmin = self.Tmin
        Tmax = self.Tmax
        grainSize = self.grainSize
        grainCount = self.grainCount
        
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        logging.info('Calculating densities of states for {0}...'.format(self))

        # Choose the energies used to compute the densities of states
        # Use Tmin to select the minimum energy and grain size
        Elist0 = self.selectEnergyGrains(Tmin, grainSize, grainCount)
        Emin0 = numpy.min(Elist0)
        grainSize0 = Elist0[1] - Elist0[0]
        # Use Tmax to select the maximum energy and grain count
        Elist0 = self.selectEnergyGrains(Tmax, grainSize, grainCount)
        grainCount0 = len(Elist0)
        Emax0 = numpy.max(Elist0)
        
        Elist = self.__getEnergyGrains(Emin0, Emax0, grainSize0, grainCount0)
        Ngrains = len(Elist)
        dE = Elist[1] - Elist[0]
        logging.info('Using {0:d} grains from {1:.2f} to {2:.2f} kJ/mol in steps of {3:.2f} kJ/mol to compute densities of states'.format(
            Ngrains, Elist[0] / 1000, Elist[-1] / 1000, dE / 1000))
        
        # Shift the energy grains so that the minimum grain is zero
        Elist -= Elist[0]
        
        # Actually compute the densities of states
        densStates = numpy.zeros((Nisom+Nreac+Nprod, Ngrains), numpy.float64)
        # Densities of states for isomers
        for i in range(Nisom):
            logging.debug('Calculating density of states for isomer "{0}"'.format(self.isomers[i]))
            densStates[i,:] = self.isomers[i].states.getDensityOfStates(Elist)
        # Densities of states for reactant channels
        for n in range(Nreac):
            if self.reactants[n][0].states is not None and self.reactants[n][1].states is not None:
                logging.debug('Calculating density of states for reactant channel "{0}"'.format(' + '.join([str(spec) for spec in self.reactants[n]])))
                densStates0 = self.reactants[n][0].states.getDensityOfStates(Elist)
                densStates1 = self.reactants[n][1].states.getDensityOfStates(Elist)
                densStates[n+Nisom,:] = states.convolve(densStates0, densStates1, Elist)
            elif self.reactants[n][0].states is not None:
                logging.debug('Calculating density of states for reactant channel "{0}"'.format(' + '.join([str(spec) for spec in self.reactants[n]])))
                densStates[n+Nisom,:] = self.reactants[n][0].states.getDensityOfStates(Elist)
            elif self.reactants[n][1].states is not None:
                logging.debug('Calculating density of states for reactant channel "{0}"'.format(' + '.join([str(spec) for spec in self.reactants[n]])))
                densStates[n+Nisom,:] = self.reactants[n][1].states.getDensityOfStates(Elist)
            else:
                logging.debug('NOT calculating density of states for reactant channel "{0}"'.format(' + '.join([str(spec) for spec in self.reactants[n]])))
        # Densities of states for product channels
        for n in range(Nprod):
            if len(self.products[n]) == 1 and self.products[n][0].states is not None:
                logging.debug('Calculating density of states for product channel "{0}"'.format(' + '.join([str(spec) for spec in self.products[n]])))
                densStates[n+Nisom+Nreac,:] = self.isomers[i].states.getDensityOfStates(Elist)
            elif len(self.products[n]) == 2 and self.products[n][0].states is not None and self.products[n][1].states is not None:
                logging.debug('Calculating density of states for product channel "{0}"'.format(' + '.join([str(spec) for spec in self.products[n]])))
                densStates0 = self.products[n][0].states.getDensityOfStates(Elist)
                densStates1 = self.products[n][1].states.getDensityOfStates(Elist)
                densStates[n+Nisom+Nreac,:] = states.convolve(densStates0, densStates1, Elist)
            elif len(self.products[n]) == 2 and self.products[n][0].states is not None:
                logging.debug('Calculating density of states for product channel "{0}"'.format(' + '.join([str(spec) for spec in self.products[n]])))
                densStates[n+Nisom+Nreac,:] = self.products[n][0].states.getDensityOfStates(Elist)
            elif len(self.products[n]) == 2 and self.products[n][1].states is not None:
                logging.debug('Calculating density of states for product channel "{0}"'.format(' + '.join([str(spec) for spec in self.products[n]])))
                densStates[n+Nisom+Nreac,:] = self.products[n][1].states.getDensityOfStates(Elist)
            else:
                logging.debug('NOT calculating density of states for product channel "{0}"'.format(' + '.join([str(spec) for spec in self.products[n]])))
        logging.debug('')
        
        # Store the computed densities of states on this object
        self.Elist_full = Elist
        self.densStates_full = densStates
  
        return densStates
    
    def calculateGroundStateEnergies(self):
        """
        Return an array of the ground-state energies in J/mol of each isomer,
        reactant channel, and product channel, in that order.
        """
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        E0 = numpy.zeros((Nisom+Nreac+Nprod), numpy.float64)
        for i in range(Nisom):
            E0[i] = self.isomers[i].E0.value_si
        for n in range(Nreac):
            E0[n+Nisom] = sum([spec.E0.value_si for spec in self.reactants[n]])
        for n in range(Nprod):
            E0[n+Nisom+Nreac] = sum([spec.E0.value_si for spec in self.products[n]])
        self.E0 = E0
        return E0
    
    def calculateFirstReactiveEnergies(self):
        """
        Return an array containing the lowest reactive energy for each isomer.
        """
        Nisom = len(self.isomers)
        Ereac = numpy.ones(Nisom, numpy.float64) * 1e20
        for i in range(Nisom):
            for rxn in self.pathReactions:
                if rxn.reactants[0] == self.isomers[i] or rxn.products[0] == self.isomers[i]:
                    if rxn.transitionState.E0.value_si < Ereac[i]:
                        Ereac[i] = rxn.transitionState.E0.value_si
        self.Ereac = Ereac
        return Ereac

    def mapDensitiesOfStates(self):
        """
        Map the overall densities of states to the current energy grains.
        Semi-logarithmic interpolation will be used if the grain sizes of 
        `Elist0` and `Elist` do not match; this should not be a significant
        source of error as long as the grain sizes are sufficiently small.
        """
        from _network import mapDensitiesOfStates
        self.densStates = mapDensitiesOfStates(
            self.Elist,
            self.Elist_full,
            self.densStates_full,
            self.T,
            self.E0,
            len(self.isomers), 
            len(self.reactants), 
            len(self.products),
        )
        return self.densStates

    def calculateMicrocanonicalRates(self, check=True):
        """
        Calculate and return arrays containing the microcanonical rate
        coefficients :math:`k(E)` for the isomerization, dissociation, and
        association path reactions in the network. If `check` is ``True``, the
        computed :math:`k(E)` values are checked that they integrate to give
        the expected :math:`k(T)` values (with some allowance for numerical
        error).
        """

        T = self.T
        Elist = self.Elist
        densStates = self.densStates
        Ngrains = len(Elist)
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)

        Kij = numpy.zeros([Nisom,Nisom,Ngrains], numpy.float64)
        Gnj = numpy.zeros([Nreac+Nprod,Nisom,Ngrains], numpy.float64)
        Fim = numpy.zeros([Nisom,Nreac,Ngrains], numpy.float64)

        for rxn in self.pathReactions:
            kf0 = 0.0; kr0 = 0.0; kf = 0.0; kr = 0.0
            if rxn.reactants[0] in self.isomers and rxn.products[0] in self.isomers:
                # Isomerization
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.isomers.index(rxn.products[0])
                Kij[prod,reac,:], Kij[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], densStates[prod,:], T)
                if check:
                    kf0 = numpy.sum(Kij[prod,reac,:] * densStates[reac,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[reac,:] * numpy.exp(-Elist / constants.R / T))
                    kr0 = numpy.sum(Kij[reac,prod,:] * densStates[prod,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[prod,:] * numpy.exp(-Elist / constants.R / T))
            elif rxn.reactants[0] in self.isomers and rxn.products in self.reactants:
                # Dissociation (reversible)
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.reactants.index(rxn.products)
                Gnj[prod,reac,:], Fim[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], densStates[prod+Nisom,:], T)
                if check:
                    kf0 = numpy.sum(Gnj[prod,reac,:] * densStates[reac,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[reac,:] * numpy.exp(-Elist / constants.R / T))
                    kr0 = numpy.sum(Fim[reac,prod,:])
            elif rxn.reactants[0] in self.isomers and rxn.products in self.products:
                # Dissociation (irreversible)
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.products.index(rxn.products) + Nreac
                Gnj[prod,reac,:], dummy = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], densStates[prod+Nisom,:], T)
                if check:
                    kf0 = numpy.sum(Gnj[prod,reac,:] * densStates[reac,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[reac,:] * numpy.exp(-Elist / constants.R / T))
                    if rxn.isIsomerization():
                        if not densStates[prod+Nisom,:].any():
                            kr0 = 0.0
                        else:
                            kr0 = numpy.sum(dummy * densStates[prod+Nisom,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[prod+Nisom,:] * numpy.exp(-Elist / constants.R / T))
                    else:
                        kr0 = numpy.sum(dummy)
            elif rxn.reactants in self.reactants and rxn.products[0] in self.isomers:
                # Association (reversible)
                reac = self.reactants.index(rxn.reactants)
                prod = self.isomers.index(rxn.products[0])
                Fim[prod,reac,:], Gnj[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac+Nisom,:], densStates[prod,:], T)
                if check:
                    kf0 = numpy.sum(Fim[prod,reac,:])
                    kr0 = numpy.sum(Gnj[reac,prod,:] * densStates[prod,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[prod,:] * numpy.exp(-Elist / constants.R / T))
            elif rxn.reactants in self.products and rxn.products[0] in self.isomers:
                # Association (irreversible)
                reac = self.products.index(rxn.reactants) + Nreac
                prod = self.isomers.index(rxn.products[0])
                dummy, Gnj[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac+Nisom,:], densStates[prod,:], T)
                if check:
                    if rxn.isIsomerization():
                        if not densStates[reac+Nisom,:].any():
                            kf0 = 0.0
                        else:
                            kf0 = numpy.sum(dummy * densStates[reac+Nisom,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[reac+Nisom,:] * numpy.exp(-Elist / constants.R / T))
                    else:
                        kf0 = numpy.sum(dummy)
                    kr0 = numpy.sum(Gnj[reac,prod,:] * densStates[prod,:] * numpy.exp(-Elist / constants.R / T)) / numpy.sum(densStates[prod,:] * numpy.exp(-Elist / constants.R / T))
            else:
                raise NetworkError('Unexpected type of path reaction "{0}"'.format(rxn))
        
            if check and rxn.kinetics is not None:
                # Check that the computed k(E) values integrate to give the 
                # expected k(T) values (allowing for some numerical error)
                error = False
                kf = rxn.kinetics.getRateCoefficient(T)
                if kr0 != 0:
                    # Check both kf(T) and Keq(T)
                    Keq0 = kf0 / kr0
                    Keq = rxn.getEquilibriumConstant(T)
                    if kf0 / kf < 0.1 or kf0 / kf > 10.0 or Keq0 / Keq < 0.5 or Keq0 / Keq > 2.0:
                        error = True
                        logging.error('For path reaction {0!s}:'.format(rxn))
                        logging.error('    Expected kf({0:g} K) = {1:g}'.format(T, kf))
                        logging.error('      Actual kf({0:g} K) = {1:g}'.format(T, kf0))
                        logging.error('    Expected Keq({0:g} K) = {1:g}'.format(T, Keq))
                        logging.error('      Actual Keq({0:g} K) = {1:g}'.format(T, Keq0))
                else:
                    # Check only kf(T)
                    if kf0 / kf < 0.1 or kf0 / kf > 10.0:
                        error = True
                        logging.error('For path reaction {0!s}:'.format(rxn))
                        logging.error('    Expected kf({0:g} K) = {1:g}'.format(T, kf))
                        logging.error('      Actual kf({0:g} K) = {1:g}'.format(T, kf0))
                # If the k(E) values are invalid (in that they give the wrong 
                # kf(T) or kr(T) when integrated), then raise an exception
                if error:
                    raise NetworkError('Invalid k(E) values computed for path reaction "{0}" do not satisfy'.format(rxn))
                
        # In the past, we have occasionally encountered k(E) values that are NaN
        # Just to be safe, let's check to be sure this isn't happening with this network
        if numpy.isnan(Kij).any() or numpy.isnan(Gnj).any() or numpy.isnan(Fim).any():
            raise NetworkError('One or more k(E) values is NaN for network {0:d}.'.format(self.index))
        
        self.Kij = Kij
        self.Gnj = Gnj
        self.Fim = Fim
        
        return Kij, Gnj, Fim

    def calculateDeltaEDown(self):
        """
        Return an array containing the average energy transferred in a 
        deactivating collision for each isomer at the current temperature. 
        First, the bath gas value is computed as a weighted sum of the values
        for each bath gas component. This value is then averaged arithmetically 
        with that from the isomer itself if the latter is provided.
        """
        T = self.T
        Nisom = len(self.isomers)
        dEdown = numpy.zeros(Nisom, numpy.float64)
        for i in range(Nisom):
            # First compute dEdown as a weighted sum of the values from each of the bath gas components
            totalFrac = 0
            for species, frac in self.bathGas.iteritems():
                if species.collisionModel is not None:
                    dEdown[i] += frac * species.collisionModel.getAlpha(T)
                    totalFrac += frac
            dEdown[i] /= totalFrac
            # If the isomer also has a collision model, then average its dEdown with that of the bath gas
            if self.isomers[i].collisionModel is not None:
                dEdown[i] = 0.5 * (dEdown[i] + self.isomers[i].collisionModel.getAlpha(T))
        self.dEdown = dEdown
        return dEdown
    
    def calculateEquilibriumRatios(self):
        """
        Return an array containing the fraction of each isomer and reactant
        channel present at equilibrium, as determined from the Gibbs free 
        energy and using the concentration equilibrium constant 
        :math:`K_\\mathrm{c}`. These values are ratios, and the absolute
        magnitude is not guaranteed; however, the implementation scales the 
        elements of the array so that they sum to unity.
        """
        T = self.T
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        eqRatios = numpy.zeros(Nisom+Nreac, numpy.float64)
        conc = 1e5 / constants.R / T
        for i in range(Nisom):
            G = self.isomers[i].thermo.getFreeEnergy(T)
            eqRatios[i] = math.exp(-G / constants.R / T)
        for i in range(Nreac):
            G = sum([spec.thermo.getFreeEnergy(T) for spec in self.reactants[i]])
            eqRatios[Nisom+i] = math.exp(-G / constants.R / T) * conc ** (len(self.reactants[i]) - 1)
        self.eqRatios = eqRatios
        return eqRatios / numpy.sum(eqRatios)

    def calculateCollisionModel(self):
        """
        Calculate the matrix of first-order rate coefficients for collisional
        population transfer between grains for each isomer, including the
        corresponding collision frequencies.
        """
        T = self.T
        P = self.P
        Elist = self.Elist
        densStates = self.densStates
        dEdown = self.dEdown
        Nisom = len(self.isomers)
        Ngrains = len(Elist)
        
        # Calculate collision frequencies
        collFreq = numpy.zeros(Nisom, numpy.float64)
        for i in range(Nisom):
            collFreq[i] = calculateCollisionFrequency(self.isomers[i], T, P, self.bathGas)
        self.collFreq = collFreq
        
        # Generate the full collision matrix for each isomer
        Mcoll = numpy.zeros((Nisom,Ngrains,Ngrains), numpy.float64)
        for i in range(Nisom):
            Mcoll[i,:,:] = collFreq[i] * SingleExponentialDown().generateCollisionMatrix(Elist, T, densStates[i,:], dEdown[i])
        self.Mcoll = Mcoll
        
        return Mcoll
    
    def initialize(self, Tmin, Tmax, Pmin, Pmax, grainSize=0.0, grainCount=0):
        """
        Initialize a pressure dependence calculation by computing several
        quantities that are independent of the conditions. You must specify
        the temperature and pressure ranges of interesting using `Tmin` and
        `Tmax` in K and `Pmin` and `Pmax` in Pa. You must also specify the
        maximum energy grain size `grainSize` in J/mol and/or the minimum
        number of grains `grainCount`.
        """
        if grainSize == 0.0 and grainCount == 0:
            raise NetworkError('Must provide either grainSize or Ngrains parameter to Network.determineEnergyGrains().')

        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.grainSize = grainSize
        self.grainCount = grainCount
        
        self.calculateGroundStateEnergies() 
        self.calculateFirstReactiveEnergies()
        self.calculateDensitiesOfStates()
        
    def setConditions(self, T, P, check=True):
        """
        Set the current network conditions to the temperature `T` in K and
        pressure `P` in Pa. All of the internal variables are updated 
        accordingly if they are out of date. For example, those variables that
        depend only on temperature will not be recomputed if the temperature
        is the same.
        """
        
        temperatureChanged = (self.T != T)
        pressureChanged = (self.P != P)
        self.T = T
        self.P = P
        
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        E0 = self.E0
        Ereac = self.Ereac
        Elist0 = self.Elist_full
        densStates0 = self.densStates_full
        grainSize = self.grainSize
        grainCount = self.grainCount
        
        # Update parameters that depend on temperature only if necessary
        if temperatureChanged:
            # Choose the energy grains to use to compute k(T,P) values at this temperature
            Elist = self.Elist = self.selectEnergyGrains(T, grainSize, grainCount)
            Ngrains = len(Elist)
            logging.info('Using {0:d} grains from {1:.2f} to {2:.2f} kJ/mol in steps of {3:.2f} kJ/mol to compute the k(T,P) values at {4:g} K'.format(
                Ngrains, numpy.min(Elist) / 1000, numpy.max(Elist) / 1000, (Elist[1] - Elist[0]) / 1000, T))
            # Map the densities of states onto this set of energies
            # Also shift each density of states to a common zero of energy
            densStates = self.mapDensitiesOfStates()
            # Calculate microcanonical rate coefficients for each path reaction
            self.calculateMicrocanonicalRates(check=check)
            # Rescale densities of states such that, when they are integrated
            # using the Boltzmann factor as a weighting factor, the result is unity
            for i in range(Nisom+Nreac):
                densStates[i,:] = densStates[i,:] / numpy.sum(densStates[i,:] * numpy.exp(-Elist / constants.R / T))
            # Use free energy to determine equilibrium ratios of each isomer and product channel
            self.calculateEquilibriumRatios()
            # Compute average energy transferred in a deactivating collision
            self.calculateDeltaEDown()
            
        # Update parameters that depend on temperature and pressure if necessary
        if temperatureChanged or pressureChanged:
            self.calculateCollisionModel()

    def applyModifiedStrongCollisionMethod(self, efficiencyModel='default'):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the modified strong collision method.
        """
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        import msc
        logging.debug('Applying modified strong collision method at {0:g} K, {1:g} bar...'.format(self.T, self.P/1e5))
        self.K, self.p0 = msc.applyModifiedStrongCollisionMethod(self.T, self.P, self.Elist, self.densStates, self.collFreq, self.dEdown, self.Kij, self.Fim, self.Gnj, self.E0, self.Ereac, efficiencyModel, Nisom, Nreac, Nprod)
        return self.K, self.p0
    
    def applyReservoirStateMethod(self):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the reservoir state method.
        """
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        import rs
        logging.debug('Applying reservoir state method at {0:g} K, {1:g} bar...'.format(self.T, self.P/1e5))
        self.K, self.p0 = rs.applyReservoirStateMethod(self.T, self.P, self.Elist, self.densStates, self.Mcoll, self.Kij, self.Fim, self.Gnj, self.Ereac, Nisom, Nreac, Nprod)
        return self.K, self.p0
    
    def applyChemicallySignificantEigenvaluesMethod(self):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the chemically-significant eigenvalues method.
        """
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        import cse
        logging.debug('Applying chemically-significant eigenvalues method at {0:g} K, {1:g} bar...'.format(self.T, self.P/1e5))
        self.K, self.p0 = cse.applyChemicallySignificantEigenvaluesMethod(self.T, self.P, self.Elist, self.densStates, self.Mcoll, self.Kij, self.Fim, self.Gnj, self.eqRatios, Nisom, Nreac, Nprod)
        return self.K, self.p0
    
    def calculateRateCoefficients(self, Tlist, Plist, method, grainSize=None, grainCount=None, errorCheck=True):
        """
        Calculate the phenomenological rate coefficients :math:`k(T,P)` for the
        network at the given temperatures `Tlist` in K and pressures `Plist` in
        Pa using the energy grains ``Elist`` in J/mol. The `method` string is
        used to indicate the method to use, and should be one of ``"modified
        strong collision"``, ``"reservoir state"``, or
        ``"chemically-significant eigenvalues"``.
        """

        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        self.initialize(
            Tmin = numpy.min(Tlist),
            Tmax = numpy.max(Tlist),
            Pmin = numpy.min(Plist),
            Pmax = numpy.max(Plist),
            grainSize = grainSize,
            grainCount = grainCount,
        )
        
        for rxn in self.pathReactions:
            if rxn.transitionState.states is not None:
                logging.debug('Using RRKM theory to compute k(E) for path reaction {0}.'.format(rxn))
            elif rxn.kinetics is not None:
                logging.debug('Using ILT method to compute k(E) for path reaction {0}.'.format(rxn))
        logging.debug('')
        
        logging.info('Calculating phenomenological rate coefficients for {0}...'.format(self))
        K = numpy.zeros((len(Tlist),len(Plist),Nisom+Nreac+Nprod,Nisom+Nreac+Nprod), numpy.float64)
        
        for t, T in enumerate(Tlist):
            for p, P in enumerate(Plist):
                self.setConditions(T, P)
                
                # Apply method
                if method.lower() == 'modified strong collision':
                    self.applyModifiedStrongCollisionMethod()
                elif method.lower() == 'reservoir state':
                    self.applyReservoirStateMethod()
                elif method.lower() == 'chemically-significant eigenvalues':
                    self.applyChemicallySignificantEigenvaluesMethod()
                else:
                    raise NetworkError('Unknown method "{0}".'.format(method))

                K[t,p,:,:] = self.K

                # Check that the k(T,P) values satisfy macroscopic equilibrium
                eqRatios = self.eqRatios
                for i in range(Nisom+Nreac):
                    for j in range(i):
                        Keq0 = K[t,p,i,j] / K[t,p,j,i]
                        Keq = eqRatios[i] / eqRatios[j]
                        if Keq0 / Keq < 0.5 or Keq0 / Keq > 2.0:
                            if i < Nisom:
                                reactants = [self.isomers[i]]
                            elif i < Nisom+Nreac:
                                reactants = self.reactants[i-Nisom]
                            else:
                                reactants = self.products[i-Nisom-Nreac]
                            if j < Nisom:
                                products = [self.isomers[j]]
                            elif j < Nisom+Nreac:
                                products = self.reactants[j-Nisom]
                            else:
                                products = self.products[j-Nisom-Nreac]
                            reaction = Reaction(reactants=reactants, products=products)
                            logging.error('For net reaction {0!s}:'.format(reaction))
                            logging.error('Expected Keq({1:g} K, {2:g} bar) = {0:11.3e}'.format(Keq, T, P/1e5))
                            logging.error('  Actual Keq({1:g} K, {2:g} bar) = {0:11.3e}'.format(Keq0, T, P/1e5))
                            raise NetworkError('MEASURE computed k(T,P) values for reaction {0!s} do not satisfy macroscopic equilibrium.'.format(reaction))
                            
                # Compute k(T,P) values from the return p0
                # This should be identical to the k(T,P) values returned by each method
                #if method.lower() != 'modified strong collision':
                #    import me
                #    K[t,p,:,:] = me.computeRateCoefficients(Mcoll, Kij, Fim, Gnj, p0[t,p,:,:,:], Nisom, Nreac, Nprod)

                # Reject if any rate coefficients are negative
                if errorCheck:
                    negativeRate = False
                    for i in range(Nisom+Nreac+Nprod):
                        for j in range(i):
                            if (K[t,p,i,j] < 0 or K[t,p,j,i] < 0) and not negativeRate:
                                negativeRate = True
                                logging.error('Negative rate coefficient generated; rejecting result.')
                                logging.info(K[t,p,0:Nisom+Nreac+Nprod,0:Nisom+Nreac])
                                K[t,p,:,:] = 0 * K[t,p,:,:]
                                p0[t,p,:,:,:] = 0 * p0[t,p,:,:,:]
                            #elif K[t,p,i,j] < 0 and i < Nisom+Nreac and j < Nisom+Nreac:
                                #K[t,p,i,j] = K[t,p,j,i] * eqRatios[j] / eqRatios[i]
                            #elif K[t,p,j,i] < 0 and i < Nisom+Nreac and j < Nisom+Nreac:
                                #K[t,p,j,i] = K[t,p,i,j] * eqRatios[i] / eqRatios[j]
                                
                            
                logging.log(0, K[t,p,0:Nisom+Nreac+Nprod,0:Nisom+Nreac])
                
        logging.debug('')

        # Mark network as valid
        self.valid = True

        return K

    def generateFullMEMatrix(self, T, P, grainSize=None, grainCount=None):
        """
        Generate the full master equation matrix for the network at the
        specified temperature `T` in K and pressure `P` in Pa and for the
        set of energy grains `Elist` in J/mol. Returns the full master equation
        matrix, the indices used to place the isomers into the rows of the
        matrix, and the computed densities of states in mol/J.
        """

        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        
        self.initialize(Tmin=T, Tmax=T, Pmin=P, Pmax=P, grainSize=grainSize, grainCount=grainCount)
        self.setConditions(T, P)
        
        Mcoll = self.Mcoll
        Kij = self.Kij
        Gnj = self.Gnj
        Fim = self.Fim
        Elist = self.Elist
        densStates = self.densStates
        Ngrains = len(Elist)
        
        # Construct accounting matrix
        # Row is grain number, column is well number, value is index into matrix
        indices = -numpy.ones([Ngrains,Nisom], numpy.int)
        Nrows = 0
        for r in range(Ngrains):
            for i in range(Nisom):
                if densStates[i,r] > 0:
                    indices[r,i] = Nrows
                    Nrows += 1
        Nrows += Nreac + Nprod
        
        # Construct full ME matrix
        M = numpy.zeros([Nrows,Nrows], numpy.float64)
        # Collision terms
        for i in range(Nisom):
            for r in range(Ngrains):
                if indices[r,i] > -1:
                    for s in range(r, Ngrains):
                        M[indices[r,i], indices[s,i]] = Mcoll[i,r,s]
                        M[indices[s,i], indices[r,i]] = Mcoll[i,s,r]
        # Isomerization terms
        for i in range(Nisom):
            for j in range(i):
                if Kij[i,j,Ngrains-1] > 0 or Kij[j,i,Ngrains-1] > 0:
                    for r in range(Ngrains):
                        u = indices[r,i]; v = indices[r,j]
                        if u > -1 and v > -1:
                            M[v,u] = Kij[j,i,r]
                            M[u,u] -= Kij[j,i,r]
                            M[u,v] = Kij[i,j,r]
                            M[v,v] -= Kij[i,j,r]
        # Association/dissociation terms
        for i in range(Nisom):
            for n in range(Nreac+Nprod):
                if Gnj[n,i,Ngrains-1] > 0:
                    for r in range(Ngrains):
                        u = indices[r,i]; v = Nrows - Nreac - Nprod + n
                        if u > -1:
                            M[u,u] -= Gnj[n,i,r]
                            M[v,u] = Gnj[n,i,r]
                            if n < Nreac:
                                M[u,v] = Fim[i,n,r]
                                M[v,v] -= Fim[i,n,r]

        return M, Elist, indices, densStates

    def __createNewSurfaceAndContext(self, ext, fstr='', width=800, height=600):
        import cairo
        if ext == '.svg':
            surface = cairo.SVGSurface(fstr, width, height)
        elif ext == '.pdf':
            surface = cairo.PDFSurface(fstr, width, height)
        elif ext == '.ps':
            surface = cairo.PSSurface(fstr, width, height)
        elif ext == '.png':
            surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(width), int(height))
        else:
            logging.warning('Unknown format for target "{0}"; not drawing potential energy surface.'.format(fstr))
            return
        cr = cairo.Context(surface)
        return surface, cr

    def drawPotentialEnergySurface(self, fstr, Eunits='kJ/mol', drawStructures=False):
        """
        Generates a file containing a rendering of the current potential
        energy surface for this reaction network. The file is saved to at
        location `fstr` on disk. The type of file saved is determined
        automatically from the extension on `fstr`; valid file types are PDF,
        SVG, PS, and PNG. The units to use for energy values can be specified
        using the `Eunits` option; allowed values are ``'J/mol'``, ``'kJ/mol'``,
        ``'cal/mol'``, ``'kcal/mol'``, and ``'cm^-1'``. The drawing is performed
        using the Cairo 2D graphics package; both Cairo and its Python wrapper
        must be installed.
        """

        try:
            import cairo
        except ImportError:
            logging.warning('Cairo not found; potential energy surface will not be drawn.')
            return

        # The order of wells is as follows:
        #   - Reactant channels come first (to the left)
        #   - Isomers are in the middle
        #   - Product channels come last (to the right)
        # This is done because most people will read the PES from left to right
        wells = []
        wells.extend(self.reactants)
        for isomer in self.isomers:
            wells.append([isomer])
        wells.extend(self.products)
        
        # Get minimum and maximum energy on surface (so we can compute bounding box)
        E0min = sum([spec.E0.value_si for spec in wells[0]]); E0max = E0min
        for i, well in enumerate(wells):
            E0 = sum([spec.E0.value_si for spec in well])
            if E0 < E0min: E0min = E0
            if E0 > E0max: E0max = E0
        for rxn in self.pathReactions:
            E0 = rxn.transitionState.E0.value_si
            if E0 < E0min: E0min = E0
            if E0 > E0max: E0max = E0
        
        # Drawing parameters
        padding_left = 96.0
        padding_right = padding_left
        padding_top = padding_left / 2.0
        padding_bottom = padding_left / 2.0
        wellWidth = 64.0; wellSpacing = 64.0; Eslope = 5.0; TSwidth = 16.0
        E0 = [sum([spec.E0.value_si for spec in well]) / 4184 for well in wells]
        E0.extend([rxn.transitionState.E0.value_si / 4184 for rxn in self.pathReactions])
        y_E0 = (max(E0) - 0.0) * Eslope + padding_top

        ext = os.path.splitext(fstr)[1].lower()

        # Draw labels for each well right away
        # We need their size information to ensure that wells don't overlap
        labels = []
        for well in wells:
            surfaces0, width0, height0, boundingRects0 = self.__drawLabel(well, ext)
            labels.append([surfaces0, width0, height0, boundingRects0])

        # Determine naive position of each well (one per column)
        coordinates = numpy.zeros((len(wells), 2), numpy.float64)
        x = padding_left + wellWidth / 2.0
        for i, well in enumerate(wells):
            E0 = sum([spec.E0.value_si for spec in well]) / 4184
            y = y_E0 - E0 * Eslope
            coordinates[i] = [x, y]
            x += wellWidth + wellSpacing

        # Squish columns together from the left where possible until an isomer is encountered
        Nleft = wells.index([self.isomers[0]])
        for i in range(Nleft-1, -1, -1):
            newX = float(coordinates[i,0])
            for j in range(i+1, Nleft):
                spacing = labels[i][2] if coordinates[i,1] < coordinates[j,1] else labels[j][2]
                spacing += 24
                if abs(coordinates[i,1] - coordinates[j,1]) < spacing:
                    newX = float(coordinates[j,0]) - (wellWidth + wellSpacing)
                    break
                else:
                    newX = float(coordinates[j,0])
            coordinates[i,0] = newX
        # Squish columns together from the right where possible until an isomer is encountered
        Nright = wells.index([self.isomers[-1]])
        for i in range(Nright+2, len(wells)):
            newX = float(coordinates[i,0])
            for j in range(i-1, Nright, -1):
                if abs(coordinates[i,1] - coordinates[j,1]) < 72:
                    newX = float(coordinates[j,0]) + (wellWidth + wellSpacing)
                    break
                else:
                    newX = float(coordinates[j,0])
            coordinates[i,0] = newX

        coordinates[:,0] -= numpy.min(coordinates[:,0]) - padding_left - wellWidth/2.0

        # Determine required size of diagram
        width = numpy.max(coordinates[:,0]) - numpy.min(coordinates[:,0]) + wellWidth + padding_left + padding_right
        height = (E0max - E0min) / 4184. * Eslope + 32.0 + padding_top + padding_bottom

        # Choose multiplier to convert energies to desired units
        if Eunits == 'J/mol':      Emult = 1.0
        elif Eunits == 'kJ/mol':   Emult = 1.0 / 1000
        elif Eunits == 'cal/mol':  Emult = 1.0 / 4.184
        elif Eunits == 'kcal/mol': Emult = 1.0 / 4184
        elif Eunits == 'cm^-1':    Emult = 1.0 / 11.96
        else:
            logging.warning('Invalid value "{0}" for Eunits parameter. Setting to "kJ/mol".'.format(Eunits))
            Emult = 1.0 / 1000

        # Initialize Cairo surface and context
        surface, cr = self.__createNewSurfaceAndContext(ext, fstr, width, height)

        # Some global settings
        cr.select_font_face("sans")
        cr.set_font_size(10)

        # If we are drawing a PNG, first fill the background with white
        ext = os.path.splitext(fstr)[1].lower()
        if ext == '.png':
            cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
            cr.rectangle(0, 0, width, height)
            cr.fill()
        
        # Draw path reactions
        for rxn in self.pathReactions:
            reac = wells.index(rxn.reactants)
            prod = wells.index(rxn.products)
            E0_reac = sum([spec.E0.value_si for spec in wells[reac]]) / 4184
            E0_prod = sum([spec.E0.value_si for spec in wells[prod]]) / 4184
            E0_TS = rxn.transitionState.E0.value_si / 4184
            if reac < prod:
                x1, y1 = coordinates[reac,:]
                x2, y2 = coordinates[prod,:]
            else:
                x1, y1 = coordinates[prod,:]
                x2, y2 = coordinates[reac,:]
            x1 += wellSpacing / 2.0; x2 -= wellSpacing / 2.0
            if abs(E0_TS - E0_reac) > 0.1 and abs(E0_TS - E0_prod) > 0.1:
                if len(rxn.reactants) == 2:
                    if reac < prod: x0 = x1 + wellSpacing * 0.5
                    else:           x0 = x2 - wellSpacing * 0.5
                elif len(rxn.products) == 2:
                    if reac < prod: x0 = x2 - wellSpacing * 0.5
                    else:           x0 = x1 + wellSpacing * 0.5
                else:
                    x0 = 0.5 * (x1 + x2)
                y0 = y_E0 - E0_TS * Eslope
                width1 = (x0 - x1)
                width2 = (x2 - x0)
                # Draw horizontal line for TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.set_line_width(2.0)
                cr.move_to(x0 - TSwidth/2.0, y0)
                cr.line_to(x0+TSwidth/2.0, y0)
                cr.stroke()
                # Add background and text for energy
                E0 = "{0:.1f}".format(rxn.transitionState.E0.value_si * Emult)
                extents = cr.text_extents(E0)
                x = x0 - extents[2] / 2.0; y = y0 - 6.0
                cr.rectangle(x + extents[0] - 2.0, y + extents[1] - 2.0, extents[2] + 4.0, extents[3] + 4.0)
                cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
                cr.fill()
                cr.move_to(x, y)
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.show_text(E0)
                # Draw Bezier curve connecting reactants and products through TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
                cr.set_line_width(1.0)
                cr.move_to(x1, y1)
                cr.curve_to(x1 + width1/8.0, y1,   x0 - width1/8.0 - TSwidth/2.0, y0,   x0 - TSwidth/2.0, y0)
                cr.move_to(x0 + TSwidth/2.0, y0)
                cr.curve_to(x0 + width2/8.0 + TSwidth/2.0, y0,   x2 - width2/8.0, y2,   x2, y2)
                cr.stroke()
            else:
                width = (x2 - x1)
                # Draw Bezier curve connecting reactants and products through TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
                cr.set_line_width(1.0)
                cr.move_to(x1, y1)
                cr.curve_to(x1 + width/4.0, y1,   x2 - width/4.0, y2,   x2, y2)
                cr.stroke()

        # Draw wells (after path reactions so that they are on top)
        for i, well in enumerate(wells):
            x0, y0 = coordinates[i,:]
            # Draw horizontal line for well
            cr.set_line_width(4.0)
            cr.move_to(x0 - wellWidth/2.0, y0)
            cr.line_to(x0 + wellWidth/2.0, y0)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.stroke()
            # Add background and text for energy
            E0 = "{0:.1f}".format(sum([spec.E0.value_si for spec in well]) * Emult)
            extents = cr.text_extents(E0)
            x = x0 - extents[2] / 2.0; y = y0 - 6.0
            cr.rectangle(x + extents[0] - 2.0, y + extents[1] - 2.0, extents[2] + 4.0, extents[3] + 4.0)
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            cr.move_to(x, y)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.show_text(E0)
            # Add background for label and label itself
            surfaces0, width0, height0, boundingRects0 = labels[i]
            left = x0 - width0/2.0; top = y0 + 4
            for j, surf in enumerate(surfaces0):
                cr.save()
                cr.rectangle(left + boundingRects0[j][0], top + boundingRects0[j][1], boundingRects0[j][2], boundingRects0[j][3])
                cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
                cr.fill()
                cr.set_source_surface(surf, left + boundingRects0[j][0], top + boundingRects0[j][1])
                cr.paint()
                cr.restore()

        # Finish Cairo drawing
        ext = os.path.splitext(fstr)[1].lower()
        if ext == '.png':
            surface.write_to_png(fstr)
        else:
            surface.finish()

    def __drawText(self, text, ext, padding=0):
        """
        Create and return a temporary Cairo surface containing the string
        `text` with an optional amount of `padding` on all sides. The type of
        surface is dictated by the `ext` parameter.
        """
        from rmgpy.molecule.draw import MoleculeDrawer, createNewSurface
        import cairo
        
        fontSizeNormal = MoleculeDrawer().options['fontSizeNormal']

        surface0 = createNewSurface(type=ext[1:])
        cr0 = cairo.Context(surface0)
        cr0.set_font_size(fontSizeNormal)
        extents = cr0.text_extents(text)
        width = extents[2] + 2 * padding; height = extents[3] + 2 * padding

        surface = createNewSurface(type=ext[1:], width=width, height=height)
        cr = cairo.Context(surface)
        cr.set_font_size(fontSizeNormal)
        cr.move_to(padding - extents[0], padding - extents[1])
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(text)
        
        return surface, cr, [0, 0, width, height]

    def __drawLabel(self, speciesList, ext, useLabels=False):
        """
        For a local minima on the potential energy surface composed of a list
        of species `speciesList`, generate a Cairo surface containing the label
        by which to mark the local minima. If molecular structure data is
        available for all species, structures will be used; otherwise labels
        will be used.
        """

        from rmgpy.molecule.draw import MoleculeDrawer, createNewSurface
        import cairo
        
        # Determine whether or not to use the molecular structures in the
        # label
        # All species must have structure data in order to use the structures
        # Otherwise the labels are used
        for spec in speciesList:
            if spec.molecule is None or len(spec.molecule) == 0:
                useLabels = True
                break

        plusSurface, cr0, plusBoundingRect = self.__drawText('+', ext, padding=2)
        plusWidth = plusBoundingRect[2]; plusHeight = plusBoundingRect[3]
        
        # Render
        surfaces = []; boundingRects = []; width = 0.0; height = 0.0
        # Draw the molecular structures or labels on temporary Cairo surfaces
        for spec in speciesList:
            if useLabels:
                surface0, cr0, boundingRect0 = self.__drawText(spec.label, ext, padding=2)
            else:
                surface0, cr0, boundingRect0 = MoleculeDrawer().draw(spec.molecule[0], format=ext[1:])
            surfaces.append(surface0)
            boundingRects.append(list(boundingRect0))
            if width < boundingRect0[2]: width = boundingRect0[2]
            height += boundingRect0[3]

        # Sort the structures from widest to narrowest
        widths = numpy.array([rect[2] for rect in boundingRects], numpy.float64)
        indices = list(widths.argsort())
        indices.reverse()
        surfaces = [surfaces[i] for i in indices]
        boundingRects = [boundingRects[i] for i in indices]

        # Insert plus between each structure
        for i in range(1, len(speciesList)):
            surfaces.insert(len(speciesList)-i, plusSurface)
            boundingRects.insert(len(speciesList)-i, [0,0,plusWidth,plusHeight])
            height += plusHeight

        # Set left and top of each bounding rectangle
        left = 0; top = 0
        for i in range(len(surfaces)):
            boundingRects[i][0] = (width - boundingRects[i][2]) / 2.0
            boundingRects[i][1] = top
            top += boundingRects[i][3]

        return surfaces, width, height, boundingRects
        
