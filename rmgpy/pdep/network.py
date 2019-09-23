#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This module contains the :class:`Network` class, a representation of a 
pressure-dependent unimolecular reaction network
"""
from __future__ import division

import logging
import math

import numpy as np

import rmgpy.constants as constants
from rmgpy.exceptions import NetworkError, InvalidMicrocanonicalRateError
from rmgpy.reaction import Reaction


################################################################################

class Network(object):
    """
    A representation of a unimolecular reaction network. The attributes are:

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `isomers`               A list of the unimolecular isomers in the network
    `reactants`             A list of the bimolecular reactant channels (Configuration objects) in the network
    `products`              A list of the bimolecular product channels (Configuration objects) in the network
    `pathReactions`         A list of "path" reaction objects that connect adjacent isomers (the high-pressure-limit)
    `bathGas`               A dictionary of the bath gas species (keys) and their mole fractions (values)
    `netReactions`          A list of "net" reaction objects that connect any pair of isomers
    ----------------------- ----------------------------------------------------
    `T`                     The current temperature in K
    `P`                     The current pressure in bar
    `Elist`                 The current array of energy grains in kJ/mol
    `Jlist`                 The current array of total angular momentum quantum numbers
    ----------------------- ----------------------------------------------------
    `Nisom`                 The number of unimolecular isomers in the network
    `Nreac`                 The number of bimolecular reactant channels in the network
    `Nprod`                 The number of bimolecular product channels in the network
    `Ngrains`               The number of energy grains
    `NJ`                    The number of angular momentum grains
    ----------------------- ----------------------------------------------------
     `grainSize`            Maximum size of separation between energies
     `grainCount`           Minimum number of descrete energies separated
     `E0`                   A list of ground state energies of isomers, reactants, and products
    `activeKRotor`          ``True`` if the K-rotor is treated as active, ``False`` if treated as adiabatic
    `activeJRotor`          ``True`` if the J-rotor is treated as active, ``False`` if treated as adiabatic
    `rmgmode`               ``True`` if in RMG mode, ``False`` otherwise
    ----------------------- ----------------------------------------------------
    `eqRatios`              An array containing concentration of each isomer and reactant channel present at equilibrium
    `collFreq`              An array of the frequency of collision between
    `Mcoll`                 Matrix of first-order rate coefficients for collisional population transfer between grains for each isomer
    `densStates`            3D np array of stable configurations, number of grains, and number of J
    ======================= ====================================================
    
    """

    def __init__(self, label='', isomers=None, reactants=None, products=None,
                 pathReactions=None, bathGas=None, netReactions=None, T=0.0, P=0.0,
                 Elist=None, Jlist=None, Ngrains=0, NJ=0, activeKRotor=True,
                 activeJRotor=True, grainSize=0.0, grainCount=0, E0=None):
        """
        To initialize a Network object for running a pressure dependent job,
        only label, isomers, reactants, products pathReactions and bathGas are useful,
        since the other attributes will be created during the run.

        The other attributes are used to reinstantiate the created network object
        for debugging and testing.
        """
        self.label = label
        self.isomers = isomers or []
        self.reactants = reactants or []
        self.products = products or []
        self.pathReactions = pathReactions or []
        self.bathGas = bathGas or {}
        self.netReactions = netReactions or []

        self.T = T
        self.P = P
        self.Elist = Elist
        self.Jlist = Jlist

        self.Nisom = len(self.isomers)
        self.Nreac = len(self.reactants)
        self.Nprod = len(self.products)
        self.Ngrains = Ngrains
        self.NJ = NJ

        self.activeKRotor = activeKRotor
        self.activeJRotor = activeJRotor

        self.grainSize = grainSize
        self.grainCount = grainCount
        self.E0 = E0

        self.valid = False

    def __repr__(self):
        string = 'Network('
        if self.label != '': string += 'label="{0}", '.format(self.label)
        if self.isomers: string += 'isomers="{0!r}", '.format(self.isomers)
        if self.reactants: string += 'reactants="{0!r}", '.format(self.reactants)
        if self.products: string += 'products="{0!r}", '.format(self.products)
        if self.pathReactions: string += 'pathReactions="{0!r}", '.format(self.pathReactions)
        if self.bathGas: string += 'bathGas="{0!r}", '.format(self.bathGas)
        if self.netReactions: string += 'netReactions="{0!r}", '.format(self.netReactions)
        if self.T != 0.0: string += 'T="{0}", '.format(self.T)
        if self.P != 0.0: string += 'P="{0}", '.format(self.P)
        if self.Elist is not None: string += 'Elist="{0}", '.format(self.Elist)
        if self.Jlist is not None: string += 'Jlist="{0}", '.format(self.Jlist)
        if self.Ngrains != 0: string += 'Ngrains="{0}", '.format(self.Ngrains)
        if self.NJ != 0: string += 'NJ="{0}", '.format(self.NJ)
        string += 'activeKRotor="{0}", '.format(self.activeKRotor)
        string += 'activeJRotor="{0}", '.format(self.activeJRotor)
        if self.grainSize != 0.0: string += 'grainSize="{0}", '.format(self.grainSize)
        if self.grainCount != 0: string += 'grainCount="{0}", '.format(self.grainCount)
        if self.E0 is not None: string += 'E0="{0}", '.format(self.E0)
        string += ')'
        return string

    def __str__(self):
        """return Network like it would be seen in an Arkane input file"""
        string = "Network(label = '{0}', isomers = {1}, reactants = {2}, products = {3}, " \
                 "pathReactions = {4}, bathGas = {5}, netReactions = {6})".format(
                        self.label, [i.species[0].label for i in self.isomers],
                        [[r.label for r in pair.species] for pair in self.reactants],
                        [[p.label for p in pair.species] for pair in self.products],
                        [r.label for r in self.pathReactions],
                        dict([(s.label, value) for s, value in self.bathGas.items()]),
                        [r.toLabeledStr() for r in self.netReactions]
                    )
        return string

    def invalidate(self):
        """
        Mark the network as in need of a new calculation to determine the
        pressure-dependent rate coefficients
        """
        self.valid = False

    def getAllSpecies(self):
        """
        Return a list of all unique species in the network, including all
        isomers, reactant and product channels, and bath gas species.
        """
        species_list = []
        for isomer in self.isomers:
            for spec in isomer.species:
                if spec not in species_list:
                    species_list.append(spec)
        for reactant in self.reactants:
            for spec in reactant.species:
                if spec not in species_list:
                    species_list.append(spec)
        for product in self.products:
            for spec in product.species:
                if spec not in species_list:
                    species_list.append(spec)
        for spec in self.bathGas:
            if spec not in species_list:
                species_list.append(spec)
        return species_list

    def initialize(self, Tmin, Tmax, Pmin, Pmax, maximumGrainSize=0.0, minimumGrainCount=0, activeJRotor=True,
                   activeKRotor=True, rmgmode=False):
        """
        Initialize a pressure dependence calculation by computing several
        quantities that are independent of the conditions. You must specify
        the temperature and pressure ranges of interesting using `Tmin` and
        `Tmax` in K and `Pmin` and `Pmax` in Pa. You must also specify the
        maximum energy grain size `grainSize` in J/mol and/or the minimum
        number of grains `grainCount`.
        """

        logging.debug("initializing network")
        if maximumGrainSize == 0.0 and minimumGrainCount == 0:
            raise NetworkError('Must provide either grainSize or Ngrains parameter to Network.determineEnergyGrains().')

        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.grainSize = maximumGrainSize
        self.grainCount = minimumGrainCount

        self.Nisom = len(self.isomers)
        self.Nreac = len(self.reactants)
        self.Nprod = len(self.products)
        self.Ngrains = 0
        self.NJ = 0

        # Calculate ground-state energies
        self.E0 = np.zeros((self.Nisom + self.Nreac + self.Nprod), np.float64)
        for i in range(self.Nisom):
            self.E0[i] = self.isomers[i].E0
        for n in range(self.Nreac):
            self.E0[n + self.Nisom] = self.reactants[n].E0
        for n in range(self.Nprod):
            self.E0[n + self.Nisom + self.Nreac] = self.products[n].E0

        # Calculate densities of states
        self.activeJRotor = activeJRotor
        self.activeKRotor = activeKRotor
        self.rmgmode = rmgmode

        self.calculateDensitiesOfStates()
        logging.debug('Finished initialization for network {0}.'.format(self.label))
        logging.debug('The network now has values of {0}'.format(repr(self)))

    def calculateRateCoefficients(self, Tlist, Plist, method, errorCheck=True):

        n_isom = len(self.isomers)
        n_reac = len(self.reactants)
        n_prod = len(self.products)

        for rxn in self.pathReactions:
            if len(rxn.transition_state.conformer.modes) > 0:
                logging.debug('Using RRKM theory to compute k(E) for path reaction {0}.'.format(rxn))
            elif rxn.kinetics is not None:
                logging.debug('Using ILT method to compute k(E) for path reaction {0}.'.format(rxn))
        logging.debug('')

        logging.info('Calculating phenomenological rate coefficients for {0}...'.format(rxn))
        K = np.zeros((len(Tlist), len(Plist), n_isom + n_reac + n_prod, n_isom + n_reac + n_prod), np.float64)

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
                    raise NetworkError('Unknown method "{0}". Valid options are "modified strong collision", '
                                       '"reservoir state", or "chemically-significant eigenvalues"'.format(method))

                K[t, p, :, :] = self.K

                # Check that the k(T,P) values satisfy macroscopic equilibrium
                eq_ratios = self.eqRatios
                for i in range(n_isom + n_reac):
                    for j in range(i):
                        Keq0 = K[t, p, j, i] / K[t, p, i, j]
                        Keq = eq_ratios[j] / eq_ratios[i]
                        if Keq0 / Keq < 0.5 or Keq0 / Keq > 2.0:
                            if i < n_isom:
                                reactants = self.isomers[i]
                            elif i < n_isom + n_reac:
                                reactants = self.reactants[i - n_isom]
                            else:
                                reactants = self.products[i - n_isom - n_reac]
                            if j < n_isom:
                                products = self.isomers[j]
                            elif j < n_isom + n_reac:
                                products = self.reactants[j - n_isom]
                            else:
                                products = self.products[j - n_isom - n_reac]
                            reaction = Reaction(reactants=reactants.species[:], products=products.species[:])
                            logging.error('For net reaction {0!s}:'.format(reaction))
                            logging.error('Expected Keq({1:g} K, {2:g} bar) = {0:11.3e}'.format(Keq, T, P * 1e-5))
                            logging.error('  Actual Keq({1:g} K, {2:g} bar) = {0:11.3e}'.format(Keq0, T, P * 1e-5))
                            raise NetworkError('Computed k(T,P) values for reaction {0!s} do not satisfy macroscopic '
                                               'equilibrium.'.format(reaction))

                # Reject if any rate coefficients are negative
                if errorCheck:
                    negative_rate = False
                    for i in range(n_isom + n_reac + n_prod):
                        for j in range(i):
                            if (K[t, p, i, j] < 0 or K[t, p, j, i] < 0) and not negative_rate:
                                negative_rate = True
                                logging.error('Negative rate coefficient generated; rejecting result.')
                                logging.info(K[t, p, 0:n_isom + n_reac + n_prod, 0:n_isom + n_reac])
                                K[t, p, :, :] = 0 * K[t, p, :, :]
                                self.K = 0 * self.K
        logging.debug('Finished calculating rate coefficients for network {0}.'.format(self.label))
        logging.debug('The network now has values of {0}'.format(repr(self)))
        logging.debug('Master equation matrix found for network {0} is {1}'.format(self.label, K))
        return K

    def setConditions(self, T, P, ymB=None):
        """
        Set the current network conditions to the temperature `T` in K and
        pressure `P` in Pa. All of the internal variables are updated 
        accordingly if they are out of date. For example, those variables that
        depend only on temperature will not be recomputed if the temperature
        is the same.
        """

        temperature_changed = (self.T != T)
        pressure_changed = (self.P != P)
        self.T = T
        self.P = P
        self.ymB = ymB

        n_isom = self.Nisom
        n_reac = self.Nreac
        n_proc = self.Nprod

        e0 = self.E0
        grain_size = self.grainSize
        grain_count = self.grainCount

        success = False
        previous_error = None
        while not success:
            success = True  # (set it to false again later if necessary)
            # Update parameters that depend on temperature only if necessary
            if temperature_changed:

                # Choose the energy grains to use to compute k(T,P) values at this temperature
                e_list = self.Elist = self.selectEnergyGrains(T, grain_size, grain_count)
                n_grains = self.Ngrains = len(e_list)
                logging.info('Using {0:d} grains from {1:.2f} to {2:.2f} kJ/mol in steps of {3:.2f} kJ/mol to compute '
                             'the k(T,P) values at {4:g} K'.format(n_grains, np.min(e_list) * 0.001,
                                                                   np.max(e_list) * 0.001,
                                                                   (e_list[1] - e_list[0]) * 0.001, T))

                # Choose the angular momenta to use to compute k(T,P) values at this temperature
                # (This only applies if the J-rotor is adiabatic
                if not self.activeJRotor:
                    j_list = self.Jlist = np.arange(0, 20, 1, np.int)
                    n_j = self.NJ = len(j_list)
                else:
                    j_list = self.Jlist = np.array([0], np.int)
                    n_j = self.NJ = 1

                # Map the densities of states onto this set of energies
                # Also shift each density of states to a common zero of energy
                self.mapDensitiesOfStates()

                # Use free energy to determine equilibrium ratios of each isomer and product channel
                self.calculateEquilibriumRatios()

                # Calculate microcanonical rate coefficients for each path reaction
                try:
                    self.calculateMicrocanonicalRates()
                except InvalidMicrocanonicalRateError as error:
                    badness = error.badness()
                    if previous_error and (previous_error.message == error.message):
                        # only compare badness if same reaction is causing problem
                        improvement = previous_error.badness() / badness
                        if improvement < 0.2 or (n_grains > 1e4 and improvement < 1.1) or (n_grains > 1.5e6):
                            # allow it to get worse at first
                            logging.error(error.message)
                            logging.error("Increasing number of grains did not decrease error enough "
                                          "(Current badness: {0:.1f}, previous {1:.1f}). Something must be wrong with "
                                          "network {2}".format(badness, previous_error.badness(), self.label))
                            raise error
                    previous_error = error
                    success = False
                    grain_size *= 0.5
                    grain_count *= 2
                    logging.warning("Increasing number of grains, decreasing grain size and trying again. "
                                    "(Current badness: {0:.1f})".format(badness))
                    continue
                else:
                    success = True

                # Rescale densities of states such that, when they are integrated
                # using the Boltzmann factor as a weighting factor, the result is unity
                for i in range(n_isom + n_reac):
                    Q = 0.0
                    for s in range(n_j):
                        Q += np.sum(
                            self.densStates[i, :, s] * (2 * j_list[s] + 1) * np.exp(-e_list / constants.R / T))
                    self.densStates[i, :, :] /= Q

            # Update parameters that depend on temperature and pressure if necessary
            if temperature_changed or pressure_changed:
                self.calculateCollisionModel()
        logging.debug('Finished setting conditions for network {0}.'.format(self.label))
        logging.debug('The network now has values of {0}'.format(repr(self)))

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
            use_grain_size = True
        elif grainCount > 0 and grainSize <= 0.0:
            # Only number of grains was specified, so we must use it
            use_grain_size = False
        else:
            # Both were specified, so we choose the tighter constraint
            # (i.e. the one that will give more grains, and so better accuracy)
            grain_size0 = (Emax - Emin) / (grainCount - 1)
            use_grain_size = (grain_size0 > grainSize)

        # Generate the array of energies
        if use_grain_size:
            e_list = np.arange(Emin, Emax + grainSize, grainSize, dtype=np.float64)
        else:
            e_list = np.linspace(Emin, Emax, grainCount, dtype=np.float64)

        return e_list

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
        e_min = np.min(self.E0)
        e_min = math.floor(e_min)  # Round to nearest whole number

        # Use the highest energy on the PES as the initial guess for Emax0
        e_max = np.max(self.E0)
        for rxn in self.pathReactions:
            E0 = float(rxn.transition_state.conformer.E0.value_si)
            if E0 > e_max: e_max = E0

        # Choose the actual e_max as many kB * T above the maximum energy on the PES
        # You should check that this is high enough so that the Boltzmann distributions have trailed off to negligible values
        e_max += 40. * constants.R * T

        return self.__getEnergyGrains(e_min, e_max, grainSize, grainCount)

    def calculateDensitiesOfStates(self):
        """
        Calculate the densities of states of each configuration that has states
        data. The densities of states are computed such that they can be
        applied to each temperature in the range of interest by interpolation.
        """

        Tmin = self.Tmin
        Tmax = self.Tmax
        grain_size = self.grainSize
        grain_count = self.grainCount

        n_isom = self.Nisom
        n_reac = self.Nreac
        n_prod = self.Nprod

        logging.info('Calculating densities of states for {0} network...'.format(self.label))

        # Choose the energies used to compute the densities of states
        # Use Tmin to select the minimum energy and grain size
        e_list0 = self.selectEnergyGrains(Tmin, grain_size, grain_count)
        e_min0 = np.min(e_list0)
        grain_size0 = e_list0[1] - e_list0[0]
        # Use Tmax to select the maximum energy and grain count
        e_list0 = self.selectEnergyGrains(Tmax, grain_size, grain_count)
        grain_count0 = len(e_list0)
        e_max0 = np.max(e_list0)

        e_list = self.__getEnergyGrains(e_min0, e_max0, grain_size0, grain_count0)
        n_grains = len(e_list)
        de = e_list[1] - e_list[0]
        logging.info('Using {0:d} grains from {1:.2f} to {2:.2f} kJ/mol in steps of {3:.2f} kJ/mol to compute '
                     'densities of states'.format(n_grains, e_list[0] * 0.001, e_list[-1] * 0.001, de * 0.001))

        # Shift the energy grains so that the minimum grain is zero
        e_list -= e_list[0]

        dens_states = np.zeros((n_isom + n_reac + n_prod, n_grains), np.float64)

        # Densities of states for isomers
        for i in range(n_isom):
            logging.debug('Calculating density of states for isomer "{0}"'.format(self.isomers[i]))
            self.isomers[i].calculateDensityOfStates(e_list, activeKRotor=self.activeKRotor,
                                                     activeJRotor=self.activeJRotor, rmgmode=self.rmgmode)

        # Densities of states for reactant channels
        for n in range(n_reac):
            if self.reactants[n].hasStatMech():
                logging.debug('Calculating density of states for reactant channel "{0}"'.format(self.reactants[n]))
                self.reactants[n].calculateDensityOfStates(e_list, activeKRotor=self.activeKRotor,
                                                           activeJRotor=self.activeJRotor, rmgmode=self.rmgmode)
            else:
                logging.warning(
                    'NOT calculating density of states for reactant channel "{0}". Missing Statmech.'.format(
                        self.reactants[n]))
                logging.warning('Reactants: {}'.format(repr(self.reactants[n])))
        # Densities of states for product channels
        if not self.rmgmode:
            for n in range(n_prod):
                if self.products[n].hasStatMech():
                    logging.debug('Calculating density of states for product channel "{0}"'.format(self.products[n]))
                    self.products[n].calculateDensityOfStates(e_list, activeKRotor=self.activeKRotor,
                                                              activeJRotor=self.activeJRotor, rmgmode=self.rmgmode)
                else:
                    logging.warning(
                        'NOT calculating density of states for product channel "{0}" Missing Statmech.'.format(
                            self.products[n]))
                    logging.warning('Products: {}'.format(repr(self.products[n])))
        logging.debug('')

        # import pylab
        # for i in range(n_isom):
        #    pylab.semilogy(e_list*0.001, self.isomers[i].densStates)
        # for n in range(n_reac):
        #    if self.reactants[n].densStates is not None:
        #        pylab.semilogy(e_list*0.001, self.reactants[n].densStates)
        # for n in range(n_prod):
        #    if self.products[n].densStates is not None:
        #        pylab.semilogy(e_list*0.001, self.products[n].densStates)
        # pylab.show()

    def mapDensitiesOfStates(self):
        """
        Map the overall densities of states to the current energy grains.
        Semi-logarithmic interpolation will be used if the grain sizes of 
        `Elist0` and `Elist` do not match; this should not be a significant
        source of error as long as the grain sizes are sufficiently small.
        """
        n_isom = len(self.isomers)
        n_reac = len(self.reactants)
        n_prod = len(self.products)
        n_grains = len(self.Elist)
        n_j = len(self.Jlist)

        self.densStates = np.zeros((n_isom + n_reac + n_prod, n_grains, n_j))
        # Densities of states for isomers
        for i in range(n_isom):
            logging.debug('Mapping density of states for isomer "{0}"'.format(self.isomers[i]))
            self.densStates[i, :, :] = self.isomers[i].mapDensityOfStates(self.Elist, self.Jlist)
        # Densities of states for reactant channels
        for n in range(n_reac):
            if self.reactants[n].densStates is not None:
                logging.debug('Mapping density of states for reactant channel "{0}"'.format(self.reactants[n]))
                self.densStates[n + n_isom, :, :] = self.reactants[n].mapDensityOfStates(self.Elist, self.Jlist)
        # Densities of states for product channels
        for n in range(n_prod):
            if self.products[n].densStates is not None:
                logging.debug('Mapping density of states for product channel "{0}"'.format(self.products[n]))
                self.densStates[n + n_isom + n_reac, :, :] = self.products[n].mapDensityOfStates(self.Elist, self.Jlist)

        # import pylab
        # for i in range(n_isom + n_reac + n_prod):
        #    pylab.semilogy(self.Elist*0.001, self.densStates[i,:])
        # pylab.show()

    def calculateMicrocanonicalRates(self):
        """
        Calculate and return arrays containing the microcanonical rate
        coefficients :math:`k(E)` for the isomerization, dissociation, and
        association path reactions in the network.
        """

        temperature = self.T
        e_list = self.Elist
        j_list = self.Jlist
        dens_states = self.densStates
        n_grains = len(e_list)
        n_isom = len(self.isomers)
        n_reac = len(self.reactants)
        n_prod = len(self.products)
        n_j = 1 if self.activeJRotor else len(j_list)

        self.Kij = np.zeros([n_isom, n_isom, n_grains, n_j], np.float64)
        self.Gnj = np.zeros([n_reac + n_prod, n_isom, n_grains, n_j], np.float64)
        self.Fim = np.zeros([n_isom, n_reac, n_grains, n_j], np.float64)

        isomers = [isomer.species[0] for isomer in self.isomers]
        reactants = [reactant.species for reactant in self.reactants]
        products = [product.species for product in self.products]

        for rxn in self.pathReactions:
            if rxn.reactants[0] in isomers and rxn.products[0] in isomers:
                # Isomerization
                reac = isomers.index(rxn.reactants[0])
                prod = isomers.index(rxn.products[0])
            elif rxn.reactants[0] in isomers and rxn.products in reactants:
                # Dissociation (reversible)
                reac = isomers.index(rxn.reactants[0])
                prod = reactants.index(rxn.products) + n_isom
            elif rxn.reactants[0] in isomers and rxn.products in products:
                # Dissociation (irreversible)
                reac = isomers.index(rxn.reactants[0])
                prod = products.index(rxn.products) + n_isom + n_reac
            elif rxn.reactants in reactants and rxn.products[0] in isomers:
                # Association (reversible)
                reac = reactants.index(rxn.reactants) + n_isom
                prod = isomers.index(rxn.products[0])
            elif rxn.reactants in products and rxn.products[0] in isomers:
                # Association (irreversible)
                reac = products.index(rxn.reactants) + n_isom + n_reac
                prod = isomers.index(rxn.products[0])
            else:
                logging.info('\nUnexpected type of path reaction.')
                logging.info('\nnetwork reactants:')
                for index, reacts in enumerate(reactants, 1):
                    logging.info('reactants {0}:'.format(index))
                    for subindex, react in enumerate(reacts, 1):
                        logging.info('reactant {0}:'.format(subindex))
                        for mol in react.molecule:
                            logging.info(mol.toAdjacencyList())
                            logging.info('reactive = {0}\n'.format(mol.reactive))
                logging.info('network products:')
                for index, prods in enumerate(products, 1):
                    logging.info('products {0}:'.format(index))
                    for subindex, pro in enumerate(prods, 1):
                        for mol in pro.molecule:
                            logging.info(mol.toAdjacencyList())
                            logging.info('reactive = {0}\n'.format(mol.reactive))
                logging.info('network isomers:')
                for index, iso in enumerate(isomers, 1):
                    logging.info('isomer {0}:'.format(index))
                    for mol in iso.molecule:
                        logging.info(mol.toAdjacencyList())
                        logging.info('reactive = {0}\n'.format(mol.reactive))
                logging.info('rxn reactants:')
                for index, react in enumerate(rxn.reactants, 1):
                    logging.info('reactant {0}:'.format(index))
                    for mol in react.molecule:
                        logging.info(mol.toAdjacencyList())
                        logging.info('reactive = {0}\n'.format(mol.reactive))
                logging.info('rxn products:')
                for index, pro in enumerate(rxn.products, 1):
                    logging.info('product {0}:'.format(index))
                    for mol in pro.molecule:
                        logging.info(mol.toAdjacencyList())
                        logging.info('reactive = {0}\n'.format(mol.reactive))
                logging.info('Path reaction {0} not found in reaction network {1}'.format(rxn, self.label))
                continue

            # Compute the microcanonical rate coefficient k(E)
            reac_dens_states = dens_states[reac, :, :]
            prod_dens_states = dens_states[prod, :, :]
            kf, kr = rxn.calculateMicrocanonicalRateCoefficient(self.Elist, self.Jlist,
                                                                reac_dens_states, prod_dens_states,
                                                                temperature)

            # Check for NaN (just to be safe)
            if np.isnan(kf).any() or np.isnan(kr).any():
                raise NetworkError('One or more k(E) values is NaN for path reaction "{0}".'.format(rxn))

            # Determine the expected value of the rate coefficient k(T)
            if rxn.canTST():
                # RRKM theory was used to compute k(E), so use TST to compute k(T)
                logging.debug('Using RRKM rate for Expected kf')
                kf_expected = rxn.calculateTSTRateCoefficient(temperature)
            else:
                # ILT was used to compute k(E), so use high-P kinetics to compute k(T)
                logging.debug('Using high pressure rate coefficient rate for Expected kf')
                kf_expected = rxn.kinetics.getRateCoefficient(temperature) if rxn.network_kinetics is None else \
                    rxn.network_kinetics.getRateCoefficient(temperature)

            # Determine the expected value of the equilibrium constant (Kc)
            Keq_expected = self.eqRatios[prod] / self.eqRatios[reac]

            # Determine the actual values of k(T) and Keq
            C0 = 1e5 / (constants.R * temperature)
            kf0 = 0.0
            kr0 = 0.0
            Qreac = 0.0
            Qprod = 0.0
            for s in range(n_j):
                kf0 += np.sum(kf[:, s] * reac_dens_states[:, s] * (2 * j_list[s] + 1)
                                 * np.exp(-e_list / constants.R / temperature))
                kr0 += np.sum(kr[:, s] * prod_dens_states[:, s] * (2 * j_list[s] + 1)
                                 * np.exp(-e_list / constants.R / temperature))
                Qreac += np.sum(reac_dens_states[:, s] * (2 * j_list[s] + 1)
                                   * np.exp(-e_list / constants.R / temperature))
                Qprod += np.sum(prod_dens_states[:, s] * (2 * j_list[s] + 1)
                                   * np.exp(-e_list / constants.R / temperature))
            kr0 *= C0 ** (len(rxn.products) - len(rxn.reactants))
            Qprod *= C0 ** (len(rxn.products) - len(rxn.reactants))
            kf_actual = kf0 / Qreac if Qreac > 0 else 0
            kr_actual = kr0 / Qprod if Qprod > 0 else 0
            Keq_actual = kf_actual / kr_actual if kr_actual > 0 else 0

            error = False
            warning = False
            k_ratio = 1.0
            Keq_ratio = 1.0
            # Check that the forward rate coefficient is correct
            if kf_actual > 0:
                k_ratio = kf_expected / kf_actual
                # Rescale kf and kr so that we get kf_expected
                kf *= k_ratio
                kr *= k_ratio
                # Decide if the disagreement warrants a warning or error
                if 0.8 < k_ratio < 1.25:
                    # The difference is probably just due to numerical error
                    pass
                elif 0.5 < k_ratio < 2.0:
                    # Might be numerical error, but is pretty large, so warn
                    warning = True
                else:
                    # Disagreement is too large, so raise exception
                    error = True

            # Check that the equilibrium constant is correct
            if Keq_actual > 0:
                Keq_ratio = Keq_expected / Keq_actual
                # Rescale kr so that we get Keq_expected
                kr /= Keq_ratio
                # In RMG jobs this never represents an error because we are
                # missing or using approximate degrees of freedom anyway
                if self.rmgmode:
                    pass
                # Decide if the disagreement warrants a warning or error
                elif 0.8 < Keq_ratio < 1.25:
                    # The difference is probably just due to numerical error
                    pass
                elif 0.5 < Keq_ratio < 2.0:
                    # Might be numerical error, but is pretty large, so warn
                    warning = True
                else:
                    # Disagreement is too large, so raise exception
                    error = True

            if rxn.reactants[0] in isomers and rxn.products[0] in isomers:
                # Isomerization
                self.Kij[prod, reac, :, :] = kf
                self.Kij[reac, prod, :, :] = kr
            elif rxn.reactants[0] in isomers and rxn.products in reactants:
                # Dissociation (reversible)
                self.Gnj[prod - n_isom, reac, :, :] = kf
                self.Fim[reac, prod - n_isom, :, :] = kr
            elif rxn.reactants[0] in isomers and rxn.products in products:
                # Dissociation (irreversible)
                self.Gnj[prod - n_isom, reac, :, :] = kf
            elif rxn.reactants in reactants and rxn.products[0] in isomers:
                # Association (reversible)
                self.Fim[prod, reac - n_isom, :, :] = kf
                self.Gnj[reac - n_isom, prod, :, :] = kr
            elif rxn.reactants in products and rxn.products[0] in isomers:
                # Association (irreversible)
                self.Gnj[reac - n_isom, prod, :, :] = kr
            else:
                raise NetworkError('Unexpected type of path reaction "{0}"'.format(rxn))

            # If the k(E) values are invalid (in that they give the wrong 
            # kf(T) or kr(T) when integrated), then raise an exception
            if error or warning:
                logging.warning('For path reaction {0!s}:'.format(rxn))
                logging.warning('    Expected kf({0:g} K) = {1:g}'.format(temperature, kf_expected))
                logging.warning('      Actual kf({0:g} K) = {1:g}'.format(temperature, kf_actual))
                logging.warning('    Expected Keq({0:g} K) = {1:g}'.format(temperature, Keq_expected))
                logging.warning('      Actual Keq({0:g} K) = {1:g}'.format(temperature, Keq_actual))
                if error:
                    raise InvalidMicrocanonicalRateError(
                        'Invalid k(E) values computed for path reaction "{0}".'.format(rxn), k_ratio, Keq_ratio)
                else:
                    logging.warning('Significant corrections to k(E) to be consistent with high-pressure limit for '
                                    'path reaction "{0}".'.format(rxn))

        # import pylab
        # for prod in range(n_isom):
        #    for reac in range(prod):
        #        pylab.semilogy(self.e_list*0.001, self.Kij[prod,reac,:])
        # for prod in range(n_reac+n_prod):
        #    for reac in range(n_isom):
        #        pylab.semilogy(self.e_list*0.001, self.Gnj[prod,reac,:])
        # pylab.show()

        return self.Kij, self.Gnj, self.Fim

    def calculateEquilibriumRatios(self):
        """
        Return an array containing the fraction of each isomer and reactant
        channel present at equilibrium, as determined from the Gibbs free 
        energy and using the concentration equilibrium constant 
        :math:`K_\\mathrm{c}`. These values are ratios, and the absolute
        magnitude is not guaranteed; however, the implementation scales the 
        elements of the array so that they sum to unity.
        """
        temperature = self.T
        n_isom = len(self.isomers)
        n_reac = len(self.reactants)
        n_prod = len(self.products)
        eq_ratios = np.zeros(n_isom + n_reac + n_prod, np.float64)
        conc = (1e5 / constants.R / temperature)  # [=] mol/m^3
        for i in range(n_isom):
            G = self.isomers[i].getFreeEnergy(temperature)
            logging.debug("Free energy for isomer {} is {}.".format(i, G))
            eq_ratios[i] = math.exp(-G / constants.R / temperature)
        for i in range(n_reac):
            G = self.reactants[i].getFreeEnergy(temperature)
            eq_ratios[n_isom + i] = math.exp(-G / constants.R / temperature) * conc ** (len(self.reactants[i].species) - 1)
        for i in range(n_prod):
            if self.products[i].hasStatMech() or self.products[i].hasThermo():
                G = self.products[i].getFreeEnergy(temperature)
                eq_ratios[n_isom + n_reac + i] = math.exp(-G / constants.R / temperature) * conc ** (len(self.products[i].species) - 1)
        self.eqRatios = eq_ratios
        return eq_ratios / np.sum(eq_ratios)

    def calculateCollisionModel(self):
        """
        Calculate the matrix of first-order rate coefficients for collisional
        population transfer between grains for each isomer, including the
        corresponding collision frequencies.
        """
        n_isom = len(self.isomers)
        n_grains = len(self.Elist)
        n_j = 1 if self.Jlist is None else len(self.Jlist)

        try:
            coll_freq = np.zeros(n_isom, np.float64)
            m_coll = np.zeros((n_isom, n_grains, n_j, n_grains, n_j), np.float64)
        except MemoryError:
            logging.warning('Collision matrix too large to manage')
            new_n_grains = int(n_grains / 2.0)
            logging.warning('Adjusting to use {0} grains instead of {1}'.format(new_n_grains, n_grains))
            self.Elist = self.selectEnergyGrains(self.T, grainCount=new_n_grains)
            return self.calculateCollisionModel()

        for i, isomer in enumerate(self.isomers):
            coll_freq[i] = isomer.calculateCollisionFrequency(self.T, self.P, self.bathGas)
            m_coll[i, :, :, :, :] = coll_freq[i] * isomer.generateCollisionMatrix(self.T, self.densStates[i, :, :],
                                                                                  self.Elist, self.Jlist)

        self.collFreq = coll_freq
        self.Mcoll = m_coll

        return m_coll

    def applyModifiedStrongCollisionMethod(self, efficiencyModel='default'):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the modified strong collision method.
        """
        import rmgpy.pdep.msc as msc
        logging.debug('Applying modified strong collision method at {0:g} K, {1:g} Pa...'.format(self.T, self.P))
        self.K, self.p0 = msc.applyModifiedStrongCollisionMethod(self, efficiencyModel)
        return self.K, self.p0

    def applyReservoirStateMethod(self):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the reservoir state method.
        """
        import rmgpy.pdep.rs as rs
        logging.debug('Applying reservoir state method at {0:g} K, {1:g} Pa...'.format(self.T, self.P))
        self.K, self.p0 = rs.applyReservoirStateMethod(self)
        return self.K, self.p0

    def applyChemicallySignificantEigenvaluesMethod(self, lumpingOrder=None):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the chemically-significant eigenvalues method.
        If a `lumpingOrder` is provided, the algorithm will attempt to lump the
        configurations (given by index) in the order provided, and return a
        reduced set of :math:`k(T,P)` values. 
        """
        import rmgpy.pdep.cse as cse
        logging.debug(
            'Applying chemically-significant eigenvalues method at {0:g} K, {1:g} Pa...'.format(self.T, self.P))
        self.K, self.p0 = cse.applyChemicallySignificantEigenvaluesMethod(self, lumpingOrder)
        return self.K, self.p0

    def generateFullMEMatrix(self, products=True):
        import rmgpy.pdep.me as me
        return me.generateFullMEMatrix(self, products=products)

    def solveFullME(self, tlist, x0):
        """
        Directly solve the full master equation using a stiff ODE solver. Pass the
        reaction `network` to solve, the temperature `T` in K and pressure `P` in
        Pa to solve at, the energies `Elist` in J/mol to use, the output time
        points `tlist` in s, the initial total populations `x0`, the full master
        equation matrix `M`, the accounting matrix `indices` relating isomer and
        energy grain indices to indices of the master equation matrix, and the
        densities of states `densStates` in mol/J of each isomer.
        Returns the times in s, population distributions for each isomer, and total
        population profiles for each configuration.
        """
        import scipy.integrate

        e_list = self.Elist
        j_list = self.Jlist
        dens_states = self.densStates

        n_isom = self.Nisom
        n_reac = self.Nreac
        n_prod = self.Nprod
        n_grains = len(e_list)
        n_j = len(j_list)
        n_time = len(tlist)

        def residual(t, y, K):
            return np.dot(K, y)

        def jacobian(t, y, K):
            return K

        ymB = self.P / constants.R / self.T
        M, indices = self.generateFullMEMatrix()
        n_rows = M.shape[0]
        M[:, n_rows - n_reac - n_prod:] *= ymB

        if self.ymB is not None:
            if isinstance(self.ymB, float):
                assert n_reac <= 1
                M[:, n_rows - n_reac - n_prod:] *= self.ymB
            else:
                for n in range(n_reac + n_prod):
                    M[:, n_rows - n_reac - n_prod + n] *= self.ymB[n]

        # Get equilibrium distributions
        eq_dist = np.zeros_like(dens_states)
        for i in range(n_isom):
            for s in range(n_j):
                eq_dist[i, :, s] = dens_states[i, :, s] * (2 * j_list[s] + 1) * np.exp(-e_list / constants.R / self.T)
            eq_dist[i, :, :] /= sum(eq_dist[i, :, :])

        # Set initial conditions
        p0 = np.zeros([M.shape[0]], float)
        for i in range(n_isom):
            for r in range(n_grains):
                for s in range(n_j):
                    index = indices[i, r, s]
                    if indices[i, r, s] > 0:
                        p0[index] = x0[i] * eq_dist[i, r, s]
        for i in range(n_reac + n_prod):
            p0[-n_reac - n_prod + i] = x0[i + n_isom]

        # Set up ODEs
        ode = scipy.integrate.ode(residual, jacobian).set_integrator('vode', method='bdf', with_jacobian=True,
                                                                     atol=1e-16, rtol=1e-8)
        ode.set_initial_value(p0, 0.0).set_f_params(M).set_jac_params(M)

        # Generate solution
        t = np.zeros([n_time], float)
        p = np.zeros([n_time, n_isom, n_grains, n_j], float)
        x = np.zeros([n_time, n_isom + n_reac + n_prod], float)
        for m in range(n_time):
            ode.integrate(tlist[m])
            t[m] = ode.t
            for r in range(n_grains):
                for s in range(n_j):
                    for i in range(0, n_isom):
                        index = indices[i, r, s]
                        if index > 0:
                            p[m, i, r, s] += ode.y[index]
                            x[m, i] += ode.y[index]
            for n in range(n_isom, n_isom + n_reac + n_prod):
                x[m, n] = ode.y[-(n_isom + n_reac + n_prod) + n]

        return t, p, x

    def solveReducedME(self, tlist, x0):
        """
        Directly solve the reduced master equation using a stiff ODE solver. 
        Pass the output time points `tlist` in s and the initial total
        populations `x0`. Be sure to run one of the methods for generating
        :math:`k(T,P)` values before calling this method.
        Returns the times in s, population distributions for each isomer, and total
        population profiles for each configuration.
        """
        import scipy.integrate

        e_list = self.Elist
        j_list = self.Jlist

        n_isom = self.Nisom
        n_reac = self.Nreac
        n_prod = self.Nprod
        n_grains = len(e_list)
        n_j = len(j_list)
        n_time = len(tlist)

        def residual(t, y, K):
            return np.dot(K, y)

        def jacobian(t, y, K):
            return K

        ymB = self.P / constants.R / self.T
        K = self.K.copy()
        K[:, n_isom:] *= ymB

        if self.ymB is not None:
            if isinstance(self.ymB, float):
                assert n_reac <= 1
                K[:, n_isom:] *= self.ymB
            else:
                for n in range(n_reac + n_prod):
                    K[:, n_isom + n] *= self.ymB[n]

        # Set up ODEs
        ode = scipy.integrate.ode(residual, jacobian).set_integrator('vode', method='bdf', with_jacobian=True,
                                                                     atol=1e-16, rtol=1e-8)
        ode.set_initial_value(x0, 0.0).set_f_params(K).set_jac_params(K)

        # Generate solution
        t = np.zeros([n_time], float)
        p = np.zeros([n_time, n_isom, n_grains, n_j], float)
        x = np.zeros([n_time, n_isom + n_reac + n_prod], float)
        for m in range(n_time):
            ode.integrate(tlist[m])
            t[m] = ode.t
            for i in range(n_isom):
                x[m, i] = ode.y[i]
                for j in range(n_isom):
                    for r in range(n_grains):
                        for s in range(n_j):
                            p[m, i, r, s] += ode.y[j] * self.p0[i, j, r, s]
                for j in range(n_isom, n_isom + n_reac):
                    for r in range(n_grains):
                        for s in range(n_j):
                            p[m, i, r, s] += ode.y[j] * self.p0[i, j, r, s] * ymB
            for n in range(n_reac + n_prod):
                x[m, n + n_isom] = ode.y[n + n_isom]

        return t, p, x

    def printSummary(self, level=logging.INFO):
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
        logging.log(level, '========================================================================')
        logging.log(level, '{0} network information'.format(self.label))
        logging.log(level, '-' * (len(self.label) + 20))
        logging.log(level, 'Isomers:')
        for isomer in self.isomers:
            logging.log(level, '    {0!s:<48} {1:12g} kJ/mol'.format(isomer, isomer.E0 * 0.001))
        logging.log(level, 'Reactant channels:')
        for reactants in self.reactants:
            logging.log(level, '    {0!s:<48} {1:12g} kJ/mol'.format(reactants, reactants.E0 * 0.001))
        logging.log(level, 'Product channels:')
        for products in self.products:
            logging.log(level, '    {0!s:<48} {1:12g} kJ/mol'.format(products, products.E0 * 0.001))
        logging.log(level, 'Path reactions:')
        for rxn in self.pathReactions:
            logging.log(level, '    {0!s:<48} {1:12g} kJ/mol'.format(
                rxn, float(rxn.transition_state.conformer.E0.value_si * 0.001)))
        logging.log(level, 'Net reactions:')
        for rxn in self.netReactions:
            logging.log(level, '    {0!s:<48}'.format(rxn))
        logging.log(level, '========================================================================')
        logging.log(level, '')
