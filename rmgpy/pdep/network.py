#!/usr/bin/env python
# encoding: utf-8

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
This module contains the :class:`Network` class, a representation of a 
pressure-dependent unimolecular reaction network
"""

import math
import numpy
import logging

import rmgpy.constants as constants
from rmgpy.reaction import Reaction

################################################################################

class NetworkError(Exception): 
    pass

class InvalidMicrocanonicalRateError(NetworkError):
    """Used when the k(E) calculation does not give the correct kf(T) or Kc(T)"""
    def __init__(self,message, k_ratio=1.0, Keq_ratio=1.0):
        self.message = message
        self.k_ratio = k_ratio
        self.Keq_ratio = Keq_ratio
    def badness(self):
        """
        How bad is the error?
        
        Returns the max of the absolute logarithmic errors of kf and Kc
        """
        return max(abs(math.log10(self.k_ratio)), abs(math.log10(self.Keq_ratio)))

################################################################################

class Network:
    """
    A representation of a unimolecular reaction network. The attributes are:

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `isomers`               A list of the unimolecular isomers in the network
    `reactants`             A list of the bimolecular reactant channels in the network
    `products`              A list of the bimolecular product channels in the network
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
    `activeKRotor`          ``True`` if the K-rotor is treated as active, ``False`` if treated as adiabatic
    `activeJRotor`          ``True`` if the J-rotor is treated as active, ``False`` if treated as adiabatic
    `rmgmode`               ``True`` if in RMG mode, ``False`` otherwise
    ======================= ====================================================
    
    """
    
    def __init__(self, label='', isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None):
        self.label = label
        self.isomers = isomers or []
        self.reactants = reactants or []
        self.products = products or []
        self.pathReactions = pathReactions or []
        self.bathGas = bathGas or {}
        self.netReactions = []
        
        self.T = 0.0
        self.P = 0.0
        self.Elist = None
        self.Jlist = None
        
        self.Nisom = len(self.isomers)
        self.Nreac = len(self.reactants)
        self.Nprod = len(self.products)
        self.Ngrains = 0
        self.NJ = 0

        self.activeKRotor = True
        self.activeJRotor = True

        self.grainSize = 0.0
        self.grainCount = 0
        self.E0 = None

        self.valid = False

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
        speciesList = []
        for isomer in self.isomers:
            for spec in isomer.species:
                if spec not in speciesList: speciesList.append(spec)
        for reactant in self.reactants:
            for spec in reactant.species:
                if spec not in speciesList: speciesList.append(spec)
        for product in self.products:
            for spec in product.species:
                if spec not in speciesList: speciesList.append(spec)
        for spec in self.bathGas:
            if spec not in speciesList: speciesList.append(spec)
        return speciesList

    def initialize(self, Tmin, Tmax, Pmin, Pmax, maximumGrainSize=0.0, minimumGrainCount=0, activeJRotor=True, activeKRotor=True, rmgmode=False):
        """
        Initialize a pressure dependence calculation by computing several
        quantities that are independent of the conditions. You must specify
        the temperature and pressure ranges of interesting using `Tmin` and
        `Tmax` in K and `Pmin` and `Pmax` in Pa. You must also specify the
        maximum energy grain size `grainSize` in J/mol and/or the minimum
        number of grains `grainCount`.
        """
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
        self.E0 = numpy.zeros((self.Nisom+self.Nreac+self.Nprod), numpy.float64)
        for i in range(self.Nisom):
            self.E0[i] = self.isomers[i].E0
        for n in range(self.Nreac):
            self.E0[n+self.Nisom] = self.reactants[n].E0
        for n in range(self.Nprod):
            self.E0[n+self.Nisom+self.Nreac] = self.products[n].E0
        
        # Calculate densities of states
        self.activeJRotor = activeJRotor
        self.activeKRotor = activeKRotor
        self.rmgmode = rmgmode
        
        self.calculateDensitiesOfStates()

    def calculateRateCoefficients(self, Tlist, Plist, method, errorCheck=True):
        
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)

        for rxn in self.pathReactions:
            if len(rxn.transitionState.conformer.modes) > 0:
                logging.debug('Using RRKM theory to compute k(E) for path reaction {0}.'.format(rxn))
            elif rxn.kinetics is not None:
                logging.debug('Using ILT method to compute k(E) for path reaction {0}.'.format(rxn))
        logging.debug('')
        
        logging.info('Calculating phenomenological rate coefficients for {0}...'.format(rxn))
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
                        Keq0 = K[t,p,j,i] / K[t,p,i,j]
                        Keq = eqRatios[j] / eqRatios[i]
                        if Keq0 / Keq < 0.5 or Keq0 / Keq > 2.0:
                            if i < Nisom:
                                reactants = self.isomers[i]
                            elif i < Nisom+Nreac:
                                reactants = self.reactants[i-Nisom]
                            else:
                                reactants = self.products[i-Nisom-Nreac]
                            if j < Nisom:
                                products = self.isomers[j]
                            elif j < Nisom+Nreac:
                                products = self.reactants[j-Nisom]
                            else:
                                products = self.products[j-Nisom-Nreac]
                            reaction = Reaction(reactants=reactants.species[:], products=products.species[:])
                            logging.error('For net reaction {0!s}:'.format(reaction))
                            logging.error('Expected Keq({1:g} K, {2:g} bar) = {0:11.3e}'.format(Keq, T, P*1e-5))
                            logging.error('  Actual Keq({1:g} K, {2:g} bar) = {0:11.3e}'.format(Keq0, T, P*1e-5))
                            raise NetworkError('Computed k(T,P) values for reaction {0!s} do not satisfy macroscopic equilibrium.'.format(reaction))
                            
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
                                self.K = 0 * self.K

        return K

    def setConditions(self, T, P, ymB=None):
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
        self.ymB = ymB
        
        Nisom = self.Nisom
        Nreac = self.Nreac
        Nprod = self.Nprod
        
        E0 = self.E0
        grainSize = self.grainSize
        grainCount = self.grainCount
        
        success = False
        previous_error = None
        while not success:
            success = True # (set it to false again later if necessary)
            # Update parameters that depend on temperature only if necessary
            if temperatureChanged:
                
                # Choose the energy grains to use to compute k(T,P) values at this temperature
                Elist = self.Elist = self.selectEnergyGrains(T, grainSize, grainCount)
                Ngrains = self.Ngrains = len(Elist)
                logging.info('Using {0:d} grains from {1:.2f} to {2:.2f} kJ/mol in steps of {3:.2f} kJ/mol to compute the k(T,P) values at {4:g} K'.format(
                    Ngrains, numpy.min(Elist) * 0.001, numpy.max(Elist) * 0.001, (Elist[1] - Elist[0]) * 0.001, T))
                
                # Choose the angular momenta to use to compute k(T,P) values at this temperature
                # (This only applies if the J-rotor is adiabatic
                if not self.activeJRotor:
                    Jlist = self.Jlist = numpy.arange(0, 20, 1, numpy.int)
                    NJ = self.NJ = len(Jlist)
                else:
                    Jlist = self.Jlist = numpy.array([0], numpy.int)
                    NJ = self.NJ = 1
                
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
                    if previous_error and (previous_error.message == error.message): # only compare badness if same reaction is causing problem
                        improvement = previous_error.badness()/badness
                        if improvement < 0.2 or (grainCount > 1e4 and improvement < 1.1) or (grainCount > 1.5e6): # allow it to get worse at first
                            logging.error(error.message)
                            logging.error("Increasing number of grains did not decrease error enough (Current badness: {0:.1f}, previous {1:.1f}). Something must be wrong with network {2}".format(badness,previous_error.badness(),self.label))
                            raise error
                    previous_error = error
                    success = False
                    grainSize *= 0.5
                    grainCount *= 2
                    logging.warning("Increasing number of grains, decreasing grain size and trying again. (Current badness: {0:.1f})".format(badness))
                    continue
                else:
                    success = True
                
                # Rescale densities of states such that, when they are integrated
                # using the Boltzmann factor as a weighting factor, the result is unity
                for i in range(Nisom+Nreac):
                    Q = 0.0
                    for s in range(NJ):
                        Q += numpy.sum(self.densStates[i,:,s] * (2*Jlist[s]+1) * numpy.exp(-Elist / constants.R / T))
                    self.densStates[i,:,:] /= Q
                
            # Update parameters that depend on temperature and pressure if necessary
            if temperatureChanged or pressureChanged:
                self.calculateCollisionModel()

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
        Emin = numpy.min(self.E0)
        Emin = math.floor(Emin) # Round to nearest whole number

        # Use the highest energy on the PES as the initial guess for Emax0
        Emax = numpy.max(self.E0)
        for rxn in self.pathReactions:
            E0 = float(rxn.transitionState.conformer.E0.value_si)
            if E0 > Emax: Emax = E0
        
        # Choose the actual Emax as many kB * T above the maximum energy on the PES
        # You should check that this is high enough so that the Boltzmann distributions have trailed off to negligible values
        Emax += 40. * constants.R * T

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
        
        Nisom = self.Nisom
        Nreac = self.Nreac
        Nprod = self.Nprod
        
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
            Ngrains, Elist[0] * 0.001, Elist[-1] * 0.001, dE * 0.001))
        
        # Shift the energy grains so that the minimum grain is zero
        Elist -= Elist[0]
    
        densStates = numpy.zeros((Nisom+Nreac+Nprod, Ngrains), numpy.float64)
        
        # Densities of states for isomers
        for i in range(Nisom):
            logging.debug('Calculating density of states for isomer "{0}"'.format(self.isomers[i]))
            self.isomers[i].calculateDensityOfStates(Elist, activeKRotor=self.activeKRotor, activeJRotor=self.activeJRotor, rmgmode=self.rmgmode)
        
        # Densities of states for reactant channels
        for n in range(Nreac):
            if self.reactants[n].hasStatMech():
                logging.debug('Calculating density of states for reactant channel "{0}"'.format(self.reactants[n]))
                self.reactants[n].calculateDensityOfStates(Elist, activeKRotor=self.activeKRotor, activeJRotor=self.activeJRotor, rmgmode=self.rmgmode)
            else:
                logging.debug('NOT calculating density of states for reactant channel "{0}"'.format(self.reactants[n]))
            
        # Densities of states for product channels
        if not self.rmgmode:
            for n in range(Nprod):
                if self.products[n].hasStatMech():
                    logging.debug('Calculating density of states for product channel "{0}"'.format(self.products[n]))
                    self.products[n].calculateDensityOfStates(Elist, activeKRotor=self.activeKRotor, activeJRotor=self.activeJRotor, rmgmode=self.rmgmode)
                else:
                    logging.debug('NOT calculating density of states for product channel "{0}"'.format(self.products[n]))

        logging.debug('')

#        import pylab
#        for i in range(Nisom):
#            pylab.semilogy(Elist*0.001, self.isomers[i].densStates)
#        for n in range(Nreac):
#            if self.reactants[n].densStates is not None:
#                pylab.semilogy(Elist*0.001, self.reactants[n].densStates)
#        for n in range(Nprod):
#            if self.products[n].densStates is not None:
#                pylab.semilogy(Elist*0.001, self.products[n].densStates)
#        pylab.show()

    def mapDensitiesOfStates(self):
        """
        Map the overall densities of states to the current energy grains.
        Semi-logarithmic interpolation will be used if the grain sizes of 
        `Elist0` and `Elist` do not match; this should not be a significant
        source of error as long as the grain sizes are sufficiently small.
        """
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        Ngrains = len(self.Elist)
        NJ = len(self.Jlist)
        
        self.densStates = numpy.zeros((Nisom+Nreac+Nprod, Ngrains, NJ))
        # Densities of states for isomers
        for i in range(Nisom):
            logging.debug('Mapping density of states for isomer "{0}"'.format(self.isomers[i]))
            self.densStates[i,:,:] = self.isomers[i].mapDensityOfStates(self.Elist, self.Jlist)
        # Densities of states for reactant channels
        for n in range(Nreac):
            if self.reactants[n].densStates is not None:
                logging.debug('Mapping density of states for reactant channel "{0}"'.format(self.reactants[n]))
                self.densStates[n+Nisom,:,:] = self.reactants[n].mapDensityOfStates(self.Elist, self.Jlist)
        # Densities of states for product channels
        for n in range(Nprod):
            if self.products[n].densStates is not None:
                logging.debug('Mapping density of states for product channel "{0}"'.format(self.products[n]))
                self.densStates[n+Nisom+Nreac,:,:] = self.products[n].mapDensityOfStates(self.Elist, self.Jlist)

#        import pylab
#        for i in range(Nisom+Nreac+Nprod):
#            pylab.semilogy(self.Elist*0.001, self.densStates[i,:])
#        pylab.show()

    def calculateMicrocanonicalRates(self):
        """
        Calculate and return arrays containing the microcanonical rate
        coefficients :math:`k(E)` for the isomerization, dissociation, and
        association path reactions in the network.
        """

        T = self.T
        Elist = self.Elist
        Jlist = self.Jlist
        densStates = self.densStates
        Ngrains = len(Elist)
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        NJ = 1 if self.activeJRotor else len(Jlist)

        self.Kij = numpy.zeros([Nisom,Nisom,Ngrains,NJ], numpy.float64)
        self.Gnj = numpy.zeros([Nreac+Nprod,Nisom,Ngrains,NJ], numpy.float64)
        self.Fim = numpy.zeros([Nisom,Nreac,Ngrains,NJ], numpy.float64)

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
                prod = reactants.index(rxn.products) + Nisom
            elif rxn.reactants[0] in isomers and rxn.products in products:
                # Dissociation (irreversible)
                reac = isomers.index(rxn.reactants[0])
                prod = products.index(rxn.products) + Nisom + Nreac
            elif rxn.reactants in reactants and rxn.products[0] in isomers:
                # Association (reversible)
                reac = reactants.index(rxn.reactants) + Nisom
                prod = isomers.index(rxn.products[0])
            elif rxn.reactants in products and rxn.products[0] in isomers:
                # Association (irreversible)
                reac = products.index(rxn.reactants) + Nisom + Nreac
                prod = isomers.index(rxn.products[0])
            else:
                raise NetworkError('Unexpected type of path reaction "{0}"'.format(rxn))
        
            # Compute the microcanonical rate coefficient k(E)
            reacDensStates = densStates[reac,:,:]
            prodDensStates = densStates[prod,:,:]
            kf, kr = rxn.calculateMicrocanonicalRateCoefficient(self.Elist, self.Jlist, reacDensStates, prodDensStates, T)
                        
            # Check for NaN (just to be safe)
            if numpy.isnan(kf).any() or numpy.isnan(kr).any():
                raise NetworkError('One or more k(E) values is NaN for path reaction "{0}".'.format(rxn))

            # Determine the expected value of the rate coefficient k(T)
            if rxn.canTST():
                # RRKM theory was used to compute k(E), so use TST to compute k(T)
                kf_expected = rxn.calculateTSTRateCoefficient(T)
            else:
                # ILT was used to compute k(E), so use high-P kinetics to compute k(T)
                kf_expected = rxn.kinetics.getRateCoefficient(T)
            
            # Determine the expected value of the equilibrium constant (Kc)
            Keq_expected = self.eqRatios[prod] / self.eqRatios[reac] 

            # Determine the actual values of k(T) and Keq
            C0 = 1e5 / (constants.R * T)
            kf0 = 0.0; kr0 = 0.0; Qreac = 0.0; Qprod = 0.0
            for s in range(NJ):
                kf0 += numpy.sum(kf[:,s] * reacDensStates[:,s] * (2*Jlist[s]+1) * numpy.exp(-Elist / constants.R / T)) 
                kr0 += numpy.sum(kr[:,s] * prodDensStates[:,s] * (2*Jlist[s]+1) * numpy.exp(-Elist / constants.R / T)) 
                Qreac += numpy.sum(reacDensStates[:,s] * (2*Jlist[s]+1) * numpy.exp(-Elist / constants.R / T)) 
                Qprod += numpy.sum(prodDensStates[:,s] * (2*Jlist[s]+1) * numpy.exp(-Elist / constants.R / T)) 
            kr0 *= C0 ** (len(rxn.products) - len(rxn.reactants))
            Qprod *= C0 ** (len(rxn.products) - len(rxn.reactants))
            kf_actual = kf0 / Qreac if Qreac > 0 else 0
            kr_actual = kr0 / Qprod if Qprod > 0 else 0
            Keq_actual = kf_actual / kr_actual if kr_actual > 0 else 0

            error = False; warning = False
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
                self.Kij[prod,reac,:,:] = kf
                self.Kij[reac,prod,:,:] = kr
            elif rxn.reactants[0] in isomers and rxn.products in reactants:
                # Dissociation (reversible)
                self.Gnj[prod-Nisom,reac,:,:] = kf
                self.Fim[reac,prod-Nisom,:,:] = kr
            elif rxn.reactants[0] in isomers and rxn.products in products:
                # Dissociation (irreversible)
                self.Gnj[prod-Nisom,reac,:,:] = kf
            elif rxn.reactants in reactants and rxn.products[0] in isomers:
                # Association (reversible)
                self.Fim[prod,reac-Nisom,:,:] = kf 
                self.Gnj[reac-Nisom,prod,:,:] = kr
            elif rxn.reactants in products and rxn.products[0] in isomers:
                # Association (irreversible)
                self.Gnj[reac-Nisom,prod,:,:] = kr
            else:
                raise NetworkError('Unexpected type of path reaction "{0}"'.format(rxn))

            # If the k(E) values are invalid (in that they give the wrong 
            # kf(T) or kr(T) when integrated), then raise an exception
            if error or warning:
                logging.warning('For path reaction {0!s}:'.format(rxn))
                logging.warning('    Expected kf({0:g} K) = {1:g}'.format(T, kf_expected))
                logging.warning('      Actual kf({0:g} K) = {1:g}'.format(T, kf_actual))
                logging.warning('    Expected Keq({0:g} K) = {1:g}'.format(T, Keq_expected))
                logging.warning('      Actual Keq({0:g} K) = {1:g}'.format(T, Keq_actual))
                if error:
                    raise InvalidMicrocanonicalRateError('Invalid k(E) values computed for path reaction "{0}".'.format(rxn), k_ratio, Keq_ratio)
                else:
                    logging.warning('Significant corrections to k(E) to be consistent with high-pressure limit for path reaction "{0}".'.format(rxn))

#        import pylab
#        for prod in range(Nisom):
#            for reac in range(prod):
#                pylab.semilogy(self.Elist*0.001, self.Kij[prod,reac,:])
#        for prod in range(Nreac+Nprod):
#            for reac in range(Nisom):
#                pylab.semilogy(self.Elist*0.001, self.Gnj[prod,reac,:])
#        pylab.show()

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
        T = self.T
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        eqRatios = numpy.zeros(Nisom+Nreac+Nprod, numpy.float64)
        conc = (1e5 / constants.R / T)          # [=] mol/m^3
        for i in range(Nisom):
            G = self.isomers[i].getFreeEnergy(T)
            eqRatios[i] = math.exp(-G / constants.R / T)
        for i in range(Nreac):
            G = self.reactants[i].getFreeEnergy(T)
            eqRatios[Nisom+i] = math.exp(-G / constants.R / T) * conc ** (len(self.reactants[i].species) - 1)
        for i in range(Nprod):
            if self.products[i].hasStatMech() or self.products[i].hasThermo():
                G = self.products[i].getFreeEnergy(T)
                eqRatios[Nisom+Nreac+i] = math.exp(-G / constants.R / T) * conc ** (len(self.products[i].species) - 1)
        self.eqRatios = eqRatios
        return eqRatios / numpy.sum(eqRatios)

    def calculateCollisionModel(self):
        """
        Calculate the matrix of first-order rate coefficients for collisional
        population transfer between grains for each isomer, including the
        corresponding collision frequencies.
        """
        Nisom = len(self.isomers)
        Ngrains = len(self.Elist)
        NJ = 1 if self.Jlist is None else len(self.Jlist)
        
        collFreq = numpy.zeros(Nisom, numpy.float64)
        Mcoll = numpy.zeros((Nisom,Ngrains,NJ,Ngrains,NJ), numpy.float64)
        
        for i, isomer in enumerate(self.isomers):
            collFreq[i] = isomer.calculateCollisionFrequency(self.T, self.P, self.bathGas)
            Mcoll[i,:,:,:,:] = collFreq[i] * isomer.generateCollisionMatrix(self.T, self.densStates[i,:,:], self.Elist, self.Jlist)
                        
        self.collFreq = collFreq
        self.Mcoll = Mcoll
        
        return Mcoll

    def applyModifiedStrongCollisionMethod(self, efficiencyModel='default'):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the modified strong collision method.
        """
        import rmgpy.pdep.msc as msc
        logging.debug('Applying modified strong collision method at {0:g} K, {1:g} bar...'.format(self.T, self.P))
        self.K, self.p0 = msc.applyModifiedStrongCollisionMethod(self, efficiencyModel)
        return self.K, self.p0
    
    def applyReservoirStateMethod(self):
        """
        Compute the phenomenological rate coefficients :math:`k(T,P)` at the
        current conditions using the reservoir state method.
        """
        import rmgpy.pdep.rs as rs
        logging.debug('Applying reservoir state method at {0:g} K, {1:g} bar...'.format(self.T, self.P))
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
        logging.debug('Applying chemically-significant eigenvalues method at {0:g} K, {1:g} bar...'.format(self.T, self.P))
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
    
        Elist = self.Elist
        Jlist = self.Jlist
        densStates = self.densStates
        
        Nisom = self.Nisom
        Nreac = self.Nreac
        Nprod = self.Nprod
        Ngrains = len(Elist)
        NJ = len(Jlist)
        Ntime = len(tlist)
        
        def residual(t, y, K):
            return numpy.dot(K, y)
        
        def jacobian(t, y, K):
            return K
    
        ymB = self.P / constants.R / self.T
        M, indices = self.generateFullMEMatrix()
        Nrows = M.shape[0]
        M[:,Nrows-Nreac-Nprod:] *= ymB
        
        if self.ymB is not None:
            if isinstance(self.ymB, float):
                assert Nreac <= 1
                M[:,Nrows-Nreac-Nprod:] *= self.ymB
            else:
                for n in range(Nreac+Nprod):
                    M[:,Nrows-Nreac-Nprod+n] *= self.ymB[n]
        
        # Get equilibrium distributions
        eqDist = numpy.zeros_like(densStates)
        for i in range(Nisom):
            for s in range(NJ):
                eqDist[i,:,s] = densStates[i,:,s] * (2*Jlist[s]+1) * numpy.exp(-Elist / constants.R / self.T)
            eqDist[i,:,:] /= sum(eqDist[i,:,:])
    
        # Set initial conditions
        p0 = numpy.zeros([M.shape[0]], float)
        for i in range(Nisom):
            for r in range(Ngrains):
                for s in range(NJ):
                    index = indices[i,r,s]
                    if indices[i,r,s] > 0:
                        p0[index] = x0[i] * eqDist[i,r,s]
        for i in range(Nreac+Nprod):
            p0[-Nreac-Nprod + i] = x0[i+Nisom]
    
        # Set up ODEs
        ode = scipy.integrate.ode(residual, jacobian).set_integrator('vode', method='bdf', with_jacobian=True, atol=1e-16, rtol=1e-8)
        ode.set_initial_value(p0, 0.0).set_f_params(M).set_jac_params(M)
    
        # Generate solution
        t = numpy.zeros([Ntime], float)
        p = numpy.zeros([Ntime, Nisom, Ngrains, NJ], float)
        x = numpy.zeros([Ntime, Nisom+Nreac+Nprod], float)
        for m in range(Ntime):
            ode.integrate(tlist[m])
            t[m] = ode.t
            for r in range(Ngrains):
                for s in range(NJ):
                    for i in range(0, Nisom):
                        index = indices[i,r,s]
                        if index > 0:
                            p[m,i,r,s] += ode.y[index]
                            x[m,i] += ode.y[index]
            for n in range(Nisom, Nisom+Nreac+Nprod):
                x[m,n] = ode.y[-(Nisom+Nreac+Nprod)+n]
    
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
    
        Elist = self.Elist
        Jlist = self.Jlist
        
        Nisom = self.Nisom
        Nreac = self.Nreac
        Nprod = self.Nprod
        Ngrains = len(Elist)
        NJ = len(Jlist)
        Ntime = len(tlist)
        
        def residual(t, y, K):
            return numpy.dot(K, y)
        
        def jacobian(t, y, K):
            return K
    
        ymB = self.P / constants.R / self.T
        K = self.K.copy()
        K[:,Nisom:] *= ymB
        
        if self.ymB is not None:
            if isinstance(self.ymB, float):
                assert Nreac <= 1
                K[:,Nisom:] *= self.ymB
            else:
                for n in range(Nreac+Nprod):
                    K[:,Nisom+n] *= self.ymB[n]
        
        # Set up ODEs
        ode = scipy.integrate.ode(residual, jacobian).set_integrator('vode', method='bdf', with_jacobian=True, atol=1e-16, rtol=1e-8)
        ode.set_initial_value(x0, 0.0).set_f_params(K).set_jac_params(K)
    
        # Generate solution
        t = numpy.zeros([Ntime], float)
        p = numpy.zeros([Ntime, Nisom, Ngrains, NJ], float)
        x = numpy.zeros([Ntime, Nisom+Nreac+Nprod], float)
        for m in range(Ntime):
            ode.integrate(tlist[m])
            t[m] = ode.t
            for i in range(Nisom):
                x[m,i] = ode.y[i]
                for j in range(Nisom):
                    for r in range(Ngrains):
                        for s in range(NJ):
                            p[m,i,r,s] += ode.y[j] * self.p0[i,j,r,s]
                for j in range(Nisom, Nisom+Nreac):
                    for r in range(Ngrains):
                        for s in range(NJ):
                            p[m,i,r,s] += ode.y[j] * self.p0[i,j,r,s] * ymB
            for n in range(Nreac+Nprod):
                x[m,n+Nisom] = ode.y[n+Nisom]
                
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
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(str(isomer), isomer.E0*0.001))
        logging.log(level, 'Reactant channels:')
        for reactants in self.reactants:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(str(reactants), reactants.E0*0.001))
        logging.log(level, 'Product channels:')
        for products in self.products:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(str(products), products.E0*0.001))
        logging.log(level, 'Path reactions:')
        for rxn in self.pathReactions:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(rxn, float(rxn.transitionState.conformer.E0.value_si*0.001)))
        logging.log(level, '========================================================================')
        logging.log(level, '')
