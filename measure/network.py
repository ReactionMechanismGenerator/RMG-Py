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

import chempy.constants as constants
import chempy.states as states

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
    `explored`          ``list``                A list of the unimolecular isomers whose reactions have been fully explored
    =================== ======================= ================================

    """

    def __init__(self, index=-1, isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None, collisionModel=None):
        self.index = index
        self.isomers = isomers or []
        self.reactants = reactants or []
        self.products = products or []
        self.pathReactions = pathReactions or []
        self.bathGas = bathGas or {}
        self.collisionModel = collisionModel
        self.netReactions = []
        self.valid = False
        self.explored = []

    def invalidate(self):
        """
        Mark the network as in need of a new calculation to determine the
        pressure-dependent rate coefficients
        """
        self.valid = False

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

        # Densities of states for isomers
        for i in range(Nisom):
            logging.debug('Calculating density of states for isomer "%s"' % self.isomers[i])
            densStates0 = self.isomers[i].states.getDensityOfStates(Elist)
            # Shift to common zero of energy
            r0 = int(round(E0[i] / dE))
            densStates[i,r0:] = densStates0[:-r0+len(densStates0)]

        # Densities of states for reactant channels
        # (Only if not minimizing the number of density of states calculations)
        if not settings.minimizeDensityOfStatesCalculations:
            for n in range(Nreac):
                r0 = int(round(E0[n+Nisom] / dE))
                if self.reactants[n][0].states is not None and self.reactants[n][1].states is not None:
                    logging.debug('Calculating density of states for reactant channel "%s"' % (' + '.join([str(spec) for spec in self.reactants[n]])))
                    densStates0 = self.reactants[n][0].states.getDensityOfStates(Elist)
                    densStates1 = self.reactants[n][1].states.getDensityOfStates(Elist)
                    densStates0 = states.convolve(densStates0, densStates1, Elist)
                    # Shift to common zero of energy
                    densStates[n+Nisom,r0:] = densStates0[:-r0+len(densStates0)]
                elif self.reactants[n][0].states is not None:
                    logging.debug('Calculating density of states for reactant channel "%s"' % (' + '.join([str(spec) for spec in self.reactants[n]])))
                    densStates0 = self.reactants[n][0].states.getDensityOfStates(Elist)
                    # Shift to common zero of energy
                    densStates[n+Nisom,r0:] = densStates0[:-r0+len(densStates0)]
                elif self.reactants[n][1].states is not None:
                    logging.debug('Calculating density of states for reactant channel "%s"' % (' + '.join([str(spec) for spec in self.reactants[n]])))
                    densStates0 = self.reactants[n][1].states.getDensityOfStates(Elist)
                    # Shift to common zero of energy
                    densStates[n+Nisom,r0:] = densStates0[:-r0+len(densStates0)]
                else:
                    logging.debug('NOT calculating density of states for reactant channel "%s"' % (' + '.join([str(spec) for spec in self.reactants[n]])))
            logging.debug('')

        return densStates

    def calculateMicrocanonicalRates(self, Elist, densStates, T=None):
        """
        Calculate and return arrays containing the microcanonical rate
        coefficients :math:`k(E)` for the isomerization, dissociation, and
        association path reactions in the network. `Elist` represents the
        array of energies in J/mol at which to compute each density of states,
        while `densStates` represents the density of states of each isomer and
        reactant channel in mol/J. The temperature `T` is K is used in certain
        circumstances when :math:`k(E)` cannot be determined without it, and
        in the detailed balance expression to obtain the reverse kinetics.
        """

        Ngrains = len(Elist)
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)

        Kij = numpy.zeros([Nisom,Nisom,Ngrains], numpy.float64)
        Gnj = numpy.zeros([Nreac+Nprod,Nisom,Ngrains], numpy.float64)
        Fim = numpy.zeros([Nisom,Nreac,Ngrains], numpy.float64)

        for rxn in self.pathReactions:
            if rxn.reactants[0] in self.isomers and rxn.products[0] in self.isomers:
                # Isomerization
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.isomers.index(rxn.products[0])
                Kij[prod,reac,:], Kij[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], densStates[prod,:], T)
            elif rxn.reactants[0] in self.isomers and rxn.products in self.reactants:
                # Dissociation (reversible)
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.reactants.index(rxn.products)
                Gnj[prod,reac,:], Fim[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], densStates[prod+Nisom,:], T)
            elif rxn.reactants[0] in self.isomers and rxn.products in self.products:
                # Dissociation (irreversible)
                reac = self.isomers.index(rxn.reactants[0])
                prod = self.products.index(rxn.products) + Nreac
                Gnj[prod,reac,:], dummy = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac,:], None, T)
            elif rxn.reactants in self.reactants and rxn.products[0] in self.isomers:
                # Association
                reac = self.reactants.index(rxn.reactants)
                prod = self.isomers.index(rxn.products[0])
                Fim[prod,reac,:], Gnj[reac,prod,:] = calculateMicrocanonicalRateCoefficient(rxn, Elist, densStates[reac+Nisom,:], densStates[prod,:], T)
            else:
                raise NetworkError('Unexpected type of path reaction "%s"' % rxn)
        logging.debug('')

        return Kij, Gnj, Fim

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
        logging.log(level, 'Network Information')
        logging.log(level, '-------------------')
        logging.log(level, 'Isomers:')
        for isomer in self.isomers:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(str(isomer), isomer.E0 / 1000.0))
        logging.log(level, 'Reactant channels:')
        for reactants in self.reactants:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(' + '.join([str(spec) for spec in reactants]), sum([spec.E0 for spec in reactants]) / 1000.0))
        logging.log(level, 'Product channels:')
        for products in self.products:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(' + '.join([str(spec) for spec in products]), sum([spec.E0 for spec in products]) / 1000.0))
        logging.log(level, 'Path reactions:')
        for rxn in self.pathReactions:
            logging.log(level, '    {0:<48s} {1:12g} kJ/mol'.format(rxn, rxn.transitionState.E0 / 1000.0))
        logging.log(level, '========================================================================')
        logging.log(level, '')

    def calculateRateCoefficients(self, Tlist, Plist, Elist, method):
        """
        Calculate the phenomenological rate coefficients :math:`k(T,P)` for the
        network at the given temperatures `Tlist` in K and pressures `Plist` in
        Pa using the energy grains ``Elist`` in J/mol. The `method` string is
        used to indicate the method to use, and should be one of ``"modified
        strong collision"``, ``"reservoir state"``, or
        ``"chemically-significant eigenvalues"``.
        """

        # Determine the values of some counters
        Ngrains = len(Elist)
        Nisom = len(self.isomers)
        Nreac = len(self.reactants)
        Nprod = len(self.products)
        dE = Elist[1] - Elist[0]

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
        Emin = Elist[0]
        for rxn in self.pathReactions:
            rxn.transitionState.E0 -= Emin
        E0 -= Emin
        Ereac -= Emin
        Elist -= Emin

        # Calculate density of states for each isomer and each reactant channel
        # that has the necessary parameters
        logging.info('Calculating densities of states for network %i...' % self.index)
        densStates0 = self.calculateDensitiesOfStates(Elist, E0)

        logging.info('Calculating phenomenological rate coefficients for network %i...' % self.index)
        K = numpy.zeros((len(Tlist),len(Plist),Nisom+Nreac+Nprod,Nisom+Nreac+Nprod), numpy.float64)

        for t, T in enumerate(Tlist):

            # Calculate microcanonical rate coefficients for each path reaction
            # If degree of freedom data is provided for the transition state, then RRKM theory is used
            # If high-pressure limit Arrhenius data is provided, then the inverse Laplace transform method is used
            # Otherwise an exception is raised
            # This is only dependent on temperature for the ILT method with
            # certain Arrhenius parameters
            Kij, Gnj, Fim = self.calculateMicrocanonicalRates(Elist, densStates0, T)

            # Rescale densities of states such that, when they are integrated
            # using the Boltzmann factor as a weighting factor, the result is unity
            densStates = numpy.zeros_like(densStates0)
            eqRatios = numpy.zeros(Nisom+Nreac, numpy.float64)
            for i in range(Nisom+Nreac):
                eqRatios[i] = numpy.sum(densStates0[i,:] * numpy.exp(-Elist / constants.R / T)) * dE
                densStates[i,:] = densStates0[i,:] / eqRatios[i] * dE
            # Use free energy to determine equilibrium ratios of each isomer and product channel
            eqRatios = numpy.zeros(Nisom+Nreac, numpy.float64)
            conc = 1e5 / constants.R / T
            for i in range(Nisom):
                G = self.isomers[i].thermo.getFreeEnergy(T)
                eqRatios[i] = math.exp(-G / constants.R / T)
            for i in range(Nreac):
                G = sum([spec.thermo.getFreeEnergy(T) for spec in self.reactants[i]])
                eqRatios[Nisom+i] = math.exp(-G / constants.R / T) * conc ** (len(self.reactants[i]) - 1)

            for p, P in enumerate(Plist):

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
                    # Apply chemically-significant eigenvalues method
                    import cse
                    K[t,p,:,:], p0 = cse.applyChemicallySignificantEigenvaluesMethod(T, P, Elist, densStates, Mcoll, Kij, Fim, Gnj, eqRatios, Nisom, Nreac, Nprod)
                else:
                    raise NetworkError('Unknown method "%s".' % method)

                logging.log(0, K[t,p,0:Nisom+Nreac+Nprod,0:Nisom+Nreac])

                logging.debug('')

        # Unshift energy grains
        for rxn in self.pathReactions:
            rxn.transitionState.E0 += Emin
        Elist += Emin

        # Mark network as valid
        self.valid = True

        return K

    def __createNewSurfaceAndContext(self, ext, fstr='', width=800, height=600):
        import cairo
        if ext == '.svg':
            surface = cairo.SVGSurface(fstr, width, height)
        elif ext == '.pdf':
            surface = cairo.PDFSurface(fstr, width, height)
        elif ext == '.ps':
            surface = cairo.PSSurface(fstr, width, height)
        else:
            logging.warning('Unknown format for target "%s"; not drawing potential energy surface.' % fstr)
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

        # Determine order of wells based on order of path reactions, but put
        # all the unimolecular isomer wells first
        wells = []
        for isomer in self.isomers: wells.append([isomer])
        for rxn in self.pathReactions:
            if rxn.reactants not in wells:
                if len(rxn.products) == 1 and rxn.products[0] in self.isomers:
                    if self.isomers.index(rxn.products[0]) < len(self.isomers) / 2:
                        wells.insert(0, rxn.reactants)
                    else:
                        wells.append(rxn.reactants)
                else:
                    wells.append(rxn.reactants)
            if rxn.products not in wells:
                if rxn.products in self.products:
                    wells.append(rxn.products)
                elif len(rxn.reactants) == 1 and rxn.reactants[0] in self.isomers:
                    if self.isomers.index(rxn.reactants[0]) < len(self.isomers) / 2:
                        wells.insert(0, rxn.products)
                    else:
                        wells.append(rxn.products)
                else:
                    wells.append(rxn.products)
        
        # Drawing parameters
        padding_left = 96.0
        padding_right = padding_left
        padding_top = padding_left / 2.0
        padding_bottom = padding_left / 2.0
        wellWidth = 64.0; wellSpacing = 64.0; Eslope = 5.0; TSwidth = 16.0
        E0 = [sum([spec.E0 for spec in well]) / 4184 for well in wells]
        E0.extend([rxn.transitionState.E0 / 4184 for rxn in self.pathReactions])
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
            E0 = sum([spec.E0 for spec in well]) / 4184
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
        height = numpy.max(coordinates[:,1]) - numpy.min(coordinates[:,1]) + 32.0 + padding_top + padding_bottom

        # Choose multiplier to convert energies to desired units
        if Eunits == 'J/mol':      Emult = 1.0
        elif Eunits == 'kJ/mol':   Emult = 1.0 / 1000
        elif Eunits == 'cal/mol':  Emult = 1.0 / 4.184
        elif Eunits == 'kcal/mol': Emult = 1.0 / 4184
        elif Eunits == 'cm^-1':    Emult = 1.0 / 11.96
        else:
            logging.warning('Invalid value "%s" for Eunits parameter. Setting to "kJ/mol".' % (Eunits))
            Emult = 1.0 / 1000

        # Initialize Cairo surface and context
        surface, cr = self.__createNewSurfaceAndContext(ext, fstr, width, height)

        # Some global settings
        cr.select_font_face("sans")
        cr.set_font_size(10)
        
        # Draw path reactions
        for rxn in self.pathReactions:
            reac = wells.index(rxn.reactants)
            prod = wells.index(rxn.products)
            E0_reac = sum([spec.E0 for spec in wells[reac]]) / 4184
            E0_prod = sum([spec.E0 for spec in wells[prod]]) / 4184
            E0_TS = rxn.transitionState.E0 / 4184
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
                E0 = "%.1f" % (rxn.transitionState.E0 * Emult)
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
            E0 = "%.1f" % (sum([spec.E0 for spec in well]) * Emult)
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
        surface.finish()

    def __drawText(self, text, ext, padding=0):
        """
        Create and return a temporary Cairo surface containing the string
        `text` with an optional amount of `padding` on all sides. The type of
        surface is dictated by the `ext` parameter.
        """

        from chempy.ext.molecule_draw import createNewSurface, fontSizeNormal
        import cairo

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

        from chempy.ext.molecule_draw import createNewSurface, drawMolecule
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
                surface0, cr0, boundingRect0 = drawMolecule(spec.molecule[0], surface=ext[1:])
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
        