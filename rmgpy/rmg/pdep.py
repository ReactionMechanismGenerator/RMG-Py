#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains classes for providing pressure-dependent kinetics estimation
functionality to RMG.
"""

import logging
import os.path

from rmgpy.quantity import Quantity
import rmgpy.measure.network
import rmgpy.reaction

################################################################################

class PressureDependenceError(Exception):
    """
    An exception class to use when an error involving pressure dependence is
    encountered. Pass a string describing the circumstances of the exceptional
    behavior.
    """
    pass

################################################################################

class PDepReaction(rmgpy.reaction.Reaction):

    def __init__(self, index=-1, reactants=None, products=None, network=None, kinetics=None, reversible=True, transitionState=None, thirdBody=False):
        rmgpy.reaction.Reaction.__init__(self, index, reactants, products, kinetics, reversible, transitionState, thirdBody)
        self.network = network

################################################################################

class PDepNetwork(rmgpy.measure.network.Network):
    """
    A representation of a *partial* unimolecular reaction network. Each partial
    network has a single `source` isomer or reactant channel, and is responsible
    only for :math:`k(T,P)` values for net reactions with source as the
    reactant. Multiple partial networks can have the same source, but networks
    with the same source and any explored isomers must be combined.

    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `source`            ``list``                The isomer or reactant channel that acts as the source
    `explored`          ``list``                A list of the unimolecular isomers whose reactions have been fully explored
    =================== ======================= ================================

    """

    def __init__(self, index=-1, source=None):
        rmgpy.measure.network.Network.__init__(self, index=index)
        self.source = source
        self.explored = []

    def getLeakCoefficient(self, T, P):
        """
        Return the pressure-dependent rate coefficient :math:`k(T,P)` describing
        the total rate of "leak" from this network. This is defined as the sum
        of the :math:`k(T,P)` values for all net reactions to nonexplored
        unimolecular isomers.
        """
        k = 0.0
        if len(self.netReactions) == 0 and len(self.pathReactions) == 1:
            # The network is of the form A + B -> C* (with C* nonincluded)
            # For this special case we use the high-pressure limit k(T) to
            # ensure that we're estimating the total leak flux
            rxn = self.pathReactions[0]
            if rxn.kinetics is None:
                if rxn.reverse.kinetics is not None:
                    rxn = rxn.reverse
                else:
                    raise PressureDependenceError('Path reaction {0} with no high-pressure-limit kinetics encountered in PDepNetwork #{1:d} while evaluating leak flux.'.format(rxn, self.index))
            if rxn.products is self.source:
                k = rxn.getRateCoefficient(T,P) / rxn.getEquilibriumConstant(T)
            else:
                k = rxn.getRateCoefficient(T,P)
        else:
            # The network has at least one included isomer, so we can calculate
            # the leak flux normally
            for rxn in self.netReactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    k += rxn.getRateCoefficient(T,P)
        return k

    def getMaximumLeakSpecies(self, T, P):
        """
        Get the unexplored (unimolecular) isomer with the maximum leak flux.
        Note that the leak rate coefficients vary with temperature and
        pressure, so you must provide these in order to get a meaningful result.
        """
        # Choose species with maximum leak flux
        maxK = 0.0; maxSpecies = None
        if len(self.netReactions) == 0 and len(self.pathReactions) == 1:
            maxK = self.getLeakCoefficient(T,P)
            rxn = self.pathReactions[0]
            if rxn.products == self.source:
                assert len(rxn.reactants) == 1
                maxSpecies = rxn.reactants[0]
            else:
                assert len(rxn.products) == 1
                maxSpecies = rxn.products[0]
        else:
            for rxn in self.netReactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    k = rxn.getRateCoefficient(T,P)
                    if maxSpecies is None or k > maxK:
                        maxSpecies = rxn.products[0]
                        maxK = k

        # Make sure we've identified a species
        if maxSpecies is None:
            raise UnirxnNetworkException('No unimolecular isomers left to explore!')
        # Return the species
        return maxSpecies

    def getLeakBranchingRatios(self, T, P):
        """
        Return a dict with the unexplored isomers in the partial network as the
        keys and the fraction of the total leak coefficient as the values.
        """
        ratios = {}
        if len(self.netReactions) == 0 and len(self.pathReactions) == 1:
            rxn = self.pathReactions[0]
            assert rxn.reactants == self.source or rxn.products == self.source
            if rxn.products == self.source:
                assert len(rxn.reactants) == 1
                ratios[rxn.reactants[0]] = 1.0
            else:
                assert len(rxn.products) == 1
                ratios[rxn.products[0]] = 1.0
        else:
            for rxn in self.netReactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    ratios[rxn.products[0]] = rxn.getRateCoefficient(T,P)

        kleak = sum(ratios.values())
        for spec in ratios:
            ratios[spec] /= kleak

        return ratios

    def exploreIsomer(self, isomer, reactionModel, database):
        """
        Explore a previously-unexplored unimolecular `isomer` in this partial
        network using the provided core-edge reaction model `reactionModel`,
        returning the new reactions and new species.
        """
        assert isomer not in self.explored
        assert isomer not in self.isomers
        assert isomer not in self.source

        logging.info('Exploring isomer {0} in pressure-dependent network #{1:d}'.format(isomer, self.index))
        self.explored.append(isomer)
        # Find reactions involving the found species as unimolecular
        # reactant or product (e.g. A <---> products)
        newReactionList = reactionModel.react(database, isomer)
        # Don't find reactions involving the new species as bimolecular
        # reactants or products with itself (e.g. A + A <---> products)
        # Don't find reactions involving the new species as bimolecular
        # reactants or products with other core species (e.g. A + B <---> products)

        return newReactionList

    def addPathReaction(self, newReaction, newSpecies):
        """
        Add a path reaction to the network. If the path reaction already exists,
        no action is taken.
        """
        # Add this reaction to that network if not already present
        found = False
        for rxn in self.pathReactions:
            if newReaction.reactants == rxn.reactants and newReaction.products == rxn.products:
                found = True
                break
            elif newReaction.products == rxn.reactants and newReaction.reactants == rxn.products:
                found = True
                break
        if not found:
            self.pathReactions.append(newReaction)
            self.invalidate()

    def merge(self, other):
        """
        Merge the partial network `other` into this network.
        """
        # Make sure the two partial networks have the same source configuration
        assert self.source == other.source

        # Merge isomers
        for isomer in other.isomers:
            if isomer not in self.isomers:
                self.isomers.append(isomer)
        # Merge explored
        for isomer in other.explored:
            if isomer not in self.explored:
                self.explored.append(isomer)
        # Merge reactants
        for reactants in other.reactants:
            if reactants not in self.reactants:
                self.reactants.append(reactants)
        # Merge products
        for products in other.products:
            if products not in self.products:
                self.products.append(products)
        
        # However, products that have been explored are actually isomers
        # These should be removed from the list of products!
        productsToRemove = []
        for products in self.products:
            if len(products) == 1 and products[0] in self.isomers:
                productsToRemove.append(products)
        for products in productsToRemove:
            self.products.remove(products)

        # Merge path reactions
        for reaction in other.pathReactions:
            found = False
            for rxn in self.pathReactions:
                if reaction.reactants == rxn.reactants and reaction.products == rxn.products:
                    # NB the isEquivalent() method that used to be on the previous line also checked reverse direction.
                    # I am not sure which is appropriate 
                    found = True
                    break
            if not found:
                self.pathReactions.append(reaction)

        # Also merge net reactions (so that when we update the network in the
        # future, we update the existing net reactions rather than making new ones)
        # Q: What to do when a net reaction exists in both networks being merged?
        for reaction in other.netReactions:
            found = False
            for rxn in self.netReactions:
                if reaction.reactants == rxn.reactants and reaction.products == rxn.products:
                    # NB the isEquivalent() method that used to be on the previous line also checked reverse direction.
                    # I am not sure which is appropriate 
                    found = True
                    break
            if not found:
                self.netReactions.append(reaction)

        # Mark this network as invalid
        self.valid = False

    def updateConfigurations(self):
        """
        Regenerate the :math:`k(T,P)` values for this partial network if the
        network is marked as invalid.
        """

        # Figure out which configurations are isomers, reactant channels, and product channels
        self.products = []
        if len(self.source) == 1:
            # The source is a unimolecular isomer
            self.isomers = self.explored[:]
            if self.source[0] not in self.isomers: self.isomers.insert(0, self.source[0])
            self.reactants = []
        else:
            # The source is a bimolecular reactant channel
            self.source.sort()
            self.isomers = self.explored[:]
            self.reactants = [self.source]
        for rxn in self.pathReactions:
            # Sort bimolecular configurations so that we always encounter them in the
            # same order
            # The actual order doesn't matter, as long as it is consistent
            rxn.reactants.sort()
            rxn.products.sort()
            # All reactants and products not already assigned as an isomer or
            # a reactant channel are product channels
            if rxn.reactants not in self.products:
                if len(rxn.reactants) == 1 and rxn.reactants[0] not in self.isomers:
                    self.products.append(rxn.reactants)
                elif len(rxn.reactants) > 1 and rxn.reactants not in self.reactants:
                    self.products.append(rxn.reactants)
            if rxn.products not in self.products:
                if len(rxn.products) == 1 and rxn.products[0] not in self.isomers:
                    self.products.append(rxn.products)
                elif len(rxn.products) > 1 and rxn.products not in self.reactants:
                    self.products.appende(rxn.products)

    def update(self, reactionModel, database, pdepSettings):
        """
        Regenerate the :math:`k(T,P)` values for this partial network if the
        network is marked as invalid.
        """
        from rmgpy.kinetics import Arrhenius, KineticsData
        from rmgpy.measure.collision import SingleExponentialDown
        from rmgpy.measure.reaction import fitInterpolationModel
        from rmgpy.measure.main import MEASURE
        import rmgpy.measure.settings
        
        # Get the parameters for the pressure dependence calculation
        # method, Tmin, Tmax, Tlist, Pmin, Pmax, Plist, grainSize, numGrains, model = settings.pressureDependence
        method = pdepSettings['method']
        Tmin = pdepSettings['Tmin']
        Tmax = pdepSettings['Tmax']
        Tlist = pdepSettings['Tlist']
        Pmin = pdepSettings['Pmin']
        Pmax = pdepSettings['Pmax']
        Plist = pdepSettings['Plist']
        grainSize = pdepSettings['minimumGrainSize']
        numGrains = pdepSettings['minimumNumberOfGrains']
        model = pdepSettings['interpolation']
        outputDirectory = pdepSettings['outputDirectory']
        
        measure = MEASURE()
        measure.method = method
        measure.Tmin = Quantity(Tmin, "K")
        measure.Tmax = Quantity(Tmax, "K")
        measure.Tlist = Quantity(Tlist, "K")
        measure.Pmin = Quantity(Pmin/1e5, "bar")
        measure.Pmax = Quantity(Pmax/1e5, "bar")
        measure.Plist = Quantity(Plist/1e5, "bar")
        measure.grainSize = Quantity(grainSize/1000., "kJ/mol")
        measure.grainCount = numGrains
        measure.model = model
        measure.network = self

        # Figure out which configurations are isomers, reactant channels, and product channels
        self.updateConfigurations()

        # Make sure we have high-P kinetics for all path reactions
        for rxn in self.pathReactions:
            if rxn.kinetics is None and rxn.reverse.kinetics is None:
                raise PressureDependenceError('Path reaction {0} with no high-pressure-limit kinetics encountered in PDepNetwork #{1:d}.'.format(rxn, self.index))
            elif rxn.kinetics is not None and rxn.kinetics.isPressureDependent():
                raise PressureDependenceError('Pressure-dependent kinetics encountered for path reaction {0} in PDepNetwork #{1:d}.'.format(rxn, self.index))
        
        # Do nothing if the network is already valid
        if self.valid: return
        # Do nothing if there are no explored wells
        if len(self.explored) == 0 and len(self.source) > 1: return

        # Generate states data for unimolecular isomers and reactants if necessary
        for spec in self.isomers:
            if spec.states is None: spec.generateStatesData(database)
        for reactants in self.reactants:
            for spec in reactants:
                if spec.states is None: spec.generateStatesData(database)

        # Determine transition state energies on potential energy surface
        # In the absence of any better information, we simply set it to
        # be the reactant ground-state energy + the activation energy
        # Note that we need Arrhenius kinetics in order to do this
        for rxn in self.pathReactions:
            if rxn.kinetics is None:
                raise Exception('Path reaction "{0}" in PDepNetwork #{1:d} has no kinetics!'.format(rxn, self.index))
            elif isinstance(rxn.kinetics, KineticsData):
                if len(rxn.reactants) == 1:
                    kunits = 's^-1'
                elif len(rxn.reactants) == 2:
                    kunits = 'm^3/(mol*s)'
                elif len(rxn.reactants) == 3:
                    kunits = 'm^6/(mol^2*s)'
                else:
                    kunits = ''
                rxn.kinetics = Arrhenius().fitToData(Tlist=rxn.kinetics.Tdata.values, klist=rxn.kinetics.kdata.values, kunits=kunits)
            elif not isinstance(rxn.kinetics, Arrhenius):
                raise Exception('Path reaction "{0}" in PDepNetwork #{1:d} has invalid kinetics type "{2!s}".'.format(rxn, self.index, rxn.kinetics.__class__))
            rxn.transitionState = rmgpy.species.TransitionState(
                E0=((sum([spec.E0.value for spec in rxn.reactants]) + rxn.kinetics.Ea.value)/1000.,"kJ/mol"),
            )

        # Set collision model
        bathGas = [spec for spec in reactionModel.core.species if not spec.reactive]
        self.bathGas = {}
        for spec in bathGas:
            self.bathGas[spec] = 1.0 / len(bathGas)
            spec.collisionModel = SingleExponentialDown(alpha0=4.86 * 4184)

        # Save input file
        measure.saveInput(os.path.join(outputDirectory, 'pdep', 'network{0:d}_{1:d}.py'.format(self.index, len(self.isomers))))
        
        self.printSummary(level=logging.INFO)

        # Automatically choose a suitable set of energy grains if they were not
        # explicitly specified in the input file
        Elist = self.autoGenerateEnergyGrains(Tmax=Tmax, grainSize=grainSize, Ngrains=numGrains)
        
        # Calculate the rate coefficients
        K, p0 = self.calculateRateCoefficients(Tlist, Plist, Elist, method)

        # Generate PDepReaction objects
        configurations = []
        configurations.extend([[isom] for isom in self.isomers])
        configurations.extend([reactants for reactants in self.reactants])
        configurations.extend([products for products in self.products])
        j = configurations.index(self.source)

        for i in range(K.shape[2]):
            if i != j:
                # Find the path reaction
                netReaction = None
                for r in self.netReactions:
                    if r.hasTemplate(configurations[j], configurations[i]):
                        netReaction = r
                # If net reaction does not already exist, make a new one
                if netReaction is None:
                    netReaction = PDepReaction(
                        reactants=configurations[j],
                        products=configurations[i],
                        network=self,
                        kinetics=None
                    )
                    netReaction = reactionModel.makeNewPDepReaction(netReaction)
                    self.netReactions.append(netReaction)

                    # Place the net reaction in the core or edge if necessary
                    # Note that leak reactions are not placed in the edge
                    if all([s in reactionModel.core.species for s in netReaction.reactants]) and all([s in reactionModel.core.species for s in netReaction.products]):
                        reactionModel.addReactionToCore(netReaction)
                    else:
                        reactionModel.addReactionToEdge(netReaction)

                # Set/update the net reaction kinetics using interpolation model
                netReaction.kinetics = fitInterpolationModel(netReaction, Tlist, Plist, K[:,:,i,j], model, Tmin, Tmax, Pmin, Pmax)

        # We're done processing this network, so mark it as valid
        self.valid = True
