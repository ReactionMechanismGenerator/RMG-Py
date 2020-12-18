#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
Contains classes for providing pressure-dependent kinetics estimation
functionality to RMG.
"""

import logging
import os.path

import mpmath as mp
import numpy as np
import scipy.optimize as opt

import rmgpy.pdep.network
import rmgpy.reaction
from rmgpy.constants import R
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.exceptions import PressureDependenceError, NetworkError
from rmgpy.pdep import Configuration
from rmgpy.rmg.react import react_species
from rmgpy.statmech import Conformer


################################################################################

class PDepReaction(rmgpy.reaction.Reaction):

    def __init__(self,
                 index=-1,
                 label='',
                 reactants=None,
                 products=None,
                 specific_collider=None,
                 network=None,
                 kinetics=None,
                 network_kinetics=None,
                 reversible=True,
                 transition_state=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None
                 ):
        rmgpy.reaction.Reaction.__init__(self,
                                         index=index,
                                         label=label,
                                         reactants=reactants,
                                         products=products,
                                         specific_collider=specific_collider,
                                         kinetics=kinetics,
                                         network_kinetics=network_kinetics,
                                         reversible=reversible,
                                         transition_state=transition_state,
                                         duplicate=duplicate,
                                         degeneracy=degeneracy,
                                         pairs=pairs
                                         )
        self.network = network

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (PDepReaction, (self.index,
                               self.label,
                               self.reactants,
                               self.products,
                               self.specific_collider,
                               self.network,
                               self.kinetics,
                               self.reversible,
                               self.transition_state,
                               self.duplicate,
                               self.degeneracy,
                               self.pairs
                               ))

    def get_source(self):
        """
        Get the source of this PDepReaction
        """
        return str(self.network)


################################################################################

class PDepNetwork(rmgpy.pdep.network.Network):
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
        rmgpy.pdep.network.Network.__init__(self, label="PDepNetwork #{0}".format(index))
        self.index = index
        self.source = source
        self.energy_correction = None
        self.explored = []

    def __str__(self):
        return "PDepNetwork #{0}".format(self.index)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (PDepNetwork, (self.index, self.source), self.__dict__)

    def __setstate__(self, dict):
        self.__dict__.update(dict)

    def cleanup(self):
        """
        Delete intermedate arrays used to compute k(T,P) values.
        """
        for isomer in self.isomers:
            isomer.cleanup()
        for reactant in self.reactants:
            reactant.cleanup()
        for product in self.products:
            product.cleanup()

        self.e_list = None
        self.j_list = None
        self.dens_states = None
        self.coll_freq = None
        self.Mcoll = None
        self.Kij = None
        self.Fim = None
        self.Gnj = None
        self.E0 = None
        self.n_grains = 0
        self.n_j = 0

        self.K = None
        self.p0 = None

    def get_leak_coefficient(self, T, P):
        """
        Return the pressure-dependent rate coefficient :math:`k(T,P)` describing
        the total rate of "leak" from this network. This is defined as the sum
        of the :math:`k(T,P)` values for all net reactions to nonexplored
        unimolecular isomers.
        """
        k = 0.0
        if len(self.net_reactions) == 0 and len(self.path_reactions) == 1:
            # The network is of the form A + B -> C* (with C* nonincluded)
            # For this special case we use the high-pressure limit k(T) to
            # ensure that we're estimating the total leak flux
            rxn = self.path_reactions[0]
            if rxn.kinetics is None:
                if rxn.reverse.kinetics is not None:
                    rxn = rxn.reverse
                else:
                    raise PressureDependenceError('Path reaction {0} with no high-pressure-limit kinetics encountered '
                                                  'in PDepNetwork #{1:d} while evaluating leak flux.'.format(rxn, self.index))
            if rxn.products is self.source:
                k = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
            else:
                k = rxn.get_rate_coefficient(T, P)
        else:
            # The network has at least one included isomer, so we can calculate
            # the leak flux normally
            for rxn in self.net_reactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    k += rxn.get_rate_coefficient(T, P)
        return k

    def get_maximum_leak_species(self, T, P):
        """
        Get the unexplored (unimolecular) isomer with the maximum leak flux.
        Note that the leak rate coefficients vary with temperature and
        pressure, so you must provide these in order to get a meaningful result.
        """
        # Choose species with maximum leak flux
        max_k = 0.0
        max_species = None
        if len(self.net_reactions) == 0 and len(self.path_reactions) == 1:
            max_k = self.get_leak_coefficient(T, P)
            rxn = self.path_reactions[0]
            if rxn.products == self.source:
                assert len(rxn.reactants) == 1
                max_species = rxn.reactants[0]
            else:
                assert len(rxn.products) == 1
                max_species = rxn.products[0]
        else:
            for rxn in self.net_reactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    k = rxn.get_rate_coefficient(T, P)
                    if max_species is None or k > max_k:
                        max_species = rxn.products[0]
                        max_k = k

        # Make sure we've identified a species
        if max_species is None:
            raise NetworkError('No unimolecular isomers left to explore!')
        # Return the species
        return max_species

    def get_leak_branching_ratios(self, T, P):
        """
        Return a dict with the unexplored isomers in the partial network as the
        keys and the fraction of the total leak coefficient as the values.
        """
        ratios = {}
        if len(self.net_reactions) == 0 and len(self.path_reactions) == 1:
            rxn = self.path_reactions[0]
            assert rxn.reactants == self.source or rxn.products == self.source
            if rxn.products == self.source:
                assert len(rxn.reactants) == 1
                ratios[rxn.reactants[0]] = 1.0
            else:
                assert len(rxn.products) == 1
                ratios[rxn.products[0]] = 1.0
        else:
            for rxn in self.net_reactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    ratios[rxn.products[0]] = rxn.get_rate_coefficient(T, P)

        kleak = sum(ratios.values())
        for spec in ratios:
            ratios[spec] /= kleak

        return ratios

    def explore_isomer(self, isomer):
        """
        Explore a previously-unexplored unimolecular `isomer` in this partial
        network using the provided core-edge reaction model `reaction_model`,
        returning the new reactions and new species.
        """

        if isomer in self.explored:
            logging.warning('Already explored isomer {0} in pressure-dependent network #{1:d}'.format(isomer,
                                                                                                      self.index))
            return []

        assert isomer not in self.source, "Attempted to explore isomer {0}, but that is the source configuration for this network.".format(isomer)

        for product in self.products:
            if product.species == [isomer]:
                break
        else:
            raise Exception('Attempted to explore isomer {0}, but that species not found in product channels.'.format(isomer))

        logging.info('Exploring isomer {0} in pressure-dependent network #{1:d}'.format(isomer, self.index))

        for mol in isomer.molecule:
            mol.update()

        self.explored.append(isomer)
        self.isomers.append(product)
        self.products.remove(product)
        # Find reactions involving the found species as unimolecular
        # reactant or product (e.g. A <---> products)

        # Don't find reactions involving the new species as bimolecular
        # reactants or products with itself (e.g. A + A <---> products)
        # Don't find reactions involving the new species as bimolecular
        # reactants or products with other core species (e.g. A + B <---> products)

        new_reactions = react_species((isomer,))

        return new_reactions

    def add_path_reaction(self, newReaction):
        """
        Add a path reaction to the network. If the path reaction already exists,
        no action is taken.
        """
        # Add this reaction to that network if not already present
        found = False
        for rxn in self.path_reactions:
            if newReaction.reactants == rxn.reactants and newReaction.products == rxn.products:
                found = True
                break
            elif newReaction.products == rxn.reactants and newReaction.reactants == rxn.products:
                found = True
                break
        if not found:
            self.path_reactions.append(newReaction)
            self.invalidate()

    def get_energy_filtered_reactions(self, T, tol):
        """
        Returns a list of products and isomers that are greater in Free Energy
        than a*R*T + Gfsource(T)
        """
        dE = tol * R * T
        for conf in self.isomers + self.products + self.reactants:
            if len(conf.species) == len(self.source):
                if len(self.source) == 1:
                    if self.source[0].is_isomorphic(conf.species[0]):
                        E0source = conf.E0
                        break
                elif len(self.source) == 2:
                    boo00 = self.source[0].is_isomorphic(conf.species[0])
                    boo01 = self.source[0].is_isomorphic(conf.species[1])
                    if boo00 or boo01:  # if we found source[0]
                        boo10 = self.source[1].is_isomorphic(conf.species[0])
                        boo11 = self.source[1].is_isomorphic(conf.species[1])
                        if (boo00 and boo11) or (boo01 and boo10):
                            E0source = conf.E0
                            break
        else:
            raise ValueError('No isomer, product or reactant channel is isomorphic to the source')

        filtered_rxns = []
        for rxn in self.path_reactions:
            E0 = rxn.transition_state.conformer.E0.value_si
            if E0 - E0source > dE:
                filtered_rxns.append(rxn)

        return filtered_rxns

    def get_rate_filtered_products(self, T, P, tol):
        """
        determines the set of path_reactions that have fluxes less than
        tol at steady state where all A => B + C reactions are irreversible
        and there is a constant flux from/to the source configuration of 1.0
        """
        c = self.solve_ss_network(T, P)
        isomer_spcs = [iso.species[0] for iso in self.isomers]
        filtered_prod = []
        if c is not None:
            for rxn in self.net_reactions:
                val = 0.0
                val2 = 0.0
                if rxn.reactants[0] in isomer_spcs:
                    ind = isomer_spcs.index(rxn.reactants[0])
                    kf = rxn.get_rate_coefficient(T, P)
                    val = kf * c[ind]
                if rxn.products[0] in isomer_spcs:
                    ind2 = isomer_spcs.index(rxn.products[0])
                    kr = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
                    val2 = kr * c[ind2]

                if max(val, val2) < tol:
                    filtered_prod.append(rxn.products)

            return filtered_prod

        else:
            logging.warning("Falling back flux reduction from Steady State analysis to rate coefficient analysis")
            ks = np.array([rxn.get_rate_coefficient(T, P) for rxn in self.net_reactions])
            frs = ks / ks.sum()
            inds = [i for i in range(len(frs)) if frs[i] < tol]
            filtered_prod = [self.net_reactions[i].products for i in inds]
            return filtered_prod

    def solve_ss_network(self, T, P):
        """
        calculates the steady state concentrations if all A => B + C
        reactions are irreversible and the flux from/to the source
        configuration is 1.0
        """
        A = np.zeros((len(self.isomers), len(self.isomers)))
        b = np.zeros(len(self.isomers))
        bimolecular = len(self.source) > 1

        isomer_spcs = [iso.species[0] for iso in self.isomers]

        for rxn in self.net_reactions:
            if rxn.reactants[0] in isomer_spcs:
                ind = isomer_spcs.index(rxn.reactants[0])
                kf = rxn.get_rate_coefficient(T, P)
                A[ind, ind] -= kf
            else:
                ind = None
            if rxn.products[0] in isomer_spcs:
                ind2 = isomer_spcs.index(rxn.products[0])
                kr = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
                A[ind2, ind2] -= kr
            else:
                ind2 = None

            if ind is not None and ind2 is not None:
                A[ind, ind2] += kr
                A[ind2, ind] += kf

            if bimolecular:
                if rxn.reactants[0] == self.source:
                    kf = rxn.get_rate_coefficient(T, P)
                    b[ind2] += kf
                elif rxn.products[0] == self.source:
                    kr = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
                    b[ind] += kr

        if not bimolecular:
            ind = isomer_spcs.index(self.source[0])
            b[ind] = -1.0  # flux at source
        else:
            b = -b / b.sum()  # 1.0 flux from source

        if len(b) == 1:
            return np.array([b[0] / A[0, 0]])

        con = np.linalg.cond(A)

        if np.log10(con) < 15:
            c = np.linalg.solve(A, b)
        else:
            logging.warning("Matrix Ill-conditioned, attempting to use Arbitrary Precision Arithmetic")
            mp.dps = 30 + int(np.log10(con))
            Amp = mp.matrix(A.tolist())
            bmp = mp.matrix(b.tolist())

            try:
                c = mp.qr_solve(Amp, bmp)

                c = np.array(list(c[0]))

                if any(c <= 0.0):
                    c, rnorm = opt.nnls(A, b)

                c = c.astype(np.float64)
            except:  # fall back to raw flux analysis rather than solve steady state problem
                return None
        
        if np.isnan(c).any():
            return None
        
        return c

    def remove_disconnected_reactions(self):
        """
        gets rid of reactions/isomers/products not connected to the source by a reaction sequence
        """
        kept_reactions = []
        kept_products = [self.source]
        incomplete = True
        while incomplete:
            s = len(kept_reactions)
            for rxn in self.path_reactions:
                if not rxn in kept_reactions:
                    if rxn.reactants in kept_products:
                        kept_products.append(rxn.products)
                        kept_reactions.append(rxn)
                    elif rxn.products in kept_products:
                        kept_products.append(rxn.reactants)
                        kept_reactions.append(rxn)

            incomplete = s != len(kept_reactions)

        logging.info('Removing disconnected items')
        for rxn in self.path_reactions:
            if rxn not in kept_reactions:
                logging.info('Removing rxn: {}'.format(rxn))
                self.path_reactions.remove(rxn)

        nrxns = []
        for nrxn in self.net_reactions:
            if nrxn.products not in kept_products or nrxn.reactants not in kept_products:
                logging.info('Removing net rxn: {}'.format(nrxn))
            else:
                logging.info('Keeping net rxn: {}'.format(nrxn))
                nrxns.append(nrxn)
        self.net_reactions = nrxns

        prods = []
        for prod in self.products:
            if prod.species not in kept_products:
                logging.info('Removing product: {}'.format(prod))
            else:
                logging.info("Keeping product: {}".format(prod))
                prods.append(prod)

        self.products = prods

        rcts = []
        for rct in self.reactants:
            if rct.species not in kept_products:
                logging.info('Removing product: {}'.format(rct))
            else:
                logging.info("Keeping product: {}".format(rct))
                rcts.append(rct)
        self.reactants = rcts

        isos = []
        for iso in self.isomers:
            if iso.species not in kept_products:
                logging.info('Removing isomer: {}'.format(iso))
            else:
                logging.info("Keeping isomer: {}".format(iso))
                isos.append(iso)

        self.isomers = isos
        self.explored = [iso.species[0] for iso in isos]

        self.n_isom = len(self.isomers)
        self.n_reac = len(self.reactants)
        self.n_prod = len(self.products)

    def remove_reactions(self, reaction_model, rxns=None, prods=None):
        """
        removes a list of reactions from the network and all reactions/products
        left disconnected by removing those reactions
        """
        if rxns:
            for rxn in rxns:
                self.path_reactions.remove(rxn)

        if prods:
            isomers = [x.species[0] for x in self.isomers]

            for prod in prods:
                prod = [x for x in prod]
                if prod[0] in isomers:  # skip isomers
                    continue
                for rxn in self.path_reactions:
                    if rxn.products == prod or rxn.reactants == prod:
                        self.path_reactions.remove(rxn)

            prodspc = [x[0] for x in prods]
            for prod in prods:
                prod = [x for x in prod]
                if prod[0] in isomers:  # deal with isomers
                    for rxn in self.path_reactions:
                        if rxn.reactants == prod and rxn.products[0] not in isomers and rxn.products[0] not in prodspc:
                            break
                        if rxn.products == prod and rxn.reactants[0] not in isomers and rxn.reactants not in prodspc:
                            break
                    else:
                        for rxn in self.path_reactions:
                            if rxn.reactants == prod or rxn.products == prod:
                                self.path_reactions.remove(rxn)

        self.remove_disconnected_reactions()

        self.cleanup()

        self.invalidate()

        assert self.path_reactions != [], 'Reduction process removed all reactions, cannot update network with no reactions'

        reaction_model.update_unimolecular_reaction_networks()

        if reaction_model.pressure_dependence.output_file:
            path = os.path.join(reaction_model.pressure_dependence.output_file, 'pdep')

            for name in os.listdir(path):  # remove the old reduced file
                if name.endswith('reduced.py'):
                    os.remove(os.path.join(path, name))

            for name in os.listdir(path):  # find the new file and name it network_reduced.py
                if not "full" in name:
                    os.rename(os.path.join(path, name), os.path.join(path, name.split("_")[0]+"_reduced.py"))

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
        products_to_remove = []
        for products in self.products:
            if len(products.species) == 1 and products.species[0] in self.isomers:
                products_to_remove.append(products)
        for products in products_to_remove:
            self.products.remove(products)

        # Merge path reactions
        for reaction in other.path_reactions:
            found = False
            for rxn in self.path_reactions:
                if reaction.reactants == rxn.reactants and reaction.products == rxn.products:
                    # NB the isEquivalent() method that used to be on the previous line also checked reverse direction.
                    # I am not sure which is appropriate 
                    found = True
                    break
            if not found:
                self.path_reactions.append(reaction)

        # Also merge net reactions (so that when we update the network in the
        # future, we update the existing net reactions rather than making new ones)
        # Q: What to do when a net reaction exists in both networks being merged?
        for reaction in other.net_reactions:
            found = False
            for rxn in self.net_reactions:
                if reaction.reactants == rxn.reactants and reaction.products == rxn.products:
                    # NB the isEquivalent() method that used to be on the previous line also checked reverse direction.
                    # I am not sure which is appropriate 
                    found = True
                    break
            if not found:
                self.net_reactions.append(reaction)

        # Mark this network as invalid
        self.valid = False

    def update_configurations(self, reaction_model):
        """
        Sort the reactants and products of each of the network's path reactions
        into isomers, reactant channels, and product channels. You must pass 
        the current `reaction_model` because some decisions on sorting are made
        based on which species are in the model core. 
        """
        reactants = []
        products = []

        # All explored species are isomers
        isomers = self.explored[:]

        # The source configuration is an isomer (if unimolecular) or a reactant channel (if bimolecular)
        if len(self.source) == 1:
            # The source is a unimolecular isomer
            if self.source[0] not in isomers: isomers.insert(0, self.source[0])
        else:
            # The source is a bimolecular reactant channel
            self.source.sort()
            reactants.append(self.source)

        # Iterate over path reactions and make sure each set of reactants and products is classified
        for rxn in self.path_reactions:
            # Sort bimolecular configurations so that we always encounter them in the
            # same order
            # The actual order doesn't matter, as long as it is consistent
            rxn.reactants.sort()
            rxn.products.sort()
            # Reactants of the path reaction
            if len(rxn.reactants) == 1 and rxn.reactants[0] not in isomers and rxn.reactants not in products:
                # We've encountered a unimolecular reactant that is not classified
                # These are always product channels (since they would be in source or explored otherwise)
                products.append(rxn.reactants)
            elif len(rxn.reactants) > 1 and rxn.reactants not in reactants and rxn.reactants not in products:
                # We've encountered bimolecular reactants that are not classified
                if all([reactant in reaction_model.core.species for reactant in rxn.reactants]):
                    # Both reactants are in the core, so treat as reactant channel
                    reactants.append(rxn.reactants)
                else:
                    # One or more reactants is an edge species, so treat as product channel
                    products.append(rxn.reactants)
            # Products of the path reaction
            if len(rxn.products) == 1 and rxn.products[0] not in isomers and rxn.products not in products:
                # We've encountered a unimolecular product that is not classified
                # These are always product channels (since they would be in source or explored otherwise)
                products.append(rxn.products)
            elif len(rxn.products) > 1 and rxn.products not in reactants and rxn.products not in products:
                # We've encountered bimolecular products that are not classified
                if all([product in reaction_model.core.species for product in rxn.products]):
                    # Both products are in the core, so treat as reactant channel
                    reactants.append(rxn.products)
                else:
                    # One or more reactants is an edge species, so treat as product channel
                    products.append(rxn.products)

        # Clear existing configurations
        self.isomers = []
        self.reactants = []
        self.products = []

        # Make a configuration object for each
        for isomer in isomers:
            self.isomers.append(Configuration(isomer))
        for reactant in reactants:
            self.reactants.append(Configuration(*reactant))
        for product in products:
            self.products.append(Configuration(*product))
        if self.energy_correction:
            for spec in self.reactants + self.products + self.isomers:
                spec.energy_correction = self.energy_correction

    def update(self, reaction_model, pdep_settings):
        """
        Regenerate the :math:`k(T,P)` values for this partial network if the
        network is marked as invalid.
        """
        from rmgpy.kinetics import Arrhenius, KineticsData, MultiArrhenius

        # Get the parameters for the pressure dependence calculation
        job = pdep_settings
        job.network = self
        output_directory = pdep_settings.output_file

        Tmin = job.Tmin.value_si
        Tmax = job.Tmax.value_si
        Pmin = job.Pmin.value_si
        Pmax = job.Pmax.value_si
        Tlist = job.Tlist.value_si
        Plist = job.Plist.value_si
        maximum_grain_size = job.maximum_grain_size.value_si if job.maximum_grain_size is not None else 0.0
        minimum_grain_count = job.minimum_grain_count
        method = job.method
        interpolation_model = job.interpolation_model
        active_j_rotor = job.active_j_rotor
        active_k_rotor = job.active_k_rotor
        rmgmode = job.rmgmode

        # Figure out which configurations are isomers, reactant channels, and product channels
        self.update_configurations(reaction_model)

        # Make sure we have high-P kinetics for all path reactions
        for rxn in self.path_reactions:
            if rxn.kinetics is None and rxn.reverse.kinetics is None:
                raise PressureDependenceError('Path reaction {0} with no high-pressure-limit kinetics encountered in '
                                              'PDepNetwork #{1:d}.'.format(rxn, self.index))
            elif rxn.kinetics is not None and rxn.kinetics.is_pressure_dependent() and rxn.network_kinetics is None:
                raise PressureDependenceError('Pressure-dependent kinetics encountered for path reaction {0} in '
                                              'PDepNetwork #{1:d}.'.format(rxn, self.index))

        # Do nothing if the network is already valid
        if self.valid:
            return
        # Do nothing if there are no explored wells
        if len(self.explored) == 0 and len(self.source) > 1:
            return
        # Log the network being updated
        logging.info("Updating {0!s}".format(self))

        E0 = []
        # Generate states data for unimolecular isomers and reactants if necessary
        for isomer in self.isomers:
            spec = isomer.species[0]
            if not spec.has_statmech():
                spec.generate_statmech()
            E0.append(spec.conformer.E0.value_si)
        for reactants in self.reactants:
            for spec in reactants.species:
                if not spec.has_statmech():
                    spec.generate_statmech()
                E0.append(spec.conformer.E0.value_si)
        # Also generate states data for any path reaction reactants, so we can
        # always apply the ILT method in the direction the kinetics are known
        for reaction in self.path_reactions:
            for spec in reaction.reactants:
                if not spec.has_statmech():
                    spec.generate_statmech()
                E0.append(spec.conformer.E0.value_si)
        # While we don't need the frequencies for product channels, we do need
        # the E0, so create a conformer object with the E0 for the product
        # channel species if necessary
        for products in self.products:
            for spec in products.species:
                if spec.conformer is None:
                    spec.conformer = Conformer(E0=spec.get_thermo_data().E0)
                E0.append(spec.conformer.E0.value_si)

        # Use the average E0 as the reference energy (`energy_correction`) for the network
        # The `energy_correction` will be added to the free energies and enthalpies for each
        # configuration in the network.
        energy_correction = -np.array(E0).mean()
        for spec in self.reactants + self.products + self.isomers:
            spec.energy_correction = energy_correction
        self.energy_correction = energy_correction

        # Determine transition state energies on potential energy surface
        # In the absence of any better information, we simply set it to
        # be the reactant ground-state energy + the activation energy
        # Note that we need Arrhenius kinetics in order to do this
        for rxn in self.path_reactions:
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
                rxn.kinetics = Arrhenius().fit_to_data(Tlist=rxn.kinetics.Tdata.value_si,
                                                       klist=rxn.kinetics.kdata.value_si, kunits=kunits)
            elif isinstance(rxn.kinetics, MultiArrhenius):
                logging.info('Converting multiple kinetics to a single Arrhenius expression for reaction {rxn}'.format(
                    rxn=rxn))
                rxn.kinetics = rxn.kinetics.to_arrhenius(Tmin=Tmin, Tmax=Tmax)
            elif not isinstance(rxn.kinetics, Arrhenius) and rxn.network_kinetics is None:
                raise Exception('Path reaction "{0}" in PDepNetwork #{1:d} has invalid kinetics '
                                'type "{2!s}".'.format(rxn, self.index, rxn.kinetics.__class__))
            rxn.fix_barrier_height(force_positive=True)
            if rxn.network_kinetics is None:
                E0 = sum([spec.conformer.E0.value_si for spec in rxn.reactants]) + rxn.kinetics.Ea.value_si + energy_correction
            else:
                E0 = sum([spec.conformer.E0.value_si for spec in rxn.reactants]) + rxn.network_kinetics.Ea.value_si + energy_correction
            rxn.transition_state = rmgpy.species.TransitionState(conformer=Conformer(E0=(E0 * 0.001, "kJ/mol")))

        # Set collision model
        bath_gas = [spec for spec in reaction_model.core.species if not spec.reactive]
        assert len(bath_gas) > 0, 'No unreactive species to identify as bath gas'

        self.bath_gas = {}
        for spec in bath_gas:
            # is this really the only/best way to weight them?
            self.bath_gas[spec] = 1.0 / len(bath_gas)

        # Save input file
        if not self.label:
            self.label = str(self.index)

        if output_directory:
            job.save_input_file(
                os.path.join(output_directory, 'pdep', 'network{0:d}_{1:d}.py'.format(self.index, len(self.isomers))))

        # Calculate the rate coefficients
        self.initialize(Tmin, Tmax, Pmin, Pmax, maximum_grain_size, minimum_grain_count, active_j_rotor, active_k_rotor,
                        rmgmode)
        K = self.calculate_rate_coefficients(Tlist, Plist, method)

        # Generate PDepReaction objects
        configurations = []
        configurations.extend([isom.species[:] for isom in self.isomers])
        configurations.extend([reactant.species[:] for reactant in self.reactants])
        configurations.extend([product.species[:] for product in self.products])
        j = configurations.index(self.source)

        for i in range(K.shape[2]):
            if i != j:
                # Find the path reaction
                net_reaction = None
                for r in self.net_reactions:
                    if r.has_template(configurations[j], configurations[i]):
                        net_reaction = r
                # If net reaction does not already exist, make a new one
                if net_reaction is None:
                    net_reaction = PDepReaction(
                        reactants=configurations[j],
                        products=configurations[i],
                        network=self,
                        kinetics=None
                    )
                    net_reaction = reaction_model.make_new_pdep_reaction(net_reaction)
                    self.net_reactions.append(net_reaction)

                    # Place the net reaction in the core or edge if necessary
                    # Note that leak reactions are not placed in the edge
                    if all([s in reaction_model.core.species for s in net_reaction.reactants]) \
                            and all([s in reaction_model.core.species for s in net_reaction.products]):
                        # Check whether netReaction already exists in the core as a LibraryReaction
                        for rxn in reaction_model.core.reactions:
                            if isinstance(rxn, LibraryReaction) \
                                    and rxn.is_isomorphic(net_reaction, either_direction=True) \
                                    and not rxn.allow_pdep_route and not rxn.elementary_high_p:
                                logging.info('Network reaction {0} matched an existing core reaction {1}'
                                             ' from the {2} library, and was not added to the model'.format(
                                    str(net_reaction), str(rxn), rxn.library))
                                break
                        else:
                            reaction_model.add_reaction_to_core(net_reaction)
                    else:
                        # Check whether netReaction already exists in the edge as a LibraryReaction
                        for rxn in reaction_model.edge.reactions:
                            if isinstance(rxn, LibraryReaction) \
                                    and rxn.is_isomorphic(net_reaction, either_direction=True) \
                                    and not rxn.allow_pdep_route and not rxn.elementary_high_p:
                                logging.info('Network reaction {0} matched an existing edge reaction {1}'
                                             ' from the {2} library, and was not added to the model'.format(
                                    str(net_reaction), str(rxn), rxn.library))
                                break
                        else:
                            reaction_model.add_reaction_to_edge(net_reaction)

                # Set/update the net reaction kinetics using interpolation model
                kdata = K[:, :, i, j].copy()
                order = len(net_reaction.reactants)
                kdata *= 1e6 ** (order - 1)
                kunits = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]
                net_reaction.kinetics = job.fit_interpolation_model(Tlist, Plist, kdata, kunits)

                # Check: For each net reaction that has a path reaction, make
                # sure the k(T,P) values for the net reaction do not exceed
                # the k(T) values of the path reaction
                # Only check the k(T,P) value at the highest P and lowest T,
                # as this is the one most likely to be in the high-pressure 
                # limit
                t = 0
                p = len(Plist) - 1
                for pathReaction in self.path_reactions:
                    if pathReaction.is_isomerization():
                        # Don't check isomerization reactions, since their
                        # k(T,P) values potentially contain both direct and
                        # well-skipping contributions, and therefore could be
                        # significantly larger than the direct k(T) value
                        # (This can also happen for association/dissociation
                        # reactions, but the effect is generally not too large)
                        continue
                    if pathReaction.reactants == net_reaction.reactants and pathReaction.products == net_reaction.products:
                        if pathReaction.network_kinetics is not None:
                            kinf = pathReaction.network_kinetics.get_rate_coefficient(Tlist[t])
                        else:
                            kinf = pathReaction.kinetics.get_rate_coefficient(Tlist[t])
                        if K[t, p, i, j] > 2 * kinf:  # To allow for a small discretization error
                            logging.warning('k(T,P) for net reaction {0} exceeds high-P k(T) by {1:g} at {2:g} K, '
                                            '{3:g} bar'.format(net_reaction, K[t, p, i, j] / kinf, Tlist[t], Plist[p] / 1e5))
                            logging.info('    k(T,P) = {0:9.2e}    k(T) = {1:9.2e}'.format(K[t, p, i, j], kinf))
                        break
                    elif pathReaction.products == net_reaction.reactants and pathReaction.reactants == net_reaction.products:
                        if pathReaction.network_kinetics is not None:
                            kinf = pathReaction.network_kinetics.get_rate_coefficient(
                                Tlist[t]) / pathReaction.get_equilibrium_constant(Tlist[t])
                        else:
                            kinf = pathReaction.kinetics.get_rate_coefficient(
                                Tlist[t]) / pathReaction.get_equilibrium_constant(Tlist[t])
                        if K[t, p, i, j] > 2 * kinf:  # To allow for a small discretization error
                            logging.warning('k(T,P) for net reaction {0} exceeds high-P k(T) by {1:g} at {2:g} K, '
                                            '{3:g} bar'.format(net_reaction, K[t, p, i, j] / kinf, Tlist[t], Plist[p] / 1e5))
                            logging.info('    k(T,P) = {0:9.2e}    k(T) = {1:9.2e}'.format(K[t, p, i, j], kinf))
                        break

        self.log_summary(level=logging.INFO)

        # Delete intermediate arrays to conserve memory
        self.cleanup()

        # We're done processing this network, so mark it as valid
        self.valid = True
