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
The Arkane Explorer module
"""

import logging
import os
import shutil
from copy import deepcopy

import numpy as np

import rmgpy
from rmgpy.data.rmg import get_db
from rmgpy.exceptions import InputError
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.statmech.conformer import Conformer

################################################################################


class ExplorerJob(object):
    """
    A representation of an Arkane explorer job. This job is used to explore a potential energy surface (PES).
    """

    def __init__(self, source, pdepjob, explore_tol, energy_tol=np.inf, flux_tol=0.0,
                 bath_gas=None, maximum_radical_electrons=np.inf):
        self.source = source
        self.explore_tol = explore_tol
        self.energy_tol = energy_tol
        self.flux_tol = flux_tol
        self.maximum_radical_electrons = maximum_radical_electrons
        self.job_rxns = None
        self.networks = None

        self.pdepjob = pdepjob

        if not hasattr(self.pdepjob, 'output_file'):
            self.pdepjob.output_file = None

        if bath_gas:
            self.bath_gas = bath_gas
        elif self.pdepjob.network and self.pdepjob.network.bath_gas:
            self.bath_gas = self.pdepjob.network.bath_gas
        else:
            raise InputError('bathGas not specified in explorer block')

    def copy(self):
        """
        Return a copy of the explorer job.
        """
        return ExplorerJob(
            source=deepcopy(self.source),
            pdepjob=self.pdepjob,
            explore_tol=self.explore_tol,
            energy_tol=self.energy_tol,
            flux_tol=self.flux_tol
        )

    def execute(self, output_file, plot, file_format='pdf', print_summary=True, species_list=None,
                thermo_library=None, kinetics_library=None):
        """Execute an ExplorerJob"""
        logging.info('Exploring network...')

        rmg = RMG()

        rmg.species_constraints = {'allowed': ['input species', 'seed mechanisms', 'reaction libraries'],
                                   'maximumRadicalElectrons': self.maximum_radical_electrons,
                                   'explicitlyAllowedMolecules': []}

        rmgpy.rmg.input.rmg = rmg

        reaction_model = CoreEdgeReactionModel()

        reaction_model.pressure_dependence = self.pdepjob

        reaction_model.pressure_dependence.rmgmode = True

        if output_file:
            reaction_model.pressure_dependence.output_file = os.path.dirname(output_file)

        kinetics_database = get_db('kinetics')
        thermo_database = get_db('thermo')

        thermo_database.libraries['thermojobs'] = thermo_library
        thermo_database.library_order.insert(0, 'thermojobs')

        kinetics_database.libraries['kineticsjobs'] = kinetics_library
        kinetics_database.library_order.insert(0, ('kineticsjobs', 'Reaction Library'))

        self.job_rxns = [rxn for rxn in reaction_model.core.reactions]

        if output_file is not None:
            if not os.path.exists(os.path.join(reaction_model.pressure_dependence.output_file, 'pdep')):
                os.mkdir(os.path.join(reaction_model.pressure_dependence.output_file, 'pdep'))
            else:
                shutil.rmtree(os.path.join(reaction_model.pressure_dependence.output_file, 'pdep'))
                os.mkdir(os.path.join(reaction_model.pressure_dependence.output_file, 'pdep'))

        # get the molecular formula for the network
        mmol = None
        for spc in self.source:
            if mmol:
                mmol = mmol.merge(spc.molecule[0])
            else:
                mmol = spc.molecule[0].copy(deep=True)

        form = mmol.get_formula()

        for spec in list(self.bath_gas.keys()) + self.source:
            nspec, is_new = reaction_model.make_new_species(spec, reactive=False)
            flags = np.array([s.molecule[0].get_formula() == form for s in reaction_model.core.species])
            reaction_model.enlarge(nspec, react_edge=False, unimolecular_react=flags,
                                   bimolecular_react=np.zeros((len(reaction_model.core.species),
                                                              len(reaction_model.core.species))))

        reaction_model.add_seed_mechanism_to_core('kineticsjobs')

        for lib in kinetics_database.library_order:
            if lib[0] != 'kineticsjobs':
                reaction_model.add_reaction_library_to_edge(lib[0])

        for spc in reaction_model.core.species:
            for i, item in enumerate(self.source):
                if spc.is_isomorphic(item):
                    self.source[i] = spc

        # react initial species
        if len(self.source) == 1:
            flags = np.array([s.molecule[0].get_formula() == form for s in reaction_model.core.species])
            biflags = np.zeros((len(reaction_model.core.species), len(reaction_model.core.species)))
        elif len(self.source) == 2:
            flags = np.array([False for s in reaction_model.core.species])
            biflags = np.array([[False for i in range(len(reaction_model.core.species))]
                                for j in range(len(reaction_model.core.species))])
            biflags[reaction_model.core.species.index(self.source[0]), reaction_model.core.species.index(
                self.source[1])] = True
        else:
            raise ValueError("Reactant channels with greater than 2 reactants not supported")

        reaction_model.enlarge(react_edge=True, unimolecular_react=flags,
                               bimolecular_react=biflags)

        # find the networks we're interested in
        networks = []
        for nwk in reaction_model.network_list:
            if set(nwk.source) == set(self.source):
                self.source = nwk.source
                networks.append(nwk)

        if len(networks) == 0:
            raise ValueError('Did not generate a network with the requested source. This usually means no unimolecular '
                             'reactions were generated for the source. Note that library reactions that are not '
                             'properly flagged as elementary_high_p can replace RMG generated reactions that would '
                             'otherwise be part of networks.')
        for network in networks:
            network.bath_gas = self.bath_gas

        self.networks = networks

        # determine T and P combinations

        if self.pdepjob.Tlist:
            t_list = self.pdepjob.Tlist.value_si
        else:
            t_list = np.linspace(self.pdepjob.Tmin.value_si, self.pdepjob.Tmax.value_si, self.pdepjob.Tcount)

        if self.pdepjob.Plist:
            p_list = self.pdepjob.Plist.value_si
        else:
            p_list = np.linspace(self.pdepjob.Pmin.value_si, self.pdepjob.Pmax.value_si, self.pdepjob.Pcount)

        # generate the network

        forbidden_structures = get_db('forbidden')
        incomplete = True
        checked_species = []

        while incomplete:
            incomplete = False
            for temperature in t_list:
                for pressure in p_list:
                    for network in self.networks:
                        # compute the characteristic rate coefficient by summing all rate coefficients
                        # from the reactant channel
                        for spc in reaction_model.edge.species:
                            if spc in checked_species:
                                continue
                            if forbidden_structures.is_molecule_forbidden(spc.molecule[0]):
                                reaction_model.remove_species_from_edge(reaction_model.reaction_systems, spc)
                                reaction_model.remove_empty_pdep_networks()
                            else:
                                checked_species.append(spc)

                        kchar = 0.0
                        for rxn in network.net_reactions:  # reaction_model.core.reactions+reaction_model.edge.reactions
                            if (set(rxn.reactants) == set(self.source)
                                    and rxn.products[0].molecule[0].get_formula() == form):
                                kchar += rxn.kinetics.get_rate_coefficient(T=temperature, P=pressure)
                            elif (set(rxn.products) == set(self.source)
                                  and rxn.reactants[0].molecule[0].get_formula() == form):
                                kchar += rxn.generate_reverse_rate_coefficient(
                                    network_kinetics=True).get_rate_coefficient(T=temperature, P=pressure)

                        if network.get_leak_coefficient(T=temperature, P=pressure) > self.explore_tol * kchar:
                            incomplete = True
                            spc = network.get_maximum_leak_species(T=temperature, P=pressure)
                            logging.info('adding new isomer {0} to network'.format(spc))
                            flags = np.array([s.molecule[0].get_formula() == form for s in reaction_model.core.species])
                            reaction_model.enlarge((network, spc), react_edge=False, unimolecular_react=flags,
                                                   bimolecular_react=np.zeros((len(reaction_model.core.species),
                                                                              len(reaction_model.core.species))))

                            flags = np.array([s.molecule[0].get_formula() == form for s in reaction_model.core.species])
                            reaction_model.enlarge(react_edge=True, unimolecular_react=flags,
                                                   bimolecular_react=np.zeros((len(reaction_model.core.species),
                                                                              len(reaction_model.core.species))))
        
        for network in self.networks:
            for rxn in network.path_reactions:
                if rxn.transition_state is None:
                    if rxn.network_kinetics is None:
                        E0 = sum([spec.conformer.E0.value_si for spec in rxn.reactants]) + rxn.kinetics.Ea.value_si
                    else:
                        E0 = sum([spec.conformer.E0.value_si for spec in rxn.reactants]) + rxn.network_kinetics.Ea.value_si
                    rxn.transition_state = rmgpy.species.TransitionState(conformer=Conformer(E0=(E0 * 0.001, "kJ/mol")))
                
        for network in self.networks:
            rm_rxns = []
            for rxn in network.path_reactions:  # remove reactions with forbidden species
                for r in rxn.reactants + rxn.products:
                    if forbidden_structures.is_molecule_forbidden(r.molecule[0]):
                        rm_rxns.append(rxn)

            for rxn in rm_rxns:
                logging.info('Removing forbidden reaction: {0}'.format(rxn))
                network.path_reactions.remove(rxn)

            # clean up output files
            if output_file is not None:
                path = os.path.join(reaction_model.pressure_dependence.output_file, 'pdep')
                for name in os.listdir(path):
                    if name.endswith('.py') and '_' in name:
                        if name.split('_')[-1].split('.')[0] != str(len(network.isomers)):
                            os.remove(os.path.join(path, name))
                        else:
                            os.rename(os.path.join(path, name),
                                      os.path.join(path, 'network_full{}.py'.format(self.networks.index(network))))

        warns = []

        for rxn in self.job_rxns:
            if rxn not in network.path_reactions:
                warns.append('Reaction {0} in the input file was not explored during network expansion and was '
                             'not included in the full network.  This is likely because your explore_tol value is '
                             'too high.'.format(rxn))

        # reduction process
        for network in self.networks:
            if self.energy_tol != np.inf or self.flux_tol != 0.0:

                rxn_set = None
                product_set = None

                for temperature in t_list:
                    if self.energy_tol != np.inf:
                        rxns = network.get_energy_filtered_reactions(temperature, self.energy_tol)
                        if rxn_set is not None:
                            rxn_set &= set(rxns)
                        else:
                            rxn_set = set(rxns)

                    for pressure in p_list:
                        if self.flux_tol != 0.0:
                            products = network.get_rate_filtered_products(temperature, pressure, self.flux_tol)
                            products = [tuple(x) for x in products]
                            if product_set is not None:
                                product_set &= set(products)
                            else:
                                product_set = set(products)

                if rxn_set:
                    logging.info('removing reactions during reduction:')
                    for rxn in rxn_set:
                        logging.info(rxn)
                    rxn_set = list(rxn_set)
                if product_set:
                    logging.info('removing products during reduction:')
                    for prod in product_set:
                        logging.info([x.label for x in prod])
                    product_set = list(product_set)

                network.remove_reactions(reaction_model, rxns=rxn_set, prods=product_set)

                for rxn in self.job_rxns:
                    if rxn not in network.path_reactions:
                        warns.append(
                            'Reaction {0} in the input file was not included in the reduced model.'.format(rxn))

        self.networks = networks
        for p, network in enumerate(self.networks):
            self.pdepjob.network = network

            if len(self.networks) > 1:
                root, file_name = os.path.split(output_file)
                s1, s2 = file_name.split(".")
                ind = str(self.networks.index(network))
                stot = os.path.join(root, s1 + "{}.".format(ind) + s2)
            else:
                stot = output_file

            self.pdepjob.execute(stot, plot, file_format='pdf', print_summary=True)
            if os.path.isfile('network.pdf'):
                os.rename('network.pdf', 'network' + str(p) + '.pdf')

            if warns:
                logging.info('\nOUTPUT WARNINGS:\n')
                for w in warns:
                    logging.warning(w)
