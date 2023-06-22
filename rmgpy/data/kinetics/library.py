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
This module contains functionality for working with kinetics libraries.
"""

import codecs
import logging
import os.path
import re
from collections import OrderedDict

import numpy as np

from rmgpy.data.base import DatabaseError, Database, Entry
from rmgpy.data.kinetics.common import save_entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics import Arrhenius, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, Chebyshev 
from rmgpy.kinetics.surface import StickingCoefficient
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species


################################################################################

class LibraryReaction(Reaction):
    """
    A Reaction object generated from a reaction library. In addition to the
    usual attributes, this class includes `library` and `entry` attributes to
    store the library and the entry in that library that it was created from.
    """

    def __init__(self,
                 index=-1,
                 reactants=None,
                 products=None,
                 specific_collider=None,
                 kinetics=None,
                 network_kinetics=None,
                 reversible=True,
                 transition_state=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 library=None,
                 allow_pdep_route=False,
                 elementary_high_p=False,
                 allow_max_rate_violation=False,
                 entry=None,
                 ):
        Reaction.__init__(self,
                          index=index,
                          reactants=reactants,
                          products=products,
                          specific_collider=specific_collider,
                          kinetics=kinetics,
                          network_kinetics=network_kinetics,
                          reversible=reversible,
                          transition_state=transition_state,
                          duplicate=duplicate,
                          degeneracy=degeneracy,
                          pairs=pairs,
                          allow_pdep_route=allow_pdep_route,
                          elementary_high_p=elementary_high_p,
                          allow_max_rate_violation=allow_max_rate_violation,
                          )
        self.library = library
        self.family = library
        self.entry = entry

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (LibraryReaction, (self.index,
                                  self.reactants,
                                  self.products,
                                  self.specific_collider,
                                  self.kinetics,
                                  self.network_kinetics,
                                  self.reversible,
                                  self.transition_state,
                                  self.duplicate,
                                  self.degeneracy,
                                  self.pairs,
                                  self.library,
                                  self.allow_pdep_route,
                                  self.elementary_high_p,
                                  self.allow_max_rate_violation,
                                  self.entry
                                  ))

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'LibraryReaction('
        if self.index != -1: string += 'index={0:d}, '.format(self.index)
        if self.reactants is not None: string += 'reactants={0!r}, '.format(self.reactants)
        if self.products is not None: string += 'products={0!r}, '.format(self.products)
        if self.specific_collider is not None: string += 'specific_collider={0!r}, '.format(self.specific_collider)
        if self.kinetics is not None: string += 'kinetics={0!r}, '.format(self.kinetics)
        if not self.reversible: string += 'reversible={0}, '.format(self.reversible)
        if self.transition_state is not None: string += 'transition_state={0!r}, '.format(self.transition_state)
        if self.duplicate: string += 'duplicate={0}, '.format(self.duplicate)
        if self.degeneracy != 1: string += 'degeneracy={0:.1f}, '.format(self.degeneracy)
        if self.pairs is not None: string += 'pairs={0}, '.format(self.pairs)
        if self.allow_pdep_route: string += "allow_pdep_route={}, ".format(self.allow_pdep_route)
        if self.elementary_high_p: string += "elementary_high_p={}, ".format(self.elementary_high_p)
        if self.allow_max_rate_violation: string += "allow_max_rate_violation={}, ".format(self.allow_max_rate_violation)
        string = string[:-2] + ')'
        return string

    def get_source(self):
        """
        Return the database that was the source of this reaction. For a
        LibraryReaction this should be a KineticsLibrary object.
        """
        return self.library

    def generate_high_p_limit_kinetics(self):
        """
        If the LibraryReactions represented by `self` has pressure dependent kinetics,
        try extracting the high pressure limit rate from it.
        Used for incorporating library reactions with pressure-dependent kinetics in PDep networks.
        Only reactions flagged as `elementary_high_p=True` should be processed here.
        If the kinetics is a :class:Lindemann or a :class:Troe, simply get the high pressure limit rate.
        If the kinetics is a :class:PDepArrhenius or a :class:Chebyshev, generate a :class:Arrhenius kinetics entry
        that represents the high pressure limit if Pmax >= 90 bar .
        This high pressure limit Arrhenius kinetics is assigned to the reaction network_kinetics attribute.
        If this method successfully generated the high pressure limit kinetics, return ``True``, otherwise ``False``.
        """
        logging.debug("Generating high pressure limit kinetics for {0}...".format(self))
        if not self.is_unimolecular():
            return False
        if isinstance(self.kinetics, Arrhenius):
            return self.elementary_high_p
        if self.network_kinetics is not None:
            return True
        if self.elementary_high_p:
            if isinstance(self.kinetics, (Lindemann, Troe)):
                self.network_kinetics = self.kinetics.arrheniusHigh
                self.network_kinetics.comment = self.kinetics.comment
                self.network_kinetics.comment = "Kinetics taken from the arrheniusHigh attribute of a" \
                    " Troe/Lindemann exprssion. Originally from reaction library {0}".format(self.library)
                return True
            if isinstance(self.kinetics, PDepArrhenius):
                if self.kinetics.pressures.value_si[-1] >= 9000000:  # Pa units
                    if isinstance(self.kinetics.arrhenius[-1], Arrhenius):
                        self.network_kinetics = self.kinetics.arrhenius[-1]
                        return True
                    else:
                        # This is probably MultiArrhenius entries inside a PDepArrhenius kinetics entry. Don't process
                        return False
            if isinstance(self.kinetics, Chebyshev):
                if self.kinetics.Pmax.value_si >= 9000000:  # Pa units
                    if len(self.reactants) == 1:
                        kunits = 's^-1'
                    elif len(self.reactants) == 2:
                        kunits = 'm^3/(mol*s)'
                    elif len(self.reactants) == 3:
                        kunits = 'm^6/(mol^2*s)'
                    else:
                        kunits = ''
                    t_step = (self.kinetics.Tmax.value_si - self.kinetics.Tmin.value_si) / 20
                    t_list = np.arange(int(self.kinetics.Tmin.value_si), int(self.kinetics.Tmax.value_si), int(t_step))
                    if t_list[-1] < int(self.kinetics.Tmax.value_si):
                        t_list = np.insert(t_list, -1, [int(self.kinetics.Tmax.value_si)])
                    k_list = []
                    for t in t_list:
                        k_list.append(self.kinetics.get_rate_coefficient(t, self.kinetics.Pmax.value_si))
                    k_list = np.array(k_list)
                    self.network_kinetics = Arrhenius().fit_to_data(Tlist=t_list, klist=k_list, kunits=kunits)
                    return True
            logging.info("NOT processing reaction {0} in a pressure-dependent reaction network.\n"
                         "Although it is marked with the `elementary_high_p=True` flag,"
                         " it doesn't answer either of the following criteria:\n1. Has a Lindemann or Troe"
                         " kinetics type; 2. Has a PDepArrhenius or Chebyshev kinetics type and has valid"
                         " kinetics at P >= 100 bar.\n".format(self))
        return False

    def get_sticking_coefficient(self, T):
        """
        Helper function to get sticking coefficient
        """
        if isinstance(self.kinetics, StickingCoefficient):
            stick = self.kinetics.get_sticking_coefficient(T)
            return stick
            
        return False


################################################################################

class KineticsLibrary(Database):
    """
    A class for working with an RMG kinetics library.
    """

    def __init__(self, label='', name='', solvent=None, short_desc='', long_desc='', auto_generated=False):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)
        self.auto_generated = auto_generated

    def __str__(self):
        return 'Kinetics Library {0}'.format(self.label)

    def __repr__(self):
        return '<KineticsLibrary "{0}">'.format(self.label)

    def get_library_reactions(self):
        """
        makes library and template reactions as appropriate from the library comments
        and returns at list of all of these LibraryReaction and TemplateReaction objects
        """
        rxns = []
        for entry in self.entries.values():
            if self.auto_generated and entry.long_desc and 'Originally from reaction library: ' in entry.long_desc:
                lib = [line for line in entry.long_desc.split('\n') if 'Originally from reaction library: ' in line]
                lib = lib[0].replace('Originally from reaction library: ', '')
                lib = lib.replace('\n', '')
                # Clean up any leading indents in kinetics comment
                entry.data.comment = '\n'.join([line.strip() for line in entry.data.comment.split('\n')])
                rxn = LibraryReaction(reactants=entry.item.reactants[:], products=entry.item.products[:],
                                      library=lib, specific_collider=entry.item.specific_collider, kinetics=entry.data,
                                      duplicate=entry.item.duplicate, reversible=entry.item.reversible,
                                      allow_pdep_route=entry.item.allow_pdep_route,
                                      elementary_high_p=entry.item.elementary_high_p)
                rxn.family = self.label  # the library the reaction was loaded from (opposed to originally from)
            elif self.auto_generated and entry.long_desc and 'rate rule' in entry.long_desc:  # template reaction
                family = ''
                template = ''
                for line in entry.long_desc.split('\n'):
                    line = line.strip()
                    if 'rate rule' in line:
                        rule_part = line.split('rate rule')[1].strip()
                        template = re.split(r'\[(.*)\]', rule_part)[1]  # Remove outer brackets and trailing text
                    elif 'family:' in line:
                        family = line.split('family:')[1].strip()
                # Clean up any leading indents in kinetics comment
                entry.data.comment = '\n'.join([line.strip() for line in entry.data.comment.split('\n')])
                rxn = TemplateReaction(reactants=entry.item.reactants[:], products=entry.item.products[:],
                                       specific_collider=entry.item.specific_collider, kinetics=entry.data,
                                       duplicate=entry.item.duplicate, reversible=entry.item.reversible,
                                       family=family, template=template, degeneracy=entry.item.degeneracy)
            else:  # pdep or standard library reaction
                rxn = LibraryReaction(reactants=entry.item.reactants[:], products=entry.item.products[:],
                                      library=self.label, specific_collider=entry.item.specific_collider,
                                      kinetics=entry.data, duplicate=entry.item.duplicate,
                                      reversible=entry.item.reversible, allow_pdep_route=entry.item.allow_pdep_route,
                                      elementary_high_p=entry.item.elementary_high_p)
            rxns.append(rxn)

        return rxns

    def mark_valid_duplicates(self, reactions1, reactions2):
        """
        Check for reactions that appear in both lists,
        and mark them as (valid) duplicates.
        """
        for r1 in reactions1:
            for r2 in reactions2:
                if (r1.reactants == r2.reactants and
                        r1.products == r2.products and
                        r1.specific_collider == r2.specific_collider and
                        r1.reversible == r2.reversible):
                    r1.duplicate = True
                    r2.duplicate = True

    def check_for_duplicates(self, mark_duplicates=False):
        """
        Check that all duplicate reactions in the kinetics library are
        properly marked (i.e. with their ``duplicate`` attribute set to 
        ``True``).
        If ``mark_duplicates`` is set to ``True``, then ignore and
        mark all duplicate reactions as duplicate.
        """
        for entry0 in self.entries.values():
            reaction0 = entry0.item
            if not reaction0.duplicate:
                # This reaction is not marked as a duplicate reaction
                # This means that if we find any duplicate reactions, it is an error
                for entry in self.entries.values():
                    reaction = entry.item
                    if reaction0 is not reaction and reaction0.is_isomorphic(reaction):
                        # We found a duplicate reaction that wasn't marked!
                        # RMG requires all duplicate reactions to be marked, unlike CHEMKIN
                        if mark_duplicates:
                            reaction0.duplicate = reaction.duplicate = True
                            logging.warning('Reaction indices {0} and {1} were marked as '
                                            'duplicate.'.format(entry0.index, entry.index))
                            continue

                        raise DatabaseError('Unexpected duplicate reaction {0} in kinetics library {1}. '
                                            'Reaction index {2} matches index {3}.'
                                            .format(reaction0, self.label, entry.index, entry0.index))

    def convert_duplicates_to_multi(self):
        """
        Merge all marked duplicate reactions in the kinetics library
        into single reactions with multiple kinetics.
        """
        logging.debug("Searching for duplicate reactions...")
        entries_to_remove = []
        for entry0 in self.entries.values():
            if entry0 in entries_to_remove:
                continue
            reaction0 = entry0.item
            if not reaction0.duplicate:
                continue
            logging.debug("Found a duplicate reaction: {0}".format(reaction0))
            duplicates = [entry0]
            for entry in self.entries.values():
                reaction = entry.item
                if reaction0 is reaction:
                    continue
                if reaction0.is_isomorphic(reaction, either_direction=False):
                    if reaction0.reversible != reaction.reversible:
                        logging.debug("Reactions isomorphic but with different reversibilities.")
                        continue
                    duplicates.append(entry)
            if len(duplicates) <= 1:
                continue
            kinetics_list = []
            long_desc = ''

            for entry in duplicates:
                kinetics = entry.data
                kinetics_list.append(kinetics)
                Tmin = kinetics.Tmin
                Tmax = kinetics.Tmax
                if kinetics.is_pressure_dependent():
                    Pmin = kinetics.Pmin
                    Pmax = kinetics.Pmax
                else:
                    Pmin = None
                    Pmax = None
                long_desc += entry.long_desc + '\n'

            if (len(kinetics_list) == 2 and
                    isinstance(kinetics_list[0], ThirdBody) and
                    not isinstance(kinetics_list[1], ThirdBody)):
                continue
            elif (len(kinetics_list) == 2 and
                    isinstance(kinetics_list[1], ThirdBody) and
                    not isinstance(kinetics_list[0], ThirdBody)):
                continue

            if all([isinstance(k, Arrhenius) for k in kinetics_list]):
                entry0.data = MultiArrhenius(arrhenius=kinetics_list, Tmin=Tmin, Tmax=Tmax)
            elif all([isinstance(k, PDepArrhenius) for k in kinetics_list]):
                entry0.data = MultiPDepArrhenius(arrhenius=kinetics_list, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax)
            else:
                logging.warning('Only Arrhenius and PDepArrhenius kinetics supported for duplicate reactions.')
                continue
            entry0.long_desc = long_desc
            entries_to_remove.extend(duplicates[1:])
        for entry in entries_to_remove:
            logging.debug("Removing duplicate reaction with index {0}.".format(entry.index))
            del (self.entries[entry.index])
        logging.debug("NB. the entries have not been renumbered, so these indices are missing.")

    def load(self, path, local_context=None, global_context=None):
        # Clear any previously-loaded data
        self.entries = OrderedDict()
        self.top = []

        # Set up global and local context
        if global_context is None:
            global_context = {}
        global_context['__builtins__'] = None
        global_context['True'] = True
        global_context['False'] = False
        if local_context is None:
            local_context = {}
        local_context['__builtins__'] = None
        local_context['entry'] = self.load_entry
        # local_context['tree'] = self._load_tree
        local_context['name'] = self.name
        local_context['solvent'] = self.solvent
        local_context['shortDesc'] = self.short_desc
        local_context['longDesc'] = self.long_desc
        local_context['metal'] = self.metal
        local_context['site'] = self.site
        local_context['facet'] = self.facet
        # add in anything from the Class level dictionary.
        for key, value in Database.local_context.items():
            local_context[key] = value

        # Process the file
        f = open(path, 'r')
        try:
            exec(f.read(), global_context, local_context)
        except Exception:
            logging.error('Error while reading database {0!r}.'.format(path))
            raise
        f.close()

        # Extract the database metadata
        self.name = local_context['name']
        self.solvent = local_context['solvent']
        self.short_desc = local_context['shortDesc']
        self.long_desc = local_context['longDesc'].strip()
        self.metal = local_context['metal']
        self.site = local_context['site']
        self.facet = local_context['facet']

        if 'autoGenerated' in local_context.keys():
            self.auto_generated = local_context['autoGenerated']

        # Load a unique set of the species in the kinetics library
        species_dict = self.get_species(os.path.join(os.path.dirname(path), 'dictionary.txt'))
        # Make sure all of the reactions draw from only this set
        entries = self.entries.values()
        for entry in entries:
            # Create a new reaction per entry
            rxn = entry.item
            rxn_string = entry.label

            # Convert coverage dependence label to species label
            if hasattr(entry.data, 'coverage_dependence'):
                if entry.data.coverage_dependence:
                    coverage_dependence_update = {}
                    for key, value in entry.data.coverage_dependence.items():
                        coverage_dependence_update[species_dict[key]] = value
                    entry.data.coverage_dependence = coverage_dependence_update

            # Convert the reactants and products to Species objects using the species_dict
            reactants, products = rxn_string.split('=>')
            reversible = True
            if '<=>' in rxn_string:
                reactants = reactants[:-1]
                products = products[1:]
            elif '=>' in rxn_string:
                products = products[1:]
                reversible = False
            if reversible != rxn.reversible:
                raise DatabaseError('Reaction string reversibility ({0}) and entry attribute `reversible` ({1}) '
                                    'must agree if reaction is irreversible.'.format(rxn.reversible, reversible))

            collider = re.search(r'\(\+[^\)]+\)', reactants)
            if collider is not None:
                collider = collider.group(0)  # save string value rather than the object
                if collider != re.search(r'\(\+[^\)]+\)', products).group(0):
                    raise ValueError('Third body colliders in reaction {0} in kinetics library {1} are missing from '
                                     'products or are not identical!'.format(rxn_string, self.label))
                extra_parenthesis = collider.count('(') - 1
                for i in range(extra_parenthesis):
                    # allow for species like N2(5) or CH2(T)(15) to be read as specific colliders,
                    # although currently not implemented in Chemkin. See RMG-Py #1070
                    collider += ')'
                reactants = reactants.replace(collider, '', 1)
                products = products.replace(collider, '', 1)
                if collider.upper().strip() != "(+M)":  # the collider is a specific species, not (+M) or (+m)
                    if collider.strip()[2:-1] not in species_dict:  # stripping spaces, '(+' and ')'
                        raise DatabaseError('Collider species {0} in kinetics library {1} is missing from its '
                                            'dictionary.'.format(collider.strip()[2:-1], self.label))
                    rxn.specific_collider = species_dict[collider.strip()[2:-1]]
            # verify there's no more than one specific_collider:
            collider = re.search(r'\(\+[^\)]+\)', reactants)
            if collider is not None:
                raise DatabaseError("Found TWO specific third body colliders, {0} and {1}, in reaction {2} in "
                                    "kinetics library {3), expecting no more than one!"
                                    .format(rxn.specific_collider, collider.group(0), rxn_string, self.label))

            for reactant in reactants.split('+'):
                reactant = reactant.strip()
                if reactant not in species_dict:
                    raise DatabaseError('Species {0} in kinetics library {1} is missing from its '
                                        'dictionary.'.format(reactant, self.label))
                rxn.reactants.append(species_dict[reactant])
            for product in products.split('+'):
                product = product.strip()
                if product not in species_dict:
                    raise DatabaseError('Species {0} in kinetics library {1} is missing from its '
                                        'dictionary.'.format(product, self.label))
                rxn.products.append(species_dict[product])

            if not rxn.is_balanced():
                raise DatabaseError('Reaction {0} in kinetics library {1} was not balanced! '
                                    'Please reformulate.'.format(rxn, self.label))

            if len(rxn.reactants) > 3:
                raise DatabaseError('RMG does not accept reactions with more than 3 reactants in its solver. '
                                    'Reaction {0} in kinetics library {1} has {2} reactants.'
                                    .format(rxn, self.label, len(rxn.reactants)))
            if len(rxn.products) > 3:
                raise DatabaseError('RMG does not accept reactions with more than 3 products in its solver. '
                                    'Reaction {0} in kinetics library {1} has {2} reactants.'
                                    .format(rxn, self.label, len(rxn.products)))

        if not self.auto_generated:
            self.check_for_duplicates()
            self.convert_duplicates_to_multi()

    def load_entry(self,
                   index,
                   label,
                   kinetics,
                   degeneracy=1,
                   duplicate=False,
                   reversible=True,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   allow_pdep_route=False,
                   elementary_high_p=False,
                   allow_max_rate_violation=False,
                   metal=None,
                   site=None,
                   facet=None,
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """

        # reactants = [Species(label=reactant1.strip().splitlines()[0].strip(), molecule=[Molecule().from_adjacency_list(reactant1)])]
        # if reactant2 is not None: reactants.append(Species(label=reactant2.strip().splitlines()[0].strip(), molecule=[Molecule().from_adjacency_list(reactant2)]))
        # if reactant3 is not None: reactants.append(Species(label=reactant3.strip().splitlines()[0].strip(), molecule=[Molecule().from_adjacency_list(reactant3)]))
        #
        # products = [Species(label=product1.strip().splitlines()[0].strip(), molecule=[Molecule().from_adjacency_list(product1)])]
        # if product2 is not None: products.append(Species(label=product2.strip().splitlines()[0].strip(), molecule=[Molecule().from_adjacency_list(product2)]))
        # if product3 is not None: products.append(Species(label=product3.strip().splitlines()[0].strip(), molecule=[Molecule().from_adjacency_list(product3)]))

        # Make a blank reaction
        rxn = Reaction(reactants=[], products=[], degeneracy=degeneracy, duplicate=duplicate, reversible=reversible,
                       allow_pdep_route=allow_pdep_route, elementary_high_p=elementary_high_p,
                       allow_max_rate_violation=allow_max_rate_violation)
        # if not rxn.is_balanced():
        #    raise DatabaseError('Reaction {0} in kinetics library {1} was not balanced! Please reformulate.'.format(rxn, self.label))
        # label = str(rxn)
        assert index not in self.entries, "Index of reaction {0} is not unique!".format(label)
        self.entries[index] = Entry(
            index=index,
            label=label,
            item=rxn,
            data=kinetics,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            metal=metal,
            site=site,
            facet=facet,
        )

    def save(self, path):
        """
        Save the current database to the file at location `path` on disk. 
        """
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass
        entries = self.get_entries_to_save()

        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}"\n'.format(self.name))
        f.write('shortDesc = "{0}"\n'.format(self.short_desc))
        f.write('longDesc = """\n')
        f.write(self.long_desc.strip() + '\n')
        f.write('"""\n')
        if self.metal:
            f.write('metal = "{0}"\n'.format(self.metal))
        if self.facet:
            f.write('facet = "{0}"\n'.format(self.facet))
        if self.site:
            f.write('site = "{0}"\n'.format(self.site))
        f.write('autoGenerated={0}\n'.format(self.auto_generated))

        for entry in entries:
            self.save_entry(f, entry)

        f.close()

        self.save_dictionary(os.path.join(os.path.split(path)[0], 'dictionary.txt'))

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the kinetics library to the file object `f`.
        """
        return save_entry(f, entry)

    def load_old(self, path):
        """
        Load an old-style RMG kinetics library from the location `path`.
        """
        path = os.path.abspath(path)

        self.load_old_dictionary(os.path.join(path, 'species.txt'), pattern=False)
        species = dict([(label, Species(label=label, molecule=[entry.item])) for label, entry in self.entries.items()])

        # Add common bath gases (Ar, Ne, He, N2) if not already present
        for label, smiles in [('Ar', '[Ar]'), ('He', '[He]'), ('Ne', '[Ne]'), ('N2', 'N#N')]:
            if label not in species:
                molecule = Molecule().from_smiles(smiles)
                spec = Species(label=label, molecule=[molecule])
                species[label] = spec

        reactions = []
        reactions.extend(self._load_old_reactions(os.path.join(path, 'reactions.txt'), species))
        if os.path.exists(os.path.join(path, 'pdepreactions.txt')):
            pdep_reactions = self._load_old_reactions(os.path.join(path, 'pdepreactions.txt'), species)
            # RMG-Py likes otherwise equivalent PDep and non-pdep reactions to be marked as duplicates
            self.mark_valid_duplicates(reactions, pdep_reactions)
            reactions.extend(pdep_reactions)

        self.entries = {}
        for index, reaction in enumerate(reactions):
            entry = Entry(
                index=index + 1,
                item=reaction,
                data=reaction.kinetics,
                label=str(reaction)
            )
            entry.long_desc = reaction.kinetics.comment
            reaction.kinetics.comment = ''
            self.entries[index + 1] = entry
            reaction.kinetics = None

        self.check_for_duplicates()
        self.convert_duplicates_to_multi()

    def _load_old_reactions(self, path, species):
        """
        Load an old-style reaction library from `path`. This algorithm can
        handle both the pressure-independent and pressure-dependent reaction
        files. If the pressure-dependent file is read, the extra pressure-
        dependent kinetics information is ignored unless the kinetics database
        is a seed mechanism.
        """
        from rmgpy.chemkin import read_reactions_block
        f = open(path, 'r')
        reaction_list = read_reactions_block(f, species_dict=species)
        f.close()
        return reaction_list

    def save_old(self, path):
        """
        Save an old-style reaction library to `path`. This creates files named
        ``species.txt``, ``reactions.txt``, and ``pdepreactions.txt`` in the
        given directory; these contain the species dictionary, high-pressure
        limit reactions and kinetics, and pressure-dependent reactions and
        kinetics, respectively.
        """
        try:
            os.makedirs(path)
        except OSError:
            pass

        def write_arrhenius(f, arrhenius):
            f.write(' {0:<12.3E} {1:>7.3f} {2:>11.2f}    {3}{4:g} {5:g} {6:g}\n'.format(
                arrhenius.A.value_si,
                arrhenius.n.value_si,
                arrhenius.Ea.value_si / 4.184,
                '*' if arrhenius.A.is_uncertainty_multiplicative() else '',
                arrhenius.A.uncertainty,
                arrhenius.n.uncertainty,
                arrhenius.Ea.uncertainty / 4.184,
            ))

        # Gather all of the species used in this kinetics library
        species_dict = self.get_species()
        # Also include colliders in the above
        for entry in self.entries.values():
            if isinstance(entry.data, ThirdBody):
                for molecule in entry.data.efficiencies:
                    formula = molecule.get_formula()
                    if formula in ['He', 'Ar', 'N2', 'Ne']:
                        pass
                    else:
                        found = False
                        for species in species_dict.values():
                            for mol in species.molecule:
                                if mol.is_isomorphic(molecule):
                                    found = True
                                    break
                        if not found:
                            species_dict[formula] = Species(label=formula, molecule=[molecule])

        entries = list(self.entries.values())
        entries.sort(key=lambda x: x.index)

        # Save the species dictionary
        species_list = list(species_dict.values())
        species_list.sort(key=lambda x: x.label)
        f = open(os.path.join(path, 'species.txt'), 'w')
        for species in species_list:
            f.write(species.molecule[0].to_adjacency_list(label=species.label, remove_h=False) + "\n")
        f.close()

        # Save the high-pressure limit reactions
        # Currently only Arrhenius kinetics are allowed
        f = open(os.path.join(path, 'reactions.txt'), 'w')
        f.write('Unit:\n')
        f.write('A: mol/m3/s\n')
        f.write('E: cal/mol\n\n')
        f.write('Reactions:\n')
        for entry in entries:
            kinetics = entry.data
            rate_list = []
            if isinstance(kinetics, MultiArrhenius):
                entry.item.duplicate = True
                rate_list = kinetics.arrhenius[:]
            else:
                if not kinetics.is_pressure_dependent():
                    rate_list.append(kinetics)
            for rate in rate_list:
                # Write reaction equation
                f.write('{0:<59}'.format(entry.item))
                # Write kinetics
                if isinstance(rate, Arrhenius):
                    write_arrhenius(f, rate)
                else:
                    raise DatabaseError('Unexpected kinetics type "{0}" encountered while saving old kinetics library '
                                        '(reactions.txt).'.format(rate.__class__))
                # Mark as duplicate if needed
                if entry.item.duplicate:
                    f.write(' DUPLICATE\n')
        f.close()

        # Save the pressure-dependent reactions
        # Currently only ThirdBody, Lindemann, Troe, and PDepArrhenius kinetics are allowed
        f = open(os.path.join(path, 'pdepreactions.txt'), 'w')
        f.write('Unit:\n')
        f.write('A: mol/m3/s\n')
        f.write('E: cal/mol\n\n')
        f.write('Reactions:\n')
        for entry in entries:
            kinetics = entry.data
            if not kinetics.is_pressure_dependent():
                continue
            rate_list = []
            if isinstance(kinetics, MultiPDepArrhenius):
                entry.item.duplicate = True
                rate_list = kinetics.arrhenius[:]
            else:
                rate_list.append(kinetics)
            for rate in rate_list:
                # Write reaction equation
                equation = str(entry.item)
                if entry.item.reversible:
                    index = equation.find('<=>')
                else:
                    index = equation.find('=>')
                if isinstance(rate, ThirdBody) and not isinstance(rate, Lindemann):
                    equation = '{0}+ M {1} + M'.format(equation[0:index], equation[index:])
                elif isinstance(rate, PDepArrhenius):
                    pass
                else:
                    equation = '{0}(+M) {1} (+M)'.format(equation[0:index], equation[index:])
                f.write('{0:<59}'.format(equation))
                # Write kinetics
                if isinstance(rate, (ThirdBody, Lindemann, Troe)):
                    if isinstance(rate, Lindemann):
                        # Lindemann (and Troe) fall-off have the High-P as default, and Low-P labeled LOW
                        write_arrhenius(f, rate.arrheniusHigh)
                    else:
                        # Non-falloff ThirdBody reactions are always in the Low-P limit
                        write_arrhenius(f, rate.arrheniusLow)
                    if len(rate.efficiencies) > 0:
                        eff_line = ''
                        for molecule, efficiency in rate.efficiencies.items():
                            for spec in species_dict.values():
                                if molecule in spec.molecule:
                                    mol_label = spec.label
                                    break
                            else:
                                mol_label = molecule.get_formula().upper()
                            eff_line += '{0}/{1:g}/  '.format(mol_label, efficiency)
                        f.write(eff_line.strip() + '\n')
                    if isinstance(rate, Lindemann):
                        f.write('     LOW  /  {0:10.3e} {1:9.3f} {2:10.2f}/\n'.format(
                            rate.arrheniusLow.A.value_si,
                            rate.arrheniusLow.n.value_si,
                            rate.arrheniusLow.Ea.value_si / 4.184,
                        ))
                    if isinstance(rate, Troe):
                        if rate.T2 is not None:
                            f.write('     TROE /  {0:10.4f} {1:10.2g} {2:10.2g} {3:10.2g}/\n'.format(
                                rate.alpha,
                                rate.T3.value_si,
                                rate.T1.value_si,
                                rate.T2.value_si,
                            ))
                        else:
                            f.write('     TROE /  {0:10.4f} {1:10.2g} {2:10.2g}/\n'.format(
                                rate.alpha,
                                rate.T3.value_si,
                                rate.T1.value_si,
                            ))

                elif isinstance(rate, PDepArrhenius):
                    write_arrhenius(f, rate.arrhenius[-1])
                    for pressure, arrhenius in zip(rate.pressures.value_si, rate.arrhenius):
                        f.write('     PLOG /  {0:10g} {1:10.3e} {2:9.3f} {3:10.2f} /\n'.format(
                            pressure / 1e5,
                            arrhenius.A.value_si,
                            arrhenius.n.value_si,
                            arrhenius.Ea.value_si / 4.184,
                        ))
                else:
                    raise DatabaseError('Unexpected kinetics type "{0}" encountered while saving old kinetics library '
                                        '(reactions.txt).'.format(rate.__class__))
                # Mark as duplicate if needed
                if entry.item.duplicate:
                    f.write(' DUPLICATE\n')
                f.write('\n')
        f.close()
