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
This module contains classes and functions that are used by multiple modules
in this subpackage.
"""
import itertools
import logging

from rmgpy.data.base import LogicNode
from rmgpy.exceptions import DatabaseError
from rmgpy.molecule import Group, Molecule
from rmgpy.molecule.fragment import Fragment
from rmgpy.reaction import Reaction
from rmgpy.species import Species


################################################################################


def save_entry(f, entry):
    """
    Save an `entry` in the kinetics database by writing a string to
    the given file object `f`.
    """

    def sort_efficiencies(efficiencies0):
        efficiencies = {}
        for mol, eff in efficiencies0.items():
            if isinstance(mol, str):
                # already in SMILES string format
                smiles = mol
            else:
                smiles = mol.to_smiles()

            efficiencies[smiles] = eff
        keys = list(efficiencies.keys())
        keys.sort()
        return [(key, efficiencies[key]) for key in keys]

    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    if entry.label != '':
        f.write('    label = "{0}",\n'.format(entry.label))

    # Entries for kinetic rules, libraries, training reactions
    # and depositories will have a Reaction object for its item
    if isinstance(entry.item, Reaction):
        # Write out additional data if depository or library
        # kinetic rules would have a Group object for its reactants instead of Species
        if isinstance(entry.item.reactants[0], Species):
            # Add degeneracy if the reaction is coming from a depository or kinetics library
            f.write('    degeneracy = {0:.1f},\n'.format(entry.item.degeneracy))
            if entry.item.duplicate:
                f.write('    duplicate = {0!r},\n'.format(entry.item.duplicate))
            if not entry.item.reversible:
                f.write('    reversible = {0!r},\n'.format(entry.item.reversible))
            if entry.item.allow_pdep_route:
                f.write('    allow_pdep_route = {0!r},\n'.format(entry.item.allow_pdep_route))
            if entry.item.elementary_high_p:
                f.write('    elementary_high_p = {0!r},\n'.format(entry.item.elementary_high_p))
            if entry.item.allow_max_rate_violation:
                f.write('    allow_max_rate_violation = {0!r},\n'.format(entry.item.allow_max_rate_violation))
    # Entries for groups with have a group or logicNode for its item
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.to_adjacency_list())
        f.write('""",\n')
    elif isinstance(entry.item, LogicNode):
        f.write('    group = "{0}",\n'.format(entry.item))
    else:
        raise DatabaseError("Encountered unexpected item of type {0} while "
                            "saving database.".format(entry.item.__class__))

    # Write kinetics
    if isinstance(entry.data, str):
        f.write('    kinetics = "{0}",\n'.format(entry.data))
    elif entry.data is not None:
        efficiencies = None
        if hasattr(entry.data, 'efficiencies'):
            efficiencies = entry.data.efficiencies
            entry.data.efficiencies = dict(sort_efficiencies(entry.data.efficiencies))
        kinetics = repr(entry.data)  # todo prettify currently does not support uncertainty attribute
        kinetics = '    kinetics = {0},\n'.format(kinetics.replace('\n', '\n    '))
        f.write(kinetics)
        if hasattr(entry.data, 'efficiencies'):
            entry.data.efficiencies = efficiencies
    else:
        f.write('    kinetics = None,\n')

    # Write reference
    if entry.reference is not None:
        reference = entry.reference.to_pretty_repr()
        lines = reference.splitlines()
        f.write('    reference = {0}\n'.format(lines[0]))
        for line in lines[1:-1]:
            f.write('    {0}\n'.format(line))
        f.write('    ),\n'.format(lines[0]))

    if entry.reference_type != "":
        f.write('    referenceType = "{0}",\n'.format(entry.reference_type))
    if entry.rank is not None:
        f.write('    rank = {0},\n'.format(entry.rank))

    if entry.short_desc.strip() != '':
        f.write(f'    shortDesc = """{entry.short_desc.strip()}""",\n')
    if entry.long_desc.strip() != '':
        f.write(f'    longDesc = \n"""\n{entry.long_desc.strip()}\n""",\n')

    # write metal attributes
    if entry.metal:
        f.write('    metal = "{0}",\n'.format(entry.metal))
    if entry.facet:
        f.write('    facet = "{0}",\n'.format(entry.facet))
    if entry.site:
        f.write('    site = "{0}",\n'.format(entry.site))

    f.write(')\n\n')


def ensure_species(input_list, resonance=False, keep_isomorphic=False):
    """
    The input list of :class:`Species` or :class:`Molecule` objects is modified
    in place to only have :class:`Species` objects. Returns None.
    """
    for index, item in enumerate(input_list):
        if isinstance(item, Molecule) or isinstance(item, Fragment):
            new_item = Species(molecule=[item])
        elif isinstance(item, Species):
            new_item = item
        else:
            raise TypeError('Only Molecule or Species objects can be handled.')
        if resonance:
            if not any([mol.reactive for mol in new_item.molecule]):
                # if generating a reaction containing a Molecule with a reactive=False flag (e.g., for degeneracy
                # calculations), that was now converted into a Species, first mark as reactive=True
                new_item.molecule[0].reactive = True
            new_item.generate_resonance_structures(keep_isomorphic=keep_isomorphic)
        input_list[index] = new_item


def generate_molecule_combos(input_species):
    """
    Generate combinations of molecules from the given species objects.
    """
    if len(input_species) == 1:
        combos = [(mol,) for mol in input_species[0].molecule]
    elif len(input_species) == 2:
        combos = itertools.product(input_species[0].molecule, input_species[1].molecule)
    elif len(input_species) == 3:
        combos = itertools.product(input_species[0].molecule, input_species[1].molecule, input_species[2].molecule)
    else:
        raise ValueError('Reaction generation can be done for 1, 2, or 3 species, not {0}.'.format(len(input_species)))

    return combos


def ensure_independent_atom_ids(input_species, resonance=True):
    """
    Given a list or tuple of :class:`Species` or :class:`Molecule` objects,
    ensure that atom ids are independent.
    The `resonance` argument can be set to False to not generate
    resonance structures.

    Modifies the list in place (replacing :class:`Molecule` with :class:`Species`).
    Returns None.
    """
    ensure_species(input_species)  # do not generate resonance structures since we do so below

    # Method to check that all species' atom ids are different
    def independent_ids():
        num_atoms = 0
        ids = []
        for spcs in input_species:
            num_atoms += len(spcs.molecule[0].atoms)
            ids.extend([atom.id for atom in spcs.molecule[0].atoms])
        num_id = len(set(ids))
        return num_id == num_atoms

    # If they are not all different, reassign ids and remake resonance structures
    if not independent_ids():
        for species in input_species:
            unreactive_mol_list = [mol for mol in species.molecule if not mol.reactive]
            mol = [mol for mol in species.molecule if mol.reactive][0]  # Choose first reactive molecule
            mol.assign_atom_ids()
            species.molecule = [mol]
            # Remake resonance structures with new labels
            if resonance:
                species.generate_resonance_structures(keep_isomorphic=True)
            if len(unreactive_mol_list):
                species.molecule.extend(unreactive_mol_list)
    elif resonance:
        # IDs are already independent, generate resonance structures if needed
        for species in input_species:
            species.generate_resonance_structures(keep_isomorphic=True)


def find_degenerate_reactions(rxn_list, same_reactants=None, template=None, kinetics_database=None,
                              kinetics_family=None, save_order=False):
    """
    Given a list of Reaction objects, this method combines degenerate
    reactions and increments the reaction degeneracy value. For multiple
    transition states, this method keeps them as duplicate reactions.

    If a template is specified, then the reaction list will be filtered
    to leave only reactions which match the specified template, then the
    degeneracy will be calculated as usual.

    A KineticsDatabase or KineticsFamily instance can also be provided to
    calculate the degeneracy for reactions generated in the reverse direction.
    If not provided, then it will be retrieved from the global database.

    This algorithm used to exist in family._generate_reactions, but was moved
    here so it could operate across reaction families.

    This method returns an updated list with degenerate reactions removed.

    Args:
        rxn_list (list):                                reactions to be analyzed
        same_reactants (bool, optional):                indicate whether the reactants are identical
        template (list, optional):                      specify a specific template to filter by
        kinetics_database (KineticsDatabase, optional): provide a KineticsDatabase instance for calculating degeneracy
        kinetics_family (KineticsFamily, optional):     provide a KineticsFamily instance for calculating degeneracy
        save_order (bool, optional):                    reset atom order after performing atom isomorphism

    Returns:
        Reaction list with degenerate reactions combined with proper degeneracy values
    """
    # If a specific reaction template is requested, filter by that template
    if template is not None:
        selected_rxns = []
        template = frozenset(template)
        for rxn in rxn_list:
            if template == frozenset(rxn.template):
                selected_rxns.append(rxn)
        if not selected_rxns:
            # Only log a warning here. If a non-empty output is expected, then the caller should raise an exception
            logging.warning('No reactions matched the specified template, {0}'.format(template))
            return []
    else:
        selected_rxns = rxn_list

    # We want to sort all the reactions into sublists composed of isomorphic reactions
    # with degenerate transition states
    sorted_rxns = []
    for rxn0 in selected_rxns:
        rxn0.ensure_species(save_order=save_order)
        if len(sorted_rxns) == 0:
            # This is the first reaction, so create a new sublist
            sorted_rxns.append([rxn0])
        else:
            # Loop through each sublist, which represents a unique reaction
            for sub_list in sorted_rxns:
                # Try to determine if the current rxn0 is identical or isomorphic to any reactions in the sublist
                isomorphic = False
                identical = False
                same_template = True
                for rxn in sub_list:
                    isomorphic = rxn0.is_isomorphic(rxn, check_identical=False, strict=False,
                                                    check_template_rxn_products=True, save_order=save_order)
                    if isomorphic:
                        identical = rxn0.is_isomorphic(rxn, check_identical=True, strict=False,
                                                       check_template_rxn_products=True, save_order=save_order)
                        if identical:
                            # An exact copy of rxn0 is already in our list, so we can move on
                            break
                        same_template = frozenset(rxn.template) == frozenset(rxn0.template)
                    else:
                        # This sublist contains a different product
                        break

                # Process the reaction depending on the results of the comparisons
                if identical:
                    # This reaction does not contribute to degeneracy
                    break
                elif isomorphic:
                    if same_template:
                        # We found the right sublist, and there is no identical reaction
                        # We should add rxn0 to the sublist as a degenerate rxn, and move on to the next rxn
                        sub_list.append(rxn0)
                        break
                    else:
                        # We found an isomorphic sublist, but the reaction templates are different
                        # We need to mark this as a duplicate and continue searching the remaining sublists
                        rxn0.duplicate = True
                        sub_list[0].duplicate = True
                        continue
                else:
                    # This is not an isomorphic sublist, so we need to continue searching the remaining sublists
                    # Note: This else statement is not technically necessary but is included for clarity
                    continue
            else:
                # We did not break, which means that there was no isomorphic sublist, so create a new one
                sorted_rxns.append([rxn0])

    rxn_list = []
    for sub_list in sorted_rxns:
        # Collapse our sorted reaction list by taking one reaction from each sublist
        rxn = sub_list[0]
        # The degeneracy of each reaction is the number of reactions that were in the sublist
        rxn.degeneracy = sum([reaction0.degeneracy for reaction0 in sub_list])
        rxn_list.append(rxn)

    for rxn in rxn_list:
        if rxn.is_forward:
            reduce_same_reactant_degeneracy(rxn, same_reactants)
        else:
            # fix the degeneracy of (not ownReverse) reactions found in the backwards direction
            try:
                family = kinetics_family or kinetics_database.families[rxn.family]
            except AttributeError:
                from rmgpy.data.rmg import get_db
                family = get_db('kinetics').families[rxn.family]
            if not family.own_reverse:
                rxn.degeneracy = family.calculate_degeneracy(rxn)

    return rxn_list


def reduce_same_reactant_degeneracy(reaction, same_reactants=None):
    """
    This method reduces the degeneracy of reactions with identical reactants,
    since translational component of the transition states are already taken
    into account (so swapping the same reactant is not valid)

    same_reactants can be None or an integer. If it is None, then isomorphism
    checks will be done to determine if the reactions are the same. If it is an
    integer, that integer denotes the number of reactants that are isomorphic.

    This comes from work by Bishop and Laidler in 1965
    """
    if not (same_reactants == 0 or same_reactants == 1):
        if len(reaction.reactants) == 2:
            if ((reaction.is_forward and same_reactants == 2) or
                    reaction.reactants[0].is_isomorphic(reaction.reactants[1])):
                reaction.degeneracy *= 0.5
                logging.debug(
                    'Degeneracy of reaction {} was decreased by 50% to {} since the reactants are identical'.format(
                        reaction, reaction.degeneracy)
                )
        elif len(reaction.reactants) == 3:
            if reaction.is_forward:
                if same_reactants == 3:
                    reaction.degeneracy /= 6.0
                    logging.debug(
                        'Degeneracy of reaction {} was divided by 6 to give {} since all of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
                elif same_reactants == 2:
                    reaction.degeneracy *= 0.5
                    logging.debug(
                        'Degeneracy of reaction {} was decreased by 50% to {} since two of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
            else:
                same_01 = reaction.reactants[0].is_isomorphic(reaction.reactants[1])
                same_02 = reaction.reactants[0].is_isomorphic(reaction.reactants[2])
                if same_01 and same_02:
                    reaction.degeneracy /= 6.0
                    logging.debug(
                        'Degeneracy of reaction {} was divided by 6 to give {} since all of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
                elif same_01 or same_02:
                    reaction.degeneracy *= 0.5
                    logging.debug(
                        'Degeneracy of reaction {} was decreased by 50% to {} since two of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
                elif reaction.reactants[1].is_isomorphic(reaction.reactants[2]):
                    reaction.degeneracy *= 0.5
                    logging.debug(
                        'Degeneracy of reaction {} was decreased by 50% to {} since two of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
