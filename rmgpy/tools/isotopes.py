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
This module contains functionality for generating mechanisms with isotopes.
"""

from __future__ import division, print_function

import logging
import math
import os
import os.path
from copy import copy, deepcopy

import numpy as np
import pandas as pd

import rmgpy.constants as constants
import rmgpy.molecule.element
from rmgpy.chemkin import ChemkinWriter
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import get_db
from rmgpy.data.thermo import find_cp0_and_cpinf
from rmgpy.kinetics.arrhenius import MultiArrhenius
from rmgpy.molecule import Molecule
from rmgpy.molecule.element import get_element
from rmgpy.reaction import Reaction, same_species_lists
from rmgpy.rmg.main import RMG, initializeLog
from rmgpy.species import Species
from rmgpy.thermo.thermoengine import processThermoData
from rmgpy.tools.loader import loadRMGJob


def initialize_isotope_model(rmg, isotopes):
    """
    Initialize the RMG object by using the parameter species list
    as initial species instead of the species from the RMG input file.
    """
    # Read input file
    rmg.loadInput(rmg.inputFile)

    # Check input file
    rmg.checkInput()

    # Load databases
    rmg.loadDatabase()

    logging.info("isotope: Adding the isotopomers into the RMG model")
    for isotopomers in isotopes:
        for spc in isotopomers:
            spec, isNew = rmg.reactionModel.makeNewSpecies(spc)
            spec.thermo = spc.thermo
            if isNew:
                rmg.reactionModel.addSpeciesToEdge(spec)
                rmg.initialSpecies.append(spec)
    logging.info("isotope: Adding standard species into the model")
    for spec in rmg.initialSpecies:
        spec.thermo = processThermoData(spec, spec.thermo)
        if not spec.reactive:
            rmg.reactionModel.enlarge(spec)
    for spec in rmg.initialSpecies:
        if spec.reactive:
            rmg.reactionModel.enlarge(spec)
    logging.info("isotope: Finalizing the species additions")
    rmg.initializeReactionThresholdAndReactFlags()
    rmg.reactionModel.initializeIndexSpeciesDict()


def generate_isotope_model(outputDirectory, rmg0, isotopes, useOriginalReactions=False,
                           kineticIsotopeEffect=None):
    """
    Replace the core species of the rmg model with the parameter list
    of species.

    Generate all reactions between new list of core species.

    Returns created RMG object.
    """
    logging.debug("isotope: called generateIsotopeModel")
    rmg = RMG(inputFile=rmg0.inputFile, outputDirectory=outputDirectory)
    rmg.attach(ChemkinWriter(outputDirectory))

    logging.info("isotope: making the isotope model for with all species")
    initialize_isotope_model(rmg, isotopes)

    if useOriginalReactions:
        logging.info("isotope: finding reactions from the original reactions")
        rxns = generate_isotope_reactions(rmg0.reactionModel.core.reactions, isotopes)
        rmg.reactionModel.processNewReactions(rxns, newSpecies=[])

    else:
        logging.info("isotope: enlarging the isotope model")
        rmg.reactionModel.enlarge(reactEdge=True,
                                  unimolecularReact=rmg.unimolecularReact,
                                  bimolecularReact=rmg.bimolecularReact)

    logging.info("isotope: clustering reactions")
    clusters = cluster(rmg.reactionModel.core.reactions)
    logging.info('isotope: fixing the directions of every reaction to a standard')
    for isotopomerRxnList in clusters:
        ensure_reaction_direction(isotopomerRxnList)

    consistent = True
    logging.info("isotope: checking symmetry is consistent among isotopomers")
    for species_list in cluster(rmg.reactionModel.core.species):
        if not ensure_correct_symmetry(species_list):
            logging.info("isotopomers of {} with index {} may have wrong "
                         "symmetry".format(species_list[0], species_list[0].index))
            consistent = False
    logging.info("isotope: checking that reaction degeneracy is consistent among isotopomers")
    for rxn_list in clusters:
        if not ensure_correct_degeneracies(rxn_list):
            logging.info("isotopomers of {} with index {} may have incorrect "
                         "degeneracy.".format(rxn_list[0], rxn_list[0].index))
            consistent = False
    if not consistent:
        logging.warning("isotope: non-consistent degeneracy and/or symmetry was detected. This may lead to "
                        "unrealistic deviations in enrichment. check log for more details")

    if kineticIsotopeEffect:
        logging.info('isotope: modifying reaction rates using kinetic isotope effect '
                     'method "{0}"'.format(kineticIsotopeEffect))
        if kineticIsotopeEffect == 'simple':
            apply_kinetic_isotope_effect_simple(clusters, rmg.database.kinetics)
        else:
            logging.warning('isotope: kinetic isotope effect {0} is not supported. '
                            'skipping adding kinetic isotope effects.')
    else:
        logging.info('isotope: not adding kinetic isotope effects since no method was supplied.')
    logging.info("isotope: saving files")
    rmg.saveEverything()

    rmg.finish()

    return rmg


def generate_isotope_reactions(isotopeless_reactions, isotopes):
    """
    Find the list of isotope reactions based on the reactions in the isotopeless
    reaction.

    uses the reactSpecies method to find reactions with proper degeneracies and
    then filters out those that don't match products. the proper reactions are
    given kinetics of the previous reaction modified for the degeneracy difference.
    """
    # make sure all isotopeless reactions have templates and are TemplateReaction objects
    for rxn in isotopeless_reactions:
        if not isinstance(rxn, TemplateReaction):
            raise TypeError('reactions sent to generate_isotope_reactions must be a TemplateReaction object')
        if rxn.template is None:
            raise AttributeError('isotope reaction {0} does not have a template attribute. '
                                 'The object is:\n\n{1}'.format(str(rxn), repr(rxn)))

    found_reactions = []
    rxn_index = 0
    while rxn_index < len(isotopeless_reactions):
        rxn = isotopeless_reactions[rxn_index]
        # find all reactions involving same reactants
        rxns_w_same_reactants = [rxn]
        rxn_index2 = rxn_index + 1
        while rxn_index2 < len(isotopeless_reactions):
            if same_species_lists(isotopeless_reactions[rxn_index].reactants,
                                  isotopeless_reactions[rxn_index2].reactants, ):
                rxns_w_same_reactants.append(isotopeless_reactions[rxn_index2])
                del isotopeless_reactions[rxn_index2]
            else:
                rxn_index2 += 1
        ##### find all pairs of reacting isotoper species #####
        # find the lists of reactants that have identical isotopomers
        reactants = []
        for reactant in rxn.reactants:
            for iso_index, isotopomers in enumerate(isotopes):
                if compare_isotopomers(reactant, isotopomers[0]):
                    reactants.append(iso_index)
                    break
        # find pairs of all reactants to react together
        reactant_pairs = []
        if len(rxn.reactants) == 1:
            reactant_pairs = [[spec] for spec in isotopes[reactants[0]]]
        elif len(rxn.reactants) == 2:
            for spec1 in isotopes[reactants[0]]:
                for spec2 in isotopes[reactants[1]]:
                    reactant_pairs.append([spec1, spec2])
        else:
            raise ValueError('Cannot process reactions with over 2 reactants')

        # remove identical pairs
        rxn_index3 = 0
        while rxn_index3 < len(reactant_pairs):
            rxn_index4 = rxn_index3 + 1
            while rxn_index4 < len(reactant_pairs):
                if same_species_lists(reactant_pairs[rxn_index3],
                                      reactant_pairs[rxn_index4]):
                    del reactant_pairs[rxn_index4]
                else:
                    rxn_index4 += 1
            rxn_index3 += 1

        # make reaction objects
        for pair in reactant_pairs:
            # copy species so they don't get modified
            species_tuple = tuple([spc.copy(deep=True) for spc in pair])
            unfiltered_rxns = get_db('kinetics').generate_reactions_from_families(species_tuple,
                                                                                  only_families=[rxn.family])
            # remove reactions whose products don't match the original reactions
            rxn_index5 = 0
            while rxn_index5 < len(unfiltered_rxns):
                for isotopeless_reaction in rxns_w_same_reactants:
                    isotopeless_kinetics = isotopeless_reaction.kinetics
                    isotopeless_degeneracy = isotopeless_reaction.degeneracy
                    if compare_isotopomers(isotopeless_reaction, unfiltered_rxns[rxn_index5], eitherDirection=False) \
                            and isotopeless_reaction.family == unfiltered_rxns[rxn_index5].family \
                            and frozenset(isotopeless_reaction.template) == \
                            frozenset(unfiltered_rxns[rxn_index5].template):
                        # apply kinetics to new reaction & modify for degeneracy
                        unfiltered_rxns[rxn_index5].kinetics = deepcopy(isotopeless_kinetics)
                        unfiltered_rxns[rxn_index5].kinetics.change_rate(
                            unfiltered_rxns[rxn_index5].degeneracy / isotopeless_degeneracy)
                        rxn_index5 += 1
                        break
                else:  # did not find same prodcuts
                    del unfiltered_rxns[rxn_index5]
            found_reactions.extend(unfiltered_rxns)
        rxn_index += 1
    return found_reactions


def generate_isotopomers(spc, N=1):
    """
    Generate all isotopomers of the parameter species by adding max. N carbon isotopes to the
    atoms of the species.
    """

    mol = spc.molecule[0]
    isotope = get_element(6, 13)

    mols = []
    add_isotope(0, N, mol, mols, isotope)

    spcs = []
    for isomol in mols:
        isotopomer = Species(molecule=[isomol], thermo=deepcopy(spc.thermo), transportData=spc.transportData,
                             reactive=spc.reactive)
        isotopomer.generate_resonance_structures(keep_isomorphic=True)
        spcs.append(isotopomer)

    # do not retain identical species:
    filtered = []
    while spcs:
        candidate = spcs.pop()
        unique = True
        for isotopomer in filtered:
            if isotopomer.is_isomorphic(candidate):
                unique = False
                break
        if unique: filtered.append(candidate)

    if spc.thermo:
        for isotopomer in filtered:
            correct_entropy(isotopomer, spc)

    return filtered


def add_isotope(i, N, mol, mols, element):
    """
    Iterate over the atoms of the molecule, and changes the element object
    of the atom by the provide parameter element object. Add the newly created
    isotopomer to the list of Molecule objects. For each created isotopomer,
    recursively call the method, until the maximum number of isotopes per molecule
    (N) is reached.

    """
    if i == N:
        return
    else:
        atoms = [at for at in mol.atoms if at.symbol == element.symbol]
        for at in atoms:
            if at.element == element:
                continue
            else:
                isotopomer = mol.copy(deep=True)
                isotopomer.atoms[mol.atoms.index(at)].element = element
                mols.append(isotopomer)
                add_isotope(i + 1, N, isotopomer, mols, element)


def cluster(objList):
    """
    Creates subcollections of isotopomers/reactions that
    only differ in their isotopic labeling.

    This method works for either species or reactions.

    It is O(n^2) efficient
    """

    unclustered = copy(objList)

    # [[list of Species objs]]
    clusters = []

    while unclustered:
        candidate = unclustered.pop()
        for cluster in clusters:
            if compare_isotopomers(cluster[0], candidate):
                cluster.append(candidate)
                break
        else:
            clusters.append([candidate])

    return clusters


def remove_isotope(labeledObj, inplace=False):
    """
    Create a deep copy of the first molecule of the species object and replace
    non-normal Element objects (of special isotopes) by the
    expected isotope.

    If the boolean `inplace` is True, the method remove the isotopic atoms of
    the Species/Reaction
    inplace and returns a list of atom objects & element pairs for adding back
    to the oritinal object. This should significantly improve speed of this method.

    If successful, the non-inplace parts should be removed
    """

    if isinstance(labeledObj, Species):
        if inplace:
            modified_atoms = []
            for mol in labeledObj.molecule:
                for atom in mol.atoms:
                    if atom.element.isotope != -1:
                        modified_atoms.append((atom, atom.element))
                        atom.element = get_element(atom.element.symbol)
            return modified_atoms
        else:
            stripped = labeledObj.copy(deep=True)

            for atom in stripped.molecule[0].atoms:
                if atom.element.isotope != -1:
                    atom.element = get_element(atom.element.symbol)

            # only do it for the first molecule, generate the other resonance isomers.
            stripped.molecule = [stripped.molecule[0]]
            stripped.generate_resonance_structures(keep_isomorphic=True)

        return stripped

    elif isinstance(labeledObj, Reaction):

        if inplace:

            atom_list = []
            for reactant in labeledObj.reactants:
                removed = remove_isotope(reactant, inplace)
                if removed:
                    atom_list += removed
            for product in labeledObj.products:
                removed = remove_isotope(product, inplace)
                if removed:
                    atom_list += removed

            return atom_list
        else:
            stripped_rxn = labeledObj.copy()

            stripped_reactants = []
            for reactant in stripped_rxn.reactants:
                stripped_reactants.append(remove_isotope(reactant, inplace))
            stripped_rxn.reactants = stripped_reactants

            stripped_products = []
            for product in stripped_rxn.products:
                stripped_products.append(remove_isotope(product, inplace))
            stripped_rxn.products = stripped_products

            return stripped_rxn
    elif isinstance(labeledObj, Molecule):
        if inplace:
            modified_atoms = []
            for atom in labeledObj.atoms:
                if atom.element.isotope != -1:
                    modified_atoms.append((atom, atom.element))
                    atom.element = get_element(atom.element.symbol)
            return modified_atoms
        else:
            stripped = labeledObj.copy(deep=True)

            for atom in stripped.atoms:
                if atom.element.isotope != -1:
                    atom.element = get_element(atom.element.symbol)

            return stripped
    else:
        raise TypeError('Only Reaction, Species, and Molecule objects are supported')


def ensure_reaction_direction(isotopomerRxns):
    """
    given a list of reactions with varying isotope labels but identical structure,
    obtained from the `cluster` method, this method remakes the kinetics so that
    they all face the same direction.
    """

    # find isotopeless reaction as standard
    reference = isotopomerRxns[0]
    family = get_db('kinetics').families[reference.family]
    if family.ownReverse:
        for rxn in isotopomerRxns:
            if not compare_isotopomers(rxn, reference, eitherDirection=False):
                # the reaction is in the oposite direction
                logging.info('isotope: identified flipped reaction direction in reaction number {} of reaction {}. '
                             'Altering the direction.'.format( rxn.index, str(rxn)))
                # obtain reverse attribute with template and degeneracy
                family.add_reverse_attribute(rxn)
                if frozenset(rxn.reverse.template) != frozenset(reference.template):
                    logging.warning("Reaction {} did not find proper reverse template, might cause "
                                    "degeneracy error.".format(str(rxn)))
                # reverse reactants and products of original reaction
                rxn.reactants, rxn.products = rxn.products, rxn.reactants
                rxn.pairs = [(p, r) for r, p in rxn.pairs]
                # set degeneracy to isotopeless reaction
                rxn.degeneracy = reference.degeneracy
                # make this reaction have kinetics of isotopeless reaction
                new_kinetics = deepcopy(reference.kinetics)
                rxn.kinetics = new_kinetics
                rxn.template = reference.template
                # set degeneracy to new reaction
                rxn.degeneracy = rxn.reverse.degeneracy
                # delete reverse attribute
                rxn.reverse = None


def redo_isotope(atomList):
    """
    This takes a list of zipped atoms with their isotopes removed, from
    and elements.
    """
    for atom, element in atomList:
        atom.element = element


def compare_isotopomers(obj1, obj2, eitherDirection=True):
    """
    This method takes two species or reaction objects and returns true if
    they only differ in isotopic labeling, and false if they have other
    differences. This also compares templates and families for TemplateReactions.

    The remove_isotope method can be slow, especially when comparing molecules
    and reactions. This was due to many copying of objects.

    This method avoid copying by storing the isotope and atom objects,
    removing them, doing the comparison, and rewriting them when
    finished the comparison.
    """

    atomlist = remove_isotope(obj1, inplace=True) + remove_isotope(obj2, inplace=True)
    if isinstance(obj1, Reaction):
        # make sure isotomorphic
        comparison_bool = obj1.is_isomorphic(obj2, eitherDirection)
        if comparison_bool and isinstance(obj1, TemplateReaction):
            # ensure families are the same
            comparison_bool = obj1.family == obj2.family
            if comparison_bool and not eitherDirection:
                # make sure templates are identical if in the same direction
                comparison_bool = frozenset(obj1.template) == frozenset(obj2.template)
    elif isinstance(obj1, Species):
        comparison_bool = obj1.is_isomorphic(obj2)
    else:
        raise TypeError('Only Reaction and Speicies Objects are supported in compareIsotopomers')
    redo_isotope(atomlist)
    return comparison_bool


def generate_RMG_model(inputFile, outputDirectory):
    """
    Generate the RMG-Py model NOT containing any non-normal isotopomers.

    Returns created RMG object.
    """
    initializeLog(logging.INFO, os.path.join(outputDirectory, 'RMG.log'))
    # generate mechanism:
    rmg = RMG(inputFile=os.path.abspath(inputFile),
              outputDirectory=os.path.abspath(outputDirectory))
    rmg.execute()

    return rmg


def correct_entropy(isotopomer, isotopeless):
    """
    Correct the entropy of the isotopomer by the following correction for symmetry:

    S(corrected) = S(original) + R*ln(sigma(isotopeless)) - R*ln(sigma(isotopomer))

    This method also copies the Enthalpy, Cp and other thermo parameters from isotopeless
    """

    # calculate -R ln (sigma) in SI units (J/K/mol)
    s_isotopeless = - constants.R * math.log(isotopeless.get_symmetry_number())
    s_isotopomer = - constants.R * math.log(isotopomer.get_symmetry_number())

    # convert species thermo to ThermoData object:
    nasa = isotopomer.thermo

    # apply correction to entropy at 298K
    delta_s = s_isotopomer - s_isotopeless
    nasa = nasa.changeBaseEntropy(delta_s)

    # put the corrected thermo back as a species attribute:
    isotopomer.thermo = nasa


def apply_kinetic_isotope_effect_simple(rxn_clusters, kinetics_database):
    """
    This method modifies reaction rates in place for the implementation of
    kinetic isotope effects using method 'simple', which is described in the
    publication: Goldman, MJ, Vandewiele, NM, Ono, S, Green WH [in preparation]

    input:
            rxn_clusters - list of list of reactions obtained from `cluster`
            kinetics_database - KineticsDatabase object for finding atom labels
    output:
            None (clusters are modified in place)

    This method has a number of dependent methods, which appear at the start of
    this method
    """

    # now for the start of applyKineticIsotopeEffectSimple
    for index, cluster in enumerate(rxn_clusters):
        # hardcoded family values determine what type of transition state
        # approximation is used.
        family = kinetics_database.families[cluster[0].family]
        if cluster[0].family.lower() == 'r_recombination':
            labels = ['*']
            three_member_ts = False
        elif cluster[0].family.lower() == 'r_addition_multiplebond':
            labels = ['*1', '*3']
            three_member_ts = False
        elif cluster[0].family.lower() == 'intra_r_add_endocyclic':
            labels = ['*1', '*3']
            three_member_ts = False
        elif cluster[0].family.lower() == 'intra_r_add_exocyclic':
            labels = ['*1', '*2']
            three_member_ts = False
        elif cluster[0].family.lower() == 'h_abstraction':
            labels = ['*1', '*3']
            three_member_ts = True
        elif cluster[0].family.lower() == 'intra_h_migration':
            labels = ['*1', '*2']
            three_member_ts = True
        elif cluster[0].family.lower() == 'disproportionation':
            labels = ['*1', '*2']
            three_member_ts = True
        else:
            logging.warning('isotope: kinetic isotope effect of family {0} not encoded into RMG. '
                            'Ignoring KIE of reaction {1}'.format(cluster[0].family, cluster[-1]))
            continue
        logging.debug('modifying reaction rate for cluster {0} for family {1}'.format(index, family.name))
        # get base reduced mass
        reaction = cluster[-1]  # set unlabeled reaction as the standard to compare
        labeled_reactants = get_labeled_reactants(reaction, family)
        base_reduced_mass = get_reduced_mass(labeled_reactants, labels, three_member_ts)
        for reaction in cluster[:-1]:
            labeled_reactants = get_labeled_reactants(reaction, family)
            reduced_mass = get_reduced_mass(labeled_reactants, labels, three_member_ts)
            reaction.kinetics.change_rate(math.sqrt(base_reduced_mass / reduced_mass))


def get_labeled_reactants(reaction, family):
    """
    Returns a list of labeled molecule objects given a species-based reaction object and it's
    corresponding kinetic family.

    Used for KIE method 'simple'
    """
    if reaction.family != family.name:
        raise AttributeError(
            "The reaction must come from the family specified: {0} != {1}".format(reaction.family, family.name))
    # save the reactants and products to replace in the reaction object
    reactants = list(reaction.reactants)
    products = list(reaction.products)

    family.add_atom_labels_for_reaction(reaction, output_with_resonance=True)
    labeled_reactants = [species.molecule[0] for species in reaction.reactants]

    # replace the original reactants and products
    reaction.reactants = reactants
    reaction.products = products

    return labeled_reactants


def get_reduced_mass(labeled_molecules, labels, three_member_ts):
    """
    Returns the reduced mass of the labeled elements
    within the labeled molecules. Used for kinetic isotope effect.
    The equations are based on Melander & Saunders, Reaction rate of isotopic molecules, 1980

    input: labeled_molecules - list of molecules with labels for the reaction
           labels - list of strings to search for labels in the product
           three_member_ts - boolean describing number of atoms in transition state
                             If True, reactions involving the movement of hydrogen between
                             two carbons is used. The labels should not include the hydrogen atom
                             If False, only a 2 member transition state is used.
    output: float

    Used for KIE method 'simple'
    """
    reduced_mass = 0.
    combined_mass = 0.
    for labeled_mol in labeled_molecules:
        for atom in labeled_mol.atoms:
            if any([atom.label == label for label in labels]):
                if three_member_ts:
                    combined_mass += atom.element.mass
                else:
                    reduced_mass += 1. / atom.element.mass
    if reduced_mass == 0. and combined_mass == 0:
        from rmgpy.exceptions import KineticsError
        raise KineticsError(
            "Did not find a labeled atom in molecules {}".format([mol.to_adjacency_list() for mol in labeled_molecules]))
    if three_member_ts:  # actually convert to reduced mass using the mass of hydrogen
        reduced_mass = 1 / rmgpy.molecule.element.H.mass + 1 / combined_mass
    return 1. / reduced_mass


def is_enriched(obj):
    """
    Returns True if the species or reaction object has any enriched isotopes.
    """

    if isinstance(obj, Species):
        for atom in obj.molecule[0].atoms:
            if atom.element.isotope != -1 and not np.allclose(atom.element.mass, get_element(atom.element.symbol).mass):
                return True
        return False
    elif isinstance(obj, Reaction):
        enriched = []
        for spec in obj.reactants:
            enriched.append(is_enriched(spec))
        for spec in obj.products:
            enriched.append(is_enriched(spec))
        return any(enriched)
    else:
        raise TypeError('is_enriched only takes species and reaction objects. {} was sent'.format(str(type(obj))))


def ensure_correct_symmetry(isotopmoper_list, isotopic_element='C'):
    """
    given a list of isotopomers (species' objects) and the element that is labeled,
    returns True if the correct symmetry is detected. False if not detected

    This uses the observation that if there are less than 2^n isotopomors, where n is the number of
    isotopic elements, then there should be an 2^n - m isotopomers where
    m is related to the amount of entropy increased by some isotopomers since they
    lost symmetry.
    """
    number_elements = 0
    for atom in isotopmoper_list[0].molecule[0].atoms:
        if atom.element.symbol == isotopic_element:
            number_elements += 1

    minimum_entropy = min([spec.getEntropy(298) for spec in isotopmoper_list])

    count = 0.
    for spec in isotopmoper_list:
        entropy_diff = spec.getEntropy(298) - minimum_entropy
        count += math.exp(entropy_diff / constants.R)
    return abs(count - 2 ** number_elements) < 0.01


def ensure_correct_degeneracies(reaction_isotopomer_list, print_data=False, r_tol_small_flux=1e-5,
                                r_tol_deviation=0.0001):
    """
    given a list of isotopomers (reaction objects), returns True if the correct
    degeneracy values are detected. False if incorrect degeneracy values
    exist.

    This method assuumes an equilibrium distribution of compounds and finds
    the fluxes created by the set of reactions (both forward and
    reverse). It then checks that the fluxes of each compounds are proportional to their
    symmetry numbers (since lower symmetry numbers have higher entropy and should have a
    larger concentration.)

    inputs:

    reaction_isotopomer_list - a list of reactions that differ based on isotope placement
    print_data - output the table of fluxes obtained. useful for debugging incorrect values
    r_tol_small_flux - fraction of maximum flux to count as 'zero'
    r_tol_deviation - allowable numerical deviation in answers

    This method has a few datastructures it utilizes:

    product_structures - a list of tuples, (index, isotopomer_structure), containing reactants and products
    product_list - a pandas.DataFrame that sotres the fluxes and symmetry values
    """

    def store_flux_info(species, flux, product_list, product_structures):
        """
        input:
            species - The desired species that you'd like to modify the flux value of
            flux - the amount to modify the flux of the species
            product_list - a list of current flux and symmetry values
            product_structures - a list of unlabled structures un the reaction class

        modifies the product list by either adding on new row with the species' structure, flux
        symmetry ratio, etc. or adds the flux to the species' current row in `product_list`

        returns: modified product_list
        """
        # find the corresponding unlabeled structure
        for struc_index, product_structure in enumerate(product_structures):
            if compare_isotopomers(product_structure, species):
                structure_index = struc_index
                break
        # store product flux and symmetry info
        for index in product_list.index:
            spec = product_list.at[index, 'product']
            # if species already listed, add to its flux
            if product_list.at[index, 'product_struc_index'] == structure_index \
                    and spec.is_isomorphic(species):
                product_list.at[index, 'flux'] += flux
                return product_list
        # add product to list
        symmetry_ratio = product_structures[structure_index].get_symmetry_number() / float(species.get_symmetry_number())
        return product_list.append({'product': species,
                                    'flux': flux,
                                    'product_struc_index': structure_index,
                                    'symmetry_ratio': symmetry_ratio},
                                   ignore_index=True)

    # copy list in case it is called upon
    reaction_isotopomer_list = copy(reaction_isotopomer_list)
    product_list = pd.DataFrame(columns=['product', 'flux', 'product_struc_index', 'symmetry_ratio'])
    product_structures = []
    unlabeled_rxn = None

    # check if first reaction is unlabeled. If not, make unlabeled reaction
    if not is_enriched(reaction_isotopomer_list[0]):
        unlabeled_rxn = reaction_isotopomer_list[0]
    else:
        unlabeled_rxn = reaction_isotopomer_list[0].copy()
        for mol in unlabeled_rxn.reactants:
            remove_isotope(mol, inplace=True)
        for mol in unlabeled_rxn.products:
            remove_isotope(mol, inplace=True)
    unlabeled_symmetry_reactants = np.prod([mol.get_symmetry_number() for mol in unlabeled_rxn.reactants])
    unlabeled_symmetry_products = np.prod([mol.get_symmetry_number() for mol in unlabeled_rxn.products])

    # prepare index of structures (product_structures)
    for struc in unlabeled_rxn.reactants + unlabeled_rxn.products:
        new_structure = struc.copy(deep=True)
        product_structures.append(new_structure)

    # go through each rxn cataloging fluxes to each product
    for rxn in reaction_isotopomer_list:
        # find characteristic flux for the forward direction
        reactant_conc = 1.
        reactant_conc /= np.prod([mol.get_symmetry_number() for mol in rxn.reactants])
        reactant_conc *= unlabeled_symmetry_reactants
        if isinstance(rxn.kinetics, MultiArrhenius):
            rate = sum([arr.A.value_si for arr in rxn.kinetics.arrhenius])
        else:
            rate = rxn.kinetics.A.value_si
        product_flux = reactant_conc * rate

        # modify fluxes for forward direction
        for rxn_product in rxn.products:
            product_list = store_flux_info(rxn_product, product_flux, product_list, product_structures)
        for rxn_reactant in rxn.reactants:
            product_list = store_flux_info(rxn_reactant, -product_flux, product_list, product_structures)

        # now find characteristic flux of reverse direction.
        if isinstance(rxn.kinetics, MultiArrhenius):
            reverse_A_factor = 0
            for arr in rxn.kinetics.arrhenius:
                reverse_A_factor += arr.A.value_si / rxn.getEquilibriumConstant(298)
        else:
            reverse_A_factor = rxn.kinetics.A.value_si / rxn.getEquilibriumConstant(298)

        # get reverse flux using product symmetries
        product_conc = 1.
        product_conc /= np.prod([mol.get_symmetry_number() for mol in rxn.products])
        product_conc *= unlabeled_symmetry_products
        reactant_flux = product_conc * reverse_A_factor

        # modify reverse fluxes
        for rxn_product in rxn.products:
            product_list = store_flux_info(rxn_product, -reactant_flux, product_list, product_structures)
        for rxn_reactant in rxn.reactants:
            product_list = store_flux_info(rxn_reactant, reactant_flux, product_list, product_structures)

    if print_data:
        print(product_list.sort_values(['product_struc_index', 'symmetry_ratio']))

    # now ensure the fluxes are correct or cancel out & throw error if not.
    pass_species = []
    for index in range(len(product_structures)):
        products = product_list[product_list.product_struc_index == index]
        fluxes = products.flux / products.flux.sum()
        symmetries = products.symmetry_ratio / products.symmetry_ratio.sum()
        # the two decisison criteria
        accurate = np.allclose(fluxes, symmetries, rtol=r_tol_deviation)
        low_fluxes = all(products.flux.abs() < max(reactant_flux, product_flux) * r_tol_small_flux)
        pass_species.append(low_fluxes or accurate)
    return all(pass_species)


def run(inputFile, outputDir, original=None, maximumIsotopicAtoms=1,
        useOriginalReactions=False,
        kineticIsotopeEffect=None):
    """
    Accepts one input file with the RMG-Py model to generate.

    Firstly, generates the RMG model for the first input file. Takes the core species of that mechanism
    and generates all isotopomers of those core species. Next, generates all reactions between the
    generated pool of isotopomers, and writes it to file.
    """
    logging.info("isotope: Starting the RMG isotope generation method 'run'")
    if not original:
        logging.info("isotope: original model not found, generating new one in directory `rmg`")
        logging.info("isotope: check `rmg/RMG.log` for the rest of the logging info.")

        outputdir_rmg = os.path.join(outputDir, 'rmg')
        os.mkdir(outputdir_rmg)

        rmg = generate_RMG_model(inputFile, outputdir_rmg)
    else:
        logging.info("isotope: original model being copied from previous RMG job in folder {}".format(original))
        outputdir_rmg = original
        chemkin_file = os.path.join(outputdir_rmg, 'chemkin', 'chem_annotated.inp')
        dict_file = os.path.join(outputdir_rmg, 'chemkin', 'species_dictionary.txt')
        rmg = loadRMGJob(inputFile, chemkin_file, dict_file, generateImages=False, useChemkinNames=True)

    logging.info("isotope: generating isotope model")
    logging.info('Generating isotopomers for the core species in {}'.format(outputdir_rmg))
    isotopes = []

    logging.info("isotope: adding all the new and old isotopomers")
    for spc in rmg.reactionModel.core.species:
        find_cp0_and_cpinf(spc, spc.thermo)
        isotopes.append([spc] + generate_isotopomers(spc, maximumIsotopicAtoms))

    logging.info('isotope: number of isotopomers: {}'.format(
        sum([len(isotopomer) for isotopomer in isotopes if isotopomer])))

    outputdir_iso = os.path.join(outputDir, 'iso')
    os.mkdir(outputdir_iso)

    logging.info('isotope: Generating RMG isotope model in {}'.format(outputdir_iso))
    generate_isotope_model(outputdir_iso, rmg, isotopes, useOriginalReactions=useOriginalReactions,
                           kineticIsotopeEffect=kineticIsotopeEffect)
