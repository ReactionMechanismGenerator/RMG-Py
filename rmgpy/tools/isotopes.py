#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import os
import os.path
import logging
import numpy as np
from copy import copy, deepcopy
import pandas as pd
import math

import rmgpy.constants as constants
from rmgpy.molecule import Molecule
from rmgpy.molecule.element import getElement
from rmgpy.tools.loader import loadRMGJob
from rmgpy.chemkin import ChemkinWriter
from rmgpy.rmg.main import RMG, initializeLog
from rmgpy.rmg.model import Species
from rmgpy.species import Species as Species2
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.thermo.thermoengine import processThermoData
from rmgpy.data.thermo import findCp0andCpInf
from rmgpy.data.rmg import getDB
import rmgpy.molecule.element
from rmgpy.kinetics.arrhenius import MultiArrhenius
from rmgpy.reaction import isomorphic_species_lists


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


def generate_isotope_model(outputDirectory, rmg0, isotopes, useOriginalReactions = False,
                         kineticIsotopeEffect = None):
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
        rmg.reactionModel.processNewReactions(rxns,newSpecies=[])

    else:
        logging.info("isotope: enlarging the isotope model")
        rmg.reactionModel.enlarge(reactEdge=True,
            unimolecularReact=rmg.unimolecularReact,
            bimolecularReact=rmg.bimolecularReact)

    logging.info("isotope: clustering reactions")
    clusters = cluster(rmg.reactionModel.core.reactions)
    logging.info('isotope: fixing the directions of every reaction to a standard')
    for isotopomerRxnList in clusters:
        ensureReactionDirection(isotopomerRxnList)

    consistent = True
    logging.info("isotope: checking symmetry is consistent among isotopomers")
    for species_list in cluster(rmg.reactionModel.core.species):
        if not ensure_correct_symmetry(species_list):
            logging.info("isotopomers of {} with index {} may have wrong symmetry".format(species_list[0], species_list[0].index))
            consistent = False
    logging.info("isotope: checking that reaction degeneracy is consistent among isotopomers")
    for rxn_list in clusters:
        if not ensure_correct_degeneracies(rxn_list):
            logging.info("isotopomers of {} with index {} may have incorrect degeneracy.".format(rxn_list[0], rxn_list[0].index))
            consistent = False
    if not consistent:
        logging.warning("isotope: non-consistent degeneracy and/or symmetry was detected. This may lead to unrealistic deviations in enrichment. check log for more details")

    if kineticIsotopeEffect:
        logging.info('isotope: modifying reaction rates using kinetic isotope effect method "{0}"'.format(kineticIsotopeEffect))
        if kineticIsotopeEffect == 'simple':
            apply_kinetic_isotope_effect_simple(clusters,rmg.database.kinetics)
        else:
            logging.warning('isotope: kinetic isotope effect {0} is not supported. skipping adding kinetic isotope effects.')
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
        assert isinstance(rxn,TemplateReaction)
        assert rxn.template is not None, 'isotope reaction {0} does not have a template attribute. Full details :\n\n{1}'.format(str(rxn),repr(rxn))

    from rmgpy.reaction import _isomorphicSpeciesList

    found_reactions = []
    rxn_index = 0
    while rxn_index < len(isotopeless_reactions):
        rxn = isotopeless_reactions[rxn_index]
        # find all reactions involving same reactants
        rxns_w_same_reactants = [rxn]
        rxn_index2 = rxn_index + 1
        while rxn_index2 < len(isotopeless_reactions):
            if isomorphic_species_lists(isotopeless_reactions[rxn_index].reactants,
                                      isotopeless_reactions[rxn_index2].reactants,
                                     ):
                rxns_w_same_reactants.append(isotopeless_reactions[rxn_index2])
                del isotopeless_reactions[rxn_index2]
            else:
                rxn_index2 += 1
        ##### find all pairs of reacting isotoper species #####
        # find the lists of reactants that have identical isotopomers
        reactants = []
        for reactant in rxn.reactants:
            for iso_index, isotopomers in enumerate(isotopes):
                if compare_isotopomers(reactant,isotopomers[0]):
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
                if isomorphic_species_lists(reactant_pairs[rxn_index3],
                                          reactant_pairs[rxn_index4]):
                    del reactant_pairs[rxn_index4]
                else:
                    rxn_index4 += 1
            rxn_index3 += 1

        # make reaction objects
        for pair in reactant_pairs:
            # copy species so they don't get modified
            speciesTuple = tuple([spc.copy(deep=True) for spc in pair])
            unfiltered_rxns = getDB('kinetics').generate_reactions_from_families(speciesTuple)
            # remove reactions whose products don't match the original reactions
            rxn_index5 = 0
            while rxn_index5 < len(unfiltered_rxns):
                for isotopeless_reaction in rxns_w_same_reactants:
                    isotopeless_kinetics = isotopeless_reaction.kinetics
                    isotopeless_degeneracy = isotopeless_reaction.degeneracy
                    if compare_isotopomers(isotopeless_reaction, unfiltered_rxns[rxn_index5],eitherDirection = False)\
                        and isotopeless_reaction.family == unfiltered_rxns[rxn_index5].family\
                        and frozenset(isotopeless_reaction.template) == \
                                     frozenset(unfiltered_rxns[rxn_index5].template):
                        # apply kinetics to new reaction & modify for degeneracy
                        unfiltered_rxns[rxn_index5].kinetics = deepcopy(isotopeless_kinetics)
                        unfiltered_rxns[rxn_index5].kinetics.changeRate(unfiltered_rxns[rxn_index5].degeneracy / isotopeless_degeneracy)
                        rxn_index5 += 1
                        break
                else: # did not find same prodcuts
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
    isotope = getElement(6, 13)

    mols = []
    add_isotope(0, N, mol, mols, isotope)

    spcs = []
    for isomol in mols:
        isotopomer = Species(molecule=[isomol], thermo=deepcopy(spc.thermo), transportData=spc.transportData, reactive=spc.reactive)
        isotopomer.generate_resonance_structures(keep_isomorphic=True)
        spcs.append(isotopomer)

    # do not retain identical species:
    filtered = []
    while spcs:
        candidate = spcs.pop()
        unique = True
        for isotopomer in filtered:
            if isotopomer.isIsomorphic(candidate):
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
    if i == N: return
    else:
        atoms = filter(lambda at: at.symbol == element.symbol, mol.atoms)
        for at in atoms:
            if at.element == element: continue
            else:
                isotopomer = mol.copy(deep=True)
                isotopomer.atoms[mol.atoms.index(at)].element = element
                mols.append(isotopomer)
                add_isotope(i+1, N, isotopomer, mols, element)

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
            if compare_isotopomers(cluster[0],candidate):
                cluster.append(candidate)
                break
        else:
            clusters.append([candidate])

    return clusters

def remove_isotope(labeledObj, inplace = False):
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

    if isinstance(labeledObj,Species2):
        if inplace:
            modifiedAtoms = []
            for mol in labeledObj.molecule:
                for atom in mol.atoms:
                    if atom.element.isotope != -1:
                        modifiedAtoms.append((atom,atom.element))
                        atom.element = getElement(atom.element.symbol)
            return modifiedAtoms
        else:
            stripped = labeledObj.copy(deep=True)
    
            for atom in stripped.molecule[0].atoms:
                if atom.element.isotope != -1:
                    atom.element = getElement(atom.element.symbol)

        # only do it for the first molecule, generate the other resonance isomers.
            stripped.molecule = [stripped.molecule[0]]
            stripped.generate_resonance_structures(keep_isomorphic=True)

        return stripped

    elif isinstance(labeledObj,Reaction):

        if inplace:

            atomList = []
            for reactant in  labeledObj.reactants:
                removed = remove_isotope(reactant,inplace)
                if removed:
                    atomList += removed
            for product in labeledObj.products:
                removed = remove_isotope(product,inplace)
                if removed:
                    atomList += removed

            return atomList
        else:
            strippedRxn = labeledObj.copy()

            strippedReactants = []
            for reactant in  strippedRxn.reactants:
                strippedReactants.append(remove_isotope(reactant,inplace))
            strippedRxn.reactants = strippedReactants

            strippedProducts = []
            for product in  strippedRxn.products:
                strippedProducts.append(remove_isotope(product,inplace))
            strippedRxn.products = strippedProducts

            return strippedRxn
    elif isinstance(labeledObj,Molecule):
        if inplace:
            modifiedAtoms = []
            for atom in labeledObj.atoms:
                if atom.element.isotope != -1:
                    modifiedAtoms.append((atom,atom.element))
                    atom.element = getElement(atom.element.symbol)
            return modifiedAtoms
        else:
            stripped = labeledObj.copy(deep=True)

            for atom in stripped.atoms:
                if atom.element.isotope != -1:
                    atom.element = getElement(atom.element.symbol)

            return stripped
    else:
        raise TypeError('Only Reaction, Species, and Molecule objects are supported')


def redo_isotope(atomList):
    """
    This takes a list of zipped atoms with their isotopes removed, from 
    and elements.
    """
    for atom, element in atomList:
        atom.element = element

def compare_isotopomers(obj1, obj2, eitherDirection = True):
    """
    This method takes two species or reaction objects and returns true if
    they only differ in isotopic labeling, and false if they have other
    differences.

    The remove_isotope method can be slow, especially when comparing molecules
    and reactions. This was due to many copying of objects.

    This method avoid copying by storing the isotope and atom objects,
    removing them, doing the comparison, and rewriting them when
    finished the comparison.
    """

    atomlist = remove_isotope(obj1,inplace=True) + remove_isotope(obj2,inplace=True)
    if isinstance(obj1,Reaction):
        comparisonBool = obj1.isIsomorphic(obj2, eitherDirection)
        if comparisonBool and isinstance(obj1, TemplateReaction):
            # ensure families are the same
            comparisonBool = obj1.family == obj2.family
            if comparisonBool and not eitherDirection:
                # make sure templates are identical if in the same direction
                comparisonBool = frozenset(obj1.template) == frozenset(obj2.template)
    elif isinstance(obj1,Species2):
        comparisonBool = obj1.isIsomorphic(obj2)
    else:
        raise TypeError('Only Reaction and Speicies Objects are supported in compareIsotopomers')
    redo_isotope(atomlist)
    return comparisonBool

def generate_RMG_model(inputFile, outputDirectory):
    """
    Generate the RMG-Py model NOT containing any non-normal isotopomers.

    Returns created RMG object.
    """
    initializeLog(logging.INFO, os.path.join(outputDirectory, 'RMG.log'))
    # generate mechanism:
    rmg = RMG(inputFile = os.path.abspath(inputFile),
            outputDirectory = os.path.abspath(outputDirectory)
        )
    rmg.execute()

    return rmg

def is_enriched(obj):
    """
    Returns True if the species or reaction object has any enriched isotopes.
    """

    if isinstance(obj,Species):
        for atom in obj.molecule[0].atoms:
            if atom.element.isotope != -1 and not np.allclose(atom.element.mass, getElement(atom.element.symbol).mass):
                return True
        return False
    elif isinstance(obj,Reaction):
        enriched = []
        for spec in obj.reactants:
            enriched.append(is_enriched(spec))
        for spec in obj.products:
            enriched.append(is_enriched(spec))
        return any(enriched)
    else:
        raise TypeError('is_enriched only takes species and reaction objects. {} was sent'.format(str(type(obj))))

def run(inputFile, outputDir, original=None, maximumIsotopicAtoms = 1,
                            useOriginalReactions = False,
                            kineticIsotopeEffect = None):
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

        outputdirRMG = os.path.join(outputDir, 'rmg')
        os.mkdir(outputdirRMG)

        rmg = generate_RMG_model(inputFile, outputdirRMG)
    else:
        logging.info("isotope: original model being copied from previous RMG job in folder {}".format(original))
        outputdirRMG = original
        chemkinFile = os.path.join(outputdirRMG, 'chemkin', 'chem_annotated.inp')
        dictFile = os.path.join(outputdirRMG, 'chemkin', 'species_dictionary.txt')
        rmg = loadRMGJob(inputFile, chemkinFile, dictFile, generateImages=False, useChemkinNames=True)

    logging.info("isotope: generating isotope model")
    logging.info('Generating isotopomers for the core species in {}'.format(outputdirRMG))
    isotopes = []

    logging.info("isotope: adding all the new and old isotopomers")
    for spc in rmg.reactionModel.core.species:
        findCp0andCpInf(spc, spc.thermo)
        isotopes.append([spc] + generate_isotopomers(spc, maximumIsotopicAtoms))
    logging.info('isotope: number of isotopomers: {}'.format(sum([len(isotopomer) for isotopomer in isotopes if isotopomer])))

    outputdirIso = os.path.join(outputDir, 'iso')
    os.mkdir(outputdirIso)

    logging.info('isotope: Generating RMG isotope model in {}'.format(outputdirIso))
    generate_isotope_model(outputdirIso, rmg, isotopes, useOriginalReactions = useOriginalReactions, kineticIsotopeEffect = kineticIsotopeEffect)
