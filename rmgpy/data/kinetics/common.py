#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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
This module contains classes and functions that are used by multiple modules
in this subpackage.
"""
import itertools
import logging

from rmgpy.data.base import LogicNode
from rmgpy.reaction import Reaction
from rmgpy.molecule import Group, Molecule
from rmgpy.species import Species
from rmgpy.exceptions import DatabaseError, KineticsError

################################################################################

def saveEntry(f, entry):
    """
    Save an `entry` in the kinetics database by writing a string to
    the given file object `f`.
    """
    from rmgpy.cantherm.output import prettify

    def sortEfficiencies(efficiencies0):
        efficiencies = {}
        for mol, eff in efficiencies0.iteritems():
            if isinstance(mol, str):
                # already in SMILES string format
                smiles = mol
            else:
                smiles = mol.toSMILES()
                
            efficiencies[smiles] = eff
        keys = efficiencies.keys()
        keys.sort()
        return [(key, efficiencies[key]) for key in keys]

    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    if entry.label != '':
        f.write('    label = "{0}",\n'.format(entry.label))


    #Entries for kinetic rules, libraries, training reactions
    #and depositories will have an Reaction object for its item
    if isinstance(entry.item, Reaction):
        #Write out additional data if depository or library
        #kinetic rules would have a Group object for its reactants instead of Species
        if isinstance(entry.item.reactants[0], Species):
            # Add degeneracy if the reaction is coming from a depository or kinetics library
            f.write('    degeneracy = {0:.1f},\n'.format(entry.item.degeneracy))
            if entry.item.duplicate:
                f.write('    duplicate = {0!r},\n'.format(entry.item.duplicate))
            if not entry.item.reversible:
                f.write('    reversible = {0!r},\n'.format(entry.item.reversible))
    #Entries for groups with have a group or logicNode for its item
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList())
        f.write('""",\n')
    elif isinstance(entry.item, LogicNode):
        f.write('    group = "{0}",\n'.format(entry.item))
    else:
        raise DatabaseError("Encountered unexpected item of type {0} while saving database.".format(entry.item.__class__))

    # Write kinetics
    if isinstance(entry.data, str):
        f.write('    kinetics = "{0}",\n'.format(entry.data))
    elif entry.data is not None:
        efficiencies = None
        if hasattr(entry.data, 'efficiencies'):
            efficiencies = entry.data.efficiencies
            entry.data.efficiencies = dict(sortEfficiencies(entry.data.efficiencies))
        kinetics = prettify(repr(entry.data))
        kinetics = '    kinetics = {0},\n'.format(kinetics.replace('\n', '\n    '))
        f.write(kinetics)
        if hasattr(entry.data, 'efficiencies'):
            entry.data.efficiencies = efficiencies
    else:
        f.write('    kinetics = None,\n')
            
    # Write reference
    if entry.reference is not None:
        reference = entry.reference.toPrettyRepr()
        lines = reference.splitlines()
        f.write('    reference = {0}\n'.format(lines[0]))
        for line in lines[1:-1]:
            f.write('    {0}\n'.format(line))
        f.write('    ),\n'.format(lines[0]))
    
    if entry.referenceType != "":
        f.write('    referenceType = "{0}",\n'.format(entry.referenceType))
    if entry.rank is not None:
        f.write('    rank = {0},\n'.format(entry.rank))
        
    if entry.shortDesc.strip() !='':
        f.write('    shortDesc = u"""')
        try:
            f.write(entry.shortDesc.encode('utf-8'))
        except:
            f.write(entry.shortDesc.strip().encode('ascii', 'ignore')+ "\n")
        f.write('""",\n')
    
    if entry.longDesc.strip() !='':
        f.write('    longDesc = \n')
        f.write('u"""\n')
        try:
            f.write(entry.longDesc.strip().encode('utf-8') + "\n")
        except:
            f.write(entry.longDesc.strip().encode('ascii', 'ignore')+ "\n")
        f.write('""",\n')

    f.write(')\n\n')


def filterReactions(reactants, products, reactionList):
    """
    Remove any reactions from the given `reactionList` whose reactants do
    not involve all the given `reactants` or whose products do not involve 
    all the given `products`. This method checks both forward and reverse
    directions, and only filters out reactions that don't match either.
    
    reactants and products can be either molecule or species objects
    """
    
    # Convert from molecules to species and generate resonance isomers.
    reactants = ensure_species(reactants, resonance=True)
    products = ensure_species(products, resonance=True)
    
    reactions = reactionList[:]
    
    for reaction in reactionList:
        # Forward direction
        reactants0 = [r for r in reaction.reactants]
        for reactant in reactants:
            for reactant0 in reactants0:
                if reactant.isIsomorphic(reactant0):
                    reactants0.remove(reactant0)
                    break
        products0 = [p for p in reaction.products]
        for product in products:
            for product0 in products0:
                if product.isIsomorphic(product0):
                    products0.remove(product0)
                    break
        forward = not (len(reactants0) != 0 or len(products0) != 0)
        # Reverse direction
        reactants0 = [r for r in reaction.products]
        for reactant in reactants:
            for reactant0 in reactants0:
                if reactant.isIsomorphic(reactant0):
                    reactants0.remove(reactant0)
                    break
        products0 = [p for p in reaction.reactants]
        for product in products:
            for product0 in products0:
                if product.isIsomorphic(product0):
                    products0.remove(product0)
                    break
        reverse = not (len(reactants0) != 0 or len(products0) != 0)
        if not forward and not reverse:
            reactions.remove(reaction)
    return reactions


def ensure_species(input_list, resonance=False, keepIsomorphic=False):
    """
    Given an input list of molecules or species, return a list with only
    species objects.
    """
    output_list = []
    for item in input_list:
        if isinstance(item, Molecule):
            new_item = Species(molecule=[item])
        elif isinstance(item, Species):
            new_item = item
        else:
            raise TypeError('Only Molecule or Species objects can be handled.')
        if resonance:
            new_item.generateResonanceIsomers(keepIsomorphic=keepIsomorphic)
        output_list.append(new_item)

    return output_list


def generate_molecule_combos(input_species):
    """
    Generate combinations of molecules from the given species objects.
    """
    if len(input_species) == 1:
        combos = [(mol,) for mol in input_species[0].molecule]
    elif len(input_species) == 2:
        combos = itertools.product(input_species[0].molecule, input_species[1].molecule)
    else:
        raise ValueError('Reaction generation can be done for 1 or 2 species, not {0}.'.format(len(input_species)))

    return combos


def _labelListOfSpecies(speciesTuple):
    """
    given a list or tuple of species' objects, ensure all their atoms' id are
    independent.

    Modifies the speciesTuple in place, nothing returned.
    """
    # assert that all species' atomlabels are different
    def independentIDs():
        num_atoms = 0
        IDs = []
        for species in speciesTuple:
            num_atoms += len(species.molecule[0].atoms)
            IDs.extend([atom.id for atom in species.molecule[0].atoms ])
        num_ID = len(set(IDs))
        return num_ID == num_atoms
    # if they are different, relabel and remake atomIDs
    if not independentIDs():
        logging.debug('identical atom ids found between species. regenerating')
        for species in speciesTuple:
            mol = species.molecule[0]
            mol.assignAtomIDs()
            # remake resonance isomers with new labeles
            species.molecule = [mol]
            species.generateResonanceIsomers(keepIsomorphic = True)


def findDegeneracies(rxnList, sameReactants=None):
    """
    given a list of Reaction object with Molecule objects, this method
    removes degenerate reactions and increments the degeneracy of the
    reaction object. For multiple transition states, this method adds
    them as separate duplicate reactions. This method modifies
    rxnList in place and does not return anything.

    This algorithm used to exist in family.__generateReactions, but was moved
    here because it didn't have any family dependence.
    """

    # We want to sort all the reactions into sublists composed of isomorphic reactions
    # with degenerate transition states
    rxnSorted = []
    for rxn0 in rxnList:
        # find resonance structures for rxn0
        convertToSpeciesObjects(rxn0)
        if len(rxnSorted) == 0:
            # This is the first reaction, so create a new sublist
            rxnSorted.append([rxn0])
        else:
            # Loop through each sublist, which represents a unique reaction
            for rxnList1 in rxnSorted:
                # Try to determine if the current rxn0 is identical or isomorphic to any reactions in the sublist
                isomorphic = False
                identical = False
                sameTemplate = False
                for rxn in rxnList1:
                    isomorphic = rxn0.isIsomorphic(rxn, checkIdentical=False, checkTemplateRxnProducts=True)
                    if not isomorphic:
                        identical = False
                    else:
                        identical = rxn0.isIsomorphic(rxn, checkIdentical=True, checkTemplateRxnProducts=True)
                    sameTemplate = frozenset(rxn.template) == frozenset(rxn0.template)
                    if not isomorphic:
                        # a different product was found, go to next list
                        break
                    elif not sameTemplate:
                        # a different transition state was found, mark as duplicate and
                        # go to the next sublist
                        rxn.duplicate = True
                        rxn0.duplicate = True
                        break
                    elif identical:
                        # An exact copy of rxn0 is already in our list, so we can move on to the next rxn
                        break
                    else: # sameTemplate and isomorphic but not identical
                        # This is the right sublist for rxn0, but continue to see if there is an identical rxn
                        continue
                else:
                    # We did not break, so this is the right sublist, but there is no identical reaction
                    # This means that we should add rxn0 to the sublist as a degenerate rxn
                    rxnList1.append(rxn0)
                if isomorphic and sameTemplate:
                    # We already found the right sublist, so we can move on to the next rxn
                    break
            else:
                # We did not break, which means that there was no isomorphic sublist, so create a new one
                rxnSorted.append([rxn0])

    rxnList = []
    for rxnList1 in rxnSorted:
        # Collapse our sorted reaction list by taking one reaction from each sublist
        rxn = rxnList1[0]
        # The degeneracy of each reaction is the number of reactions that were in the sublist
        rxn.degeneracy = sum([reaction0.degeneracy for reaction0 in rxnList1])
        rxnList.append(rxn)

    for rxn in rxnList:
        if rxn.isForward:
            reduceSameReactantDegeneracy(rxn, sameReactants)
        else:
            # fix the degeneracy of (not ownReverse) reactions found in the backwards direction
            correctDegeneracyOfReverseReaction(rxn)

    return rxnList


def convertToSpeciesObjects(reaction):
    """
    modifies a reaction holding Molecule objects to a reaction holding
    Species objects, with generated resonance isomers.
    """
    # if already species' objects, return none
    if isinstance(reaction.reactants[0],Species):
        return None
    # obtain species with all resonance isomers
    for i, mol in enumerate(reaction.reactants):
        spec = Species(molecule = [mol])
        if not reaction.isForward:
            spec.generateResonanceIsomers(keepIsomorphic=True)
        reaction.reactants[i] = spec
    for i, mol in enumerate(reaction.products):
        spec = Species(molecule = [mol])
        if reaction.isForward:
            spec.generateResonanceIsomers(keepIsomorphic=True)
        reaction.products[i] = spec

    # convert reaction.pairs object to species
    newPairs=[]
    for reactant, product in reaction.pairs:
        newPair = []
        for reactant0 in reaction.reactants:
            if reactant0.isIsomorphic(reactant):
                newPair.append(reactant0)
                break
        for product0 in reaction.products:
            if product0.isIsomorphic(product):
                newPair.append(product0)
                break
        newPairs.append(newPair)
    reaction.pairs = newPairs

    try:
        convertToSpeciesObjects(reaction.reverse)
    except AttributeError:
        pass


def reduceSameReactantDegeneracy(reaction, sameReactants=None):
    """
    This method reduces the degeneracy of reactions with identical reactants,
    since translational component of the transition states are already taken
    into account (so swapping the same reactant is not valid)

    This comes from work by Bishop and Laidler in 1965
    """
    if len(reaction.reactants) == 2 and (
                (reaction.isForward and sameReactants) or
                reaction.reactants[0].isIsomorphic(reaction.reactants[1])
            ):
        reaction.degeneracy *= 0.5
        logging.debug('Degeneracy of reaction {} was decreased by 50% to {} since the reactants are identical'.format(reaction, reaction.degeneracy))


def correctDegeneracyOfReverseReaction(reaction):
    """
    This method corrects the degeneracy of reactions found when the backwards
    template is used. Given the following parameters:

        reaction - list of reactions with their degeneracies already counted

    This method modifies reaction in place and returns nothing

    This does not adjust for identical reactants, you need to use `reduceSameReactantDegeneracy`
    to adjust for that.
    """
    from rmgpy.data.rmg import getDB
    family = getDB('kinetics').families[reaction.family]
    if not family.ownReverse:
        reaction.degeneracy = family.calculateDegeneracy(reaction)
