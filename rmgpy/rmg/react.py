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
Contains functions for generating reactions.
"""
import logging
import itertools

from rmgpy.data.rmg import getDB
from rmgpy.scoop_framework.util import map_
from rmgpy.species import Species
        
def react(*spcTuples):
    """
    Generate reactions between the species in the 
    list of species tuples for all the reaction families available.

    For each tuple of one or more Species objects [(spc1,), (spc2, spc3), ...]
    the following is done:

    A list of tuples is created for each resonance isomer of the species.
    Each tuple consists of (Molecule, index) with the index the species index of the Species object.

    Possible combinations between the first spc in the tuple, and the second species in the tuple
    is obtained by taking the combinatorial product of the two generated [(Molecule, index)] lists.

    Returns a flat generator object containing the generated Reaction objects.
    """

    results = map_(
                reactSpecies,
                spcTuples)

    reactions = itertools.chain.from_iterable(results)

    return reactions

def reactSpecies(speciesTuple):
    """
    given one species tuple, will find the reactions and remove degeneracy
    from them.
    """
    speciesTuple = tuple([spc.copy(deep=True) for spc in speciesTuple])

    _labelListOfSpecies(speciesTuple)

    combos = getMoleculeTuples(speciesTuple)

    reactions = map(reactMolecules,combos)
    reactions = list(itertools.chain.from_iterable(reactions))
    # remove reverse reaction
    reactions = findDegeneracies(reactions)
    # add reverse attribute to families with ownReverse
    for rxn in reactions:
        family = getDB('kinetics').families[rxn.family]
        if family.ownReverse:
            family.addReverseAttribute(rxn)
    # fix the degneracy of (not ownReverse) reactions found in the backwards
    # direction
    correctDegeneracyOfReverseReactions(reactions)
    reduceSameReactantDegeneracy(reactions)
    # get a molecule list with species indexes
    zippedList = []
    for spec in speciesTuple:
        for mol in spec.molecule:
            zippedList.append((mol,spec.index))

    molecules, reactantIndices = zip(*zippedList)

    deflate(reactions,
            [spec for spec in speciesTuple],
            [spec.index for spec in speciesTuple])

    return reactions

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

def getMoleculeTuples(speciesTuple):
    """
    returns a list of molule tuples from given speciesTuples.

    The species objects should already have resonance isomers
    generated for the function to work
    """
    combos = []
    if len(speciesTuple) == 1:#unimolecular reaction
        spc, = speciesTuple
        mols = [(mol, spc.index) for mol in spc.molecule]
        combos.extend([(combo,) for combo in mols])
    elif len(speciesTuple) == 2:#bimolecular reaction
        spcA, spcB = speciesTuple
        molsA = [(mol, spcA.index) for mol in spcA.molecule]
        molsB = [(mol, spcB.index) for mol in spcB.molecule]
        combos.extend(itertools.product(molsA, molsB))
    return combos

def reactMolecules(moleculeTuples):
    """
    Performs a reaction between
    the resonance isomers.

    The parameter contains a list of tuples with each tuple:
    (Molecule, index of the core species it belongs to)
    """

    families = getDB('kinetics').families
    
    molecules, reactantIndices = zip(*moleculeTuples)

    reactionList = []
    for _, family in families.iteritems():
        rxns = family.generateReactions(molecules)
        reactionList.extend(rxns)

    for reactant in molecules:
        reactant.clearLabeledAtoms()

    return reactionList

def findDegeneracies(rxnList, useSpeciesReaction = True):
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
                    isomorphic = rxn0.isIsomorphic(rxn,checkIdentical=False)
                    if not isomorphic:
                        identical = False
                    else:
                        identical = rxn0.isIsomorphic(rxn,checkIdentical=True)
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

def reduceSameReactantDegeneracy(rxnList):
    """
    This method reduces the degeneracy of reactions with identical reactants,
    since translational component of the transition states are already taken
    into account (so swapping the same reactant is not valid)
    
    This comes from work by Bishop and Laidler in 1965
    """
    for reaction in rxnList:
        if len(reaction.reactants) == 2 and reaction.reactants[0].isIsomorphic(reaction.reactants[1]):
            reaction.degeneracy *= 0.5
            logging.debug('Degeneracy of reaction {} was decreased by 50% to {} since the reactants are identical'.format(reaction,reaction.degeneracy))
        if reaction.reverse and len(reaction.reverse.reactants) == 2 and \
                                   reaction.reverse.reactants[0].isIsomorphic(reaction.reverse.reactants[1]):
            reaction.reverse.degeneracy *= 0.5
            logging.debug('Degeneracy of reaction {} was decreased by 50% to {} since the reactants are identical'.format(reaction.reverse,reaction.reverse.degeneracy))

def correctDegeneracyOfReverseReactions(reactionList):
    """
    This method corrects the degeneracy of reactions found when the backwards
    template is used. Given the following parameters:

        reactionList - list of reactions with their degeneracies already counted
        reactants - list/tuple of species used in the generateReactions method

    This method modifies reactionList in place and returns nothing
    
    This does not adjust for identical reactants, you need to use `reduceSameReactantDegeneracy`
    to adjust for that.
    """
    for rxn in reactionList:
        if not rxn.isForward:
            # was reverse reaction so should find degeneracy
            family = getDB('kinetics').families[rxn.family]
            if not family.ownReverse: 
                rxn.degeneracy = family.calculateDegeneracy(rxn, ignoreSameReactants=True)

def deflate(rxns, species, reactantIndices):
    """
    The purpose of this function is to replace the reactants and
    products of a reaction, stored as Molecule objects by 
    integer indices, corresponding to the species core index.

    Creates a dictionary with Molecule objects as keys and newly 
    created Species objects as values.

    It iterates over the reactantIndices array, with elements in this array
    corresponding to the indices of the core species. It creates a 
    Molecule -> index entry in the previously created dictionary.

    It iterates over the reaction list, and iteratively updates the
    created dictionary as more reactions are processed.    
    """    

    molDict = {}

    for i, coreIndex in enumerate(reactantIndices):
        if coreIndex != -1:
            for mol in species[i].molecule:
                molDict[mol] = coreIndex

    for rxn in rxns:
        deflateReaction(rxn, molDict)
        try:
            deflateReaction(rxn.reverse, molDict) 
        except AttributeError:
            pass



def reactAll(coreSpcList, numOldCoreSpecies, unimolecularReact, bimolecularReact):
    """
    Reacts the core species list via uni- and bimolecular reactions.
    """

    # Select reactive species that can undergo unimolecular reactions:
    spcTuples = [(coreSpcList[i],)
     for i in xrange(numOldCoreSpecies) if (unimolecularReact[i] and coreSpcList[i].reactive)]

    for i in xrange(numOldCoreSpecies):
        for j in xrange(i, numOldCoreSpecies):
            # Find reactions involving the species that are bimolecular
            # This includes a species reacting with itself (if its own concentration is high enough)
            if bimolecularReact[i,j]:
                if coreSpcList[i].reactive and coreSpcList[j].reactive:
                    spcTuples.append((coreSpcList[i], coreSpcList[j]))

    rxns = list(react(*spcTuples))
    return rxns

def deflateReaction(rxn, molDict):
    """
    This function deflates a single reaction holding speices objects, and uses the provided 
    dictionary to populate reactants/products/pairs with integer indices,
    if possible.

    If the Molecule object could not be found in the dictionary, a new
    dictionary entry is created, using the Species object as the value
    for the entry.

    The reactants/products/pairs of both the forward and reverse reaction 
    object are populated with the value of the dictionary, either an
    integer index, or either a Species object.
    """
    for spec in itertools.chain(rxn.reactants, rxn.products):
        if not spec.molecule[0] in molDict:
            molDict[spec.molecule[0]] = spec

    rxn.reactants = [molDict[spec.molecule[0]] for spec in rxn.reactants]
    rxn.products = [molDict[spec.molecule[0]] for spec in rxn.products]
    rxn.pairs = [(molDict[reactant.molecule[0]], molDict[product.molecule[0]]) for reactant, product in rxn.pairs]
