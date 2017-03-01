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
Contains functions for generating reactions.
"""
import logging
import itertools

from rmgpy.molecule.molecule import Molecule
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

    reactions = map(reactMolecules,combos)
    reactions = list(itertools.chain.from_iterable(reactions))
    findDegeneracies(reactions)
    reduceSameReactantDegeneracy(reactions)

    deflate(reactions,
            [spec.molecule[0] for spec in speciesTuple],
            [spec.index for spec in speciesTuple])

    return reactions

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

def findDegeneracies(rxnList):
    """
    given a list of reaction object with Molecule objects, this method 
    removes degenerate reactions and
    increments the degeneracy of the reaction object. This method modifies
    rxnList in place and does not return anything.

    This algorithm used to exist in family.__generateReactions, but was moved
    here because it didn't have any family dependence.
    """

    index0 = 0
    while index0 < len(rxnList):
        reaction0 = rxnList[index0]

        # obtain species with all resonance isomers
        storeReactants = reaction0.reactants[:]
        storeProducts = reaction0.products[:]
        for i, mol in enumerate(reaction0.reactants):
            spec = Species(molecule = [mol])
            spec.generateResonanceIsomers()
            reaction0.reactants[i] = spec
        for i, mol in enumerate(reaction0.products):
            spec = Species(molecule = [mol])
            spec.generateResonanceIsomers()
            reaction0.products[i] = spec

        # Remove duplicates from the reaction list
        index = index0 + 1
        while index < len(rxnList):
            reaction = rxnList[index]
            # ensure same family and structure
            # (template can be added later after multiple TS's are enabled with all([template0 in reaction.template for template0 in reaction0.template]))
            if reaction.family == reaction0.family and \
                 reaction0.isIsomorphic(reaction): 
                reaction0.degeneracy += reaction.degeneracy
                try:
                    reaction0.reverse.degeneracy += reaction.reverse.degeneracy
                except AttributeError:
                    pass
                rxnList.remove(reaction)
                
            else:
                index += 1

        index0 += 1

        # reconvert the reaction back to molecule objects (for deflating)
        reaction0.reactants = storeReactants
        reaction0.products = storeProducts


def reduceSameReactantDegeneracy(rxnList):
    """
    This method reduces the degeneracy of reactions with identical reactants,
    since translational component of the transition states are already taken
    into account (so swapping the same reactant is not valid)
    
    This comes from work by Bishop and Laidler in 1965
    """
    for reaction in rxnList:
        if len(reaction.reactants) == 2 and reaction.reactants[0].isIsomorphic(reaction.reactants[1]):
            assert reaction.degeneracy % 2 == 0, 'Reaction with isomorphic reactants had an odd degeneracy of {0}. The reaction is {1}'.format(reaction.degeneracy, str(reaction))
            reaction.degeneracy /= 2

    

def deflate(rxns, molecules, reactantIndices):
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
            molDict[molecules[i]] = coreIndex 

    for rxn in rxns:
        deflateReaction(rxn, molDict)
        try:
            deflateReaction(rxn.reverse, molDict) 
        except AttributeError, e:
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
    This function deflates a single reaction, and uses the provided 
    dictionary to populate reactants/products/pairs with integer indices,
    if possible.

    If the Molecule object could not be found in the dictionary, a new
    dictionary entry is created, creating a new Species object as the value
    for the entry.

    The reactants/products/pairs of both the forward and reverse reaction 
    object are populated with the value of the dictionary, either an
    integer index, or either a Species object.
    """

    for mol in itertools.chain(rxn.reactants, rxn.products):
        if not mol in molDict:
            molDict[mol] = Species(molecule=[mol])

    rxn.reactants = [molDict[mol] for mol in rxn.reactants]
    rxn.products = [molDict[mol] for mol in rxn.products]
    rxn.pairs = [(molDict[reactant], molDict[product]) for reactant, product in rxn.pairs]