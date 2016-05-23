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
from rmgpy.scoop_framework.util import map_, WorkerWrapper
from rmgpy.species import Species
        
def react(spcA, speciesList=[]):
    """
    Generate reactions between spcA and the list of 
    species for all the reaction families available.

    Returns an empty list if the spcA is non-reactive.

    For the spcA, a list of tuples is created for each
    resonance isomer of the species. Each tuple consists of (Molecule, index)
    with the index the species index of the Species object.

    For each Species in the the speciesList, its corresponding
    resonance isomers are stored in a tuple ([Molecule], index) with the index 
    the species index of the Species object.

    Each tuple ([Molecule], index) is expanded into a list of tuples (Molecule, index)
    resulting in one large list [(Molecule, index)].

    Possible combinations between the spcA, and a species from the 
    speciesList is obtained by taking the combinatorial product of the
    two generated [(Molecule, index)] lists.
    """
    if not spcA.reactive: return []
    
    molsA = [(mol, spcA.index) for mol in spcA.molecule]

    molsB = [(spcB.molecule, spcB.index) for spcB in speciesList if spcB.reactive]

    temp = []
    for mols, index in molsB:
        for molB in mols:
            temp.append((molB, index))
    molsB = temp

    if not molsB:
        combos = [(t,) for t in molsA]
    else:
        combos = list(itertools.product(molsA, molsB))

    results = map_(
                WorkerWrapper(reactMolecules),
                combos
            )

    reactionList = itertools.chain.from_iterable(results)
    return reactionList

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

    for rxn in reactionList:
        deflate(rxn, molecules, reactantIndices)

    return reactionList

def deflate(rxn, reactants, reactantIndices):
    """
    Creates a dictionary with Molecule objects as keys and newly 
    creatd Species objects as values.

    Iterates over the reactantIndices array.
    The elements in this array correspond to the indices of the
    core species.

    The Species object, stored as a value of the newly created dictionary,
    is retrieved by using a Molecule object from the reactants array 
    as a dictionary key. The elements in this array correspond to the 
    Molecule objects that were used as reactants to generate reactions.

    The Species object is replaced by the core species index in the newly
    created dictionary.

    The reactants, products, and pairs objects of the reaction
    are replaced by index integers for those cases the Molecule object
    exists as a key in the newly created dictionary.
    """    

    molDict = {}
    for mol in itertools.chain(rxn.reactants, rxn.products):
        molDict[mol] = Species(molecule=[mol])

    for i, coreIndex in enumerate(reactantIndices):
        if coreIndex != -1:
            molDict[reactants[i]] = coreIndex 

    rxn.reactants = [molDict[mol] for mol in rxn.reactants]
    rxn.products = [molDict[mol] for mol in rxn.products]
    rxn.pairs = [(molDict[reactant],molDict[product]) for reactant, product in rxn.pairs]