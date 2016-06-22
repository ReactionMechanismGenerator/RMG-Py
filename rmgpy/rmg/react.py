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
    
    combos = []

    for t in spcTuples:
        t = tuple([spc.copy(deep=True) for spc in t])
        if len(t) == 1:#unimolecular reaction
            spc, = t
            mols = [(mol, spc.index) for mol in spc.molecule]
            combos.extend([(combo,) for combo in mols])
        elif len(t) == 2:#bimolecular reaction
            spcA, spcB = t
            molsA = [(mol, spcA.index) for mol in spcA.molecule]
            molsB = [(mol, spcB.index) for mol in spcB.molecule]
            combos.extend(itertools.product(molsA, molsB))

    results = map_(
                reactMolecules,
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