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
    """
    if not spcA.reactive: return []
    
    molsA = [(mol, spcA.index) for mol in spcA.molecule]

    molsB = molsB = [spcB.molecule for spcB in speciesList if spcB.reactive]
    molsB = list(itertools.chain.from_iterable(molsB))
    molsB = [(mol, spcB.index) for mol in molsB]

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

    deflate(reactionList, molecules, reactantIndices)

    return reactionList

def deflate(reactionList, reactants, reactantIndices):
    """
    Return the reactions' reactants, products and pairs
    as containing Species objects, not Molecule objects.
    
    Replace the species objects by unique integers for molecule 
    objects that are found the parameter 'reactants' list.
    """

    for rxn in reactionList:
        molDict = {}

        for mol in rxn.reactants:
            molDict[mol] = Species(molecule=[mol])

        for mol in rxn.products:
            molDict[mol] = Species(molecule=[mol])

        for i, coreIndex in enumerate(reactantIndices):
            if coreIndex != -1:
                molDict[reactants[i]] = coreIndex 

        rxn.reactants = [molDict[mol] for mol in rxn.reactants]
        rxn.products = [molDict[mol] for mol in rxn.products]
        rxn.pairs = [(molDict[reactant],molDict[product]) for reactant, product in rxn.pairs]

    return reactionList