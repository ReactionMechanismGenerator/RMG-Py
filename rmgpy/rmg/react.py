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

from rmgpy.data.rmg import getDB
from rmgpy.scoop_framework.util import map_, WorkerWrapper
        

def reactFamilies(spcA, speciesList=[]):
    """
    Generate reactions between spcA and the list of 
    species for all the reaction families available.
    """
    if not spcA.reactive: return []
    
    families = getDB('kinetics').families
    familyKeys = families.keys()
    familieCount = len(familyKeys)
    results = map_(
                WorkerWrapper(reactFamily),
                familyKeys,
                [spcA.copy(deep=True)] * familieCount,
                [speciesList] * familieCount
            )
    reactionList = itertools.chain.from_iterable(results)
    return reactionList

def reactFamily(familyKey, spcA, speciesList):
    """
    Generate uni and bimolecular reactions for one specific family.
    :return: a list of new reactions
    """

    if not speciesList:
        combos = [[spcA]]
    else:
        reactive_species = [spc for spc in speciesList if spc.reactive]
        if reactive_species:
            combos = list(itertools.product(reactive_species, [spcA]))
        else:
            return []

    results = map_(
                WorkerWrapper(reactSpecies),
                combos,
                [familyKey] * len(combos),
                )

    # flatten list of lists:
    reactionList = list(itertools.chain.from_iterable(results))
    
    return reactionList

def reactSpecies(speciesList, familyKey):
    """
    Performs a reaction between the
    species in the list for the given family key.
    """

    molList = [spc.molecule for spc in speciesList]

    combos = list(itertools.product(*molList))

    results = map_(
                WorkerWrapper(reactMolecules),
                combos,
                [familyKey] * len(combos),
                )

    # flatten list of lists:
    reactionList = list(itertools.chain.from_iterable(results))

    return reactionList

def reactMolecules(molecules, familyKey):
    """
    Performs a reaction between
    the resonance isomers for the given family key.
    """

    families = getDB('kinetics').families
    family = families[familyKey]

    reactionList = family.generateReactions(list(molecules))

    for reactant in molecules:
        reactant.clearLabeledAtoms()

    return reactionList
