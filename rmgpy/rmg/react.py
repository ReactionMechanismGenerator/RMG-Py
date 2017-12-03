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
Contains functions for generating reactions.
"""
import itertools

from rmgpy.data.rmg import getDB
from rmgpy.scoop_framework.util import map_


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
    Given a tuple of Species objects, generates all possible reactions
    from the loaded reaction families and combines degenerate reactions.

    The generated reactions are deflated.
    """
    speciesTuple = tuple([spc.copy(deep=True) for spc in speciesTuple])

    reactions = getDB('kinetics').generate_reactions_from_families(speciesTuple)

    deflate(reactions,
            [spec for spec in speciesTuple],
            [spec.index for spec in speciesTuple])

    return reactions


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
    This function deflates a single reaction holding species objects, and uses the provided
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
    try:
        rxn.pairs = [(molDict[reactant.molecule[0]], molDict[product.molecule[0]]) for reactant, product in rxn.pairs]
    except ValueError:
        rxn.pairs = None
