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
        
def react(spcA, speciesList=[]):
    """
    Generate reactions between spcA and the list of 
    species for all the reaction families available.
    """
    if not spcA.reactive: return []
    
    molsB = [spcB.molecule for spcB in speciesList if spcB.reactive]
    molsB = list(itertools.chain.from_iterable(molsB))
    
    if not molsB:
        combos = [(mol,) for mol in spcA.molecule]
    else:
        combos = list(itertools.product(spcA.molecule, molsB))

    results = map_(
                WorkerWrapper(reactMolecules),
                combos
            )

    reactionList = itertools.chain.from_iterable(results)
    return reactionList

def reactMolecules(molecules):
    """
    Performs a reaction between
    the resonance isomers.
    """
    
    families = getDB('kinetics').families
    
    reactionList = []
    for _, family in families.iteritems():
        rxns = family.generateReactions(molecules)
        reactionList.extend(rxns)

    for reactant in molecules:
        reactant.clearLabeledAtoms()

    return reactionList
