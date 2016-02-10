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

    smiA = [mol.toSMILES() for mol in spcA.molecule]

    smiB = [[mol.toSMILES() for mol in spcB.molecule] for spcB in speciesList if spcB.reactive]
    smiB = list(itertools.chain.from_iterable(smiB))
    
    if not list(smiB):
        combos = [(smi,) for smi in smiA]
    else:
        combos = list(itertools.product(smiA, smiB))

    results = map_(
                WorkerWrapper(reactMolecules),
                combos
            )

    reactionList = itertools.chain.from_iterable(results)
    return reactionList

def reactMolecules(smis):
    """
    Performs a reaction between
    the resonance isomers for the given family key.
    """
    molecules = [Molecule(SMILES=smi) for smi in smis]

    rxns = []
    families = getDB('kinetics').families
    for key, family in families.iteritems():

        
        reactionList = family.generateReactions(molecules)

        for reactant in molecules:
            reactant.clearLabeledAtoms()

        rxns.extend(reactionList)

    return rxns
