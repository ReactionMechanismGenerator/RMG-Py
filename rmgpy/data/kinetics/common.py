#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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

from rmgpy.data.base import DatabaseError, LogicNode
from rmgpy.reaction import Reaction, ReactionError
from rmgpy.molecule import Molecule, Group
from rmgpy.species import Species

################################################################################

# The names of all of the RMG reaction families that are bimolecular
BIMOLECULAR_KINETICS_FAMILIES = [
    'H_Abstraction',
    'R_Addition_MultipleBond',
    'R_Recombination',
    'Disproportionation',
    '1+2_Cycloaddition',
    '2+2_cycloaddition_Cd',
    '2+2_cycloaddition_CO',
    '2+2_cycloaddition_CCO',
    'Diels_alder_addition',
    '1,2_Insertion',
    '1,3_Insertion_CO2',
    '1,3_Insertion_ROR',
    'R_Addition_COm',
    'Oa_R_Recombination',
    'Substitution_O',
    'SubstitutionS',
    'R_Addition_CSm',
    '1,3_Insertion_RSR',
]

# The names of all of the RMG reaction families that are unimolecular
UNIMOLECULAR_KINETICS_FAMILIES = [
    'intra_H_migration',
    'Birad_recombination',
    'intra_OH_migration',
    'HO2_Elimination_from_PeroxyRadical',
    'Cyclic_Ether_Formation',
    'Enol_keto_tautomerism',
    'Intra_R_Add_Exocyclic',
    'Intra_R_Add_Endocyclic',
    '1,2-Birad_to_alkene',
    'Intra_Disproportionation',
    'Korcek_step1',
    'Korcek_step2',
    '1,2_shiftS',
    'intra_substitutionCS_cyclization',
    'intra_substitutionCS_isomerization',
    'intra_substitutionS_cyclization',
    'intra_substitutionS_isomerization',
]

################################################################################

class KineticsError(Exception):
    """
    An exception for problems with kinetics. Pass a string describing the problem.
    """
    pass

class UndeterminableKineticsError(ReactionError):
    """
    An exception raised when attempts to estimate appropriate kinetic parameters
    for a chemical reaction are unsuccessful.
    """
    def __init__(self, reaction, message=''):
        new_message = 'Kinetics could not be determined. '+message
        ReactionError.__init__(self,reaction,new_message)

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

    if isinstance(entry.item, Reaction):
        for i, reactant in enumerate(entry.item.reactants):
            if isinstance(reactant, Molecule):
                f.write('    reactant{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(reactant.toAdjacencyList(removeH=True))
                f.write('""",\n')
            elif isinstance(reactant, Species):
                f.write('    reactant{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(reactant.molecule[0].toAdjacencyList(label=reactant.label, removeH=True))
                f.write('""",\n')
            elif isinstance(reactant, Group):
                f.write('    group{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(reactant.toAdjacencyList())
                f.write('""",\n')
            elif isinstance(reactant, LogicNode):
                f.write('    group{0:d} = "{1}",\n'.format(i+1, reactant))
        for i, product in enumerate(entry.item.products):
            if isinstance(product, Molecule):
                f.write('    product{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(product.toAdjacencyList(removeH=True))
                f.write('""",\n')
            elif isinstance(reactant, Species):
                f.write('    product{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(product.molecule[0].toAdjacencyList(label=product.label, removeH=True))
                f.write('""",\n')
        if not isinstance(entry.item.reactants[0], Group) and not isinstance(entry.item.reactants[0], LogicNode):
            f.write('    degeneracy = {0:d},\n'.format(entry.item.degeneracy))
        if entry.item.duplicate: 
            f.write('    duplicate = {0!r},\n'.format(entry.item.duplicate))
        if not entry.item.reversible:
            f.write('    reversible = {0!r},\n'.format(entry.item.reversible))
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
    else:
        f.write('    reference = None,\n')
    
    f.write('    referenceType = "{0}",\n'.format(entry.referenceType))
    if entry.rank is not None:
        f.write('    rank = {0},\n'.format(entry.rank))
    f.write('    shortDesc = u"""')
    f.write(entry.shortDesc)
    f.write('""",\n')
    f.write('    longDesc = \n')
    f.write('u"""\n')
    f.write(entry.longDesc.strip() + "\n")
    f.write('""",\n')

    f.write('    history = [\n')
    for time, user, action, description in entry.history:
        f.write('        ("{0}","{1}","{2}","""{3}"""),\n'.format(time, user, action, description))
    f.write('    ],\n')

    f.write(')\n\n')
