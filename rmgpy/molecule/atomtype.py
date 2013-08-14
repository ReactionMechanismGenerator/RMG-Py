#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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
This module defines the atom types that are available for representing
molecular functional groups and substructure patterns. Each available atom type
is defined as an instance of the :class:`AtomType` class. The atom types 
themselves are available in the ``atomTypes`` module-level variable, or as
the return value from the :meth:`getAtomType()` method.

If you want to change which atom types are available in RMG and/or what they
represent, this should be the only module you need to change to do so.
"""

import cython

################################################################################

class AtomTypeError(Exception):
    """
    An exception to be raised when an error occurs while working with atom
    types. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

################################################################################

class AtomType:
    """
    A class for internal representation of atom types. Using unique objects
    rather than strings allows us to use fast pointer comparisons instead of
    slow string comparisons, as well as store extra metadata. In particular,
    we store metadata describing the atom type's hierarchy with regard to other
    atom types, and the atom types that can result when various actions
    involving this atom type are taken. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `label`             ``str``             A unique label for the atom type
    `generic`           ``list``            The atom types that are more generic than this one
    `specific`          ``list``            The atom types that are more specific than this one
    `incrementBond`     ``list``            The atom type(s) that result when an adjacent bond's order is incremented
    `decrementBond`     ``list``            The atom type(s) that result when an adjacent bond's order is decremented
    `formBond`          ``list``            The atom type(s) that result when a new single bond is formed to this atom type
    `breakBond`         ``list``            The atom type(s) that result when an existing single bond to this atom type is broken
    `incrementRadical`  ``list``            The atom type(s) that result when the number of radical electrons is incremented
    `decrementRadical`  ``list``            The atom type(s) that result when the number of radical electrons is decremented
    =================== =================== ====================================

    """

    def __init__(self, label='', generic=None, specific=None):
        self.label = label
        self.generic = generic or []
        self.specific = specific or []
        self.incrementBond = []
        self.decrementBond = []
        self.formBond = []
        self.breakBond = []
        self.incrementRadical = []
        self.decrementRadical = []

    def __repr__(self):
        return '<AtomType "%s">' % self.label

    def __reduce__(self):
        """
        A helper function used when pickling an AtomType object.
        """
        d = {
            'label': self.label,
            'generic': self.generic,
            'specific': self.specific,
            'incrementBond': self.incrementBond,
            'decrementBond': self.decrementBond,
            'formBond': self.formBond,
            'breakBond': self.breakBond,
            'incrementRadical': self.incrementRadical,
            'decrementRadical': self.decrementRadical,
        }
        return (AtomType, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an AtomType object.
        """
        self.label = d['label']
        self.generic = d['generic']
        self.specific = d['specific']
        self.incrementBond = d['incrementBond']
        self.decrementBond = d['decrementBond']
        self.formBond = d['formBond']
        self.breakBond = d['breakBond']
        self.incrementRadical = d['incrementRadical']
        self.decrementRadical = d['decrementRadical']

    def setActions(self, incrementBond, decrementBond, formBond, breakBond, incrementRadical, decrementRadical):
        self.incrementBond = incrementBond
        self.decrementBond = decrementBond
        self.formBond = formBond
        self.breakBond = breakBond
        self.incrementRadical = incrementRadical
        self.decrementRadical = decrementRadical

    def equivalent(self, other):
        """
        Returns ``True`` if two atom types `atomType1` and `atomType2` are
        equivalent or ``False``  otherwise. This function respects wildcards,
        e.g. ``R!H`` is equivalent to ``C``.
        """
        return self is other or self in other.specific or other in self.specific

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if atom type `atomType1` is a specific case of
        atom type `atomType2` or ``False``  otherwise.
        """
        return self is other or self in other.specific

################################################################################

atomTypes = {}
atomTypes['R']    = AtomType(label='R', generic=[], specific=[
    'R!H',
    'C','Cs','Cd','Cdd','Ct','CO','Cb','Cbf', 'CS',
    'H',
    'O','Os','Od','Oa',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Ss','Sd','Sa']
)
atomTypes['R!H']  = AtomType(label='R!H', generic=['R'], specific=[
    'C','Cs','Cd','Cdd','Ct','CO','Cb','Cbf', 'CS',
    'O','Os','Od','Oa',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Ss','Sd','Sa']
)
atomTypes['C']    = AtomType('C', generic=['R','R!H'], specific=['Cs','Cd','Cdd','Ct','CO','Cb','Cbf','CS'])
atomTypes['Cs']   = AtomType('Cs', generic=['R','R!H', 'C'], specific=[])
atomTypes['Cd']   = AtomType('Cd', generic=['R','R!H', 'C'], specific=[])
atomTypes['Cdd']  = AtomType('Cdd', generic=['R','R!H', 'C'], specific=[])
atomTypes['Ct']   = AtomType('Ct', generic=['R','R!H', 'C'], specific=[])
atomTypes['CO']   = AtomType('CO', generic=['R','R!H', 'C'], specific=[])
atomTypes['Cb']   = AtomType('Cb', generic=['R','R!H', 'C'], specific=[])
atomTypes['Cbf']  = AtomType('Cbf', generic=['R','R!H', 'C'], specific=[])
atomTypes['CS']   = AtomType('CS', generic=['R','R!H', 'C'], specific=[])
atomTypes['H']    = AtomType('H', generic=['R'], specific=[])
atomTypes['O']    = AtomType('O', generic=['R','R!H'], specific=['Os','Od','Oa'])
atomTypes['Os']   = AtomType('Os', generic=['R','R!H','O'], specific=[])
atomTypes['Od']   = AtomType('Od', generic=['R','R!H','O'], specific=[])
atomTypes['Oa']   = AtomType('Oa', generic=['R','R!H','O'], specific=[])
atomTypes['Si']   = AtomType('Si', generic=['R','R!H'], specific=['Sis','Sid','Sidd','Sit','SiO','Sib','Sibf'])
atomTypes['Sis']  = AtomType('Sis', generic=['R','R!H','Si'], specific=[])
atomTypes['Sid']  = AtomType('Sid', generic=['R','R!H','Si'], specific=[])
atomTypes['Sidd'] = AtomType('Sidd', generic=['R','R!H','Si'], specific=[])
atomTypes['Sit']  = AtomType('Sit', generic=['R','R!H','Si'], specific=[])
atomTypes['SiO']  = AtomType('SiO', generic=['R','R!H','Si'], specific=[])
atomTypes['Sib']  = AtomType('Sib', generic=['R','R!H','Si'], specific=[])
atomTypes['Sibf'] = AtomType('Sibf', generic=['R','R!H','Si'], specific=[])
atomTypes['S']    = AtomType('S', generic=['R','R!H'], specific=['Ss','Sd','Sa'])
atomTypes['Ss']   = AtomType('Ss', generic=['R','R!H','S'], specific=[])
atomTypes['Sd']   = AtomType('Sd', generic=['R','R!H','S'], specific=[])
atomTypes['Sa']   = AtomType('Sa', generic=['R','R!H','S'], specific=[])

atomTypes['R'   ].setActions(incrementBond=['R'],            decrementBond=['R'],         formBond=['R'],     breakBond=['R'],     incrementRadical=['R'],     decrementRadical=['R'])
atomTypes['R!H' ].setActions(incrementBond=['R!H'],          decrementBond=['R!H'],       formBond=['R!H'],   breakBond=['R!H'],   incrementRadical=['R!H'],   decrementRadical=['R!H'])

atomTypes['C'   ].setActions(incrementBond=['C'],            decrementBond=['C'],         formBond=['C'],     breakBond=['C'],     incrementRadical=['C'],     decrementRadical=['C'])
atomTypes['Cs'  ].setActions(incrementBond=['Cd','CO','Cs'], decrementBond=[],            formBond=['Cs'],    breakBond=['Cs'],    incrementRadical=['Cs'],    decrementRadical=['Cs'])
atomTypes['Cd'  ].setActions(incrementBond=['Cdd','Ct'],     decrementBond=['Cs'],        formBond=['Cd'],    breakBond=['Cd'],    incrementRadical=['Cd'],    decrementRadical=['Cd'])
atomTypes['Cdd' ].setActions(incrementBond=[],               decrementBond=['Cd','CO','CS'], formBond=[],     breakBond=[],        incrementRadical=[],        decrementRadical=[])
atomTypes['Ct'  ].setActions(incrementBond=[],               decrementBond=['Cd'],        formBond=['Ct'],    breakBond=['Ct'],    incrementRadical=['Ct'],    decrementRadical=['Ct'])
atomTypes['CO'  ].setActions(incrementBond=['Cdd'],          decrementBond=['Cs'],        formBond=['CO'],    breakBond=['CO'],    incrementRadical=['CO'],    decrementRadical=['CO'])
atomTypes['CS'  ].setActions(incrementBond=['Cdd'],          decrementBond=['Cs'],        formBond=['CS'],    breakBond=['CS'],    incrementRadical=['CS'],    decrementRadical=['CS'])
atomTypes['Cb'  ].setActions(incrementBond=[],               decrementBond=[],            formBond=['Cb'],    breakBond=['Cb'],    incrementRadical=['Cb'],    decrementRadical=['Cb'])
atomTypes['Cbf' ].setActions(incrementBond=[],               decrementBond=[],            formBond=[],        breakBond=[],        incrementRadical=[],        decrementRadical=[])

atomTypes['H'   ].setActions(incrementBond=[],               decrementBond=[],            formBond=['H'],     breakBond=['H'],     incrementRadical=['H'],     decrementRadical=['H'])

atomTypes['O'   ].setActions(incrementBond=['O'],            decrementBond=['O'],         formBond=['O'],     breakBond=['O'],     incrementRadical=['O'],     decrementRadical=['O'])
atomTypes['Os'  ].setActions(incrementBond=['Od'],           decrementBond=[],            formBond=['Os'],    breakBond=['Os'],    incrementRadical=['Os'],    decrementRadical=['Os'])
atomTypes['Od'  ].setActions(incrementBond=[],               decrementBond=['Os'],        formBond=[],        breakBond=[],        incrementRadical=[],        decrementRadical=[])
atomTypes['Oa'  ].setActions(incrementBond=[],               decrementBond=[],            formBond=[],        breakBond=[],        incrementRadical=[],        decrementRadical=[])

atomTypes['Si'   ].setActions(incrementBond=['Si'],          decrementBond=['Si'],        formBond=['Si'],    breakBond=['Si'],    incrementRadical=['Si'],    decrementRadical=['Si'])
atomTypes['Sis'  ].setActions(incrementBond=['Sid','SiO'],   decrementBond=[],            formBond=['Sis'],   breakBond=['Sis'],   incrementRadical=['Sis'],   decrementRadical=['Sis'])
atomTypes['Sid'  ].setActions(incrementBond=['Sidd','Sit'],  decrementBond=['Sis'],       formBond=['Sid'],   breakBond=['Sid'],   incrementRadical=['Sid'],   decrementRadical=['Sid'])
atomTypes['Sidd' ].setActions(incrementBond=[],              decrementBond=['Sid','SiO'], formBond=[],        breakBond=[],        incrementRadical=[],        decrementRadical=[])
atomTypes['Sit'  ].setActions(incrementBond=[],              decrementBond=['Sid'],       formBond=['Sit'],   breakBond=['Sit'],   incrementRadical=['Sit'],   decrementRadical=['Sit'])
atomTypes['SiO'  ].setActions(incrementBond=['Sidd'],        decrementBond=['Sis'],       formBond=['SiO'],   breakBond=['SiO'],   incrementRadical=['SiO'],   decrementRadical=['SiO'])
atomTypes['Sib'  ].setActions(incrementBond=[],              decrementBond=[],            formBond=['Sib'],   breakBond=['Sib'],   incrementRadical=['Sib'],   decrementRadical=['Sib'])
atomTypes['Sibf' ].setActions(incrementBond=[],              decrementBond=[],            formBond=[],        breakBond=[],        incrementRadical=[],        decrementRadical=[])

atomTypes['S'   ].setActions(incrementBond=['S'],            decrementBond=['S'],         formBond=['S'],     breakBond=['S'],     incrementRadical=['S'],     decrementRadical=['S'])
atomTypes['Ss'  ].setActions(incrementBond=['Sd'],           decrementBond=[],            formBond=['Ss'],    breakBond=['Ss'],    incrementRadical=['Ss'],    decrementRadical=['Ss'])
atomTypes['Sd'  ].setActions(incrementBond=[],               decrementBond=['Ss'],        formBond=[],        breakBond=[],        incrementRadical=[],        decrementRadical=[])
atomTypes['Sa'  ].setActions(incrementBond=[],               decrementBond=[],            formBond=[],        breakBond=[],        incrementRadical=[],        decrementRadical=[])

for atomType in atomTypes.values():
    for items in [atomType.generic, atomType.specific,
      atomType.incrementBond, atomType.decrementBond, atomType.formBond,
      atomType.breakBond, atomType.incrementRadical, atomType.decrementRadical]:
        for index in range(len(items)):
            items[index] = atomTypes[items[index]]

def getAtomType(atom, bonds):
    """
    Determine the appropriate atom type for an :class:`Atom` object `atom`
    with local bond structure `bonds`, a ``dict`` containing atom-bond pairs.
    """

    cython.declare(atomType=str)
    cython.declare(double=cython.int, double0=cython.int, triple=cython.int, benzene=cython.int)
    
    atomType = ''
    
    # Count numbers of each higher-order bond type
    double = 0; doubleO = 0; triple = 0; benzene = 0
    for atom2, bond12 in bonds.iteritems():
        if bond12.isDouble():
            if atom2.isOxygen(): doubleO +=1
            else:                double += 1
        elif bond12.isTriple(): triple += 1
        elif bond12.isBenzene(): benzene += 1

    # Use element and counts to determine proper atom type
    if atom.symbol == 'C':
        if   double == 0 and doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Cs'
        elif double == 1 and doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Cd'
        elif double + doubleO == 2        and triple == 0 and benzene == 0: atomType = 'Cdd'
        elif double == 0 and doubleO == 0 and triple == 1 and benzene == 0: atomType = 'Ct'
        elif double == 0 and doubleO == 1 and triple == 0 and benzene == 0: atomType = 'CO'
        elif double == 0 and doubleO == 0 and triple == 0 and benzene == 2: atomType = 'Cb'
        elif double == 0 and doubleO == 0 and triple == 0 and benzene == 3: atomType = 'Cbf'
    elif atom.symbol == 'H':
        atomType = 'H'
    elif atom.symbol == 'O':
        if   double + doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Os'
        elif double + doubleO == 1 and triple == 0 and benzene == 0: atomType = 'Od'
        elif len(bonds) == 0:                                        atomType = 'Oa'
    elif atom.symbol == 'Si':
        if   double == 0 and doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Sis'
        elif double == 1 and doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Sid'
        elif double + doubleO == 2        and triple == 0 and benzene == 0: atomType = 'Sidd'
        elif double == 0 and doubleO == 0 and triple == 1 and benzene == 0: atomType = 'Sit'
        elif double == 0 and doubleO == 1 and triple == 0 and benzene == 0: atomType = 'SiO'
        elif double == 0 and doubleO == 0 and triple == 0 and benzene == 2: atomType = 'Sib'
        elif double == 0 and doubleO == 0 and triple == 0 and benzene == 3: atomType = 'Sibf'
    elif atom.symbol == 'S':
        if   double + doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Ss'
        elif double + doubleO == 1 and triple == 0 and benzene == 0: atomType = 'Sd'
        elif len(bonds) == 0:                                        atomType = 'Sa'
    elif atom.symbol == 'N' or atom.symbol == 'Ar' or atom.symbol == 'He' or atom.symbol == 'Ne':
        return None

    # Raise exception if we could not identify the proper atom type
    if atomType == '':
        raise AtomTypeError('Unable to determine atom type for atom {0}, which has {1:d} double bonds to C, {2:d} double bonds to O, {3:d} triple bonds, and {4:d} benzene bonds.'.format(atom, double, doubleO, triple, benzene))

    return atomTypes[atomType]
