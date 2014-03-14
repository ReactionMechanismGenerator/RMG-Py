#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a 
#   copy of this software and associated documentation files (the "Software"), 
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#   and/or sell copies of the Software, and to permit persons to whom the 
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
#   DEALINGS IN THE SOFTWARE. 
#
################################################################################

"""
This module contains functionality for reading from and writing to the
adjacency list format used by Reaction Mechanism Generator (RMG).
"""

import re
from .molecule import Atom, Bond
from .group import GroupAtom, GroupBond
#import chempy.molecule.atomtype as atomtypes

################################################################################

class InvalidAdjacencyListError(Exception):
    """
    An exception used to indicate that an RMG-style adjacency list is invalid.
    Pass a string describing the reason the adjacency list is invalid
    """
    pass

################################################################################

def fromAdjacencyList(adjlist, group=False, saturateH=False):
    """
    Convert a string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects.
    """
    atoms = []
    atomdict = {}
    bonds = {}
    multiplicity = -1
    
    try:
        
        adjlist = adjlist.strip()
        lines = adjlist.splitlines()
        if adjlist == '' or len(lines) == 0:
            raise InvalidAdjacencyListError('Empty adjacency list.')

        # Skip the first line if it contains a label
        if len(lines[0].split()) == 1:
            label = lines.pop(0)
            if len(lines) == 0:
                raise InvalidAdjacencyListError('No atoms specified in adjacency list.')
            
        # Skip the second line if it contains a multiplicity
        if len(lines[0].split()) == 1:
            multiplicity = int(lines.pop(0))
            if len(lines) == 0:
                raise InvalidAdjacencyListError('No atoms specified in adjacency list.')
        
        mistake1 = re.compile('\{[^}]*\s+[^}]*\}')
        # Iterate over the remaining lines, generating Atom or GroupAtom objects
        for line in lines:

            # Sometimes people put spaces after commas, which messes up the
            # parse-by-whitespace. Examples include '{Cd, Ct}'.
            if mistake1.search(line):
                raise InvalidAdjacencyListError(
                    "Shouldn't have spaces inside braces: {0}".format(mistake1.search(line).group())
                    )

            # Sometimes commas are used to delimit bonds in the bond list,
            # so replace them just in case
            line = line.replace('},{', '} {')
            
            data = line.split()

            # Skip if blank line
            if len(data) == 0: continue

            # First item is index for atom
            # Sometimes these have a trailing period (as if in a numbered list),
            # so remove it just in case
            aid = int(data[0].strip('.'))

            # If second item starts with '*', then atom is labeled
            label = ''; index = 1
            if data[1][0] == '*':
                label = data[1]
                index += 1

            # Next is the element or functional group element
            # A list can be specified with the {,} syntax
            atomType = data[index]
            if atomType[0] == '{':
                atomType = atomType[1:-1].split(',')
            else:
                atomType = [atomType]
            index += 1
            
            # Next is the electron state
            radicalElectrons = []
            elecState = data[index].upper()
            if elecState[0] == '{':
                elecState = elecState[1:-1].split(',')
            else:
                elecState = [elecState]
            for e in elecState:
                if e == '0':
                    radicalElectrons.append(0)
                elif e == '1':
                    radicalElectrons.append(1)
                elif e == '2':
                    radicalElectrons.append(2)
                elif e == '3':
                    radicalElectrons.append(3)
                elif e == '4':
                    radicalElectrons.append(4)
                elif e == 'X':
                    radicalElectrons.extend([0,1,2,2])
                else:
                    raise InvalidAdjacencyListError('Invalid number of radicals.')
            index += 1
            
            # Next number defines the number of lone electron pairs (if provided)
            lonePairElectrons = -1
         #   if not group and len(data) > index:
            if len(data) > index:
                lpState = data[index]
                if lpState[0] != '{':
                    if lpState == '0':
                        lonePairElectrons = 0
                    elif lpState == '1':
                        lonePairElectrons = 1
                    elif lpState == '2':
                        lonePairElectrons = 2
                    elif lpState == '3':
                        lonePairElectrons = 3
                    elif lpState == '4':
                        lonePairElectrons = 4
                    else:
                        raise InvalidAdjacencyListError('Invalid number of lone electron pairs.')
                    index += 1
                else:
                    lonePairElectrons = -1
            else:
                lonePairElectrons = -1
            
            # Create a new atom based on the above information
            if group:
                atom = GroupAtom(atomType, radicalElectrons, [0 for e in radicalElectrons], label, [lonePairElectrons])
            else:
                atom = Atom(atomType[0], radicalElectrons[0], 0, label, lonePairElectrons)

            # Add the atom to the list
            atoms.append(atom)
            atomdict[aid] = atom
            
            # Process list of bonds
            bonds[aid] = {}
            for datum in data[index:]:

                # Sometimes commas are used to delimit bonds in the bond list,
                # so strip them just in case
                datum = datum.strip(',')
                
                aid2, comma, order = datum[1:-1].partition(',')
                aid2 = int(aid2)
                if aid == aid2:
                    raise InvalidAdjacencyListError('Attempted to create a bond between atom {0:d} and itself.'.format(aid))
                
                if order[0] == '{':
                    order = order[1:-1].split(',')
                else:
                    order = [order]

                bonds[aid][aid2] = order

        # Check consistency using bonddict
        for atom1 in bonds:
            for atom2 in bonds[atom1]:
                if atom2 not in bonds:
                    raise InvalidAdjacencyListError('Atom {0:d} not in bond dictionary.'.format(atom2))
                elif atom1 not in bonds[atom2]:
                    raise InvalidAdjacencyListError('Found bond between {0:d} and {1:d}, but not the reverse.'.format(atom1, atom2))
                elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                    raise InvalidAdjacencyListError('Found bonds between {0:d} and {1:d}, but of different orders "{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1]))

        # Convert bonddict to use Atom[group] and Bond[group] objects
        atomkeys = atomdict.keys()
        atomkeys.sort()
        for aid1 in atomkeys:
            atomkeys2 = bonds[aid1].keys()
            atomkeys2.sort()
            for aid2 in atomkeys2:
                if aid1 < aid2:
                    atom1 = atomdict[aid1]
                    atom2 = atomdict[aid2]
                    order = bonds[aid1][aid2]
                    if group:
                        bond = GroupBond(atom1, atom2, order)
                    elif len(order) == 1:
                        bond = Bond(atom1, atom2, order[0])
                    else:
                        raise InvalidAdjacencyListError('Multiple bond orders specified for an atom in a Molecule.')
                    atom1.edges[atom2] = bond
                    atom2.edges[atom1] = bond
        
        if saturateH:
            # Add explicit hydrogen atoms to complete structure if desired
            if not group:
                valences = {'H': 1, 'C': 4, 'O': 2, 'N': 3, 'S': 2, 'Si': 4, 'Cl': 1, 'He': 0, 'Ne': 0, 'Ar': 0}
                orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}
                newAtoms = []
                for atom in atoms:
                    try:
                        valence = valences[atom.symbol]
                    except KeyError:
                        raise InvalidAdjacencyListError('Cannot add hydrogens to adjacency list: Unknown valence for atom "{0}".'.format(atom.symbol))
                    radical = atom.radicalElectrons
                    order = 0
                    for atom2, bond in atom.bonds.items():
                        order += orders[bond.order]
                    count = valence - radical - int(order)
                    for i in range(count):
                        a = Atom('H', 0, 1, 0, '')
                        b = Bond(atom, a, 'S')
                        newAtoms.append(a)
                        atom.bonds[a] = b
                        a.bonds[atom] = b
                atoms.extend(newAtoms)
        
        # Calculate the number of lone pair electrons requiring molecule with all hydrogen atoms present
        if not group and lonePairElectrons == -1:
            valences = {'H': 1, 'C': 4, 'O': 2, 'N': 3, 'S': 2, 'Si': 4, 'He': 0, 'Ne': 0, 'Ar': 0, 'Cl': 1}
            orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}
            for atom in atoms:
                if not atom.isHydrogen():
                    try:
                        valence = valences[atom.symbol]
                    except KeyError:
                        raise InvalidAdjacencyListError('Cannot add hydrogens to adjacency list: Unknown valence for atom "{0}".'.format(atom.symbol))
                    radical = atom.radicalElectrons
                    order = 0
                    for atom2, bond in atom.bonds.items():
                        order += orders[bond.order]
                    lonePairs = 4 - order - radical
                    charge = 8 - valence - order - radical - 2*lonePairs
                    atom.setLonePairs(lonePairs)
                    atom.updateCharge()
                else:
                    valence = valences[atom.symbol]
                    radical = atom.radicalElectrons
                    order = 0
                    for atom2, bond in atom.bonds.items():
                        order += orders[bond.order]
                    lonePairs = 1 - order - radical
                    charge = 2 - valence - order - radical - 2*lonePairs
                    atom.setLonePairs(lonePairs)
                    atom.updateCharge()
        elif not group:
            for atom in atoms:
                atom.updateCharge()
                    
    except InvalidAdjacencyListError:
        print adjlist
        raise
    
    if not group:
        if multiplicity == -1:
            nRad = 0
            for atom in atoms:
                nRad += atom.radicalElectrons
            multiplicity = nRad + 1
            
        return atoms, multiplicity
    
    else:
        
        return atoms

################################################################################

def getElectronState(radicalElectrons, spinMultiplicity):
    """
    Return the electron state corresponding to the given number of radical
    electrons `radicalElectrons` and spin multiplicity `spinMultiplicity`.
    Raise a :class:`ValueError` if the electron state cannot be determined.
    """
    electronState = ''
    if radicalElectrons == 0: 
        electronState = '0'
    elif radicalElectrons == 1:
        electronState = '1'
    elif radicalElectrons == 2 and spinMultiplicity == 1:
        electronState = '2S'
    elif radicalElectrons == 2 and spinMultiplicity == 3: 
        electronState = '2T'
    elif radicalElectrons == 3 and spinMultiplicity == 2: 
        electronState = '3D'
    elif radicalElectrons == 3 and spinMultiplicity == 4: 
        electronState = '3Q'
    elif radicalElectrons == 4 and spinMultiplicity == 1: 
        electronState = '4S'
    elif radicalElectrons == 4 and spinMultiplicity == 3: 
        electronState = '4T'
    elif radicalElectrons == 4 and spinMultiplicity == 5: 
        electronState = '4V'
    else:
        raise ValueError('Unable to determine electron state for {0:d} radical electrons with spin multiplicity of {1:d}.'.format(radicalElectrons, spinMultiplicity))
    return electronState

def toAdjacencyList(atoms, multiplicity=0, label=None, group=False, removeH=False, removeLonePairs=False, printMultiplicity=False):
    """
    Convert a chemical graph defined by a list of `atoms` into a string
    adjacency list.
    """
    adjlist = ''

    # Don't remove hydrogen atoms if the molecule consists only of hydrogen atoms
    try:
        if removeH and all([atom.element.symbol == 'H' for atom in atoms]): removeH = False
    except AttributeError:
        pass

    if label: adjlist += label + '\n'
    
    if printMultiplicity: adjlist += str(multiplicity) + '\n'

    # Determine the numbers to use for each atom
    atomNumbers = {}
    index = 0
    for atom in atoms:
        if removeH and atom.element.symbol == 'H' and atom.label == '': continue
        atomNumbers[atom] = '{0:d}'.format(index + 1)
        index += 1
    
    atomLabels = dict([(atom, '{0}'.format(atom.label)) for atom in atomNumbers])
    
    atomTypes = {}
    atomElectronStates = {}
    atomLonePairs = {}
    if group:
        for atom in atomNumbers:
            # Atom type(s)
            if len(atom.atomType) == 1:
                atomTypes[atom] = atom.atomType[0].label
            else:
                atomTypes[atom] = '{{{0}}}'.format(','.join([a.label for a in atom.atomType]))
            # Electron state(s)
            if len(atom.radicalElectrons) == 1: 
                atomElectronStates[atom] = str(atom.radicalElectrons[0])
            else:
                atomElectronStates[atom] = '{{{0}}}'.format(','.join([str(radical) for radical in atom.radicalElectrons]))  
    else:
        for atom in atomNumbers:
            # Atom type
            atomTypes[atom] = '{0}'.format(atom.element.symbol)
            # Electron state(s)
            atomElectronStates[atom] = '{0}'.format(str(atom.radicalElectrons))    
            if not removeLonePairs:
                # Lone Pair(s)
                atomLonePairs[atom] = atom.lonePairs
    
    # Determine field widths
    atomNumberWidth = max([len(s) for s in atomNumbers.values()]) + 1
    atomLabelWidth = max([len(s) for s in atomLabels.values()])
    if atomLabelWidth > 0: atomLabelWidth += 1
    atomTypeWidth = max([len(s) for s in atomTypes.values()]) + 1
    atomElectronStateWidth = max([len(s) for s in atomElectronStates.values()])
    atomLonePairWidth = 2
    
    # Assemble the adjacency list
    for atom in atoms:
        if atom not in atomNumbers: continue

        # Atom number
        adjlist += '{0:<{1:d}}'.format(atomNumbers[atom], atomNumberWidth)
        # Atom label
        adjlist += '{0:<{1:d}}'.format(atomLabels[atom], atomLabelWidth)
        # Atom type(s)
        adjlist += '{0:<{1:d}}'.format(atomTypes[atom], atomTypeWidth)
        # Electron state(s)
        adjlist += '{0:<{1:d}}'.format(atomElectronStates[atom], atomElectronStateWidth)
        if group == False and not removeLonePairs:
            # Lone Pair(s)
            adjlist += '{0:>{1:d}}'.format(atomLonePairs[atom], atomLonePairWidth)
        
        # Bonds list
        atoms2 = atom.bonds.keys()
        # sort them the same way as the atoms
        atoms2.sort(key=atoms.index)

        for atom2 in atoms2:
            if atom2 not in atomNumbers: continue

            bond = atom.bonds[atom2]
            adjlist += ' {{{0},'.format(atomNumbers[atom2])

            # Bond type(s)
            if group:
                if len(bond.order) == 1:
                    adjlist += bond.order[0]
                else:
                    adjlist += '{{{0}}}'.format(','.join(bond.order))
            else:
                adjlist += bond.order
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist
