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
import logging
import re
from .molecule import Atom, Bond
from .group import GroupAtom, GroupBond
#import chempy.molecule.atomtype as atomtypes

bond_orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}

class PeriodicSystem(object):
    valence_electrons_first_period_elements  = {'H':1, 'He':2}
        
    valence_electrons_second_period_elements = {'C':4, 'N':5, 'O':6, 'Ne':8}
        
    valence_electrons_third_period_elements  = {'Si':4, 'S':6, 'Cl':7, 'Ar':8}
        
    valence_electrons = {}
    valence_electrons.update(valence_electrons_first_period_elements)
    valence_electrons.update(valence_electrons_second_period_elements)
    valence_electrons.update(valence_electrons_third_period_elements)
    
    lone_pairs         = {'H': 0, 'C': 0, 'N': 1, 'O': 2, 'Si':0, 'S': 2, 'Cl':3 }
    
class Saturator(object):
    @staticmethod
    def saturate(atoms):
            '''
            Returns a list of atoms that is extended 
            (and bond attributes) by saturating the valency of the non-hydrogen atoms with an 
            appropriate number of hydrogen atoms.
            
            The required number of hydrogen atoms per heavy atom is determined as follows:
            H's =     max number of valence electrons - atom.radicalElectrons
                        - 2* atom.lonePairs - order - atom.charge
            
            '''
            global bond_orders
            newAtoms = []
            for atom in atoms:
                try:
                    max_number_of_valence_electrons = PeriodicSystem.valence_electrons[atom.symbol]
                except KeyError:
                    raise InvalidAdjacencyListError('Cannot add hydrogens to adjacency list: Unknown orbital for atom "{0}".'.format(atom.symbol))
                
                order = 0
                for _, bond in atom.bonds.items():
                    order += bond_orders[bond.order]
                    
                number_of_H_to_be_added = max_number_of_valence_electrons - atom.radicalElectrons - 2* atom.lonePairs - int(order) - atom.charge
                
                if number_of_H_to_be_added < 0:
                    raise InvalidAdjacencyListError('Incorrect electron configuration on atom.')
                    
                for _ in range(number_of_H_to_be_added):
                    a = Atom(element='H', radicalElectrons=0, charge=0, label='', lonePairs=0)
                    b = Bond(atom, a, 'S')
                    newAtoms.append(a)
                    atom.bonds[a] = b
                    a.bonds[atom] = b
            atoms.extend(newAtoms)  

class ConsistencyChecker(object):
        
    @staticmethod
    def check_partial_charge(atom):
            '''
            Checks whether the partial charge attribute of the atom checks out with 
            the theoretical one:
            
            '''
            global bond_orders
            valence = PeriodicSystem.valence_electrons[atom.symbol]
            order = 0
            for _, bond in atom.bonds.items():
                order += bond_orders[bond.order]
                
            theoretical = valence - order - atom.radicalElectrons - 2*atom.lonePairs

            if atom.charge != theoretical:
                raise InvalidAdjacencyListError('Invalid valency for atom {symbol} with {radicals} unpaired electrons, {lonePairs} pairs of electrons, and {charge} charge.'
                                                .format(symbol=atom.symbol, radicals=atom.radicalElectrons, lonePairs=atom.lonePairs, charge=atom.charge))

    @staticmethod
    def check_multiplicity(nRad, multiplicity):
        '''
        Check if the parameter multiplicity is an odd number,
        since the multiplicity should comply with the formula
        
        m = 2s + 1, with s the sum of the spin [+/- 1/2) ] of the unpaired electrons
        
        For a simple radical (nRad = 1): 
        s = +1/2 , m = 2 (doublet)
        
        For a biradical, s can be either 0 [+0.5 + (-0.5) ] or 1 [+0.5 + (+0.5) ] 
        and m = 1 (singlet) or m = 3 (triplet).
        '''
        if nRad in [0,1]:
            if multiplicity != (2*nRad/2 + 1):
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of radicals {1}.'.format(multiplicity, nRad))
        elif nRad == 2:
            if not int(multiplicity) in [1,3]: raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of radicals {1}.'.format(multiplicity, nRad))
        else: logging.info("Consistency checking of multiplicity of molecules with more than 2 unpaired electrons is not implemented yet!")
    
    @staticmethod
    def check_hund_rule(atom, multiplicity):
        '''
        It is checked whether atoms with 2 unpaired electrons on the same atom
        result in a multiplicity of 3, and not 1. 
        
        Unpaired electrons in 2 different orbitals belonging to the same atom
        should have the same spin, and hence, should result in a multiplicity of 3. 
        '''
        if atom.radicalElectrons == 2 and multiplicity == 1:
            raise InvalidAdjacencyListError("Violation of hund's rule. Invalid multiplicity of {0} because there is an atom with {1} unpaired electrons"
                                            .format(multiplicity, atom.radicalElectrons))
            
################################################################################

class InvalidAdjacencyListError(Exception):
    """
    An exception used to indicate that an RMG-style adjacency list is invalid.
    Pass a string describing the reason the adjacency list is invalid
    """
    pass

################################################################################

def fromOldAdjacencyList(adjlist, group=False, saturateH=False):
    """
    Convert a pre-June-2014 string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects. 
    It can read both "old style" that existed for years, an the "intermediate style" that
    existed for a few months in 2014, with the extra column of integers for lone pairs.
    """
    atoms = []
    atomdict = {}
    bonds = {}
    
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
        
        mistake1 = re.compile('\{[^}]*\s+[^}]*\}')
        atomicMultiplicities = {} # these are no longer stored on atoms, so we make a separate dictionary
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
            radicalElectrons = []; 
            additionalLonePairs = []
            elecState = data[index].upper()
            if elecState[0] == '{':
                elecState = elecState[1:-1].split(',')
            else:
                elecState = [elecState]
            if len(elecState) == 0:
                raise InvalidAdjacencyListError("There must be some electronic state defined for an old adjlist.")
            for e in elecState:
                if e == '0':
                    radicalElectrons.append(0); additionalLonePairs.append(0)
                elif e == '1':
                    radicalElectrons.append(1); additionalLonePairs.append(0)
                elif e == '2': 
                    if not group:
                        raise InvalidAdjacencyListError("Number of radical electrons = 2 is not specific enough for a molecule adjlist.  Please use 2S or 2T.")
                    # includes 2S and 2T
                    radicalElectrons.append(0); additionalLonePairs.append(1)
                    radicalElectrons.append(2); additionalLonePairs.append(0)
                elif e == '2S':
                    radicalElectrons.append(0); additionalLonePairs.append(1)
                elif e == '2T':
                    radicalElectrons.append(2); additionalLonePairs.append(0)
                elif e == '3':
                    if not group:
                        raise InvalidAdjacencyListError("Number of radical electrons = 3 is not specific enough for a molecule adjlist.  Please use 3D or 3Q.")
                    # includes 3D and 3Q
                    radicalElectrons.append(1); additionalLonePairs.append(1)
                    radicalElectrons.append(3); additionalLonePairs.append(0)
                elif e == '3D':
                    radicalElectrons.append(1); additionalLonePairs.append(1)
                elif e == '3Q':
                    radicalElectrons.append(3); additionalLonePairs.append(0)
                elif e == '4':
                    if not group:
                        raise InvalidAdjacencyListError("Number of radical electrons = 4 is not specific enough for a molecule adjlist.  Please use 4S, 4T, or 4V.")
                    # includes 4S, 4T, and 4V
                    radicalElectrons.append(0); additionalLonePairs.append(2)
                    radicalElectrons.append(2); additionalLonePairs.append(1)
                    radicalElectrons.append(4); additionalLonePairs.append(0)
                elif e == '4S':
                    radicalElectrons.append(0); additionalLonePairs.append(2)
                elif e == '4T':
                    radicalElectrons.append(2); additionalLonePairs.append(1)
                elif e == '4V':
                    radicalElectrons.append(4); additionalLonePairs.append(0)
                elif e == 'X':
                    if not group:
                        raise InvalidAdjacencyListError("Number of radical electrons = X is not specific enough for a molecule adjlist.  Wildcards should only be used for groups.")
                    radicalElectrons = []
            index += 1
            
            # Next number defines the number of lone electron pairs (if provided)
            lonePairsOfElectrons = None
            if len(data) > index:
                lpState = data[index]
                if lpState[0] == '{':
                    # this is the start of the chemical bonds - no lone pair info was provided
                    lonePairsOfElectrons = None
                else:
                    if lpState == '0':
                        lonePairsOfElectrons = 0
                    if lpState == '1':
                        lonePairsOfElectrons = 1
                    if lpState == '2':
                        lonePairsOfElectrons = 2
                    if lpState == '3':
                        lonePairsOfElectrons = 3
                    if lpState == '4':
                        lonePairsOfElectrons = 4
                    index += 1
            else: # no bonds or lone pair info provided.
                lonePairsOfElectrons = None

            # Create a new atom based on the above information
            if group:
                if lonePairsOfElectrons is not None:
                    lonePairsOfElectrons = [additional + lonePairsOfElectrons for additional in additionalLonePairs]
                atom = GroupAtom(atomType=atomType,
                                 radicalElectrons=sorted(set(radicalElectrons)),
                                 charge=None,
                                 label=label,
                                 lonePairs=lonePairsOfElectrons,  # Assign lonePairsOfElectrons as None if it is not explicitly provided
                                 )
                
            else:
                standardLonePairs = {'H': 0, 'C': 0, 'O': 2, 'S': 2, 'Si': 0, 'Cl': 3, 'He': 1, 'Ne': 4, 'Ar': 4}
                if lonePairsOfElectrons is not None:
                    # Intermediate adjlist representation
                    lonePairsOfElectrons = lonePairsOfElectrons + additionalLonePairs[0]
                else:
                    # Add the standard number of lone pairs with the additional lone pairs
                    lonePairsOfElectrons = standardLonePairs[atomType[0]] + additionalLonePairs[0]
                    
                atom = Atom(element=atomType[0],
                        radicalElectrons=radicalElectrons[0],
                        charge=0,
                        label=label,
                        lonePairs=lonePairsOfElectrons,
                        )
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
        
        if not group:
            if saturateH:
                # Add explicit hydrogen atoms to complete structure if desired
                orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}
                standardLonePairs = {'H': 0, 'C': 0, 'O': 2, 'S': 2, 'Si': 0, 'Cl': 3, 'He': 1, 'Ne': 4, 'Ar': 4}
                valences = {'H': 1, 'C': 4, 'O': 2, 'N': 3, 'S': 2, 'Si': 4, 'Cl': 1, 'He': 0, 'Ne': 0, 'Ar': 0}
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
                    count = valence - radical - int(order) - 2*(atom.lonePairs-standardLonePairs[atom.symbol])
                    for i in range(count):
                        a = Atom(element='H', radicalElectrons=0, charge=0, label='', lonePairs=0)
                        b = Bond(atom, a, 'S')
                        newAtoms.append(a)
                        atom.bonds[a] = b
                        a.bonds[atom] = b
                atoms.extend(newAtoms)
        
            # Calculate the multiplicity for the molecule and update the charges on each atom
            nRad = 0   # total number of radical electrons
            for atom in atoms:
                atom.updateCharge()
                nRad += atom.radicalElectrons
            multiplicity = nRad + 1     # 2 s + 1, where s is the combined spin of unpaired electrons (s = 1/2 per unpaired electron)
        
        else:
            # Don't set a multiplicity for groups when converting from an old adjlist
            multiplicity = None
                    
    except InvalidAdjacencyListError:
        logging.error("Troublesome adjacency list:\n" + adjlist)
        raise
    
    return atoms, multiplicity
###############################

re_IntermediateAdjList = re.compile('^\s*(\d*)\s+' +  # atom number digit
                          '(?P<label>\*\d*\s+)?' +  # optional label eg * or *2
                          '(?P<atomtype>\{?[A-Z]\S*)\s+' +  # atomtype eg R!H or {Cb,Cd}
                          '(?P<radicals>X|\d[STDQV]?|\{?\d[^}]*\})\s+' +  #radicals eg. X or 2T or {1,2,2T}
                          '(?P<lonepairs>\d)' +  # lone pairs eg. 0
                          '(?P<bonds>(\s+\{\d+\,(?:[SDTB]|\{.+?\})\},?)*)' +  # bonds, eg {2,S} {4,{S,D}}
                          '\s*$')  # the end!

re_OldAdjList = re.compile('^\s*(\d*)\s+' +  # atom number digit
                          '(?P<label>\*\d*\s+)?' +  # optional label eg * or *2
                          '(?P<atomtype>\{?[A-Z]\S*)\s+' +  # atomtype eg R!H or {Cb,Cd}
                          '(?P<radicals>X|\d[STDQV]?|\{?\d[^}]*\})' +  #radicals eg. X or 2T or {1,2,2T}
                          '(?P<bonds>(\s+\{\d+\,(?:[SDTB]|\{.+?\})\},?)*)' +  # bonds, eg {2,S} {4,{S,D}}
                          '\s*$')  # the end!

def fromAdjacencyList(adjlist, group=False, saturateH=False):
    """
    Convert a string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects.
    """
    atoms = []
    atomdict = {}
    bonds = {}
    multiplicity = None
    
    adjlist = adjlist.strip()
    lines = adjlist.splitlines()
    if adjlist == '' or len(lines) == 0:
        raise InvalidAdjacencyListError('Empty adjacency list.')

    # Detect old-style adjacency lists by looking at the last line's syntax
    lastLine = lines[-1].strip()
    while not lastLine:  # Remove any empty lines from the end
        lines.pop()
        lastLine = lines[-1].strip()
    if re_IntermediateAdjList.match(lastLine):
        logging.debug("Adjacency list line '{0}' looks like an intermediate style adjacency list".format(lastLine))
        return fromOldAdjacencyList(adjlist, group=group, saturateH=saturateH)
    if re_OldAdjList.match(lastLine):
        logging.debug("Adjacency list line '{0}' looks like an old style adjacency list".format(lastLine))
        if not group:
            logging.debug("Will assume implicit H atoms")
        return fromOldAdjacencyList(adjlist, group=group, saturateH=(not group))

    # Interpret the first line if it contains a label
    if len(lines[0].split()) == 1:
        label = lines.pop(0)
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list.')
        
    # Interpret the second line if it contains a multiplicity
    if lines[0].split()[0] == 'multiplicity':
        line = lines.pop(0)
        if group:
            match = re.match('\s*multiplicity\s+\[\s*(\d(?:,\s*\d)*)\s*\]\s*$', line)
            # should match "multiplicity [1]" or " multiplicity   [ 1, 2, 3 ]" or " multiplicity [1,2,3]"
            # and whatever's inside the [] (excluding leading and trailing spaces) should be captured as group 1.
            # Multiplicities must be only one digit (i.e. less than 10)
            # The (?:,\s*\d)* matches patters like ", 2" 0 or more times, but doesn't capture them (because of the leading ?:)
            assert match, "Invalid multiplicity line '{0}'. Should be a list like 'multiplicity [1,2,3]'".format(line)
            multiplicities = match.group(1).split(',')
            multiplicity = [int(i) for i in multiplicities]
        else:
            match = re.match('\s*multiplicity\s+\d+\s*$', line)
            assert match, "Invalid multiplicity line '{0}'. Should be an integer like 'multiplicity 2'".format(line)
            multiplicity = int(line.split()[1])
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list.')
    
    mistake1 = re.compile('\{[^}]*\s+[^}]*\}')
    # Iterate over the remaining lines, generating Atom or GroupAtom objects
    for line in lines:

        # Sometimes people put spaces after commas, which messes up the
        # parse-by-whitespace. Examples include '[Cd, Ct]'.
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
        if atomType[0] == '[':
            if not group:
                raise InvalidAdjacencyListError("A molecule should not assign more than one atomtype per atom.")
            atomType = atomType[1:-1].split(',')
        else:
            atomType = [atomType]
        index += 1
        
        # Next the number of unpaired electrons
        unpairedElectrons = []
        uState = data[index]
        if uState[0] == 'u':
            if uState[1] == '[':
                uState = uState[2:-1].split(',')
            else:
                uState = [uState[1]]
            for u in uState:
                if u == '0':
                    unpairedElectrons.append(0)
                elif u == '1':
                    unpairedElectrons.append(1)
                elif u == '2':
                    unpairedElectrons.append(2)
                elif u == '3':
                    unpairedElectrons.append(3)
                elif u == '4':
                    unpairedElectrons.append(4)
                elif u == 'x':
                    if not group:
                        raise InvalidAdjacencyListError("A molecule should not assign a wildcard to number of unpaired electrons.")
                else:
                    raise InvalidAdjacencyListError('Number of unpaired electrons not recognized.')
            index += 1
        else:
            raise InvalidAdjacencyListError('Number of unpaired electrons not defined.')
        
        # Next the number of lone electron pairs (if provided)
        lonePairs = []
        if len(data) > index:
            lpState = data[index]
            if lpState[0] == 'p':
                if lpState[1] == '[':
                    lpState = lpState[2:-1].split(',')
                else:
                    lpState = [lpState[1]]
                for l in lpState:
                    if l == '0':
                        lonePairs.append(0)
                    elif l == '1':
                        lonePairs.append(1)
                    elif l == '2':
                        lonePairs.append(2)
                    elif l == '3':
                        lonePairs.append(3)
                    elif l == '4':
                        lonePairs.append(4)
                    elif l == 'x':
                        if not group:
                            raise InvalidAdjacencyListError("A molecule should not have a wildcard assigned to number of lone pairs.")
                    else:
                        raise InvalidAdjacencyListError('Number of lone electron pairs not recognized.')
                index += 1
            else:
                if not group:
                    lonePairs.append(0)
        else:
            if not group:
                lonePairs.append(0)
            
        # Next the number of partial charges (if provided)
        partialCharges = []
        if len(data) > index:
            eState = data[index]
            if eState[0] == 'c':
                if eState[1] == '[':
                    eState = eState[2:-1].split(',')
                else:
                    eState = [eState[1:]]
                for e in eState:
                    if e == '0':
                        partialCharges.append(0)
                    elif e == '+1':
                        partialCharges.append(1)
                    elif e == '+2':
                        partialCharges.append(2)
                    elif e == '+3':
                        partialCharges.append(3)
                    elif e == '+4':
                        partialCharges.append(4)
                    elif e == '-1':
                        partialCharges.append(-1)
                    elif e == '-2':
                        partialCharges.append(-2)
                    elif e == '-3':
                        partialCharges.append(-3)
                    elif e == '-4':
                        partialCharges.append(-4)
                    elif e == 'x':
                        if not group:
                            raise InvalidAdjacencyListError("A molecule should not have a wildcard assigned to number of charges.")
                    else:
                        raise InvalidAdjacencyListError('Number of partial charges not recognized.')
                index += 1
            else:
                if not group:
                    partialCharges.append(0)
        else:
            if not group:
                partialCharges.append(0)
        
        # Create a new atom based on the above information
        if group:
            atom = GroupAtom(atomType, unpairedElectrons, partialCharges, label, lonePairs)
        else:
            atom = Atom(atomType[0], unpairedElectrons[0], partialCharges[0], label, lonePairs[0])

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
            
            if order[0] == '[':
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
            Saturator.saturate(atoms)

    
    # Consistency checks
    if not group:
        # Molecule consistency check
        # Electron and valency consistency check for each atom
        for atom in atoms: ConsistencyChecker.check_partial_charge(atom)

        nRad = sum([atom.radicalElectrons for atom in atoms])
        absolute_spin_per_electron = 1/2.
        if multiplicity == None: multiplicity = 2* (nRad * absolute_spin_per_electron) + 1
            
        ConsistencyChecker.check_multiplicity(nRad, multiplicity)
        for atom in atoms: ConsistencyChecker.check_hund_rule(atom, multiplicity)
        return atoms, multiplicity
    else:
        # Currently no group consistency check
        return atoms, multiplicity


def toAdjacencyList(atoms, multiplicity, label=None, group=False, removeH=False, removeLonePairs=False):
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
    
    if group:
        if multiplicity is not None:
            assert isinstance(multiplicity, list), "Functional group should have a list of possible multiplicities"
            if multiplicity != [1,2,3,4,5]:
                adjlist += 'multiplicity [{0!s}]\n'.format(','.join(str(i) for i in multiplicity))
    else:
        assert isinstance(multiplicity, int), "Molecule should have an integer multiplicity"
        if multiplicity != 1 or any( atom.radicalElectrons for atom in atoms ):
            adjlist += 'multiplicity {0!r}\n'.format(multiplicity)

    # Determine the numbers to use for each atom
    atomNumbers = {}
    index = 0
    for atom in atoms:
        if removeH and atom.element.symbol == 'H' and atom.label == '': continue
        atomNumbers[atom] = '{0:d}'.format(index + 1)
        index += 1
    
    atomLabels = dict([(atom, '{0}'.format(atom.label)) for atom in atomNumbers])
    
    atomTypes = {}
    atomUnpairedElectrons = {}
    atomLonePairs = {}
    atomCharge = {}
    if group:
        for atom in atomNumbers:
            # Atom type(s)
            if len(atom.atomType) == 1: 
                atomTypes[atom] = atom.atomType[0].label
            else:
                atomTypes[atom] = '[{0}]'.format(','.join([a.label for a in atom.atomType]))
            # Unpaired Electron(s)
            if len(atom.radicalElectrons) == 1: 
                atomUnpairedElectrons[atom] = str(atom.radicalElectrons[0])
            elif len(atom.radicalElectrons) == 0:
                atomUnpairedElectrons[atom] = 'x'  # the wildcard represents [0,1,2,3,4]
            else:
                atomUnpairedElectrons[atom] = '[{0}]'.format(','.join([str(radical) for radical in atom.radicalElectrons]))

            # Lone Electron Pair(s)
            if atom.lonePairs is None or not atom.lonePairs: # if None or []
                atomLonePairs[atom] = None
            elif len(atom.lonePairs) == 1: 
                atomLonePairs[atom] = str(atom.lonePairs[0])
            elif set(atom.lonePairs) == set(range(5)):
                atomLonePairs[atom] = 'x'  # the wildcard represents [0,1,2,3,4]
            else:
                atomLonePairs[atom] = '[{0}]'.format(','.join([str(pair) for pair in atom.lonePairs]))
    else:
        for atom in atomNumbers:
            # Atom type
            atomTypes[atom] = '{0}'.format(atom.element.symbol)
            # Unpaired Electron(s)
            atomUnpairedElectrons[atom] = '{0}'.format(atom.radicalElectrons)
            # Lone Electron Pair(s)
            atomLonePairs[atom] = atom.lonePairs
            # Partial Charge(s)
            atomCharge[atom] = atom.charge
    
    # Determine field widths
    atomNumberWidth = max([len(s) for s in atomNumbers.values()]) + 1
    atomLabelWidth = max([len(s) for s in atomLabels.values()])
    if atomLabelWidth > 0: atomLabelWidth += 1
    atomTypeWidth = max([len(s) for s in atomTypes.values()]) + 1
    atomUnpairedElectronsWidth = max([len(s) for s in atomUnpairedElectrons.values()])
    atomLonePairWidth = 1
    atomChargeWidth = 1
    
    # Assemble the adjacency list
    for atom in atoms:
        if atom not in atomNumbers: continue

        # Atom number
        adjlist += '{0:<{1:d}}'.format(atomNumbers[atom], atomNumberWidth)
        # Atom label
        adjlist += '{0:<{1:d}}'.format(atomLabels[atom], atomLabelWidth)
        # Atom type(s)
        adjlist += '{0:<{1:d}}'.format(atomTypes[atom], atomTypeWidth)
        # Unpaired Electron(s)
        adjlist += 'u{0:<{1:d}}'.format(atomUnpairedElectrons[atom], atomUnpairedElectronsWidth)
        # Lone Electron Pair(s)
        if str(atomLonePairs[atom]) != 'None':
            adjlist += ' p{0:>{1:d}}'.format(atomLonePairs[atom], atomLonePairWidth)
        if not group:
            # Partial Charge(s)
            if atomCharge[atom] > 0:
                adjlist += ' c+{0:>{1:d}}'.format(atomCharge[atom], atomChargeWidth)
            elif atomCharge[atom] < 0:
                adjlist += ' c{0:>{1:d}}'.format(atomCharge[atom], atomChargeWidth)
            else:
                adjlist += ' c{0:>{1:d}} '.format(atomCharge[atom], atomChargeWidth)
        
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
                    adjlist += '[{0}]'.format(','.join(bond.order))
            else:
                adjlist += bond.order
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist
