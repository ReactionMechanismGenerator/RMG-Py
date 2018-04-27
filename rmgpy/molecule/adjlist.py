#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains functionality for reading from and writing to the
adjacency list format used by Reaction Mechanism Generator (RMG).
"""
import logging
import warnings
import re
import numpy as np
from .molecule import Atom, Bond, getAtomType
from .group import GroupAtom, GroupBond
from .element import getElement, PeriodicSystem
from rmgpy.exceptions import InvalidAdjacencyListError

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
            newAtoms = []
            for atom in atoms:
                try:
                    max_number_of_valence_electrons = PeriodicSystem.valence_electrons[atom.symbol]
                except KeyError:
                    raise InvalidAdjacencyListError('Cannot add hydrogens to adjacency list: Unknown orbital for atom "{0}".'.format(atom.symbol))
                
                order = atom.getBondOrdersForAtom()
                    
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
            if atom.symbol == 'X':
                return  # because we can't check it.
        
            valence = PeriodicSystem.valence_electrons[atom.symbol]
            order = atom.getBondOrdersForAtom()
                
            theoretical = valence - order - atom.radicalElectrons - 2*atom.lonePairs

            if not (-0.301 < atom.charge - theoretical < 0.301):
                # It should be 0, but -0.1 is caused by a Hydrogen bond
                raise InvalidAdjacencyListError(
                    ('Invalid valency for atom {symbol} ({type}) with {radicals} unpaired electrons, '
                    '{lonePairs} pairs of electrons, {charge} charge, and bonds [{bonds}].'
                    ).format(symbol=atom.symbol,
                             type=getAtomType(atom, atom.edges).label,
                             radicals=atom.radicalElectrons,
                             lonePairs=atom.lonePairs,
                             charge=atom.charge,
                             bonds=','.join([str(bond.order) for bond in atom.bonds.values()])
                            ))

    @staticmethod
    def check_multiplicity(nRad, multiplicity):
        '''
        Check that the multiplicity complies with the formula: m = 2s + 1,
        where s is the sum of the spin [+/- (1/2) ] of the unpaired electrons
        
        For a simple radical (nRad = 1): 
        s = +1/2 , m = 2 (doublet)
        
        For a biradical, s can be either 0 [+0.5 + (-0.5) ] or 1 [+0.5 + (+0.5) ] 
        and m = 1 (singlet) or m = 3 (triplet).
        '''
        if nRad in [0,1]:
            if multiplicity != (nRad + 1):
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of radicals {1}.'.format(multiplicity, nRad))
        elif nRad == 2:
            if not int(multiplicity) in [1,3]: raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of radicals {1}.'.format(multiplicity, nRad))
        elif nRad == 3:
            if not int(multiplicity) in [4,2]: raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of radicals {1}.'.format(multiplicity, nRad))
        elif nRad == 4:
            if not int(multiplicity) in [5,3,1]: raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of radicals {1}.'.format(multiplicity, nRad))
        else: logging.warning("Consistency checking of multiplicity of molecules with more than 4 unpaired electrons is not implemented yet!")
    
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
                raise InvalidAdjacencyListError("""Error in adjacency list\n{0}\nNo atoms specified.""".format(adjlist))
        
        mistake1 = re.compile('\{[^}]*\s+[^}]*\}')
        atomicMultiplicities = {} # these are no longer stored on atoms, so we make a separate dictionary
        # Iterate over the remaining lines, generating Atom or GroupAtom objects
        for line in lines:

            # Sometimes people put spaces after commas, which messes up the
            # parse-by-whitespace. Examples include '{Cd, Ct}'.
            if mistake1.search(line):
                raise InvalidAdjacencyListError(
                    "Error in adjacency list: \n{1}\nspecies shouldn't have spaces inside braces: {0}".format(mistake1.search(line).group(), adjlist)
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
                raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nThere must be some electronic state defined for an old adjlist".format(adjlist))
            for e in elecState:
                if e == '0':
                    radicalElectrons.append(0); additionalLonePairs.append(0)
                elif e == '1':
                    radicalElectrons.append(1); additionalLonePairs.append(0)
                elif e == '2': 
                    if not group:
                        raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nNumber of radical electrons = 2 is not specific enough.  Please use 2S or 2T.".format(adjlist))
                    # includes 2S and 2T
                    radicalElectrons.append(0); additionalLonePairs.append(1)
                    radicalElectrons.append(2); additionalLonePairs.append(0)
                elif e == '2S':
                    radicalElectrons.append(0); additionalLonePairs.append(1)
                elif e == '2T':
                    radicalElectrons.append(2); additionalLonePairs.append(0)
                elif e == '3':
                    if not group:
                        raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nNumber of radical electrons = 3 is not specific enough.  Please use 3D or 3Q.".format(adjlist))
                    # includes 3D and 3Q
                    radicalElectrons.append(1); additionalLonePairs.append(1)
                    radicalElectrons.append(3); additionalLonePairs.append(0)
                elif e == '3D':
                    radicalElectrons.append(1); additionalLonePairs.append(1)
                elif e == '3Q':
                    radicalElectrons.append(3); additionalLonePairs.append(0)
                elif e == '4':
                    if not group:
                        raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nNumber of radical electrons = 4 is not specific enough. Please use 4S, 4T, or 4V.".format(adjlist))
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
                        raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nNumber of radical electrons = X is not specific enough.  Wildcards should only be used for groups.".format(adjlist))
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
                if lonePairsOfElectrons is not None:
                    # Intermediate adjlist representation
                    lonePairsOfElectrons = lonePairsOfElectrons + additionalLonePairs[0]
                else:
                    # Add the standard number of lone pairs with the additional lone pairs
                    lonePairsOfElectrons = PeriodicSystem.lone_pairs[atomType[0]] + additionalLonePairs[0]
                    
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
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAttempted to create a bond between atom {0:d} and itself.'.format(aid,adjlist))
                
                if order[0] == '{':
                    order = order[1:-1].split(',')
                else:
                    order = [order]

                bonds[aid][aid2] = order

        # Check consistency using bonddict
        for atom1 in bonds:
            for atom2 in bonds[atom1]:
                if atom2 not in bonds:
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAtom {0:d} not in bond dictionary.'.format(atom2,adjlist))
                elif atom1 not in bonds[atom2]:
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{2}\nFound bond between {0:d} and {1:d}, but not the reverse'.format(atom1, atom2, adjlist))
                elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                    raise InvalidAdjacencyListError('Error in adjacency list: \n{4}\nFound bonds between {0:d} and {1:d}, but of different orders "{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1], adjlist))

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
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nMultiple bond orders specified for an atom.'.format(adjlist))
                    atom1.edges[atom2] = bond
                    atom2.edges[atom1] = bond
        
        if not group:
            if saturateH:
                # Add explicit hydrogen atoms to complete structure if desired
                newAtoms = []
                for atom in atoms:
                    try:
                        valence = PeriodicSystem.valences[atom.symbol]
                    except KeyError:
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nCannot add hydrogens: Unknown valence for atom "{0}".'.format(atom.symbol, adjlist))
                    radical = atom.radicalElectrons
                    order = atom.getBondOrdersForAtom()
                    count = valence - radical - int(order) - 2*(atom.lonePairs-PeriodicSystem.lone_pairs[atom.symbol])
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
        logging.debug("adjacency list:\n{1}\nline '{0}' looks like an intermediate style adjacency list".format(lastLine, adjlist))
        return fromOldAdjacencyList(adjlist, group=group, saturateH=saturateH)
    if re_OldAdjList.match(lastLine):
        logging.debug("Adjacency list:\n{1}\nline '{0}' looks like an old style adjacency list".format(lastLine, adjlist))
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
            if not match:
                rematch = re.match('\s*multiplicity\s+x\s*$', line)
                assert rematch, "Invalid multiplicity line '{0}'. Should be a list like 'multiplicity [1,2,3]' or a wildcard 'multiplicity x'".format(line)
            else:
            # should match "multiplicity [1]" or " multiplicity   [ 1, 2, 3 ]" or " multiplicity [1,2,3]"
            # and whatever's inside the [] (excluding leading and trailing spaces) should be captured as group 1.
            # If a wildcard is desired, this line can be omitted or replaced with 'multiplicity x'
            # Multiplicities must be only one digit (i.e. less than 10)
            # The (?:,\s*\d)* matches patters like ", 2" 0 or more times, but doesn't capture them (because of the leading ?:)            
                multiplicities = match.group(1).split(',')
                multiplicity = [int(i) for i in multiplicities]
        else:
            match = re.match('\s*multiplicity\s+\d+\s*$', line)
            assert match, "Invalid multiplicity line '{0}'. Should be an integer like 'multiplicity 2'".format(line)
            multiplicity = int(line.split()[1])
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list: \n{0}'.format(adjlist))
    
    mistake1 = re.compile('\{[^}]*\s+[^}]*\}')
    # Iterate over the remaining lines, generating Atom or GroupAtom objects
    for line in lines:

        # Sometimes people put spaces after commas, which messes up the
        # parse-by-whitespace. Examples include '[Cd, Ct]'.
        if mistake1.search(line):
            raise InvalidAdjacencyListError(
                "{1} Shouldn't have spaces inside braces:\n{0}".format(mistake1.search(line).group(), adjlist)
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
                raise InvalidAdjacencyListError("Error on:\n{0}\nA molecule should not assign more than one atomtype per atom.".format(adjlist))
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
                        raise InvalidAdjacencyListError("Error on:\n{0}\nA molecule should not assign a wildcard to number of unpaired electrons.".format(adjlist))
                else:
                    raise InvalidAdjacencyListError('Number of unpaired electrons not recognized on\n{0}.'.format(adjlist))
            index += 1
        else:
            raise InvalidAdjacencyListError('Number of unpaired electrons not defined on\n{0}.'.format(adjlist))
        
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
                            raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nA molecule should not have a wildcard assigned to number of lone pairs.".format(adjlist))
                    else:
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nNumber of lone electron pairs not recognized.'.format(adjlist))
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
                            raise InvalidAdjacencyListError("Error on adjacency list:\n{0}\nA molecule should not have a wildcard assigned to number of charges.".format(adjlist))
                    else:
                        raise InvalidAdjacencyListError('Error on adjacency list:\n{0}\nNumber of partial charges not recognized.'.format(adjlist))
                index += 1
            else:
                if not group:
                    partialCharges.append(0)
        else:
            if not group:
                partialCharges.append(0)
        

        # Next the isotope (if provided)
        isotope = -1
        if len(data) > index:
            iState = data[index]
            if iState[0] == 'i':
                isotope = int(iState[1:])
                index += 1

        # Next ring membership info (if provided)
        props = {}
        if len(data) > index:
            rState = data[index]
            if rState[0] == 'r':
                props['inRing'] = bool(int(rState[1]))
                index += 1

        # Create a new atom based on the above information
        if group:
            atom = GroupAtom(atomType, unpairedElectrons, partialCharges, label, lonePairs, props)
        else:
            atom = Atom(atomType[0], unpairedElectrons[0], partialCharges[0], label, lonePairs[0])
            if isotope != -1:
                atom.element = getElement(atom.number, isotope)

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
                raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAttempted to create a bond between atom {0:d} and itself.'.format(aid, adjlist))
            
            if order[0] == '[':
                order = order[1:-1].split(',')
            else:
                order = [order]

            bonds[aid][aid2] = order

    # Check consistency using bonddict
    for atom1 in bonds:
        for atom2 in bonds[atom1]:
            if atom2 not in bonds:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAtom {0:d} not in bond dictionary.'.format(atom2, adjlist))
            elif atom1 not in bonds[atom2]:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{2}\nFound bond between {0:d} and {1:d}, but not the reverse.'.format(atom1, atom2, adjlist))
            elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{4}\nFound bonds between {0:d} and {1:d}, but of different orders "{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1], adjlist))

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
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nMultiple bond orders specified for an atom in a Molecule.'.format(adjlist))
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


def toAdjacencyList(atoms, multiplicity, label=None, group=False, removeH=False, removeLonePairs=False, oldStyle=False):
    """
    Convert a chemical graph defined by a list of `atoms` into a string
    adjacency list.
    """
    if not atoms:
        return ''

    if oldStyle:
        return toOldAdjacencyList(atoms, multiplicity, label, group, removeH)
    
    adjlist = ''

    # Don't remove hydrogen atoms if the molecule consists only of hydrogen atoms
    try:
        if removeH and all([atom.element.symbol == 'H' for atom in atoms]): removeH = False
    except AttributeError:
        pass

    if label: adjlist += label + '\n'
    
    if group:
        if multiplicity:
            # Functional group should have a list of possible multiplicities.  
            # If the list is empty, then it does not need to be written
            adjlist += 'multiplicity [{0!s}]\n'.format(','.join(str(i) for i in multiplicity))
    else:
        assert isinstance(multiplicity, int), "Molecule should have an integer multiplicity"
        if multiplicity != 1 or any( atom.radicalElectrons for atom in atoms ):
            adjlist += 'multiplicity {0!r}\n'.format(multiplicity)

    # Determine the numbers to use for each atom
    atomNumbers = {}
    index = 0
    for atom in atoms:
        if removeH and atom.symbol == 'H' and atom.label == '': continue
        atomNumbers[atom] = '{0:d}'.format(index + 1)
        index += 1
    
    atomLabels = dict([(atom, '{0}'.format(atom.label)) for atom in atomNumbers])
    
    atomTypes = {}
    atomUnpairedElectrons = {}
    atomLonePairs = {}
    atomCharge = {}
    atomIsotope = {}
    atomProps = {}
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
                atomUnpairedElectrons[atom] = 'x'  # Empty list indicates wildcard
            else:
                atomUnpairedElectrons[atom] = '[{0}]'.format(','.join([str(radical) for radical in atom.radicalElectrons]))

            # Lone Electron Pair(s)
            if len(atom.lonePairs) == 1: 
                atomLonePairs[atom] = str(atom.lonePairs[0])
            elif len(atom.lonePairs) == 0:  
                atomLonePairs[atom] = None   # Empty list indicates wildcard
            else:
                atomLonePairs[atom] = '[{0}]'.format(','.join([str(pair) for pair in atom.lonePairs]))
                
            # Charges
            if len(atom.charge) == 1: 
                atomCharge[atom] = '+' + str(atom.charge[0]) if atom.charge[0] > 0 else str(atom.charge[0])
            elif len(atom.charge) == 0:  
                atomCharge[atom] = None   # Empty list indicates wildcard
            else:
                atomCharge[atom] = '[{0}]'.format(','.join(['+'+str(charge) if charge > 0 else ''+str(charge) for charge in atom.charge]))

            # Isotopes
            atomIsotope[atom] = -1

            # Other props
            props = []
            if 'inRing' in atom.props:
                props.append(' r{0}'.format(int(atom.props['inRing'])))
            atomProps[atom] = props
    else:
        for atom in atomNumbers:
            # Atom type
            atomTypes[atom] = '{0}'.format(atom.symbol)
            # Unpaired Electron(s)
            atomUnpairedElectrons[atom] = '{0}'.format(atom.radicalElectrons)
            # Lone Electron Pair(s)
            atomLonePairs[atom] = str(atom.lonePairs)
            # Partial Charge(s)
            atomCharge[atom] = '+'+str(atom.charge) if atom.charge > 0 else '' + str(atom.charge)
            # Isotopes
            if isinstance(atom, Atom):
                atomIsotope[atom] = atom.element.isotope
            else:
                # cutting labels in 
                # fragment cases
                atomIsotope[atom] = atom.isotope

    
    # Determine field widths
    atomNumberWidth = max([len(s) for s in atomNumbers.values()]) + 1
    atomLabelWidth = max([len(s) for s in atomLabels.values()])
    if atomLabelWidth > 0: atomLabelWidth += 1
    atomTypeWidth = max([len(s) for s in atomTypes.values()]) + 1
    atomUnpairedElectronsWidth = max([len(s) for s in atomUnpairedElectrons.values()])
    #atomLonePairWidth = max([len(s) for s in atomLonePairs.values()])
    #atomChargeWidth = max([len(s) for s in atomCharge.values()])
    
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
        if atomLonePairs[atom] != None:
            adjlist += ' p{0}'.format(atomLonePairs[atom])
        # Partial charges
        if atomCharge[atom] != None:
            adjlist += ' c{0}'.format(atomCharge[atom])
        # Isotopes
        if atomIsotope[atom] != -1:
            adjlist += ' i{0}'.format(atomIsotope[atom])
        if group and len(atomProps[atom]) > 0:
            for prop in atomProps[atom]:
                adjlist += prop

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
                code = '[{0}]'
                if len(bond.order) == 1:
                    code = '{0}'
                # preference is for string representation, backs down to number
                # numbers if doesn't work
                try:
                    adjlist += code.format(','.join(bond.getOrderStr()))
                except ValueError:
                    adjlist += code.format(','.join(str(bond.getOrderNum())))
            else:
                # preference is for string representation, backs down to number
                # numbers if doesn't work
                try:
                    adjlist += bond.getOrderStr()
                except ValueError:
                    adjlist += str(bond.getOrderNum())
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist

def getOldElectronState(atom):
    """
    Get the old adjacency list format electronic state
    """
    additionalLonePairs = atom.lonePairs - PeriodicSystem.lone_pairs[atom.element.symbol]
    electrons = atom.radicalElectrons + additionalLonePairs * 2
    if electrons == 0:
        electronState = '0'
    elif electrons == 1:
        electronState = '1'
    elif electrons == 2:
        if additionalLonePairs == 0:
            electronState = '2T'
        elif additionalLonePairs == 1:
            electronState = '2S'
        else:
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    elif electrons == 3:
        if additionalLonePairs == 0: 
            electronState = '3Q'
        elif additionalLonePairs == 1:
            electronState = '3D'
        else: 
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    elif electrons == 4:
        if additionalLonePairs == 0:
            electronState = '4V'
        elif additionalLonePairs == 1:
            electronState = '4T'
        elif additionalLonePairs == 2:
            electronState = '4S'
        else: 
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    else: 
        raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    return electronState



def toOldAdjacencyList(atoms, multiplicity=None, label=None, group=False, removeH=False):
    """
    Convert a chemical graph defined by a list of `atoms` into a string old-style 
    adjacency list that can be used in RMG-Java.  Currently not working for groups.
    """
    warnings.warn("The old adjacency lists are no longer supported and may be"
                  " removed in version 2.3.", DeprecationWarning)
    adjlist = ''
    
    if group:
        raise InvalidAdjacencyListError("Not yet implemented.")
    # Filter out all non-valid atoms
    if not group:
        for atom in atoms:
            if atom.element.symbol in ['He','Ne','Ar','N']:
                raise InvalidAdjacencyListError("Old-style adjacency list does not accept He, Ne, Ar, N elements.")

    # Don't remove hydrogen atoms if the molecule consists only of hydrogen atoms
    try:
        if removeH and all([atom.element.symbol == 'H' for atom in atoms]): removeH = False
    except AttributeError:
        pass

    if label: adjlist += label + '\n'

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
    if group:
        raise InvalidAdjacencyListError("Not yet implemented.")
    else:
        for atom in atomNumbers:
            # Atom type
            atomTypes[atom] = '{0}'.format(atom.element.symbol)
            # Electron state(s)
            atomElectronStates[atom] = '{0}'.format(getOldElectronState(atom))    
    
    # Determine field widths
    atomNumberWidth = max([len(s) for s in atomNumbers.values()]) + 1
    atomLabelWidth = max([len(s) for s in atomLabels.values()])
    if atomLabelWidth > 0: atomLabelWidth += 1
    atomTypeWidth = max([len(s) for s in atomTypes.values()]) + 1
    atomElectronStateWidth = max([len(s) for s in atomElectronStates.values()])
    
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
                    adjlist += bond.getOrderStr()[0]
                else:
                    adjlist += '{{{0}}}'.format(','.join(bond.getOrderStr()))
            else:
                adjlist += bond.getOrderStr()
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist
