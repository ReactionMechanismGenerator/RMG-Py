#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
import re
import warnings

from rmgpy.exceptions import InvalidAdjacencyListError
from rmgpy.molecule.atomtype import get_atomtype
from rmgpy.molecule.element import get_element, PeriodicSystem
from rmgpy.molecule.group import GroupAtom, GroupBond
from rmgpy.molecule.molecule import Atom, Bond
from rmgpy.molecule.fragment import Fragment, CuttingLabel


class Saturator(object):
    @staticmethod
    def saturate(atoms):
        """
        Returns a list of atoms that is extended
        (and bond attributes) by saturating the valency of the non-hydrogen atoms with an
        appropriate number of hydrogen atoms.

        The required number of hydrogen atoms per heavy atom is determined as follows:
        H's =     max number of valence electrons - atom.radical_electrons
                    - 2* atom.lone_pairs - order - atom.charge

        """
        new_atoms = []
        for atom in atoms:
            if not isinstance(atom, Atom): continue
            try:
                max_number_of_valence_electrons = PeriodicSystem.valence_electrons[atom.symbol]
            except KeyError:
                raise InvalidAdjacencyListError(
                    'Cannot add hydrogens to adjacency list: Unknown orbital for atom "{0}".'.format(atom.symbol))

            order = atom.get_total_bond_order()

            number_of_h_to_be_added = max_number_of_valence_electrons - atom.radical_electrons - 2 * atom.lone_pairs - int(
                order) - atom.charge

            if number_of_h_to_be_added < 0:
                raise InvalidAdjacencyListError('Incorrect electron configuration on atom.')

            for _ in range(number_of_h_to_be_added):
                a = Atom(element='H', radical_electrons=0, charge=0, label='', lone_pairs=0)
                b = Bond(atom, a, 'S')
                new_atoms.append(a)
                atom.bonds[a] = b
                a.bonds[atom] = b
        atoms.extend(new_atoms)


class ConsistencyChecker(object):

    @staticmethod
    def check_partial_charge(atom):
        """
        Checks whether the partial charge attribute of the atom checks out with
        the theoretical one:

        """
        if atom.symbol in {'X','L','R'}:
            return  # because we can't check it.

        valence = PeriodicSystem.valence_electrons[atom.symbol]
        order = atom.get_total_bond_order()

        theoretical = valence - order - atom.radical_electrons - 2 * atom.lone_pairs

        if not (-0.301 < atom.charge - theoretical < 0.301):
            # It should be 0, but -0.1 is caused by a Hydrogen bond
            raise InvalidAdjacencyListError(
                'Invalid valency for atom {symbol} ({type}) with {radicals} unpaired electrons, '
                '{lone_pairs} pairs of electrons, {charge} charge, and bonds [{bonds}].'.format(
                    symbol=atom.symbol,
                    type=get_atomtype(atom, atom.edges).label,
                    radicals=atom.radical_electrons,
                    lone_pairs=atom.lone_pairs,
                    charge=atom.charge,
                    bonds=','.join([str(bond.order) for bond in atom.bonds.values()])
                )
            )

    @staticmethod
    def check_multiplicity(n_rad, multiplicity):
        """
        Check that the multiplicity complies with the formula: m = 2s + 1,
        where s is the sum of the spin [+/- (1/2) ] of the unpaired electrons

        For a simple radical (n_rad = 1):
        s = +1/2 , m = 2 (doublet)

        For a biradical, s can be either 0 [+0.5 + (-0.5) ] or 1 [+0.5 + (+0.5) ]
        and m = 1 (singlet) or m = 3 (triplet).
        """
        if n_rad in [0, 1]:
            if multiplicity != (n_rad + 1):
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        elif n_rad == 2:
            if not int(multiplicity) in [1, 3]:
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        elif n_rad == 3:
            if not int(multiplicity) in [4, 2]:
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        elif n_rad == 4:
            if not int(multiplicity) in [5, 3, 1]:
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        else:
            logging.warning("Consistency checking of multiplicity of molecules with "
                            "more than 4 unpaired electrons is not implemented yet!")

    @staticmethod
    def check_hund_rule(atom, multiplicity):
        """
        It is checked whether atoms with 2 unpaired electrons on the same atom
        result in a multiplicity of 3, and not 1.

        Unpaired electrons in 2 different orbitals belonging to the same atom
        should have the same spin, and hence, should result in a multiplicity of 3.
        """
        if atom.radical_electrons == 2 and multiplicity == 1:
            raise InvalidAdjacencyListError(
                "Violation of hund's rule. Invalid multiplicity of {0} because there is an "
                "atom with {1} unpaired electrons".format(multiplicity, atom.radical_electrons))


################################################################################

def from_old_adjacency_list(adjlist, group=False, saturate_h=False):
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

        mistake1 = re.compile(r'\{[^}]*\s+[^}]*\}')
        atomic_multiplicities = {}  # these are no longer stored on atoms, so we make a separate dictionary
        # Iterate over the remaining lines, generating Atom or GroupAtom objects
        for line in lines:

            # Sometimes people put spaces after commas, which messes up the
            # parse-by-whitespace. Examples include '{Cd, Ct}'.
            if mistake1.search(line):
                raise InvalidAdjacencyListError("Error in adjacency list: \n{1}\nspecies shouldn't have spaces inside "
                                                "braces: {0}".format(mistake1.search(line).group(), adjlist))

            # Sometimes commas are used to delimit bonds in the bond list,
            # so replace them just in case
            line = line.replace('},{', '} {')

            data = line.split()

            # Skip if blank line
            if len(data) == 0:
                continue

            # First item is index for atom
            # Sometimes these have a trailing period (as if in a numbered list),
            # so remove it just in case
            aid = int(data[0].strip('.'))

            # If second item starts with '*', then atom is labeled
            label = ''
            index = 1
            if data[1][0] == '*':
                label = data[1]
                index += 1

            # Next is the element or functional group element
            # A list can be specified with the {,} syntax
            atom_type = data[index]
            if atom_type[0] == '{':
                atom_type = atom_type[1:-1].split(',')
            else:
                atom_type = [atom_type]
            index += 1

            # Next is the electron state
            radical_electrons = []
            additional_lone_pairs = []
            elec_state = data[index].upper()
            if elec_state[0] == '{':
                elec_state = elec_state[1:-1].split(',')
            else:
                elec_state = [elec_state]
            if len(elec_state) == 0:
                raise InvalidAdjacencyListError(
                    "Error in adjacency list:\n{0}\nThere must be some electronic state defined for an "
                    "old adjlist".format(adjlist))
            for e in elec_state:
                if e == '0':
                    radical_electrons.append(0)
                    additional_lone_pairs.append(0)
                elif e == '1':
                    radical_electrons.append(1)
                    additional_lone_pairs.append(0)
                elif e == '2': 
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = 2 is not specific enough. "
                            "Please use 2S or 2T.".format(adjlist))
                    # includes 2S and 2T
                    radical_electrons.append(0); additional_lone_pairs.append(1)
                    radical_electrons.append(2); additional_lone_pairs.append(0)
                elif e == '2S':
                    radical_electrons.append(0); additional_lone_pairs.append(1)
                elif e == '2T':
                    radical_electrons.append(2); additional_lone_pairs.append(0)
                elif e == '3':
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = 3 is not specific enough. "
                            "Please use 3D or 3Q.".format(adjlist))
                    # includes 3D and 3Q
                    radical_electrons.append(1); additional_lone_pairs.append(1)
                    radical_electrons.append(3); additional_lone_pairs.append(0)
                elif e == '3D':
                    radical_electrons.append(1); additional_lone_pairs.append(1)
                elif e == '3Q':
                    radical_electrons.append(3); additional_lone_pairs.append(0)
                elif e == '4':
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = 4 is not specific enough. "
                            "Please use 4S, 4T, or 4V.".format(adjlist))
                    # includes 4S, 4T, and 4V
                    radical_electrons.append(0); additional_lone_pairs.append(2)
                    radical_electrons.append(2); additional_lone_pairs.append(1)
                    radical_electrons.append(4); additional_lone_pairs.append(0)
                elif e == '4S':
                    radical_electrons.append(0); additional_lone_pairs.append(2)
                elif e == '4T':
                    radical_electrons.append(2); additional_lone_pairs.append(1)
                elif e == '4V':
                    radical_electrons.append(4); additional_lone_pairs.append(0)
                elif e == 'X':
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = X is not specific enough. "
                            "Wildcards should only be used for groups.".format(adjlist))
                    radical_electrons = []
            index += 1

            # Next number defines the number of lone electron pairs (if provided)
            lone_pairs_of_electrons = None
            if len(data) > index:
                lp_state = data[index]
                if lp_state[0] == '{':
                    # this is the start of the chemical bonds - no lone pair info was provided
                    lone_pairs_of_electrons = None
                else:
                    if lp_state == '0':
                        lone_pairs_of_electrons = 0
                    if lp_state == '1':
                        lone_pairs_of_electrons = 1
                    if lp_state == '2':
                        lone_pairs_of_electrons = 2
                    if lp_state == '3':
                        lone_pairs_of_electrons = 3
                    if lp_state == '4':
                        lone_pairs_of_electrons = 4
                    index += 1
            else:  # no bonds or lone pair info provided.
                lone_pairs_of_electrons = None

            # Create a new atom based on the above information
            if group:
                if lone_pairs_of_electrons is not None:
                    lone_pairs_of_electrons = [additional + lone_pairs_of_electrons for additional in additional_lone_pairs]
                atom = GroupAtom(atomtype=atom_type,
                                 radical_electrons=sorted(set(radical_electrons)),
                                 charge=None,
                                 label=label,
                                 lone_pairs=lone_pairs_of_electrons,
                                 # Assign lone_pairs_of_electrons as None if it is not explicitly provided
                                 )

            else:
                if lone_pairs_of_electrons is not None:
                    # Intermediate adjlist representation
                    lone_pairs_of_electrons = lone_pairs_of_electrons + additional_lone_pairs[0]
                else:
                    # Add the standard number of lone pairs with the additional lone pairs
                    lone_pairs_of_electrons = PeriodicSystem.lone_pairs[atom_type[0]] + additional_lone_pairs[0]

                atom = Atom(element=atom_type[0],
                            radical_electrons=radical_electrons[0],
                            charge=0,
                            label=label,
                            lone_pairs=lone_pairs_of_electrons,
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
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAttempted to create a bond between '
                                                    'atom {0:d} and itself.'.format(aid, adjlist))

                if order[0] == '{':
                    order = order[1:-1].split(',')
                else:
                    order = [order]

                bonds[aid][aid2] = order

        # Check consistency using bonddict
        for atom1 in bonds:
            for atom2 in bonds[atom1]:
                if atom2 not in bonds:
                    raise InvalidAdjacencyListError(
                        'Error in adjacency list:\n{1}\nAtom {0:d} not in bond dictionary.'.format(atom2, adjlist))
                elif atom1 not in bonds[atom2]:
                    raise InvalidAdjacencyListError(
                        'Error in adjacency list:\n{2}\nFound bond between {0:d} and {1:d}, '
                        'but not the reverse'.format(atom1, atom2, adjlist))
                elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                    raise InvalidAdjacencyListError(
                        'Error in adjacency list: \n{4}\nFound bonds between {0:d} and {1:d}, but of different orders '
                        '"{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1], adjlist))

        # Convert bonddict to use Atom[group] and Bond[group] objects
        atomkeys = list(atomdict.keys())
        atomkeys.sort()
        for aid1 in atomkeys:
            atomkeys2 = list(bonds[aid1].keys())
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
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nMultiple bond orders specified '
                                                        'for an atom.'.format(adjlist))
                    atom1.edges[atom2] = bond
                    atom2.edges[atom1] = bond

        if not group:
            if saturate_h:
                # Add explicit hydrogen atoms to complete structure if desired
                new_atoms = []
                for atom in atoms:
                    try:
                        valence = PeriodicSystem.valences[atom.symbol]
                    except KeyError:
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nCannot add hydrogens: Unknown '
                                                        'valence for atom "{0}".'.format(atom.symbol, adjlist))
                    radical = atom.radical_electrons
                    order = atom.get_total_bond_order()
                    count = valence - radical - int(order) - 2 * (
                            atom.lone_pairs - PeriodicSystem.lone_pairs[atom.symbol])
                    for i in range(count):
                        a = Atom(element='H', radical_electrons=0, charge=0, label='', lone_pairs=0)
                        b = Bond(atom, a, 'S')
                        new_atoms.append(a)
                        atom.bonds[a] = b
                        a.bonds[atom] = b
                atoms.extend(new_atoms)

            # Calculate the multiplicity for the molecule and update the charges on each atom
            n_rad = 0  # total number of radical electrons
            for atom in atoms:
                atom.update_charge()
                n_rad += atom.radical_electrons
            multiplicity = n_rad + 1  # 2 s + 1, where s is the combined spin of unpaired electrons (s = 1/2 per unpaired electron)

        else:
            # Don't set a multiplicity for groups when converting from an old adjlist
            multiplicity = None

    except InvalidAdjacencyListError:
        logging.error("Troublesome adjacency list:\n" + adjlist)
        raise
    if group:
        return atoms, multiplicity, [], []
    else:
        return atoms, multiplicity, '', ''


###############################

re_intermediate_adjlist = re.compile(r'^\s*(\d*)\s+' +  # atom number digit
                                     r'(?P<label>\*\d*\s+)?' +  # optional label eg * or *2
                                     r'(?P<atomtype>\{?[A-Z]\S*)\s+' +  # atomtype eg R!H or {Cb,Cd}
                                     r'(?P<radicals>X|\d[STDQV]?|\{?\d[^}]*\})\s+' +  # radicals eg. X or 2T or {1,2,2T}
                                     r'(?P<lonepairs>\d)' +  # lone pairs eg. 0
                                     r'(?P<bonds>(\s+\{\d+\,(?:[SDTB]|\{.+?\})\},?)*)' +  # bonds, eg {2,S} {4,{S,D}}
                                     r'\s*$')  # the end!

re_old_adjlist = re.compile(r'^\s*(\d*)\s+' +  # atom number digit
                            r'(?P<label>\*\d*\s+)?' +  # optional label eg * or *2
                            r'(?P<atomtype>\{?[A-Z]\S*)\s+' +  # atomtype eg R!H or {Cb,Cd}
                            r'(?P<radicals>X|\d[STDQV]?|\{?\d[^}]*\})' +  # radicals eg. X or 2T or {1,2,2T}
                            r'(?P<bonds>(\s+\{\d+\,(?:[SDTB]|\{.+?\})\},?)*)' +  # bonds, eg {2,S} {4,{S,D}}
                            r'\s*$')  # the end!


def from_adjacency_list(adjlist, group=False, saturate_h=False, check_consistency=True):
    """
    Convert a string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects.
    """
    atoms = []
    atom_dict = {}
    bonds = {}
    multiplicity = None

    adjlist = adjlist.strip()
    lines = adjlist.splitlines()
    if adjlist == '' or len(lines) == 0:
        raise InvalidAdjacencyListError('Empty adjacency list.')

    # Detect old-style adjacency lists by looking at the last line's syntax
    last_line = lines[-1].strip()
    while not last_line:  # Remove any empty lines from the end
        lines.pop()
        last_line = lines[-1].strip()
    if re_intermediate_adjlist.match(last_line):
        logging.debug(
            "adjacency list:\n{1}\nline '{0}' looks like an intermediate style "
            "adjacency list".format(last_line, adjlist))
        return from_old_adjacency_list(adjlist, group=group, saturate_h=saturate_h)
    if re_old_adjlist.match(last_line):
        logging.debug(
            "Adjacency list:\n{1}\nline '{0}' looks like an old style adjacency list".format(last_line, adjlist))
        if not group:
            logging.debug("Will assume implicit H atoms")
        return from_old_adjacency_list(adjlist, group=group, saturate_h=(not group))

    # Interpret the first line if it contains a label
    if len(lines[0].split()) == 1:
        label = lines.pop(0)
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list.')

    # Interpret the second line if it contains a multiplicity
    if lines[0].split()[0] == 'multiplicity':
        line = lines.pop(0)
        if group:
            match = re.match(r'\s*multiplicity\s+\[\s*(\d(?:,\s*\d)*)\s*\]\s*$', line)
            if not match:
                rematch = re.match(r'\s*multiplicity\s+x\s*$', line)
                if not rematch:
                    raise InvalidAdjacencyListError("Invalid multiplicity line '{0}'. Should be a list like "
                                                    "'multiplicity [1,2,3]' or a wildcard 'multiplicity x'".format(line))
            else:
                # should match "multiplicity [1]" or " multiplicity   [ 1, 2, 3 ]" or " multiplicity [1,2,3]"
                # and whatever's inside the [] (excluding leading and trailing spaces) should be captured as group 1.
                # If a wildcard is desired, this line can be omitted or replaced with 'multiplicity x'
                # Multiplicities must be only one digit (i.e. less than 10)
                # The (?:,\s*\d)* matches patters like ", 2" 0 or more times, but doesn't capture them (because of the leading ?:)
                multiplicities = match.group(1).split(',')
                multiplicity = [int(i) for i in multiplicities]
        else:
            match = re.match(r'\s*multiplicity\s+\d+\s*$', line)
            if not match:
                raise InvalidAdjacencyListError("Invalid multiplicity line '{0}'. Should be an integer like "
                                                "'multiplicity 2'".format(line))
            multiplicity = int(line.split()[1])
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list: \n{0}'.format(adjlist))
        
    mistake1 = re.compile(r'\{[^}]*\s+[^}]*\}')
    if group:
        metal = []
        facet = []
    else:
        metal = ''
        facet = ''
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
        if len(data) == 0:
            continue
        
        if line.split()[0] == 'metal':
            if group:
                match = re.match(r'\s*metal\s+\[\s*[\w,\s]+\s*\]\s*$', line)
                if not match:
                    rematch = re.match(r'\s*metal\s+x\s*$', line)
                    if not rematch:
                        raise InvalidAdjacencyListError("Invalid metal line '{0}'. Should be a list like "
                                                        "'metal [Cu,Fe,Ag]' or a wildcard 'metal x'".format(line))
                else:
                    out = line.split('[')[1][:-1]
                    metal = [x.strip() for x in out.split(',') if x.strip() != '']
        
            else:
                match = re.match(r'\s*metal\s+\w+\s*$', line)
                if not match:
                    raise InvalidAdjacencyListError("Invalid metal line '{0}'. Should be a string like "
                                                "'metal Fe'".format(line))
                metal = line.split()[1].strip()
            
            continue
                
        if line.split()[0] == 'facet':
            if group:
                match = re.match(r'\s*facet\s+\[\s*[\w,\s]+\s*\]\s*$', line)
                if not match:
                    rematch = re.match(r'\s*facet\s+x\s*$', line)
                    if not rematch:
                        raise InvalidAdjacencyListError("Invalid facet line '{0}'. Should be a list like "
                                                        "'facet [111,211,110]' or a wildcard 'facet x'".format(line))
                else:
                    out = line.split('[')[1][:-1]
                    facet = [x.strip() for x in out.split(',') if x.strip() != '']
        
            else:
                match = re.match(r'\s*facet\s+\w+\s*$', line)
                if not match:
                    raise InvalidAdjacencyListError("Invalid facet line '{0}'. Should be a string like "
                                                "'facet 111'".format(line))
                facet = line.split()[1].strip()
            
            continue
        
        # First item is index for atom
        # Sometimes these have a trailing period (as if in a numbered list),
        # so remove it just in case
        aid = int(data[0].strip('.'))

        # If second item starts with '*', then atom is labeled
        label = ''
        index = 1
        if data[1][0] == '*':
            label = data[1]
            index += 1

        # Next is the element or functional group element
        # A list can be specified with the {,} syntax
        atom_type = data[index]
        if atom_type[0] == '[':
            if not group:
                raise InvalidAdjacencyListError("Error on:\n{0}\nA molecule should not assign more than one "
                                                "atomtype per atom.".format(adjlist))
            atom_type = atom_type[1:-1].split(',')
        else:
            atom_type = [atom_type]
        index += 1

        # Next the number of unpaired electrons
        unpaired_electrons = []
        u_state = data[index]
        if u_state[0] == 'u':
            if u_state[1] == '[':
                u_state = u_state[2:-1].split(',')
            else:
                u_state = [u_state[1]]
            for u in u_state:
                if u == '0':
                    unpaired_electrons.append(0)
                elif u == '1':
                    unpaired_electrons.append(1)
                elif u == '2':
                    unpaired_electrons.append(2)
                elif u == '3':
                    unpaired_electrons.append(3)
                elif u == '4':
                    unpaired_electrons.append(4)
                elif u == 'x':
                    if not group:
                        raise InvalidAdjacencyListError("Error on:\n{0}\nA molecule should not assign a wildcard to "
                                                        "number of unpaired electrons.".format(adjlist))
                else:
                    raise InvalidAdjacencyListError('Number of unpaired electrons not recognized on\n{0}.'.format(adjlist))
            index += 1
        else:
            raise InvalidAdjacencyListError('Number of unpaired electrons not defined on\n{0}.'.format(adjlist))

        # Next the number of lone electron pairs (if provided)
        lone_pairs = []
        if len(data) > index:
            lp_state = data[index]
            if lp_state[0] == 'p':
                if lp_state[1] == '[':
                    lp_state = lp_state[2:-1].split(',')
                else:
                    lp_state = [lp_state[1]]
                for lp in lp_state:
                    if lp == '0':
                        lone_pairs.append(0)
                    elif lp == '1':
                        lone_pairs.append(1)
                    elif lp == '2':
                        lone_pairs.append(2)
                    elif lp == '3':
                        lone_pairs.append(3)
                    elif lp == '4':
                        lone_pairs.append(4)
                    elif lp == 'x':
                        if not group:
                            raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nA molecule should not have "
                                                            "a wildcard assigned to number of lone pairs.".format(adjlist))
                    else:
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nNumber of lone electron pairs '
                                                        'not recognized.'.format(adjlist))
                index += 1
            else:
                if not group:
                    lone_pairs.append(0)
        else:
            if not group:
                lone_pairs.append(0)

        # Next the number of partial charges (if provided)
        partial_charges = []
        if len(data) > index:
            e_state = data[index]
            if e_state[0] == 'c':
                if e_state[1] == '[':
                    e_state = e_state[2:-1].split(',')
                else:
                    e_state = [e_state[1:]]
                for e in e_state:
                    if e == '0':
                        partial_charges.append(0)
                    elif e == '+1':
                        partial_charges.append(1)
                    elif e == '+2':
                        partial_charges.append(2)
                    elif e == '+3':
                        partial_charges.append(3)
                    elif e == '+4':
                        partial_charges.append(4)
                    elif e == '-1':
                        partial_charges.append(-1)
                    elif e == '-2':
                        partial_charges.append(-2)
                    elif e == '-3':
                        partial_charges.append(-3)
                    elif e == '-4':
                        partial_charges.append(-4)
                    elif e == 'x':
                        if not group:
                            raise InvalidAdjacencyListError("Error on adjacency list:\n{0}\nA molecule should not have "
                                                            "a wildcard assigned to number of charges.".format(adjlist))
                    else:
                        raise InvalidAdjacencyListError('Error on adjacency list:\n{0}\nNumber of partial charges '
                                                        'not recognized.'.format(adjlist))
                index += 1
            else:
                if not group:
                    partial_charges.append(0)
        else:
            if not group:
                partial_charges.append(0)

        # Next the sites (if provided)
        sites = []
        if len(data) > index:
            s_state = data[index]
            if s_state[0] == 's':
                if s_state[1] == '[':
                    s_state = s_state[2:-1].split(',')
                else:
                    s_state = [s_state[1:]]
                for s in s_state:
                    sites.append(s[1:-1])
                index += 1
            
        # Next the morphologys (if provided)
        morphologies = []
        if len(data) > index:
            m_state = data[index]
            if m_state[0] == 'm':
                if m_state[1] == '[':
                    m_state = m_state[2:-1].split(',')
                else:
                    m_state = [m_state[1:]]
                for m in m_state:
                    morphologies.append(m[1:-1])
                index += 1
        
        # Next the isotope (if provided)
        isotope = -1
        if len(data) > index:
            i_state = data[index]
            if i_state[0] == 'i':
                isotope = int(i_state[1:])
                index += 1

        # Next ring membership info (if provided)
        props = {}
        if len(data) > index:
            r_state = data[index]
            if r_state[0] == 'r':
                props['inRing'] = bool(int(r_state[1]))
                index += 1

        # Create a new atom based on the above information
        if group:
            atom = GroupAtom(atom_type, unpaired_electrons, partial_charges, label, lone_pairs, sites, morphologies, props)
        else:
            # detect if this is cutting label or atom
            _ , cutting_label_list = Fragment().detect_cutting_label(atom_type[0])
            if cutting_label_list == []:
                if sites == []:
                    site = ''
                else:
                    site = sites[0]
                if morphologies == []:
                    morphology = ''
                else:
                    morphology = morphologies[0]
                atom = Atom(atom_type[0], unpaired_electrons[0], partial_charges[0], label, lone_pairs[0], site, morphology)
                if isotope != -1:
                    atom.element = get_element(atom.number, isotope)
            else:
                atom = CuttingLabel(name=atom_type[0], label=label)

        # Add the atom to the list
        atoms.append(atom)
        atom_dict[aid] = atom

        # Process list of bonds
        bonds[aid] = {}
        for datum in data[index:]:

            # Sometimes commas are used to delimit bonds in the bond list,
            # so strip them just in case
            datum = datum.strip(',')

            aid2, comma, order = datum[1:-1].partition(',')
            aid2 = int(aid2)
            if aid == aid2:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAttempted to create a bond between '
                                                'atom {0:d} and itself.'.format(aid, adjlist))

            if order[0] == '[':
                order = order[1:-1].split(',')
            else:
                order = [order]

            bonds[aid][aid2] = order

    # Check consistency using bonddict
    for atom1 in bonds:
        for atom2 in bonds[atom1]:
            if atom2 not in bonds:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAtom {0:d} not in bond '
                                                'dictionary.'.format(atom2, adjlist))
            elif atom1 not in bonds[atom2]:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{2}\nFound bond between {0:d} and {1:d}, '
                                                'but not the reverse.'.format(atom1, atom2, adjlist))
            elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                raise InvalidAdjacencyListError(
                    'Error in adjacency list:\n{4}\nFound bonds between {0:d} and {1:d}, but of different orders '
                    '"{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1], adjlist))

    # Convert bonddict to use Atom[group] and Bond[group] objects
    atomkeys = list(atom_dict.keys())
    atomkeys.sort()
    for aid1 in atomkeys:
        atomkeys2 = list(bonds[aid1].keys())
        atomkeys2.sort()
        for aid2 in atomkeys2:
            if aid1 < aid2:
                atom1 = atom_dict[aid1]
                atom2 = atom_dict[aid2]
                order = bonds[aid1][aid2]
                if group:
                    bond = GroupBond(atom1, atom2, order)
                elif len(order) == 1:
                    bond = Bond(atom1, atom2, order[0])
                else:
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nMultiple bond orders specified for '
                                                    'an atom in a Molecule.'.format(adjlist))
                atom1.edges[atom2] = bond
                atom2.edges[atom1] = bond

    if saturate_h:
        # Add explicit hydrogen atoms to complete structure if desired
        if not group:
            Saturator.saturate(atoms)

    # Consistency checks
    if not group and check_consistency:
        # Molecule consistency check
        # Electron and valency consistency check for each atom
        for atom in atoms:
            if isinstance(atom, Atom):
                ConsistencyChecker.check_partial_charge(atom)

        n_rad = sum([atom.radical_electrons for atom in atoms])
        absolute_spin_per_electron = 1 / 2.
        if multiplicity is None:
            multiplicity = 2 * (n_rad * absolute_spin_per_electron) + 1

        ConsistencyChecker.check_multiplicity(n_rad, multiplicity)
        for atom in atoms:
            ConsistencyChecker.check_hund_rule(atom, multiplicity)
        return atoms, multiplicity, metal, facet
    else:
        # Currently no group consistency check
        if not group:
            if multiplicity is None:
                n_rad = sum([atom.radical_electrons for atom in atoms])
                multiplicity = n_rad + 1

        return atoms, multiplicity, metal, facet



def to_adjacency_list(atoms, multiplicity, metal='', facet='', label=None, group=False, remove_h=False, remove_lone_pairs=False,
                      old_style=False):
    """
    Convert a chemical graph defined by a list of `atoms` into a string
    adjacency list.
    """
    if old_style:
        warnings.warn("Support for writing old style adjacency lists has been removed in RMG-Py v3.", RuntimeWarning)
    if not atoms:
        return ''

    adjlist = ''

    # Don't remove hydrogen atoms if the molecule consists only of hydrogen atoms
    try:
        if remove_h and all([atom.element.symbol == 'H' for atom in atoms]): remove_h = False
    except AttributeError:
        pass

    if label:
        adjlist += label + '\n'

    if group:
        if multiplicity:
            # Functional group should have a list of possible multiplicities.  
            # If the list is empty, then it does not need to be written
            adjlist += 'multiplicity [{0!s}]\n'.format(','.join(str(i) for i in multiplicity))
        if metal:
            adjlist += 'metal [{0!s}]\n'.format(','.join(i for i in metal))
        if facet:
            adjlist += 'facet [{0!s}]\n'.format(','.join(i for i in facet))
    else:
        assert isinstance(multiplicity, int), "Molecule should have an integer multiplicity"
        if multiplicity != 1 or any(atom.radical_electrons for atom in atoms):
            adjlist += 'multiplicity {0!r}\n'.format(multiplicity)
        if metal:
            adjlist += f"metal {metal}\n"
        if facet:
            adjlist += f"facet {facet}\n"

    # Determine the numbers to use for each atom
    atom_numbers = {}
    index = 0
    for atom in atoms:
        if remove_h and atom.symbol == 'H' and atom.label == '':
            continue
        atom_numbers[atom] = '{0:d}'.format(index + 1)
        index += 1

    atom_labels = dict([(atom, '{0}'.format(atom.label)) for atom in atom_numbers])

    atom_types = {}
    atom_unpaired_electrons = {}
    atom_lone_pairs = {}
    atom_charge = {}
    atom_isotope = {}
    atom_props = {}
    atom_site = {}
    atom_morphology = {}
    if group:
        for atom in atom_numbers:
            # Atom type(s)
            if len(atom.atomtype) == 1:
                atom_types[atom] = atom.atomtype[0].label
            else:
                atom_types[atom] = '[{0}]'.format(','.join([a.label for a in atom.atomtype]))
            # Unpaired Electron(s)
            if len(atom.radical_electrons) == 1:
                atom_unpaired_electrons[atom] = str(atom.radical_electrons[0])
            elif len(atom.radical_electrons) == 0:
                atom_unpaired_electrons[atom] = 'x'  # Empty list indicates wildcard
            else:
                atom_unpaired_electrons[atom] = '[{0}]'.format(','.join([str(radical) for radical in atom.radical_electrons]))

            # Lone Electron Pair(s)
            if len(atom.lone_pairs) == 1:
                atom_lone_pairs[atom] = str(atom.lone_pairs[0])
            elif len(atom.lone_pairs) == 0:
                atom_lone_pairs[atom] = None  # Empty list indicates wildcard
            else:
                atom_lone_pairs[atom] = '[{0}]'.format(','.join([str(pair) for pair in atom.lone_pairs]))

            # Charges
            if len(atom.charge) == 1:
                atom_charge[atom] = '+' + str(atom.charge[0]) if atom.charge[0] > 0 else str(atom.charge[0])
            elif len(atom.charge) == 0:
                atom_charge[atom] = None  # Empty list indicates wildcard
            else:
                atom_charge[atom] = '[{0}]'.format(','.join(['+'+str(charge) if charge > 0 else ''+str(charge) for charge in atom.charge]))

            # Sites
            if len(atom.site) == 1:
                atom_site[atom] = "\"" + atom.site[0] + "\""
            elif len(atom.site) == 0:
                atom_site[atom] = None  # Empty list indicates wildcard
            else:
                atom_site[atom] = '["{0}"]'.format('","'.join(s for s in atom.site))
            
            # Morphologies
            if len(atom.morphology) == 1:
                atom_morphology[atom] = "\"" + atom.morphology[0] + "\""
            elif len(atom.morphology) == 0:
                atom_morphology[atom] = None  # Empty list indicates wildcard
            else:
                atom_morphology[atom] = '["{0}"]'.format('","'.join(s for s in atom.morphology))
                
            # Isotopes
            atom_isotope[atom] = -1

            # Other props
            props = []
            if 'inRing' in atom.props:
                props.append(' r{0}'.format(int(atom.props['inRing'])))
            atom_props[atom] = props
    else:
        for atom in atom_numbers:
            # Atom type
            atom_types[atom] = '{0}'.format(atom.symbol)
            # Unpaired Electron(s)
            atom_unpaired_electrons[atom] = '{0}'.format(atom.radical_electrons)
            # Lone Electron Pair(s)
            atom_lone_pairs[atom] = str(atom.lone_pairs)
            # Partial Charge(s)
            atom_charge[atom] = '+' + str(atom.charge) if atom.charge > 0 else '' + str(atom.charge)
            # Sites
            if atom.site:
                atom_site[atom] = "\"" + atom.site + "\""
            else:
                atom_site[atom] = None
            # Morphology
            if atom.morphology:
                atom_morphology[atom] = "\"" + atom.morphology + "\""
            else:
                atom_morphology[atom] = None
            # Isotopes
            if isinstance(atom, Atom):
                atom_isotope[atom] = atom.element.isotope
            else:
                # cutting labels in fragment cases
                atom_isotope[atom] = atom.isotope

    # Determine field widths
    atom_number_width = max([len(s) for s in atom_numbers.values()]) + 1
    atom_label_width = max([len(s) for s in atom_labels.values()])
    if atom_label_width > 0:
        atom_label_width += 1
    atom_type_width = max([len(s) for s in atom_types.values()]) + 1
    atom_unpaired_electrons_width = max([len(s) for s in atom_unpaired_electrons.values()])

    # Assemble the adjacency list
    for atom in atoms:
        if atom not in atom_numbers:
            continue

        # Atom number
        adjlist += '{0:<{1:d}}'.format(atom_numbers[atom], atom_number_width)
        # Atom label
        adjlist += '{0:<{1:d}}'.format(atom_labels[atom], atom_label_width)
        # Atom type(s)
        adjlist += '{0:<{1:d}}'.format(atom_types[atom], atom_type_width)
        # Unpaired Electron(s)
        adjlist += 'u{0:<{1:d}}'.format(atom_unpaired_electrons[atom], atom_unpaired_electrons_width)
        # Lone Electron Pair(s)
        if atom_lone_pairs[atom] is not None:
            adjlist += ' p{0}'.format(atom_lone_pairs[atom])
        # Partial charges
        if atom_charge[atom] is not None:
            adjlist += ' c{0}'.format(atom_charge[atom])
        # Sites
        if atom_site[atom]:
            adjlist += ' s{0}'.format(atom_site[atom])
        # Morphologies
        if atom_morphology[atom]:
            adjlist += ' m{0}'.format(atom_morphology[atom])
        # Isotopes
        if atom_isotope[atom] != -1:
            adjlist += ' i{0}'.format(atom_isotope[atom])
        if group and len(atom_props[atom]) > 0:
            for prop in atom_props[atom]:
                adjlist += prop

        # Bonds list
        atoms2 = list(atom.bonds.keys())
        # sort them the same way as the atoms
        atoms2.sort(key=atoms.index)

        for atom2 in atoms2:
            if atom2 not in atom_numbers:
                continue

            bond = atom.bonds[atom2]
            adjlist += ' {{{0},'.format(atom_numbers[atom2])

            # Bond type(s)
            if group:
                code = '[{0}]'
                if len(bond.order) == 1:
                    code = '{0}'
                # preference is for string representation, backs down to number
                # numbers if doesn't work
                try:
                    adjlist += code.format(','.join(bond.get_order_str()))
                except ValueError:
                    adjlist += code.format(','.join(str(bond.get_order_num())))
            else:
                # preference is for string representation, backs down to number
                # numbers if doesn't work
                try:
                    adjlist += bond.get_order_str()
                except ValueError:
                    adjlist += str(bond.get_order_num())
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist


def get_old_electron_state(atom):
    """
    Get the old adjacency list format electronic state
    """
    additional_lone_pairs = atom.lone_pairs - PeriodicSystem.lone_pairs[atom.element.symbol]
    electrons = atom.radical_electrons + additional_lone_pairs * 2
    if electrons == 0:
        electron_state = '0'
    elif electrons == 1:
        electron_state = '1'
    elif electrons == 2:
        if additional_lone_pairs == 0:
            electron_state = '2T'
        elif additional_lone_pairs == 1:
            electron_state = '2S'
        else:
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    elif electrons == 3:
        if additional_lone_pairs == 0:
            electron_state = '3Q'
        elif additional_lone_pairs == 1:
            electron_state = '3D'
        else:
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    elif electrons == 4:
        if additional_lone_pairs == 0:
            electron_state = '4V'
        elif additional_lone_pairs == 1:
            electron_state = '4T'
        elif additional_lone_pairs == 2:
            electron_state = '4S'
        else:
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    else:
        raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    return electron_state
