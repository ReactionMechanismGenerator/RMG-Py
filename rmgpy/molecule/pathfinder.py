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
This module provides functions for searching paths within a molecule.
The paths generally consist of alternating atoms and bonds.
"""
import cython
import itertools

from Queue import Queue


def find_butadiene(start, end):
    """
    Search for a path between start and end atom that consists of 
    alternating non-single and single bonds.

    Returns a list with atom and bond elements from start to end, or
    None if nothing was found.
    """
    
    q = Queue()#FIFO queue of paths that need to be analyzed
    q.put([start])

    while not q.empty():
        path = q.get()
        # search for end atom among the neighbors of the terminal atom of the path:
        terminal = path[-1]
        assert isinstance(terminal, Atom)
        for atom4, bond34 in terminal.bonds.iteritems():
            if atom4 == end and not bond34.isSingle():# we have found the path we are looking for
                #add the final bond and atom and return
                path.append(bond34)
                path.append(atom4)
                return path
        else:#none of the neighbors is the end atom.
            # Add a new allyl path and try again:
            new_paths = add_allyls(path)
            [q.put(p) if p else '' for p in new_paths]

    # Could not find a resonance path from start atom to end atom
    return None  

def find_butadiene_end_with_charge(start):
    """
    Search for a (4-atom, 3-bond) path between start and end atom that consists of 
    alternating non-single and single bonds and ends with a charged atom.

    Returns a list with atom and bond elements from start to end, or
    None if nothing was found.
    """

    q = Queue()#FIFO queue of paths that need to be analyzed
    q.put([start])

    while not q.empty():
        path = q.get()
        # search for end atom among the neighbors of the terminal atom of the path:
        terminal = path[-1]
        assert isinstance(terminal, Atom)
        for atom4, bond34 in terminal.bonds.iteritems():
            if atom4.charge != 0 and not bond34.isSingle() and not atom4 in path:# we have found the path we are looking for
                #add the final bond and atom and return
                path.append(bond34)
                path.append(atom4)
                return path
        else:#none of the neighbors is the end atom.
            # Add a new allyl path and try again:
            new_paths = add_allyls(path)
            [q.put(p) if p else '' for p in new_paths]

    # Could not find a resonance path from start atom to end atom
    return None

def find_allyl_end_with_charge(start):
    """
    Search for a (3-atom, 2-bond) path between start and end atom that consists of 
    alternating non-single and single bonds and ends with a charged atom.

    Returns a list with atom and bond elements from start to end, or
    an empty list if nothing was found.
    """
    paths = []

    q = Queue()#FIFO queue of paths that need to be analyzed
    unsaturated_bonds = add_unsaturated_bonds([start])
    
    if not unsaturated_bonds:
        return []
    
    [q.put(path) for path in unsaturated_bonds]

    while not q.empty():
        path = q.get()
        # search for end atom among the neighbors of the terminal atom of the path:
        terminal = path[-1]
        assert isinstance(terminal, Atom)

        path_copy = path[:]
        for atom3, bond23 in terminal.bonds.iteritems():
            if atom3.charge != 0 and not atom3 in path_copy:# we have found the path we are looking for
                #add the final bond and atom and return
                path_copy_copy = path_copy[:]
                path_copy_copy.extend([bond23, atom3])
                paths.append(path_copy_copy)
        else:#none of the neighbors is the end atom.
            # Add a new inverse allyl path and try again:
            new_paths = add_inverse_allyls(path)
            [q.put(p) if p else '' for p in new_paths]

    # Could not find a resonance path from start atom to end atom
    return paths

def find_shortest_path(start, end, path=None):
    path = path if path else []
    path = path + [start]
    if start == end:
        return path

    shortest = None
    for node,_ in start.edges.iteritems():
        if node not in path:
            newpath = find_shortest_path(node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest

def add_unsaturated_bonds(path):
    """
    Find all the (2-atom, 1-bond) patterns "X=X" starting from the 
    last atom of the existing path.

    The bond attached to the starting atom should be non single.
    """
    paths = []
    start = path[-1]
    assert isinstance(start, Atom)

    for atom2, bond12 in start.bonds.iteritems():
        if not bond12.isSingle() and not atom2 in path and atom2.number!= 1:
            new_path = path[:]
            new_path.extend((bond12, atom2))
            paths.append(new_path)
    return paths 

def add_allyls(path):
    """
    Find all the (3-atom, 2-bond) patterns "X=X-X" starting from the 
    last atom of the existing path.

    The bond attached to the starting atom should be non single.
    The second bond should be single.
    """
    paths = []
    start = path[-1]
    assert isinstance(start, Atom)

    for atom2, bond12 in start.bonds.iteritems():
        if not bond12.isSingle() and not atom2 in path:
            for atom3, bond23 in atom2.bonds.iteritems():
                if start is not atom3 and atom3.number!= 1:
                    new_path = path[:]
                    new_path.extend((bond12, atom2, bond23, atom3))
                    paths.append(new_path)
    return paths

def add_inverse_allyls(path):
    """
    Find all the (3-atom, 2-bond) patterns "start~atom2=atom3" starting from the 
    last atom of the existing path.

    The second bond should be non-single.
    """
    paths = []
    start = path[-1]
    assert isinstance(start, Atom)

    for atom2, bond12 in start.bonds.iteritems():
        if not atom2 in path:
            for atom3, bond23 in atom2.bonds.iteritems():
                if not atom3 in path and atom3.number!= 1 and not bond23.isSingle():
                    new_path = path[:]
                    new_path.extend((bond12, atom2, bond23, atom3))
                    paths.append(new_path)
    return paths

def compute_atom_distance(atom_indices, mol):
    """
    Compute the distances between each pair of atoms in the atom_indices.

    The distance between two atoms is defined as the length of the shortest path
    between the two atoms minus 1, because the start atom is part of the path.

    The distance between multiple atoms is defined by generating all possible
    combinations between two atoms and storing the distance between each combination
    of atoms in a dictionary.

    The parameter 'atom_indices' is a  list of 1-based atom indices.

    """
    if len(atom_indices) == 1: return {(atom_indices[0],): 0}

    distances = {}
    combos = [sorted(tup) for tup in itertools.combinations(atom_indices, 2)]
    
    for i1, i2 in combos:
        start, end = mol.atoms[i1 - 1], mol.atoms[i2 - 1]
        path = find_shortest_path(start, end)
        distances[(i1, i2)] = len(path) - 1  

    return distances  


def find_allyl_delocalization_paths(atom1):
    """
    Find all the delocalization paths allyl to the radical center indicated
    by `atom1`. Used to generate resonance isomers.
    """
    cython.declare(paths=list, atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)
    
    # No paths if atom1 is not a radical
    if atom1.radicalElectrons <= 0:
        return []

    # Find all delocalization paths
    paths = []
    for atom2, bond12 in atom1.edges.items():
        # Vinyl bond must be capable of gaining an order
        if (bond12.isSingle() or bond12.isDouble()) and atom1.radicalElectrons:
            for atom3, bond23 in atom2.edges.items():
                # Allyl bond must be capable of losing an order without breaking
                if atom1 is not atom3 and (bond23.isDouble() or bond23.isTriple()):
                    paths.append([atom1, atom2, atom3, bond12, bond23])
    return paths


def find_lone_pair_radical_delocalization_paths(atom1):
    """
    Find all the delocalization paths of lone electron pairs next to the radical center indicated
    by `atom1`. Used to generate resonance isomers in adjacent N/O/S atoms.
    Two adjacent O atoms are not allowed since (a) currently RMG has no good thermo/kinetics for R[:O+.][:::O-] which
    could have been generated as a resonance structure of R[::O][::O.].
    The radical site (atom1) could be either:
    - `N u1 p0`, eg O=[N.+][:::O-]
    - `N u1 p1`, eg R[:NH][:NH.]
    - `O u1 p1`, eg [:O.+]=[::N-]; not allowed when adjacent to another O atom
    - `O u1 p2`, eg O=N[::O.]; not allowed when adjacent to another O atom
    - `S u1 p0`, eg O[S.+]([O-])=O
    - `S u1 p1`, eg O[:S.+][O-]
    - `S u1 p2`, eg O=N[::S.]
    - any of the above with more than 1 radical where possible
    The non-radical site (atom2) could respectively be:
    - `N u0 p1`
    - `N u0 p2`
    - `O u0 p2`
    - `O u0 p3`
    - `S u0 p1`
    - `S u0 p2`
    - `S u0 p3`
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    The bond between the sites does not have to be single, e.g.: [:O.+]=[::N-] <=> [::O]=[:N.]
    """
    cython.declare(paths=list, atom2=Atom, bond12=Bond)

    paths = []
    if atom1.isNOS():
        if ((atom1.isNitrogen() and atom1.radicalElectrons >= 1 and atom1.lonePairs in [0, 1])
                or (atom1.isOxygen() and atom1.radicalElectrons >= 1 and atom1.lonePairs in [1, 2])
                or (atom1.isSulfur() and atom1.radicalElectrons >= 1 and atom1.lonePairs in [0, 1, 2])):
            for atom2 in atom1.edges.keys():
                if ((atom2.isNitrogen() and atom2.lonePairs in [1, 2])
                        or (atom2.isOxygen() and atom2.lonePairs in [2, 3] and not atom1.isOxygen())  #avoid RO[::O.] <-> R[:O.+][:::O-], see RMG-Py #1223
                        or (atom2.isSulfur() and atom2.lonePairs in [1, 2, 3])):
                    paths.append([atom1, atom2])
    return paths


def find_lone_pair_multiple_bond_delocalization_paths(atom1):
    """
    Find all the delocalization paths of a N/O/S atom1 which either:
    - Has a lonePair and is bonded by a single/double bond (e.g., [::NH-]-[CH2+], [::N-]=[CH+]) -- direction 1
    - Can obtain a lonePair and is bonded by a double/triple bond (e.g., [:NH]=[CH2], [:N]#[CH]) -- direction 2
    Giving the following resonance transitions, for example:
    - [::NH-]-[CH2+] <=> [:NH]=[CH2]
    - [:N]#[CH] <=> [::N-]=[CH+]
    - N#[N+][O-] <=> <=> [N-]=[N+]=O
    - C[N+](=O)[NH-] <=> <=> C[N+]([O-])=[NH]
    - other examples: S#N, N#[S], O=S([O])=O, [NH]=[N+]=[N-]
    Direction "1" is the direction <increasing> the bond order as in [::NH-]-[CH2+] <=> [:NH]=[CH2]
    Direction "2" is the direction <decreasing> the bond order as in [:NH]=[CH2] <=> [::NH-]-[CH2+]
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    (In direction 1 atom1 <losses> a lone pair, in direction 2 atom1 <gains> a lone pair)
    """
    cython.declare(paths=list, atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)

    paths = []
    if atom1.isNOS():
        for atom2, bond12 in atom1.edges.items():
            if atom2.isNonHydrogen():  # don't bother with hydrogen atoms.
                # Find paths in the direction <increasing> the bond order,
                # atom1 must posses at least one lone pair to loose it
                if ((bond12.isSingle() or bond12.isDouble())
                        and is_NOS_able_to_lose_lone_pair(atom1)):
                    paths.append([atom1, atom2, bond12, 1])  # direction = 1
                # Find paths in the direction <decreasing> the bond order,
                # atom1 gains a lone pair, hence cannot already have more than two lone pairs
                if ((bond12.isDouble() or bond12.isTriple())
                        and is_NOS_able_to_gain_lone_pair(atom1)):
                    paths.append([atom1, atom2, bond12, 2])  # direction = 2
    return paths


def find_lone_pair_radical_multiple_bond_delocalization_paths(atom1):
    """
    Find all the delocalization paths of a N/O/S atom1 which either:
    - Has a lonePair and is bonded by a single/double bond to a radical atom (e.g., [::N]-[.CH2])
    - Can obtain a lonePair, has a radical, and is bonded by a double/triple bond (e.g., [:N.]=[CH2])
    Giving the following resonance transitions, for example:
    - [::N]-[.CH2] <=> [:N.]=[CH2]
    - O[:S](=O)[::O.] <=> O[S.](=O)=[::O]
    Direction "1" is the direction <increasing> the bond order as in [::N]-[.CH2] <=> [:N.]=[CH2]
    Direction "2" is the direction <decreasing> the bond order as in [:N.]=[CH2] <=> [::N]-[.CH2]
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    (In direction 1 atom1 <losses> a lone pair, gains a radical, and atom2 looses a radical.
    In direction 2 atom1 <gains> a lone pair, looses a radical, and atom2 gains a radical)
    """
    cython.declare(paths=list, atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)

    paths = []
    if atom1.isNOS():
        for atom2, bond12 in atom1.edges.items():
            # Find paths in the direction <increasing> the bond order
            # atom1 must posses at least one lone pair to loose it, atom2 must be a radical
            if (atom2.radicalElectrons and (bond12.isSingle() or bond12.isDouble())
                    and is_NOS_able_to_lose_lone_pair(atom1)):
                paths.append([atom1, atom2, bond12, 1])  # direction = 1
            # Find paths in the direction <decreasing> the bond order
            # atom1 gains a lone pair, hence cannot already have more than two lone pairs, and is also a radical
            if (atom1.radicalElectrons and (bond12.isDouble() or bond12.isTriple())
                    and is_NOS_able_to_gain_lone_pair(atom1)):
                paths.append([atom1, atom2, bond12, 2])  # direction = 2
    return paths


def find_N5ddc_N5tc_delocalization_paths(atom1):
    """
    Find all the resonance structures of nitrogen atoms with two double bonds (atomType N5ddc)
    and nitrogen atoms with one triple and one single bond (atomType N5tc).
    Examples:
    - N2O (N#[N+][O-] <-> [N-]=[N+]=O)
    - Azide (N#[N+][NH-] <-> [N-]=[N+]=N <-> [N-2][N+]#[NH+])
    - N#N group on sulfur (O[S-](O)[N+]#N <-> OS(O)=[N+]=[N-] <-> O[S+](O)#[N+][N-2])
    In this transition atom1 is the middle N+ (N5ddc or N5tc)
    A "if atom.atomType.label in ['N5ddc','N5tc']" check should be done before calling this function
    """
    cython.declare(paths=list, atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)

    paths = []

    for atom2, bond12 in atom1.edges.items():
        if atom2.isNOS():
            for atom3, bond23 in atom1.edges.items():
                if atom2 is not atom3:
                    # Find this transitions from the N5tc side.
                    # atom2 is single-bonded to atom1 and looses a lone pair
                    # atom3 is triple-bonded to atom1 and gains a lone pair (atom3 cannot be oxygen)
                    # bond12 is incremented, bond23 is decremented
                    if (((atom2.isNitrogen() and atom2.lonePairs in [2, 3])
                            or (atom2.isOxygen() and atom2.lonePairs == 3)
                            or (atom2.isSulfur() and atom2.lonePairs in [1, 2, 3]))
                            and ((atom3.isNitrogen() or atom3.isSulfur()) and atom3.lonePairs in [0, 1])
                            and bond12.isSingle() and bond23.isTriple()):
                        paths.append([atom1, atom2, atom3, bond12, bond23])
                    # Find this transitions from the N5ddc side.
                    # atom2 looses a lone pair and ends up with a triple bond to atom1 (atom2 cannot be oxygen)
                    # atom3 gains a lone pair and ends up with a single bond to atom1
                    # bond12 is incremented, bond23 is decremented (same actions as above under these definitions)
                    elif (((atom2.isNitrogen() or atom2.isSulfur()) and atom2.lonePairs in [1, 2])
                            and ((atom3.isNitrogen() and atom3.lonePairs in [1, 2])
                            or (atom3.isOxygen() and atom3.lonePairs == 2)
                            or (atom3.isSulfur() and atom3.lonePairs in [0, 1, 2]))
                            and bond12.isDouble() and bond23.isDouble()):
                        paths.append([atom1, atom2, atom3, bond12, bond23])
    return paths


def find_N5dc_delocalization_paths(atom1):
    """
    Find all the resonance structures of an N5dc nitrogen atom with a single bond to a radical N/O/S site, another
    single bond to a negatively charged N/O/S site, and one double bond (not participating in this transformation)
    Example:
    - N=[N+]([O])([O-]) <=> N=[N+]([O-])([O]), these structures are isomorphic but not identical, the transition is
    important for correct degeneracy calculations
    In this transition atom1 is the middle N+ (N5dc), atom2 is the radical site, and atom3 is negatively charged
    A "if atom1.atomType.label == 'N5dc'" check should be done before calling this function
    """
    cython.declare(paths=list, atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)

    paths = []

    for atom2, bond12 in atom1.edges.items():
        if atom2.radicalElectrons and bond12.isSingle() and not atom2.charge and is_NOS_able_to_gain_lone_pair(atom2):
            for atom3, bond23 in atom1.edges.items():
                if (atom2 is not atom3 and bond23.isSingle() and atom3.charge < 0
                        and is_NOS_able_to_lose_lone_pair(atom3)):
                    paths.append([atom2, atom3])
                    return paths  # there could only be one such path per atom1, return if found
    return paths


def is_NOS_able_to_gain_lone_pair(atom):
    """
    Helper function
    Returns True if atom is N/O/S and is able to <gain> an additional lone pair, False otherwise
    We don't allow O to remain with no lone pairs
    """
    return (((atom.isNitrogen() or atom.isSulfur()) and atom.lonePairs in [0, 1, 2])
                        or (atom.isOxygen() and atom.lonePairs in [1, 2]))


def is_NOS_able_to_lose_lone_pair(atom):
    """
    Helper function
    Returns True if atom is N/O/S and is able to <loose> a lone pair, False otherwise
    We don't allow O to remain with no lone pairs
    """
    return (((atom.isNitrogen() or atom.isSulfur()) and atom.lonePairs in [1, 2, 3])
                        or (atom.isOxygen() and atom.lonePairs in [2, 3]))
