

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

def find_shortest_path(start, end, path=[]):
    path = path + [start]
    if start == end:
        return path

    shortest = None
    for node,_ in start.bonds.iteritems():
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


def findAllDelocalizationPaths(atom1):
    """
    Find all the delocalization paths allyl to the radical center indicated
    by `atom1`. Used to generate resonance isomers.
    """
    cython.declare(paths=list)
    cython.declare(atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)
    
    # No paths if atom1 is not a radical
    if atom1.radicalElectrons <= 0:
        return []

    # Find all delocalization paths
    paths = []
    for atom2, bond12 in atom1.edges.items():
        # Vinyl bond must be capable of gaining an order
        if (bond12.isSingle() or bond12.isDouble()) and atom1.radicalElectrons == 1:
            for atom3, bond23 in atom2.edges.items():
                # Allyl bond must be capable of losing an order without breaking
                if atom1 is not atom3 and (bond23.isDouble() or bond23.isTriple()):
                    paths.append([atom1, atom2, atom3, bond12, bond23])
    return paths

def findAllDelocalizationPathsLonePairRadical(atom1):
    """
    Find all the delocalization paths of lone electron pairs next to the radical center indicated
    by `atom1`. Used to generate resonance isomers.
    """
    cython.declare(paths=list)
    cython.declare(atom2=Atom, bond12=Bond)
    
    # No paths if atom1 is not a radical
    if atom1.radicalElectrons <= 0:
        return []
    
    # In a first step we only consider nitrogen and oxygen atoms as possible radical centers
    if not ((atom1.lonePairs == 0 and atom1.isNitrogen()) or(atom1.lonePairs == 2 and atom1.isOxygen())):
        return []
    
    # Find all delocalization paths
    paths = []
    for atom2, bond12 in atom1.edges.items():
        # Only single bonds are considered
        if bond12.isSingle():
            # Neighboring atom must posses a lone electron pair to loose it
            if ((atom2.lonePairs == 1 and atom2.isNitrogen()) or (atom2.lonePairs == 3 and atom2.isOxygen())) and (atom2.radicalElectrons == 0):
                paths.append([atom1, atom2])
                
    return paths

def findAllDelocalizationPathsN5dd_N5ts(atom1):
    """
    Find all the resonance structures of nitrogen atoms with two double bonds (N5dd)
    and nitrogen atoms with one triple and one single bond (N5ts)
    """
    cython.declare(paths=list)
    cython.declare(atom2=Atom, bond12=Bond)
    
    # No paths if atom1 is not nitrogen
    if not (atom1.isNitrogen()):
        return []
    
    # Find all delocalization paths
    paths = []
    index_atom_2 = 0
    index_atom_3 = 0
    
    for atom2, bond12 in atom1.edges.items():
        index_atom_2 = index_atom_2 + 1
        # Only double bonds are considered
        if bond12.isDouble():
            for atom3, bond13 in atom1.edges.items():
                index_atom_3 = index_atom_3 + 1
                # Only double bonds are considered, at the moment we only consider non-radical nitrogen and oxygen atoms
                if (bond13.isDouble() and atom3.radicalElectrons == 0 and atom3.lonePairs > 0 and not atom3.isOxygen() and not atom3.isCarbon() and (index_atom_2 != index_atom_3)):
                    paths.append([atom1, atom2, atom3, bond12, bond13, 1])
    
    for atom2, bond12 in atom1.edges.items():
        # Only triple bonds are considered
        if bond12.isTriple():
            for atom3, bond13 in atom1.edges.items():
                # Only single bonds are considered, at the moment we only consider negatively charged nitrogen and oxygen
                if (bond13.isSingle() and ((atom3.isNitrogen() and atom3.lonePairs >= 2) or (atom3.isOxygen() and atom3.lonePairs >= 3))):
                    paths.append([atom1, atom2, atom3, bond12, bond13, 2])
    
    return paths    