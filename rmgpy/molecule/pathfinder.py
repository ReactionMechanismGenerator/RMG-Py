

"""
This module provides functions for searching paths within a molecule.
The paths generally consist of alternating atoms and bonds.
"""

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