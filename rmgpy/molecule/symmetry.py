#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This module provides functionality for estimating the symmetry number of a
molecule from its chemical graph representation.
"""
import itertools


def calculate_atom_symmetry_number(molecule, atom):
    """
    Return the symmetry number centered at `atom` in the structure. The
    `atom` of interest must not be in a cycle.
    """
    symmetry_number = 1

    single = double = triple = benzene = num_neighbors = 0  # note that 0 is immutable
    for bond in atom.edges.values():
        if bond.is_single():
            single += 1
        elif bond.is_double():
            double += 1
        elif bond.is_triple():
            triple += 1
        elif bond.is_benzene():
            benzene += 1
        num_neighbors += 1

    # If atom has zero or one neighbors, the symmetry number is 1
    if num_neighbors < 2:
        return symmetry_number

    # Create temporary structures for each functional group attached to atom
    molecule0 = molecule
    molecule = molecule0.copy(True)
    atom = molecule.vertices[molecule0.vertices.index(atom)]
    molecule.remove_atom(atom)
    groups = molecule.split()

    # Determine equivalence of functional groups around atom
    group_isomorphism = dict([(group, dict()) for group in groups])
    for group1 in groups:
        for group2 in groups:
            if group1 is not group2 and group2 not in group_isomorphism[group1]:
                group_isomorphism[group1][group2] = group1.is_isomorphic(group2)
                group_isomorphism[group2][group1] = group_isomorphism[group1][group2]
            elif group1 is group2:
                group_isomorphism[group1][group1] = True
    count = [sum([int(group_isomorphism[group1][group2]) for group2 in groups]) for group1 in groups]
    for i in range(count.count(2) // 2):
        count.remove(2)
    for i in range(count.count(3) // 3):
        count.remove(3)
        count.remove(3)
    for i in range(count.count(4) // 4):
        count.remove(4)
        count.remove(4)
        count.remove(4)
    count.sort()
    count.reverse()

    if atom.radical_electrons == 0:
        if single == 4:
            # Four single bonds
            if count == [4]:
                symmetry_number *= 12
            elif count == [3, 1]:
                symmetry_number *= 3
            elif count == [2, 2]:
                symmetry_number *= 2
            elif count == [2, 1, 1]:
                symmetry_number *= 1
            elif count == [1, 1, 1, 1]:
                symmetry_number *= 0.5  # found chirality
        elif single == 3:
            # Three single bonds
            if count == [3]:
                symmetry_number *= 3
            elif count == [2, 1]:
                symmetry_number *= 1
            elif count == [1, 1, 1]:
                symmetry_number *= 1
        elif single == 2:
            # Two single bonds
            if count == [2]:
                symmetry_number *= 2
        # for resonance hybrids
        elif single == 1:
            if count == [2, 1]:
                symmetry_number *= 2
        elif double == 2:
            # Two double bonds
            if count == [2]:
                symmetry_number *= 2
        # for nitrogen resonance hybrids
        elif single == 0:
            if count == [2]:
                symmetry_number *= 2
    elif atom.radical_electrons == 1:
        if single == 3:
            # Three single bonds
            if count == [3]:
                symmetry_number *= 6
            elif count == [2, 1]:
                symmetry_number *= 2
            elif count == [1, 1, 1]:
                symmetry_number *= 1
        elif single == 1:
            if count == [2, 1]:
                symmetry_number *= 2
            elif count == [1, 1, 1]:
                symmetry_number *= 1
    elif atom.radical_electrons == 2:
        if single == 2:
            # Two single bonds
            if count == [2]:
                symmetry_number *= 2

    return symmetry_number


################################################################################

def calculate_bond_symmetry_number(molecule, atom1, atom2):
    """
    Return the symmetry number centered at `bond` in the structure.
    """
    bond = atom1.edges[atom2]
    symmetry_number = 1
    if atom1.equivalent(atom2):
        # An O-O bond is considered to be an "optical isomer" and so no
        # symmetry correction will be applied
        if (atom1.atomtype.label == 'O2s' and atom2.atomtype.label == 'O2s' and
                atom1.radical_electrons == atom2.radical_electrons == 0):
            return symmetry_number
        # If the molecule is diatomic, then we don't have to check the
        # ligands on the two atoms in this bond (since we know there
        # aren't any)
        elif len(molecule.vertices) == 2:
            symmetry_number = 2
        else:
            molecule.remove_bond(bond)
            structure = molecule.copy(True)
            molecule.add_bond(bond)

            atom1 = structure.atoms[molecule.atoms.index(atom1)]
            atom2 = structure.atoms[molecule.atoms.index(atom2)]
            fragments = structure.split()

            if len(fragments) != 2:
                return symmetry_number

            fragment1, fragment2 = fragments
            if atom1 in fragment1.atoms:
                fragment1.remove_atom(atom1)
            if atom2 in fragment1.atoms:
                fragment1.remove_atom(atom2)
            if atom1 in fragment2.atoms:
                fragment2.remove_atom(atom1)
            if atom2 in fragment2.atoms:
                fragment2.remove_atom(atom2)
            groups1 = fragment1.split()
            groups2 = fragment2.split()

            # Test functional groups for symmetry
            if len(groups1) == len(groups2) == 1:
                if groups1[0].is_isomorphic(groups2[0]):
                    symmetry_number *= 2
            elif len(groups1) == len(groups2) == 2:
                if groups1[0].is_isomorphic(groups2[0]) and groups1[1].is_isomorphic(groups2[1]):
                    symmetry_number *= 2
                elif groups1[1].is_isomorphic(groups2[0]) and groups1[0].is_isomorphic(groups2[1]):
                    symmetry_number *= 2
            elif len(groups1) == len(groups2) == 3:
                if groups1[0].is_isomorphic(groups2[0]):
                    if groups1[1].is_isomorphic(groups2[1]) and groups1[2].is_isomorphic(groups2[2]):
                        symmetry_number *= 2
                    elif groups1[1].is_isomorphic(groups2[2]) and groups1[2].is_isomorphic(groups2[1]):
                        symmetry_number *= 2
                elif groups1[0].is_isomorphic(groups2[1]):
                    if groups1[1].is_isomorphic(groups2[2]) and groups1[2].is_isomorphic(groups2[0]):
                        symmetry_number *= 2
                    elif groups1[1].is_isomorphic(groups2[0]) and groups1[2].is_isomorphic(groups2[2]):
                        symmetry_number *= 2
                elif groups1[0].is_isomorphic(groups2[2]):
                    if groups1[1].is_isomorphic(groups2[0]) and groups1[2].is_isomorphic(groups2[1]):
                        symmetry_number *= 2
                    elif groups1[1].is_isomorphic(groups2[1]) and groups1[2].is_isomorphic(groups2[0]):
                        symmetry_number *= 2

    return symmetry_number


################################################################################

def calculate_axis_symmetry_number(molecule):
    r"""
    Get the axis symmetry number correction. The "axis" refers to a series
    of two or more cumulated double bonds (e.g. C=C=C, etc.). Corrections
    for single C=C bonds are handled in getBondSymmetryNumber().
    
    Each axis (C=C=C) has the potential to double the symmetry number.
    If an end has 0 or 1 groups (eg. =C=CJJ or =C=C-R) then it cannot 
    alter the axis symmetry and is disregarded::
    
        A=C=C=C..        A-C=C=C=C-A
        
          s=1                s=1
    
    If an end has 2 groups that are different then it breaks the symmetry 
    and the symmetry for that axis is 1, no matter what's at the other end::
    
        A\               A\         /A
          T=C=C=C=C-A      T=C=C=C=T
        B/               A/         \B
              s=1             s=1
    
    If you have one or more ends with 2 groups, and neither end breaks the 
    symmetry, then you have an axis symmetry number of 2::
    
        A\         /B      A\         
          C=C=C=C=C          C=C=C=C-B
        A/         \B      A/         
              s=2                s=2
    """

    symmetry_number = 1

    # List all double bonds in the structure
    double_bonds = []
    for atom1 in molecule.vertices:
        for atom2 in atom1.edges:
            if (atom1.edges[atom2].is_double() or atom1.edges[atom2].order > 2) \
                    and molecule.vertices.index(atom1) < molecule.vertices.index(atom2):
                double_bonds.append((atom1, atom2))

    # Search for adjacent double bonds
    cumulated_bonds = []
    for i, bond1 in enumerate(double_bonds):
        atom11, atom12 = bond1
        for bond2 in double_bonds[i + 1:]:
            atom21, atom22 = bond2
            if atom11 is atom21 or atom11 is atom22 or atom12 is atom21 or atom12 is atom22:
                list_to_add_to = None
                for cumBonds in cumulated_bonds:
                    if (atom11, atom12) in cumBonds or (atom21, atom22) in cumBonds:
                        list_to_add_to = cumBonds
                if list_to_add_to is not None:
                    if (atom11, atom12) not in list_to_add_to:
                        list_to_add_to.append((atom11, atom12))
                    if (atom21, atom22) not in list_to_add_to:
                        list_to_add_to.append((atom21, atom22))
                else:
                    cumulated_bonds.append([(atom11, atom12), (atom21, atom22)])

    # Also keep isolated double bonds
    for bond1 in double_bonds:
        for bonds in cumulated_bonds:
            if bond1 in bonds:
                break
        else:
            cumulated_bonds.append([bond1])

    # For each set of adjacent double bonds, check for axis symmetry
    for bonds in cumulated_bonds:
        # Do nothing if axis is in cycle
        found = False
        for atom1, atom2 in bonds:
            if molecule.is_bond_in_cycle(atom1.edges[atom2]):
                found = True
        if found:
            continue

        # Find terminal atoms in axis
        # Terminal atoms labelled T:  T=C=C=C=T
        axis = []
        for bond in bonds:
            axis.extend(bond)
        terminal_atoms = []
        for atom in axis:
            if axis.count(atom) == 1:
                terminal_atoms.append(atom)
        if len(terminal_atoms) != 2:
            continue

        # Remove axis from (copy of) structure
        bond_list = []
        for atom1, atom2 in bonds:
            bond = atom1.edges[atom2]
            bond_list.append(bond)
            molecule.remove_bond(bond)
        structure = molecule.copy(True)
        terminal_atoms = [structure.vertices[molecule.vertices.index(atom)] for atom in terminal_atoms]
        for bond in bond_list:
            molecule.add_bond(bond)

        atoms_to_remove = []
        for atom in structure.vertices:
            if len(atom.edges) == 0 and atom not in terminal_atoms:  # it's not bonded to anything
                atoms_to_remove.append(atom)
        for atom in atoms_to_remove:
            structure.remove_atom(atom)

        # Split remaining fragments of structure
        end_fragments = structure.split()

        # 
        # there can be two groups at each end     A\         /B
        #                                           T=C=C=C=T
        #                                         A/         \B

        # to start with nothing has broken symmetry about the axis
        symmetry_broken = False
        end_fragments_to_remove = []
        for fragment in end_fragments:  # a fragment is one end of the axis
            # remove the atom that was at the end of the axis and split what's left into groups
            terminal_atom = None
            for atom in terminal_atoms:
                if atom in fragment.atoms:
                    terminal_atom = atom
                    fragment.remove_atom(atom)
                    break
            else:
                continue

            groups = []
            if len(fragment.atoms) > 0:
                groups = fragment.split()

            # If end has only one group then it can't contribute to (nor break) axial symmetry
            #   Eg. this has no axis symmetry:   A-T=C=C=C=T-A
            # so we remove this end from the list of interesting end fragments
            if len(groups) == 0:
                end_fragments_to_remove.append(fragment)
                continue  # next end fragment
            elif len(groups) == 1 and terminal_atom.radical_electrons == 0:
                if terminal_atom.atomtype.label == 'N3d':
                    symmetry_broken = True
                else:
                    end_fragments_to_remove.append(fragment)
                    continue  # next end fragment
            elif len(groups) == 1 and terminal_atom.radical_electrons != 0:
                symmetry_broken = True
            elif len(groups) == 2:
                if not groups[0].is_isomorphic(groups[1]):
                    # this end has broken the symmetry of the axis
                    symmetry_broken = True

        for fragment in end_fragments_to_remove:
            end_fragments.remove(fragment)

        # If there are end fragments left that can contribute to symmetry,
        # and none of them broke it, then double the symmetry number
        # NB>> This assumes coordination number of 4 (eg. Carbon).
        #      And would be wrong if we had /B
        #                         =C=C=C=C=T-B
        #                                   \B
        #      (for some T with coordination number 5).
        if end_fragments and not symmetry_broken:
            symmetry_number *= 2

    return symmetry_number


################################################################################

def calculate_cyclic_symmetry_number(molecule):
    """
    Get the symmetry number correction for cyclic regions of a molecule.
    For complicated fused rings the smallest set of smallest rings is used.
    """
    symmetry_number = 1

    # for polycyclics, We should be getting the largest ring, not the smallest
    single_rings, polycyclic_rings = molecule.get_disparate_cycles()
    for ring in polycyclic_rings:
        single_rings.append(molecule.get_largest_ring(ring[0]))
    # Get symmetry number for each ring in structure & multiply
    for ring in single_rings:
        ring = molecule.sort_cyclic_vertices(ring)
        size = len(ring)

        # look for twisting rotation
        # go through each possible number of symmetrical sets
        for num_sections in range(size, 0, -1):
            # only go through if it can give symmetry (only factors of size)
            if size % num_sections == 0:
                # only check the minimum number of sections necessary
                num_rotations = size // num_sections
                # check rotation around the ring
                all_the_same = True
                starting_index = 0
                while all_the_same and starting_index < size // 2:
                    for atom_index in range(num_rotations, size, num_rotations):
                        if not _indistinguishable(ring[starting_index], ring[(starting_index + atom_index) % size]):
                            all_the_same = False
                            break
                    starting_index += 1
                if all_the_same:
                    symmetry_number *= num_sections
                    break

        # look for flipping rotation. 
        if size % 2 == 0:  # even length only has to go through half
            flipping_atom_indexes = list(range(size // 2))
        else:
            flipping_atom_indexes = list(range(size))

        for flipping_atom_index in flipping_atom_indexes:

            # check for flipping with across an axis containing atoms
            all_the_same = True
            min_index = flipping_atom_index + 1
            max_index = flipping_atom_index + size - 1
            while min_index <= max_index:
                # ensure the two atoms are different. use mod size to loop to the start
                # of the list when index out of bounds
                if not _indistinguishable(ring[min_index % size], ring[max_index % size]):
                    all_the_same = False
                    break
                min_index += 1
                max_index += -1
            # check to make sure that the groups are identical on centers of flipping
            if all_the_same:
                ringed_atom_ids = [id(_atom) for _atom in ring]
                if size % 2 == 0:  # look at two atoms
                    atom1 = ring[flipping_atom_index]
                    atom2 = ring[flipping_atom_index + size // 2]
                    non_ring_bonded_atoms = [bonded_atom for bonded_atom in itertools.chain(atom1.bonds.keys(), atom2.bonds.keys())
                                             if id(bonded_atom) not in ringed_atom_ids]
                    if len(non_ring_bonded_atoms) < 3:
                        pass  # all_the_same still true
                    elif len(non_ring_bonded_atoms) == 3:
                        # at least one of these much be a match for flipping to happen
                        identical = _indistinguishable(non_ring_bonded_atoms[0], non_ring_bonded_atoms[1])
                        identical2 = _indistinguishable(non_ring_bonded_atoms[0], non_ring_bonded_atoms[2])
                        identical3 = _indistinguishable(non_ring_bonded_atoms[1], non_ring_bonded_atoms[2])
                        if not (identical or identical2 or identical3):
                            all_the_same = False
                    elif len(non_ring_bonded_atoms) == 4:
                        same_sides = _indistinguishable(non_ring_bonded_atoms[0], non_ring_bonded_atoms[1]) and \
                                     _indistinguishable(non_ring_bonded_atoms[2], non_ring_bonded_atoms[3])
                        if not same_sides:
                            atom0_matching = _indistinguishable(non_ring_bonded_atoms[0], non_ring_bonded_atoms[2]) or \
                                             _indistinguishable(non_ring_bonded_atoms[0], non_ring_bonded_atoms[3])
                            atom1_matching = _indistinguishable(non_ring_bonded_atoms[1], non_ring_bonded_atoms[2]) or \
                                             _indistinguishable(non_ring_bonded_atoms[1], non_ring_bonded_atoms[3])
                            if not (atom0_matching and atom1_matching):
                                all_the_same = False
                else:
                    atom = ring[flipping_atom_index]
                    non_ring_bonded_atoms = [bonded_atom for bonded_atom in atom.bonds.keys()
                                             if id(bonded_atom) not in ringed_atom_ids]
                    if len(non_ring_bonded_atoms) < 2:
                        pass  # all_the_same still true
                    elif len(non_ring_bonded_atoms) == 2:
                        identical = _indistinguishable(non_ring_bonded_atoms[0], non_ring_bonded_atoms[1])
                        if not identical:
                            # flipping a tetrahedral will not work
                            all_the_same = False
                        else:
                            # having 5+ bonnds is not modeled here
                            pass
            # for even rings, check for flipping accross bonds too
            if not all_the_same and size % 2 == 0:
                all_the_same = True
                min_index = flipping_atom_index
                max_index = flipping_atom_index + size - 1
                while min_index < max_index:
                    # ensure the two atoms are different. use mod size to loop to the start
                    # of the list when index out of bounds
                    if not _indistinguishable(ring[min_index % size], ring[max_index % size]):
                        all_the_same = False
                        break
                    min_index += 1
                    max_index += -1

            if all_the_same:
                symmetry_number *= 2
                break
    return symmetry_number


################################################################################


def _indistinguishable(atom1, atom2):
    """
    Determine if two atoms are feasibly indistinguishable based on connections
    to nearest neighbors.
    """
    if (not atom1.equivalent(atom2)
            or atom1.connectivity1 != atom2.connectivity1
            or atom1.connectivity2 != atom2.connectivity2
            or atom1.connectivity3 != atom2.connectivity3):
        return False

    bond_orders_1 = [bond.order for bond in atom1.bonds.values()].sort()
    bond_orders_2 = [bond.order for bond in atom2.bonds.values()].sort()

    if bond_orders_1 != bond_orders_2:
        return False

    bonds_1 = list(atom1.bonds.items())
    bonds_2 = list(atom2.bonds.items())

    for i, (neighbor1, bond1) in enumerate(bonds_1):
        for j, (neighbor2, bond2) in enumerate(bonds_2):
            if bond1.equivalent(bond2) and neighbor1.equivalent(neighbor2):
                del bonds_2[j]
                break
        else:
            return False

    # We were able to match up all neighbors
    return True


def calculate_symmetry_number(molecule):
    """
    Return the symmetry number for the structure. The symmetry number
    includes both external and internal modes.
    """
    symmetry_number = 1

    for atom in molecule.vertices:
        if not molecule.is_atom_in_cycle(atom):
            symmetry_number *= calculate_atom_symmetry_number(molecule, atom)

    for atom1 in molecule.vertices:
        for atom2 in list(atom1.edges):  # Make a copy of the list of neighbors since we modify the dictionary
            if (molecule.vertices.index(atom1) < molecule.vertices.index(atom2) and
                    not molecule.is_bond_in_cycle(atom1.edges[atom2])):
                symmetry_number *= calculate_bond_symmetry_number(molecule, atom1, atom2)

    symmetry_number *= calculate_axis_symmetry_number(molecule)

    if molecule.is_cyclic():
        symmetry_number *= calculate_cyclic_symmetry_number(molecule)

    return symmetry_number
