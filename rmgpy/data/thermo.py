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

"""
import itertools
import logging
import math
import os.path
import re
import time
from copy import deepcopy

import numpy as np

import rmgpy.constants as constants
import rmgpy.molecule
import rmgpy.quantity
from rmgpy.data.base import Database, Entry, make_logic_node, DatabaseError
from rmgpy.ml.estimator import MLEstimator
from rmgpy.molecule import Molecule, Bond, Group
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.data.surface import MetalDatabase
from rmgpy import settings
from rmgpy.molecule.fragment import Fragment
from rmgpy.data.surface import MetalDatabase
from rmgpy import settings

#: This dictionary is used to add multiplicity to species label
_multiplicity_labels = {1: 'S', 2: 'D', 3: 'T', 4: 'Q', 5: 'V'}


################################################################################

def save_entry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the thermo
    database to the file object `f`.
    """

    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    f.write('    label = "{0}",\n'.format(entry.label))

    if isinstance(entry.item, Molecule):
        f.write('    molecule = \n')
        f.write('"""\n')
        f.write(entry.item.to_adjacency_list(remove_h=False))
        f.write('""",\n')
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.to_adjacency_list())
        f.write('""",\n')
    else:
        f.write('    group = "{0}",\n'.format(entry.item))

    if isinstance(entry.data, ThermoData):
        f.write('    thermo = ThermoData(\n')
        f.write('        Tdata = {0!r},\n'.format(entry.data.Tdata))
        f.write('        Cpdata = {0!r},\n'.format(entry.data.Cpdata))
        f.write('        H298 = {0!r},\n'.format(entry.data.H298))
        f.write('        S298 = {0!r},\n'.format(entry.data.S298))
        if entry.data.Tmin is not None:
            f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None:
            f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, Wilhoit):
        f.write('    thermo = Wilhoit(\n')
        f.write('        cp0 = {0!r},\n'.format(entry.data.cp0))
        f.write('        cpInf = {0!r},\n'.format(entry.data.cpInf))
        f.write('        a0 = {0:g},\n'.format(entry.data.a0))
        f.write('        a1 = {0:g},\n'.format(entry.data.a1))
        f.write('        a2 = {0:g},\n'.format(entry.data.a2))
        f.write('        a3 = {0:g},\n'.format(entry.data.a3))
        f.write('        B = {0!r},\n'.format(entry.data.B))
        f.write('        H0 = {0!r},\n'.format(entry.data.H0))
        f.write('        S0 = {0!r},\n'.format(entry.data.S0))
        if entry.data.Tmin is not None:
            f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None:
            f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, NASA):
        f.write('    thermo = NASA(\n')
        f.write('        polynomials = [\n')
        for poly in entry.data.polynomials:
            f.write('            {0!r},\n'.format(poly))
        f.write('        ],\n')
        if entry.data.Tmin is not None:
            f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None:
            f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.E0 is not None:
            f.write('        E0 = {0!r},\n'.format(entry.data.E0))
        if entry.data.Cp0 is not None:
            f.write('        Cp0 = {0!r},\n'.format(entry.data.Cp0))
        if entry.data.CpInf is not None:
            f.write('        CpInf = {0!r},\n'.format(entry.data.CpInf))
        f.write('    ),\n')
    else:
        f.write('    thermo = {0!r},\n'.format(entry.data))

    if entry.reference is not None:
        f.write('    reference = {0!r},\n'.format(entry.reference))
    if entry.reference_type != "":
        f.write('    referenceType = "{0}",\n'.format(entry.reference_type))
    f.write(f'    shortDesc = """{entry.short_desc.strip()}""",\n')
    f.write(f'    longDesc = \n"""\n{entry.long_desc.strip()}\n""",\n')
    if entry.rank:
        f.write("    rank = {0},\n".format(entry.rank))

    if entry.metal:
        f.write('    metal = "{0}",\n'.format(entry.metal))
    if entry.facet:
        f.write('    facet = "{0}",\n'.format(entry.facet))
    if entry.site:
        f.write('    site = "{0}",\n'.format(entry.site))

    f.write(')\n\n')


def generate_old_library_entry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    if isinstance(data, ThermoData):
        return '{0:9g} {1:9g} {2:9g} {3:9g} {4:9g} {5:9g} {6:9g} {7:9g} {8:9g} {9:9g} {10:9g} {11:9g}'.format(
            data.H298.value_si / 4184.,
            data.S298.value_si / 4.184,
            data.Cpdata.value_si[0] / 4.184,
            data.Cpdata.value_si[1] / 4.184,
            data.Cpdata.value_si[2] / 4.184,
            data.Cpdata.value_si[3] / 4.184,
            data.Cpdata.value_si[4] / 4.184,
            data.Cpdata.value_si[5] / 4.184,
            data.Cpdata.value_si[6] / 4.184,
            data.H298.uncertainty / 4184.,
            data.S298.uncertainty / 4.184,
            max(data.Cpdata.uncertainty) / 4.184,
        )
    elif isinstance(data, str):
        return data
    else:
        return '{0:9g} {1:9g} {2:9g} {3:9g} {4:9g} {5:9g} {6:9g} {7:9g} {8:9g} {9:9g} {10:9g} {11:9g}'.format(
            data.get_enthalpy(298) / 4184.,
            data.get_entropy(298) / 4.184,
            data.get_heat_capacity(300) / 4.184,
            data.get_heat_capacity(400) / 4.184,
            data.get_heat_capacity(500) / 4.184,
            data.get_heat_capacity(600) / 4.184,
            data.get_heat_capacity(800) / 4.184,
            data.get_heat_capacity(1000) / 4.184,
            data.get_heat_capacity(1500) / 4.184,
            0,
            0,
            0,
        )


def process_old_library_entry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    return ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([float(d) for d in data[2:9]], "cal/(mol*K)", "+|-", float(data[11])),
        H298=(float(data[0]), "kcal/mol", "+|-", float(data[9])),
        S298=(float(data[1]), "cal/(mol*K)", "+|-", float(data[10])),
    )


def add_thermo_data(thermo_data1, thermo_data2, group_additivity=False, verbose=False):
    """
    Add the thermodynamic data `thermo_data2` to the data `thermo_data1`,
    and return `thermo_data1`.

    If `group_additivity` is True, append comments related to group additivity estimation
    If `verbose` is False, omit the comments from a "zero entry", whose H298, S298, and Cp are all 0.
    If `verbose` is True, or thermo_data2 is not a zero entry, add thermo_data2.comment to thermo_data1.comment.
    """
    if (len(thermo_data1.Tdata.value_si) != len(thermo_data2.Tdata.value_si) or
            any([T1 != T2 for T1, T2 in zip(thermo_data1.Tdata.value_si, thermo_data2.Tdata.value_si)])):
        raise ValueError('Cannot add these ThermoData objects due to their having different temperature points.')

    for i in range(thermo_data1.Tdata.value_si.shape[0]):
        thermo_data1.Cpdata.value_si[i] += thermo_data2.Cpdata.value_si[i]
    thermo_data1.H298.value_si += thermo_data2.H298.value_si
    thermo_data1.S298.value_si += thermo_data2.S298.value_si

    test_zero = sum(abs(value) for value in
                    [thermo_data2.H298.value_si, thermo_data2.S298.value_si] + thermo_data2.Cpdata.value_si.tolist())
    # Used to check if all of the entries in thermo_data2 are zero

    if group_additivity:
        if verbose or test_zero != 0:
            # If verbose==True or test_zero!=0, add thermo_data2.comment to thermo_data1.comment.
            if thermo_data1.comment:
                thermo_data1.comment += ' + {0}'.format(thermo_data2.comment)
            else:
                thermo_data1.comment = 'Thermo group additivity estimation: ' + thermo_data2.comment

    return thermo_data1


def remove_thermo_data(thermo_data1, thermo_data2, group_additivity=False, verbose=False):
    """
    Remove the thermodynamic data `thermo_data2` from the data `thermo_data1`,
    and return `thermo_data1`.
    If `verbose` is True, append ' - thermo_data2.comment' to the thermo_data1.comment.
    If `verbose` is False, remove the thermo_data2.comment from the thermo_data1.comment.
    """
    if (len(thermo_data1.Tdata.value_si) != len(thermo_data2.Tdata.value_si) or
            any([T1 != T2 for T1, T2 in zip(thermo_data1.Tdata.value_si, thermo_data2.Tdata.value_si)])):
        raise ValueError('Cannot take the difference between these ThermoData objects due to their having different '
                         'temperature points.')

    for i in range(thermo_data1.Tdata.value_si.shape[0]):
        thermo_data1.Cpdata.value_si[i] -= thermo_data2.Cpdata.value_si[i]
    thermo_data1.H298.value_si -= thermo_data2.H298.value_si
    thermo_data1.S298.value_si -= thermo_data2.S298.value_si

    if group_additivity:
        if verbose:
            thermo_data1.comment += ' - {0}'.format(thermo_data2.comment)
        else:
            thermo_data1.comment = re.sub(re.escape(' + ' + thermo_data2.comment), '', thermo_data1.comment, 1)
    return thermo_data1


def average_thermo_data(thermo_data_list=None):
    """
    Average a list of ThermoData values together.
    Sets uncertainty values to be the approximately the 95% confidence interval, equivalent to
    2 standard deviations calculated using the sample standard variance:
    
    Uncertainty = 2s
    s = sqrt( sum(abs(x - x.mean())^2) / N - 1) where N is the number of values averaged
    
    Note that uncertainties are only computed when number of values is greater than 1.
    """
    if thermo_data_list is None:
        thermo_data_list = []

    num_values = len(thermo_data_list)

    if num_values == 0:
        raise ValueError('No thermo data values were inputted to be averaged.')
    else:
        logging.debug('Averaging thermo data over {0} value(s).'.format(num_values))

        if num_values == 1:
            return deepcopy(thermo_data_list[0])

        else:
            averaged_thermo_data = deepcopy(thermo_data_list[0])
            for thermo_data in thermo_data_list[1:]:
                averaged_thermo_data = add_thermo_data(averaged_thermo_data, thermo_data)

            for i in range(averaged_thermo_data.Tdata.value_si.shape[0]):
                averaged_thermo_data.Cpdata.value_si[i] /= num_values

                cp_data = [thermo_data.Cpdata.value_si[i] for thermo_data in thermo_data_list]
                averaged_thermo_data.Cpdata.uncertainty[i] = 2 * np.std(cp_data, ddof=1)

            h_data = [thermo_data.H298.value_si for thermo_data in thermo_data_list]
            averaged_thermo_data.H298.value_si /= num_values
            averaged_thermo_data.H298.uncertainty_si = 2 * np.std(h_data, ddof=1)

            s_data = [thermo_data.S298.value_si for thermo_data in thermo_data_list]
            averaged_thermo_data.S298.value_si /= num_values
            averaged_thermo_data.S298.uncertainty_si = 2 * np.std(s_data, ddof=1)
            return averaged_thermo_data


def common_atoms(cycle1, cycle2):
    """
    INPUT: two cycles with type: list of atoms
    OUTPUT: a set of common atoms
    """
    set1 = set(cycle1)
    set2 = set(cycle2)
    return set1.intersection(set2)


def combine_cycles(cycle1, cycle2):
    """
    INPUT: two cycles with type: list of atoms
    OUTPUT: a combined cycle with type: list of atoms
    """
    set1 = set(cycle1)
    set2 = set(cycle2)
    return list(set1.union(set2))


def is_aromatic_ring(submol):
    """
    This method takes a monoring submol (Molecule initialized with a list of atoms containing just 
    the ring), and check if it is a aromatic ring.
    """
    ring_size = len(submol.atoms)
    if ring_size not in [5, 6]:
        return False
    for ring_atom in submol.atoms:
        for bonded_atom, bond in ring_atom.edges.items():
            if bonded_atom in submol.atoms:
                if not bond.is_benzene():
                    return False
    return True


def is_bicyclic(polyring):
    """
    Given a polyring (a list of `Atom`s)
    returns True if it's a bicyclic, False otherwise
    """
    submol, _ = convert_ring_to_sub_molecule(polyring)
    sssr = submol.get_smallest_set_of_smallest_rings()

    return len(sssr) == 2


def find_aromatic_bonds_from_sub_molecule(submol):
    """
    This method finds all the aromatic bonds within a input submolecule and 
    returns a set of unique aromatic bonds
    """

    aromatic_bonds = []
    for atom in submol.atoms:
        bonds = submol.get_bonds(atom)
        for atom_j in bonds:
            if atom_j in submol.atoms:
                bond = bonds[atom_j]
                if bond.is_benzene():
                    aromatic_bonds.append(bond)
    return set(aromatic_bonds)


def convert_ring_to_sub_molecule(ring):
    """
    This function takes a ring structure (can either be monoring or polyring) to create a new 
    submolecule with newly deep copied atoms

    Outputted submolecules may have incomplete valence and may cause errors with some Molecule.methods(), such
    as update_atomtypes() or update(). In the future we may consider using groups for the sub-molecules.
    """

    atoms_mapping = {}
    for atom in ring:
        atoms_mapping[atom] = atom.copy()  # this copy is deep copy of origin atom with empty edges

    mol0 = Molecule(atoms=list(atoms_mapping.values()))

    for atom in ring:
        for bonded_atom, bond in atom.edges.items():
            if bonded_atom in ring:
                if not mol0.has_bond(atoms_mapping[atom], atoms_mapping[bonded_atom]):
                    mol0.add_bond(Bond(atoms_mapping[atom], atoms_mapping[bonded_atom], order=bond.order))

    mol0.update_multiplicity()
    mol0.update_connectivity_values()
    return mol0, atoms_mapping


def combine_two_rings_into_sub_molecule(ring1, ring2):
    """
    This function combines 2 rings (with common atoms) to create a new 
    submolecule with newly deep copied atoms
    """

    assert len(common_atoms(ring1, ring2)) > 0, "The two input rings don't have common atoms."

    atoms_mapping = {}
    for atom in ring1 + ring2:
        if atom not in atoms_mapping:
            atoms_mapping[atom] = atom.copy()

    mol0 = Molecule(atoms=list(atoms_mapping.values()))

    for atom in ring1:
        for bonded_atom, bond in atom.edges.items():
            if bonded_atom in ring1:
                if not mol0.has_bond(atoms_mapping[atom], atoms_mapping[bonded_atom]):
                    mol0.add_bond(Bond(atoms_mapping[atom], atoms_mapping[bonded_atom], order=bond.order))

    for atom in ring2:
        for bonded_atom, bond in atom.edges.items():
            if bonded_atom in ring2:
                if not mol0.has_bond(atoms_mapping[atom], atoms_mapping[bonded_atom]):
                    mol0.add_bond(Bond(atoms_mapping[atom], atoms_mapping[bonded_atom], order=bond.order))

    mol0.update_multiplicity()
    mol0.update_connectivity_values()

    return mol0, atoms_mapping


def get_copy_for_one_ring(ring):
    """
    Make a copy of a single ring from a molecule.

    Returns a list of atoms.
    """
    atoms_mapping = convert_ring_to_sub_molecule(ring)[1]

    ring_copy = [atoms_mapping[atom] for atom in ring]

    return ring_copy


def get_copy_from_two_rings_with_common_atoms(ring1, ring2):
    """
    Make a copy of a two rings from a molecule and also generates the merged ring.

    Returns a copy of ring1, a copy of ring2, and the merged rings, each as a list of atoms.
    """
    merged_ring, atoms_mapping = combine_two_rings_into_sub_molecule(ring1, ring2)

    ring1_copy = [atoms_mapping[atom] for atom in ring1]
    ring2_copy = [atoms_mapping[atom] for atom in ring2]

    return ring1_copy, ring2_copy, merged_ring


def is_ring_partial_matched(ring, matched_group):
    """
    An example of ring partial match is tricyclic ring is matched by a bicyclic group
    usually because of not enough data in polycyclic tree. The method takes a matched group 
    returned from descend_tree and the ring (a list of non-hydrogen atoms in the ring)
    """
    # if matched group has less atoms than the target ring
    # it's surely a partial match
    if len(ring) > len(matched_group.atoms):
        return True
    else:
        submol_ring, _ = convert_ring_to_sub_molecule(ring)
        sssr = submol_ring.get_smallest_set_of_smallest_rings()
        sssr_grp = matched_group.get_smallest_set_of_smallest_rings()
        if sorted([len(sr) for sr in sssr]) == sorted([len(sr_grp) for sr_grp in sssr_grp]):
            return False
        else:
            return True


def bicyclic_decomposition_for_polyring(polyring):
    """
    Decompose a polycyclic ring into all possible bicyclic combinations: `bicyclics_merged_from_ring_pair`
    and return a `ring_occurances_dict` that contains all single ring tuples as keys and the number of times
    they appear each bicyclic submolecule.  These bicyclic and single rings are used 
    later in the heuristic polycyclic thermo algorithm.
    """

    submol, _ = convert_ring_to_sub_molecule(polyring)
    sssr = submol.get_deterministic_sssr()

    ring_pair_with_common_atoms_list = []
    ring_occurances_dict = {}

    # Initialize ringOccuranceDict
    for ring in sssr:
        ring_occurances_dict[tuple(ring)] = 0

    ring_num = len(sssr)
    for i in range(ring_num):
        for j in range(i + 1, ring_num):
            if common_atoms(sssr[i], sssr[j]):
                # Copy the SSSR's again because these ones are going to be merged into bicyclics
                # and manipulated (aromatic bonds have to be screened and changed to single if needed)
                sssr_i, sssr_j, merged_ring = get_copy_from_two_rings_with_common_atoms(sssr[i], sssr[j])
                ring_pair_with_common_atoms_list.append([sssr_i, sssr_j, merged_ring])
                # Save the single ring SSSRs that appear in bicyclics using the original copy
                # because they will be manipulated (differently) in _add_poly_ring_correction_thermo_data_from_heuristic
                ring_occurances_dict[tuple(sssr[i])] += 1
                ring_occurances_dict[tuple(sssr[j])] += 1

    bicyclics_merged_from_ring_pair = []
    # pre-process 2-ring cores
    for ringA, ringB, merged_ring in ring_pair_with_common_atoms_list:
        submol_a = Molecule(atoms=ringA)
        submol_b = Molecule(atoms=ringB)
        is_a_aromatic = is_aromatic_ring(submol_a)
        is_b_aromatic = is_aromatic_ring(submol_b)
        # if ringA and ringB are both aromatic or not aromatic
        # don't need to do anything extra
        if is_a_aromatic and is_b_aromatic:
            pass
        elif not is_a_aromatic and not is_b_aromatic:
            aromatic_bonds_in_a = find_aromatic_bonds_from_sub_molecule(submol_a)
            for aromaticBond_inA in aromatic_bonds_in_a:
                aromaticBond_inA.set_order_num(1)

            aromatic_bonds_in_b = find_aromatic_bonds_from_sub_molecule(submol_b)
            for aromaticBond_inB in aromatic_bonds_in_b:
                aromaticBond_inB.set_order_num(1)
        elif is_a_aromatic:
            aromatic_bonds_in_b = find_aromatic_bonds_from_sub_molecule(submol_b)
            for aromaticBond_inB in aromatic_bonds_in_b:
                # Make sure the aromatic bond in ringB is in ringA, and both ringB atoms are in ringA 
                # If so, preserve the B bond status, otherwise change to single bond order
                if ((aromaticBond_inB.atom1 in submol_a.atoms) and
                        (aromaticBond_inB.atom2 in submol_a.atoms) and
                        (submol_a.has_bond(aromaticBond_inB.atom1, aromaticBond_inB.atom2))):
                    pass
                else:
                    aromaticBond_inB.set_order_num(1)
        else:
            aromatic_bonds_in_a = find_aromatic_bonds_from_sub_molecule(submol_a)
            for aromaticBond_inA in aromatic_bonds_in_a:
                if ((aromaticBond_inA.atom1 in submol_b.atoms) and
                        (aromaticBond_inA.atom2 in submol_b.atoms) and
                        (submol_b.has_bond(aromaticBond_inA.atom1, aromaticBond_inA.atom2))):
                    pass
                else:
                    aromaticBond_inA.set_order_num(1)
        merged_ring.saturate_unfilled_valence(update=True)
        bicyclics_merged_from_ring_pair.append(merged_ring)

    return bicyclics_merged_from_ring_pair, ring_occurances_dict


def split_bicyclic_into_single_rings(bicyclic_submol):
    """
    Splits a given bicyclic submolecule into two individual single 
    ring submolecules (a list of `Molecule`s ).
    """
    sssr = bicyclic_submol.get_deterministic_sssr()

    return [convert_ring_to_sub_molecule(sssr[0])[0],
            convert_ring_to_sub_molecule(sssr[1])[0]]


def saturate_ring_bonds(ring_submol):
    """
    Given a ring submolelcule (`Molecule`), makes a deep copy and converts non-single bonds 
    into single bonds, returns a new saturated submolecule (`Molecule`)
    """
    atoms_mapping = {}
    for atom in ring_submol.atoms:
        if atom not in atoms_mapping:
            atoms_mapping[atom] = atom.copy()

    mol0 = Molecule(atoms=list(atoms_mapping.values()))

    already_saturated = True
    for atom in ring_submol.atoms:
        for bonded_atom, bond in atom.edges.items():
            if bonded_atom in ring_submol.atoms:
                if bond.order > 1.0 and not bond.is_benzene():
                    already_saturated = False
                if not mol0.has_bond(atoms_mapping[atom], atoms_mapping[bonded_atom]):
                    bond_order = 1.0
                    if bond.is_benzene():
                        bond_order = 1.5
                    mol0.add_bond(Bond(atoms_mapping[atom], atoms_mapping[bonded_atom], order=bond_order))

    mol0.saturate_unfilled_valence()
    mol0.update_atomtypes()
    mol0.update_multiplicity()
    mol0.update_connectivity_values()
    return mol0, already_saturated


################################################################################

class ThermoDepository(Database):
    """
    A class for working with the RMG thermodynamics depository.
    """

    def __init__(self, label='', name='', short_desc='', long_desc='', metal=None, site=None, facet=None):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc, metal=metal, site=site, facet=facet)

    def load_entry(self, index, label, molecule, thermo, reference=None, referenceType='', shortDesc='', longDesc='',
                   rank=None, metal=None, site=None, facet=None):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """
        entry = Entry(
            index=index,
            label=label,
            item=Molecule().from_adjacency_list(molecule),
            data=thermo,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            rank=rank,
            metal=metal,
            site=site,
            facet=facet,
        )
        self.entries[label] = entry
        return entry

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)


################################################################################

class ThermoLibrary(Database):
    """
    A class for working with a RMG thermodynamics library.
    """

    def __init__(self, label='', name='', solvent=None, short_desc='', long_desc='', metal=None, site=None, facet=None):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc,
                          metal=metal, site=site, facet=facet)

    def load_entry(self,
                   index,
                   label,
                   molecule,
                   thermo,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   rank=None,
                   metal=None,
                   facet=None,
                   site=None,
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """

        try:
            molecule = Molecule().from_adjacency_list(molecule)
        except TypeError:
            molecule = Fragment().from_adjacency_list(molecule)

        # Internal checks for adding entry to the thermo library
        if label in list(self.entries.keys()):
            raise DatabaseError('Found a duplicate molecule with label {0} in the thermo library {1}. '
                                'Please correct your library.'.format(label, self.name))

        for entry in self.entries.values():
            if molecule.is_isomorphic(entry.item):
                if molecule.multiplicity == entry.item.multiplicity:
                    raise DatabaseError('Adjacency list and multiplicity of {0} matches that of '
                                        'existing molecule {1} in thermo library {2}. Please '
                                        'correct your library.'.format(label, entry.label, self.name))

        self.entries[label] = Entry(
            index=index,
            label=label,
            item=molecule,
            data=thermo,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            rank=rank,
            metal=metal,
            facet=facet,
            site=site,
        )

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)

    def generate_old_library_entry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        return generate_old_library_entry(data)

    def process_old_library_entry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return process_old_library_entry(data)


################################################################################

class ThermoGroups(Database):
    """
    A class for working with an RMG thermodynamics group additivity database.
    """

    def __init__(self, label='', name='', short_desc='', long_desc='', metal=None, site=None, facet=None):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc,
                          metal=metal, site=site, facet=facet)

    def load_entry(self,
                   index,
                   label,
                   group,
                   thermo,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   rank=None,
                   metal=None,
                   facet=None,
                   site=None,
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """

        if (group[0:3].upper() == 'OR{' or
                group[0:4].upper() == 'AND{' or
                group[0:7].upper() == 'NOT OR{' or
                group[0:8].upper() == 'NOT AND{'):
            item = make_logic_node(group)
        else:
            item = Group().from_adjacency_list(group)
        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=thermo,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            rank=rank,
            metal=metal,
            facet=facet,
            site=site,
        )

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)

    def generate_old_library_entry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """

        return generate_old_library_entry(data)

    def process_old_library_entry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return process_old_library_entry(data)

    def copy_data(self, source, destination):
        """
        This method copys the ThermoData object and all meta data
        from source to destination
        Args:
            source: The entry for which data is being copied
            destination: The entry for which data is being overwritten

        """
        destination.data = source.data
        destination.reference = source.reference
        destination.short_desc = source.short_desc
        destination.long_desc = source.long_desc
        destination.rank = source.rank
        destination.reference_type = source.reference_type
        destination.metal = source.metal
        destination.facet = source.facet
        destination.site = source.site

    def remove_group(self, group_to_remove):
        """
        Removes a group that is in a tree from the database. For thermo
        groups we also, need to re-point any unicode thermo_data that may
        have pointed to the entry.

        Returns the removed group
        """

        # First call base class method
        Database.remove_group(self, group_to_remove)

        parent_r = group_to_remove.parent

        # look for other pointers that point toward entry
        for entry in self.entries.values():
            if isinstance(entry.data, str):
                if entry.data == group_to_remove.label:
                    # if the entryToRemove.data is also a pointer, then copy
                    if isinstance(group_to_remove.data, str):
                        entry.data = group_to_remove.data
                    # if the parent points toward entry and the data is
                    # not a base string, we need to copy the data to the parent
                    elif entry is parent_r:
                        self.copy_data(group_to_remove, parent_r)
                    # otherwise, point toward entryToRemove's parent
                    else:
                        entry.data = str(parent_r.label)

        return group_to_remove


################################################################################

class ThermoDatabase(object):
    """
    A class for working with the RMG thermodynamics database.
    """

    def __init__(self):
        self.depository = {}
        self.libraries = {}
        self.surface = {}
        self.groups = {}
        self.adsorption_groups = "adsorptionPt111"
        self.library_order = []
        self.local_context = {
            'ThermoData': ThermoData,
            'Wilhoit': Wilhoit,
            'NASAPolynomial': NASAPolynomial,
            'NASA': NASA,
        }
        self.global_context = {}

        # Use Pt111 binding energies as default
        self.binding_energies = {
            'H': rmgpy.quantity.Energy(-2.75368,'eV/molecule'),
            'C': rmgpy.quantity.Energy(-7.02516,'eV/molecule'),
            'N': rmgpy.quantity.Energy(-4.63225,'eV/molecule'),
            'O': rmgpy.quantity.Energy(-3.81153,'eV/molecule'),
        }

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoDatabase object.
        """
        d = {
            'depository': self.depository,
            'libraries': self.libraries,
            'groups': self.groups,
            'library_order': self.library_order,
            'surface' : self.surface,
        }
        return ThermoDatabase, (), d

    def __setstate__(self, d):
        """
        A helper function used when unpickling a ThermoDatabase object.
        """
        self.depository = d['depository']
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.library_order = d['library_order']
        self.surface = d['surface']

    def load(self, path, libraries=None, depository=True, surface=False):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        if depository:
            self.load_depository(os.path.join(path, 'depository'))
        else:
            self.depository = {}
        self.load_libraries(os.path.join(path, 'libraries'), libraries)
        self.load_groups(os.path.join(path, 'groups'))
        if surface:
            self.load_surface()

    def load_depository(self, path):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.depository = {
            'stable': ThermoDepository().load(os.path.join(path, 'stable.py'),
                                              self.local_context, self.global_context),
            'radical': ThermoDepository().load(os.path.join(path, 'radical.py'),
                                               self.local_context, self.global_context)
        }

    def load_libraries(self, path, libraries=None):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        
        If no libraries are given, all are loaded.
        """
        self.libraries = {}
        self.library_order = []
        if libraries is None:
            for (root, dirs, files) in os.walk(os.path.join(path)):
                for f in files:
                    name, ext = os.path.splitext(f)
                    if ext.lower() == '.py':
                        logging.info('Loading thermodynamics library from {0} in {1}...'.format(f, root))
                        library = ThermoLibrary()
                        library.load(os.path.join(root, f), self.local_context, self.global_context)
                        library.label = os.path.splitext(f)[0]
                        self.libraries[library.label] = library
                        self.library_order.append(library.label)

        else:
            for libraryName in libraries:
                f = f'{libraryName}.py'
                if os.path.isfile(libraryName):
                    logging.info(f'Loading thermodynamics library from an external location: {libraryName}..')
                    library = ThermoLibrary()
                    library.load(libraryName, self.local_context, self.global_context)
                    library.label = os.path.splitext(os.path.split(libraryName)[-1])[0]
                    self.libraries[library.label] = library
                    self.library_order.append(library.label)
                elif os.path.exists(os.path.join(path, f)):
                    logging.info(f'Loading thermodynamics library from {f} in {path}...')
                    library = ThermoLibrary()
                    library.load(os.path.join(path, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.library_order.append(library.label)
                else:
                    raise DatabaseError('Library {} not found in {}... Please check if your library is '
                                        'correctly placed'.format(libraryName, path))

    def load_surface(self):
        """
        Load the metal database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        MetalDB = MetalDatabase()
        MetalDB.load(os.path.join(settings['database.directory'], 'surface'))

        self.surface = {
            'metal': MetalDB
        }

    def load_groups(self, path):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        logging.info('Loading thermodynamics group database from {0}...'.format(path))
        categories = [
            'group',
            'ring',
            'radical',
            'polycyclic',
            'other',
            'longDistanceInteraction_cyclic',
            'longDistanceInteraction_noncyclic',
            'adsorptionPt111',
            'adsorptionLi'
        ]
        # categories.append(self.adsorption_groups)
        self.groups = {
            category: ThermoGroups(label=category).load(os.path.join(path, category + '.py'),
                                                        self.local_context, self.global_context)
            for category in categories
        }

        self.record_ring_generic_nodes()
        self.record_polycylic_generic_nodes()

    def save(self, path):
        """
        Save the thermo database to the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path):
            os.mkdir(path)
        self.save_depository(os.path.join(path, 'depository'))
        self.save_libraries(os.path.join(path, 'libraries'))
        self.save_groups(os.path.join(path, 'groups'))
        self.save_surface(os.path.join(path, 'surface'))

    def save_depository(self, path):
        """
        Save the thermo depository to the given `path` on disk, where `path`
        points to the top-level folder of the thermo depository.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for depo in self.depository.keys():
            self.depository[depo].save(os.path.join(path, depo + '.py'))

    def save_libraries(self, path):
        """
        Save the thermo libraries to the given `path` on disk, where `path`
        points to the top-level folder of the thermo libraries.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for library in self.libraries.values():
            library.save(os.path.join(path, '{0}.py'.format(library.label)))

    def save_groups(self, path):
        """
        Save the thermo groups to the given `path` on disk, where `path`
        points to the top-level folder of the thermo groups.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for group in self.groups.keys():
            self.groups[group].save(os.path.join(path, group + '.py'))

    def save_surface(self, path):
        """
        Save the metal library to the given `path` on disk, where `path`
        points to the top-level folder of the metal library.
        """

        if not os.path.exists(path):
            os.mkdir(path)
        for library in self.surface.keys():
            self.surface[library].save(os.path.join(path, library + '.py'))

    def load_old(self, path):
        """
        Load the old RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # The old database does not have a depository, so create an empty one
        self.depository = {}
        self.depository['stable'] = ThermoDepository(label='stable', name='Stable Molecules')
        self.depository['radical'] = ThermoDepository(label='radical', name='Radical Molecules')

        for (root, dirs, files) in os.walk(os.path.join(path, 'thermo_libraries')):
            if (os.path.exists(os.path.join(root, 'Dictionary.txt')) and
                    os.path.exists(os.path.join(root, 'Library.txt'))):
                library = ThermoLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.load_old(
                    dictstr=os.path.join(root, 'Dictionary.txt'),
                    treestr='',
                    libstr=os.path.join(root, 'Library.txt'),
                    num_parameters=12,
                    num_labels=1,
                    pattern=False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups = {}
        self.groups['group'] = ThermoGroups(label='group', name='Functional Group Additivity Values').load_old(
            dictstr=os.path.join(path, 'thermo_groups', 'Group_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', 'Group_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', 'Group_Library.txt'),
            num_parameters=12,
            num_labels=1,
            pattern=True,
        )
        self.groups['gauche'] = ThermoGroups(label='gauche', name='Gauche Interaction Corrections').load_old(
            dictstr=os.path.join(path, 'thermo_groups', 'Gauche_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', 'Gauche_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', 'Gauche_Library.txt'),
            num_parameters=12,
            num_labels=1,
            pattern=True,
        )
        self.groups['int15'] = ThermoGroups(label='int15', name='1,5-Interaction Corrections').load_old(
            dictstr=os.path.join(path, 'thermo_groups', '15_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', '15_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', '15_Library.txt'),
            num_parameters=12,
            num_labels=1,
            pattern=True,
        )
        self.groups['radical'] = ThermoGroups(label='radical', name='Radical Corrections').load_old(
            dictstr=os.path.join(path, 'thermo_groups', 'Radical_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', 'Radical_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', 'Radical_Library.txt'),
            num_parameters=12,
            num_labels=1,
            pattern=True,
        )
        self.groups['ring'] = ThermoGroups(label='ring', name='Ring Corrections').load_old(
            dictstr=os.path.join(path, 'thermo_groups', 'Ring_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', 'Ring_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', 'Ring_Library.txt'),
            num_parameters=12,
            num_labels=1,
            pattern=True,
        )
        self.groups['polycyclic'] = ThermoGroups(label='other', name='Polycyclic Ring Corrections').load_old(
            dictstr=os.path.join(path, 'thermo_groups', 'Polycyclic_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', 'Polycyclic_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', 'Polycyclic_Library.txt'),
            num_parameters=12,
            num_labels=1,
            pattern=True,
        )
        self.groups['other'] = ThermoGroups(label='other', name='Other Corrections').load_old(
            dictstr=os.path.join(path, 'thermo_groups', 'Other_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', 'Other_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', 'Other_Library.txt'),
            num_parameters=12,
            num_labels=1,
            pattern=True,
        )

    def prune_heteroatoms(self, allowed=None):
        """
        Remove all species from thermo libraries that contain atoms other than those allowed.
        
        This is useful before saving the database for use in RMG-Java
        """
        if allowed is None:
            allowed = ['C', 'H', 'O', 'S']
        allowed_elements = [rmgpy.molecule.element.get_element(label) for label in allowed]
        for library in self.libraries.values():
            logging.info("Removing hetoroatoms from thermo library '{0}'".format(library.name))
            to_delete = []
            for entry in library.entries.values():
                for atom in entry.item.atoms:
                    if atom.element not in allowed_elements:
                        to_delete.append(entry.label)
                        break
            for label in to_delete:
                logging.info(" {0}".format(label))
                library.entries.pop(label)

    def save_old(self, path):
        """
        Save the old RMG thermo database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """

        # Depository not used in old database, so it is not saved

        libraries_path = os.path.join(path, 'thermo_libraries')
        if not os.path.exists(libraries_path):
            os.mkdir(libraries_path)
        for library in self.libraries.values():
            library_path = os.path.join(libraries_path, library.label)
            if not os.path.exists(library_path):
                os.mkdir(library_path)
            library.save_old(
                dictstr=os.path.join(library_path, 'Dictionary.txt'),
                treestr='',
                libstr=os.path.join(library_path, 'Library.txt'),
            )

        groups_path = os.path.join(path, 'thermo_groups')
        if not os.path.exists(groups_path):
            os.mkdir(groups_path)
        self.groups['group'].save_old(
            dictstr=os.path.join(groups_path, 'Group_Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Group_Tree.txt'),
            libstr=os.path.join(groups_path, 'Group_Library.txt'),
        )
        self.groups['gauche'].save_old(
            dictstr=os.path.join(groups_path, 'Gauche_Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Gauche_Tree.txt'),
            libstr=os.path.join(groups_path, 'Gauche_Library.txt'),
        )
        self.groups['int15'].save_old(
            dictstr=os.path.join(groups_path, '15_Dictionary.txt'),
            treestr=os.path.join(groups_path, '15_Tree.txt'),
            libstr=os.path.join(groups_path, '15_Library.txt'),
        )
        self.groups['radical'].save_old(
            dictstr=os.path.join(groups_path, 'Radical_Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Radical_Tree.txt'),
            libstr=os.path.join(groups_path, 'Radical_Library.txt'),
        )
        self.groups['ring'].save_old(
            dictstr=os.path.join(groups_path, 'Ring_Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Ring_Tree.txt'),
            libstr=os.path.join(groups_path, 'Ring_Library.txt'),
        )
        self.groups['polycyclic'].save_old(
            dictstr=os.path.join(groups_path, 'Polycyclic_Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Polycyclic_Tree.txt'),
            libstr=os.path.join(groups_path, 'Polycyclic_Library.txt'),
        )
        self.groups['other'].save_old(
            dictstr=os.path.join(groups_path, 'Other_Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Other_Tree.txt'),
            libstr=os.path.join(groups_path, 'Other_Library.txt'),
        )

    def record_polycylic_generic_nodes(self):
        """
        Identify generic nodes in tree for polycyclic groups.
        Saves them as a list in the `generic_nodes` attribute
        in the polycyclic :class:`ThermoGroups` object, which
        must be pre-loaded.

        Necessary for polycyclic heuristic.
        """
        self.groups['polycyclic'].generic_nodes = ['PolycyclicRing']
        for label, entry in self.groups['polycyclic'].entries.items():
            if isinstance(entry.data, ThermoData):
                continue
            self.groups['polycyclic'].generic_nodes.append(label)

    def record_ring_generic_nodes(self):
        """
        Identify generic nodes in tree for ring groups.
        Saves them as a list in the `generic_nodes` attribute
        in the ring :class:`ThermoGroups` object, which
        must be pre-loaded.

        Necessary for polycyclic heuristic.
        """
        self.groups['ring'].generic_nodes = ['Ring']
        for label, entry in self.groups['ring'].entries.items():
            if isinstance(entry.data, ThermoData):
                continue
            self.groups['ring'].generic_nodes.append(label)

    def get_thermo_data(self, species, metal_to_scale_to=None, training_set=None):
        """
        Return the thermodynamic parameters for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via machine learning and then group additivity.
        
        The method corrects for symmetry when the molecule uses machine
        learning or group additivity. Libraries and direct QM calculations
        are already corrected.

        If either metal to scale to or from is not specified, assume the binding energies given in the input file
        
        Returns: ThermoData
        """
        from rmgpy.rmg.input import get_input

        thermo0 = self.get_thermo_data_from_libraries(species)

        if thermo0 is not None:  # was able to find thermodata in the loaded libraries
            if len(thermo0) != 3:
                raise RuntimeError("thermo0 should be a tuple (thermo_data, library, entry), not {0}".format(thermo0))
            entry = thermo0[2]
            thermo0 = thermo0[0]

            if species.contains_surface_site():
                if entry.metal is not None:
                    if entry.facet is not None:
                        db_label = entry.metal + entry.facet
                        thermo0 = self.correct_binding_energy(thermo0, species, metal_to_scale_from=db_label,
                                                              metal_to_scale_to=metal_to_scale_to)
                    else:  # no facet was given
                        thermo0 = self.correct_binding_energy(thermo0, species, metal_to_scale_from=entry.metal, metal_to_scale_to=metal_to_scale_to)
                else:  # assume the thermo came from pt 111
                    thermo0 = self.correct_binding_energy(thermo0, species, metal_to_scale_from=None, metal_to_scale_to=metal_to_scale_to)
            return thermo0

        if species.contains_surface_site():
            try:
                thermo0 = self.get_thermo_data_for_surface_species(species)
                metal_to_scale_from = self.adsorption_groups.split('adsorption')[-1]
                if metal_to_scale_from != metal_to_scale_to:
                    thermo0 = self.correct_binding_energy(thermo0, species, metal_to_scale_from=metal_to_scale_from, metal_to_scale_to=metal_to_scale_to)  # group adsorption values come from Pt111
                return thermo0
            except:
                logging.error("Error attempting to get thermo for species %s with structure \n%s", 
                    species, species.molecule[0].to_adjacency_list())
                raise

        try:
            quantum_mechanics = get_input('quantum_mechanics')
        except Exception:
            logging.debug('Quantum Mechanics DB could not be found.')
            quantum_mechanics = None

        try:
            ml_estimator, ml_settings = get_input('ml_estimator')
        except Exception:
            logging.debug('ML estimator could not be found.')
            ml_estimator, ml_settings = None, None

        if quantum_mechanics:
            try:
                original_molecule = species.molecule[0]
                if quantum_mechanics.settings.onlyCyclics and not original_molecule.is_cyclic():
                    pass
                else:  # try a QM calculation
                    if original_molecule.get_radical_count() > quantum_mechanics.settings.maxRadicalNumber:
                        # Too many radicals for direct calculation: use HBI.
                        logging.info("{0} radicals on {1} exceeds limit of {2}. Using HBI method.".format(
                            original_molecule.get_radical_count(),
                            species.label,
                            quantum_mechanics.settings.maxRadicalNumber,
                        ))

                        # Need to estimate thermo via each resonance isomer
                        thermo = []
                        for molecule in species.molecule:
                            molecule.clear_labeled_atoms()
                            # Try to see if the saturated molecule can be found in the libraries
                            tdata = self.estimate_radical_thermo_via_hbi(molecule, self.get_thermo_data_from_libraries)
                            priority = 1
                            if tdata is None:
                                # Then attempt quantum mechanics job on the saturated molecule
                                tdata = self.estimate_radical_thermo_via_hbi(molecule, quantum_mechanics.get_thermo_data)
                                priority = 2
                            if tdata is None:
                                # Fall back to group additivity
                                tdata = self.estimate_thermo_via_group_additivity(molecule)
                                priority = 3

                            thermo.append((priority, tdata.get_enthalpy(298.), molecule, tdata))

                        if len(thermo) > 1:
                            # Sort thermo first by the priority, then by the most stable H298 value
                            thermo = sorted(thermo, key=lambda x: (x[0], x[1]))
                            for i, therm in enumerate(thermo):
                                logging.debug("Resonance isomer {0} {1} gives H298={2:.0f} J/mol"
                                              "".format(i + 1, therm[2].to_smiles(), therm[1]))
                            # Save resonance isomers reordered by their thermo
                            species.molecule = [item[2] for item in thermo]
                            original_molecule = species.molecule[0]
                        thermo0 = thermo[0][3]

                        # update entropy by symmetry correction
                        thermo0.S298.value_si -= constants.R * math.log(species.get_symmetry_number())

                    else:  # Not too many radicals: do a direct calculation.
                        thermo0 = quantum_mechanics.get_thermo_data(original_molecule)  # returns None if it fails
            except ValueError as e: #rdkit fails to generate conformers 
                logging.error("Quantum Mechanics calculation failed for species: %s with ValueError: %s", species.label, e.args[0])
                logging.error("Falling back to ML (If turned on) or GAV (If not)")
                    
        if thermo0 is None:
            # First try finding stable species in libraries and using HBI
            for mol in species.molecule:
                if mol.reactive:
                    original_molecule = mol
                    break
            else:
                for mol in species.molecule:
                    logging.info(mol.to_adjacency_list())
                    logging.info('reactive = {0}'.format(mol.reactive))
                    logging.info('\n')
                raise ValueError('Could not process a species with no reactive structures')
            if original_molecule.get_radical_count() > 0:
                # If the molecule is a radical, check if any of the saturated forms are in the libraries
                # first and perform an HBI correction on them
                thermo = []
                for molecule in species.molecule:
                    if molecule.reactive:
                        molecule.clear_labeled_atoms()
                        # First see if the saturated molecule is in the libaries
                        tdata = self.estimate_radical_thermo_via_hbi(molecule, self.get_thermo_data_from_libraries)
                        if tdata:
                            thermo.append((tdata.get_enthalpy(298.), molecule, tdata))

                if thermo:
                    # Sort thermo by the most stable H298 value when choosing between thermoLibrary values
                    thermo = sorted(thermo, key=lambda x: x[0])
                    # Sort thermo by the structure reactive attribute, with `reactive=True` structures first
                    thermo.sort(key=lambda x: x[1].reactive, reverse=True)
                    for i, therm in enumerate(thermo):
                        if therm[1].reactive:
                            logging.debug("Resonance isomer {0} {1} gives H298={2:.0f} J/mol"
                                          "".format(i + 1, therm[1].to_smiles(), therm[0]))
                        else:
                            logging.debug("Non-reactive resonance isomer {0} {1} gives H298={2:.0f} J/mol"
                                          "".format(i + 1, therm[1].to_smiles(), therm[0]))
                    # Save resonance isomers reordered by their thermo
                    new_mol_list = [item[1] for item in thermo]
                    if len(new_mol_list) < len(species.molecule):
                        new_mol_list.extend([mol for mol in species.molecule if mol not in new_mol_list])
                    species.molecule = new_mol_list
                    thermo0 = thermo[0][2]

            if thermo0 is None:
                # If we still don't have thermo, use ML to estimate it, but
                # only if the molecule is made up of H, C, N, and O atoms and
                # is not a singlet carbene. ML settings are checked in
                # `self.get_thermo_data_from_ml`.
                if (ml_estimator is not None
                        and all(a.element.number in {1, 6, 7, 8} for a in species.molecule[0].atoms)
                        and species.molecule[0].get_singlet_carbene_count() == 0):
                    thermo0 = self.get_thermo_data_from_ml(species,
                                                           ml_estimator,
                                                           ml_settings)

            if thermo0 is None:
                # And lastly, resort back to group additivity to determine thermo for molecule
                thermo0 = self.get_thermo_data_from_groups(species)

            # Update entropy by symmetry correction (not included in trained ML model)
            thermo0.S298.value_si -= constants.R * math.log(species.get_symmetry_number())

        # Make sure to calculate Cp0 and CpInf if it wasn't done already
        find_cp0_and_cpinf(species, thermo0)

        # Return the resulting thermo parameters
        return thermo0

    def set_binding_energies(self, binding_energies='Pt111'):
        """
        Sets and stores the atomic binding energies specified in the input file.

        All adsorbates will be scaled to use these elemental binding energies.

        Args:
            binding_energies (dict, optional): the desired binding energies with
                elements as keys and binding energy/unit tuples (or Energy 
                quantities) as values

        Returns:
            None, stores result in self.binding_energies
        """
        
        if isinstance(binding_energies, str):
            if not self.surface:
                self.load_surface()
            binding_energies = self.surface['metal'].find_binding_energies(binding_energies)

        for element, energy in binding_energies.items():
            binding_energies[element] = rmgpy.quantity.Energy(energy)

        self.binding_energies = binding_energies

    def correct_binding_energy(self, thermo, species, metal_to_scale_from=None, metal_to_scale_to=None):
        """
        Changes the provided thermo, by applying a linear scaling relation
        to correct the adsorption energy.

        :param thermo: starting thermo data
        :param species: the species (which is an adsorbate)
        :param metal_to_scale_from: the metal you want to scale from (string eg. 'Pt111' or None)
        :param metal_to_scale_to: the metal you want to scale to (string e.g 'Pt111' or None)
        :return: corrected thermo
        """

        if metal_to_scale_from == metal_to_scale_to:
            return thermo

        if metal_to_scale_to is None:
            metal_to_scale_to_binding_energies = self.binding_energies
        else:
            metal_to_scale_to_binding_energies = self.surface['metal'].find_binding_energies(metal_to_scale_to)

        if metal_to_scale_from is None:
            metal_to_scale_from_binding_energies = self.binding_energies
        else:
            metal_to_scale_from_binding_energies = self.surface['metal'].find_binding_energies(metal_to_scale_from)

        delta_atomic_adsorption_energy = {
            'C': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
            'H': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
            'O': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
            'N': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
            'F': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
        }

        for element, delta_energy in delta_atomic_adsorption_energy.items():
            try:
                delta_energy.value_si = metal_to_scale_to_binding_energies[element].value_si - metal_to_scale_from_binding_energies[element].value_si
            except KeyError:
                pass

        if all(-0.01 < v.value_si < 0.01 for v in delta_atomic_adsorption_energy.values()):
            return thermo

        molecule = species.molecule[0]
        # only want/need to do one resonance structure
        surface_sites = []
        for atom in molecule.atoms:
            if atom.is_surface_site():
                surface_sites.append(atom)
        normalized_bonds = {'C': 0., 'O': 0., 'N': 0., 'H': 0., 'F': 0., 'Li': 0.}
        max_bond_order = {'C': 4., 'O': 2., 'N': 3., 'H': 1., 'F': 1, 'Li': 1.}
        for site in surface_sites:
            numbonds = len(site.bonds)
            if numbonds == 0:
                # vanDerWaals
                pass
            else:
                assert len(site.bonds) == 1, "Each surface site can only be bonded to 1 atom"
                bonded_atom = list(site.bonds.keys())[0]
                bond = site.bonds[bonded_atom]
                if bond.is_single():
                    bond_order = 1.
                elif bond.is_double():
                    bond_order = 2.
                elif bond.is_triple():
                    bond_order = 3.
                elif bond.is_quadruple():
                    bond_order = 4.
                else:
                    raise NotImplementedError("Unsupported bond order {0} for binding energy "
                                              "correction.".format(bond.order))

                normalized_bonds[bonded_atom.symbol] += bond_order / max_bond_order[bonded_atom.symbol]

        if not isinstance(thermo, ThermoData):
            thermo = thermo.to_thermo_data()
            find_cp0_and_cpinf(species, thermo)

        # now edit the adsorptionThermo using LSR
        comments = []
        for element,bond in normalized_bonds.items():
            if bond:
                try:
                    change_in_binding_energy = delta_atomic_adsorption_energy[element].value_si * bond
                except KeyError:
                    continue
                thermo.H298.value_si += change_in_binding_energy
                comments.append(f'{bond:.2f}{element}')
        thermo.comment += " Binding energy corrected by LSR ({}) from {}".format('+'.join(comments), metal_to_scale_from)
        return thermo

    def get_thermo_data_for_surface_species(self, species):
        """
        Get the thermo data for an adsorbed species,
        by desorbing it, finding the thermo of the gas-phase
        species, then adding an adsorption correction that
        is found from the groups/adsorption tree.
        Does not apply linear scaling relationship.
        
        Returns a :class:`ThermoData` object, with no Cp0 or CpInf
        """

        if species.is_surface_site():
            raise DatabaseError("Can't estimate thermo of vacant site. Should be in library (and should be 0).")

        logging.debug("Trying to generate thermo for surface species using first of %d resonance isomer(s):",
                      len(species.molecule))
        molecule = species.molecule[0]
        # store any labeled atoms to reapply at the end
        labeled_atoms = molecule.get_all_labeled_atoms()
        molecule.clear_labeled_atoms()
        logging.debug("Before removing from surface:\n" + molecule.to_adjacency_list())
        # only want/need to do one resonance structure,
        # because will need to regenerate others in gas phase
        dummy_molecules = molecule.get_desorbed_molecules()
        for mol in dummy_molecules:
            mol.clear_labeled_atoms()
        if len(dummy_molecules) == 0:
            raise RuntimeError(f"Cannot get thermo for gas-phase molecule. No valid dummy molecules from original molecule:\n{molecule.to_adjacency_list()}")
        
        # if len(molecule) > 1, it will assume all resonance structures have already been generated when it tries to generate them, so evaluate each configuration separately and pick the lowest energy one by H298 value
        gas_phase_species_from_libraries = []
        gas_phase_species_estimates = []
        for dummy_molecule in dummy_molecules:
            dummy_species = Species()
            dummy_species.molecule = [dummy_molecule]
            dummy_species.generate_resonance_structures()
            dummy_species.thermo = self.get_thermo_data(dummy_species)
            if dummy_species.thermo.label:
                gas_phase_species_from_libraries.append(dummy_species)
            else:
                gas_phase_species_estimates.append(dummy_species)

        # define the comparison function to find the lowest energy
        def lowest_energy(species):
            if hasattr(species.thermo, 'H298'):
                return species.thermo.H298.value_si
            else:
                return species.thermo.get_enthalpy(298.0)

        if gas_phase_species_from_libraries:
            species = min(gas_phase_species_from_libraries, key=lowest_energy)
        else:
            species = min(gas_phase_species_estimates, key=lowest_energy)

        thermo = species.thermo
        thermo.comment = f"Gas phase thermo for {thermo.label or species.molecule[0].to_smiles()} from {thermo.comment}. Adsorption correction:"
        logging.debug("Using thermo from gas phase for species {}\n".format(species.label) + repr(thermo))

        if not isinstance(thermo, ThermoData):
            thermo = thermo.to_thermo_data()
            find_cp0_and_cpinf(species, thermo)

        # Get the adsorption energy
        # Create the ThermoData object
        adsorption_thermo = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
            H298=(0.0, "kJ/mol"),
            S298=(0.0, "J/(mol*K)"),
        )

        surface_sites = molecule.get_surface_sites()
        try:
            self._add_adsorption_correction(adsorption_thermo, self.groups[self.adsorption_groups], molecule, surface_sites)
        except (KeyError, DatabaseError):
            logging.error("Couldn't find in adsorption thermo database:")
            logging.error(molecule)
            logging.error(molecule.to_adjacency_list())
            raise

        # (group_additivity=True means it appends the comments)
        add_thermo_data(thermo, adsorption_thermo, group_additivity=True)

        if thermo.label:
            thermo.label += 'X' * len(surface_sites)

        find_cp0_and_cpinf(species, thermo)

        # if the molecule had labels, reapply them 
        for label,atom in labeled_atoms.items():
            if isinstance(atom,list):
                for a in atom:
                    a.label = label
            else:
                atom.label = label

        return thermo

    def _add_adsorption_correction(self, adsorption_thermo, adsorption_groups, molecule, surface_sites):
        """Add thermo adsorption correction(s) to estimate adsorbate thermo from gas phase.
        If the molecule is multidentate, multiple adsoption corrections may be applied if 
        there does not exist a multidentate adsorption group with the same number of sites.
        In this case, only the enthalpy correction (H298) will be used for subsequent groups
        to avoid over-correcting the entropy and heat capacity due to the loss of translational 
        and rotational degrees of freedom from the gas phase.

        Args:
            adsorption_thermo ([ThermoData]): the ThermoData object to add the correction(s)
            adsorption_groups ([database]): the groups database (adsorptionPt111)
            molecule ([Molecule]): the molecule to apply the thermo correction
            surface_sites ([list([Atom])]): a list of the surface site atoms in the molecule
        """

        number_of_surface_sites = len(surface_sites)

        matches = []
        for atom in surface_sites:
            labeled_atoms = {'*': atom}
            node = adsorption_groups.descend_tree(molecule, labeled_atoms)
            if node is None: 
                # no match, so try the next surface site
                continue
            while node is not None and node.data is None:
                node = node.parent
            if node is None:
                # no data, so try the next surface site
                continue
            data = node.data
            comment = node.label
            loop_count = 0
            while isinstance(data, str):
                loop_count += 1
                if loop_count > 100:
                    raise DatabaseError("Maximum iterations reached while following thermo group data pointers. A circular"
                                    f" reference may exist. Last node was {node.label} pointing to group called {data} in "
                                    f"database {adsorption_groups.label}")

                for entry in adsorption_groups.entries.values():
                    if entry.label == data:
                        data = entry.data
                        comment = entry.label
                        break
                else:
                    raise DatabaseError(f"Node {node.label} points to a non-existing group called {data} "
                                    f"in database {adsorption_groups.label}")
            data.comment = f'{adsorption_groups.label}({comment})'
            group_surface_sites = node.item.get_surface_sites()
            if len(group_surface_sites) == number_of_surface_sites:
                # all the surface sites are accounted for so add the adsorption group and return
                add_thermo_data(adsorption_thermo, data, group_additivity=True)
                return True
            else:
                # we have not found a full match yet, so append and keep looking
                matches.append((len(group_surface_sites),data))
        
        if len(matches) == 0:
            raise DatabaseError(f"Could not find an adsorption correction in {adsorption_groups.label} for {molecule}")
        matches.sort(key = lambda x: -x[0])
        # sort the matches by descending number of surface sites
        corrections_applied = 0
        # start a counter for the number of corrections applied
        for number_of_group_sites, data in matches:
            if number_of_surface_sites - number_of_group_sites < 0:
                # too many sites in this group, skip to the next one
                continue
            if not corrections_applied:
                # this is the first correction, so add H298, S298, and Cp
                add_thermo_data(adsorption_thermo, data, group_additivity=True)
            else:
                # We have already corrected S298 and Cp, so we only want to correct H298
                adsorption_thermo.H298.value_si += data.H298.value_si
                adsorption_thermo.comment += ' + H298({0})'.format(data.comment)
            corrections_applied += 1
            number_of_surface_sites -= number_of_group_sites
            if number_of_surface_sites <= 0:
                # we have corrected for all the sites
                if number_of_surface_sites < 0:
                    adsorption_thermo.comment += ' WARNING(Too many adsorption corrections were added to the thermo!'
                    adsorption_thermo.comment += 'The H298 is very likely understimated as a result!)'
                break
        
        if number_of_surface_sites > 0:
            adsorption_thermo.comment += ' WARNING({} surface sites were unaccounted for with adsorption corrections!'.format(number_of_surface_sites)
            adsorption_thermo.comment += 'The H298 is very likely overestimated as a result!)'

        return True

    def get_thermo_data_from_libraries(self, species, training_set=None):
        """
        Return the thermodynamic parameters for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before failing and returning None.
        `training_set` is used to identify if function is called during training set or not.
        During training set calculation we want to use gas phase thermo to not affect reverse
        rate calculation.
        
        Returns: ThermoData or None
        """
        import rmgpy.rmg.main
        thermo_data = None

        # chatelak 11/15/14: modification to introduce liquid phase thermo libraries
        library_list = deepcopy(self.library_order)  # copy the value to not affect initial object

        if rmgpy.rmg.main.solvent is not None:
            liq_libraries = []
            # Liquid phase simulation part:
            # This bloc "for": Identify liquid phase libraries and store them in liq_libraries
            for iterLib in library_list:
                if self.libraries[iterLib].solvent:
                    liq_libraries.append(iterLib)
            # Check in liq_libraries if thermo for species exists and return the first match.
            # Only if function not called by training_set
            if liq_libraries and training_set is None:
                for label in liq_libraries:
                    thermo_data = self.get_thermo_data_from_library(species, self.libraries[label])
                    if thermo_data is not None:
                        if len(thermo_data) != 3:
                            raise RuntimeError("thermo_data should be a tuple (thermo_data, library, entry), "
                                               "not {0}".format(thermo_data))
                        # Watch out comments changed: this is used later to apply solvation or not on
                        # species matching thermo. If required, Modify this carefully.
                        thermo_data[0].comment += 'Liquid thermo library: ' + label
                        return thermo_data
            # Remove liq_libraries from library_list if:
            #     called by training set (training_set=True) or if no thermo found in liqLibrairies
            # if no liquid library found this does nothing.
            for libIter in liq_libraries:
                library_list.remove(libIter)

        # Condition to execute this part: gas phase simulation or training set or liquid phase simulation with:
        #     noliquid libraries found or no matching species found in liquid libraries
        # If gas phase simulation library_list = self.library_order (just like before modifications) and they are
        # all gas phase, already checked by checkLibrairies function in database.load()
        # Check the libraries in order; return the first successful match
        for label in library_list:
            thermo_data = self.get_thermo_data_from_library(species, self.libraries[label])
            if thermo_data is not None:
                if len(thermo_data) != 3:
                    raise RuntimeError("thermo_data should be a tuple (thermo_data, library, entry), "
                                       "not {0}".format(thermo_data))
                if rmgpy.rmg.main.solvent is not None and training_set is None:
                    thermo_data[0].comment += 'Thermo library corrected for liquid phase: ' + label
                else:
                    thermo_data[0].comment += 'Thermo library: ' + label
                return thermo_data

        return None

    def get_all_thermo_data(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species`. The hits from the depository come
        first, then the libraries (in order), and then the group additivity
        estimate. This method is useful for a generic search job.
        
        Returns: a list of tuples (ThermoData, source, entry) 
        (Source is a library or depository, or None)
        """
        thermo_data_list = []
        # Data from depository comes first
        thermo_data_list.extend(self.get_thermo_data_from_depository(species))
        # Data from libraries comes second
        for label in self.library_order:
            data = self.get_thermo_data_from_library(species, self.libraries[label])
            if data:
                if len(data) != 3:
                    raise RuntimeError("data should be a tuple (thermo_data, library, entry), "
                                       "not {0}".format(data))
                data[0].comment += label
                thermo_data_list.append(data)

        # Last entry is always the estimate from group additivity
        # Make it a tuple
        # Distinguish surface species, as orignial get_thermo_data_from_groups does
        # not work for surface sites or surface species
        if species.is_surface_site():
            # Cannot estimate thermo of vacant site. Thermo stores in library
            pass
        elif species.contains_surface_site():
            try:
                # Estimate thermo of surface species based on modfied GA method
                data = (self.get_thermo_data_for_surface_species(species), None, None)
            except DatabaseError:
                # We don't have a GAV estimate, e.g. unsupported element
                pass
            else:
                thermo_data_list.append(data)
        else:
            try:
                data = (self.get_thermo_data_from_groups(species), None, None)
            except DatabaseError:
                # We don't have a GAV estimate, e.g. unsupported element
                pass
            else:
                # update group activity for symmetry
                data[0].S298.value_si -= constants.R * math.log(species.get_symmetry_number())
                thermo_data_list.append(data)
        # Return all of the resulting thermo parameters
        return thermo_data_list

    def get_thermo_data_from_depository(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species` from the depository. If no
        depository is loaded, a :class:`DatabaseError` is raised.
        
        Returns: a list of tuples (thermo_data, depository, entry) without any Cp0 or CpInf data.
        """
        items = []
        for entry in self.depository['stable'].entries.values():
            for molecule in species.molecule:
                if molecule.is_isomorphic(entry.item):
                    items.append((deepcopy(entry.data), self.depository['stable'], entry))
                    break
        for entry in self.depository['radical'].entries.values():
            for molecule in species.molecule:
                if molecule.is_isomorphic(entry.item):
                    items.append((deepcopy(entry.data), self.depository['radical'], entry))
                    break
        return items

    def get_thermo_data_from_library(self, species, library):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` from the specified thermodynamics
        `library`. If `library` is a string, the list of libraries is searched
        for a library with that name. If no match is found in that library,
        ``None`` is returned. If no corresponding library is found, a
        :class:`DatabaseError` is raised.
        
        Returns a tuple: (ThermoData, library, entry)  or None.
        """
        match = None
        for entry in library.entries.values():
            for molecule in species.molecule:
                if molecule.is_isomorphic(entry.item) and entry.data is not None:
                    thermo_data = deepcopy(entry.data)
                    thermo_data.label = entry.label
                    find_cp0_and_cpinf(species, thermo_data)
                    match = (thermo_data, library, entry)
                    break
            if match is not None:
                break
        if match is not None:
            # Move the matched molecule to the first position in the list
            species.molecule.remove(molecule)
            species.molecule.insert(0, molecule)
        return match

    def get_thermo_data_from_groups(self, species):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        The resonance isomer (molecule) with the lowest H298 is used, and as a side-effect
        the resonance isomers (items in `species.molecule` list) are sorted in ascending order.
        
        This does not account for symmetry. The method calling this sould correct for it.
        
        Returns: ThermoData
        """
        thermo = []
        for molecule in species.molecule:
            molecule.clear_labeled_atoms()
            molecule.update_atomtypes()
            tdata = self.estimate_thermo_via_group_additivity(molecule)
            thermo.append(tdata)

        indices = self.prioritize_thermo(species, thermo)

        species.molecule = [species.molecule[ind] for ind in indices]

        thermo_data = thermo[indices[0]]
        find_cp0_and_cpinf(species, thermo_data)
        return thermo_data

    def get_thermo_data_from_ml(self, species, ml_estimator, ml_settings):
        """
        Return the set of thermodynamic parameters corresponding to a
        given :class:`Species` object `species` by estimation using the
        ML estimator. Also compare the estimated uncertainties to the
        user-defined cutoffs. If any of the uncertainties are larger
        than their corresponding cutoffs, return None. Also check all
        other options in `ml_settings`.

        For HBI, the resonance isomer with the lowest H298 is used and
        the resonance isomers in species are sorted in ascending order.

        The entropy is not corrected for the symmetry of the molecule.
        This should be done later by the calling function.
        """
        molecule = species.molecule[0]

        min_heavy = ml_settings['min_heavy_atoms'] or 1
        max_heavy = ml_settings['max_heavy_atoms'] or np.inf
        min_carbon = ml_settings['min_carbon_atoms'] or 0
        max_carbon = ml_settings['max_carbon_atoms'] or np.inf
        min_oxygen = ml_settings['min_oxygen_atoms'] or 0
        max_oxygen = ml_settings['max_oxygen_atoms'] or np.inf
        min_nitrogen = ml_settings['min_nitrogen_atoms'] or 0
        max_nitrogen = ml_settings['max_nitrogen_atoms'] or np.inf

        element_count = molecule.get_element_count()
        n_heavy = sum(count for element, count in element_count.items() if element != 'H')

        if not (min_heavy <= n_heavy <= max_heavy):
            return None
        if not (min_carbon <= element_count.get('C', 0) <= max_carbon):
            return None
        if not (min_oxygen <= element_count.get('O', 0) <= max_oxygen):
            return None
        if not (min_nitrogen <= element_count.get('N', 0) <= max_nitrogen):
            return None
        if ml_settings['only_heterocyclics'] and not molecule.is_heterocyclic():
            return None
        if ml_settings['only_cyclics'] and not molecule.is_cyclic():
            return None
        min_cycle_overlap = ml_settings['min_cycle_overlap']
        if min_cycle_overlap > 0 and molecule.get_max_cycle_overlap() < min_cycle_overlap:
            return None

        if molecule.is_radical():
            thermo = [self.estimate_radical_thermo_via_hbi(mol, ml_estimator.get_thermo_data) for mol in species.molecule]
            H298 = np.array([tdata.H298.value_si for tdata in thermo])
            indices = H298.argsort()
            species.molecule = [species.molecule[ind] for ind in indices]
            thermo0 = thermo[indices[0]]
        else:
            thermo0 = ml_estimator.get_thermo_data_for_species(species)

        # The keys for this dictionary should match the keys in
        # `RMG.ml_uncertainty_cutoffs`. Use a temperature-weighted
        # average to estimate uncertainty for Cp.
        uncertainties = dict(
            H298=thermo0.H298.uncertainty_si,
            S298=thermo0.S298.uncertainty_si,
            Cp=np.average(thermo0.Cpdata.uncertainty_si, weights=thermo0.Tdata.value_si)
        )

        ml_uncertainty_cutoffs = ml_settings['uncertainty_cutoffs']
        if any(uncertainties[p] > ml_uncertainty_cutoffs[p].value_si for p in ml_uncertainty_cutoffs):
            return None
        else:
            return thermo0

    def prioritize_thermo(self, species, thermo_data_list):
        """
        Use some metrics to reorder a list of thermo data from best to worst.
        Return a list of indices with the desired order associated with the index of thermo from the data list.
        """
        if len(species.molecule) > 1:
            # Go further only if there is more than one isomer
            if species.molecule[0].is_cyclic():
                # Special treatment for cyclic compounds
                entries = []
                for i, thermo in enumerate(thermo_data_list):
                    ring_groups, polycyclic_groups = self.get_ring_groups_from_comments(thermo)

                    # Use rank as a metric for prioritizing thermo. 
                    # The smaller the rank, the better.
                    sum_rank = np.sum(
                        [3 if entry.rank is None else entry.rank for entry in ring_groups + polycyclic_groups])

                    # Also use number of aromatic rings as a metric, more aromatic rings is better
                    # Group values are generally fitted to the most aromatic resonance structure
                    num_arom_rings = species.molecule[i].count_aromatic_rings()

                    entries.append((i, thermo, sum_rank, -num_arom_rings))

                # Sort first by number of aromatic rings, then rank, then by enthalpy at 298 K
                entries.sort(key=lambda entry: (entry[3], entry[2], entry[1].get_enthalpy(298.)))
                indices = [entry[0] for entry in entries]
            else:
                # For noncyclics, default to original algorithm of ordering thermo based on the most stable enthalpy
                H298 = np.array([t.get_enthalpy(298.) for t in thermo_data_list])
                indices = H298.argsort().tolist()
            # Sort indices again by the Molecule.has_charge()
            indices.sort(key=lambda index: species.molecule[index].has_charge(), reverse=False)
            # Sort indices again by the Molecule.reactive flag
            indices.sort(key=lambda index: species.molecule[index].reactive, reverse=True)
        else:
            indices = [0]
        return indices

    def estimate_radical_thermo_via_hbi(self, molecule, stable_thermo_estimator):
        """
        Estimate the thermodynamics of a radical by saturating it,
        applying the provided stable_thermo_estimator method on the saturated species,
        then applying hydrogen bond increment corrections for the radical
        site(s) and correcting for the symmetry.
        
        No entropy is included in the returning term.
        This should be done later by the calling function.
        """
        if not molecule.is_radical():
            raise ValueError("Method only valid for radicals.")

        saturated_struct = molecule.copy(deep=True)
        added = saturated_struct.saturate_radicals()
        saturated_struct.props['saturated'] = True

        # Get thermo estimate for saturated form of structure
        if stable_thermo_estimator == self.get_thermo_data_from_libraries:
            # Get data from libraries
            saturated_spec = Species(molecule=[saturated_struct])
            thermo_data_sat = stable_thermo_estimator(saturated_spec)
            if thermo_data_sat:
                if len(thermo_data_sat) != 3:
                    raise RuntimeError("thermo_data should be a tuple (thermo_data, library, entry), "
                                       "not {0}".format(thermo_data_sat))
                thermo_data_sat = thermo_data_sat[0]
        else:
            thermo_data_sat = stable_thermo_estimator(saturated_struct)

        if thermo_data_sat is None:
            # We couldn't get thermo for the saturated species from libraries, ml, or qm
            # However, if we were trying group additivity, this could be a problem
            if stable_thermo_estimator == self.compute_group_additivity_thermo:
                logging.info("Thermo data of saturated {0} of molecule {1} is None.".format(saturated_struct, molecule))
            return None

        # Convert to ThermoData object if necessary in order to add and subtract from enthalpy and entropy values
        if not isinstance(thermo_data_sat, ThermoData):
            thermo_data_sat = thermo_data_sat.to_thermo_data()

        if not (stable_thermo_estimator == self.compute_group_additivity_thermo
                or isinstance(stable_thermo_estimator.__self__, MLEstimator)):
            # remove the symmetry contribution to the entropy of the saturated molecule
            # assumes that the thermo data comes from QMTP or from a thermolibrary
            thermo_data_sat.S298.value_si += constants.R * math.log(saturated_struct.get_symmetry_number())

        thermo_data = thermo_data_sat

        # For each radical site, get radical correction
        # Only one radical site should be considered at a time; all others
        # should be saturated with hydrogen atoms
        for atom in added:
            # Remove the added hydrogen atoms and bond and restore the radical
            for H, bond in added[atom]:
                saturated_struct.remove_bond(bond)
                saturated_struct.remove_atom(H)
                atom.increment_radical()
            saturated_struct.update()
            try:
                self._add_group_thermo_data(thermo_data, self.groups['radical'], saturated_struct, {'*': atom})
            except KeyError:
                logging.error("Couldn't find in radical thermo database:")
                logging.error(molecule)
                logging.error(molecule.to_adjacency_list())
                raise
            # Re-saturate
            for H, bond in added[atom]:
                saturated_struct.add_atom(H)
                saturated_struct.add_bond(bond)
                atom.decrement_radical()
            # Subtract the enthalpy of the added hydrogens
            for H, bond in added[atom]:
                thermo_data.H298.value_si -= 52.103 * 4184

        # Remove all of the interactions of the saturated structure. Then add the interactions of the radical.
        # Take C1=CC=C([O])C(O)=C1 as an example, we need to remove the interation of OH-OH, then add the interaction of Oj-OH.
        # For now, we only apply this part to cyclic structure because we only have radical interaction data for aromatic radical.
        if saturated_struct.is_cyclic():
            sssr = saturated_struct.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._remove_group_thermo_data(thermo_data, self.groups['longDistanceInteraction_cyclic'],
                                                       saturated_struct, {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass
            sssr = molecule.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._add_group_thermo_data(thermo_data, self.groups['longDistanceInteraction_cyclic'], molecule,
                                                    {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass

        # prevents the original thermo species name being used for the HBI corrected radical in species generation
        thermo_data.label = ''

        return thermo_data

    def estimate_thermo_via_group_additivity(self, molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Molecule` object ``molecule`` using the group additivity values
        method. If no group additivity values are loaded, a :class:`DatabaseError`
        is raised.

        The entropy is not corrected for the symmetry of the molecule,
        this should be done later by the calling function.
        """
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong.
        molecule.sort_atoms()

        if molecule.is_radical():
            thermo_data = self.estimate_radical_thermo_via_hbi(molecule, self.compute_group_additivity_thermo)
        else:
            thermo_data = self.compute_group_additivity_thermo(molecule)
        return thermo_data

    def compute_group_additivity_thermo(self, molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Molecule` object ``molecule`` using the group additivity values
        method. If no group additivity values are loaded, a :class:`DatabaseError`
        is raised.

        The entropy is not corrected for the symmetry of the molecule,
        this should be done later by the calling function.
        """

        assert not molecule.is_radical(), "This method is only for saturated non-radical species."
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong.
        molecule.sort_atoms()

        # Create the ThermoData object
        thermo_data = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
            H298=(0.0, "kJ/mol"),
            S298=(0.0, "J/(mol*K)"),
        )

        cyclic = molecule.is_cyclic()
        # Generate estimates of the thermodynamics parameters
        for atom in molecule.atoms:
            # Iterate over atoms and skip hydogens and halogens (since there are no groups centered on these atomtypes)
            if atom.is_non_hydrogen() and not atom.is_halogen():
                # Get initial thermo estimate from main group database
                data_added = False
                try:
                    data_added = self._add_group_thermo_data(thermo_data, self.groups['group'], molecule, {'*': atom})[1]
                except KeyError:
                    logging.error("Couldn't find in main thermo database:")
                    logging.error(molecule)
                    logging.error(molecule.to_adjacency_list())
                    raise
                if not data_added:
                    neighbors = ''.join(sorted([atom2.atomtype.label for atom2 in atom.edges.keys()
                                                if atom2.atomtype.label != 'H']))
                    neighbors += 'H' * len(['H' for atom2 in atom.edges.keys() if atom2.atomtype.label == 'H'])
                    if atom.atomtype.label == 'Cb':
                        neighbors = neighbors.replace('Cb', '')
                    group_str = f'{atom.atomtype.label}-{neighbors}'
                    if group_str not in ['O2d-CO', 'S2d-CS']:
                        if thermo_data.comment:
                            thermo_data.comment += f' + missing({group_str})'
                        else:
                            thermo_data.comment = f'Thermo group additivity estimation: missing({group_str})'
                # Correct for gauche and 1,5- interactions
                # Pair atom with its 1st and 2nd nonHydrogen neighbors, 
                # Then match the pair with the entries in the database longDistanceInteraction_noncyclic.py
                # Currently we only have gauche(1,4) and 1,5 interactions in that file. 
                # If you want to add more corrections for longer distance, please call get_nth_neighbor() method accordingly.
                # Potentially we could include other.py in this database, but it's a little confusing how to label atoms for the entries in other.py
                if not molecule.is_atom_in_cycle(atom):
                    for atom_2 in molecule.get_nth_neighbor([atom], [1, 2]):
                        if molecule.is_atom_in_cycle(atom_2) and not atom_2.is_bonded_to_halogen():
                            continue
                        # This is the correction for noncyclic structure. 
                        # If `atom_2` is bonded to a halogen, we apply noncyclic corrections regardless if `atom_2` is in a cycle or not.
                        # If `atom_2` is not bonded to a halogen, and `atom` or `atom_2` is in a cycle, do not apply this correction.
                        # Note that previously we do not do gauche for cyclic molecule, which is unreasonable for cyclic molecule with a long tail.
                        try:
                            self._add_group_thermo_data(thermo_data, self.groups['longDistanceInteraction_noncyclic'],
                                                        molecule, {'*1': atom, '*2': atom_2})
                        except KeyError:
                            pass
                try:
                    self._add_group_thermo_data(thermo_data, self.groups['other'], molecule, {'*': atom})
                except KeyError:
                    pass

        # Do long distance interaction correction for cyclic molecule. 
        # First get smallest set of smallest rings. 
        # Then for every single ring, generate the atom pairs by itertools.permutation.
        # Finally match the atom pair with the database.
        # WIPWIPWIPWIPWIPWIPWIP         #########################################         WIPWIPWIPWIPWIPWIPWIP
        # WIP: For now, in the database, if an entry describes the interaction between same groups, 
        # it will be halved because it will be counted twice here. 
        # Alternatively we could keep all the entries as their full values by using combinations instead of permutations here.
        # In that case, we need to add more lines to match from reverse side when we didn't hit the most specific level from the forward side.
        # PS: by saying 'forward side', I mean {'*1':atomPair[0], '*2':atomPair[1]}. So the following is the reverse side '{'*1':atomPair[1], '*2':atomPair[0]}'
        # In my opinion, it's cleaner to do it in the current way.
        # WIPWIPWIPWIPWIPWIPWIP         #########################################         WIPWIPWIPWIPWIPWIPWIP
        if cyclic:
            sssr = molecule.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._add_group_thermo_data(thermo_data, self.groups['longDistanceInteraction_cyclic'], molecule,
                                                    {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass

        # Do ring corrections separately because we only want to match
        # each ring one time

        if cyclic:
            monorings, polyrings = molecule.get_disparate_cycles()
            for ring in monorings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self._add_ring_correction_thermo_data_from_tree(thermo_data, self.groups['ring'], molecule, ring)
                except KeyError:
                    logging.error("Couldn't find a match in the monocyclic ring database even though "
                                  "monocyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.to_adjacency_list())
                    raise
            for polyring in polyrings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self._add_polycyclic_correction_thermo_data(thermo_data, molecule, polyring)
                except KeyError:
                    logging.error("Couldn't find a match in the polycyclic ring database even though "
                                  "polycyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.to_adjacency_list())
                    raise

        return thermo_data

    def _add_polycyclic_correction_thermo_data(self, thermo_data, molecule, polyring):
        """
        INPUT: `polyring` as a list of `Atom` forming a polycyclic ring
        OUTPUT: if the input `polyring` can be fully matched in polycyclic database, the correction
        will be directly added to `thermo_data`; otherwise, a heuristic approach will
        be applied.
        """
        # look up polycylic tree directly
        matched_group_thermodata, matched_group, is_partial_match = self._add_ring_correction_thermo_data_from_tree(
            None, self.groups['polycyclic'], molecule, polyring)

        # if partial match (non-H atoms number same between 
        # polycylic ring in molecule and match group)
        # otherwise, apply heuristic algorithm
        if not is_partial_match:
            if is_bicyclic(polyring) and matched_group.label in self.groups['polycyclic'].generic_nodes:
                # apply secondary decompostion formula
                # to get a estimated_group_thermodata
                estimated_bicyclic_thermodata = self.get_bicyclic_correction_thermo_data_from_heuristic(polyring)
                if not estimated_bicyclic_thermodata:
                    estimated_bicyclic_thermodata = matched_group_thermodata
                thermo_data = add_thermo_data(thermo_data, estimated_bicyclic_thermodata, group_additivity=True,
                                              verbose=True)
            else:
                # keep matched_group_thermodata as is
                thermo_data = add_thermo_data(thermo_data, matched_group_thermodata, group_additivity=True, verbose=True)
                # By setting verbose=True, we turn on the comments of polycyclic correction to pass the unittest.
                # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.
        else:
            self._add_poly_ring_correction_thermo_data_from_heuristic(thermo_data, polyring)

    def _add_poly_ring_correction_thermo_data_from_heuristic(self, thermo_data, polyring):
        """
        INPUT: `polyring` as a list of `Atom` forming a polycyclic ring, which can 
        only be partially matched.
        OUTPUT: `polyring` will be decomposed into a combination of 2-ring polycyclics
        and each one will be looked up from polycyclic database. The heuristic formula 
        is "polyring thermo correction = sum of correction of all 2-ring sub-polycyclics - 
        overlapped single-ring correction"; the calculated polyring thermo correction 
        will be finally added to input `thermo_data`.
        """

        # polyring decomposition
        bicyclics_merged_from_ring_pair, ring_occurrences_dict = bicyclic_decomposition_for_polyring(polyring)

        # loop over 2-ring cores
        for bicyclic in bicyclics_merged_from_ring_pair:
            matched_group_thermodata, matched_group, _ = self._add_ring_correction_thermo_data_from_tree(
                None, self.groups['polycyclic'], bicyclic, bicyclic.atoms)

            if matched_group.label in self.groups['polycyclic'].generic_nodes:
                # apply secondary decompostion formula
                # to get a estimated_group_thermodata
                estimated_bicyclic_thermodata = self.get_bicyclic_correction_thermo_data_from_heuristic(bicyclic.atoms)
                if not estimated_bicyclic_thermodata:
                    estimated_bicyclic_thermodata = matched_group_thermodata
                thermo_data = add_thermo_data(thermo_data, estimated_bicyclic_thermodata, group_additivity=True,
                                             verbose=True)
            else:
                # keep matched_group_thermodata as is
                thermo_data = add_thermo_data(thermo_data, matched_group_thermodata, group_additivity=True, verbose=True)

        # loop over 1-ring 
        for singleRingTuple, occurrence in ring_occurrences_dict.items():
            single_ring = list(singleRingTuple)

            if occurrence >= 2:
                submol, _ = convert_ring_to_sub_molecule(single_ring)

                if not is_aromatic_ring(submol):
                    aromatic_bonds = find_aromatic_bonds_from_sub_molecule(submol)
                    for aromaticBond in aromatic_bonds:
                        aromaticBond.set_order_num(1)

                    submol.saturate_unfilled_valence()
                    single_ring_thermodata = self._add_ring_correction_thermo_data_from_tree(
                        None, self.groups['ring'], submol, submol.atoms)[0]

                else:
                    submol.update()
                    single_ring_thermodata = self._add_ring_correction_thermo_data_from_tree(
                        None, self.groups['ring'], submol, submol.atoms)[0]
            for _ in range(occurrence - 1):
                thermo_data = remove_thermo_data(thermo_data, single_ring_thermodata, True, True)
                # By setting verbose=True, we turn on the comments of polycyclic correction to pass the unittest.
                # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.

    def get_bicyclic_correction_thermo_data_from_heuristic(self, bicyclic):

        # saturate if the bicyclic has unsaturated bonds
        # otherwise return None
        bicyclic_submol = convert_ring_to_sub_molecule(bicyclic)[0]
        saturated_bicyclic_submol, already_saturated = saturate_ring_bonds(bicyclic_submol)

        if already_saturated:
            return None
        # split bicyclic into two single ring submols
        single_ring_submols = split_bicyclic_into_single_rings(bicyclic_submol)

        # split saturated bicyclic into two single ring submols
        saturated_single_ring_submols = split_bicyclic_into_single_rings(saturated_bicyclic_submol)

        # apply formula: 
        # bicyclic correction ~= saturated bicyclic correction - 
        # saturated single ring corrections + single ring corrections

        estimated_bicyclic_thermo_data = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
            H298=(0.0, "kJ/mol"),
            S298=(0.0, "J/(mol*K)")
        )

        saturated_bicyclic_thermo_data = self._add_ring_correction_thermo_data_from_tree(
            None, self.groups['polycyclic'], saturated_bicyclic_submol, saturated_bicyclic_submol.atoms)[0]

        estimated_bicyclic_thermo_data = add_thermo_data(estimated_bicyclic_thermo_data,
                                                         saturated_bicyclic_thermo_data,
                                                         group_additivity=True)

        estimated_bicyclic_thermo_data.comment = "Estimated bicyclic component: " + \
                                                 saturated_bicyclic_thermo_data.comment

        for submol in saturated_single_ring_submols:

            if not is_aromatic_ring(submol):
                aromatic_bonds = find_aromatic_bonds_from_sub_molecule(submol)
                for aromatic_bond in aromatic_bonds:
                    aromatic_bond.set_order_num(1)

                submol.saturate_unfilled_valence()
                single_ring_thermo_data = self._add_ring_correction_thermo_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]

            else:
                submol.update()
                single_ring_thermo_data = self._add_ring_correction_thermo_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]
            estimated_bicyclic_thermo_data = remove_thermo_data(estimated_bicyclic_thermo_data,
                                                                single_ring_thermo_data,
                                                                group_additivity=True, verbose=True)

        for submol in single_ring_submols:

            if not is_aromatic_ring(submol):
                aromatic_bonds = find_aromatic_bonds_from_sub_molecule(submol)
                for aromatic_bond in aromatic_bonds:
                    aromatic_bond.set_order_num(1)

                submol.saturate_unfilled_valence()
                single_ring_thermo_data = self._add_ring_correction_thermo_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]

            else:
                submol.update()
                single_ring_thermo_data = self._add_ring_correction_thermo_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]

            estimated_bicyclic_thermo_data = add_thermo_data(estimated_bicyclic_thermo_data,
                                                             single_ring_thermo_data, group_additivity=True, verbose=True)

        return estimated_bicyclic_thermo_data

    def _add_ring_correction_thermo_data_from_tree(self, thermo_data, ring_database, molecule, ring):
        """
        Determine the ring correction group additivity thermodynamic data for the given
         `ring` in the `molecule`, and add it to the existing thermo data
        `thermo_data`.
        Also returns the matched ring group from the database from which the data originated.
        """
        matched_ring_entries = []
        # label each atom in the ring individually to try to match the group
        # for each ring, save only the ring that is matches the most specific leaf in the tree.
        for atom in ring:
            atoms = {'*': atom}
            entry = ring_database.descend_tree(molecule, atoms)
            matched_ring_entries.append(entry)

        if matched_ring_entries is []:
            raise KeyError('Node not found in database.')
        # Decide which group to keep
        is_partial_match = True
        complete_matched_groups = [entry for entry in matched_ring_entries
                                   if not is_ring_partial_matched(ring, entry.item)]

        if complete_matched_groups:
            is_partial_match = False
            matched_ring_entries = complete_matched_groups

        depth_list = [len(ring_database.ancestors(entry)) for entry in matched_ring_entries]
        most_specific_match_indices = [i for i, x in enumerate(depth_list) if x == max(depth_list)]

        most_specific_matched_entries = [matched_ring_entries[idx] for idx in most_specific_match_indices]
        if len(set(most_specific_matched_entries)) != 1:
            logging.debug('More than one type of node was found to be most specific for this ring.')
            logging.debug('This is either due to a database error in the ring or polycyclic groups, '
                          'or a partial match between the group and the full ring.')
            logging.debug(most_specific_matched_entries)

        # Condense the number of most specific groups down to one
        most_specific_matched_entry = matched_ring_entries[most_specific_match_indices[0]]

        node = most_specific_matched_entry

        if node is None:
            raise DatabaseError('Unable to determine thermo parameters for {0}: no data for {1} or '
                                'any of its ancestors.'.format(molecule, mostSpecificGroup))

        while node is not None and node.data is None:
            # do average of its children
            success, averaged_thermo_data = self._average_children_thermo(node, ring_database)
            if success:
                node.data = averaged_thermo_data
            else:
                node = node.parent

        data = node.data
        comment = node.label
        while isinstance(data, str) and data is not None:
            for entry in ring_database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    node = entry
                    break
        data.comment = '{0}({1})'.format(ring_database.label, comment)

        if thermo_data is None:
            return data, node, is_partial_match
        else:
            return add_thermo_data(thermo_data, data, group_additivity=True, verbose=True), node, is_partial_match
            # By setting verbose=True, we turn on the comments of ring correction to pass the unittest.
            # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.

    def _average_children_thermo(self, node, database):
        """
        Use children's thermo data to guess thermo data of parent `node` 
        that doesn't have thermo data built-in in tree yet. 
        For `node` has children that have thermo data, return success flag 
        `True` and the average thermo data.
        For `node` whose children that all have no thermo data, return flag
        `False` and None for the thermo data.
        """
        if not node.children:
            if node.data is None:
                return False, None
            else:
                return True, node.data
        else:
            children_thermo_data_list = []
            for child in node.children:
                if child.data is None:
                    success, child_thermo_data_average = self._average_children_thermo(child, database)
                    if success:
                        children_thermo_data_list.append(child_thermo_data_average)
                else:
                    data = child.data
                    while isinstance(data, str):
                        data = database.entries[data].data
                    children_thermo_data_list.append(data)
            if children_thermo_data_list:
                return True, average_thermo_data(children_thermo_data_list)
            else:
                return False, None

    def _add_group_thermo_data(self, thermo_data, database, molecule, atom):
        """
        Determine the group additivity thermodynamic data for the atom ``atom``
        in the structure ``molecule``, and add it to the existing thermo data
        ``thermo_data``.
        The parameter ``atom`` is a dictionary of label-atom pairs like {'*',atom}

        Returns:
            tuple: The combined ThermoData object and a bool flag indicating whether new data was added to it.
        """
        node0 = database.descend_tree(molecule, atom, None)
        if node0 is None:
            raise KeyError(f'Node not found for atom {atom} in molecule {molecule} in thermo database {database.label}.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node is not None and node.data is None:
            node = node.parent
        if node is None:
            raise DatabaseError(f'Unable to determine thermo parameters for atom {atom} in molecule {molecule}: '
                                f'no data for node {node0} or any of its ancestors in database {database.label}.')

        data = node.data
        comment = node.label
        loop_count = 0
        while isinstance(data, str):
            loop_count += 1
            if loop_count > 100:
                raise DatabaseError("Maximum iterations reached while following thermo group data pointers. A circular"
                                    f" reference may exist. Last node was {node.label} pointing to group called {data} in "
                                    f"database {database.label}")

            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
            else:
                raise DatabaseError(f"Node {node.label} points to a non-existing group called {data} "
                                    f"in database {database.label}")
        data.comment = f'{database.label}({comment})'

        # This code prints the hierarchy of the found node; useful for debugging
        # result = ''
        # while node is not None:
        #   result = ' -> ' + node.label + result
        #   node = node.parent
        # print result[4:]

        if thermo_data is None:
            return data, False
        else:
            if data.is_all_zeros():
                return thermo_data, False
            return add_thermo_data(thermo_data, data, group_additivity=True), True

    def _remove_group_thermo_data(self, thermo_data, database, molecule, atom):
        """
        Based on the _add_group_thermo_data method. Just replace the last line with 'return remove_thermo_data()'.
        Determine the group additivity thermodynamic data for the atom `atom` in the structure `structure`,
        and REMOVE it from the existing thermo data `thermo_data`.
        """
        node0 = database.descend_tree(molecule, atom, None)
        if node0 is None:
            raise KeyError(f'Node not found for atom {atom} in molecule {molecule} in thermo database {database.label}.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node.data is None and node is not None:
            node = node.parent
        if node is None:
            raise DatabaseError(f'Unable to determine thermo parameters for atom {atom} in molecule {molecule}: '
                                f'no data for node {node0} or any of its ancestors in database {database.label}.')

        data = node.data
        comment = node.label
        loop_count = 0
        while isinstance(data, str):
            loop_count += 1
            if loop_count > 100:
                raise DatabaseError("Maximum iterations reached while following thermo group data pointers. A circular"
                                    f" reference may exist. Last node was {node.label} pointing to group called {data} in "
                                    f"database {database.label}")
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
            else:
                raise DatabaseError(f"Node {node.label} points to a non-existing group called {data} "
                                    f"in database {database.label}")
        data.comment = f'{database.label}({comment})'

        # This code prints the hierarchy of the found node; useful for debugging
        # result = ''
        # while node is not None:
        #    result = ' -> ' + node.label + result
        #    node = node.parent
        # print result[4:]

        if thermo_data is None:
            return data
        else:
            return remove_thermo_data(thermo_data, data, True)

    def get_ring_groups_from_comments(self, thermo_data):
        """
        Takes a string of comments from group additivity estimation, and extracts the ring and polycyclic ring groups
        from them, returning them as lists.
        """
        tokens = thermo_data.comment.split()
        ring_groups = []
        polycyclic_groups = []
        regex = r"\((.*)\)"  # only hit outermost parentheses
        for token in tokens:
            if token.startswith('ring'):
                split_tokens = re.split(regex, token)
                assert len(split_tokens) == 3, 'token: {}'.format(token)
                group_label = split_tokens[1]
                ring_groups.append(self.groups['ring'].entries[group_label])
            if token.startswith('polycyclic'):
                split_tokens = re.split(regex, token)
                assert len(split_tokens) == 3, 'token: {}'.format(token)
                group_label = split_tokens[1]
                polycyclic_groups.append(self.groups['polycyclic'].entries[group_label])

        return ring_groups, polycyclic_groups

    def extract_source_from_comments(self, species):
        """
        `species`: A species object containing thermo data and thermo data comments
        
        Parses the verbose string of comments from the thermo data of the species object,
        and extracts the thermo sources.

        Returns a dictionary with keys of either 'Library', 'QM', and/or 'GAV'.
        Commonly, species thermo are estimated using only one of these sources.
        However, a radical can be estimated with more than one type of source, for 
        instance a saturated library value and a GAV HBI correction, or a QM saturated value
        and a GAV HBI correction.  
        
        source = {'Library': String_Name_of_Library_Used,
                  'QM': String_of_Method_Used,
                  'GAV': Dictionary_of_Groups_Used 
                  }
                  
        The Dictionary_of_Groups_Used looks like 
        {'groupType':[List of tuples containing (Entry, Weight)]
        """
        comment = species.thermo.comment
        tokens = comment.split()

        source = {}

        if comment.startswith('Thermo library'):
            # Store name of the library source, which is the 3rd token in the comments
            source['Library'] = tokens[2]

        elif comment.startswith('QM'):
            # Store the level of the calculation, which is the 2nd token in the comments
            source['QM'] = tokens[1]

        # Check for group additivity contributions to the thermo in this species            

        # The contribution of the groups can be either additive or substracting
        # after changes to the polycyclic algorithm

        comment = comment.replace(' + ', ' +')
        comment = comment.replace(' - ', ' -')
        tokens = comment.split()

        groups = {}
        group_types = list(self.groups.keys())

        regex = r"\((.*)\)"  # only hit outermost parentheses
        for token in tokens:
            weight = 1  # default contribution is additive
            if token.startswith('+'):
                token = token[1:]
            elif token.startswith('-'):
                weight = -1
                token = token[1:]
            for groupType in group_types:
                if token.startswith(groupType + '(') and token.endswith(')'):
                    split_tokens = re.split(regex, token)
                    group_label = split_tokens[1]
                    group_entry = self.groups[groupType].entries[group_label]
                    # Use dictionary to combine into weights when necessary
                    if groupType not in groups:
                        groups[groupType] = {group_entry: weight}
                    else:
                        if group_entry in groups[groupType]:
                            groups[groupType][group_entry] += weight
                        else:
                            groups[groupType][group_entry] = weight
                    break

        if groups:
            # Indicate that group additivity is used when it is either an HBI correction
            # onto a  thermo library or QM value, or if the entire molecule is estimated using group additivity
            # Save the groups into the source dictionary

            # Convert groups back into tuples 
            for groupType, groupDict in groups.items():
                groups[groupType] = list(groupDict.items())

            source['GAV'] = groups

        # Perform a sanity check that this molecule is estimated by at least one method
        if not list(source.keys()):
            raise ValueError('Species {0} thermo appears to not be estimated using any methods.'.format(species))

        return source


class ThermoCentralDatabaseInterface(object):
    """
    A class for interfacing with RMG online thermo central database.
    """

    def __init__(self, host, port, username, password, application):
        self.host = host
        self.port = port
        self.username = username
        self.password = password
        self.application = application
        self.client = self.connect()

    def connect(self):

        import pymongo

        remote_address = 'mongodb://{0}:{1}@{2}/thermoCentralDB'.format(self.username,
                                                                        self.password,
                                                                        self.host)
        client = pymongo.MongoClient(remote_address,
                                     self.port,
                                     serverSelectionTimeoutMS=2000)
        try:
            client.server_info()
            logging.info("\nConnection success to RMG Thermo Central Database!\n")
            return client

        except (pymongo.errors.ServerSelectionTimeoutError,
                pymongo.errors.OperationFailure):
            logging.info("\nConnection failure to RMG Thermo Central Database...")
            logging.info("This RMG job still can run but cannot utilize data from central database.\n")
            return None

    def satisfy_registration_requirements(self, species, thermo, thermodb):
        """
        Given a species, check if it's allowed to register in 
        central thermo database.

        Requirements for now: 
        cyclic, 
        its thermo is estimated by GAV and no exact match/use heuristics
        """
        if not species.molecule[0].is_cyclic():
            return False

        gav_keywords = 'Thermo group additivity estimation'
        if isinstance(thermo, ThermoData) and thermo.comment.startswith(gav_keywords):
            ring_groups, polycyclic_groups = thermodb.get_ring_groups_from_comments(thermo)

            # use GAV generic node to estimate thermo
            for group in ring_groups + polycyclic_groups:
                if group.label in thermodb.groups['ring'].generic_nodes + thermodb.groups['polycyclic'].generic_nodes:
                    return True

            # used some heuristic way to estimate thermo
            if ") - ring(" in thermo.comment:
                return True
            else:
                return False
        else:
            return False

    def register_in_central_thermo_db(self, species):

        # choose registration table
        db = getattr(self.client, 'thermoCentralDB')
        registration_table = getattr(db, 'registration_table')
        results_table = getattr(db, 'results_table')

        # prepare registration entry
        try:
            aug_inchi = species.get_augmented_inchi()

            # check if it's registered before or
            # already have available data in results_table
            registered_entries = list(registration_table.find({"aug_inchi": aug_inchi}))
            finished_entries = list(results_table.find({"aug_inchi": aug_inchi}))

            if len(registered_entries) + len(finished_entries) > 0:
                return

            smiles_input = species.molecule[0].to_smiles()
            status = 'pending'
            species_registration_entry = {'aug_inchi': aug_inchi,
                                          'SMILES_input': smiles_input,
                                          'radical_number': species.molecule[0].get_radical_count(),
                                          'status': status,
                                          'user': self.username,
                                          'application': self.application,
                                          'timestamp': time.time()
                                          }

            registration_table.insert(species_registration_entry)

        except ValueError:
            logging.info('Fail to generate inchi/smiles for species below:\n{0}'.format(species.to_adjacency_list()))


def find_cp0_and_cpinf(species, heat_capacity):
    """
    Calculate the Cp0 and CpInf values, and add them to the HeatCapacityModel object.
    """
    if heat_capacity.Cp0 is None:
        cp_0 = species.calculate_cp0()
        heat_capacity.Cp0 = (cp_0, "J/(mol*K)")
    if heat_capacity.CpInf is None:
        cp_inf = species.calculate_cpinf()
        heat_capacity.CpInf = (cp_inf, "J/(mol*K)")
