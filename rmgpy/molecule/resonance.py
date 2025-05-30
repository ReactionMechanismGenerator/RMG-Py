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
This module contains methods for generation of resonance structures of molecules.

The main function to generate all relevant resonance structures for a given
Molecule object is ``generate_resonance_structures``. It calls the necessary
functions for generating each type of resonance structure.

Currently supported resonance types:

- All species:
    - ``generate_allyl_delocalization_resonance_structures``: single radical shift with double or triple bond
    - ``generate_lone_pair_multiple_bond_resonance_structures``: lone pair shift with double or triple bond in a 3-atom system (between nonadjacent atoms)
    - ``generate_adj_lone_pair_radical_resonance_structures``: single radical shift with lone pair between adjacent atoms
    - ``generate_adj_lone_pair_multiple_bond_resonance_structures``: multiple bond shift with lone pair between adjacent atoms
    - ``generate_adj_lone_pair_radical_multiple_bond_resonance_structures``: multiple bond and radical shift with lone pair and radical  between adjacent atoms
    - ``generate_N5dc_radical_resonance_structures``: shift between radical and lone pair mediated by an N5dc atom
    - ``generate_aryne_resonance_structures``: shift between cumulene and alkyne forms of arynes, which are not considered aromatic in RMG
- Aromatic species only:
    - ``generate_optimal_aromatic_resonance_structures``: fully delocalized structure, where all aromatic rings have benzene bonds
    - ``generate_kekule_structure``: generate a single Kekule structure for an aromatic compound (single/double bond form)
    - ``generate_opposite_kekule_structure``: for monocyclic aromatic species, rotate the double bond assignment
    - ``generate_clar_structures``: generate all structures with the maximum number of pi-sextet assignments
- Multidentate adsorbates only
    - ``generate_adsorbate_shift_down_resonance_structures``: shift 2 electrons from a C=/#C bond to the X-C bond
    - ``generate_adsorbate_shift_up_resonance_structures``: shift 2 electrons from a X=/#C bond to a C-C bond
    - ``generate_adsorbate_conjugate_resonance_structures``: shift 2 electrons in a conjugate pi system for bridged X-C-C-C-X adsorbates
"""

import logging
from operator import attrgetter

import cython
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

import rmgpy.molecule.filtration as filtration
import rmgpy.molecule.pathfinder as pathfinder
from rmgpy.exceptions import AtomTypeError, KekulizationError, ResonanceError
from rmgpy.molecule.adjlist import Saturator
from rmgpy.molecule.fragment import CuttingLabel
from rmgpy.molecule.graph import Vertex
from rmgpy.molecule.kekulize import kekulize
from rmgpy.molecule.molecule import Atom, Bond, Molecule


def populate_resonance_algorithms(features=None):
    """
    Generate list of resonance structure algorithms relevant to the current molecule.

    Takes a dictionary of features generated by analyze_molecule().
    Returns a list of resonance algorithms.
    """
    cython.declare(method_list=list)
    method_list = []

    if features is None:
        method_list = [
            generate_allyl_delocalization_resonance_structures,
            generate_lone_pair_multiple_bond_resonance_structures,
            generate_adj_lone_pair_radical_resonance_structures,
            generate_adj_lone_pair_multiple_bond_resonance_structures,
            generate_adj_lone_pair_radical_multiple_bond_resonance_structures,
            generate_N5dc_radical_resonance_structures,
            generate_optimal_aromatic_resonance_structures,
            generate_aryne_resonance_structures,
            generate_kekule_structure,
            generate_clar_structures,
            generate_adsorbate_shift_down_resonance_structures,
            generate_adsorbate_shift_up_resonance_structures,
            generate_adsorbate_conjugate_resonance_structures
        ]
    else:
        # If the molecule is aromatic, then radical resonance has already been considered
        # If the molecule was falsely identified as aromatic, then is_aryl_radical will still accurately capture
        # cases where the radical is in an orbital that is orthogonal to the pi orbitals.
        if features['is_radical'] and not features['is_aromatic'] and not features['is_aryl_radical']:
            method_list.append(generate_allyl_delocalization_resonance_structures)
        if features['is_cyclic']:
            method_list.append(generate_aryne_resonance_structures)
        if features['hasNitrogenVal5']:
            method_list.append(generate_N5dc_radical_resonance_structures)
        if features['hasLonePairs']:
            method_list.append(generate_adj_lone_pair_radical_resonance_structures)
            method_list.append(generate_adj_lone_pair_multiple_bond_resonance_structures)
            method_list.append(generate_adj_lone_pair_radical_multiple_bond_resonance_structures)
            if not features['is_aromatic']:
                # The generate_lone_pair_multiple_bond_resonance_structures method may perturb the electronic
                # configuration of a conjugated aromatic system, causing a major slow-down (two orders of magnitude
                # slower in one observed case), and it doesn't necessarily result in new representative localized
                # structures. Here we forbid it for all structures bearing at least one aromatic ring as a "good enough"
                # solution. A more holistic approach would be to identify these cases in generate_resonance_structures,
                # and pass a list of forbidden atom ID's to find_lone_pair_multiple_bond_paths.
                method_list.append(generate_lone_pair_multiple_bond_resonance_structures)
        if features['is_multidentate']:
            method_list.append(generate_adsorbate_shift_down_resonance_structures)
            method_list.append(generate_adsorbate_shift_up_resonance_structures)
            method_list.append(generate_adsorbate_conjugate_resonance_structures)
    return method_list


def analyze_molecule(mol, save_order=False):
    """
    Identify key features of molecule important for resonance structure generation.
    `save_order` is used to maintain the atom order, when analyzing the molecule, defaults to False.

    Returns a dictionary of features.
    """
    cython.declare(features=dict)

    features = {'is_radical': mol.is_radical(),
                'is_cyclic': mol.is_cyclic(),
                'is_aromatic': False,
                'isPolycyclicAromatic': False,
                'is_aryl_radical': False,
                'hasNitrogenVal5': False,
                'hasLonePairs': False,
                'is_multidentate': mol.is_multidentate(),
                }

    if features['is_cyclic']:
        aromatic_rings = mol.get_aromatic_rings(save_order=save_order)[0]
        if len(aromatic_rings) > 0:
            features['is_aromatic'] = True
        if len(aromatic_rings) > 1:
            features['isPolycyclicAromatic'] = True
        if features['is_radical'] and features['is_aromatic']:
            features['is_aryl_radical'] = mol.is_aryl_radical(aromatic_rings)
    for atom in mol.vertices:
        if atom.is_nitrogen() and atom.lone_pairs == 0:
            features['hasNitrogenVal5'] = True
        if atom.lone_pairs > 0:
            features['hasLonePairs'] = True

    return features


def generate_resonance_structures(mol, clar_structures=True, keep_isomorphic=False,
                                  filter_structures=True, save_order=False):
    """
    Generate and return all of the resonance structures for the input molecule.

    Most of the complexity of this method goes into handling aromatic species, particularly to generate an accurate
    set of resonance structures that is consistent regardless of the input structure. The following considerations
    are made:

    1. False positives from RDKit aromaticity detection can occur if a molecule has exocyclic double bonds
    2. False negatives from RDKit aromaticity detection can occur if a radical is delocalized into an aromatic ring
    3. sp2 hybridized radicals in the plane of an aromatic ring do not participate in hyperconjugation
    4. Non-aromatic resonance structures of PAHs are not important resonance contributors (assumption)

    Aromatic species are broken into the following categories for resonance treatment:

    - Radical polycyclic aromatic species: Kekule structures are generated in order to generate adjacent resonance
      structures. The resulting structures are then used for Clar structure generation. After all three steps, any
      non-aromatic structures are removed, under the assumption that they are not important resonance contributors.
    - Radical monocyclic aromatic species: Kekule structures are generated along with adjacent resonance structures.
      All are kept regardless of aromaticity because the radical is more likely to delocalize into the ring.
    - Stable polycyclic aromatic species: Clar structures are generated
    - Stable monocyclic aromatic species: Kekule structures are generated
    """
    cython.declare(mol_list=list, new_mol_list=list, features=dict, method_list=list)

    # Check that mol is a valid structure in terms of atomTypes and net charge. Since SMILES with hypervalance
    # heteroatoms are not always read correctly, print a suggestion to input the structure using an adjList.
    try:
        mol.update(sort_atoms=not save_order)
    except AtomTypeError:
        logging.error("The following molecule has at least one atom with an undefined atomtype:\n{0}"
                      "\nIf this structure was entered in SMILES, try using the adjacencyList format for an unambiguous"
                      " definition.".format(mol.to_adjacency_list()))
        raise
    if mol.get_net_charge() != 0:
        # logging.info("Got the following structure:\nSMILES: {0}\nAdjacencyList:\n{1}\nNet charge: {2}\n\n"
        #                  "Currently RMG cannot process charged species correctly."
        #                  "\nIf this structure was entered in SMILES, try using the adjacencyList format for an"
        #                  " unambiguous definition. "
        #                  "Returning the input mol".format(mol.to_smiles(), mol.to_adjacency_list(), mol.get_net_charge()))
        return [mol]
    if not mol.reactive:
        raise ResonanceError('Can only generate resonance structures for reactive molecules! Got the following '
                             'unreactive structure:\n{0}Reactive = {1}'.format(mol.to_adjacency_list(), mol.reactive))

    mol_list = [mol]

    # Analyze molecule
    features = analyze_molecule(mol, save_order=save_order)

    # Use generate_optimal_aromatic_resonance_structures to check for false positives and negatives
    if features['is_aromatic'] or (features['is_cyclic'] and features['is_radical'] and not features['is_aryl_radical']):
        new_mol_list = generate_optimal_aromatic_resonance_structures(mol, features, save_order=save_order)
        if len(new_mol_list) == 0:
            # Encountered false positive, ie. the molecule is not actually aromatic
            features['is_aromatic'] = False
            features['isPolycyclicAromatic'] = False
        else:
            features['is_aromatic'] = True
            if len(new_mol_list[0].get_aromatic_rings(save_order=save_order)[0]) > 1:
                features['isPolycyclicAromatic'] = True
            for new_mol in new_mol_list:
                # Append to structure list if unique
                if not keep_isomorphic and mol.is_isomorphic(new_mol,
                                                             initial_map=None,
                                                             generate_initial_map=False,
                                                             save_order=save_order):
                    # Note: `initial_map` and `generate_initial_map` is using default values.
                    # They are required in compilation before assigning `save_order`.
                    continue
                elif keep_isomorphic and mol.is_identical(new_mol):
                    continue
                else:
                    mol_list.append(new_mol)

    # Special handling for aromatic species
    if features['is_aromatic']:
        if features['is_radical'] and not features['is_aryl_radical']:
            _generate_resonance_structures(mol_list, [generate_kekule_structure],
                                           keep_isomorphic=keep_isomorphic,
                                           save_order=save_order)
            _generate_resonance_structures(mol_list, [generate_allyl_delocalization_resonance_structures],
                                           keep_isomorphic=keep_isomorphic,
                                           save_order=save_order)
        if features['isPolycyclicAromatic'] and clar_structures:
            _generate_resonance_structures(mol_list, [generate_clar_structures],
                                           keep_isomorphic=keep_isomorphic,
                                           save_order=save_order)
        else:
            _generate_resonance_structures(mol_list, [generate_aromatic_resonance_structure],
                                           keep_isomorphic=keep_isomorphic,
                                           save_order=save_order)

    # Generate remaining resonance structures
    method_list = populate_resonance_algorithms(features)
    _generate_resonance_structures(mol_list, method_list, keep_isomorphic=keep_isomorphic,
                                   save_order=save_order)

    if filter_structures:
        return filtration.filter_structures(mol_list, features=features, save_order=save_order)

    return mol_list


def _generate_resonance_structures(mol_list, method_list, keep_isomorphic=False, copy=False,
                                   save_order=False):
    """
    Iteratively generate all resonance structures for a list of starting molecules using the specified methods.

    Args:
        mol_list             starting list of molecules
        method_list          list of resonance structure algorithms
        keep_isomorphic      if False, removes any structures that give is_isomorphic=True (default)
                            if True, only remove structures that give is_identical=True
        copy                if False, append new resonance structures to input list (default)
                            if True, make a new list with all of the resonance structures
    """
    cython.declare(index=cython.int, molecule=Graph, new_mol_list=list, new_mol=Graph, mol=Graph,
                   input_charge=cython.int, x=Vertex)

    if copy:
        # Make a copy of the list so we don't modify the input list
        mol_list = mol_list[:]

    min_octet_deviation = min(filtration.get_octet_deviation_list(mol_list))
    min_charge_span = min(filtration.get_charge_span_list(mol_list))

    # Iterate over resonance structures
    index = 0
    while index < len(mol_list):
        molecule = mol_list[index]
        new_mol_list = []

        # On-the-fly filtration: Extend methods only for molecule that don't deviate too much from the octet rule
        # (a +2 distance from the minimal deviation is used, octet deviations per species are in increments of 2)
        # Sometimes rearranging the structure requires an additional higher charge span structure, so allow
        # structures with a +1 higher charge span compared to the minimum, e.g., [O-]S#S[N+]#N
        # Filtration is always called.
        octet_deviation = filtration.get_octet_deviation(molecule)
        charge_span = molecule.get_charge_span()
        if octet_deviation <= min_octet_deviation + 2 and charge_span <= min_charge_span + 1:
            for method in method_list:
                new_mol_list.extend(method(molecule))
            if octet_deviation < min_octet_deviation:
                # update min_octet_deviation to make this criterion tighter
                min_octet_deviation = octet_deviation
            if charge_span < min_charge_span:
                # update min_charge_span to make this criterion tighter
                min_charge_span = charge_span

        for new_mol in new_mol_list:
            # Append to structure list if unique
            for mol in mol_list:
                if not keep_isomorphic and mol.is_isomorphic(new_mol,
                                                             initial_map=None,
                                                             generate_initial_map=False,
                                                             save_order=save_order):
                    # Note: `initial_map` and `generate_initial_map` is using default values.
                    # They are required in compilation before assigning `save_order`.
                    break
                elif keep_isomorphic and mol.is_identical(new_mol):
                    break
            else:
                mol_list.append(new_mol)

        # Move to the next resonance structure
        index += 1

    # check net charge
    input_charge = mol_list[0].get_net_charge()

    for mol in mol_list[1:]:
        if mol.get_net_charge() != input_charge:
            mol_list.remove(mol)
            logging.debug('Resonance generation created a molecule %s with a net charge of %d '
                          'which does not match the input mol charge of %d.\n'
                          'Removing it from resonance structures', mol.smiles, mol.get_net_charge(), input_charge)
        if mol.contains_surface_site():
            for x in [atom for atom in mol.atoms if atom.is_surface_site()]:
                if x.radical_electrons != 0:
                    mol_list.remove(mol)
                    logging.debug('Resonance generation created a molecule %s with %d radicals on %s.\n'
                                  'Removing it from resonance structures', mol.smiles, x.radical_electrons, x.symbol)
                elif x.lone_pairs != 0:
                    mol_list.remove(mol)
                    logging.debug('Resonance generation created a molecule %s with %d lone pairs on %s.\n'
                                  'Removing it from resonance structures', mol.smiles, x.lone_pairs, x.symbol)
                elif x.charge != 0:
                    mol_list.remove(mol)
                    logging.debug('Resonance generation created a molecule %s with a charge of %d on %s.\n'
                                  'Removing it from resonance structures', mol.smiles, x.charge, x.symbol)

    return mol_list


def generate_allyl_delocalization_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by one allyl radical shift.

    Biradicals on a single atom are not supported.
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex, atom3=Vertex, bond12=Edge, bond23=Edge)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    if mol.is_radical():  # Iterate over radicals in structure
        for atom in mol.vertices:
            paths = pathfinder.find_allyl_delocalization_paths(atom)
            for atom1, atom2, atom3, bond12, bond23 in paths:
                # Adjust to (potentially) new resonance structure
                atom1.decrement_radical()
                atom3.increment_radical()
                bond12.increment_order()
                bond23.decrement_order()
                # Make a copy of structure
                structure = mol.copy(deep=True)
                # Restore current structure
                atom1.increment_radical()
                atom3.decrement_radical()
                bond12.decrement_order()
                bond23.increment_order()
                try:
                    structure.update_atomtypes(log_species=False)
                except AtomTypeError:
                    pass  # Don't append resonance structure if it creates an undefined atomtype
                else:
                    structures.append(structure)
    return structures


def generate_lone_pair_multiple_bond_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by lone electron pair - multiple bond shifts in 3-atom systems.
    Examples: aniline (Nc1ccccc1), azide, [:NH2]C=[::O] <=> [NH2+]=C[:::O-]
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex, atom3=Vertex, bond12=Edge, bond23=Edge)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    for atom in mol.vertices:
        if atom.lone_pairs >= 1:
            paths = pathfinder.find_lone_pair_multiple_bond_paths(atom)
            for atom1, atom2, atom3, bond12, bond23 in paths:
                # Adjust to (potentially) new resonance structure
                atom1.decrement_lone_pairs()
                atom3.increment_lone_pairs()
                bond12.increment_order()
                bond23.decrement_order()
                atom1.update_charge()
                atom3.update_charge()
                # Make a copy of structure
                structure = mol.copy(deep=True)
                # Restore current structure
                atom1.increment_lone_pairs()
                atom3.decrement_lone_pairs()
                bond12.decrement_order()
                bond23.increment_order()
                atom1.update_charge()
                atom3.update_charge()
                try:
                    structure.update_atomtypes(log_species=False)
                except AtomTypeError:
                    pass  # Don't append resonance structure if it creates an undefined atomtype
                else:
                    structures.append(structure)
    return structures


def generate_adj_lone_pair_radical_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by lone electron pair - radical shifts between adjacent atoms.
    These resonance transformations do not involve changing bond orders.
    NO2 example: O=[:N]-[::O.] <=> O=[N.+]-[:::O-]
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    if mol.is_radical():  # Iterate over radicals in structure
        for atom in mol.vertices:
            paths = pathfinder.find_adj_lone_pair_radical_delocalization_paths(atom)
            for atom1, atom2 in paths:
                # Adjust to (potentially) new resonance structure
                atom1.decrement_radical()
                atom1.increment_lone_pairs()
                atom1.update_charge()
                atom2.increment_radical()
                atom2.decrement_lone_pairs()
                atom2.update_charge()
                # Make a copy of structure
                structure = mol.copy(deep=True)
                # Restore current structure
                atom1.increment_radical()
                atom1.decrement_lone_pairs()
                atom1.update_charge()
                atom2.decrement_radical()
                atom2.increment_lone_pairs()
                atom2.update_charge()
                try:
                    structure.update_atomtypes(log_species=False)
                except AtomTypeError:
                    pass  # Don't append resonance structure if it creates an undefined atomtype
                else:
                    structures.append(structure)
    return structures


def generate_adj_lone_pair_multiple_bond_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by lone electron pair - multiple bond shifts between adjacent atoms.
    Example: [:NH]=[CH2] <=> [::NH-]-[CH2+]
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    Here atom1 refers to the N/S/O atom, atom 2 refers to the any R!H (atom2's lone_pairs aren't affected)
    (In direction 1 atom1 <losses> a lone pair, in direction 2 atom1 <gains> a lone pair)
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph, direction=cython.int)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex, bond12=Edge)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    for atom in mol.vertices:
        paths = pathfinder.find_adj_lone_pair_multiple_bond_delocalization_paths(atom)
        for atom1, atom2, bond12, direction in paths:
            if direction == 1:  # The direction <increasing> the bond order
                atom1.decrement_lone_pairs()
                bond12.increment_order()
            elif direction == 2:  # The direction <decreasing> the bond order
                atom1.increment_lone_pairs()
                bond12.decrement_order()
            atom1.update_charge()
            atom2.update_charge()
            # Make a copy of structure
            structure = mol.copy(deep=True)
            # Restore current structure
            if direction == 1:  # The direction <increasing> the bond order
                atom1.increment_lone_pairs()
                bond12.decrement_order()
            elif direction == 2:  # The direction <decreasing> the bond order
                atom1.decrement_lone_pairs()
                bond12.increment_order()
            atom1.update_charge()
            atom2.update_charge()
            try:
                structure.update_atomtypes(log_species=False)
            except AtomTypeError:
                pass  # Don't append resonance structure if it creates an undefined atomtype
            else:
                if not (structure.get_net_charge() and structure.contains_surface_site()):
                    structures.append(structure)
    return structures


def generate_adj_lone_pair_radical_multiple_bond_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by lone electron pair - radical - multiple bond shifts between adjacent atoms.
    Example: [:N.]=[CH2] <=> [::N]-[.CH2]
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    Here atom1 refers to the N/S/O atom, atom 2 refers to the any R!H (atom2's lone_pairs aren't affected)
    This function is similar to generate_adj_lone_pair_multiple_bond_resonance_structures() except for dealing with the
    radical transformations.
    (In direction 1 atom1 <losses> a lone pair, gains a radical, and atom2 looses a radical.
    In direction 2 atom1 <gains> a lone pair, looses a radical, and atom2 gains a radical)
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph, direction=cython.int)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex, bond12=Edge)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    if mol.is_radical():  # Iterate over radicals in structure
        for atom in mol.vertices:
            paths = pathfinder.find_adj_lone_pair_radical_multiple_bond_delocalization_paths(atom)
            for atom1, atom2, bond12, direction in paths:
                if direction == 1:  # The direction <increasing> the bond order
                    atom1.decrement_lone_pairs()
                    bond12.increment_order()
                    atom1.increment_radical()
                    atom2.decrement_radical()
                elif direction == 2:  # The direction <decreasing> the bond order
                    atom1.increment_lone_pairs()
                    bond12.decrement_order()
                    atom1.decrement_radical()
                    atom2.increment_radical()
                atom1.update_charge()
                atom2.update_charge()
                # Make a copy of structure
                structure = mol.copy(deep=True)
                # Restore current structure
                if direction == 1:  # The direction <increasing> the bond order
                    atom1.increment_lone_pairs()
                    bond12.decrement_order()
                    atom1.decrement_radical()
                    atom2.increment_radical()
                elif direction == 2:  # The direction <decreasing> the bond order
                    atom1.decrement_lone_pairs()
                    bond12.increment_order()
                    atom1.increment_radical()
                    atom2.decrement_radical()
                atom1.update_charge()
                atom2.update_charge()
                try:
                    structure.update_atomtypes(log_species=False)
                except AtomTypeError:
                    pass  # Don't append resonance structure if it creates an undefined atomtype
                else:
                    structures.append(structure)
    return structures


def generate_N5dc_radical_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by radical and lone pair shifts mediated by an N5dc atom.
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Molecule)
    cython.declare(atom=Atom, atom2=Atom, atom3=Atom)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    for atom in mol.vertices:
        if atom.atomtype.label == 'N5dc' and atom.radical_electrons == 0 and len(atom.edges) == 3:
            paths = pathfinder.find_N5dc_radical_delocalization_paths(atom)
            for atom2, atom3 in paths:
                atom2.decrement_radical()
                atom2.increment_lone_pairs()
                atom3.decrement_lone_pairs()
                atom3.increment_radical()
                atom2.update_charge()
                atom3.update_charge()
                # Make a copy of structure
                structure = mol.copy(deep=True)
                # Restore current structure
                atom2.increment_radical()
                atom2.decrement_lone_pairs()
                atom3.increment_lone_pairs()
                atom3.decrement_radical()
                atom2.update_charge()
                atom3.update_charge()
                try:
                    structure.update_atomtypes(log_species=False)
                except AtomTypeError:
                    pass  # Don't append resonance structure if it creates an undefined atomtype
                else:
                    structures.append(structure)
    return structures


def generate_optimal_aromatic_resonance_structures(mol, features=None, save_order=False):
    """
    Generate the aromatic form of the molecule. For radicals, generates the form with the most aromatic rings.

    Returns result as a list.
    In most cases, only one structure will be returned.
    In certain cases where multiple forms have the same number of aromatic rings, multiple structures will be returned.
    If there's an error (eg. in RDKit) it just returns an empty list.
    """
    cython.declare(molecule=Graph, rings=list, aromaticBonds=list, kekuleList=list, maxNum=cython.int, mol_list=list,
                   new_mol_list=list, ring=list, bond=Bond, order=float, originalBonds=list, originalOrder=list,
                   i=cython.int, counter=cython.int)

    if features is None:
        features = analyze_molecule(mol, save_order=save_order)

    if not features['is_cyclic']:
        return []

    # Copy the molecule so we don't affect the original
    molecule = mol.copy(deep=True)

    # Attempt to rearrange electrons to obtain a structure with the most aromatic rings
    # Possible rearrangements include aryne resonance and allyl resonance
    res_list = [generate_aryne_resonance_structures]
    if features['is_radical'] and not features['is_aryl_radical']:
        res_list.append(generate_allyl_delocalization_resonance_structures)

    if molecule.is_aromatic():
        kekule_list = generate_kekule_structure(molecule)
    else:
        kekule_list = [molecule]

    _generate_resonance_structures(kekule_list, res_list, save_order=save_order)

    # Sort all of the generated structures by number of perceived aromatic rings
    mol_dict = {}
    for mol0 in kekule_list:
        aromatic_bonds = mol0.get_aromatic_rings(save_order=save_order)[1]
        num_aromatic = len(aromatic_bonds)
        mol_dict.setdefault(num_aromatic, []).append((mol0, aromatic_bonds))

    # List of potential number of aromatic rings, sorted from largest to smallest
    arom_options = sorted(mol_dict.keys(), reverse=True)

    new_mol_list = []
    for num in arom_options:
        mol_list = mol_dict[num]
        # Generate the aromatic resonance structure(s)
        for mol0, aromatic_bonds in mol_list:
            # Aromatize the molecule in place
            result = generate_aromatic_resonance_structure(mol0, aromatic_bonds, copy=False, save_order=save_order)
            if not result:
                # We failed to aromatize this molecule
                # This could be due to incorrect aromaticity perception by RDKit
                continue

            for mol1 in new_mol_list:
                if mol1.is_isomorphic(mol0, initial_map=None,
                                      generate_initial_map=False, save_order=save_order):
                    # Note: `initial_map` and `generagenerate_initial_map` is using default values.
                    # They are required in compilation before assigning `save_order`.
                    break
            else:
                new_mol_list.append(mol0)

        if new_mol_list:
            # We found the most aromatic resonance structures so there's no need to try smaller numbers
            break

    return new_mol_list


def generate_aromatic_resonance_structure(mol, aromatic_bonds=None, copy=True, save_order=False):
    """
    Generate the aromatic form of the molecule in place without considering other resonance.

    Args:
        mol: :class:`Molecule` object to modify
        aromatic_bonds (optional): list of previously identified aromatic bonds
        copy (optional): copy the molecule if ``True``, otherwise modify in place

    Returns:
        List of one molecule if successful, empty list otherwise
    """
    if copy:
        molecule = mol.copy(deep=True)
    else:
        molecule = mol

    if aromatic_bonds is None:
        aromatic_bonds = molecule.get_aromatic_rings(save_order=save_order)[1]
    if len(aromatic_bonds) == 0:
        return []

    # Save original bond orders in case this doesn't work out
    original_bonds = []
    for ring in aromatic_bonds:
        original_order = []
        for bond in ring:
            original_order.append(bond.order)
        original_bonds.append(original_order)
    # Change bond types to benzene bonds for all aromatic rings
    for ring in aromatic_bonds:
        for bond in ring:
            bond.order = 1.5

    try:
        molecule.update_atomtypes(log_species=False)
    except AtomTypeError:
        # If this didn't work the first time, then there might be a ring that is not actually aromatic
        # Reset our changes
        for ring, original_order in zip(aromatic_bonds, original_bonds):
            for bond, order in zip(ring, original_order):
                bond.order = order
        # Try to make each ring aromatic, one by one
        i = 0  # Track how many rings are aromatic
        counter = 0  # Track total number of attempts to avoid infinite loops
        while i < len(aromatic_bonds) and counter < 2 * len(aromatic_bonds):
            counter += 1
            original_order = []
            for bond in aromatic_bonds[i]:
                original_order.append(bond.order)
                bond.order = 1.5
            try:
                molecule.update_atomtypes(log_species=False)
            except AtomTypeError:
                # This ring could not be made aromatic, possibly because it depends on other rings
                # Undo changes
                for bond, order in zip(aromatic_bonds[i], original_order):
                    bond.order = order
                # Move it to the end of the list, and go on to the next ring
                aromatic_bonds.append(aromatic_bonds.pop(i))
                molecule.update_atomtypes(log_species=False)
                continue
            else:
                # We're done with this ring, so go on to the next ring
                i += 1
        # If we didn't end up making any of the rings aromatic, then this molecule is not actually aromatic
        if i == 0:
            # Move onto next molecule in the list
            return []

    return [molecule]


def generate_aryne_resonance_structures(mol):
    """
    Generate aryne resonance structures, including the cumulene and alkyne forms.

    For all 6-membered rings, check for the following bond patterns:

      - DDDSDS
      - STSDSD

    This does NOT cover all possible aryne resonance forms, only the simplest ones.
    Especially for polycyclic arynes, enumeration of all resonance forms is
    related to enumeration of all Kekule structures, which is very difficult.
    """
    cython.declare(rings=list, ring=list, new_mol_list=list, bond_list=list,
                   i=cython.int, j=cython.int, bond_orders=str, new_orders=str,
                   ind=cython.int, bond=Edge, new_mol=Graph)

    rings = mol.get_relevant_cycles()
    rings = [ring for ring in rings if len(ring) == 6]

    new_mol_list = []
    for ring in rings:
        # Get bond orders
        bond_list = mol.get_edges_in_cycle(ring)
        bond_orders = ''.join([bond.get_order_str() for bond in bond_list])
        new_orders = None
        # Check for expected bond patterns
        if bond_orders.count('T') == 1:
            # Reorder the list so the triple bond is first
            ind = bond_orders.index('T')
            bond_orders = bond_orders[ind:] + bond_orders[:ind]
            bond_list = bond_list[ind:] + bond_list[:ind]
            # Check for patterns
            if bond_orders == 'TSDSDS':
                new_orders = 'DDSDSD'
        elif bond_orders.count('D') == 4:
            # Search for DDD and reorder the list so that it comes first
            if 'DDD' in bond_orders:
                ind = bond_orders.index('DDD')
                bond_orders = bond_orders[ind:] + bond_orders[:ind]
                bond_list = bond_list[ind:] + bond_list[:ind]
            elif bond_orders.startswith('DD') and bond_orders.endswith('D'):
                bond_orders = bond_orders[-1:] + bond_orders[:-1]
                bond_list = bond_list[-1:] + bond_list[:-1]
            elif bond_orders.startswith('D') and bond_orders.endswith('DD'):
                bond_orders = bond_orders[-2:] + bond_orders[:-2]
                bond_list = bond_list[-2:] + bond_list[:-2]
            # Check for patterns
            if bond_orders == 'DDDSDS':
                new_orders = 'STSDSD'

        if new_orders is not None:
            # We matched one of our patterns, so we can now change the bonds
            for i, bond in enumerate(bond_list):
                bond.set_order_str(new_orders[i])
            # Make a copy of the molecule
            new_mol = mol.copy(deep=True)
            # Undo the changes to the current molecule
            for i, bond in enumerate(bond_list):
                bond.set_order_str(bond_orders[i])
            # Try to update atom types
            try:
                new_mol.update_atomtypes(log_species=False)
            except AtomTypeError:
                pass  # Don't append resonance structure if it creates an undefined atomtype
            else:
                new_mol_list.append(new_mol)

    return new_mol_list


def generate_kekule_structure(mol):
    """
    Generate a kekulized (single-double bond) form of the molecule.
    The specific arrangement of double bonds is non-deterministic, and depends on RDKit.

    Returns a single Kekule structure as an element of a list of length 1.
    If there's an error (eg. in RDKit) then it just returns an empty list.
    """
    cython.declare(atom=Vertex, molecule=Graph)

    for atom in mol.atoms:
        if isinstance(atom,CuttingLabel):
            continue
        if atom.atomtype.label == 'Cb' or atom.atomtype.label == 'Cbf':
            break
    else:
        return []

    molecule = mol.copy(deep=True)

    try:
        kekulize(molecule)
    except KekulizationError:
        return []

    return [molecule]


def generate_isomorphic_resonance_structures(mol, saturate_h=False):
    """
    Select the resonance isomer that is isomorphic to the parameter isomer, with the lowest unpaired
    electrons descriptor.

    We generate over all resonance isomers (non-isomorphic as well as isomorphic) and retain isomorphic
    isomers.

    If `saturate_h` is `True`, then saturate `mol` with hydrogens before generating the resonance structures,
    and remove the hydrogens before returning `isomorphic_isomers`. This is useful when resonance structures are
    generated for molecules in which all hydrogens were intentionally removed as in generating augInChI. Otherwise,
    RMG will probably get many of the lone_pairs and partial charges in a molecule wrong.

    WIP: do not generate aromatic resonance isomers.
    """

    cython.declare(isomorphic_isomers=list, isomers=list, index=int, max_val_e=int, order=float, num_h_to_add=int,
                   isomer=Molecule, newIsomer=Molecule, isom=Molecule, atom=Atom, a=Atom, b=Bond, newAtoms=list)

    if saturate_h:  # Add explicit hydrogen atoms to complete structure if desired
        Saturator.saturate(mol.vertices)

    isomorphic_isomers = [mol]  # resonance isomers that are isomorphic to the parameter isomer.

    isomers = [mol]

    # Iterate over resonance isomers
    index = 0
    while index < len(isomers):
        isomer = isomers[index]

        new_isomers = []
        for algo in populate_resonance_algorithms():
            new_isomers.extend(algo(isomer))

        for newIsomer in new_isomers:
            # Append to isomer list if unique
            for isom in isomers:
                if isom.copy(deep=True).is_isomorphic(newIsomer.copy(deep=True)):
                    isomorphic_isomers.append(newIsomer)
                    break
            else:
                isomers.append(newIsomer)

        # Move to next resonance isomer
        index += 1

    if saturate_h:  # remove hydrogens before returning isomorphic_isomers
        for isomer in isomorphic_isomers:
            isomer.delete_hydrogens()

    return isomorphic_isomers


def generate_clar_structures(mol, save_order=False):
    """
    Generate Clar structures for a given molecule.

    Returns a list of :class:`Molecule` objects corresponding to the Clar structures.
    """
    cython.declare(output=list, mol_list=list, new_mol=Graph, aromatic_rings=list, bonds=list, index=cython.int, bond=Edge, ring=list)

    if not mol.is_cyclic():
        return []

    # Atom IDs are necessary in order to maintain consistent matrices between iterations
    if not mol.atom_ids_valid():
        mol.assign_atom_ids()

    try:
        solutions = _clar_optimization(mol, save_order=save_order)
    except (RuntimeError, ValueError):  # either a crash during optimization or the result was an empty tuple
        # The optimization algorithm did not work on the first iteration
        return []

    mol_list = []

    for new_mol, aromatic_rings, bonds, solution in solutions:
        # The solution includes a part corresponding to rings, y, and a part corresponding to bonds, x, using
        # nomenclature from the paper. In y, 1 means the ring as a sextet, 0 means it does not.
        # In x, 1 corresponds to a double bond, 0 either means a single bond or the bond is part of a sextet.
        y = solution[:len(aromatic_rings)]
        x = solution[len(aromatic_rings):]

        # Apply results to molecule - double bond locations first
        for index, bond in enumerate(bonds):
            if x[index] == 0:
                bond.order = 1  # single
            elif x[index] == 1:
                bond.order = 2  # double
            else:
                raise ValueError(f'Unaccepted bond value {x[index]} obtained from optimization.')

        # Then apply locations of aromatic sextets by converting to benzene bonds
        for index, ring in enumerate(aromatic_rings):
            if y[index] == 1:
                for i, atom_1 in enumerate(ring):
                    for j, atom_2 in enumerate(ring):
                        if new_mol.has_bond(atom_1, atom_2):
                            new_mol.get_bond(atom_1, atom_2).order = 1.5

        try:
            new_mol.update_atomtypes()
        except AtomTypeError:
            pass
        finally:
            mol_list.append(new_mol)

    return mol_list

# helper functions for sorting
def _sum_atom_ids(atom_list):
    return sum(atom.id for atom in atom_list)

def _tuplize_bond(bond):
    return (bond.atom1.id, bond.atom2.id)


def _clar_optimization(mol, recursion_constraints=None, clar_number=-1, save_order=False):
    """
    Implements Mixed Integer Linear Programming for finding Clar structures.

    First finds the Clar number (and an arbitrary structure with that number), then recursively
    calls itself to enumerate more structures with that Clar number. No guarantees about
    which structures will be found, or how many (we stop solving once solution would require
    branching, which typically happens after at least a couple have already been found).

    Returns a list of valid Clar solutions in the form of a tuple, with the following entries:
        [0] Copy of mol
        [1] List of aromatic rings
        [2] List of bonds
        [3] Solution vector

    The solution vector is a binary integer list of length (# aromatic rings + # bonds) indicating
    if each ring is aromatic and if each bond is double.

    This implementation follows the original implementation very closely (Hansen, P.; Zheng, M. The
    Clar Number of a Benzenoid Hydrocarbon and Linear Programming. J. Math. Chem. 1994, 15 (1), 93–107.)
    with only slight modifications to prevent changing bond orders for exocyclic atoms.
    """
    cython.declare(
        molecule=Graph,
        aromatic_rings=list,
        exo_info=list,
        n_ring=cython.int,
        n_atom=cython.int,
        n_bond=cython.int,
        A=list,
        solutions=list,
    )

    # after we find all the Clar structures we will modify the molecule bond orders to be aromatic,
    # so we make explicit copies to avoid overwrites.
    molecule = mol.copy(deep=True)

    aromatic_rings = molecule.get_aromatic_rings(save_order=save_order)[0]

    # doesn't work for system without aromatic rings, just exit
    if len(aromatic_rings) == 0:
        return []

    # stability between multiple runs
    aromatic_rings.sort(key=_sum_atom_ids)

    # Cython doesn't allow mutable defaults, so we just set it here instead
    if recursion_constraints is None:
        recursion_constraints = []

    # Get set of atoms that are in rings
    atoms = set()
    for ring in aromatic_rings:
        atoms.update(ring)
    atoms = sorted(atoms, key=attrgetter('id'))

    # Get set of bonds involving the ring atoms, ignoring bonds to hydrogen
    bonds = set()
    for atom in atoms:
        bonds.update([atom.bonds[key] for key in atom.bonds.keys() if key.is_non_hydrogen()])
    bonds = sorted(bonds, key=_tuplize_bond)

    # identify exocyclic bonds, and save their order if exocyclic
    exo_info = []
    for bond in bonds:
        if bond.atom1 not in atoms or bond.atom2 not in atoms:
            # save order for exocyclic
            if bond.is_double():
                exo_info.append(1)
            else:
                exo_info.append(0)
        else:  # not exocyclic
            exo_info.append(None)

    # Dimensions
    n_ring = len(aromatic_rings)
    n_atom = len(atoms)
    n_bond = len(bonds)

    # The aromaticity assignment problem is formulated as an MILP problem
    # minimize:
    #       c @ x
    # such that
    #       b_l <= A @ x <= b_u
    #       l <= x <= u
    #       x are integers

    # Connectivity matrix which indicates which rings and bonds each atom is in
    # Part of equality constraint A_eq @ x = b_eq
    A = []
    for atom in atoms:
        in_ring = [1 if atom in ring else 0 for ring in aromatic_rings]
        in_bond = [1 if atom in [bond.atom1, bond.atom2] else 0 for bond in bonds]
        A.append(in_ring + in_bond)
    initial_constraint = [LinearConstraint(  # i.e. an equality constraint
        A=np.array(A, dtype=int), lb=np.ones(n_atom, dtype=int), ub=np.ones(n_atom, dtype=int)
    )]

    # on recursive calls we already know the Clar number, so we can additionally constrain the system
    # to only find structures with this number. Allows detecting when all formulae have been found and,
    # in theory, should solve faster
    if clar_number != -1:
        initial_constraint += [
            LinearConstraint(
                A=np.array([1] * n_ring + [0] * n_bond),
                ub=clar_number,
                lb=clar_number,
            ),
        ]

    # Objective vector for optimization: sextets have a weight of 1, double bonds have a weight of 0
    # we negate because the original problem is formulated as maximization, whereas SciPy only does min
    c = -np.array([1] * n_ring + [0] * n_bond, dtype=int)

    # variable bounds
    # - binary problem, so everything must be either zero or one
    # - rings are simply 0 or 1
    # - bonds are also 0 or 1, except exocyclic bonds which must remain unchanged
    bounds = Bounds(
        lb=np.array([0] * n_ring + [0 if b is None else b for b in exo_info], dtype=int),
        ub=np.array([1] * n_ring + [1 if b is None else b for b in exo_info], dtype=int),
    )

    result = milp(
        c=c,
        integrality=1,
        bounds=bounds,
        constraints=initial_constraint + recursion_constraints,
        options={'time_limit': 10},
    )

    if (status := result.status) != 0:
        if status == 2:  # infeasible
            raise RuntimeError("All valid Clar formulae have been enumerated!")
        else:
            raise RuntimeError(f"Optimization failed (Exit Code {status}) for an unexpected reason '{result.message}'")

    _clar_number, solution = -result.fun, result.x

    # optimization may have reached a bad local minimum - this case is rare
    if _clar_number == 0:
        return []

    # on later runs, non-integer solutions occur - branching might be able to find actual solutions,
    # but we just call it good enough here
    if any([x != 1 and x != 0 for x in solution]):
        raise RuntimeError("Optimization obtained a non-integer solution - no more formulae will be enumerated.")

    # first solution, so the result should be an upper limit
    if clar_number == -1:
        clar_number = _clar_number

    selected_sextets = list(solution[:n_ring])
    # restrict those rings which were selected for the current clar structure
    # from forming it again by requiring that their sum is 1 less than it used
    # to be
    recursion_constraints += [
        LinearConstraint(
            A=np.array(selected_sextets + [0] * n_bond),
            ub=clar_number - 1,
            lb=0,
        ),
    ]

    # Run optimization with additional constraints
    try:
        inner_solutions = _clar_optimization(mol, recursion_constraints=recursion_constraints, clar_number=clar_number, save_order=save_order)
    except RuntimeError as e:
        logging.debug(f"Clar Optimization stopped: {e}")
        inner_solutions = []

    return inner_solutions + [(molecule, aromatic_rings, bonds, solution)]


def generate_adsorbate_shift_down_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by the shift a pi bond between two C-C atoms to both X-C bonds.
    Example XCHXCH: [X]C=C[X] <=> [X]=CC=[X]
    (where '=' denotes a double bond)
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex, atom3=Vertex, atom4=Vertex, bond12=Edge, bond23=Edge, bond34=Edge)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    if mol.is_multidentate():
        for atom in mol.vertices:
            paths = pathfinder.find_adsorbate_delocalization_paths(atom)
            for atom1, atom2, atom3, atom4, bond12, bond23, bond34 in paths:
                if bond23.is_single():
                    continue
                else:
                    bond12.increment_order()
                    bond23.decrement_order()
                    bond34.increment_order()
                    structure = mol.copy(deep=True)
                    bond12.decrement_order()
                    bond23.increment_order()
                    bond34.decrement_order()
                    try:
                        structure.update_atomtypes(log_species=False)
                    except AtomTypeError:
                        pass
                    else:
                        structures.append(structure)
    return structures


def generate_adsorbate_shift_up_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by the shift of two electrons from X-C bonds to increase the bond
    order between two C-C atoms by 1.
    Example XCHXCH: [X]=CC=[X] <=> [X]C=C[X]
    (where '=' denotes a double bond, '#' denotes a triple bond)
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex, atom3=Vertex, atom4=Vertex, bond12=Edge, bond23=Edge, bond34=Edge)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    if mol.is_multidentate():
        for atom in mol.vertices:
            paths = pathfinder.find_adsorbate_delocalization_paths(atom)
            for atom1, atom2, atom3, atom4, bond12, bond23, bond34 in paths:
                if ((bond12.is_double_or_triple() and bond23.is_single() and bond34.is_double_or_triple()) or
                    (bond12.is_double() and bond23.is_double() and bond34.is_double())):
                        bond12.decrement_order()
                        bond23.increment_order()
                        bond34.decrement_order()
                        structure = mol.copy(deep=True)
                        bond12.increment_order()
                        bond23.decrement_order()
                        bond34.increment_order()
                        try:
                            structure.update_atomtypes(log_species=False)
                        except AtomTypeError:
                            pass
                        else:
                            structures.append(structure)
    return structures


def generate_adsorbate_conjugate_resonance_structures(mol):
    """
    Generate all of the resonance structures formed by the shift of two
    electrons in a conjugated pi bond system of a bidentate adsorbate
    with a bridging atom in between.

    Example XCHCHXC: [X]#CC=C[X] <=> [X]=C=CC=[X]
    (where '#' denotes a triple bond, '=' denotes a double bond)
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Graph)
    cython.declare(atom=Vertex, atom1=Vertex, atom2=Vertex, atom3=Vertex, atom4=Vertex, atom5=Vertex, bond12=Edge, bond23=Edge, bond34=Edge, bond45=Edge)
    cython.declare(v1=Vertex, v2=Vertex)

    structures = []
    if mol.is_multidentate():
        for atom in mol.vertices:
            paths = pathfinder.find_adsorbate_conjugate_delocalization_paths(atom)
            for atom1, atom2, atom3, atom4, atom45, bond12, bond23, bond34, bond45 in paths:
                    if (bond12.is_double_or_triple() and
                        (bond23.is_single() or bond23.is_double()) and
                        bond34.is_double_or_triple() and
                        (bond45.is_single() or bond45.is_double())):
                                bond12.decrement_order()
                                bond23.increment_order()
                                bond34.decrement_order()
                                bond45.increment_order()
                                structure = mol.copy(deep=True)
                                bond12.increment_order()
                                bond23.decrement_order()
                                bond34.increment_order()
                                bond45.decrement_order()
                                try:
                                    structure.update_atomtypes(log_species=False)
                                except AtomTypeError:
                                    pass
                                else:
                                    structures.append(structure)
    return structures
