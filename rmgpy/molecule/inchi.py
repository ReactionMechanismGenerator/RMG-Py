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

import itertools
import re

import cython
from rdkit import Chem

import rmgpy.molecule.element as elements
import rmgpy.molecule.pathfinder as pathfinder
import rmgpy.molecule.resonance as resonance
from rmgpy.exceptions import InchiException
from rmgpy.molecule.adjlist import ConsistencyChecker
from rmgpy.molecule.converter import toRDKitMol
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.util import agglomerate, partition, generate_combo, swap

# search for (*) PARENTHESES
PARENTHESES = re.compile(r'\((.[^\(\)]*)\)')

INCHI_PREFIX = 'InChI=1'

"""
The prefix with the information on the distribution of unpaired electrons across the atoms.

For example, a triplet biradical with unpaired electrons on atom 1 and atom 3 
will have the following unpaired electron layer in the augmented InChI:

InChI=1/.../u1,3

The indices refer to the 1-based indices in the InChI string (NOT the 0-based
    indices of the Molecule container!)

"""
U_LAYER_PREFIX = '/u'

"""The separator that separates the indices of the atoms that bear unpaired electrons."""
U_LAYER_SEPARATOR = ','

"""
The prefix with the information on the distribution of the atoms 
with an unexpected number of lone pairs.

For example, a singlet methylene with a lone pair on atom 1
will have the following lone pair layer in the augmented InChI:

InChI=1/.../lp1

The indices refer to the 1-based indices in the InChI string (NOT the 0-based
    indices of the Molecule container!)

"""
P_LAYER_PREFIX = '/lp'

"""The separator that separates the indices of the atoms that bear unpaired electrons."""
P_LAYER_SEPARATOR = ','

ulayer_pattern = re.compile(U_LAYER_PREFIX + r'(.*)')
player_pattern = re.compile(P_LAYER_PREFIX + r'(.*)')


def decompose_aug_inchi(string):
    """
    Converts an augmented inchi into 
    - an inchi, 
    - indices array for the atoms bearing unpaired electrons.
    - indices array for the atoms bearing (unexpected) lone pairs.

    Atoms that unexpectedly bear zero lone pairs will be mentioned as follows:
    "x(0)", with x the index of the atom, and (0) denoting that the atom does not 
    bear any lone pairs.

    The "x(0)" will be parsed into a tuple (x, 0).
    """
    cython.declare(
        inchi=str,
        u_indices=list,
        p_indices=list,
    )

    if U_LAYER_PREFIX in string:
        inchi = string.split(U_LAYER_PREFIX)[0]
    elif P_LAYER_PREFIX in string:
        inchi = string.split(P_LAYER_PREFIX)[0]
    else:
        inchi = string

    u_indices, p_indices = [], []
    matches = re.findall(ulayer_pattern, string)
    if matches:
        u_indices = map(int, matches.pop().split('/')[0].split(U_LAYER_SEPARATOR))

    matches = re.findall(player_pattern, string)
    if matches:
        dummy = matches.pop().split('/')[0].split(P_LAYER_SEPARATOR)
        for index in dummy:
            if '(0)' in str(index):
                index = int(str(index).split('(0)')[0])
                p_indices.append((index, 0))
            else:
                p_indices.append(int(index))

    return inchi, u_indices, p_indices


def remove_inchi_prefix(string):
    """
    Splits off the 'InChI=1S' or 'InChI=1' layer of an InChI
    and returns the last part.
    """

    if not INCHI_PREFIX in string:
        raise InchiException('Not a valid InChI: {}'.format(string))

    return re.split(r"(InChI=1+)(S*)/", string)[-1]


def compose_aug_inchi(inchi, ulayer=None, player=None):
    """
    Composes an augmented InChI by concatenating the different pieces
    as follows:

    InChI=1S/XXXX.../c.../h.../ux,x,/...
    """
    cython.declare(
        temp=str,
    )

    aug_inchi = INCHI_PREFIX + '/' if not INCHI_PREFIX in inchi else ''
    aug_inchi += inchi

    for layer in filter(None, [ulayer, player]):
        aug_inchi += layer

    return aug_inchi


def compose_aug_inchi_key(inchi_key, ulayer=None, player=None):
    """
    Composes an augmented InChI Key by concatenating the different pieces
    as follows:

    XXXXXXXXXXXXXX-XXXXXXXXXX-ux,x,xxx

    Uses hyphens rather than forward slashes to avoid messing up file paths.
    """

    aug_inchi_key = inchi_key

    for layer in filter(None, [ulayer, player]):
        aug_inchi_key += '-' + layer[1:]  # cut off the '/'

    return aug_inchi_key


##########################################################
# Methods for parsing layers in InChI and auxiliary info #
##########################################################

def _parse_H_layer(inchi):
    """
    Converts the Mobile H layer of an inchi string into a 
    list of atom indices couples that carry a mobile hydrogen.

    Example:
    The hydrogen of the hydroxyl group can migrate to the carbonyl 
    oxygen.

    O=C-O
    InChI=1S/CH2O2/c2-1-3/h1H,(H,2,3)

    The atoms bearing a mobile hydrogen will be found in
    the hydrogen layer, within the parentheses.

    An empty list will be returned when there are no mobile hydrogens.

    """

    cython.declare(
        pieces=list,
        h_layer=str,
        piece=str,
        couples=list,
        match=str,
        mobile_h_atoms=list,
    )

    pieces = inchi.split('/')
    h_layer = None
    for piece in pieces:
        if piece.startswith('h'):
            h_layer = piece
            break
    else:
        raise Exception('Could not find the hydrogen layer in the inchi: {}'.format(inchi))

    couples = []
    for match in re.findall(PARENTHESES, h_layer):
        mobile_h_atoms = map(int, match[2:].split(','))
        couples.append(mobile_h_atoms)

    return couples


def _parse_E_layer(auxinfo):
    """
    Converts the layer with equivalence information (E-layer) 
    on atoms into a list of lists of equivalent atom indices.

    Example:
    Auxiliary info of InChI=1S/C8H14/c1-5-7(3)8(4)6-2/h5-8H,1-2H2,3-4H3:
    AuxInfo=1/0/N:1,8,4,6,2,7,3,5/E:(1,2)(3,4)(5,6)(7,8)/rA:8C.2C.2CCCCCC/rB:s1;s2;s3;s3;s5;s5;d7;/rC:;;;;;;;;
    E-layer: 

    /E:(1,2)(3,4)(5,6)(7,8)/

    denotes that atoms (1,2), (3,4), (5,6), (7,8) are equivalent and cannot be distinguished based on the 
    implemented canonicalization algorithm.

    Returned object:

    [[1,2],[3,4],[5,6],[7,8]]

    Returns an empty list of the E-layer could not be found.

    """

    cython.declare(
        pieces=list,
        e_layer=str,
        piece=str,
        equivalent_atoms=list,
        atomtuple=str,
        indices=list,
    )

    pieces = auxinfo.split('/')
    e_layer = None
    for piece in pieces:
        if piece.startswith('E'):
            e_layer = piece[2:]  # cut off /E:
            break
    else:
        return []

    equivalent_atoms = []
    for atomtuple in re.findall(PARENTHESES, e_layer):
        indices = list(map(int, atomtuple.split(',')))
        equivalent_atoms.append(indices)

    return equivalent_atoms


def _parse_N_layer(auxinfo):
    """
    Parses the layer with atom ordering information (N-layer)
    and returns a list of atom indices that reflect how the atoms of the original
    molecule should be ordered according to the InChI algorithm.


    Example:
    Auxiliary info of SMILES OCCC (InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3):
    AuxInfo=1/0/N:4,3,2,1/rA:4OCCC/rB:s1;s2;s3;/rC:;;;;

    N-layer:
    /N:4,3,2,1

    The original number of an atom with identification number n is given as the
    n-th member of this list for a component; the lists are separated with ";".

    Raises an exception when the N-layer could not be found.
    """

    cython.declare(
        pieces=list,
        atom_numbers=str,
        piece=str,
        indices=list,
    )

    pieces = auxinfo.split('/')
    atom_numbers = None
    for piece in pieces:
        if piece.startswith('N'):
            atom_numbers = piece[2:]  # cut off N:
            break
    else:
        raise Exception('Could not find the N-layer in the auxiliary info: {}'.format(auxinfo))

    indices = map(int, atom_numbers.split(','))

    return indices


###########################################
# Methods for generating augmented layers #
###########################################

def _has_unexpected_lone_pairs(mol):
    """
    Iterates over the atoms of the Molecule and returns whether
    at least one atom bears an unexpected number of lone pairs.

    E.g.
    carbon with > 0 lone pairs
    nitrogen with > 1 lone pairs
    oxygen with > 2 lone pairs

    The expected number of lone pairs of an element is equal to
    """

    for at in mol.atoms:
        try:
            exp = elements.PeriodicSystem.lone_pairs[at.symbol]
        except KeyError:
            raise Exception("Unrecognized element: {}".format(at.symbol))
        else:
            if at.lonePairs != exp: return True

    return False


def _get_unpaired_electrons(mol):
    """
    Returns a sorted list of the indices of the atoms that bear one or more
    unpaired electrons.
    """

    cython.declare(
        locations=list,
        index=int,
        at=Atom,
    )
    locations = []
    for index, at in enumerate(mol.atoms):
        if at.radicalElectrons >= 1:
            locations.append(index)

    return sorted(locations)


def _generate_minimum_resonance_isomer(mol):
    """
    Select the resonance isomer that is isomorphic to the parameter isomer, with the lowest unpaired
    electrons descriptor.

    First, we generate all isomorphic resonance isomers.
    Next, we return the candidate with the lowest unpaired electrons metric.

    The metric is a sorted list with indices of the atoms that bear an unpaired electron
    """

    cython.declare(
        candidates=list,
        sel=Molecule,
        cand=Molecule,
        metric_sel=list,
        metric_cand=list,
    )

    candidates = resonance.generate_isomorphic_resonance_structures(mol, saturate_h=True)

    sel = candidates[0]
    metric_sel = _get_unpaired_electrons(sel)
    for cand in candidates[1:]:
        metric_cand = _get_unpaired_electrons(cand)
        if metric_cand < metric_sel:
            sel = cand
            metric_sel = metric_cand

    return sel


def _compute_agglomerate_distance(agglomerates, mol):
    """
    Iterates over a list of lists containing atom indices.
    For each list the distances between the atoms is computed.
    A list of distances is returned.

    """

    cython.declare(
        distances=list,
        agglomerate=list,
        dist=dict,
    )

    distances = []
    for agglomerate in agglomerates:
        dist = pathfinder.compute_atom_distance(agglomerate, mol)
        distances.append(dist)

    return distances


def _is_valid_combo(combo, mol, distances):
    """
    Check if the combination of atom indices refers to
    atoms that are adjacent in the molecule.
    """
    cython.declare(
        agglomerates=list,
        new_distances=list,
        orig_dist=dict,
        new_dist=dict,
    )

    # compute shortest path between atoms
    agglomerates = agglomerate(combo)
    new_distances = _compute_agglomerate_distance(agglomerates, mol)

    # combo is valid if the distance is equal to the parameter distance

    if len(distances) != len(new_distances): return False

    for orig_dist, new_dist in zip(distances, new_distances):
        # only compare the values of the dictionaries:
        if sorted(orig_dist.values()) != sorted(new_dist.values()):
            return False

    return True


def _find_lowest_u_layer(mol, u_layer, equivalent_atoms):
    """
    Searches for the "minimum" combination of indices of atoms that bear unpaired electrons.

    It does so by using the information on equivalent atoms to permute equivalent atoms to
    obtain a combination of atoms that is the (numerically) lowest possible combination.

    Each possible combination is valid if and only if the distances between the atoms of the
    combination is identical to the distances between the original combination.

    First, the algorithm partitions equivalent atoms that bear an unpaired electron.
    Next, the combinations are generated, and for each combination it is verified whether
    it pertains to a "valid" combination.

    Returns a list of indices corresponding to the lowest combination of atom indices bearing
    unpaired electrons.
    """

    cython.declare(
        new_u_layer=list,
        grouped_electrons=list,
        corresponding_E_layers=list,
        group=list,
        e_layer=list,
        combos=list,
        orig_agglomerates=list,
        orig_distances=list,
        selected_group=list,
        combo=list,
    )
    if not equivalent_atoms:
        return u_layer

    new_u_layer = []

    grouped_electrons, corresponding_E_layers = partition(u_layer, equivalent_atoms)

    # don't process atoms that do not belong to an equivalence layer
    for group, e_layer in zip(grouped_electrons[:], corresponding_E_layers[:]):
        if not e_layer:
            new_u_layer.extend(group)
            grouped_electrons.remove(group)
            corresponding_E_layers.remove(e_layer)

    combos = generate_combo(grouped_electrons, corresponding_E_layers)
    # compute original distance:
    orig_agglomerates = agglomerate(grouped_electrons)
    orig_distances = _compute_agglomerate_distance(orig_agglomerates, mol)

    # deflate the list of lists to be able to numerically compare them
    selected_group = sorted(itertools.chain.from_iterable(grouped_electrons))

    # see if any of the combos is valid and results in a lower numerical combination than the original
    for combo in combos:
        if _is_valid_combo(combo, mol, orig_distances):
            combo = sorted(itertools.chain.from_iterable(combo))
            if combo < selected_group:
                selected_group = combo

    # add the minimized unpaired electron positions to the u-layer:
    new_u_layer.extend(selected_group)

    return sorted(new_u_layer)


def _create_U_layer(mol, auxinfo):
    """
    Creates a string with the positions of the atoms that bear unpaired electrons. The string
    can be used to complement the InChI with an additional layer that allows for the differentiation
    between structures with multiple unpaired electrons.

    The string is composed of a prefix ('u') followed by the positions of each of the unpaired electrons,
    sorted in numerical order.

    Example:
    - methyl radical ([CH3]) : u1
    - triplet methylene biradical ([CH2]) : u1,1
    - ethane-1,2-diyl biradical ([CH2][CH2]): u1,2

    When the molecule does not bear any unpaired electrons, None is returned.

    """

    cython.declare(
        minmol=Molecule,
        # rdkitmol=,
        u_layer=list,
        i=int,
        at=Atom,
        equivalent_atoms=list,
    )

    if mol.getRadicalCount() == 0:
        return None
    elif mol.getFormula() == 'H':
        return U_LAYER_PREFIX + '1'

    # create preliminary u-layer:
    u_layer = []
    for i, at in enumerate(mol.atoms):
        u_layer.extend([i + 1] * at.radicalElectrons)

    # extract equivalent atom pairs from E-layer of auxiliary info:
    equivalent_atoms = _parse_E_layer(auxinfo)
    if equivalent_atoms:
        # select lowest u-layer:
        u_layer = _find_lowest_u_layer(mol, u_layer, equivalent_atoms)

    return (U_LAYER_PREFIX + ','.join(map(str, u_layer)))


def _find_lowest_p_layer(minmol, p_layer, equivalent_atoms):
    """
    Permute the equivalent atoms and return the combination with the
    lowest p-layer.

    TODO: The presence of unpaired electrons complicates stuff.
    """
    return p_layer


def _create_P_layer(mol, auxinfo):
    """
    Creates a string with the positions of the atoms that bear an unexpected number of lone pairs. The string
    can be used to complement the InChI with an additional layer that allows for the differentiation
    between structures with lone pairs.

    The string is composed of a prefix ('P_LAYER_PREFIX') followed by the positions of each of the atoms with an
    unexpected number of lone pairs, sorted in numerical order.

    Example:
    - singlet methylene biradical ([CH2]) : 'P_LAYER_PREFIX'1

    When the molecule does not bear any atoms with an unexpected number of lone pairs,
    None is returned.
    """

    # create preliminary p-layer:
    p_layer = []
    for i, at in enumerate(mol.atoms):
        try:
            exp = elements.PeriodicSystem.lone_pairs[at.symbol]
        except KeyError:
            raise Exception("Unrecognized element: {}".format(at.symbol))
        else:
            if at.lonePairs != exp:
                if at.lonePairs == 0:
                    p_layer.append('{}{}'.format(i, '(0)'))
                else:
                    p_layer.extend([i + 1] * at.lonePairs)

    # extract equivalent atom pairs from E-layer of auxiliary info:
    equivalent_atoms = _parse_E_layer(auxinfo)
    if equivalent_atoms:
        # select lowest u-layer:
        p_layer = _find_lowest_p_layer(mol, p_layer, equivalent_atoms)

    if p_layer:
        return (P_LAYER_PREFIX + P_LAYER_SEPARATOR.join(map(str, p_layer)))
    else:
        return None


def create_augmented_layers(mol):
    """
    The indices in the string refer to the atom indices in the molecule, according to the atom order
    obtained by sorting the atoms using the InChI canonicalization algorithm.

    First a deep copy is created of the original molecule and hydrogen atoms are removed from the molecule.
    Next, the molecule is converted into an InChI string, and the auxiliary information of the inchification
    procedure is retrieved.

    The N-layer is parsed and used to sort the atoms of the original order according
    to the order in the InChI. In case, the molecule contains atoms that cannot be distinguished
    with the InChI algorithm ('equivalent atoms'), the position of the unpaired electrons is changed
    as to ensure the atoms with the lowest indices are used to compose the string.
    """

    if mol.getRadicalCount() == 0 and not _has_unexpected_lone_pairs(mol):
        return None, None
    elif mol.getFormula() == 'H':
        return U_LAYER_PREFIX + '1', None
    else:
        molcopy = mol.copy(deep=True)

        hydrogens = filter(lambda at: at.number == 1, molcopy.atoms)
        for h in hydrogens:
            molcopy.removeAtom(h)

        rdkitmol = toRDKitMol(molcopy)
        _, auxinfo = Chem.MolToInchiAndAuxInfo(rdkitmol, options='-SNon')  # suppress stereo warnings

        # extract the atom numbers from N-layer of auxiliary info:
        atom_indices = _parse_N_layer(auxinfo)
        atom_indices = [atom_indices.index(i + 1) for i, atom in enumerate(molcopy.atoms)]

        # sort the atoms based on the order of the atom indices
        molcopy.atoms = [x for (y, x) in sorted(zip(atom_indices, molcopy.atoms), key=lambda pair: pair[0])]

        ulayer = _create_U_layer(molcopy, auxinfo)

        player = _create_P_layer(molcopy, auxinfo)

        return ulayer, player


##################################################################
# Methods for fixing molecules generated from an augmented InChI #
##################################################################

def _fix_triplet_to_singlet(mol, p_indices):
    """
    Iterates over the atoms and checks whether atoms bearing two unpaired electrons are
    also present in the p_indices list.

    If so, convert to the two unpaired electrons into a lone pair, and remove that atom
    index from the p_indices list.
    """

    for at in mol.atoms:
        index = mol.atoms.index(at) + 1
        if mol.getRadicalCount() == 2 and index in p_indices:
            at.lonePairs += 1
            at.radicalElectrons -= 2
            p_indices.remove(index)


def _convert_charge_to_unpaired_electron(mol, u_indices):
    """
    Iterates over the atoms foundin the parameter list and
    converts a unit of charge on atoms into an unpaired electron.

    Removes treated atoms from the parameter list.
    """
    for at in mol.atoms:
        at_index = mol.atoms.index(at) + 1
        if at.charge != 0 and at_index in u_indices:
            at.charge += 1 if at.charge < 0 else -1
            at.radicalElectrons += 1
            u_indices.remove(at_index)


def _convert_4_atom_3_bond_path(start):
    """
    Searches for 4-atom-3-bond [X=X-X=X+] paths starting from the parameter atom.
    If a path is found, the starting atom receives an unpaired electron while
    the bonds in the delocalization path are "inverted". A unit of charge on the
    end atom is neutralized and a lone pair is added.
    """
    path = pathfinder.find_butadiene_end_with_charge(start)

    if path is not None:
        start.radicalElectrons += 1
        end = path[-1]
        end.charge += 1 if end.charge < 0 else -1
        end.lonePairs += 1

        # filter bonds from path and convert bond orders:
        bonds = path[1::2]  # odd
        for bond in bonds[::2]:  # even
            assert isinstance(bond, Bond)
            bond.decrementOrder()
        for bond in bonds[1::2]:  # odd bonds
            assert isinstance(bond, Bond)
            bond.incrementOrder()

        return True

    return False


def _convert_3_atom_2_bond_path(start, mol):
    """
    Searches for 3-atom-2-bond [X=X-X+] paths paths starting from the parameter atom.
    If a correct path is found, the starting atom receives an unpaired electron while
    the bonds in the delocalization path are "inverted". A unit of charge on the
    end atom is neutralized and a lone pair is added.

    If it turns out the path was invalid, the actions are reverted, and another path
    is tried instead.

    To facilitate reverting the changes, we use a reaction recipe and populate it
    with a number of actions that reflect the changes in bond orders and unpaired
    electrons that the molecule should undergo.
    """
    from rmgpy.data.kinetics.family import ReactionRecipe

    def is_valid(mol):
        """Check if total bond order of oxygen atoms is smaller than 4."""

        for at in mol.atoms:
            if at.number == 8:
                order = at.getBondOrdersForAtom()
                not_correct = order >= 4
                if not_correct:
                    return False

        return True

    paths = pathfinder.find_allyl_end_with_charge(start)

    for path in paths:
        # label atoms so that we can use the labels in the actions of the recipe
        for i, at in enumerate(path[::2]):
            at.label = str(i)
        # we have found the atom we are looking for
        recipe = ReactionRecipe()
        recipe.addAction(['GAIN_RADICAL', start.label, 1])

        end = path[-1]
        end_original_charge = end.charge

        # filter bonds from path and convert bond orders:
        bonds = path[1::2]  # odd elements
        for bond in bonds[::2]:  # even
            recipe.addAction(['CHANGE_BOND', bond.atom1.label, -1, bond.atom2.label])
        for bond in bonds[1::2]:  # odd
            recipe.addAction(['CHANGE_BOND', bond.atom1.label, 1, bond.atom2.label])

        end.charge += 1 if end.charge < 0 else -1
        recipe.applyForward(mol)

        if is_valid(mol):
            # unlabel atoms so that they never cause trouble downstream
            for i, at in enumerate(path[::2]):
                at.label = ''
            return True
        else:
            recipe.applyReverse(mol)
            end.charge = end_original_charge

            # unlabel atoms so that they never cause trouble downstream
            for i, at in enumerate(path[::2]):
                assert isinstance(at, Atom)
                at.label = ''

    return False


def _convert_delocalized_charge_to_unpaired_electron(mol, u_indices):
    """
    Iterates over the atom indices of the parameter list and searches
    a charged atom that is connected to that atom via some kind of
    delocalization path.

    """
    u_indices_copy = u_indices[:]
    for index in u_indices_copy:
        start = mol.atoms[index - 1]

        found = _convert_4_atom_3_bond_path(start)
        if found:
            u_indices.remove(index)
            continue

        found = _convert_3_atom_2_bond_path(start, mol)
        if found:
            u_indices.remove(index)
            continue


def _fix_adjacent_charges(mol):
    """
    Searches for pairs of charged atoms.
    Neutralizes one unit of charge on each atom,
    and increments the bond order of the bond in between
    the atoms.
    """
    for at in mol.atoms:
        if at.charge != 0:
            for neigh, bond in at.bonds.iteritems():
                if neigh.charge != 0:
                    bond.incrementOrder()
                    at.charge += 1 if at.charge < 0 else -1
                    neigh.charge += 1 if neigh.charge < 0 else -1


def _fix_charge(mol, u_indices):
    """
    Tries to fix a number of structural features in the molecule related to charge,
    based on the information from the parameter list of atom indices with unpaired electrons.
    """

    if not u_indices:
        return

    is_charged = sum([abs(at.charge) for at in mol.atoms]) != 0
    is_correct = mol.getRadicalCount() == (mol.multiplicity - 1)
    if mol.multiplicity < 3 or is_correct or not is_charged:
        return

    # converting charges to unpaired electrons for atoms in the u-layer
    _convert_charge_to_unpaired_electron(mol, u_indices)

    # convert neighboring atoms (or delocalized paths) to unpaired electrons
    _convert_delocalized_charge_to_unpaired_electron(mol, u_indices)

    _fix_adjacent_charges(mol)


def _reset_lone_pairs(mol, p_indices):
    """
    Iterates over the atoms of the molecule and
    resets the atom's lone pair count to the value stored in the p_indices list,
    or to the default value.

    """
    for at in mol.atoms:
        index = mol.atoms.index(at) + 1  # 1-based index
        count = p_indices.count(index)
        if count != 0:
            at.lonePairs = count
        else:
            order = at.getBondOrdersForAtom()
            at.lonePairs = (elements.PeriodicSystem.valence_electrons[
                                at.symbol] - order - at.radicalElectrons - at.charge) / 2


def _fix_oxygen_unsaturated_bond(mol, u_indices):
    """
    Searches for a radical or a charged oxygen atom connected to
    a closed-shell carbon via an unsatured bond.

    Decrements the unsatured bond,
    transfers the unpaired electron from O to C or
    converts the charge from O to an unpaired electron on C,
    increases the lone pair count of O to 2.

    Only do this once per molecule.
    """

    for at in mol.atoms:
        if at.isOxygen() and at.radicalElectrons == 1 and at.lonePairs == 1:
            bonds = mol.getBonds(at)
            oxygen = at
            for atom2, bond in bonds.iteritems():
                if bond.isTriple():
                    bond.decrementOrder()
                    oxygen.radicalElectrons -= 1
                    atom2.radicalElectrons += 1
                    oxygen.lonePairs += 1
                    return
        elif at.isOxygen() and at.charge == 1 and at.lonePairs == 1:
            bonds = mol.getBonds(at)
            oxygen = at

            start = oxygen
            # search for 3-atom-2-bond [X=X-X] paths
            paths = pathfinder.find_allyl_end_with_charge(start)
            for path in paths:
                end = path[-1]
                start.charge += 1 if start.charge < 0 else -1
                end.charge += 1 if end.charge < 0 else -1
                start.lonePairs += 1
                # filter bonds from path and convert bond orders:
                bonds = path[1::2]  # odd elements
                for bond in bonds[::2]:  # even bonds
                    assert isinstance(bond, Bond)
                    bond.decrementOrder()
                for bond in bonds[1::2]:  # odd bonds
                    assert isinstance(bond, Bond)
                    bond.incrementOrder()
                break
            else:
                for atom2, bond in bonds.iteritems():
                    if not bond.isSingle() and atom2.charge == 0:
                        oxygen.charge -= 1
                        if (mol.atoms.index(atom2) + 1) in u_indices:
                            bond.decrementOrder()
                            atom2.radicalElectrons += 1
                            u_indices.remove(mol.atoms.index(atom2) + 1)
                        oxygen.lonePairs += 1


def _is_unsaturated(mol):
    """
    Does the molecule have a bond that's not single?
    Eg. a bond that is double or triple or benzene
    """
    cython.declare(atom1=Atom,
                   atom2=Atom,
                   bonds=dict,
                   bond=Bond)
    for atom in mol.atoms:
        for bond in atom.bonds.itervalues():
            if not bond.isSingle():
                return True

    return False


def _convert_unsaturated_bond_to_triplet(bond):
    """
    Decrements the bond if it is unsatured, and adds an unpaired
    electron to each of the atoms connected by the bond.
    """
    if not bond.isSingle():
        for at in (bond.atom1, bond.atom2):
            at.radicalElectrons += 1
        bond.decrementOrder()
        return True
    return False


def _fix_mobile_h(mol, inchi, u1, u2):
    """

    Identifies a system of atoms bearing unpaired electrons and mobile hydrogens
    at the same time.

    The system will consist of a central atom that does not bear any mobile hydrogens,
    but that is bound to an atom that does bear a mobile hydrogen, called the "original atom".

    The algorithm identifies the "new partner" atom that is part of the mobile hydrogen
    system.

    Next, the mobile hydrogen is transferred from the original atom, to the new partner,
    and a bond is removed and added respectively.

    Finally, the central atom and the original atom will each receive an unpaired electron,
    and the bond between them will decrease in order.
    """

    mobile_hydrogens = _parse_H_layer(inchi)

    if mobile_hydrogens:
        # WIP: only consider the first system of mobile hydrogens:
        mobile_hydrogens = mobile_hydrogens[0]

        # find central atom:
        central, original, new_partner = swap(mobile_hydrogens, [u1, u2])

        central, original, new_partner = \
            mol.atoms[central - 1], mol.atoms[original - 1], mol.atoms[new_partner - 1]

        # search hydrogen atom and bond
        hydrogen = None
        for at, bond in original.bonds.iteritems():
            if at.number == 1:
                hydrogen = at
                mol.removeBond(bond)
                break

        new_h_bond = Bond(new_partner, hydrogen, order='S')
        mol.addBond(new_h_bond)

        mol.getBond(central, new_partner).decrementOrder()

        central.radicalElectrons += 1
        original.radicalElectrons += 1
        return True

    return False


def _fix_butadiene_path(start, end):
    """
    Searches for a 1,3-butadiene path between the start and end atom.
    Adds an unpaired electron to start and end atom, and "inverts" the bonds
    in between them.
    """
    path = pathfinder.find_butadiene(start, end)
    if path is not None:
        start.radicalElectrons += 1
        end.radicalElectrons += 1
        # filter bonds from path and convert bond orders:
        bonds = path[1::2]  # odd elements
        for bond in bonds[::2]:  # even bonds
            assert isinstance(bond, Bond)
            bond.decrementOrder()
        for bond in bonds[1::2]:  # odd bonds
            assert isinstance(bond, Bond)
            bond.incrementOrder()

        return True

    return False


def _fix_unsaturated_bond_to_biradical(mol, inchi, u_indices):
    """
    Convert an unsaturated bond (double, triple) into a bond
    with a lower bond order (single, double), and give an unpaired electron
    to each of the neighboring atoms, with indices referring to the 1-based
    index in the InChI string.
    """
    cython.declare(u1=cython.int, u2=cython.int)
    cython.declare(atom1=Atom, atom2=Atom)
    cython.declare(b=Bond)

    combos = itertools.combinations(u_indices, 2)

    isFixed = False
    for u1, u2 in combos:
        atom1 = mol.atoms[u1 - 1]  # convert to 0-based index for atoms in molecule
        atom2 = mol.atoms[u2 - 1]  # convert to 0-based index for atoms in molecule
        if mol.hasBond(atom1, atom2):
            b = mol.getBond(atom1, atom2)
            isFixed = _convert_unsaturated_bond_to_triplet(b)
            if isFixed:
                break

            else:
                isFixed = _fix_mobile_h(mol, inchi, u1, u2)
                if isFixed:
                    break
        else:
            isFixed = _fix_butadiene_path(atom1, atom2)
            if isFixed:
                break

    if isFixed:
        u_indices.remove(u1)
        u_indices.remove(u2)
        return mol
    else:
        raise Exception(
            'Could not convert an unsaturated bond into a biradical for the \
            indices {} provided in the molecule: {}.'
                .format(u_indices, mol.toAdjacencyList())
        )


def _fix_unsaturated_bond(mol, indices, aug_inchi):
    """
    Adds unpaired electrons to the molecule by converting unsaturated bonds into triplets.

    It does so by converting an unsaturated bond into a triplet, and verifying whether
    the total number of unpaired electrons matches the multiplicity.

    Finishes when all unsaturated bonds have been tried, or when there are no pairs
    of atoms that should be unpaired electrons left.
    """

    correct = mol.getRadicalCount() == (mol.multiplicity - 1)

    if not correct and not indices:
        raise Exception('Cannot correct {} based on {} by converting unsaturated bonds into unpaired electrons...' \
                        .format(mol.toAdjacencyList(), aug_inchi))

    unsaturated = _is_unsaturated(mol)

    while not correct and unsaturated and len(indices) > 1:
        mol = _fix_unsaturated_bond_to_biradical(mol, aug_inchi.inchi, indices)
        correct = mol.getRadicalCount() == (mol.multiplicity - 1)
        unsaturated = _is_unsaturated(mol)


def _check_molecule(mol, aug_inchi):
    """
    Check if the molecular structure is correct.

    Checks whether the multiplicity contained in the augmented inchi,
    corresponds to the number of unpaired electrons + 1 found in the molecule.

    Checks whether the valence of each atom is compatible with the bond order,
    number of unpaired electrons, lone pairs and charge.

    """
    cython.declare(inchi=str,
                   at=Atom
                   )

    ConsistencyChecker.check_multiplicity(mol.getRadicalCount(), mol.multiplicity)
    _, u_indices, _ = decompose_aug_inchi(str(aug_inchi))
    assert(mol.getRadicalCount() == len(u_indices))

    for at in mol.atoms:
        ConsistencyChecker.check_partial_charge(at)


def fix_molecule(mol, aug_inchi):
    """
    Fixes a number of structural features of the erroneous Molecule
    parsed by the backends, based on multiplicity and unpaired electron information
    stored in the augmented inchi.
    """

    u_indices = aug_inchi.u_indices[:] if aug_inchi.u_indices else []
    p_indices = aug_inchi.p_indices[:] if aug_inchi.p_indices else []

    # ignore atoms that bear already unpaired electrons:
    for i in set(u_indices[:]):
        atom = mol.atoms[i - 1]
        for _ in range(atom.radicalElectrons): u_indices.remove(i)

        # ignore atoms that bear already lone pairs:
    for i in set(p_indices[:]):
        atom = mol.atoms[i - 1]
        for _ in range(atom.lonePairs): p_indices.remove(i)

    _fix_triplet_to_singlet(mol, p_indices)

    _fix_charge(mol, u_indices)

    _reset_lone_pairs(mol, p_indices)

    _fix_oxygen_unsaturated_bond(mol, u_indices)

    _fix_unsaturated_bond(mol, u_indices, aug_inchi)

    _check_molecule(mol, aug_inchi)


class InChI(str):
    """InChI is a type of string in which the InChI=1 prefix is ignored."""

    def __new__(self, inchi):
        if not INCHI_PREFIX in inchi:
            raise InchiException('Not a valid InChI: {}'.format(inchi))

        return str.__new__(self, remove_inchi_prefix(inchi))


class AugmentedInChI(InChI):
    """AugmentedInChI is an InChI with inchi, and unpaired electron attributes."""

    def __init__(self, aug_inchi):
        super(AugmentedInChI, self).__init__()
        inchi, u_indices, p_indices = decompose_aug_inchi(str(self))

        self.inchi = str(inchi)

        # default to None
        self.u_indices = u_indices or None
        self.p_indices = p_indices or None
