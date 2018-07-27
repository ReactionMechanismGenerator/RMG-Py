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
This module contains functions for filtering a list of Molecules representing a single Species,
keeping only the representative structures. Relevant for filtration of negligible mesomerism contributing structures.

The rules this module follows are (by order of importance):
1. Minimum overall deviation from the Octet Rule (elaborated for Dectet for sulfur as a third row element)
2. Additional charge separation is only allowed for radicals if it makes a new radical site in the species
3. If a structure must have charge separation, negative charges will be assigned to more electronegative atoms, whereas
   positive charges will be assigned to less electronegative atoms (charge stabilization)
4. Opposite charges will be as close as possible to one another, and vice versa (charge stabilization)

(inspired by http://www.chem.ucla.edu/~harding/tutorials/resonance/imp_res_str.html)
"""

import logging

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.element import PeriodicSystem
from rmgpy.molecule.pathfinder import find_shortest_path
from rmgpy.exceptions import ResonanceError


def filter_structures(mol_list, mark_unreactive=True, allow_expanded_octet=True, features=None):
    """
    We often get too many resonance structures from the combination of all rules, particularly for species containing
    lone pairs. This function filters them out by minimizing the number of C/N/O/S atoms without a full octet.
    """
    if not all([(mol.multiplicity == mol_list[0].multiplicity) for mol in mol_list]):
        raise ValueError("Cannot filter structures with different multiplicities!")

    # Get an octet deviation list
    octet_deviation_list = get_octet_deviation_list(mol_list, allow_expanded_octet=allow_expanded_octet)

    # Filter mol_list using the octet rule and the respective octet deviation list
    filtered_list, charge_span_list = octet_filtration(mol_list, octet_deviation_list)

    # Filter by charge
    filtered_list = charge_filtration(filtered_list, charge_span_list)

    # Filter aromatic structures
    if features is not None and features['isAromatic']:
        filtered_list = aromaticity_filtration(filtered_list, features)

    if not filtered_list:
        raise ResonanceError('Could not determine representative localized structures for species {0}'.format(
            mol_list[0].toSMILES()))

    if mark_unreactive:
        # Mark selected unreactive structures if OS and/or adjacent birad unidirectional transitions were used
        mark_unreactive_structures(filtered_list, mol_list)

    # Check that there's at least one reactive structure in the list
    check_reactive(filtered_list)

    return filtered_list


def get_octet_deviation_list(mol_list, allow_expanded_octet=True):
    """
    Returns the a list of octet deviations for a respective list of :class:Molecule objects
    """
    octet_deviation_list = []
    for mol in mol_list:
        octet_deviation_list.append(get_octet_deviation(mol, allow_expanded_octet=allow_expanded_octet))

    return octet_deviation_list


def get_octet_deviation(mol, allow_expanded_octet=True):
    """
    Returns the octet deviation for a :class:Molecule object
    if `allow_expanded_octet` is ``True`` (by default), then the function also considers dectet for
    third row elements (currently sulfur is the only hypervalance third row element in RMG)
    """
    from afm.fragment import Fragment, CuttingLabel
    if not isinstance(mol, (Molecule, Fragment)):
        raise TypeError("Octet deviation could only be determined for Molecule objects.")

    octet_deviation = 0  # This is the overall "score" for the molecule, summed across all non-H atoms
    for atom in mol.vertices:
        if atom.isHydrogen() or isinstance(atom, CuttingLabel):
            continue
        val_electrons = 2 * (int(atom.getBondOrdersForAtom()) + atom.lonePairs) + atom.radicalElectrons
        if atom.isCarbon() or atom.isNitrogen() or atom.isOxygen():
            octet_deviation += abs(8 - val_electrons)  # expecting C/N/O to be near octet
        elif atom.isSulfur():
            if not allow_expanded_octet:
                # If allow_expanded_octet is False, then adhere to the octet rule for sulfur as well.
                # This is in accordance with J. Chem. Educ., 1995, 72 (7), p 583, DOI: 10.1021/ed072p583
                # This results in O=[:S+][:::O-] as a representative structure for SO2 rather than O=S=O,
                # and in C[:S+]([:::O-])C as a representative structure for DMSO rather than CS(=O)C.
                octet_deviation += abs(8 - val_electrons)
            else:
                # If allow_expanded_octet is True, then do not adhere to the octet rule for sulfur
                # and allow dectet structures (but don't prefer duedectet).
                # This is in accordance with:
                # -  J. Chem. Educ., 1972, 49 (12), p 819, DOI: 10.1021/ed049p819
                # -  J. Chem. Educ., 1986, 63 (1), p 28, DOI: 10.1021/ed063p28
                # -  J. Chem. Educ., 1992, 69 (10), p 791, DOI: 10.1021/ed069p791
                # -  J. Chem. Educ., 1999, 76 (7), p 1013, DOI: 10.1021/ed076p1013
                # This results in O=S=O as a representative structure for SO2 rather than O=[:S+][:::O-],
                # and in CS(=O)C as a representative structure for DMSO rather than C[:S+]([:::O-])C.
                if atom.lonePairs <= 1:
                    octet_deviation += min(abs(8 - val_electrons), abs(10 - val_electrons))  # octet/dectet on S p[0,1]
                    # eg [O-][S+]=O, O[S]=O, OS([O])=O, O=S(=O)(O)O
                elif atom.lonePairs >= 2:
                    octet_deviation += abs(8 - val_electrons)  # octet on S p[2,3]
                    # eg [S][S], OS[O], [NH+]#[N+][S-][O-], O[S-](O)[N+]#N, S=[O+][O-]
            for atom2, bond in atom.bonds.iteritems():
                if atom2.isSulfur() and bond.isTriple():
                    octet_deviation += 0.5  # penalty for S#S substructures. Often times sulfur can have a triple
                    # bond to another sulfur in a structure that obeys the octet rule, but probably shouldn't be a
                    # correct resonance structure. This adds to the combinatorial effect of resonance structures
                    # when generating reactions, yet probably isn't too important for reactivity. The penalty value
                    # is 0.5 since S#S substructures are captured twice (once for each S atom).
                    # Examples: CS(=O)SC <=> CS(=O)#SC;
                    # [O.]OSS[O.] <=> [O.]OS#S[O.] <=> [O.]OS#[S.]=O; N#[N+]SS[O-] <=> N#[N+]C#S[O-]
        # Penalize birad sites only if they theoretically substitute a lone pair.
        # E.g., O=[:S..] is penalized, but [C..]=C=O isn't.
        if (atom.radicalElectrons >= 2 and
                ((atom.isNitrogen() and atom.lonePairs == 0)
                 or (atom.isOxygen() and atom.lonePairs in [0,1,2])
                 or (atom.isSulfur() and atom.lonePairs in [0,1,2]))):
            octet_deviation += 3

    return octet_deviation


def octet_filtration(mol_list, octet_deviation_list):
    """
    Returns a filtered list based on the octet_deviation_list. Also computes and returns a charge_span_list.
    Filtering using the octet deviation criterion rules out most unrepresentative structures. However, since some
    charge-strained species are still kept (e.g., [NH]N=S=O <-> [NH+]#[N+][S-][O-]), we also generate during the same
    loop a charge_span_list to keep track of the charge spans. This is used for further filtering.
    """
    filtered_list = []
    charge_span_list = []
    for index, mol in enumerate(mol_list):
        if octet_deviation_list[index] == min(octet_deviation_list):
            filtered_list.append(mol)
            charge_span_list.append(mol.getChargeSpan())

    return filtered_list, charge_span_list


def get_charge_span_list(mol_list):
    """
    Returns the a list of charge spans for a respective list of :class:Molecule objects
    This is also calculated in the octet_filtration() function along with the octet filtration process
    """
    charge_span_list = []
    for mol in mol_list:
        charge_span_list.append(mol.getChargeSpan())

    return charge_span_list


def charge_filtration(filtered_list, charge_span_list):
    """
    Returns a new filtered_list, filtered based on charge_span_list, electronegativity and proximity considerations.
    If structures with an additional charge layer introduce reactive sites (i.e., radicals or multiple bonds) they will
    also be considered.
    For example:
    - Both of NO2's resonance structures will be kept: [O]N=O <=> O=[N+.][O-]
    - NCO will only have two resonance structures [N.]=C=O <=> N#C[O.], and will loose the third structure which has
      the same octet deviation, has a charge separation, but the radical site has already been considered: [N+.]#C[O-]
    - CH2NO keeps all three structures, since a new radical site is introduced: [CH2.]N=O <=> C=N[O.] <=> C=[N+.][O-]
    - NH2CHO has two structures, one of which is charged since it introduces a multiple bond: NC=O <=> [NH2+]=C[O-]
    However, if the species is not a radical, or multiple bonds do not alter, we only keep the structures with the
    minimal charge span. For example:
    - NSH will only keep the N#S form and not [N-]=[SH+]
    - The following species will loose two thirds of its resonance structures, which are charged: CS(=O)SC <=>
      CS(=O)#SC <=> C[S+]([O-]SC <=> CS([O-])=[S+]C <=> C[S+]([O-])#SC <=> C[S+](=O)=[S-]C
    - Azide is know to have three resonance structures: [NH-][N+]#N <=> N=[N+]=[N-] <=> [NH+]#[N+][N-2];
      here we filter the third one out due to the higher charge span, which does not contribute to reactivity in RMG
    """
    min_charge_span = min(charge_span_list)
    if len(set(charge_span_list)) > 1:
        # Proceed if there are structures with different charge spans
        charged_list = [filtered_mol for index, filtered_mol in enumerate(filtered_list) if
                        charge_span_list[index] == min_charge_span + 1]  # save the 2nd charge span layer
        filtered_list = [filtered_mol for index, filtered_mol in enumerate(filtered_list) if
                         charge_span_list[index] == min_charge_span]  # the minimal charge span layer
        # Find the radical and multiple bond sites in all filtered_list structures:
        rad_sorting_list = []  # sortingLabels for radical sites
        mul_bond_sorting_list = []  # sortingLabels for multiple bind sites in the form of (atom1,atom2) tuples
        for mol in filtered_list:
            for atom in mol.vertices:
                if atom.radicalElectrons and int(atom.sortingLabel) not in rad_sorting_list:
                    rad_sorting_list.append(int(atom.sortingLabel))
                for atom2, bond in atom.bonds.iteritems():
                    # check if bond is multiple, store only from one side (atom1 < atom2) for consistency
                    if atom2.sortingLabel > atom.sortingLabel and bond.isDouble() or bond.isTriple():
                        mul_bond_sorting_list.append((int(atom.sortingLabel), int(atom2.sortingLabel)))
        # Find unique radical and multiple bond sites in charged_list and append to unique_charged_list:
        unique_charged_list = []
        for mol in charged_list:
            unique_charged_list.extend(find_unique_sites_in_charged_list(mol, rad_sorting_list, mul_bond_sorting_list))

        # Charge stabilization considerations for the case where there are several charge span layers
        # are checked here for filtered_list and unique_charged_list separately.
        filtered_list = stabilize_charges_by_electronegativity(filtered_list)
        filtered_list = stabilize_charges_by_proximity(filtered_list)
        if unique_charged_list:
            unique_charged_list = stabilize_charges_by_electronegativity(unique_charged_list, allow_empty_list=True)
            unique_charged_list = stabilize_charges_by_proximity(unique_charged_list)
            filtered_list.extend(unique_charged_list)

    if min_charge_span:
        # If the species has charge separation, apply charge stability considerations.
        # These considerations should be checked regardless of the existence of radical sites.
        # They should also be checked if len(set(charge_span_list)) == 1.
        filtered_list = stabilize_charges_by_electronegativity(filtered_list)
        filtered_list = stabilize_charges_by_proximity(filtered_list)

    return filtered_list


def find_unique_sites_in_charged_list(mol, rad_sorting_list, mul_bond_sorting_list):
    """
    A helper function for reactive site discovery in charged species
    """
    for atom in mol.vertices:
        if atom.radicalElectrons and int(atom.sortingLabel) not in rad_sorting_list:
            return [mol]
        for atom2, bond in atom.bonds.iteritems():
            if atom2.sortingLabel > atom.sortingLabel and (bond.isDouble() or bond.isTriple())\
                    and (int(atom.sortingLabel), int(atom2.sortingLabel)) not in mul_bond_sorting_list\
                    and not (atom.isSulfur() and atom2.isSulfur()):
                # We check that both atoms aren't S, otherwise we get [S.-]=[S.+] as a structure of S2 triplet
                return [mol]
    return []


def stabilize_charges_by_electronegativity(mol_list, allow_empty_list=False):
    """
    Only keep structures that obey the electronegativity rule. If a structure must have charge separation, negative
    charges will be assigned to more electronegative atoms, and vice versa.
    If allow_empty_list is set to ``False`` (default), this function will not return an empty list. If it is set
    to ``True`` and all structures in `mol_list` violate the electronegativity heuristic, the original `mol_list`
    is returned (examples: [C-]#[O+], CS, [NH+]#[C-], [OH+]=[N-], [C-][S+]=C violate this heuristic).
    """
    indices_to_pop = []
    mol_list_copy = list(mol_list)
    for i, mol in enumerate(mol_list):
        electroneg_positively_charged_atoms = electroneg_negatively_charged_atoms = 0
        for atom in mol.vertices:
            if atom.charge > 0:
                electroneg_positively_charged_atoms += PeriodicSystem.electronegativity[atom.symbol] * abs(atom.charge)
                if atom.isOxygen():
                    for atom2 in atom.edges.keys():
                        if atom2.isFluorine() and atom2.charge < 0:
                            break
                    else:
                        electroneg_positively_charged_atoms += 1  # penalty for positively charged O not adjacent to F-,
                    # as in [N-2][N+]#[O+], [O-]S#[O+], OS(S)([O-])#[O+], [OH+]=S(O)(=O)[O-], [OH.+][S-]=O.
                    # [C-]#[O+] and [O-][O+]=O, which are correct structures, also get penalized here, but that's OK
                    # since they are still eventually selected as representative structures according to the rules here.
            elif atom.charge < 0:
                electroneg_negatively_charged_atoms += PeriodicSystem.electronegativity[atom.symbol] * abs(atom.charge)
        if electroneg_positively_charged_atoms > electroneg_negatively_charged_atoms:
            # Filter structures in which more electronegative atoms are positively charged.
            # This condition is NOT hermetic: It is possible to think of a situation where one structure has
            # several pairs of formally charged atoms, where one of the pairs isn't obeying the
            # electronegativity rule, while the sum of the pairs does.
            indices_to_pop.append(i)
    for i in reversed(xrange(len(mol_list))):  # pop starting from the end, so indices won't change
        if i in indices_to_pop:
            mol_list.pop(i)
    if mol_list or allow_empty_list:
        return mol_list
    else:
        return mol_list_copy


def stabilize_charges_by_proximity(mol_list):
    """
    Only keep structures that obey the charge proximity rule.
    Opposite charges will be as close as possible to one another, and vice versa.
    """
    indices_to_pop = []
    charge_distance_list = []  # indices match mol_list
    for i, mol in enumerate(mol_list):
        # Try finding well-defined pairs of formally-charged atoms to apply the proximity principle
        # (opposite charges will be as close as possible to one another, and vice versa)
        cumulative_opposite_charge_distance = cumulative_similar_charge_distance = 0
        for atom1 in mol.vertices:
            if atom1.charge:
                for atom2 in mol.vertices:
                    if atom2.charge and atom2.sortingLabel > atom1.sortingLabel:
                        # found two charged atoms
                        if (atom1.charge > 0) ^ (atom2.charge > 0):  # xor
                            # they have opposing signs when ONLY one is positive
                            cumulative_opposite_charge_distance += len(find_shortest_path(atom1, atom2))
                        else:
                            # they have similar signs
                            cumulative_similar_charge_distance += len(find_shortest_path(atom1, atom2))
        charge_distance_list.append([cumulative_opposite_charge_distance,
                                     cumulative_similar_charge_distance])
    min_cumulative_opposite_charge_distance = min([distances[0] for distances in charge_distance_list]
                                                  or [0])  # in Python 3 use `min(list, default=0)`
    for i, distances in enumerate(charge_distance_list):
        # after generating the charge_distance_list, iterate through it and mark structures to pop
        if distances[0] > min_cumulative_opposite_charge_distance:
            indices_to_pop.append(i)
    max_cumulative_similar_charge_distance = max([distances[1] for i, distances in
                                                  enumerate(charge_distance_list) if i not in indices_to_pop] or [0])
    for i, distances in enumerate(charge_distance_list):
        if distances[0] < max_cumulative_similar_charge_distance:
            indices_to_pop.append(i)
    for i in reversed(xrange(len(mol_list))):  # pop starting from the end, so indices won't change
        if i in indices_to_pop:
            mol_list.pop(i)
    return mol_list


def aromaticity_filtration(mol_list, features):
    """
    Returns a filtered list of molecules based on heuristics for determining
    representative aromatic resonance structures.

    For monocyclic aromatics, Kekule structures are removed, with the
    assumption that an equivalent aromatic structure exists. Non-aromatic
    structures are maintained if they present new radical sites. Instead of
    explicitly checking the radical sites, we only check for the SDSDSD bond
    motif since radical delocalization will disrupt that pattern.

    For polycyclic aromatics, structures without any benzene bonds are removed.
    The idea is that radical delocalization into the aromatic pi system is
    unfavorable because it disrupts aromaticity. Therefore, structures where
    the radical is delocalized so far into the molecule such that none of the
    rings are aromatic anymore are not representative. While this isn't strictly
    true, it helps reduce the number of representative structures by focusing
    on the most important ones.
    """
    # Start by selecting all aromatic resonance structures
    filtered_list = []
    other_list = []
    for mol in mol_list:
        if mol.isAromatic():
            filtered_list.append(mol)
        else:
            other_list.append(mol)

    if not features['isPolycyclicAromatic']:
        # Look for structures that don't have standard SDSDSD bond orders
        for mol in other_list:
            # Check all 6 membered rings
            rings = [ring for ring in mol.getRelevantCycles() if len(ring) == 6]
            for ring in rings:
                bond_list = mol.get_edges_in_cycle(ring)
                bond_orders = ''.join([bond.getOrderStr() for bond in bond_list])
                if bond_orders == 'SDSDSD' or bond_orders == 'DSDSDS':
                    break
            else:
                filtered_list.append(mol)

    return filtered_list


def mark_unreactive_structures(filtered_list, mol_list):
    """
    Mark selected structures in filtered_list with the Molecule.reactive flag set to `False` (it is `True` by default)
    Changes the filtered_list object, and does not return anything
    """
    # sort all structures in filtered_list so that the reactive ones are first
    filtered_list.sort(key=lambda mol: mol.reactive, reverse=True)

    # Make sure that the (first) original structure is always first in the list (unless it was filtered out).
    # Important whenever Species.molecule[0] is expected to be used (e.g., training reactions) after generating
    # resonance structures. However, if it was filtered out, it should be appended to the end of the list.
    for index, filtered in enumerate(filtered_list):
        if filtered.copy(deep=True).isIsomorphic(mol_list[0].copy(deep=True)):
            filtered_list.insert(0, filtered_list.pop(index))
            break
    else:
        # Append the original structure to list and set `reactive` to `False`.
        # This structure may very well deviate from the octet rule or have other attributes by which it should have
        # been filtered out. However, for processing reactions (e.g., degeneracy calculations) it should be kept
        # (e.g., [::N]O <=> [::N][::O.] + [H.], where [::N][::O.] should be recognized as [:N.]=[::O]).
        mol = mol_list[0]
        logging.debug("Setting the unrepresentative resonance structure {0} as unreactive in species {1}.".format(
            mol.toSMILES(),filtered_list[0].toSMILES()))
        logging.debug("Unreactive structure:\n{0}\nA representative reactive structure in this species:\n{1}\n".format(
            mol.toAdjacencyList(),filtered_list[0].toAdjacencyList()))
        mol.reactive = False
        filtered_list.append(mol)


def check_reactive(filtered_list):
    """
    Check that there's at least one reactive structure in the returned list.
    If not, raise an error (does not return anything)
    """
    if not any([mol.reactive for mol in filtered_list]):
        logging.info('\n\n')
        logging.error('No reactive structures were attributed to species {0}'.format(filtered_list[0].toSMILES()))
        for mol in filtered_list:
            logging.info('Structure: {0}\n{1}Reactive: {2}'.format(mol.toSMILES(),mol.toAdjacencyList(),mol.reactive))
        logging.info('\n')
        raise ResonanceError('Each species must have at least one reactive structure. Something probably went wrong'
                             ' when exploring resonance structures for species {0}'.format(filtered_list[0].toSMILES()))
