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
This module contains methods for filtering a list of Molecules representing a single Species,
keeping only the representative structures. Relevant for filtration of negligible mesomerism contributing structures.

The rules this module follows are (by order of importance):
1. Minimum overall deviation from the Octet Rule (elaborated for Dectet and Duodectet for hypervalance heteroatoms)
2. Additional charge separation is only important for radicals and if it makes a new radical site in the species
3. If a structure must have charge separation, negative charges will be assigned to more electronegative atoms, whereas
   positive charges will be assigned to less electronegative atoms
4. Opposite charges will be as close as possible to one another, and vice versa

(inspired by http://www.chem.ucla.edu/~harding/tutorials/resonance/imp_res_str.html)
"""

import logging

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.element import PeriodicSystem
from rmgpy.molecule.pathfinder import find_shortest_path
from rmgpy.exceptions import ResonanceError


def filter_structures(mol_list, mark_unreactive=True):
    """
    We often get too many resonance structures from the combination of all rules, particularly for species containing
    lone pairs. This method filters them out by minimizing the number of C/N/O/S atoms without a full octet.
    """
    if not all([(mol.multiplicity == mol_list[0].multiplicity) for mol in mol_list]):
        raise ValueError("Cannot filter structures with different multiplicities!")

    # Get an octet deviation list
    octet_deviation_list = get_octet_deviation_list(mol_list)

    # Filter mol_list using the octet rule and the respective octet deviation list
    filtered_list, charge_span_list = octet_filtration(mol_list, octet_deviation_list)

    # Filter by charge span
    filtered_list = charge_filtration(filtered_list, charge_span_list)

    if mark_unreactive:
        # Mark selected unreactive structures if OS and/or adjacent birad unidirectional transitions were used
        mark_unreactive_structures(filtered_list, mol_list)

    # Check that there's at least one reactive structure in the list
    check_reactive(filtered_list)

    return filtered_list


def get_octet_deviation_list(mol_list):
    """
    Returns the a list of octet deviations for a respective list of :class:Molecule objects
    """
    octet_deviation_list = []
    for mol in mol_list:
        octet_deviation_list.append(get_octet_deviation(mol))

    return octet_deviation_list


def get_octet_deviation(mol):
    """
    Returns the octet deviation for a :class:Molecule object
    """
    if not isinstance(mol, Molecule):
        raise ValueError("Octet deviation could only be determined for Molecule objects.")

    octet_deviation = 0  # This is the overall "score" for the molecule, summed across all C/N/O/S atoms
    for atom in mol.vertices:
        val_electrons = 2 * (int(atom.getBondOrdersForAtom()) + atom.lonePairs) + atom.radicalElectrons
        if atom.isCarbon():
            octet_deviation += abs(8 - val_electrons)  # expecting C to be near octet
        elif atom.isNitrogen():
            if atom.lonePairs:
                octet_deviation += abs(8 - val_electrons)  # expecting N p1/2/3 to be near octet
            else:
                octet_deviation += min(abs(10 - val_electrons), abs(8 - val_electrons))  # N p0 could be near octet or
                # dectet. N p0 could be closer to an octet rather than a dectet such as in O=[N+][O-]
            if val_electrons > 8:
                octet_deviation += 1  # penalty for N with valance greater than 8 (as in O=[N.]=O,
                # [NH2.]=[:NH.], N#[N.]O, CCN=N#N)
        elif atom.isOxygen():
            octet_deviation += abs(8 - val_electrons)  # expecting O to be near octet
            if atom.atomType.label in ['O4sc', 'O4dc', 'O4tc']:
                octet_deviation += 1  # penalty for O p1 c+1
                # as in [N-2][N+]#[O+], [O-]S#[O+], OS(S)([O-])#[O+], [OH+]=S(O)(=O)[O-], [OH.+][S-]=O.
                # [C-]#[O+] and [O-][O+]=O which are correct structures also get penalized here, but that's OK
                # since they are still eventually selected as representative structures according to the rules here.
        elif atom.isSulfur():
            if atom.lonePairs == 0:
                octet_deviation += abs(12 - val_electrons)  # duodectet on S p0, eg O=S(=O)(O)O val 12, O[S](=O)=O val 11
            elif atom.lonePairs == 1:
                octet_deviation += min(abs(8 - val_electrons), abs(10 - val_electrons))  # octet/dectet on S p1,
                # eg [O-][S+]=O val 8, O[S]=O val 9, OS([O])=O val 10
            elif atom.lonePairs == 2:
                octet_deviation += min(abs(8 - val_electrons), abs(10 - val_electrons))  # octet/dectet on S p2,
                # eg [S][S] val 7, OS[O] val 8, [NH+]#[N+][S-][O-] val 9, O[S-](O)[N+]#N val 10
            elif atom.lonePairs == 3:
                octet_deviation += abs(8 - val_electrons)  # octet on S p3, eg [S-][O+]=O
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
        # E.g. O=[:S..] is penalized, but [C]=C=O isn't.
        if (atom.radicalElectrons == 2 and
                ((atom.isNitrogen() and atom.lonePairs == 0)
                 or (atom.isOxygen() and atom.lonePairs in [0,1])
                 or (atom.isSulfur() and atom.lonePairs in [0,1]))):
            octet_deviation += 3

    return octet_deviation


def octet_filtration(mol_list, octet_deviation_list):
    """
    Returns the a filtered list based on the octet_deviation_list. Also computes and returns a charge_span_list.
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
    This is also calculated in the octet_filtration() method along with the octet filtration process
    """
    charge_span_list = []
    for mol in mol_list:
        charge_span_list.append(mol.getChargeSpan())

    return charge_span_list


def charge_filtration(filtered_list, charge_span_list):
    """
    Returns a new filtered_list, filtered based on charge_span_list, electronegativity and proximity considerations.
    If a structure must have charge separation, negative charges will be assigned to more electronegative atoms, whereas
    positive charges will be assigned to less electronegative atoms. Also, opposite charges will be as close as possible
    to one another, and vice versa.
    If the species is a radical, we first check whether keeping an extra charge span separation might be important for
    reactivity by relocating the radical site. If so, we keep these structures.
    For example:
    - Both of NO2's resonance structures will be kept: [O]N=O <=> O=[N+.][O-]
    - NCO will only have two resonance structures [N.]=C=O <=> N#C[O.], and will loose the third structure which has
      the same octet deviation, has a charge separation, but the radical site has already been considered: [N+.]#C[O-]
    - CH2NO keeps all three structures, since a new radical site is introduced: [CH2.]N=O <=> C=N[O.] <=> C=[N+.][O-]
    However, if the species is not a radical we only keep the structures with the minimal charge span.
    For example:
    - NSH will only keep N#S and not [N-]=[SH+]
    - The following species will loose two thirds of its resonance structures, which are charged: CS(=O)SC <=>
      CS(=O)#SC <=> C[S+]([O-]SC <=> CS([O-])=[S+]C <=> C[S+]([O-])#SC <=> C[S+](=O)=[S-]C
    - The azide structure is know to have three resonance structures: [NH-][N+]#N <=> N=[N+]=[N-] <=> [NH+]#[N+][N-2];
      here we'll lose the third one, which is theoretically "true", but doesn't contribute to reactivity.
    """
    min_charge_span = min(charge_span_list)
    if len(set(charge_span_list)) > 1:
        # Proceed only if there are structures with different charge spans and the species is a radical.
        if filtered_list[0].isRadical():
            charged_list = [filtered_mol for index, filtered_mol in enumerate(filtered_list) if
                            charge_span_list[index] == min_charge_span + 1]  # save the 2nd charge span layer
            filtered_list = [filtered_mol for index, filtered_mol in enumerate(filtered_list) if
                            charge_span_list[index] == min_charge_span]  # keep at least one charge span layer
            # Find the radical sites in all filtered_list structures:
            sorting_list = []
            for mol in filtered_list:
                for atom in mol.vertices:
                    if atom.radicalElectrons:
                        sorting_list.append(int(atom.sortingLabel))
            # Find unique radical sites in charged_list and append these structures to filtered_list:
            unique_charged_list = []
            for mol in charged_list:
                for atom in mol.vertices:
                    if atom.radicalElectrons and int(atom.sortingLabel) not in sorting_list:
                        unique_charged_list.append(mol)

            if unique_charged_list:
                # only keep structures that obey the electronegativity and the charge proximity rules
                indices_to_pop = []
                charge_distance_list = []  # matching i indices to unique_charged_list
                for i, mol in enumerate(unique_charged_list):
                    electroneg_positively_charged_atoms = electroneg_negatively_charged_atoms = 0
                    for atom in mol.vertices:
                        if atom.charge > 0:
                            electroneg_positively_charged_atoms += PeriodicSystem.electronegativity[atom.symbol]
                        elif atom.charge < 0:
                            electroneg_negatively_charged_atoms += PeriodicSystem.electronegativity[atom.symbol]
                    if electroneg_positively_charged_atoms > electroneg_negatively_charged_atoms:
                        # This condition is NOT hermetic: It is possible to think of a situation where one structure has
                        # several pairs of formally charged atoms, where one of the pairs isn't obeying the
                        # electronegativity rule, while the sum of the pairs does.
                        indices_to_pop.append(i)
                    # try to find well-defined pairs of formally-charged atoms to apply the proximity principle
                    # (opposite charges will be as close as possible to one another, and vice versa)
                    cumulative_opposite_charge_distance = cumulative_similar_charge_distance = 0
                    for atom1 in mol.vertices:
                        if atom1.charge:
                            for atom2 in mol.vertices:
                                if atom2.charge and atom2.sortingLabel > atom1.sortingLabel:
                                    # found two charged atoms
                                    if (atom1.charge > 0) ^ (atom2.charge > 0):  # xor
                                        # they have opposing signs when ONLY one is positive
                                        cumulative_opposite_charge_distance += len(find_shortest_path(atom1,atom2))
                                    else:
                                        # they have similar signs
                                        cumulative_similar_charge_distance += len(find_shortest_path(atom1,atom2))
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

                # pop structures that do not obey the electronegativity or charge proximity rules
                for i in reversed(xrange(len(unique_charged_list))):
                    # pop from the end, so indices won't change
                    if i in indices_to_pop:
                        unique_charged_list.pop(i)

                # append unique_charged_list to filtered_list
                filtered_list += unique_charged_list

        else:
            # There are structures with different charge spans, but the species is not a radical.
            # Additional charge span levels are removed
            filtered_list = [filtered for index, filtered in enumerate(filtered_list) if
                             charge_span_list[index] == min_charge_span]

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
