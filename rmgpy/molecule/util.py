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

import itertools
import re

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.fragment import Fragment


def get_element_count(obj):
    """Converts an (augmented) inchi or Molecule into a dictionary element -> count"""
    element_count = {}

    if isinstance(obj, str):
        assert 'InChI=1' in obj
        mf = obj.split('/')[1]
        pieces = re.findall(r'[A-Z][^A-Z]*', mf)  # split on capital letters

        for piece in pieces:
            match = re.match(r"([a-z]+)([0-9]*)", piece, re.I)
            if match:
                element, count = match.groups()
                if count is '':
                    count = 1
                if element in element_count:
                    element_count[element] += int(count)
                else:
                    element_count[element] = int(count)

        # For surface species, replace Pt with X again
        if 'Pt' in element_count:
            element_count['X'] = element_count['Pt']
            del element_count['Pt']

        return element_count

    elif isinstance(obj, Molecule) or isinstance(obj, Fragment):
        return obj.get_element_count()

    else:
        raise Exception


def partition(sample, list_of_samples):
    """
    Group indices from the parameter sample 
    that belong to the same list.

    Returns a list of lists with partitioned indices, and 
    a list of lists with the corresponding sample list they were found in.

    E.g.:
    sample : [1,3]
    list_of_samples : [[1,2,3,4],[5,6]]
    returns: [[1,3]], [[1,2,3,4]]

    Indices not part of any of the lists should be in singleton list and 
    have corresponding empty list:

    E.g.:
    sample : [7]
    list_of_samples : [[1,2,3,4],[5,6]]
    returns: [[7]], [[]]

    """

    partitions, sample_lists = [], []

    for s in sample:
        for one_sample_list in list_of_samples:
            if s in one_sample_list:
                try:
                    index = sample_lists.index(one_sample_list)
                    partitions[index].append(s)
                except ValueError:
                    partitions.append([s])
                    sample_lists.append(one_sample_list)

                break
        else:  # s does not belong to any list of samples
            partitions.append([s])
            sample_lists.append([])

    return partitions, sample_lists


def agglomerate(groups):
    """
    Iterates over the parameter list of lists, and identifies all singletons.
    A new list of lists is created in which all singletons are combined together,
    while the other lists consisting of more than 1 element are simply copied.
    The newly created collapsed list singletons is appended at the end.

    Example:

    [[1,2,3], [4], [5,6], [7]]

    Returns:
    [[1,2,3], [5,6], [4,7]]
    """

    result = [x for x in groups if len(x) > 1]
    singletons = [x for x in groups if len(x) == 1]

    # collapse singletons:
    singletons = list(itertools.chain.from_iterable(singletons))

    result.append(singletons)
    return result


def generate_combo(samples, sample_spaces):
    """
    First, generate combinations of i samples from the corresponding sample_spaces.
    Next, generate the cross product between the elements in the previously generated list.

    Finally, filter out the original combination from the result.
    """
    combos = []
    for sample, sample_space in zip(samples, sample_spaces):
        combos_one_sample = [list(tup) for tup in itertools.combinations(sample_space, len(sample))]
        combos.append(combos_one_sample)

    # cross product of the elements in the list:
    combos = [list(tup) for tup in itertools.product(*combos)]

    # don't add the original samples
    combos = [x for x in combos if x != samples]

    return combos


def swap(to_be_swapped, sample):
    """
    Identifies which index of the list samples  is present in
    the list to be swapped. 

    E.g.:
    to be swapped: [2,3]
    sample: [1,3]

    Returns: 
    1, 3, 2
    
    """

    to_be_swapped = set(to_be_swapped)
    sample = set(sample)

    original = (sample.intersection(to_be_swapped)).pop()
    central = (sample - to_be_swapped).pop()
    new_partner = (to_be_swapped - sample).pop()

    return central, original, new_partner

def generate_closed_shell_singlet(m: Molecule):
    for i in range(len(m.atoms)):
            m.atoms[i].id = i
    assert m.multiplicity == 1 and m.get_radical_count()>0
    radical_center_ids = [x.id for x in m.atoms if x.radical_electrons > 0]
    # remove radicals from radical centers (2)
    for radical_center_id in radical_center_ids:
        m.atoms[radical_center_id].decrement_radical()
    # add removed radicals (2) to one of the radical sites as a lone pair (1)
    m.atoms[radical_center_ids[0]].increment_lone_pairs()
    # pick the best resonance structure
    print([x.smiles for x in m.generate_resonance_structures()])
    return m.generate_resonance_structures()[0]

def generate_singlet_diradicals(m: Molecule):
    for i in range(len(m.atoms)):
        m.atoms[i].id = i
    
    singlet_diradicals = [] 
    for edge in m.get_all_edges():
        
        M = m.copy(deep=True) 
        if edge.get_order_num() > 1: # find a pi bond
            atom1_id = edge.atom1.id
            atom2_id = edge.atom2.id
            M.atoms[atom1_id].increment_radical()  # add a radical to each atom of the pi bond
            M.atoms[atom2_id].increment_radical()
            M.get_bond(M.atoms[atom1_id], M.atoms[atom2_id]).decrement_order() # remove 1 pi bond
            potential_singlet_diradicals = M.generate_resonance_structures()  # generate resonance structures
            
            for potential_singlet_diradical in potential_singlet_diradicals: # find all resonance structures with non-neighboring radical sites
                radical_center_ids = sorted([x.id for x in potential_singlet_diradical.atoms if x.radical_electrons==1]) 
                potential_singlet_diradical_edges = potential_singlet_diradical.get_all_edges()
                potential_singlet_diradical_edge_ids = [sorted([x.atom1.id, x.atom2.id]) for x in potential_singlet_diradical_edges]
                if radical_center_ids not in potential_singlet_diradical_edge_ids:
                    if potential_singlet_diradical not in singlet_diradicals:
                        singlet_diradicals.append(potential_singlet_diradical)
    return singlet_diradicals