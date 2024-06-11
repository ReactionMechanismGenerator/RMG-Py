from rmgpy.molecule.fragment import Fragment
from rmgpy.tools.canteramodel import Cantera
from rmgpy.chemkin import load_chemkin_file
import re
import os
import numpy as np
import matplotlib.pyplot as plt


def match_sequences(seq1, seq2, rtol=1e-6):
    '''
    Given two lists (each item is int or float):
    seq1 and seq2 with same sum, the method returns
    matched indices and values.
    Example:
    seq1 = [1, 3, 1]
    seq2 = [2, 1, 2]
    return: [[(0,0),1],
                    [(1,0),1],
                    [(1,1),1],
                    [(1,2),1],
                    [(2,2),1]]
    '''
    sum_diff = sum(seq2) - sum(seq1)
    assert (
        np.isclose(sum(seq1), sum(seq2), rtol=rtol)
    ), "seq1 has different sum (diff={0}) than seq2.".format(sum_diff)

    # force the sum to be same if the difference
    # is small enough
    if sum_diff >= 0:
        seq1[-1] = seq1[-1] + sum_diff
    else:
        seq2[-1] = seq2[-1] - sum_diff

    # make cumulative sequences
    cum_seq1 = np.cumsum(seq1)
    cum_seq2 = np.cumsum(seq2)

    # add index tags two both cumulative seqs
    pin1 = 0
    pin2 = 0
    matched_indices = []
    matched_cum_values = []
    while pin1 < len(cum_seq1) and pin2 < len(cum_seq2):
        matched_indices.append((pin1, pin2))

        if cum_seq1[pin1] > cum_seq2[pin2]:
            matched_cum_values.append(cum_seq2[pin2])
            pin2 += 1
        elif cum_seq1[pin1] < cum_seq2[pin2]:
            matched_cum_values.append(cum_seq1[pin1])
            pin1 += 1
        else:
            matched_cum_values.append(cum_seq2[pin2])
            pin1 += 1
            pin2 += 1

    # get matches
    matches = []
    for i in range(len(matched_indices)):
        matched_index_tup = matched_indices[i]
        matched_cum_value = matched_cum_values[i]
        if i == 0:
            previous_cum_value = 0
        else:
            previous_cum_value = matched_cum_values[i - 1]

        matches.append(
            [matched_index_tup, matched_cum_value - previous_cum_value])

    return matches


def match_concentrations_with_same_sums(conc1, conc2, rtol=1e-6):
    '''match_concentrations_with_same_sums
    Given two lists with each item to be a tuple
    (species label, concentration)
    conc1 and conc2 with same total concentrations,
    the method returns matched species labels and
    concentrations.
    Example:
    conc1 = [('a', 1),
                    ('b', 3),
                    ('c', 1)]
    conc2 = [('x', 2),
                    ('y', 1),
                    ('z', 2)]
    return: [(('a','x'),1),
                    (('b','x'),1),
                    (('b','y'),1),
                    (('b','z'),1),
                    (('c','z'),1)]
    '''
    labels1 = [tup[0] for tup in conc1]
    labels2 = [tup[0] for tup in conc2]

    seq1 = [tup[1] for tup in conc1]
    seq2 = [tup[1] for tup in conc2]

    matches_seq = FragList.match_sequences(seq1, seq2, rtol)

    matches_conc = []
    for match_seq in matches_seq:
        matched_label_index1 = match_seq[0][0]
        matched_label_index2 = match_seq[0][1]
        matched_value = match_seq[1]

        matched_label1 = labels1[matched_label_index1]
        matched_label2 = labels2[matched_label_index2]
        match_conc = ((matched_label1, matched_label2), matched_value)
        matches_conc.append(match_conc)
    return matches_conc


def match_concentrations_with_different_sums(conc1, conc2):
    """
    Given two lists with each item to be a tuple
    (species label, concentration)
    conc1 and conc2 with different total concentrations,
    the method returns matched species labels and
    concentrations.
    Example:
    conc1 = [('a', 1),
                    ('b', 3),
                    ('c', 1)]
    conc2 = [('x', 2),
                    ('y', 1),
                    ('z', 10)]
    return: [(('a','x', 'z', 'z'),1),
                    (('b','x', 'z', 'z'),1),
                    (('b','y', 'z', 'z'),1),
                    (('b','z', 'z'),1),
                    (('c','z', 'z'),1)]
    """
    labels1 = [tup[0] for tup in conc1]
    labels2 = [tup[0] for tup in conc2]

    seq1 = [tup[1] for tup in conc1]
    seq2 = [tup[1] for tup in conc2]

    matches_conc = []
    pin1 = 0
    pin2 = 0
    val1 = seq1[pin1]
    val2 = seq2[pin2]

    while True:
        if val1 > val2:
            match = ((labels1[pin1], labels2[pin2]), val2)
            matches_conc.append(match)
            val1 = val1 - val2
            pin2 += 1
            if pin2 == len(seq2):
                break
            val2 = seq2[pin2]
        elif val1 < val2:
            match = ((labels1[pin1], labels2[pin2]), val1)
            matches_conc.append(match)
            val2 = val2 - val1
            pin1 += 1
            if pin1 == len(seq1):
                break
            val1 = seq1[pin1]
        else:
            match = ((labels1[pin1], labels2[pin2]), val1)
            matches_conc.append(match)
            pin1 += 1
            pin2 += 1
            if pin1 == len(seq1):
                break
            val1 = seq1[pin1]
            if pin2 == len(seq2):
                break
            val2 = seq2[pin2]

    # if pin2 first reaches the end
    # append all the remaining seq1 to matches_conc
    if pin2 == len(seq2) and pin1 < len(seq1):
        remain_conc1 = [(labels1[pin1], val1)] + conc1[(pin1 + 1):]
        matches_conc.extend(remain_conc1)

    # if pin1 first reaches the end
    # let matches_conc match with remaining seq2
    elif pin1 == len(seq1) and pin2 < len(seq2):
        remain_conc2 = [(labels2[pin2], val2)] + conc2[(pin2 + 1):]
        matches_conc = FragList.match_concentrations_with_different_sums(
            matches_conc, remain_conc2
        )

    # if pin1 and pin2 reach the ends at same time
    # matches_conc is ready to return
    return matches_conc


def shuffle(conc, seed=None):
    """
    Randomly shuffle a list of fragments
    """
    idx_arr = np.arange(len(conc))

    if seed is not None:
        np.random.seed(seed)
    np.random.shuffle(idx_arr)

    return [conc[idx] for idx in idx_arr]


def flatten(combo):
    """
    Given a combo nested `tuple`, e.g.,
    ((('LY', 'XR'), ('LWL', 'RUR'))
    return a list of labels contained in
    the combo ['LY', 'XR', 'LWL', 'RUR']
    """
    return_list = []
    for i in combo:
        if isinstance(i, tuple):
            return_list.extend(FragList.flatten(i))
        else:
            return_list.append(i)
    return return_list


# label should match the desired merging l/'abel on frag2
def merge_frag_to_frag(frag1, frag2, label):
    from rmgpy.molecule import Bond
    from rmgpy.molecule.fragment import Fragment, CuttingLabel

    frag_spe1 = Fragment().from_smiles_like_string(frag1)
    frag_spe2 = Fragment().from_smiles_like_string(frag2)
    # find position of desired CuttingLabel
    # need to find CuttingLabel on frag2 first
    for vertex in frag_spe2.vertices:
        if isinstance(vertex, CuttingLabel):
            if vertex.symbol == label:
                cut2 = vertex

                atom2 = list(cut2.edges.keys())[0]
                frag_spe2.remove_atom(cut2)
                break

    if cut2.symbol[0] == 'L':
        Ctl = cut2.symbol.replace('L', 'R')
    else:  # that means this CuttingLabel is R something

        Ctl = cut2.symbol.replace('R', 'L')

    # merge to frag_spe1
    for vertex in frag_spe1.vertices:
        if isinstance(vertex, CuttingLabel):
            if vertex.symbol == Ctl:
                cut1 = vertex
                atom1 = list(cut1.edges.keys())[0]
                frag_spe1.remove_atom(cut1)
                break

    # new merged fragment
    new_frag = frag_spe1.merge(frag_spe2)
    new_frag.add_bond(Bond(atom1=atom1, atom2=atom2, order=1))
    new_frag = new_frag.copy(deep=True)
    new_frag.update()
    return new_frag  # return Fragment obtl


def merge_frag_list(to_be_merged):
    import os
    # merges fragments in list from right to left
    species_list = []
    ethylene = []
    newlist = []
    warnings = []

    while len(to_be_merged) > 1:

        # second to last fragmentin list
        frag1 = to_be_merged[-2].smiles
        frag2 = to_be_merged[-1].smiles  # last fragment in list

        if 'R' in frag1 and 'L' in frag2:
            newfrag = FragList.merge_frag_to_frag(frag1, frag2, 'L')

        elif 'L' in frag1 and 'R' in frag2:
            newfrag = FragList.merge_frag_to_frag(frag1, frag2, 'R')

        # warn user if last two fragments in list cannot be merged (no R/L
        # combo to be made)
        else:
            print('Warning! Could not merge fragments {} and {}'.format(
                frag1, frag2))

            if 'L' in frag1 and 'L' in frag2:
                newfrag = FragList.merge_frag_to_frag(
                    frag1.replace('L', 'R'), frag2, 'L')
        if len(to_be_merged) > 2:
            cut = len(to_be_merged) - 2
            newfraglist = to_be_merged[:cut]

            newfraglist.append(newfrag)
        elif len(to_be_merged) == 2:

            newfraglist = [newfrag]

        to_be_merged = newfraglist

    to_be_merged = newfraglist
    # newlist.append(newfraglist) # if done merging list, write final
    # structure to list of smiles structures

    # print('{}% of fragments fully merged...'.format(np.round(100*(i+1)/len(flattened_matches_random)),1))
# print(newfraglist)
    return newfraglist
