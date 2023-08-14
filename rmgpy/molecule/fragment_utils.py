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

import numpy as np


def matches_resolve(matches, rr_ll_list):
    """
    Sort out the pair of fragments and correct the amount. If the pair
    contains additional cutting label, it will be added into new_r_l_moles
    for further matching with other pairs.
    """
    new_matches = []
    new_r_l_moles = []
    for match in matches:
        pair = match[0]
        value = match[1]
        l_frag, r_frag = pair

        if l_frag not in rr_ll_list:
            if r_frag not in rr_ll_list:
                # cases like (L-Y, X-R)
                new_matches.append((pair, value))
            else:
                # cases like (L-Y, R-U1-R)
                new_matches.append(((l_frag, r_frag, l_frag), value / 2.0))
        else:
            if r_frag not in rr_ll_list:
                # cases like (L-W1-L, X-R)
                new_matches.append(((r_frag, l_frag, r_frag), value / 2.0))
            else:
                # cases like (L-W1-L, R-U1-R)
                new_r_l_moles.append((pair, value / 2.0))

    return new_matches, new_r_l_moles


def shuffle(conc, seed=None):
    """
    Randomly shuffle a list of fragments
    """
    idx_arr = np.arange(len(conc))

    if seed is not None:
        np.random.seed(seed)
    np.random.shuffle(idx_arr)

    return [conc[idx] for idx in idx_arr]


def grind(conc, size):
    """
    Split fragment concentrations into several repeating concentration units with specified size
    """
    grinded_conc = []
    for label, c in conc:
        times = int(c / size)
        grinded_conc.extend([(label, size)] * times)

        if c - size * times > 0:
            grinded_conc.append((label, c - size * times))

    return grinded_conc


def match_concentrations_with_same_sums(conc1, conc2, diff_tol=1e-6):
    """
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
    """
    labels1 = [tup[0] for tup in conc1]
    labels2 = [tup[0] for tup in conc2]

    seq1 = [tup[1] for tup in conc1]
    seq2 = [tup[1] for tup in conc2]

    matches_seq = match_sequences(seq1, seq2, diff_tol)

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
        remain_conc1 = [(labels1[pin1], val1)] + conc1[(pin1 + 1) :]
        matches_conc.extend(remain_conc1)

    # if pin1 first reaches the end
    # let matches_conc match with remaining seq2
    elif pin1 == len(seq1) and pin2 < len(seq2):
        remain_conc2 = [(labels2[pin2], val2)] + conc2[(pin2 + 1) :]
        matches_conc = match_concentrations_with_different_sums(
            matches_conc, remain_conc2
        )

    # if pin1 and pin2 reach the ends at same time
    # matches_conc is ready to return
    return matches_conc


def match_sequences(seq1, seq2, diff_tol=1e-6):
    """
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
    """
    # check if sums are close to same
    sum_diff = sum(seq2) - sum(seq1)
    assert (
        abs(sum_diff / 1.0 / sum(seq1)) <= diff_tol
    ), "seq1 has different sum (diff={0}) than seq2.".format(sum_diff)

    # force the sum to be same if the difference
    # is small enough
    if sum_diff >= 0:
        seq1[-1] = seq1[-1] + sum_diff
    else:
        seq2[-1] = seq2[-1] - sum_diff

    # make cumulative sequences
    cum_seq1 = [seq1[0]]
    for item1 in seq1[1:]:
        cum_seq1.append(cum_seq1[-1] + item1)

    cum_seq2 = [seq2[0]]
    for item2 in seq2[1:]:
        cum_seq2.append(cum_seq2[-1] + item2)

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

        matches.append([matched_index_tup, matched_cum_value - previous_cum_value])

    return matches


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
            return_list.extend(flatten(i))
        else:
            return_list.append(i)
    return return_list
