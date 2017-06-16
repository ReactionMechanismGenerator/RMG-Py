################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import unittest
from scipy.special import comb

from .util import *

class PartitionTest(unittest.TestCase):

    def test_singleton(self):
        """
        Test that a index not part of the parameter list, results in a key-value pair with
        an empty list.
        """
        indices = [7]
        list_of_samples = [[1,2,3,4],[5,6]]
        expected_partitions, expected_sample_lists = [[7]], [[]]

        partitions, sample_lists = partition(indices, list_of_samples)

        self.assertEquals(partitions, expected_partitions)
        self.assertEquals(sample_lists, expected_sample_lists)

    def test_2_elements_in_1_layer(self):
        indices = [1,3]
        list_of_samples = [[1,2,3,4],[5,6]]
        expected_partitions, expected_sample_lists = [[1,3]], [[1,2,3,4]]

        partitions, sample_lists = partition(indices, list_of_samples)

        self.assertEquals(partitions, expected_partitions)
        self.assertEquals(sample_lists, expected_sample_lists)

    def test_2_elements_in_2_layers(self):
        indices = [1,5]
        list_of_samples = [[1,2,3,4],[5,6]]
        expected_partitions, expected_sample_lists = [[1], [5]], [[1,2,3,4], [5,6]]

        partitions, sample_lists = partition(indices, list_of_samples)

        self.assertEquals(partitions, expected_partitions)
        self.assertEquals(sample_lists, expected_sample_lists)

    def test_3_elements_in_2_layers(self):
        indices = [1,4,5]
        list_of_samples = [[1,2,3,4],[5,6]]
        expected_partitions, expected_sample_lists = [[1,4], [5]], [[1,2,3,4], [5,6]]

        partitions, sample_lists = partition(indices, list_of_samples)

        self.assertEquals(partitions, expected_partitions)
        self.assertEquals(sample_lists, expected_sample_lists)

    def test_3_elements_in_2_layers_1_singleton(self):
        indices = [1,5,7]
        list_of_samples = [[1,2,3,4],[5,6]]
        expected_partitions, expected_sample_lists = [[1], [5], [7]], [[1,2,3,4], [5,6], []]

        partitions, sample_lists = partition(indices, list_of_samples)

        self.assertEquals(partitions, expected_partitions)
        self.assertEquals(sample_lists, expected_sample_lists)

class AgglomerateTest(unittest.TestCase):

    def test_normal(self):
        groups = [[1,2,3], [4], [5,6], [7]]
        agglomerates = agglomerate(groups)
        expected = [[1,2,3], [5,6], [4,7]]

        self.assertEquals(agglomerates, expected)

class ComboGeneratorTest(unittest.TestCase):
    def test_2_elements(self):

        samples = [[1,2,3],[6]]
        sample_spaces = [[1,2,3,4], [5,6]]

        combos = generate_combo(samples, sample_spaces)
        
        expected = 1
        for sample, sample_space in zip(samples, sample_spaces):
            expected *= comb(len(sample_space), len(sample), exact=True)

        # we leave out the original combination
        expected -= 1

        self.assertEquals(len(combos), expected)        

class SwapTest(unittest.TestCase):

    def test_2_elements_sets(self):
        to_be_swapped = [2,3]
        sample = [1,3]

        result = swap(to_be_swapped, sample)
        expected = (1,3,2)
        self.assertEquals(result, expected)
