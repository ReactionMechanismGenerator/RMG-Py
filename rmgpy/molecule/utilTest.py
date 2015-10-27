
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