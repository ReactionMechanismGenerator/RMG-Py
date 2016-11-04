
from .cnn_model import *
from .layers import GraphFP
from keras.layers.core import Dense
import unittest

class Test_CNN_Model(unittest.TestCase):

	def test_build_model(self):

		embedding_size = 300
		attribute_vector_size = 10
		hidden = 0
		test_model = build_model(embedding_size=embedding_size, 
								attribute_vector_size=attribute_vector_size,
								hidden=hidden
								)
		self.assertEqual(len(test_model.layers), 2)
		self.assertTrue(isinstance(test_model.layers[0], GraphFP))
		self.assertTrue(isinstance(test_model.layers[1], Dense))

		self.assertEqual(test_model.layers[0].inner_dim, attribute_vector_size-1)
		self.assertEqual(test_model.layers[0].output_dim, embedding_size)
