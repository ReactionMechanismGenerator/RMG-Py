from .predictor import *
from .layers import GraphFP
from keras.layers.core import Dense
import unittest
import os
import rmgpy

class Test_Predictor(unittest.TestCase):

	def setUp(self):

		self.predictor = Predictor()

	def test_model(self):

		self.predictor.build_model()
		predictor_model = self.predictor.model
		self.assertEqual(len(predictor_model.layers), 2)
		self.assertTrue(isinstance(predictor_model.layers[0], GraphFP))
		self.assertTrue(isinstance(predictor_model.layers[1], Dense))

		self.assertEqual(predictor_model.layers[0].inner_dim, 8)
		self.assertEqual(predictor_model.layers[0].output_dim, 512)

	def test_load_input(self):

		test_predictor_input = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'predictor_input.py'
											)
		self.predictor.load_input(test_predictor_input)

		predictor_model = self.predictor.model
		self.assertEqual(len(predictor_model.layers), 2)
		self.assertTrue(isinstance(predictor_model.layers[0], GraphFP))
		self.assertTrue(isinstance(predictor_model.layers[1], Dense))

		self.assertEqual(predictor_model.layers[0].inner_dim, 8)
		self.assertEqual(predictor_model.layers[0].output_dim, 512)