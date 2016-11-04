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
		self.assertEqual(len(predictor_model.layers), 3)
		self.assertTrue(isinstance(predictor_model.layers[0], GraphFP))
		self.assertTrue(isinstance(predictor_model.layers[1], Dense))

		self.assertEqual(predictor_model.layers[0].inner_dim, 32)
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
		self.assertEqual(len(predictor_model.layers), 3)
		self.assertTrue(isinstance(predictor_model.layers[0], GraphFP))
		self.assertTrue(isinstance(predictor_model.layers[1], Dense))
		self.assertTrue(isinstance(predictor_model.layers[2], Dense))

		gfp = self.predictor.model.layers[0]
		dense1 = self.predictor.model.layers[1]
		dense2 = self.predictor.model.layers[2]

		self.assertEqual(gfp.W_inner.shape.eval()[0], 6)
		self.assertEqual(gfp.W_inner.shape.eval()[1], 32)
		self.assertEqual(gfp.W_inner.shape.eval()[2], 32)
		self.assertEqual(gfp.b_inner.shape.eval()[0], 6)
		self.assertEqual(gfp.b_inner.shape.eval()[1], 1)
		self.assertEqual(gfp.b_inner.shape.eval()[2], 32)

		self.assertEqual(gfp.W_output.shape.eval()[0], 6)
		self.assertEqual(gfp.W_output.shape.eval()[1], 32)
		self.assertEqual(gfp.W_output.shape.eval()[2], 512)
		self.assertEqual(gfp.b_output.shape.eval()[0], 6)
		self.assertEqual(gfp.b_output.shape.eval()[1], 1)
		self.assertEqual(gfp.b_output.shape.eval()[2], 512)

		self.assertEqual(dense1.W.shape.eval()[0], 512)
		self.assertEqual(dense1.W.shape.eval()[1], 50)
		self.assertEqual(dense1.b.shape.eval()[0], 50)

		self.assertEqual(dense2.W.shape.eval()[0], 50)
		self.assertEqual(dense2.W.shape.eval()[1], 1)
		self.assertEqual(dense2.b.shape.eval()[0], 1)
