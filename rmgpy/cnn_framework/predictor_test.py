from .predictor import *
from .layers import GraphFP
from keras.layers.core import Dense
import unittest
import os
import rmgpy
from rmgpy.molecule.molecule import Molecule
from rmgpy.cnn_framework.molecule_tensor import get_attribute_vector_size

class Test_Predictor(unittest.TestCase):

	def setUp(self):

		self.predictor = Predictor()

	def test_model(self):

		self.predictor.build_model()
		predictor_model = self.predictor.model
		self.assertEqual(len(predictor_model.layers), 3)
		self.assertTrue(isinstance(predictor_model.layers[0], GraphFP))
		self.assertTrue(isinstance(predictor_model.layers[1], Dense))

		self.assertEqual(predictor_model.layers[0].inner_dim, get_attribute_vector_size()-1)
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

	def test_load_parameters(self):

		test_predictor_input = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'predictor_input.py'
											)
		self.predictor.load_input(test_predictor_input)

		param_path = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'weights.h5'
											)
		self.predictor.load_parameters(param_path)

		gfp = self.predictor.model.layers[0]
		dense1 = self.predictor.model.layers[1]
		dense2 = self.predictor.model.layers[2]

		self.assertAlmostEqual(gfp.W_inner.eval()[0][0][0], 1.050, 3)
		self.assertAlmostEqual(gfp.b_inner.eval()[0][0][0], 0.036, 3)
		self.assertAlmostEqual(gfp.W_output.eval()[0][0][0], 0.181, 3)
		self.assertAlmostEqual(gfp.b_output.eval()[0][0][0], 0.071, 3)

		self.assertAlmostEqual(dense1.W.eval()[0][0], 0.344, 3)
		self.assertAlmostEqual(dense1.b.eval()[0], 0.486, 3)

		self.assertAlmostEqual(dense2.W.eval()[0][0], 3.284, 3)
		self.assertAlmostEqual(dense2.b.eval()[0], 3.074, 3)
		
	def test_predict(self):
		
		test_predictor_input = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'predictor_input.py'
											)
		self.predictor.load_input(test_predictor_input)

		param_path = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'weights.h5'
											)
		self.predictor.load_parameters(param_path)

		mol_test = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')

		h298_predicted = self.predictor.predict(mol_test)

		self.assertAlmostEqual(h298_predicted, 30.16, 1)