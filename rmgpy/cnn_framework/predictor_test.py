from .predictor import *
from .layers import GraphFP
from keras.layers.core import Dense
import unittest
import os
import rmgpy
from rmgpy.molecule.molecule import Molecule

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

		self.assertEqual(gfp.W_inner.shape.eval()[0], 3)
		self.assertEqual(gfp.W_inner.shape.eval()[1], 38)
		self.assertEqual(gfp.W_inner.shape.eval()[2], 38)
		self.assertEqual(gfp.b_inner.shape.eval()[0], 3)
		self.assertEqual(gfp.b_inner.shape.eval()[1], 1)
		self.assertEqual(gfp.b_inner.shape.eval()[2], 38)

		self.assertEqual(gfp.W_output.shape.eval()[0], 3)
		self.assertEqual(gfp.W_output.shape.eval()[1], 38)
		self.assertEqual(gfp.W_output.shape.eval()[2], 512)
		self.assertEqual(gfp.b_output.shape.eval()[0], 3)
		self.assertEqual(gfp.b_output.shape.eval()[1], 1)
		self.assertEqual(gfp.b_output.shape.eval()[2], 512)

		self.assertEqual(dense1.W.shape.eval()[0], 512)
		self.assertEqual(dense1.W.shape.eval()[1], 50)
		self.assertEqual(dense1.b.shape.eval()[0], 50)

		self.assertEqual(dense2.W.shape.eval()[0], 50)
		self.assertEqual(dense2.W.shape.eval()[1], 1)
		self.assertEqual(dense2.b.shape.eval()[0], 1)

	def test_specify_datasets(self):
		"""
		Test the datasets specification is done properly
		"""
		expected_datasets = [('sdata134k', 'polycyclic_2954_table'),
							('sdata134k', 'cyclic_O_only_table')]

		self.assertEqual(self.predictor.datasets, expected_datasets)

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

		self.assertAlmostEqual(gfp.W_inner.eval()[0][0][0], 0.950, 3)
		self.assertAlmostEqual(gfp.b_inner.eval()[0][0][0], 0.028, 3)
		self.assertAlmostEqual(gfp.W_output.eval()[0][0][0], -0.203, 3)
		self.assertAlmostEqual(gfp.b_output.eval()[0][0][0], -0.014, 3)

		self.assertAlmostEqual(dense1.W.eval()[0][0], 0.191, 3)
		self.assertAlmostEqual(dense1.b.eval()[0], -0.024, 3)

		self.assertAlmostEqual(dense2.W.eval()[0][0], 2.987, 3)
		self.assertAlmostEqual(dense2.b.eval()[0], 2.874, 3)
		
	def test_predict(self):
		"""
		Test predictor is predicting within a reasonable range
		we should change weights.h5 everytime we change feature space
		"""
		
		test_predictor_input = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'predictor_input.py'
											)
		self.predictor.load_input(test_predictor_input)
		self.assertTrue(self.predictor.add_extra_atom_attribute)
		self.assertFalse(self.predictor.add_extra_bond_attribute)

		param_path = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'weights.h5'
											)
		self.predictor.load_parameters(param_path)

		mol_test = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')

		h298_predicted = self.predictor.predict(mol_test)

		self.assertAlmostEqual(h298_predicted, 7.4, 0)