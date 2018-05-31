
import unittest
import rmgpy
from .predictor import Predictor
from .input import *
from .layers import MoleculeConv
from keras.layers.core import Dense

class Test_Input(unittest.TestCase):

	def test_read_input_file(self):

		predictor_test = Predictor()

		path = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'predictor_input.py'
											)
		read_input_file(path, predictor_test)

		predictor_model = predictor_test.model
		self.assertEqual(len(predictor_model.layers), 3)
		self.assertTrue(isinstance(predictor_model.layers[0], MoleculeConv))
		self.assertTrue(isinstance(predictor_model.layers[1], Dense))

		self.assertEqual(predictor_model.layers[0].inner_dim, 38)
		self.assertEqual(predictor_model.layers[0].units, 512)

