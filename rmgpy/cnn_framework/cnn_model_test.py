
from .cnn_model import *
from .layers import MoleculeConv
from keras.layers.core import Dense
import unittest
import rmgpy
import os

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
		self.assertTrue(isinstance(test_model.layers[0], MoleculeConv))
		self.assertTrue(isinstance(test_model.layers[1], Dense))

		self.assertEqual(test_model.layers[0].inner_dim, attribute_vector_size-1)
		self.assertEqual(test_model.layers[0].output_dim, embedding_size)

	def test_save_model(self):

		embedding_size = 300
		attribute_vector_size = 10
		hidden = 0
		test_model = build_model(embedding_size=embedding_size, 
								attribute_vector_size=attribute_vector_size,
								hidden=hidden
								)

		save_model_folder = os.path.join(os.path.dirname(rmgpy.__file__), 
											'cnn_framework',
											'test_data',  
											'save_model_test')
		if not os.path.exists(save_model_folder):
			os.mkdir(save_model_folder)
		
		fpath = os.path.join(save_model_folder, 'model')

		save_model(test_model, [1.0], [1.0], 1.0, 1.0, fpath)

		import shutil
		shutil.rmtree(save_model_folder)
