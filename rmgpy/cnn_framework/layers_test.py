from rmgpy.cnn_framework import layers
import keras.backend as K
from keras.engine import Model, Input
from keras.models import Sequential
import numpy as np
from numpy.testing import assert_allclose
import inspect
from nose.plugins.attrib import attr
import unittest


class TestLayer(unittest.TestCase):

	def test_MoleculeConv1(self):
		units = 10
		attribute_size = 3
		atom_num = 9
		batch_size = 2

		layer_test(self, layers.MoleculeConv,
					kwargs={
						"units": units, 
						"inner_dim": attribute_size-1, 
					},
					input_shape=(batch_size, 
								atom_num, 
								atom_num,
								attribute_size)
					)

	def test_MoleculeConv2(self):

		units = 10

		M = np.array([
					[[ 1.,  0.,  0.],[ 2.,  2.,  1.], [ 0.,  0.,  0.]],

					[[ 1.,  2.,  1.], [ 2.,  0.,  0.], [ 3.,  1.,  1.]],

					[[ 0.,  0.,  0.], [ 2.,  1.,  1.], [ 3.,  0.,  0.]]
					])

		expected_fp = np.array([26.0]*units)

		attribute_size = M.shape[2]
		depth = 1

		layer_test(self, layers.MoleculeConv,
					kwargs={
						"units": units, 
						"inner_dim": attribute_size-1,
						"init_inner": 'identity',
						"init_output": 'one',
						"activation_inner": 'linear',
						"activation_output": 'linear',
						"depth": depth
					},
					input_data=np.array([M]),
					input_dtype=K.floatx(),
					expected_output = np.array([expected_fp])
					)

	def test_MoleculeConv3(self):

		units = 10

		M = np.array([
					[[ 1.,  0.,  0.],[ 2.,  2.,  1.], [ 0.,  0.,  0.]],

					[[ 1.,  2.,  1.], [ 2.,  0.,  0.], [ 3.,  1.,  1.]],

					[[ 0.,  0.,  0.], [ 2.,  1.,  1.], [ 3.,  0.,  0.]]
					])

		expected_fp = np.array([81.0]*units)

		attribute_size = M.shape[2]
		depth = 2

		layer_test(self, layers.MoleculeConv,
					kwargs={
						"units": units, 
						"inner_dim": attribute_size-1,
						"init_inner": 'identity',
						"init_output": 'one',
						"activation_inner": 'linear',
						"activation_output": 'linear',
						"depth": depth
					},
					input_data=np.array([M]),
					input_dtype=K.floatx(),
					expected_output = np.array([expected_fp])
					)

@attr('helper')
def layer_test(test_case, layer_cls, kwargs={}, input_shape=None, input_dtype=None,
			   input_data=None, expected_output=None,
			   expected_output_dtype=None, fixed_batch_size=False):
	"""Test routine for a layer with a single input tensor
	and single output tensor.
	"""
	if input_data is None:
		assert input_shape
		if not input_dtype:
			input_dtype = K.floatx()
		input_data_shape = list(input_shape)
		for i, e in enumerate(input_data_shape):
			if e is None:
				input_data_shape[i] = np.random.randint(1, 4)
		input_data = (10 * np.random.random(input_data_shape))
		input_data = input_data.astype(input_dtype)
	elif input_shape is None:
		input_shape = input_data.shape

	if expected_output_dtype is None:
		expected_output_dtype = input_dtype

	# instantiation
	layer = layer_cls(**kwargs)

	# test get_weights , set_weights
	weights = layer.get_weights()
	layer.set_weights(weights)

	# test and instantiation from weights
	if 'weights' in inspect.getargspec(layer_cls.__init__):
		kwargs['weights'] = weights
		layer = layer_cls(**kwargs)

	# test in functional API
	if fixed_batch_size:
		x = Input(batch_shape=input_shape, dtype=input_dtype)
	else:
		x = Input(shape=input_shape[1:], dtype=input_dtype)
	y = layer(x)
	test_case.assertEqual(K.dtype(y), expected_output_dtype)

	model = Model(input=x, output=y)
	model.compile('rmsprop', 'mse')

	expected_output_shape = layer.get_output_shape_for(input_shape)
	actual_output = model.predict(input_data)
	actual_output_shape = actual_output.shape
	for expected_dim, actual_dim in zip(expected_output_shape,
										actual_output_shape):
		if expected_dim is not None:
			test_case.assertEqual(expected_dim, actual_dim)
	if expected_output is not None:
		assert_allclose(actual_output, expected_output, rtol=1e-3)


	model = Sequential()
	model.add(layer)
	model.compile('rmsprop', 'mse')
	actual_output = model.predict(input_data)
	actual_output_shape = actual_output.shape
	for expected_dim, actual_dim in zip(expected_output_shape,
										actual_output_shape):
		if expected_dim is not None:
			test_case.assertEqual(expected_dim, actual_dim)
	if expected_output is not None:
		assert_allclose(actual_output, expected_output, rtol=1e-3)


	# for further checks in the caller function
	return actual_output