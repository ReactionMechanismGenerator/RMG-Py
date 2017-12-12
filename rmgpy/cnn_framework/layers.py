
import numpy as np
import keras.backend as K
import theano.tensor as T # should write custom back-end eventually, but this is quick fix
import theano
from keras import activations, initializations
from keras.engine.topology import Layer
from keras.layers import merge

class MoleculeConv(Layer):
	
	def __init__(self, units, inner_dim, depth = 2, init_output='uniform', 
			activation_output='softmax', init_inner='identity',
			activation_inner='linear', scale_output=0.01, padding=False, **kwargs):
		if depth < 1:
			quit('Cannot use MoleculeConv with depth zero')
		self.init_output = initializations.get(init_output)
		self.activation_output = activations.get(activation_output)
		self.init_inner = initializations.get(init_inner)
		self.activation_inner = activations.get(activation_inner)
		self.units = units
		self.inner_dim = inner_dim
		self.depth = depth
		self.scale_output = scale_output
		self.padding = padding

		self.initial_weights = None
		self.input_dim = 4 # each entry is a 3D N_atom x N_atom x N_feature tensor
		if self.input_dim:
			kwargs['input_shape'] = (None, None, None,) # 3D tensor for each input
		#self.input = K.placeholder(ndim = 4)
		super(MoleculeConv, self).__init__(**kwargs)


	def build(self, input_shape):
		'''Builds internal weights and paramer attribute'''
		# NOTE: NEED TO TILE AND EVALUATE SO THAT PARAMS CAN BE VARIABLES
		# OTHERWISE K.GET_VALUE() DOES NOT WORK

		# Define template weights for inner FxF
		W_inner = self.init_inner((self.inner_dim, self.inner_dim))
		b_inner = K.zeros((1, self.inner_dim))
		# Initialize weights tensor 
		self.W_inner = K.variable(T.tile(W_inner, (self.depth + 1, 1, 1)).eval())
		self.W_inner.name = 'T:W_inner'
		self.b_inner = K.variable(T.tile(b_inner, (self.depth + 1, 1, 1)).eval())
		self.b_inner.name = 'T:b_inner'
		# # Concatenate third dimension (depth) so different layers can have 
		# # different weights. Now, self.W_inner[#,:,:] corresponds to the 
		# # weight matrix for layer/depth #.

		# Define template weights for output FxL
		# W_output = self.init_output((self.inner_dim, self.units), scale = self.scale_output)
		W_output = self.init_output((self.inner_dim, self.units))
		b_output = K.zeros((1, self.units))
		# Initialize weights tensor
		self.W_output = K.variable(T.tile(W_output, (self.depth + 1, 1, 1)).eval())
		self.W_output.name = 'T:W_output'
		self.b_output = K.variable(T.tile(b_output, (self.depth + 1, 1, 1)).eval())
		self.b_output.name = 'T:b_output'
		# # Concatenate third dimension (depth) so different layers can have 
		# # different weights. Now, self.W_output[#,:,:] corresponds to the 
		# # weight matrix for layer/depth #.

		# Pack params
		self.trainable_weights = [self.W_inner, 
					   self.b_inner,
					   self.W_output,
					   self.b_output]

	def get_output_shape_for(self, input_shape):
		return (input_shape[0], self.units)

	def call(self, x, mask=None):
		(output, updates) = theano.scan(lambda x_one: self.get_output_singlesample(x_one), sequences = x)
		return output

	def get_output_singlesample(self, M):
		"""
		Given a molecule tensor M, calcualte its fingerprint
		"""

		# if incoming tensor M has padding 
		# remove padding first 
		if self.padding:
			rowsum = M.sum(axis = 0) 
			trim = rowsum[:, -1] 
			trim_to = T.eq(trim, 0).nonzero()[0][0] # first index with no bonds
			M = M[:trim_to, :trim_to, :] # reduced graph

		# dimshuffle to get diagonal items to
		# form atom matrix A
		(A_tmp, updates) = theano.scan(lambda x: x.diagonal(), sequences = M[:,:,:-1].dimshuffle((2, 0, 1)))
		# Now the attributes is (N_features x N_atom), so we need to transpose
		A = A_tmp.T

		# get connectivity matrix: N_atom * N_atom
		C = M[:, :, -1] + T.identity_like(M[:, :, -1])

		# get bond tensor: N_atom * N_atom * (N_features-1) 
		B_tmp = M[:, :, :-1] - A
		coeff = K.concatenate([M[:, :, -1:]]*self.inner_dim, axis = 2)
		B = merge([B_tmp, coeff], mode="mul")

		# Get initial fingerprint
		presum_fp = self.attributes_to_fp_contribution(A, 0)
		fp = K.sum(presum_fp, axis = 0) # sum across atom contributions

		# Iterate through different depths, updating atom matrix each time
		A_new = A
		for depth in range(self.depth):
			A_new = self.activation_inner(K.dot(K.dot(C, A_new) + K.sum(B, axis=1), 
												self.W_inner[depth+1, :, :]) 
										+ self.b_inner[depth+1, 0, :])

			presum_fp_new = self.attributes_to_fp_contribution(A_new, depth + 1)
			fp = fp + K.sum(presum_fp_new, axis = 0) 

		return fp


	def attributes_to_fp_contribution(self, attributes, depth):
		'''Given a 2D tensor of attributes where the first dimension corresponds to a single
		node, this method will apply the output sparsifying (often softmax) function and return
		the contribution to the fingerprint'''
		# Apply output activation function
		output_dot = K.dot(attributes, self.W_output[depth, :, :]) # ignore last attribute (bond flag)
		output_dot.name = 'output_dot'
		output_bias = self.b_output[depth, 0, :]
		output_bias.name = 'output_bias'
		output_activated = self.activation_output(output_dot + output_bias)
		output_activated.name = 'output_activated'
		return output_activated

	def get_config(self):
		config = {'units': self.units,
				  'inner_dim' : self.inner_dim,
				  'init_output' : self.init_output.__name__,
				  'init_inner' : self.init_inner.__name__,
				  'activation_inner': self.activation_inner.__name__,
				  'activation_output' : self.activation_output.__name__,
				  'input_dim': self.input_dim,
				  'depth' : self.depth}
		base_config = super(MoleculeConv, self).get_config()
		return dict(list(base_config.items()) + list(config.items()))