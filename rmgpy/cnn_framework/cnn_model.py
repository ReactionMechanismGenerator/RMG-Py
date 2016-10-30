from .layers import GraphFP
from keras.models import Sequential
from keras.layers.core import Dense
from keras.optimizers import RMSprop, Adam

import logging

def build_model(embedding_size=512, attribute_vector_size=9, depth=2, scale_output=0.05, padding=True, 
				hidden=0, hidden_activation='tanh',
				output_activation='linear', output_size=1, 
				lr=0.01, optimizer='adam', loss='mse'):
	
	model = Sequential()

	model.add(GraphFP(embedding_size, attribute_vector_size-1, 
		depth=depth,
		scale_output=scale_output,
		padding=padding,
		activation_inner='tanh'))
	
	logging.info('cnn_model: added GraphFP layer ({} -> {})'.format('mol', embedding_size))
	if hidden > 0:
		
		model.add(Dense(hidden, activation=hidden_activation))
		logging.info('cnn_model: added {} Dense layer (-> {})'.format(hidden_activation, hidden))
		
	model.add(Dense(output_size, activation=output_activation))
	logging.info('cnn_model: added lin Dense layer (-> {})'.format(output_size))

	# Compile
	if optimizer == 'adam':
		optimizer = Adam(lr=lr)
	elif optimizer == 'rmsprop':
		optimizer = RMSprop(lr=lr)
	else:
		logging.info('Can only handle adam or rmsprop optimizers currently')
		quit(1)

	# Custom loss to filter out NaN values in multi-task predictions
	if loss == 'custom':
		loss = mse_no_NaN

	logging.info('compiling cnn_model...')
	model.compile(loss=loss, optimizer=optimizer)
	logging.info('done compiling.')

	return model