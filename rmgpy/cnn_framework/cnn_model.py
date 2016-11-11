from .layers import GraphFP
from keras.models import Sequential
from keras.layers.core import Dense
from keras.optimizers import RMSprop, Adam
import numpy as np

import logging

def build_model(embedding_size=512, attribute_vector_size=33, depth=5, scale_output=0.05, padding=False, 
				hidden=50, hidden_activation='tanh',
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

def train_model(model, data, nb_epoch=0, lr_func='0.01', patience=10):
	"""
	inputs:
		model - a Keras model
		data - X_train, X_val, X_test, y_train, y_val, y_test
		nb_epoch - number of epochs to train for
		lr_func - string which is evaluated with 'epoch' to produce the learning 
				rate at each epoch 
		patience - number of epochs to wait when no progress is being made in 
				the validation loss

	outputs:
		model - a trained Keras model
		loss - list of training losses corresponding to each epoch 
		val_loss - list of validation losses corresponding to each epoch
	"""

	# Get data from helper function
	(X_train, X_val, _, y_train, y_val, _) = data

	# Create learning rate function
	lr_func_string = 'def lr(epoch):\n    return {}\n'.format(lr_func)
	exec lr_func_string

	# Fit (allows keyboard interrupts in the middle)
	# Because molecular graph tensors are different sizes based on N_atoms, can only do one at a time
	# (alternative is to pad with zeros and try to add some masking feature to GraphFP)
	try:
		loss = []
		val_loss = []

		wait = 0
		prev_best_val_loss = 99999999
		for i in range(nb_epoch):
			logging.info('Epoch {}/{}, lr = {}'.format(i + 1, nb_epoch, lr(i)))
			this_loss = []
			this_val_loss = []
			model.optimizer.lr.set_value(lr(i))
			
			# Run through training set
			logging.info('Training...')
			training_order = range(len(X_train))
			np.random.shuffle(training_order)
			for j in training_order:
				single_mol_as_array = np.array(X_train[j:j+1])
				single_y_as_array = np.reshape(y_train[j], (1, -1))
				sloss = model.train_on_batch(single_mol_as_array, single_y_as_array)
				this_loss.append(sloss)
			
			# Run through testing set
			logging.info('Validating..')
			for j in range(len(X_val)):
				single_mol_as_array = np.array(X_val[j:j+1])
				single_y_as_array = np.reshape(y_val[j], (1, -1))
				sloss = model.test_on_batch(single_mol_as_array, single_y_as_array)
				this_val_loss.append(sloss)
			
			loss.append(np.mean(this_loss))
			val_loss.append(np.mean(this_val_loss))
			logging.info('mse loss: {}\tmse val_loss: {}'.format(loss[i], val_loss[i]))

			# Check progress
			if np.mean(this_val_loss) < prev_best_val_loss:
				wait = 0
				prev_best_val_loss = np.mean(this_val_loss)
				if patience == -1:
					model.save_weights('train_cnn_results/best.h5', overwrite=True)
			else:
				wait = wait + 1
				logging.info('{} epochs without val_loss progress'.format(wait))
				if wait == patience:
					logging.info('stopping early!')
					break
		if patience == -1:
			model.load_weights('train_cnn_results/best.h5')

	except KeyboardInterrupt:
		logging.info('User terminated training early (intentionally)')

	return (model, loss, val_loss)