from __future__ import print_function
from .layers import GraphFP
from keras.models import Sequential
from keras.layers.core import Dense
from keras.optimizers import RMSprop, Adam
import numpy as np
import theano.tensor as T 
from keras import initializations
from keras.utils.visualize_util import plot
import json
import datetime
import logging
from rmgpy.cnn_framework.molecule_tensor import get_attribute_vector_size


def build_model(embedding_size=512, attribute_vector_size=None, depth=5, scale_output=0.05, padding=False, 
				hidden=50, hidden_activation='tanh',
				output_activation='linear', output_size=1, 
				lr=0.01, optimizer='adam', loss='mse'):

	"""
	build generic cnn model that takes molecule tensor and predict output 
	with size of output_size.
	"""
	
	model = Sequential()

	if attribute_vector_size is None:
		attribute_vector_size = get_attribute_vector_size()

	model.add(GraphFP(output_dim=embedding_size, 
		inner_dim=attribute_vector_size-1, 
		depth=depth,
		scale_output=scale_output,
		padding=padding,
		activation_inner='tanh'))
	
	logging.info('cnn_model: added GraphFP layer ({} -> {})'.format('mol', embedding_size))
	if hidden > 0:
		
		model.add(Dense(hidden, activation=hidden_activation))
		logging.info('cnn_model: added {} Dense layer (-> {})'.format(hidden_activation, hidden))
		
	model.add(Dense(output_size, activation=output_activation))
	logging.info('cnn_model: added {} Dense layer (-> {})'.format(output_activation, output_size))

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
	(X_train, X_val, X_test, y_train, y_val, y_test) = data

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

		# evaluate test loss upon final model
		mean_test_loss = evaluate_mean_tst_loss(model, X_test, y_test)

	except KeyboardInterrupt:
		logging.info('User terminated training early (intentionally)')

	return (model, loss, val_loss, mean_test_loss)

def evaluate_mean_tst_loss(model, X_test, y_test):

	"""
	Given final model and test examples
	returns mean test loss: a float number
	"""
	test_losses = []
	for j in range(len(X_test)):
		single_mol_as_array = np.array(X_test[j:j+1])
		single_y_as_array = np.reshape(y_test[j], (1, -1))
		sloss = model.test_on_batch(single_mol_as_array, single_y_as_array)
		test_losses.append(sloss)

	mean_test_loss = np.mean(test_losses)
	return mean_test_loss

def reset_model(model):
	'''Given a Keras model consisting only of GraphFP, Dense, and Dropout layers,
	this function will reset the trainable weights to save time for CV tests.'''

	for layer in model.layers:
		# Note: these are custom depending on the layer type
		if '.GraphFP' in str(layer):
			W_inner = layer.init_inner((layer.inner_dim, layer.inner_dim))
			b_inner = np.zeros((1, layer.inner_dim))
			# Inner weights
			layer.W_inner.set_value((T.tile(W_inner, (layer.depth + 1, 1, 1)).eval() + \
				initializations.uniform((layer.depth + 1, layer.inner_dim, layer.inner_dim)).eval()).astype(np.float32))
			layer.b_inner.set_value((T.tile(b_inner, (layer.depth + 1, 1, 1)).eval()  + \
				initializations.uniform((layer.depth + 1, 1, layer.inner_dim)).eval()).astype(np.float32))

			# Outer weights
			W_output = layer.init_output((layer.inner_dim, layer.output_dim), scale = layer.scale_output)
			b_output = np.zeros((1, layer.output_dim))
			# Initialize weights tensor
			layer.W_output.set_value((T.tile(W_output, (layer.depth + 1, 1, 1)).eval()).astype(np.float32))
			layer.b_output.set_value((T.tile(b_output, (layer.depth + 1, 1, 1)).eval()).astype(np.float32))
			logging.info('graphFP layer reset')

		elif '.Dense' in str(layer):
			layer.W.set_value((layer.init(layer.W.shape.eval()).eval()).astype(np.float32))
			layer.b.set_value(np.zeros(layer.b.shape.eval(), dtype=np.float32))
			logging.info('dense layer reset')

		elif '.Dropout' in str(layer):
			logging.info('dropout unchanged')
		else:
			raise ValueError('Unknown layer {}, cannot reset weights'.format(str(layer)))
	logging.info('Reset model weights')
	return model

def save_model(model, loss, val_loss, mean_test_loss, fpath):
	'''Saves NN model object and associated information.

	inputs:
		model - a Keras model
		loss - list of training losses 
		val_loss - list of validation losses
		mean_test_loss - test loss on test set
		fpath - root filepath to save everything to (with .json, h5, png, info 
		config - the configuration dictionary that defined this model 
		tstamp - current timestamp to log in info file'''

	# Dump data
	with open(fpath + '.json', 'w') as structure_fpath:
		json.dump(model.to_json(), structure_fpath)
	logging.info('...saved structural information')

	# Dump weights
	model.save_weights(fpath + '.h5', overwrite = True)
	logging.info('...saved weights')

	# Dump image
	plot(model, to_file = fpath + '.png')
	logging.info('...saved image')

	# Dump history
	save_model_history_manual(loss, val_loss, fpath + '.hist')

	mean_loss = loss[-1]
	mean_val_loss = val_loss[-1]
	write_loss_report(mean_loss, mean_val_loss, mean_test_loss, fpath + '_loss_report.txt')
	logging.info ('...saved history')

	logging.info('...saved model to {}.[json, h5, png]'.format(fpath))

def save_model_history_manual(loss, val_loss, fpath):
	''''This function saves the history returned by model.fit to a tab-
	delimited file, where model is a keras model'''

	# Open file
	fid = open(fpath, 'a')
	logging.info('trained at {}'.format(datetime.datetime.utcnow()))
	print('iteration\tloss\tval_loss', file=fid)

	try:
		# Iterate through
		for i in range(len(loss)):
			print('{}\t{}\t{}'.format(i + 1, 
							loss[i], val_loss[i]),
							file = fid)
	except KeyError:
		print('<no history found>', file = fid)

	# Close file
	fid.close()

def write_loss_report(mean_loss, mean_val_loss, mean_test_loss, fpath):

	"""
	Write training, validation and test mean loss
	"""
	loss_report = open(fpath, 'a')
	print("{:50} {}".format("Training loss (mse):", mean_loss), file=loss_report)
	print("{:50} {}".format("Validation loss (mse):", mean_val_loss), file=loss_report)
	print("{:50} {:.4f}".format("Test loss (mse):", mean_test_loss), file=loss_report)

	# Close file
	loss_report.close()