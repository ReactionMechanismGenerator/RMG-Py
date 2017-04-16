
from cnn_model import build_model, train_model, reset_model, save_model, write_loss_report
from .input import read_input_file
from .molecule_tensor import get_molecule_tensor, pad_molecule_tensor
import os
import rmgpy
import numpy as np
from .data import (prepare_data_one_fold, prepare_folded_data_from_multiple_datasets, 
	prepare_full_train_data_from_multiple_datasets, split_inner_val_from_train_data)
import logging
from keras.callbacks import EarlyStopping
import json

class Predictor(object):

	def __init__(self, input_file=None, datasets_file=None):

		self.model = None
		if input_file:
			self.input_file = input_file
		else:
			self.input_file = os.path.join(os.path.dirname(rmgpy.__file__),
										'cnn_framework',
										'data',  
										'predictor_input.py'
										)

		if datasets_file:
			self.datasets_file = datasets_file
		else:
			self.datasets_file = os.path.join(os.path.dirname(rmgpy.__file__),
										'cnn_framework',
										'data',  
										'datasets.txt'
										)
		self.specify_datasets(self.datasets_file)
	def build_model(self):
		"""
		This method is intended to provide a way to build default model 
		"""

		self.model = build_model()

	def load_input(self, path=None):
		"""
		This method is intended to provide a way to build model from an input file
		"""
		
		if path is None: 
			path = self.input_file
			print path
		read_input_file(path, self)

	def specify_datasets(self, datasets_file_path=None):
		"""
		This method specify which datasets to use for training
		"""
		self.datasets = []
		with open(datasets_file_path, 'r') as f_in:
			for line in f_in:
				line = line.strip()
				if line and not line.startswith('#'):
					db, table = [token.strip() for token in line.split('.')]
					self.datasets.append((db, table))

	def kfcv_train(self, folds, lr_func, save_model_path, batch_size=1):

		# prepare data for training
		folded_data = prepare_folded_data_from_multiple_datasets(self.datasets, folds, 
																self.add_extra_atom_attribute, 
																self.add_extra_bond_attribute,
																self.padding,
																self.padding_final_size)

		X_test, y_test, folded_Xs, folded_ys = folded_data

		losses = []
		inner_val_losses = []
		outer_val_losses = []
		test_losses = []
		for fold in range(folds):
			data = prepare_data_one_fold(folded_Xs, 
										folded_ys, 
										current_fold=fold, 
										shuffle_seed=4)

			# execute train_model
			X_train, X_inner_val, X_outer_val, y_train, y_inner_val, y_outer_val = data
			train_model_output = train_model(self.model, 
											X_train,
											y_train,
											X_inner_val,
											y_inner_val,
											X_test,
											y_test,
											X_outer_val,
											y_outer_val, 
											nb_epoch=150,
											batch_size=batch_size, 
											lr_func=lr_func, 
											patience=10)

			model, loss, inner_val_loss, mean_outer_val_loss, mean_test_loss = train_model_output

			# loss and inner_val_loss each is a list
			# containing loss for each epoch
			losses.append(loss)
			inner_val_losses.append(inner_val_loss)
			outer_val_losses.append(mean_outer_val_loss)
			test_losses.append(mean_test_loss)
			
			# save model and write fold report
			fpath = os.path.join(save_model_path, 'fold_{0}'.format(fold))
			self.save_model(loss, inner_val_loss, mean_outer_val_loss, mean_test_loss, fpath)

			# once finish training one fold, reset the model
			self.reset_model()

		# mean inner_val_loss and outer_val_loss used for selecting parameters, 
		# e.g., lr, epoch, attributes, etc
		full_folds_mean_loss = np.mean([l[-1] for l in losses if len(l) > 0 ])
		full_folds_mean_inner_val_loss = np.mean([l[-1] for v_l in inner_val_losses if len(l) > 0 ])
		full_folds_mean_outer_val_loss = np.mean(outer_val_losses)
		full_folds_mean_test_loss = np.mean(test_losses)

		full_folds_loss_report_path = os.path.join(save_model_path, 'full_folds_loss_report.txt')

		write_loss_report(full_folds_mean_loss, 
						full_folds_mean_inner_val_loss, 
						full_folds_mean_outer_val_loss,
						full_folds_mean_test_loss, 
						full_folds_loss_report_path)

	def full_train(self, lr_func, save_model_path, batch_size=1):

		# prepare data for training
		folded_data = prepare_full_train_data_from_multiple_datasets(self.datasets, 
																self.add_extra_atom_attribute, 
																self.add_extra_bond_attribute,
																self.padding,
																self.padding_final_size)

		X_test, y_test, X_train, y_train = folded_data

		losses = []
		inner_val_losses = []
		test_losses = []
		data = split_inner_val_from_train_data(X_train, y_train)

		X_train, X_inner_val, y_train, y_inner_val = data

		# execute train_model
		logging.info('\nStart full training...')
		logging.info('Training data: {} points'.format(len(X_train)))
		logging.info('Inner val data: {} points'.format(len(X_inner_val)))
		logging.info('Test data: {} points'.format(len(X_test)))
		train_model_output = train_model(self.model, 
										X_train,
										y_train,
										X_inner_val,
										y_inner_val,
										X_test,
										y_test,
										X_outer_val=None,
										y_outer_val=None, 
										nb_epoch=150,
										batch_size=batch_size, 
										lr_func=lr_func, 
										patience=10)

		model, loss, inner_val_loss, mean_outer_val_loss, mean_test_loss = train_model_output

		# loss and inner_val_loss each is a list
		# containing loss for each epoch
		losses.append(loss)
		inner_val_losses.append(inner_val_loss)
		test_losses.append(mean_test_loss)
		
		# save model and write report
		fpath = os.path.join(save_model_path, 'full_train')
		self.save_model(loss, inner_val_loss, mean_outer_val_loss, mean_test_loss, fpath)


	def kfcv_batch_train(self, folds, batch_size=50):

		# prepare data for training
		folded_data = prepare_folded_data_from_multiple_datasets(self.datasets, folds, 
																self.add_extra_atom_attribute, 
																self.add_extra_bond_attribute,
																self.padding,
																self.padding_final_size)

		X_test, y_test, folded_Xs, folded_ys = folded_data

		losses = []
		inner_val_losses = []
		outer_val_losses = []
		test_losses = []
		for fold in range(folds):
			data = prepare_data_one_fold(folded_Xs, 
										folded_ys, 
										current_fold=fold, 
										shuffle_seed=4)

			X_train, X_inner_val, X_outer_val, y_train, y_inner_val, y_outer_val = data

			X_train.extend(X_inner_val)
			y_train.extend(y_inner_val)

			earlyStopping = EarlyStopping(monitor='val_loss', patience=10, verbose=1, mode='auto')

			history_callback = self.model.fit(np.array(X_train), 
							np.array(y_train), 
							callbacks=[earlyStopping], 
							nb_epoch=150, 
							batch_size=batch_size, 
							validation_split=0.1)

			loss_history = history_callback.history
			with open('history.json_fold_{0}'.format(fold), 'w') as f_in:
				json.dump(loss_history, f_in, indent=2)

			# evaluate outer valiation loss
			outer_val_loss = self.model.evaluate(np.array(X_outer_val), 
											np.array(y_outer_val), 
											batch_size=50)
			logging.info("\nOuter val loss: {0}".format(outer_val_loss))

			test_loss = self.model.evaluate(np.array(X_test), np.array(y_test), batch_size=50)
			logging.info("\nTest loss: {0}".format(test_loss))

			# once finish training one fold, reset the model
			self.reset_model()


	def load_parameters(self, param_path=None):

		if not param_path:
			param_path = os.path.join(os.path.dirname(rmgpy.__file__),
									'cnn_framework',
									'data', 
									'weights', 
									'polycyclic_enthalpy_weights.h5'
									)

		self.model.load_weights(param_path)

	def reset_model(self):

		self.model = reset_model(self.model)

	def save_model(self, loss, inner_val_loss, mean_outer_val_loss, mean_test_loss, fpath):

		save_model(self.model, loss, inner_val_loss, mean_outer_val_loss, mean_test_loss, fpath)

	def predict(self, molecule):

		molecule_tensor = get_molecule_tensor(molecule, \
							self.add_extra_atom_attribute, self.add_extra_bond_attribute)
		if self.padding:
			molecule_tensor = pad_molecule_tensor(molecule_tensor, self.padding_final_size)
		molecule_tensor_array = np.array([molecule_tensor])
		return self.model.predict(molecule_tensor_array)[0][0]
