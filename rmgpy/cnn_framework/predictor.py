
from cnn_model import build_model, train_model, reset_model, save_model, write_loss_report
from .input import read_input_file
from .molecule_tensor import get_molecule_tensor
import os
import rmgpy
import numpy as np
from .data import get_data_from_db, prepare_folded_data, prepare_data_one_fold
import logging

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

	def kfcv_train(self, dataset, folds, lr_func, save_model_path):

		# prepare data for training
		[db, data_table] = dataset.split('.')
		(X, y) = get_data_from_db(db, data_table, self.add_extra_atom_attribute, self.add_extra_bond_attribute)
		(folded_Xs, folded_ys) = prepare_folded_data(X, y, folds, shuffle_seed=0)

		losses = []
		val_losses = []
		test_losses = []
		for fold in range(folds):
			data = prepare_data_one_fold(folded_Xs, folded_ys, current_fold=fold, shuffle_seed=4)

			# execute train_model
			(model, loss, val_loss, mean_test_loss) = train_model(self.model, data, 
												nb_epoch=150, lr_func=lr_func, 
												patience=10)

			# loss and val_loss each is a list
			# containing loss for each epoch
			losses.append(loss)
			val_losses.append(val_loss)
			test_losses.append(mean_test_loss)
			
			# save model and write fold report
			fpath = os.path.join(save_model_path, 'fold_{0}'.format(fold))
			self.save_model(loss, val_loss, mean_test_loss, fpath)

			# once finish training one fold, reset the model
			self.reset_model()

		# mean loss and val_loss used for selecting parameters, 
		# e.g., lr, epoch, attributes, etc
		full_folds_mean_loss = np.mean([l[-1] for l in losses if len(l) > 0 ])
		full_folds_mean_val_loss = np.mean([l[-1] for v_l in val_losses if len(l) > 0 ])
		full_folds_mean_test_loss = np.mean(test_losses)

		full_folds_loss_report_path = os.path.join(save_model_path, 'full_folds_loss_report.txt')
		write_loss_report(full_folds_mean_loss, full_folds_mean_val_loss, \
			full_folds_mean_test_loss, full_folds_loss_report_path)

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

	def save_model(self, loss, val_loss, mean_test_loss, fpath):

		save_model(self.model, loss, val_loss, mean_test_loss, fpath)

	def predict(self, molecule):

		molecule_tensor = get_molecule_tensor(molecule, \
							self.add_extra_atom_attribute, self.add_extra_bond_attribute)

		molecule_tensor_array = np.array([molecule_tensor])
		return self.model.predict(molecule_tensor_array)[0][0]
