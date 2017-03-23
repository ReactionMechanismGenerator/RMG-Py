
import logging
import numpy as np
from rmgpy.molecule.molecule import Molecule
from .molecule_tensor import get_molecule_tensor
from pymongo import MongoClient

def get_data_from_db(db_name, collection_name, 
					add_extra_atom_attribute=True, add_extra_bond_attribute=True):

	# connect to db and query
	client = MongoClient('mongodb://user:user@rmg.mit.edu/admin', 27018)
	db =  getattr(client, db_name)
	collection = getattr(db, collection_name)
	db_cursor = collection.find()

	# collect data
	logging.info('Collecting polycyclic data: {0}.{1}...'.format(db_name, collection_name))
	X = []
	y = []
	for db_mol in db_cursor:
		smile = str(db_mol["SMILES_input"])
		mol = Molecule().fromSMILES(smile)
		mol_tensor = get_molecule_tensor(mol, add_extra_atom_attribute, 
										 add_extra_bond_attribute)
		hf298_qm = float(db_mol["Hf298"])
		X.append(mol_tensor)
		y.append(hf298_qm)

	logging.info('Done collecting data: {0} points...'.format(len(X)))
	
	return (X, y)

def prepare_folded_data(X, y, folds, shuffle_seed=None):

	# Get target size of each fold
	n = len(X)
	logging.info('Total of {} input data points'.format(n))
	target_fold_size = int(np.ceil(float(n) / folds))

	# Feed shuffle seed
	if shuffle_seed is not None:
		np.random.seed(shuffle_seed)

	all_indices = range(n)
	np.random.shuffle(all_indices)

	# shuffle X and y
	X_shuffled = [X[i] for i in all_indices]
	y_shuffled = [y[i] for i in all_indices]

	# Split up data
	folded_Xs 		= [X_shuffled[i:i+target_fold_size]   for i in range(0, n, target_fold_size)]
	folded_ys 		= [y_shuffled[i:i+target_fold_size]   for i in range(0, n, target_fold_size)]

	logging.info('Split data into {} folds'.format(folds))

	# reset np random seed to avoid side-effect on other methods
	# relying on np.random
	if shuffle_seed is not None:
		np.random.seed()

	return (folded_Xs, folded_ys)


def prepare_data_one_fold(folded_Xs, folded_ys, current_fold=0, 
						shuffle_seed=None, training_ratio=0.9):

	"""
	this method prepares X_train, y_train, X_val, y_val,
	and X_test, y_test
	"""
	logging.info('...using fold {}'.format(current_fold+1))

	# Recombine into training and testing
	X_train   = [x for folded_X in (folded_Xs[:current_fold] + folded_Xs[(current_fold + 1):])  for x in folded_X]
	y_train   = [y for folded_y in (folded_ys[:current_fold] + folded_ys[(current_fold + 1):])  for y in folded_y]

	# Test is current_fold
	X_test    = folded_Xs[current_fold]
	y_test    = folded_ys[current_fold]

	# Define validation set as random 10% of training
	if shuffle_seed is not None:
		np.random.seed(shuffle_seed)

	training_indices = range(len(X_train))
	np.random.shuffle(training_indices)
	split = int(len(training_indices) * training_ratio)
	X_train,   X_val  = [X_train[i] for i in training_indices[:split]],   [X_train[i] for i in training_indices[split:]]
	y_train,   y_val  = [y_train[i] for i in training_indices[:split]],   [y_train[i] for i in training_indices[split:]]

	logging.info('Total training: {}'.format(len(X_train)))
	logging.info('Total validation: {}'.format(len(X_val)))
	logging.info('Total testing: {}'.format(len(X_test)))

	# reset np random seed to avoid side-effect on other methods
	# relying on np.random
	if shuffle_seed is not None:
		np.random.seed()

	return (X_train, X_val, X_test, y_train, y_val, y_test)