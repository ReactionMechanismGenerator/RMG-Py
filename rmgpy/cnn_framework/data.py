
import logging
import numpy as np
from rmgpy.molecule.molecule import Molecule
from .molecule_tensor import get_molecule_tensor, pad_molecule_tensor
from pymongo import MongoClient

def get_data_from_db(db_name, collection_name, 
					add_extra_atom_attribute=True, 
					add_extra_bond_attribute=True,
					padding=True,
					padding_final_size=20):

	# connect to db and query
	client = MongoClient('mongodb://user:user@rmg.mit.edu/admin', 27018)
	db =  getattr(client, db_name)
	collection = getattr(db, collection_name)
	db_cursor = collection.find()

	# collect data
	logging.info('Collecting polycyclic data: {0}.{1}...'.format(db_name, collection_name))
	X = []
	y = []
	db_mols = []
	for db_mol in db_cursor:
		db_mols.append(db_mol)

	logging.info('Done collecting data: {0} points...'.format(len(db_mols)))

	logging.info('Generating molecular tensor data...')

	for db_mol in db_mols:
		smile = str(db_mol["SMILES_input"])
		mol = Molecule().fromSMILES(smile)
		mol_tensor = get_molecule_tensor(mol, add_extra_atom_attribute, 
										 add_extra_bond_attribute)
		if padding:
			mol_tensor = pad_molecule_tensor(mol_tensor, padding_final_size)
		hf298_qm = float(db_mol["Hf298"])
		X.append(mol_tensor)
		y.append(hf298_qm)

	logging.info('Done: generated {0} tensors with padding={1}'.format(len(X), padding))

	return (X, y)

def prepare_folded_data_from_multiple_datasets(datasets, 
												folds, 
												add_extra_atom_attribute, 
												add_extra_bond_attribute,
												padding=True,
												padding_final_size=20):

	folded_datasets = []
	test_data_datasets = []
	for db, table in datasets:
		(X, y) = get_data_from_db(db, 
								table, 
								add_extra_atom_attribute, 
								add_extra_bond_attribute,
								padding,
								padding_final_size)

		(X_test, y_test, X_train_and_val, y_train_and_val) = split_tst_from_train_and_val(X, y, shuffle_seed=0)
		test_data_datasets.append((X_test, y_test))

		(folded_Xs, folded_ys) = prepare_folded_data(X_train_and_val, y_train_and_val, folds, shuffle_seed=2)
		folded_datasets.append((folded_Xs, folded_ys))

	# merge into one folded_Xs and folded_ys
	logging.info('Merging {} datasets for training...'.format(len(datasets)))
	(folded_Xs, folded_ys) = folded_datasets[0]
	if len(folded_datasets) > 1:
		for folded_Xs_1, folded_ys_1 in folded_datasets[1:]:
			folded_Xs_ext = []
			folded_ys_ext = []
			for idx, folded_X in enumerate(folded_Xs):
				folded_X.extend(folded_Xs_1[idx])
				folded_Xs_ext.append(folded_X)

				folded_y = folded_ys[idx]
				folded_y.extend(folded_ys_1[idx])
				folded_ys_ext.append(folded_y)

			folded_Xs = folded_Xs_ext
			folded_ys = folded_ys_ext

	# merge into one X_test and y_test
	(X_test, y_test) = test_data_datasets[0]
	if len(test_data_datasets) > 1:
		for X_test_1, y_test_1 in test_data_datasets[1:]:
			X_test.extend(X_test_1)
			y_test.extend(y_test_1)

	return X_test, y_test, folded_Xs, folded_ys
def split_tst_from_train_and_val(X, y, shuffle_seed=None, testing_ratio=0.1):

	n = len(X)
	# Feed shuffle seed
	if shuffle_seed is not None:
		np.random.seed(shuffle_seed)

	all_indices = range(n)
	np.random.shuffle(all_indices)

	# shuffle X and y
	X_shuffled = [X[i] for i in all_indices]
	y_shuffled = [y[i] for i in all_indices]

	split = int(len(all_indices) * testing_ratio)
	X_test, X_train_and_val = [X_shuffled[i] for i in all_indices[:split]],   [X_shuffled[i] for i in all_indices[split:]]
	y_test, y_train_and_val = [y_shuffled[i] for i in all_indices[:split]],   [y_shuffled[i] for i in all_indices[split:]]

	# reset np random seed to avoid side-effect on other methods
	# relying on np.random
	if shuffle_seed is not None:
		np.random.seed()

	return (X_test, y_test, X_train_and_val, y_train_and_val)

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
	this method prepares X_train, y_train, X_inner_val, y_inner_val,
	and X_outer_val, y_outer_val
	"""
	logging.info('...using fold {}'.format(current_fold+1))

	# Recombine into training and testing
	X_train   = [x for folded_X in (folded_Xs[:current_fold] + folded_Xs[(current_fold + 1):])  for x in folded_X]
	y_train   = [y for folded_y in (folded_ys[:current_fold] + folded_ys[(current_fold + 1):])  for y in folded_y]

	# Test is current_fold
	X_outer_val    = folded_Xs[current_fold]
	y_outer_val    = folded_ys[current_fold]

	# Define validation set as random 10% of training
	if shuffle_seed is not None:
		np.random.seed(shuffle_seed)

	training_indices = range(len(X_train))
	np.random.shuffle(training_indices)
	split = int(len(training_indices) * training_ratio)
	X_train,   X_inner_val  = [X_train[i] for i in training_indices[:split]],   [X_train[i] for i in training_indices[split:]]
	y_train,   y_inner_val  = [y_train[i] for i in training_indices[:split]],   [y_train[i] for i in training_indices[split:]]

	logging.info('Total training: {}'.format(len(X_train)))
	logging.info('Total inner validation: {}'.format(len(X_inner_val)))
	logging.info('Total outer validation: {}'.format(len(X_outer_val)))

	# reset np random seed to avoid side-effect on other methods
	# relying on np.random
	if shuffle_seed is not None:
		np.random.seed()

	return (X_train, X_inner_val, X_outer_val, y_train, y_inner_val, y_outer_val)