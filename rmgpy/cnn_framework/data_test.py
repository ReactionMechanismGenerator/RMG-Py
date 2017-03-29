import unittest
from .data import *

class Test_Data(unittest.TestCase):

	def setUp(self):

		(X, y) = get_data_from_db('sdata134k', 'polycyclic_2954_table')
		
		self.X = X
		self.y = y


	def test_get_HC_polycyclics_data_from_db(self):

		self.assertEqual(len(self.X), 2954)
		self.assertEqual(len(self.y), 2954)

	def test_prepare_folded_data(self):

		folds = 5
		shuffle_seed = 4
		(folded_Xs, folded_ys) = prepare_folded_data(self.X, self.y, folds, shuffle_seed)
		self.assertEqual(len(folded_Xs), folds)
		self.assertEqual(len(folded_ys), folds)

		# test shuffle is expected
		import numpy as np
		all_indices = range(len(self.X))
		np.random.seed(shuffle_seed)
		np.random.shuffle(all_indices)

		first_X = folded_Xs[0][0]
		expected_first_X = self.X[all_indices[0]]

		self.assertTrue(np.all(np.equal(first_X, expected_first_X)))

	def test_prepare_data_one_fold(self):

		import numpy as np

		folds = 5
		training_ratio=0.9

		n = len(self.X)
		target_fold_size = int(np.ceil(float(n) / folds))
		folded_Xs 		= [self.X[i:i+target_fold_size]   for i in range(0, n, target_fold_size)]
		folded_ys 		= [self.y[i:i+target_fold_size]   for i in range(0, n, target_fold_size)]

		shuffle_seed = 4 # seed for method `prepare_data_one_fold()`
		data = prepare_data_one_fold(folded_Xs, folded_ys, current_fold=0, training_ratio=training_ratio,\
									shuffle_seed=shuffle_seed)

		self.assertEqual(len(data), 6)

		X_train = data[0]
		X_val = data[1]
		X_test = data[2]

		self.assertAlmostEqual(len(X_train)/10.0, 
							training_ratio*int(np.ceil(1.0*len(self.X)/folds))*(folds - 1)/10.0, 
							0)
		self.assertAlmostEqual(len(X_val)/10.0, 
							(1-training_ratio)*int(np.ceil(1.0*len(self.X)/folds))*(folds - 1)/10.0, 
							0)
		self.assertAlmostEqual(len(X_test)/10.0, 
							int(np.ceil(1.0*len(self.X)/folds))/10.0, 
							0)


		# test shuffle is expected
		testset_size = len(folded_Xs[0])
		training_val_indices = range(testset_size, len(self.X))
		np.random.seed(shuffle_seed)
		np.random.shuffle(training_val_indices)

		first_X_in_train = X_train[0]
		expected_first_X_in_train = self.X[training_val_indices[0]]

		self.assertTrue(np.all(np.equal(first_X_in_train, expected_first_X_in_train)))

	def test_prepare_folded_data_from_multiple_datasets(self):

		datasets =  [
					('sdata134k', 'kh_tricyclic_table'), 
					('sdata134k', 'kh_tricyclic_table'),
					('sdata134k', 'kh_tricyclic_table')
					]


		X_test, y_test, folded_Xs, folded_ys = prepare_folded_data_from_multiple_datasets(
										datasets=datasets, 
										folds=5, add_extra_atom_attribute=True, 
										add_extra_bond_attribute=True)
		self.assertEqual(len(folded_Xs), 5)
		self.assertEqual(len(folded_ys), 5)

		self.assertEqual(len(X_test), 54)
		self.assertEqual(len(y_test), 54)
		self.assertEqual(len(folded_Xs[0]), 99)
		self.assertEqual(len(folded_Xs[0]), 99)

	def test_split_tst_from_train_and_val(self):

		X_test, y_test, X_train_and_val, y_train_and_val = split_tst_from_train_and_val(\
															self.X, self.y, 
															shuffle_seed=None, 
															testing_ratio=0.1)

		self.assertEqual(len(X_test), 295)
		self.assertEqual(len(y_test), 295)
		self.assertEqual(len(X_train_and_val), 2659)
		self.assertEqual(len(y_train_and_val), 2659)