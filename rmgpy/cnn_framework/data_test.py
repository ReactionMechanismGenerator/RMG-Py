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
		(folded_Xs, folded_ys) = prepare_folded_data(self.X, self.y, folds)
		self.assertEqual(len(folded_Xs), folds)
		self.assertEqual(len(folded_ys), folds)

	def test_prepare_data_one_fold(self):

		folds = 5
		training_ratio=0.9
		(folded_Xs, folded_ys) = prepare_folded_data(self.X, self.y, folds)
		data = prepare_data_one_fold(folded_Xs, folded_ys, current_fold=0, training_ratio=training_ratio)

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
		

