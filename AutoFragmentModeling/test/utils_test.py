import unittest

import afm.utils

class TestUtils(unittest.TestCase):

	def test_match_sequences1(self):

		seq1 = [1, 3, 1]
		seq2 = [2, 1, 2]

		matches = afm.utils.match_sequences(seq1, seq2)

		expected_matches = [[(0,0),1], 
							[(1,0),1], 
							[(1,1),1],
							[(1,2),1],
							[(2,2),1]]

		self.assertEqual(matches, expected_matches)


	def test_match_sequences2(self):
		"""
		Test match_sequences() can tolerate slight
		sum difference between two sequences, default tolerance
		is 1e-6.
		"""

		seq1 = [1, 3, 1-1e-6]
		seq2 = [2, 1, 2]

		matches = afm.utils.match_sequences(seq1, seq2)

		expected_matches = [[(0,0),1], 
							[(1,0),1], 
							[(1,1),1],
							[(1,2),1],
							[(2,2),1]]

		self.assertEqual(matches, expected_matches)

	def test_match_concentrations_with_same_sums(self):

		conc1 = [('a', 1),
				('b', 3),
				('c', 1)]

		conc2 = [('x', 2),
				('y', 1),
				('z', 2)]

		matches = afm.utils.match_concentrations_with_same_sums(conc1, conc2)

		expected_matches = [(('a','x'),1), 
							 (('b','x'),1), 
							 (('b','y'),1),
							 (('b','z'),1),
							 (('c','z'),1)]

		self.assertEqual(matches, expected_matches)

	def test_match_concentrations_with_different_sums1(self):

		conc1 = [('a', 1),
				('b', 3),
				('c', 1)]
		conc2 = [('x', 2),
				('y', 1),
				('z', 10)]

		matches = afm.utils.match_concentrations_with_different_sums(conc1, conc2)

		expected_matches = [(((('a', 'x'), 'z'), 'z'), 1), 
							(((('b', 'x'), 'z'), 'z'), 1), 
							(((('b', 'y'), 'z'), 'z'), 1),
							((('b','z'), 'z'), 1),
							((('c','z'), 'z'), 1)]
		
		self.assertEqual(matches, expected_matches)

	def test_match_concentrations_with_different_sums2(self):

		conc1 = [(('LY', 'XR'), 10),
			   	 (('XR', 'LWL', 'XR'), 2),
			   	 (('LY', 'TR'), 3)]
		conc2 = [(('LWL', 'RUR'), 3),
				 ('LQR', 3)]

		matches = afm.utils.match_concentrations_with_different_sums(conc1, conc2)

		expected_matches = [((('LY', 'XR'), ('LWL', 'RUR')), 3),
							((('LY', 'XR'), 'LQR'), 3),
							(('LY', 'XR'), 4),
							(('XR', 'LWL', 'XR'), 2),
							(('LY', 'TR'), 3)]
		
		self.assertEqual(matches, expected_matches)

	def test_match_concentrations_with_different_sums3(self):

		conc1 = [(('LY', 'XR'), 10),
			   	 (('XR', 'LWL', 'XR'), 2),
			   	 (('LY', 'TR'), 3)]
		conc2 = [(('LWL', 'RUR'), 3),
				 ('LQR', 7)]

		matches = afm.utils.match_concentrations_with_different_sums(conc1, conc2)

		expected_matches = [((('LY', 'XR'), ('LWL', 'RUR')), 3),
							((('LY', 'XR'), 'LQR'), 7),
							(('XR', 'LWL', 'XR'), 2),
							(('LY', 'TR'), 3)]
		
		self.assertEqual(matches, expected_matches)

	def test_match_concentrations_with_different_sums4(self):

		conc1 = [(('LY', 'XR'), 10),
			   	 (('XR', 'LWL', 'XR'), 2),
			   	 (('LY', 'TR'), 3)]
		conc2 = [(('LWL', 'RUR'), 3),
				 ('LQR', 10)]

		matches = afm.utils.match_concentrations_with_different_sums(conc1, conc2)

		expected_matches = [((('LY', 'XR'), ('LWL', 'RUR')), 3),
							((('LY', 'XR'), 'LQR'), 7),
							((('XR', 'LWL', 'XR'), 'LQR'), 2),
							((('LY', 'TR'), 'LQR'), 1),
							(('LY', 'TR'), 2)]
		
		self.assertEqual(matches, expected_matches)

	def test_match_concentrations_with_different_sums5(self):

		conc1 = [(('LY', 'XR'), 10),
			   	 (('XR', 'LWL', 'XR'), 2),
			   	 (('LY', 'TR'), 3)]
		conc2 = [(('LWL', 'RUR'), 3),
				 ('LQR', 12)]

		matches = afm.utils.match_concentrations_with_different_sums(conc1, conc2)

		expected_matches = [((('LY', 'XR'), ('LWL', 'RUR')), 3),
							((('LY', 'XR'), 'LQR'), 7),
							((('XR', 'LWL', 'XR'), 'LQR'), 2),
							((('LY', 'TR'), 'LQR'), 3)]
		
		self.assertEqual(matches, expected_matches)

	def test_match_concentrations_with_different_sums6(self):

		conc1 = [(('LY', 'XR'), 10),
			   	 (('XR', 'LWL', 'XR'), 2),
			   	 (('LY', 'TR'), 3)]
		conc2 = [(('LWL', 'RUR'), 3),
				 ('LQR', 30)]

		matches = afm.utils.match_concentrations_with_different_sums(conc1, conc2)

		expected_matches = [((((('LY', 'XR'), ('LWL', 'RUR')), 'LQR'), 'LQR'), 3),
							(((('LY', 'XR'), 'LQR'), 'LQR'), 7),
							(((('XR', 'LWL', 'XR'), 'LQR'), 'LQR'), 2),
							(((('LY', 'TR'), 'LQR'), 'LQR'), 3)]
		
		self.assertEqual(matches, expected_matches)

	def test_grind(self):

		conc = [('a', 1),
				('b', 3),
				('c', 1)]
		size = 0.6

		grinded_conc = afm.utils.grind(conc, size)
		expected_grinded_conc = [('a', 0.6),
								 ('a', 0.4),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('c', 0.6),
								 ('c', 0.4)]
		self.assertEqual(grinded_conc, expected_grinded_conc)

	def test_shuffle(self):

		conc = [('a', 1),
				('b', 3),
				('c', 1),
				('d', 2)]
		
		seed = 0
		shuffled_conc = afm.utils.shuffle(conc, seed)

		expected_conc = [('c', 1),
						 ('d', 2),
						 ('b', 3),
						 ('a', 1)]		 

		self.assertEqual(shuffled_conc, expected_conc)

	def test_matches_resolve(self):

		matches = [(('LY', 'XR'), 10),
				   (('LWL', 'XR'), 4),
				   (('LWL', 'RUR'), 6),
				   (('LY', 'TR'), 3)]
		rr_ll_list = ['LWL', 'RUR']

		new_matches, new_r_l_moles = afm.utils.matches_resolve(matches, rr_ll_list)

		expected_new_matches = [(('LY', 'XR'), 10),
							    (('XR', 'LWL', 'XR'), 2),
							    (('LY', 'TR'), 3)]
		expected_new_r_l_moles = [(('LWL', 'RUR'), 3)]

		self.assertEqual(new_matches, expected_new_matches)
		self.assertEqual(new_r_l_moles, expected_new_r_l_moles)

# 1 label (R label only)
# test matches resolve for 1 label


	def test_matches_resolve_1_label(self):

		matches = [(('RY', 'CCR'), 10),
				   (('ArCR', 'RCCCCR'), 4),
				   (('C_CCR', 'RCC_CCR'), 6),
				   (('RCCCCR', 'RCC_CCR'), 6),
				   (('RY', 'TR'), 3)]
		rr_list = ['RCCCCR', 'RCC_CCR']

		new_matches_1_label, new_r_l_moles_1_label = afm.utils.matches_resolve(matches, rr_list)

		expected_new_matches_1_label = [(('RY', 'CCR'), 10),
							    (('ArCR', 'RCCCCR', 'ArCR'), 2),
				   			    (('C_CCR', 'RCC_CCR', 'C_CCR'), 3),
							    (('RY', 'TR'), 3)]
		expected_new_r_l_moles_1_label = [(('RCCCCR', 'RCC_CCR'), 3)]

		self.assertEqual(new_matches_1_label, expected_new_matches_1_label)
		self.assertEqual(new_r_l_moles_1_label, expected_new_r_l_moles_1_label)
