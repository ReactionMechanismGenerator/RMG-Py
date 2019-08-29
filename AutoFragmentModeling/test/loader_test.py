import os
import unittest

import afm.loader

class TestLoader(unittest.TestCase):

	def test_load_fragment_reactions_from_chemkin(self):

		chemkin_path = os.path.join(os.path.dirname(__file__), 
									'data', 
									'loader_data',
									'chem.inp')

		dictionary_path = os.path.join(os.path.dirname(__file__), 
									'data', 
									'loader_data',
									'species_dictionary.txt')

		fragment_smiles_path = os.path.join(os.path.dirname(__file__), 
									'data', 
									'loader_data',
									'fragment_smiles.txt')

		fragments_dict, fragment_rxns = afm.loader.load_fragment_reactions_from_chemkin(chemkin_path,
                                        												dictionary_path,
                                        												fragment_smiles_path)

		self.assertEqual(40, len(fragments_dict))
		self.assertEqual(312, len(fragment_rxns))

	def test_load_pseudo_fragment_reactions(self):

		chemkin_path = os.path.join(os.path.dirname(__file__), 
									'data', 
									'loader_data',
									'chem.inp')

		dictionary_path = os.path.join(os.path.dirname(__file__), 
									'data', 
									'loader_data',
									'species_dictionary.txt')

		fragment_smiles_path = os.path.join(os.path.dirname(__file__), 
									'data', 
									'loader_data',
									'fragment_smiles.txt')

		fragments_dict, _ = afm.loader.load_fragment_reactions_from_chemkin(chemkin_path,
                                        									dictionary_path,
                                        									fragment_smiles_path)

		pseudo_fragrxns = afm.loader.load_pseudo_fragment_reactions(fragments_dict)

		# currently only one reaction is added as
		# pseudo reaction
		self.assertEqual(1, len(pseudo_fragrxns))

		pseudo_fragrxn = pseudo_fragrxns[0]

		self.assertEqual('pseudo_rxn', pseudo_fragrxn.family)
