import os
import unittest

from rmgpy import cnn_framework
from rmgpy.molecule import Molecule
from rmgpy.cnn_framework.main import *


class TestMCNNEstimator(unittest.TestCase):

	def setUp(self):

		pretrained_models_path = os.path.join(cnn_framework.__path__[0], 
											  'pretrained_models')
		Hf298_path = os.path.join(pretrained_models_path, 'Hf298')
		S298_path = os.path.join(pretrained_models_path, 'S298')
		Cp_path = os.path.join(pretrained_models_path, 'Cp')
		self.mcnn_estimator = MCNNEstimator(Hf298_path, S298_path, Cp_path)

	def test_get_thermo_data(self):

		mol = Molecule().fromSMILES('C1C2C1C2')
		thermo = self.mcnn_estimator.get_thermo_data(mol)

		self.assertTrue(thermo.comment.startswith('MCNN Estimation.'))
		self.assertAlmostEqual(thermo.Cp0.value_si, 33.3, 1)
		self.assertAlmostEqual(thermo.CpInf.value_si, 232.8, 1)
		self.assertEqual(len(thermo.Cpdata.value_si), 7)
