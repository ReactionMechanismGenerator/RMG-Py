import os

from rmgpy.molecule import Molecule
from rmgpy.cnn_framework.predictor import Predictor

def load_cnn_thermo_predictor(thermo_predictor_model='zb_biycyclics_t.r.0.5'):

	thermo_predictor_data = os.path.join(os.path.dirname(__file__),
										'cnn_model_selection')

	datasets_file = os.path.join(thermo_predictor_data,
								thermo_predictor_model,
								'datasets.txt')

	predictor_input = os.path.join(thermo_predictor_data,
								thermo_predictor_model,
								'predictor_input.py')

	h298_predictor = Predictor(datasets_file=datasets_file)
	h298_predictor.load_input(predictor_input)
	param_path = os.path.join(thermo_predictor_data,
							thermo_predictor_model,
							'saved_model',
							'full_train.h5')
	h298_predictor.load_parameters(param_path)

	return h298_predictor

molecule_smiles = 'C1CCCCC1'
h298_predictor = load_cnn_thermo_predictor()
mol = Molecule(SMILES=molecule_smiles)
thermo_result = h298_predictor.predict(mol)
print thermo_result