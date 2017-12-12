
import os
import argparse
import pandas as pd
from tqdm import tqdm

from rmgpy.cnn_framework.data import get_db_mols
from rmgpy.molecule.molecule import Molecule
from rmgpy.cnn_framework.predictor import Predictor

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to RMG Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--datasets', metavar='FILE', type=str, 
        help='path to da file specifies on which datasets to test')

    parser.add_argument('-m', '--model', type=str, 
        help='path to the testing model')


    return parser.parse_args()
################################################################################

def read_datasets_file(datasets_file_path):
    """
    This method specify which datasets to use for validation
    """
    datasets = []
    with open(datasets_file_path, 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if line and not line.startswith('#'):
                host, db, table = [token.strip() for token in line.split('.')]
                datasets.append((host, db, table))

    return datasets

def prepare_data(host, db_name, collection_name):

    # load validation data
    db_mols = get_db_mols(host, db_name, collection_name)

    smiles_list = []
    hf298s_qm = []
    for i, db_mol in enumerate(db_mols):
        smiles = str(db_mol["SMILES_input"])
        hf298_qm = float(db_mol["Hf298(kcal/mol)"])
        
        smiles_list.append(smiles)
        hf298s_qm.append(hf298_qm)

    return smiles_list, hf298s_qm

def prepare_predictor(model):

    h298_predictor = Predictor()

    predictor_input = os.path.join(model,
                                  'predictor_input.py')

    h298_predictor.load_input(predictor_input)

    param_path = os.path.join(model,
                             'saved_model',
                             'full_train.h5')
    h298_predictor.load_parameters(param_path)

    return h298_predictor

def make_predictions(predictor, smiles_list):

    hf298s_cnn = []
    for smiles in tqdm(smiles_list):
        mol = Molecule().fromSMILES(smiles)
        hf298_cnn = predictor.predict(mol)
        hf298s_cnn.append(hf298_cnn)

    return hf298s_cnn

def evaluate(smiles_list, hf298s_qm, hf298s_cnn):

    hf298_df = pd.DataFrame(index=smiles_list)

    hf298_df['Hf298_qm(kcal/mol)'] = pd.Series(hf298s_qm, index=hf298_df.index)
    hf298_df['Hf298_cnn(kcal/mol)'] = pd.Series(hf298s_cnn, index=hf298_df.index)

    qm_cnn_diff = abs(hf298_df['Hf298_qm(kcal/mol)']-hf298_df['Hf298_cnn(kcal/mol)'])
    hf298_df['H298_qm_cnn_diff(kcal/mol)'] = pd.Series(qm_cnn_diff, index=hf298_df.index)

    return hf298_df

def display_result(hf298_df):

    descr = hf298_df['H298_qm_cnn_diff(kcal/mol)'].describe()

    count = int(descr.loc['count'])
    mean = descr.loc['mean']
    std = descr.loc['std']

    display_str = 'count: {0}, error mean: {1:.02f} kcal/mol, error std: {2:.02f} kcal/mol'.format(count, mean, std) 
    print display_str 

def validate(datasets_file, model):

    # load cnn predictor
    h298_predictor = prepare_predictor(model)

    datasets = read_datasets_file(datasets_file)
    for host, db_name, collection_name in datasets:
        
        print "\nhost: {0}, db: {1}, collection: {2}".format(host, db_name, collection_name)
        
        # prepare data for testing
        smiles_list, hf298s_qm = prepare_data(host, db_name, collection_name)

        # Do the predictions
        hf298s_cnn = make_predictions(h298_predictor, smiles_list)

        # evaluate performance
        hf298_df = evaluate(smiles_list, hf298s_qm, hf298s_cnn)

        # display result
        display_result(hf298_df)


def main():

    args = parseCommandLineArguments()

    datasets_file = args.datasets
    model = args.model
    validate(datasets_file, model)

main()