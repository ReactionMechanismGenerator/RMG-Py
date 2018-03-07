
import os
import argparse
import json
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

def prepare_data(host, db_name, collection_name, prediction_task="Hf298(kcal/mol)"):

    # load validation data
    db_mols = get_db_mols(host, db_name, collection_name)

    smiles_list = []
    ys = []
    # decide what predict task is
    if prediction_task not in ["Hf298(kcal/mol)", "S298(cal/mol/K)", "Cp"]:
        raise NotImplementedError("Prediction task: {0} not supported yet!".format(prediction_task))

    for i, db_mol in enumerate(db_mols):
        smiles = str(db_mol["SMILES_input"])

        if prediction_task != "Cp":
            y = float(db_mol[prediction_task])
        else:
            try:
                y = float(db_mol["Cp298(cal/mol/K)"])
            except KeyError:
                y = float(db_mol["Cp300(cal/mol/K)"])
        
        smiles_list.append(smiles)
        ys.append(y)

    return smiles_list, ys

def prepare_predictor(model):

    predictor = Predictor()

    predictor_input = os.path.join(model,
                                  'predictor_input.py')

    predictor.load_input(predictor_input)

    param_path = os.path.join(model,
                             'saved_model',
                             'full_train.h5')
    predictor.load_parameters(param_path)

    return predictor

def make_predictions(predictor, smiles_list):

    ys_cnn = []
    for smiles in tqdm(smiles_list):
        mol = Molecule().fromSMILES(smiles)
        y_cnn = predictor.predict(mol)
        ys_cnn.append(y_cnn)

    return ys_cnn

def evaluate(smiles_list, ys, ys_pred, prediction_task="Hf298(kcal/mol)"):

    result_df = pd.DataFrame(index=smiles_list)

    result_df[prediction_task+"_true"] = pd.Series(ys, index=result_df.index)
    result_df[prediction_task+"_pred"] = pd.Series(ys_pred, index=result_df.index)

    diff = abs(result_df[prediction_task+"_true"]-result_df[prediction_task+"_pred"])
    result_df[prediction_task+"_diff"] = pd.Series(diff, index=result_df.index)

    return result_df

def display_result(result_df, prediction_task="Hf298(kcal/mol)"):

    descr = result_df[prediction_task+"_diff"].describe()

    count = int(descr.loc['count'])
    mean = descr.loc['mean']
    std = descr.loc['std']

    display_str = 'prediction task: {0}, count: {1}, error mean: {2:.02f}, error std: {3:.02f}'.format(prediction_task, 
                                                                                                       count, mean, std) 
    print display_str

    return (count, mean, std) 

def validate(datasets_file, model):

    # load cnn predictor
    predictor = prepare_predictor(model)

    datasets = read_datasets_file(datasets_file)

    evaluation_results = {}
    for host, db_name, collection_name in datasets:
        
        print "\nhost: {0}, db: {1}, collection: {2}".format(host, db_name, collection_name)
        
        # prepare data for testing
        smiles_list, ys = prepare_data(host, db_name, collection_name,
                                       prediction_task=predictor.prediction_task)

        # Do the predictions
        ys_pred = make_predictions(predictor, smiles_list)

        # evaluate performance
        result_df = evaluate(smiles_list, ys, ys_pred, prediction_task=predictor.prediction_task)

        # display result
        count, mean, std = display_result(result_df, prediction_task=predictor.prediction_task)
        
        table = '.'.join([host, db_name, collection_name])
        evaluation_results[table] = {"count": count,
                                     "MAE": mean,
                                     "MAE std": std}

    return evaluation_results

def main():

    args = parseCommandLineArguments()

    datasets_file = args.datasets
    model = args.model
    evaluation_results = validate(datasets_file, model)

    with open('evaluation_results.json', 'w') as f_out:
        json.dump(evaluation_results, f_out, indent=4, sort_keys=True)

main()