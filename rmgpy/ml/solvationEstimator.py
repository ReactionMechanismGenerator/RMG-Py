from chemprop_solvation.utils import load_args, load_checkpoint, load_scalers
from chemprop_solvation.data import MoleculeDatapoint, MoleculeDataset
from tqdm import tqdm
from chemprop_solvation.train.predict import predict
import numpy as np

def predict_solvation(smiles: list(list()) = None, prop: str = 'Gsolv'):
    """Computes the desired solvation property using chemprop_solvation.
    Supported properties are 'Gsolv'"""
    
    # define the path to the models used to make predictions
    model_path = []
    if prop == 'Gsolv':
        model_path = ['./../Yunsie/test/fold_0/model_0/model.pt', './../Yunsie/test/fold_0/model_1/model.pt']
    else:
        raise ValueError(f'Property "{prop}" not supported.')
        
    # load the arguments and scalers to normalize predictions and features
    train_args = load_args(model_path[0])
    scaler, features_scaler = load_scalers(model_path[0])
    
    # convert the smiles, remove invalid smiles from the list
    data = MoleculeDataset([MoleculeDatapoint(sm, train_args) for sm in smiles])
    valid_indices = [i for i in range(len(data)) if data[i].solvent_mol is not None and data[i].solute_mol is not None]
    data = MoleculeDataset([data[i] for i in valid_indices])
    if len(data) == 0:
        raise ValueError(f'No valid smiles are given.')
    if train_args.features_scaling:
        data.normalize_features(features_scaler)
        
    # make predictions for all models
    all_preds = []
    for checkpoint_path in tqdm(model_path, total=len(model_path)):
        # load model
        model = load_checkpoint(checkpoint_path, cuda=False)
        model_preds = predict(
            model=model,
            data=data,
            batch_size=1,
            scaler=scaler
        )
        all_preds.append(np.array(model_preds))
    
    # calculate average prediction and variance on prediction (=epistemic or model uncertainty)
    epi_unc = []
    avg_pre = []
    all_preds = np.array(all_preds)
    for i in range(len(model_path)):
        epi_unc.append(np.var(all_preds[:,i]))
        avg_pre.append(np.mean(all_preds[:,i]))
    
    return avg_pre, epi_unc
