###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import contextlib
import os
from argparse import Namespace
from typing import Callable, Union

try:
    import chemprop
except ImportError as e:
    chemprop = None
    chemprop_exception = e
import numpy as np

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


class MLEstimator:
    """
    A machine learning based estimator for thermochemistry prediction.

    The attributes are:

    ==================== ======================= =======================
    Attribute            Type                    Description
    ==================== ======================= =======================
    `hf298_estimator`    :class:`Predictor`      Hf298 estimator
    `s298_cp_estimator`  :class:`Predictor`      S298 and Cp estimator
    `temps`              ``list``                Cp temperatures
    ==================== ======================= =======================
    """

    # These should correspond to the temperatures that the ML model was
    # trained on for Cp.
    temps = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0]

    def __init__(self, hf298_path: str, s298_cp_path: str):
        self.hf298_estimator = load_estimator(hf298_path)
        self.s298_cp_estimator = load_estimator(s298_cp_path)

    def get_thermo_data(self, molecule: Union[Molecule, str]) -> ThermoData:
        """
        Return thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` or a SMILES string.

        Returns: ThermoData
        """
        molecule = Molecule(smiles=molecule) if isinstance(molecule, str) else molecule

        hf298 = self.hf298_estimator(molecule.smiles)[0][0]
        s298_cp = self.s298_cp_estimator(molecule.smiles)[0]
        s298, cp = s298_cp[0], s298_cp[1:]

        cp0 = molecule.calculate_cp0()
        cpinf = molecule.calculate_cpinf()

        # Set uncertainties to 0 because the current model cannot estimate them
        thermo = ThermoData(
            Tdata=(self.temps, 'K'),
            Cpdata=(cp, 'cal/(mol*K)', np.zeros(len(self.temps))),
            H298=(hf298, 'kcal/mol', 0),
            S298=(s298, 'cal/(mol*K)', 0),
            Cp0=(cp0, 'J/(mol*K)'),
            CpInf=(cpinf, 'J/(mol*K)'),
            Tmin=(300.0, 'K'),
            Tmax=(2000.0, 'K'),
            comment='ML Estimation'
        )

        return thermo

    def get_thermo_data_for_species(self, species: Species) -> ThermoData:
        """
        Return the set of thermodynamic parameters corresponding to a
        given :class:`Species` object `species`.

        The current ML estimator treats each resonance isomer
        identically, i.e., any of the resonance isomers can be chosen.

        Returns: ThermoData
        """
        return self.get_thermo_data(species.molecule[0])


def load_estimator(model_dir: str) -> Callable[[str], np.ndarray]:
    """
    Load chemprop model and return function for evaluating it.
    """
    if chemprop is None:
        # Delay chemprop ImportError until we actually try to use it
        # so that RMG can load successfully without chemprop.
        raise chemprop_exception

    args = Namespace()  # Simple class to hold attributes

    # Set up chemprop predict arguments
    args.checkpoint_dir = model_dir
    args.checkpoint_path = None
    chemprop.parsing.update_checkpoint_args(args)
    args.cuda = False

    scaler, features_scaler = chemprop.utils.load_scalers(args.checkpoint_paths[0])
    train_args = chemprop.utils.load_args(args.checkpoint_paths[0])

    # Update args with training arguments
    for key, value in vars(train_args).items():
        if not hasattr(args, key):
            setattr(args, key, value)

    # Load models in ensemble
    models = []
    for checkpoint_path in args.checkpoint_paths:
        models.append(chemprop.utils.load_checkpoint(checkpoint_path, cuda=args.cuda))

    # Set up estimator
    def estimator(smi: str):
        # Make dataset
        data = chemprop.data.MoleculeDataset(
            [chemprop.data.MoleculeDatapoint(line=[smi], args=args)]
        )

        # Normalize features
        if train_args.features_scaling:
            data.normalize_features(features_scaler)

        # Redirect chemprop stderr to null device so that it doesn't
        # print progress bars every time a prediction is made
        with open(os.devnull, 'w') as f, contextlib.redirect_stderr(f):
            # Predict with each model individually and sum predictions
            sum_preds = np.zeros((len(data), args.num_tasks))
            for model in models:
                model_preds = chemprop.train.predict(
                    model=model,
                    data=data,
                    batch_size=1,  # We'll only predict one molecule at a time
                    scaler=scaler
                )
                sum_preds += np.array(model_preds)

        avg_preds = sum_preds / len(models)
        return avg_preds

    return estimator
