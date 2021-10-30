#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
from enum import Enum

class FeatEnum(str, Enum):
    from_smiles = "from_smiles"
    from_rdkit_mol = "from_rdkit_mol"


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

class GNNEstimator:
    """
    A machine learning based estimator for thermochemistry prediction.  
    """

    # These should correspond to the temperatures that the ML model was
    # trained on for Cp.
    # we also train on cp0 and cpinf in 0,1 th index of prediction
    temps = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2400.0]

    def __init__(
        self,
        model_type: str = "attn_mpn",
        inference_type: str = "ensemble",
        dimenetpp_featurizer: str = "gfn1",
    ):
        # once gnns thermo is updated in upstream this will change
        # lazy import here
        from gnns_thermo.inference import (
            GNNCalculator,
            EnsembleGNNCalculator,
            EnsembleDimeNetPPCalculator,
        )
        from gnns_thermo.testing import get_chkpt
        from gnns_thermo.config.enums import ModelInferenceEnum

        self.model_type = model_type  # for logging
        self.dimenetpp_featurizer = dimenetpp_featurizer
        self.inference_type = inference_type
        h298_chkpt, h298_config = get_chkpt(
            model_type, "h298", inference_type=inference_type
        )  # get a single checkpoint or folder full of checkpoints based on inference type
        s298_chkpt, s298_config = get_chkpt(
            model_type, "s298", inference_type=inference_type
        )
        cp_chkpt, cp_config = get_chkpt(model_type, "cp", inference_type=inference_type)
        # can be offloaded to a method when everything is in upstream
        if inference_type == ModelInferenceEnum.ensemble:
            if model_type in ["attn_mpn", "mpnn", "dmpnn"]:
                calculator_type = EnsembleGNNCalculator
            elif model_type == "dimenetpp":
                calculator_type = EnsembleDimeNetPPCalculator
            else:
                raise NotImplementedError(
                    f"Given model type {model_type} is not implemented"
                )
        else:
            raise NotImplementedError(
                f"Given inference type {inference_type} is not supported"
            )
        # make calculator objects here
        # ensemble calculators will raise an error here if a sinlge checkpoint is given
        self.hf298_estimator = calculator_type(
            model_type=model_type,
            device="cpu",
            models_path=h298_chkpt,
            config_file=h298_config,
        )
        self.s298_estimator = calculator_type(
            model_type=model_type,
            device="cpu",
            models_path=s298_chkpt,
            config_file=s298_config,
        )
        self.cp_estimator = calculator_type(
            model_type=model_type,
            device="cpu",
            models_path=cp_chkpt,
            config_file=cp_config,
        )

    def get_thermo_data(
        self, molecule: Union[Molecule, str], mode: FeatEnum = FeatEnum.from_rdkit_mol
    ) -> ThermoData:
        """
        Return thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` or a SMILES string.

        Returns: ThermoData
        """
        molecule = Molecule(smiles=molecule) if isinstance(molecule, str) else molecule
        if mode == FeatEnum.from_smiles:
            # calculator takes a list so we can batch if needed
            input = [molecule.smiles]
        elif mode == FeatEnum.from_rdkit_mol:
            input = [molecule.to_rdkit_mol(remove_h=False)]
        else:
            raise NotImplementedError(
                f"Give mode {mode} is not implemented for ML stimation"
            )
        if self.model_type == "dimenetpp":
            uhf = molecule.multiplicity - 1
        else:
            uhf = None
        # need to find
        hf298_pred = self.hf298_estimator.calculate(
            input, uhf=uhf, feat=self.dimenetpp_featurizer
        )
        s298_pred = self.s298_estimator.calculate(
            input, uhf=uhf, feat=self.dimenetpp_featurizer
        )
        cp_pred = self.cp_estimator.calculate(
            input, uhf=uhf, feat=self.dimenetpp_featurizer
        )
        if self.inference_type == "ensemble":
            hf298, hf298_std = hf298_pred.mean(0), hf298_pred.std(0)
            s298, s298_std = s298_pred.mean(0), s298_pred.std(0)
            cp, cp_std = cp_pred.mean(1), cp_pred.std(1)
            # explcit cast for cp
        elif self.inference_type == "single_model":
            # some weird quantity errors and this is a hack
            hf298, hf298_std = np.array([hf298_pred]).mean(), np.array([0.0]).mean()
            s298, s298_std = np.array([s298_pred]).mean(), np.array([0.0]).mean()
            cp, cp_std = cp_pred, 0.0
        else:
            raise NotImplementedError(f"Given {self.inference_type} is not supported")
        # type conversion for quantity handling
        hf298, hf298_std = hf298.astype(np.float64), hf298_std.astype(np.float64)
        s298, s298_std = s298.astype(np.float64), s298_std.astype(np.float64)

        cp0 = molecule.calculate_cp0()
        cpinf = molecule.calculate_cpinf()

        cp = cp.astype(np.float64)

        # Set uncertainties to std
        thermo = ThermoData(
            Tdata=(self.temps, "K"),
            Cpdata=(cp, "J/(mol*K)", np.zeros(len(self.temps))),
            H298=(hf298, "kcal/mol", hf298_std),
            S298=(s298, "J/(mol*K)", s298_std),
            Cp0=(cp0, "J/(mol*K)"),
            CpInf=(cpinf, "J/(mol*K)"),
            Tmin=(300.0, "K"),
            Tmax=(2400.0, "K"),
            comment=f"ML Estimation using featurizer {mode} with model type {self.model_type} in {self.inference_type} mode",
        )

        return thermo

    def get_thermo_data_for_species(self, species: Species) -> ThermoData:
        """
        Return the set of thermodynamic parameters corresponding to a
        given :class:`Species` object `species`.

        If the `species` contains resonance structures, thermo is estimated for
        each `molecule` and the most stable thermodata (with the lowest H298)
        is returned. The `species` molecule list is reorder by ascending H298.

        Returns: ThermoData
        """

        if len(species.molecule) == 1:
            return self.get_thermo_data(species.molecule[0])
        else:  # This species has resonance structures
            thermo = []
            for molecule in species.molecule:
                tdata = self.get_thermo_data(molecule)
                thermo.append((tdata.get_enthalpy(298.0), molecule, tdata))
            thermo = sorted(thermo, key=lambda x: x[0])  # sort by ascending H298
            species.molecule = [
                item[1] for item in thermo
            ]  # reorder molecules by ascending H298
            return thermo[0][
                2
            ]  # return the thermodata of lowest energy resonance structure


