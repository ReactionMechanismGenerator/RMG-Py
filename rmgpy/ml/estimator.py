#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

from typing import Union

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
    """

    # These should correspond to the temperatures that the ML model was
    # trained on for Cp.
    # we also train on cp0 and cpinf in 0,1 th index of prediction
    temps = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2400.0]

    def __init__(self, model_type: str = "attn_mpn", inference_type: str = "ensemble"):
        # once gnns thermo is updated in upstream this will change
        # lazy import here
        from gnns_thermo.inference import GNNCalculator, EnsembleGNNCalculator
        from gnns_thermo.testing import get_chkpt
        from gnns_thermo.config import ModelInferenceEnum

        self.model_type = model_type  # for logging
        self.inference_type = inference_type
        h298_chkpt, h298_config = get_chkpt(
            model_type, "h298", inference_type=inference_type
        )  # get a single checkpoint or folder full of checkpoints based on inference type
        s298_chkpt, s298_config = get_chkpt(
            model_type, "s298", inference_type=inference_type
        )
        cp_chkpt, cp_config = get_chkpt(model_type, "cp", inference_type=inference_type)
        if inference_type == ModelInferenceEnum.ensemble:
            calculator_type = EnsembleGNNCalculator
        elif inference_type == ModelInferenceEnum.single_model:
            calculator_type = GNNCalculator

        else:
            raise NotImplementedError(f"Given inference type is not supported")
        # make calculator objects here
        self.hf298_estimator = calculator_type(
            model_type, "cpu", h298_chkpt, h298_config
        )
        self.s298_estimator = calculator_type(
            model_type, "cpu", s298_chkpt, s298_config
        )
        self.cp_estimator = calculator_type(model_type, "cpu", cp_chkpt, cp_config)

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
            input = [molecule.to_rdkit_mol()]
        else:
            raise NotImplementedError(
                f"Give mode {mode} is not implemented for ML stimation"
            )
        hf298_pred = self.hf298_estimator.calculate(input)
        s298_pred = self.s298_estimator.calculate(input)
        cp_pred = self.cp_estimator.calculate(input)
        if self.inference_type == "ensemble":
            hf298, hf298_std = hf298_pred.mean(0), hf298_pred.std(0)
            s298, s298_std = s298_pred.mean(0), s298_pred.std(0)
            cp, cp_std = cp_pred.mean(0), cp_pred.std(0)
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

        cp0, cpinf = cp[0].astype(np.float64), cp[1].astype(np.float64)
        cp = cp[2:].astype(np.float64)

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
        else: # This species has resonance structures
            thermo = []
            for molecule in species.molecule:
                tdata = self.get_thermo_data(molecule)
                thermo.append((tdata.get_enthalpy(298.0), molecule, tdata))
            thermo = sorted(thermo, key=lambda x: x[0]) # sort by ascending H298
            species.molecule = [item[1] for item in thermo] # reorder molecules by ascending H298
            return thermo[0][2] # return the thermodata of lowest energy resonance structure

