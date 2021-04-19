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
    # we also train on cp0 and cpinf in 0,1 th index of prediction
    temps = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2400.0]

    def __init__(self, model_type: str = "attn_mpn"):
        # once gnns thermo is updated in upstream this will change
        # lazy import here cause segfaults with torch
        from gnns_thermo.inference import GNNCalculator
        from gnns_thermo.testing import get_chkpt

        self.model_type = model_type
        h298_chkpt, h298_config = get_chkpt(model_type, "h298")
        s298_chkpt, s298_config = get_chkpt(model_type, "s298")
        cp_chkpt, cp_config = get_chkpt(model_type, "cp")
        self.hf298_estimator = GNNCalculator(model_type, "cpu", h298_chkpt, h298_config)
        self.s298_estimator = GNNCalculator(model_type, "cpu", s298_chkpt, s298_config)
        self.cp_estimator = GNNCalculator(model_type, "cpu", cp_chkpt, cp_config)

    def get_thermo_data(self, molecule: Union[Molecule, str]) -> ThermoData:
        """
        Return thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` or a SMILES string.

        Returns: ThermoData
        """
        molecule = Molecule(smiles=molecule) if isinstance(molecule, str) else molecule
        # convert to float from np.float64
        hf298 = self.hf298_estimator.calculate([molecule.smiles])
        s298 = self.s298_estimator.calculate([molecule.smiles])
        cp = self.cp_estimator.calculate([molecule.smiles])
        cp0, cpinf = cp[0], cp[1]
        cp = cp[2:]

        # cp0 = molecule.calculate_cp0()
        # cpinf = molecule.calculate_cpinf()

        # Set uncertainties to 0 because the current model cannot estimate them
        thermo = ThermoData(
            Tdata=(self.temps, "K"),
            Cpdata=(cp, "J/(mol*K)", np.zeros(len(self.temps))),
            H298=(hf298, "kcal/mol", 0),
            S298=(s298, "J/(mol*K)", 0),
            Cp0=(cp0, "J/(mol*K)"),
            CpInf=(cpinf, "J/(mol*K)"),
            Tmin=(300.0, "K"),
            Tmax=(2400.0, "K"),
            comment=f"ML Estimation with model type {self.model_type}",
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
