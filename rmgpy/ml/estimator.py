#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import os
import numpy as np

from dde.predictor import Predictor

from rmgpy.thermo import ThermoData


class MLEstimator():

    """
    A machine learning based estimator for thermochemistry prediction.

    The attributes are:

    ==================== ======================= =======================
    Attribute            Type                    Description
    ==================== ======================= =======================
    `hf298_estimator`    :class:`Predictor`      Hf298 estimator
    `hf298_uncertainty`  ``bool``                Hf298 uncertainty flag
    `s298_estimator`     :class:`Predictor`      S298 estimator
    `s298_uncertainty`   ``bool``                S298 uncertainty flag
    `cp_estimator`       :class:`Predictor`      Cp estimator
    `cp_uncertainty`     ``bool``                Cp uncertainty flag
    ==================== ======================= =======================

    """

    def __init__(self, hf298_path, s298_path, cp_path):
        self.hf298_estimator, self.hf298_uncertainty = load_estimator(hf298_path)
        self.s298_estimator, self.s298_uncertainty = load_estimator(s298_path)
        self.cp_estimator, self.cp_uncertainty = load_estimator(cp_path)

    def get_thermo_data(self, molecule):
        """
        Return thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule`. Also set the
        uncertainties estimated by the ML model if available.

        Returns: ThermoData
        """

        # These should correspond to the temperatures that the ML model was
        # trained on for Cp.
        temps = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0]

        hf298 = self.hf298_estimator.predict(molecule=molecule, sigma=self.hf298_uncertainty)
        s298 = self.s298_estimator.predict(molecule=molecule, sigma=self.s298_uncertainty)
        cp = self.cp_estimator.predict(molecule=molecule, sigma=self.cp_uncertainty)

        # If uncertainty is available for the ML model, a tuple of predicted
        # value and estimated uncertainty is returned. An uncertainty of None
        # gets set to a valua of 0 by :class:`Quantity`.
        hf298, hf298u = hf298 if isinstance(hf298, tuple) else (hf298, None)
        s298, s298u = s298 if isinstance(s298, tuple) else (s298, None)
        cp, cpu = cp if isinstance(cp, tuple) else (cp, None)

        cp = [np.float64(cp_i) for cp_i in cp]
        if cpu is not None:
            cpu = [np.float64(cpu_i) for cpu_i in cpu]

        cp0 = molecule.calculateCp0()
        cpinf = molecule.calculateCpInf()
        thermo = ThermoData(
            Tdata=(temps, 'K'),
            Cpdata=(cp, 'cal/(mol*K)', cpu),
            H298=(hf298, 'kcal/mol', hf298u),
            S298=(s298, 'cal/(mol*K)', s298u),
            Cp0=(cp0, 'J/(mol*K)'),
            CpInf=(cpinf, 'J/(mol*K)'),
            Tmin=(300.0, 'K'),
            Tmax=(2000.0, 'K'),
            comment='ML Estimation'
        )

        return thermo

    def get_thermo_data_for_species(self, species):
        """
        Return the set of thermodynamic parameters corresponding to a
        given :class:`Species` object `species`.

        The current ML estimator treats each resonance isomer
        identically, i.e., any of the resonance isomers can be chosen.

        Returns: ThermoData
        """
        return self.get_thermo_data(species.molecule[0])


def load_estimator(model_path):
    estimator = Predictor()

    input_file = os.path.join(model_path, 'predictor_input.py')
    weights_file = os.path.join(model_path, 'full_train.h5')
    model_file = os.path.join(model_path, 'full_train.json')
    mean_and_std_file = os.path.join(model_path, 'full_train_mean_std.npz')

    estimator.load_input(input_file)
    if os.path.exists(model_file):
        estimator.load_architecture(model_file)
        uncertainty = True
    else:
        uncertainty = False
    mean_and_std_file = mean_and_std_file if os.path.exists(mean_and_std_file) else None
    estimator.load_parameters(param_path=weights_file, mean_and_std_path=mean_and_std_file)

    return estimator, uncertainty
