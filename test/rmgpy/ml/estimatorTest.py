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

import os


from rmgpy import settings
from rmgpy.ml.estimator import MLEstimator


class TestMLEstimator:
    """
    Contains unit tests for rmgpy.ml.estimator
    """

    def setup_class(self):
        """
        Set up the MLEstimator class. This method is run once before all
        other unit tests.
        """
        models_path = os.path.join(settings["database.directory"], "thermo", "ml", "main")
        hf298_path = os.path.join(models_path, "hf298")
        s298_cp_path = os.path.join(models_path, "s298_cp")
        self.ml_estimator = MLEstimator(hf298_path, s298_cp_path)

    def test_get_thermo_data(self):
        """
        Test that we can make a prediction using MLEstimator.
        """
        smi = "C1C2C1C2"
        thermo = self.ml_estimator.get_thermo_data(smi)

        assert thermo.comment.startswith("ML Estimation")
        assert round(abs(thermo.Cp0.value_si - 33.3), 1) == 0
        assert round(abs(thermo.CpInf.value_si - 232.8), 1) == 0
        assert len(thermo.Cpdata.value_si) == 7
