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

import unittest

from rmgpy.ml.estimator import MLEstimator


class TestMLEstimator(unittest.TestCase):
    """
    Contains unit tests for rmgpy.ml.estimator
    """

    def setUp(self):
        """
        Set up the MLEstimator class. This method is run once before all
        other unit tests.
        """
        # ensemble attn_mpn (best for 2D SMILES -> prop)
        self.ml_estimator = MLEstimator("attn_mpn")
        # single attn_mpn
        self.ml_estimator_single_model = MLEstimator(
            "attn_mpn", inference_type="single_model"
        )
        # ensemble mpnn with sum pooling
        self.ml_estimator_mpnn = MLEstimator("mpnn")
        # dimenetpp with xTB1 geometry optimization
        self.ml_estimator_dimenetpp = MLEstimator(
            "dimenetpp", inference_type="single_model"
        )

    def test_get_thermo_data_from_smiles_ensemble(self):
        """
        Test that we can make a prediction using MLEstimator using gnns_thermo ensemble models and smiles as input.
        """
        smi = "C1C2C1C2"
        thermo = self.ml_estimator.get_thermo_data(smi, mode="from_smiles")
        thermo_mpnn = self.ml_estimator_mpnn.get_thermo_data(smi, mode="from_smiles")
        self.assertTrue(
            thermo.comment.startswith("ML Estimation using featurizer from_smiles")
        )
        # regression tests with qm9 all systematic datasets (april 20)
        self.assertAlmostEqual(
            thermo.Cp0.value_si, 33.15302276611328, 1, msg="Cp0 regression error"
        )
        self.assertAlmostEqual(
            thermo.CpInf.value_si, 232.6637420654297, 1, msg="CpInf regression error"
        )
        self.assertEqual(len(thermo.Cpdata.value_si), 9)
        # check we get some uncertainty values
        self.assertGreater(
            abs(thermo.H298.uncertainty), 0.0, msg="No uncertainty values"
        )
        # only check uncertainty for other models
        self.assertGreater(
            abs(thermo_mpnn.H298.uncertainty), 0.0, msg="No uncertainty values"
        )

    def test_get_thermo_data_from_rdkit_mol_ensemble(self):
        """
        Test that we can make a prediction using MLEstimator using gnns_thermo ensemble models and rdkit mol as input.
        """
        smi = "C1C2C1C2"
        thermo = self.ml_estimator.get_thermo_data(smi, mode="from_rdkit_mol")
        thermo_mpnn = self.ml_estimator_mpnn.get_thermo_data(smi, mode="from_rdkit_mol")
        self.assertTrue(
            thermo.comment.startswith("ML Estimation using featurizer from_rdkit")
        )

        # regression tests with qm9 all systematic datasets (april 20)
        self.assertAlmostEqual(
            thermo.Cp0.value_si, 33.15302276611328, 1, msg="Cp0 regression error"
        )
        self.assertAlmostEqual(
            thermo.CpInf.value_si, 232.6637420654297, 1, msg="CpInf regression error"
        )
        self.assertEqual(len(thermo.Cpdata.value_si), 9)
        # check we get some uncertainty values
        self.assertGreater(
            abs(thermo.H298.uncertainty), 0.0, msg="No uncertainty values"
        )
        # only check uncertainty for other models
        self.assertGreater(
            abs(thermo_mpnn.H298.uncertainty), 0.0, msg="No uncertainty values"
        )

    def test_get_thermo_data_from_rdkit_mol_single_model(self):
        """
        Test that we can make a prediction using MLEstimator using gnns_thermo single model.
        """
        smi = "C1C2C1C2"
        thermo = self.ml_estimator_single_model.get_thermo_data(
            smi, mode="from_rdkit_mol"
        )
        self.assertTrue(
            thermo.comment.startswith("ML Estimation using featurizer from_rdkit")
        )
        self.assertEqual(len(thermo.Cpdata.value_si), 9)

    def test_get_thermo_data_with_dimenetpp_from_smiles_single(self):
        """
        Test that we can make a prediction using MLEstimator using gnns_thermo dimenetpp single model and smiles as input.
        """
        smi = "C1C2C1C2"
        thermo = self.ml_estimator_dimenetpp.get_thermo_data(smi, mode="from_smiles")
        self.assertTrue(
            thermo.comment.startswith("ML Estimation using featurizer from_smiles")
        )
        # regression tests with qm9 all systematic datasets (april 20)
        self.assertAlmostEqual(
            thermo.Cp0.value_si, 33.74019241333008, 1, msg="Cp0 regression error"
        )
        # self.assertAlmostEqual(
        # thermo.CpInf.value_si, 232.6637420654297, 1, msg="CpInf regression error"
        # )
        self.assertEqual(len(thermo.Cpdata.value_si), 9)
        # check we get some uncertainty values
        # self.assertGreater(
        # abs(thermo.H298.uncertainty), 0.0, msg="No uncertainty values"
        # )

    def test_convert_thermo_data(self):
        """
        Test that we can make a prediction using gnns_thermo and convert to wilholt.
        """
        smi = "C1C2C1C2"
        thermo = self.ml_estimator.get_thermo_data(smi)
        wilhoit = thermo.to_wilhoit(B=1000.0)
        s_nasa = thermo.get_entropy(1000.0)

        s_nasa_to_wh = wilhoit.get_entropy(1000.0)

        self.assertAlmostEqual(
            s_nasa,
            s_nasa_to_wh,
            -1,
            msg="this will break if there are any type issues with quantities",
        )
