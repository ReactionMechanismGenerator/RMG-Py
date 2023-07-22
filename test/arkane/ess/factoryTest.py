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

"""
This module contains unit tests of the :mod:`arkane.ess.factory` module.
"""

import os


from rmgpy.exceptions import InputError

from arkane.ess import GaussianLog, MolproLog, OrcaLog, Psi4Log, QChemLog, TeraChemLog
from arkane.ess.factory import ess_factory
import pytest


class TestThermo:
    """
    Contains unit tests of the factory module.
    """

    @classmethod
    def setup_class(cls):
        """
        A method that is run before all unit tests in this class.
        """
        cls.data_path = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "..", "..", "arkane", "data"))

    def test_ess_factory(self):
        """Test identifying the electronic structure software from the log file"""
        gaussian_log_path1 = os.path.join(self.data_path, "gaussian", "ethylene_G3.log")
        gaussian_log_path2 = os.path.join(self.data_path, "gaussian", "oxygen.log")
        molpro_log_path1 = os.path.join(self.data_path, "molpro", "HOSI_ccsd_t1.out")
        molpro_log_path2 = os.path.join(self.data_path, "molpro", "molpro_mrci+q.out")
        orca_path_1 = os.path.join(self.data_path, "orca", "Orca_dlpno_test.log")
        orca_path_2 = os.path.join(self.data_path, "orca", "Orca_opt_freq_test.log")
        orca_path_3 = os.path.join(self.data_path, "orca", "Orca_TS_test.log")
        psi4_path_1 = os.path.join(self.data_path, "psi4", "opt_freq.out")
        psi4_path_2 = os.path.join(self.data_path, "psi4", "opt_freq_dft.out")
        psi4_path_3 = os.path.join(self.data_path, "psi4", "opt_freq_dft_ts.out")
        psi4_path_4 = os.path.join(self.data_path, "psi4", "opt_freq_ts.out")
        qchem_log_path1 = os.path.join(self.data_path, "qchem", "CH4_sp.out")
        qchem_log_path2 = os.path.join(self.data_path, "qchem", "co.out")
        terachem_log_path_1 = os.path.join(self.data_path, "terachem", "ethane_minimize_output.out")
        terachem_log_path_2 = os.path.join(self.data_path, "terachem", "formaldehyde_sp_terachem_output.out")
        terachem_log_path_3 = os.path.join(self.data_path, "terachem", "formaldehyde_sp_terachem_results.dat")
        terachem_log_path_4 = os.path.join(self.data_path, "terachem", "formaldehyde_coords.xyz")
        terachem_log_path_5 = os.path.join(self.data_path, "terachem", "formaldehyde_output.geometry")
        non_ess_log_path = os.path.abspath(os.path.join(self.data_path, "methoxy.py"))

        assert isinstance(ess_factory(gaussian_log_path1), GaussianLog)
        assert isinstance(ess_factory(gaussian_log_path2), GaussianLog)

        assert isinstance(ess_factory(molpro_log_path1), MolproLog)
        assert isinstance(ess_factory(molpro_log_path2), MolproLog)

        for orca_path in [orca_path_1, orca_path_2, orca_path_3]:
            assert isinstance(ess_factory(orca_path), OrcaLog)

        for psi4_path in [psi4_path_1, psi4_path_2, psi4_path_3, psi4_path_4]:
            assert isinstance(ess_factory(psi4_path), Psi4Log)

        assert isinstance(ess_factory(qchem_log_path1), QChemLog)
        assert isinstance(ess_factory(qchem_log_path2), QChemLog)

        for terachem_path in [
            terachem_log_path_1,
            terachem_log_path_2,
            terachem_log_path_3,
            terachem_log_path_4,
            terachem_log_path_5,
        ]:
            assert isinstance(ess_factory(terachem_path), TeraChemLog)

        with pytest.raises(InputError):
            ess_factory(non_ess_log_path)
