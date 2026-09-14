#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains unit tests of the :mod:`arkane.pdep` module.
"""

import logging
import math
import os
import shutil

import pytest
import yaml

from rmgpy import settings
from rmgpy.chemkin import read_reactions_block
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.species import Species

from arkane.main import Arkane


@pytest.mark.functional
class ArkaneTest:
    """
    Contains unit tests for the sensitivity module in Arkane
    """

    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.directory = os.path.join(settings["test_data.directory"], "arkane", "tst1", "")
        cls.input_file = os.path.join(cls.directory, "pdep_sa.py")

        # clean working folder from all previous test output
        dirs = [d for d in os.listdir(cls.directory) if not os.path.isfile(os.path.join(cls.directory, d))]
        for d in dirs:
            shutil.rmtree(os.path.join(settings["test_data.directory"], "arkane", "tst1", d, ""))
        files = [f for f in os.listdir(cls.directory) if os.path.isfile(os.path.join(cls.directory, f))]
        for f in files:
            if "pdep_sa" not in f:
                os.remove(os.path.join(settings["test_data.directory"], "arkane", "tst1", f))

    def test_pdep_job(self):
        """
        A general test for a PDep job in Arkane
        """
        self.tst1 = Arkane()
        self.tst1.input_file = self.input_file
        self.tst1.output_directory = self.directory
        self.tst1.verbose = logging.WARN
        self.tst1.plot = False
        self.tst1.job_list = []
        self.tst1.job_list = self.tst1.load_input_file(self.tst1.input_file)
        self.tst1.execute()

        job = self.tst1.job_list[0]
        assert job.Tmin.value_si == 300.0
        assert job.minimum_grain_count == 100
        assert not job.rmgmode
        assert job.active_j_rotor
        assert job.network.path_reactions[0].label == "acetylperoxy <=> hydroperoxylvinoxy"
        assert round(abs(job.network.path_reactions[0].transition_state.tunneling.E0_TS.value_si - -24267.2), 7) == 0
        assert round(abs(job.network.path_reactions[0].transition_state.tunneling.frequency.value_si - -1679.04), 7) == 0
        assert len(job.network.net_reactions[0].reactants[0].conformer.modes) == 6

        # test that a network pdf was generated
        files = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        assert any(f == "network.pdf" for f in files)

        # Test the generated network reaction
        dictionary = {
            "hydroperoxylvinoxy": Species().from_smiles("[CH2]C(=O)OO"),
            "acetylperoxy": Species().from_smiles("CC(=O)O[O]"),
        }
        with open(os.path.join(self.directory, "chem.inp"), "r") as chem:
            reaction_list = read_reactions_block(chem, dictionary)
        rxn = reaction_list[0]
        assert isinstance(rxn.kinetics, Chebyshev)
        # Accept a delta of 0.2, which could result from numerical discrepancies
        # See RMG-Py #1682 on GitHub for discussion
        assert abs(rxn.kinetics.get_rate_coefficient(1000.0, 1.0) - 88.88253229631246) < 0.2

        files = [
            f for f in os.listdir(os.path.join(self.directory, "sensitivity", "")) if os.path.isfile(os.path.join(self.directory, "sensitivity", f))
        ]
        assert any("hydroperoxylvinoxy.pdf" in f for f in files)

        with open(os.path.join(self.directory, "sensitivity", "network1.txt"), "r") as f:
            lines = f.readlines()
            for line in lines:
                if "1000.0" in line:
                    break
        sa_coeff = line.split()[-2]
        assert abs(float(sa_coeff) - -7.02e-07) < 0.02e-6

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        cls.directory = os.path.join(settings["test_data.directory"], "arkane", "tst1", "")
        cls.input_file = os.path.join(cls.directory, "pdep_sa.py")

        # clean working folder from all previous test output
        dirs = [d for d in os.listdir(cls.directory) if not os.path.isfile(os.path.join(cls.directory, d))]
        for d in dirs:
            shutil.rmtree(os.path.join(settings["test_data.directory"], "arkane", "tst1", d, ""))
        files = [f for f in os.listdir(cls.directory) if os.path.isfile(os.path.join(cls.directory, f))]
        for f in files:
            if "pdep_sa" not in f:
                os.remove(os.path.join(settings["test_data.directory"], "arkane", "tst1", f))


@pytest.mark.functional
class ArkaneILTSensitivityTest:
    """
    Integration test for `arkane.sensitivity.PDepSensitivity` on a network that genuinely
    contains an ILT-based ("modeless", ``not rxn.can_tst()``) path reaction alongside RRKM-based
    ones, exercising the real perturb-both-E0-and-Ea code path (see arkane/sensitivity.py) rather
    than a synthetic/fabricated one.

    The fixture (test/rmgpy/test_data/arkane/tst_ilt_mixed/pdep_sa_ilt_mixed.py) reuses the real
    species/transitionState/reaction/network thermochemistry from the official
    examples/arkane/networks/acetyl+O2_cse example (its 'entrance1' transition state is modeless,
    with real Arrhenius kinetics on the owning reaction; its 'isom1'/'exit1'/'exit2'/'exit3'
    transition states are RRKM, with full statmech modes), with two deliberate deviations: the
    `pressureDependence(...)` block's computational settings (grid density/range, solver method)
    were narrowed to keep this test fast, and the 'entrance1' TS E0 was moved to 2.0 kcal/mol so
    that it does NOT satisfy E0(TS) = sum(E0(reactants)) + Ea. The latter is the regression shape
    that matters: real RMG-written network files generally break that relation, and an earlier
    revision of the sensitivity fix gated the Ea perturbation on it -- passing on this fixture
    (where every energy was 0.0 and the relation held trivially) while producing structurally
    zero TS coefficients on essentially all real networks (24 of 6,434 resolvable path reactions
    across 400 real RMG networks satisfied the relation, i.e. 0.4%).
    """

    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.directory = os.path.join(settings["test_data.directory"], "arkane", "tst_ilt_mixed", "")
        cls.input_file = os.path.join(cls.directory, "pdep_sa_ilt_mixed.py")

        # clear the working folder from any previous test output
        dirs = [d for d in os.listdir(cls.directory) if not os.path.isfile(os.path.join(cls.directory, d))]
        for d in dirs:
            shutil.rmtree(os.path.join(cls.directory, d, ""))
        files = [f for f in os.listdir(cls.directory) if os.path.isfile(os.path.join(cls.directory, f))]
        for f in files:
            if "pdep_sa_ilt_mixed" not in f:
                os.remove(os.path.join(cls.directory, f))

    def test_ilt_and_rrkm_path_reactions_produce_finite_sensitivity_coefficients(self):
        """
        Test that a real PDep sensitivity analysis job, run on a network with both an ILT-based
        and several RRKM-based path reactions, produces finite, non-all-zero TS sensitivity
        coefficients for the ILT-based transition state (whose Ea is perturbed alongside its E0),
        while leaving the RRKM-based transition states' coefficients well-defined too (they are
        unaffected by the Ea-perturbation logic, since it only applies when `not rxn.can_tst()`).
        """
        arkane = Arkane()
        arkane.input_file = self.input_file
        arkane.output_directory = self.directory
        arkane.verbose = logging.WARN
        arkane.plot = False
        arkane.job_list = []
        arkane.job_list = arkane.load_input_file(self.input_file)
        arkane.execute()

        job = arkane.job_list[0]
        path_reactions_by_label = {rxn.label: rxn for rxn in job.network.path_reactions}
        ilt_rxn = path_reactions_by_label["entrance1"]
        assert not ilt_rxn.can_tst()  # the genuinely ILT-based reaction
        for rrkm_label in ("isom1", "exit1", "exit2", "exit3"):
            assert path_reactions_by_label[rrkm_label].can_tst()

        # Guard the fixture's regression shape: the ILT TS's E0 must NOT satisfy
        # E0(TS) = sum(E0) + Ea on either side (real RMG-written networks generally break this
        # relation; an earlier revision gated the Ea perturbation on it and thus passed on a
        # relation-satisfying fixture while failing on real networks).
        ts_e0 = ilt_rxn.transition_state.conformer.E0.value_si
        ea = ilt_rxn.kinetics.Ea.value_si
        for side in (ilt_rxn.reactants, ilt_rxn.products):
            e0_sum = sum(spc.conformer.E0.value_si for spc in side)
            assert abs((ts_e0 - ea) - e0_sum) > 5000.0  # J/mol

        yaml_path = os.path.join(self.directory, "sensitivity", "sa_coefficients.yml")
        assert os.path.isfile(yaml_path)
        with open(yaml_path, "r") as f:
            sa_data = yaml.unsafe_load(f)

        ilt_ts_values = []
        rrkm_ts_values = []
        well_values = []
        for reaction_str, conditions in sa_data.items():
            if reaction_str in ("structures", "metadata"):
                continue
            for condition_data in conditions.values():
                for entry_label, coefficient in condition_data.items():
                    if entry_label == "(TS) entrance1":
                        ilt_ts_values.append(coefficient)
                    elif entry_label.startswith("(TS) "):
                        rrkm_ts_values.append(coefficient)
                    else:
                        well_values.append(coefficient)

        assert len(ilt_ts_values) > 0
        assert all(math.isfinite(v) for v in ilt_ts_values)
        assert any(v != 0.0 for v in ilt_ts_values)

        # The YAML carries a machine-discoverable marker of which TS rows are ILT-based (combined
        # E0+Ea perturbation), so a consumer need not re-derive can_tst() per reaction.
        assert "(TS) entrance1" in sa_data["metadata"]["ilt_transition_states"]

        # The ILT TS coefficients must be genuinely meaningful, not numerical noise: without the
        # Ea perturbation they sit many orders of magnitude below the well coefficients (the E0
        # perturbation alone is a structural no-op for an ILT-based reaction), so require the
        # largest ILT TS coefficient to be within two orders of magnitude of the largest well
        # coefficient.
        assert len(well_values) > 0
        max_well = max(abs(v) for v in well_values)
        max_ilt_ts = max(abs(v) for v in ilt_ts_values)
        assert max_well > 0.0
        assert max_ilt_ts >= 1e-2 * max_well

        assert len(rrkm_ts_values) > 0
        assert all(math.isfinite(v) for v in rrkm_ts_values)
        assert any(v != 0.0 for v in rrkm_ts_values)

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        cls.directory = os.path.join(settings["test_data.directory"], "arkane", "tst_ilt_mixed", "")

        # clean working folder from all previous test output
        dirs = [d for d in os.listdir(cls.directory) if not os.path.isfile(os.path.join(cls.directory, d))]
        for d in dirs:
            shutil.rmtree(os.path.join(cls.directory, d, ""))
        files = [f for f in os.listdir(cls.directory) if os.path.isfile(os.path.join(cls.directory, f))]
        for f in files:
            if "pdep_sa_ilt_mixed" not in f:
                os.remove(os.path.join(cls.directory, f))
