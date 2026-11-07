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

from unittest.mock import patch

"""
Unit tests for how the RMS reactors in rmgpy.rmg.reactionmechanismsimulator_reactors call
ReactionMechanismSimulator's edge analysis. Julia is replaced by a mock, so these tests run
without RMS installed.
"""

import logging
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

import rmgpy.rmg.reactionmechanismsimulator_reactors as rms_reactors
from rmgpy.rmg.settings import SimulatorSettings


def _make_reactor():
    """A Reactor with generate_reactor stubbed out, skipping __init__ (which needs RMS objects)."""
    reactor = rms_reactors.Reactor.__new__(rms_reactors.Reactor)
    reactor.core_phase_system = "core"
    reactor.edge_phase_system = "edge"
    reactor.terminations = []
    reactor.generate_reactor = lambda phase_system: (phase_system + "_react", [], [], None)
    return reactor


def _mock_main(resurrected=False, invalid_objects=()):
    """A mock of juliacall's Main whose selectobjects returns the given resurrection result."""
    main = MagicMock()
    main.Symbol.side_effect = lambda name: ("Symbol", name)
    main.Sundials.CVODE_BDF.side_effect = lambda **kwargs: ("CVODE_BDF", kwargs)
    main.selectobjects.return_value = (
        False,  # terminated
        resurrected,
        list(invalid_objects),
        "uni",
        "bi",
        "tri",
        "rate_ratios",
        1.0e-11,  # t
        0.0,  # conversion
    )
    return main


class TestRMSReactorSimulate:
    @pytest.mark.parametrize(
        "linear_solver,expected",
        [("dense", "Dense"), ("GMRES", "GMRES")],
    )
    def test_linear_solver_passed_to_cvode(self, linear_solver, expected):
        """Reactor.simulate builds CVODE_BDF with the linear solver chosen in the simulator settings"""
        main = _mock_main()
        with patch.object(rms_reactors, "Main", main):
            _make_reactor().simulate(MagicMock(), SimulatorSettings(linear_solver=linear_solver), None)
        main.Sundials.CVODE_BDF.assert_called_once_with(linear_solver=("Symbol", expected))
        kwargs = main.selectobjects.call_args.kwargs
        assert kwargs["solver"] == ("CVODE_BDF", {"linear_solver": ("Symbol", expected)})
        assert kwargs["atol"] == 1e-16
        assert kwargs["rtol"] == 1e-8

    def test_default_is_dense(self):
        """With default simulator settings, CVODE uses a dense direct linear solver, not GMRES"""
        main = _mock_main()
        with patch.object(rms_reactors, "Main", main):
            _make_reactor().simulate(MagicMock(), SimulatorSettings(), None)
        main.Sundials.CVODE_BDF.assert_called_once_with(linear_solver=("Symbol", "Dense"))

    def test_resurrection_is_logged(self, caplog):
        """When RMS resurrects a failed simulation, RMG logs a warning naming the species added"""
        main = _mock_main(resurrected=True, invalid_objects=[SimpleNamespace(name="CO2X(12)")])
        with patch.object(rms_reactors, "Main", main), caplog.at_level(logging.WARNING):
            result = _make_reactor().simulate(MagicMock(), SimulatorSettings(), None)
        assert result[1] is True
        assert "Resurrecting the model" in caplog.text
        assert "CO2X(12)" in caplog.text

    def test_resurrection_with_empty_edge_is_logged(self, caplog):
        """A resurrection with no edge species to add still logs a warning, without naming a species"""
        main = _mock_main(resurrected=True, invalid_objects=[])
        with patch.object(rms_reactors, "Main", main), caplog.at_level(logging.WARNING):
            _make_reactor().simulate(MagicMock(), SimulatorSettings(), None)
        assert "Resurrecting the model" in caplog.text
        assert "no edge species to add" in caplog.text

    def test_no_warning_without_resurrection(self, caplog):
        """A successful simulation does not log the resurrection warning"""
        main = _mock_main(resurrected=False)
        with patch.object(rms_reactors, "Main", main), caplog.at_level(logging.WARNING):
            _make_reactor().simulate(MagicMock(), SimulatorSettings(), None)
        assert "Resurrecting" not in caplog.text
