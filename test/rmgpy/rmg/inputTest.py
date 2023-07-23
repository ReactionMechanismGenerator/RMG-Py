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

from unittest.mock import patch

import rmgpy.rmg.input as inp
from rmgpy.rmg.main import RMG

import pytest


def setup_module(module):
    """
    A method that is run before the class.
    """
    # set-up RMG object and get global rmg object in input.py file
    # so methods can be tested
    global rmg
    rmg = RMG()
    inp.set_global_rmg(rmg)


def teardown_module(module):
    # remove RMG object
    global rmg
    rmg = None


class TestInputDatabase:
    """
    Contains unit tests rmgpy.rmg.input.database
    """

    def teardown_class(self):
        # remove the reactionLibraries value
        global rmg
        rmg.reaction_libraries = None

    def test_importing_database_reaction_libraries_from_string(self):
        """
        Test that we can import Reaction Libraries using the non-tuple form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=["test"])
        assert isinstance(rmg.reaction_libraries[0], tuple)
        assert not rmg.reaction_libraries[0][1]

    def test_importing_database_reaction_libraries_from_false_tuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple False form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[("test", False)])
        assert isinstance(rmg.reaction_libraries[0], tuple)
        assert not rmg.reaction_libraries[0][1]

    def test_importing_database_reaction_libraries_from_true_tuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple True form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[("test", True)])
        assert isinstance(rmg.reaction_libraries[0], tuple)
        assert rmg.reaction_libraries[0][1]


class TestInputMLEstimator:
    """
    Contains unit tests rmgpy.rmg.input.mlEstimator
    """

    def teardown_class(self):
        # remove the reactionLibraries value
        global rmg
        rmg.ml_estimator = None

    def test_ml_estimator(self):
        """
        Test that we can input.
        """
        from rmgpy.ml.estimator import MLEstimator

        global rmg
        # add database properties to RMG
        inp.ml_estimator(thermo=True)
        assert isinstance(rmg.ml_estimator, MLEstimator)
        assert isinstance(rmg.ml_settings, dict)


class TestInputThemoCentralDatabase:
    """
    Contains unit tests rmgpy.rmg.input.thermo_central_database
    """

    def teardown_class(self):
        # remove the reactionLibraries value
        global rmg
        rmg.thermo_central_database = None

    def test_themo_central_database(self):
        """
        Test that we can input.
        """
        global rmg
        # add database properties to RMG
        inp.thermo_central_database(
            host="some_host",
            port=0,
            username="some_usr",
            password="some_pw",
            application="some_app",
        )
        assert rmg.thermo_central_database.host == "some_host"
        assert rmg.thermo_central_database.port == 0
        assert rmg.thermo_central_database.username == "some_usr"
        assert rmg.thermo_central_database.password == "some_pw"
        assert rmg.thermo_central_database.application == "some_app"
        assert rmg.thermo_central_database.client == None


class TestInputReactors:
    """
    Contains unit tests for reactor input classes
    """

    @pytest.fixture(autouse=True)
    def setup_rmg(self):
        """This method is run before every test in this class"""
        # Create a mock species dictionary
        # In reality, the values would be Species objects, but it doesn't matter for testing
        species_dict = {
            "A": "A",
            "B": "B",
            "C": "C",
            "X": "X",
        }

        # Assign to global variable in the input module
        inp.species_dict = species_dict

        # Initialize the rmg.reaction_systems attribute
        global rmg
        rmg.reaction_systems = []
        yield

    def teardown_class(self):
        """This method is run after every test in this class"""
        # Reset the global species_dict variable in the input module
        inp.species_dict = {}

        # Reset the rmg.reaction_systems attribute
        global rmg
        rmg.reaction_systems = []

    def test_simple_reactor_mole_fractions(self):
        """Test that SimpleReactor mole fractions are set properly"""
        inp.simple_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 0.5,
                "B": 0.3,
                "C": 0.2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

    @patch("rmgpy.rmg.input.logging")
    def test_simple_reactor_mole_fractions_normalize_1(self, mock_logging):
        """Test that SimpleReactor mole fractions are normalized properly"""
        inp.simple_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial mole fractions do not sum to one; normalizing.")

    @patch("rmgpy.rmg.input.logging")
    def test_simple_reactor_mole_fractions_normalize_2(self, mock_logging):
        """Test that SimpleReactor mole fractions are normalized properly"""
        inp.simple_reactor(
            temperature=[(1000, "K"), (2000, "K")],
            pressure=[(1, "atm"), (10, "atm")],
            initialMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial mole fractions do not sum to one; normalizing.")

    def test_simple_reactor_mole_fractions_ranged(self):
        """Test that SimpleReactor ranged mole fractions are not normalized"""
        inp.simple_reactor(
            temperature=[(1000, "K"), (2000, "K")],
            pressure=[(1, "atm"), (10, "atm")],
            initialMoleFractions={
                "A": [5, 8],
                "B": 3,
                "C": 2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == [5, 8]
        assert reactor.initial_mole_fractions["B"] == 3
        assert reactor.initial_mole_fractions["C"] == 2

    def test_liquid_reactor_concentrations(self):
        """Test that LiquidReactor concentrations are set properly"""
        inp.liquid_reactor(
            temperature=(1000, "K"),
            initialConcentrations={
                "A": (0.3, "mol/L"),
                "B": (0.2, "mol/L"),
                "C": (0.1, "mol/L"),
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]

        # Values get converted to default SI units, mol/m^3
        assert reactor.initial_concentrations["A"] == 300
        assert reactor.initial_concentrations["B"] == 200
        assert reactor.initial_concentrations["C"] == 100

    def test_surface_reactor_mole_fractions(self):
        """Test that SurfaceReactor mole fractions are set properly"""
        inp.surface_reactor(
            temperature=(1000, "K"),
            initialPressure=(1, "atm"),
            initialGasMoleFractions={
                "A": 0.5,
                "B": 0.3,
                "C": 0.2,
            },
            initialSurfaceCoverages={"X": 1.0},
            surfaceVolumeRatio=(1e1, "m^-1"),
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_gas_mole_fractions["A"] == 0.5
        assert reactor.initial_gas_mole_fractions["B"] == 0.3
        assert reactor.initial_gas_mole_fractions["C"] == 0.2

    @patch("rmgpy.rmg.input.logging")
    def test_surface_reactor_mole_fractions_normalize_1(self, mock_logging):
        """Test that SurfaceReactor mole fractions are normalized properly"""
        inp.surface_reactor(
            temperature=(1000, "K"),
            initialPressure=(1, "atm"),
            initialGasMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            initialSurfaceCoverages={"X": 1.0},
            surfaceVolumeRatio=(1e1, "m^-1"),
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_gas_mole_fractions["A"] == 0.5
        assert reactor.initial_gas_mole_fractions["B"] == 0.3
        assert reactor.initial_gas_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial gas mole fractions do not sum to one; renormalizing.")

    def test_mb_sampled_reactor_mole_fractions(self):
        """Test that MBSampledReactor mole fractions are set properly"""
        inp.mb_sampled_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 0.5,
                "B": 0.3,
                "C": 0.2,
            },
            mbsamplingRate=3500,
            terminationTime=(1, "s"),
            constantSpecies=["B", "C"],
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

    @patch("rmgpy.rmg.input.logging")
    def test_mb_sampled_reactor_mole_fractions_normalize_1(self, mock_logging):
        """Test that MBSampledReactor mole fractions are normalized properly"""
        inp.mb_sampled_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            mbsamplingRate=3500,
            terminationTime=(1, "s"),
            constantSpecies=["B", "C"],
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial mole fractions do not sum to one; normalizing.")
