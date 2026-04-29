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

import rmgpy.rmg.input as inp
from rmgpy.exceptions import InputError
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.ml.estimator import ADMONITION

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
        Test that reaction libraries given as plain strings are stored as strings.
        """
        global rmg
        inp.database(reactionLibraries=["test"])
        assert rmg.reaction_libraries == ["test"]
        assert "test" not in rmg.reaction_libraries_output_edge

    def test_importing_database_reaction_libraries_from_false_tuple(self):
        """
        Test that (name, False) tuples are stored as plain strings without output_edge.
        """
        global rmg
        inp.database(reactionLibraries=[("test", False)])
        assert rmg.reaction_libraries == ["test"]
        assert "test" not in rmg.reaction_libraries_output_edge

    def test_importing_database_reaction_libraries_from_true_tuple(self):
        """
        Test that (name, True) tuples are stored as plain strings with the name in output_edge.
        """
        global rmg
        inp.database(reactionLibraries=[("test", True)])
        assert rmg.reaction_libraries == ["test"]
        assert "test" in rmg.reaction_libraries_output_edge

class TestInputDatabaseAutoSelection:
    """Tests for 'auto' and '<PAH_libs>' token handling in database()."""

    def test_thermo_auto_string(self):
        global rmg
        inp.database(thermoLibraries='auto')
        assert rmg.thermo_libraries == 'auto'

    def test_thermo_auto_in_list(self):
        global rmg
        inp.database(thermoLibraries=['myLib', 'auto'])
        assert rmg.thermo_libraries == ['myLib', 'auto']

    def test_thermo_pah_libs_in_list(self):
        global rmg
        inp.database(thermoLibraries=['auto', '<PAH_libs>'])
        assert rmg.thermo_libraries == ['auto', '<PAH_libs>']

    def test_transport_auto_string(self):
        global rmg
        inp.database(transportLibraries='auto')
        assert rmg.transport_libraries == 'auto'

    def test_reaction_libraries_auto_string(self):
        global rmg
        inp.database(reactionLibraries='auto')
        assert rmg.reaction_libraries == 'auto'

    def test_reaction_libraries_auto_in_list(self):
        global rmg
        inp.database(reactionLibraries=['auto', '<PAH_libs>'])
        assert 'auto' in rmg.reaction_libraries
        assert '<PAH_libs>' in rmg.reaction_libraries

    def test_reaction_libraries_mixed_auto_and_tuple(self):
        global rmg
        inp.database(reactionLibraries=[('myLib', True), 'auto'])
        assert rmg.reaction_libraries == ['myLib', 'auto']
        assert 'myLib' in rmg.reaction_libraries_output_edge

    def test_seed_mechanisms_auto(self):
        global rmg
        inp.database(seedMechanisms='auto')
        assert rmg.seed_mechanisms == 'auto'

    def test_kinetics_families_auto(self):
        global rmg
        inp.database(kineticsFamilies='auto')
        assert rmg.kinetics_families == 'auto'

    def test_kinetics_families_auto_with_exclusion(self):
        global rmg
        inp.database(kineticsFamilies=['!H_Abstraction', 'auto'])
        assert rmg.kinetics_families == ['!H_Abstraction', 'auto']

    def test_default_no_auto(self):
        """Without 'auto', fields should behave as before."""
        global rmg
        inp.database(
            thermoLibraries=['primaryThermoLibrary'],
            reactionLibraries=[('lib1', False)],
            seedMechanisms=[],
            transportLibraries=None,
            kineticsFamilies='default',
        )
        assert rmg.thermo_libraries == ['primaryThermoLibrary']
        assert rmg.reaction_libraries == ['lib1']
        assert 'lib1' not in rmg.reaction_libraries_output_edge
        assert rmg.seed_mechanisms == []
        assert rmg.transport_libraries is None
        assert rmg.kinetics_families == 'default'

    def test_pah_libs_standalone_thermo_raises(self):
        with pytest.raises(InputError):
            inp.database(thermoLibraries='<PAH_libs>')

    def test_pah_libs_standalone_reaction_raises(self):
        with pytest.raises(InputError):
            inp.database(reactionLibraries='<PAH_libs>')

    def test_pah_libs_standalone_transport_raises(self):
        with pytest.raises(InputError):
            inp.database(transportLibraries='<PAH_libs>')

    def test_pah_libs_standalone_seeds_raises(self):
        with pytest.raises(InputError):
            inp.database(seedMechanisms='<PAH_libs>')

    def test_pah_libs_tuple_in_reaction_libs_raises(self):
        with pytest.raises(InputError):
            inp.database(reactionLibraries=[('<PAH_libs>', True)])


@pytest.mark.skip(reason=ADMONITION)
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


class TestInputThermoCentralDatabase:
    """
    Contains unit tests rmgpy.rmg.input.thermo_central_database
    """

    def teardown_class(self):
        # remove the reactionLibraries value
        global rmg
        rmg.thermo_central_database = None

    def test_thermo_central_database(self):
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


class TestInputPressureDependence:
    """
    Contains unit tests for pressure dependence input, including completedNetworks
    """

    def setup_method(self):
        """This method is run before every test in this class"""
        global rmg
        # Reset the completed networks set before each test
        rmg.reaction_model = CoreEdgeReactionModel()
        rmg.reaction_model.completed_pdep_networks = set()

    def test_completed_networks_single(self):
        """Test that a single completedNetwork can be added via pressure_dependence"""
        global rmg
        
        inp.pressure_dependence(
            method='modified strong collision',
            temperatures=(300, 2000, 'K', 8),
            pressures=(0.01, 100, 'bar', 5),
            maximumGrainSize=(0.5, 'kcal/mol'),
            minimumNumberOfGrains=250,
            interpolation=('Chebyshev', 6, 4),
            maximumAtoms=16,
            completedNetworks=['CH2O2'],
        )
        
        # Check that the network was added
        assert len(rmg.reaction_model.completed_pdep_networks) == 1
        # The formula CH2O2 should be converted to a sorted tuple of elements
        expected_key = (('C', 1), ('H', 2), ('O', 2))
        assert expected_key in rmg.reaction_model.completed_pdep_networks

    def test_completed_networks_multiple(self):
        """Test that multiple completedNetworks can be added via pressure_dependence"""
        global rmg
        
        inp.pressure_dependence(
            method='modified strong collision',
            temperatures=(300, 2000, 'K', 8),
            pressures=(0.01, 100, 'bar', 5),
            maximumGrainSize=(0.5, 'kcal/mol'),
            minimumNumberOfGrains=250,
            interpolation=('Chebyshev', 6, 4),
            maximumAtoms=16,
            completedNetworks=['CH2O2', 'C2H6'],
        )
        
        # Check that both networks were added
        assert len(rmg.reaction_model.completed_pdep_networks) == 2
        expected_key1 = (('C', 1), ('H', 2), ('O', 2))
        expected_key2 = (('C', 2), ('H', 6))
        assert expected_key1 in rmg.reaction_model.completed_pdep_networks
        assert expected_key2 in rmg.reaction_model.completed_pdep_networks

    def test_completed_networks_none(self):
        """Test that pressure_dependence works without completedNetworks"""
        global rmg
        
        inp.pressure_dependence(
            method='modified strong collision',
            temperatures=(300, 2000, 'K', 8),
            pressures=(0.01, 100, 'bar', 5),
            maximumGrainSize=(0.5, 'kcal/mol'),
            minimumNumberOfGrains=250,
            interpolation=('Chebyshev', 6, 4),
            maximumAtoms=16,
        )
        
        # Check that no networks were added
        assert len(rmg.reaction_model.completed_pdep_networks) == 0
