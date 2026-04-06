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


class TestWriteInputFile:
    """
    Contains unit test for writing input files for each of the reactor types:

        'simpleReactor': simple_reactor, ✅
        'constantVIdealGasReactor' : constant_V_ideal_gas_reactor, ✅
        'constantTPIdealGasReactor' : constant_TP_ideal_gas_reactor, ✅
        'liquidSurfaceReactor' : liquid_cat_reactor, ✅
        'constantTVLiquidReactor': constant_T_V_liquid_reactor,
        'liquidReactor': liquid_reactor, ✅
        'surfaceReactor': surface_reactor, ✅
        'mbsampledReactor': mb_sampled_reactor,

    """
    def setup_method(self):
        """This method is run before every test in this class"""
        global rmg
        rmg.reaction_systems = []

    def test_write_superminimal_input(self):
        """
        Test that we can write superminimal input file and read it back in with the same values.
        """

        superminimal_input_file = '../../../examples/rmg/superminimal/input.py'
        superminimal_output_file = 'temp_superminimal_input.py'

        rmg = RMG()
        inp.read_input_file(superminimal_input_file, rmg)

        # read a bunch of values in from input file to check they are the same after writing
        T = rmg.reaction_systems[0].T.value_si
        P = rmg.reaction_systems[0].P.value_si
        initialMoleFractions = {k.label: v for k, v in rmg.reaction_systems[0].initial_mole_fractions.items()}
        for term in rmg.reaction_systems[0].termination:
            if hasattr(term, 'time'):
                termination_time = term.time.value_si
            elif hasattr(term, 'conversion'):
                termination_conversion = term.conversion
                termination_converstion_species = term.species.label

        inp.save_input_file(superminimal_output_file, rmg)
        # read it back in and confirm all the values match
        rmg1 = RMG()
        inp.read_input_file(superminimal_output_file, rmg1)
        assert rmg1.reaction_systems[0].T.value_si == T
        assert rmg1.reaction_systems[0].P.value_si == P
        output_mol_fractions = {k.label: v for k, v in rmg1.reaction_systems[0].initial_mole_fractions.items()}
        assert output_mol_fractions == initialMoleFractions
        for term in rmg1.reaction_systems[0].termination:
            if hasattr(term, 'time'):
                assert term.time.value_si == termination_time
            elif hasattr(term, 'conversion'):
                assert term.conversion == termination_conversion
                assert term.species.label == termination_converstion_species

        # clean up
        import os
        os.remove(superminimal_output_file)

    @pytest.mark.skip(reason="Slow test that runs a full RMG job")
    def test_write_superminimal_and_run(self):
        """
        Test that we can write superminimal input file and then run RMG without errors
        """
        import os
        import shutil

        superminimal_input_file = '../../../examples/rmg/superminimal/input.py'
        new_run_dir = 'temp_superminimal_run'
        os.makedirs(new_run_dir, exist_ok=True)
        superminimal_output_file = os.path.join(new_run_dir, 'temp_superminimal_input.py')

        rmg = RMG()
        inp.read_input_file(superminimal_input_file, rmg)
        inp.save_input_file(superminimal_output_file, rmg)

        # run RMG with the new input file
        import subprocess
        subprocess.run(['python', '../../../rmg.py', superminimal_output_file], check=True)

        # clean up
        shutil.rmtree(new_run_dir)

    def test_write_min_surf_input(self):
        """
        Test that we can write the minimal surface input file and read it back in with the same values.
        """

        min_surf_input_file = '../../../examples/rmg/minimal_surface/input.py'
        min_surf_output_file = 'temp_min_surf_input.py'

        rmg = RMG()
        inp.read_input_file(min_surf_input_file, rmg)

        # read a bunch of values in from input file to check they are the same after writing
        T = rmg.reaction_systems[0].T.value_si
        P = rmg.reaction_systems[0].P_initial.value_si
        initialMoleFractions = {k.label: v for k, v in rmg.reaction_systems[0].initial_gas_mole_fractions.items()}
        initialSurfaceCoverages = {k.label: v for k, v in rmg.reaction_systems[0].initial_surface_coverages.items()}
        for term in rmg.reaction_systems[0].termination:
            if hasattr(term, 'time'):
                termination_time = term.time.value_si
            elif hasattr(term, 'conversion'):
                termination_conversion = term.conversion
                termination_converstion_species = term.species.label
            elif hasattr(term, 'ratio'):
                termination_ratio = term.ratio

        binding_energies = {k: v.value_si for k, v in rmg.binding_energies.items()}
        surface_site_density = rmg.surface_site_density.value_si

        inp.save_input_file(min_surf_output_file, rmg)
        # read it back in and confirm all the values match
        rmg1 = RMG()
        inp.read_input_file(min_surf_output_file, rmg1)
        assert rmg1.reaction_systems[0].T.value_si == T
        assert rmg1.reaction_systems[0].P_initial.value_si == P
        output_mol_fractions = {k.label: v for k, v in rmg1.reaction_systems[0].initial_gas_mole_fractions.items()}
        assert output_mol_fractions == initialMoleFractions
        output_surface_coverages = {k.label: v for k, v in rmg1.reaction_systems[0].initial_surface_coverages.items()}
        assert output_surface_coverages == initialSurfaceCoverages

        output_binding_energies = {k: v.value_si for k, v in rmg1.binding_energies.items()}
        assert output_binding_energies == binding_energies

        assert rmg1.surface_site_density.value_si == surface_site_density

        for term in rmg1.reaction_systems[0].termination:
            if hasattr(term, 'time'):
                assert term.time.value_si == termination_time
            elif hasattr(term, 'conversion'):
                assert term.conversion == termination_conversion
                assert term.species.label == termination_converstion_species
            elif hasattr(term, 'ratio'):
                assert term.ratio == termination_ratio

        # clean up
        import os
        os.remove(min_surf_output_file)

    @pytest.mark.skip(reason="Slow test that runs a full RMG job")
    def test_write_min_surf_and_run(self):
        """
        Test that we can write minimal surface input file and then run RMG without errors
        """
        import os
        import shutil

        min_surf_input_file = '../../../examples/rmg/minimal_surface/input.py'
        new_run_dir = 'temp_min_surf_run'
        os.makedirs(new_run_dir, exist_ok=True)
        min_surf_output_file = os.path.join(new_run_dir, 'temp_min_surf_input.py')

        rmg = RMG()
        inp.read_input_file(min_surf_input_file, rmg)
        inp.save_input_file(min_surf_output_file, rmg)

        # run RMG with the new input file
        import subprocess
        subprocess.run(['python', '../../../rmg.py', min_surf_output_file], check=True)

        # clean up
        shutil.rmtree(new_run_dir)

    @pytest.mark.skip(reason="Slow because it has to compile Julia")
    def test_write_liquid_cat_input(self):
        """
        Test that we can write liquid catalyst input file and read it back in with the same values.
        """

        liquid_cat_input_file = '../../../examples/rmg/liquid_cat/input.py'
        liquid_cat_output_file = 'temp_liquid_cat_input.py'

        rmg = RMG()
        inp.read_input_file(liquid_cat_input_file, rmg)

        # read a bunch of values in from input file to check they are the same after writing
        T = rmg.reaction_systems[0].T.value_si
        tf = rmg.reaction_systems[0].tf
        liquid_init_conditions = {k: v for k, v in rmg.reaction_systems[0].initial_conditions['liquid'].items()}
        surf_init_conditions = {k: v for k, v in rmg.reaction_systems[0].initial_conditions['surface'].items()}

        termination_species = rmg.reaction_systems[0].terminations[0][0].label
        termination_conversion = rmg.reaction_systems[0].terminations[0][1]

        inp.save_input_file(liquid_cat_output_file, rmg)
        # read it back in and confirm all the values match
        rmg1 = RMG()
        inp.read_input_file(liquid_cat_output_file, rmg1)
        assert rmg1.reaction_systems[0].T.value_si == T
        assert rmg1.reaction_systems[0].tf == tf

        liquid_init_conditions_output = {k: v for k, v in rmg1.reaction_systems[0].initial_conditions['liquid'].items()}
        assert pytest.approx(liquid_init_conditions_output) == liquid_init_conditions
        surf_init_conditions_output = {k: v for k, v in rmg1.reaction_systems[0].initial_conditions['surface'].items()}
        assert pytest.approx(surf_init_conditions_output) == surf_init_conditions

        assert rmg1.reaction_systems[0].terminations[0][0].label == termination_species
        assert rmg1.reaction_systems[0].terminations[0][1] == termination_conversion

        # clean up
        import os
        os.remove(liquid_cat_output_file)

    @pytest.mark.skip(reason="Slow test that runs a full RMG job")
    def test_write_liquid_cat_and_run(self):
        """
        Test that we can write liquid catalyst input file and then run RMG without errors
        """
        import os
        import shutil

        liquid_cat_input_file = '../../../examples/rmg/liquid_cat/input.py'
        new_run_dir = 'temp_liquid_cat_run'
        os.makedirs(new_run_dir, exist_ok=True)
        liquid_cat_output_file = os.path.join(new_run_dir, 'temp_liquid_cat_input.py')

        rmg = RMG()
        inp.read_input_file(liquid_cat_input_file, rmg)
        inp.save_input_file(liquid_cat_output_file, rmg)

        # run RMG with the new input file
        import subprocess
        subprocess.run(['python', '../../../rmg.py', '-t', '00:00:01:30', liquid_cat_output_file], check=True)

        # clean up
        shutil.rmtree(new_run_dir)

    def test_write_liquid_input(self):
        """
        Test that we can write the liquid reactor input file and read it back in with the same values.
        """

        liquid_input_file = '../../../examples/rmg/liquid_phase/input.py'
        liquid_output_file = 'temp_liquid_input.py'

        rmg = RMG()
        inp.read_input_file(liquid_input_file, rmg)

        # read a bunch of values in from input file to check they are the same after writing
        T = rmg.reaction_systems[0].T.value_si
        P = rmg.reaction_systems[0].P.value_si
        initialConcentrations = {k.label: v for k, v in rmg.reaction_systems[0].initial_concentrations.items()}
        for term in rmg.reaction_systems[0].termination:
            if hasattr(term, 'time'):
                termination_time = term.time.value_si
        solvent = rmg.solvent

        inp.save_input_file(liquid_output_file, rmg)
        # read it back in and confirm all the values match
        rmg1 = RMG()
        inp.read_input_file(liquid_output_file, rmg1)
        assert rmg1.reaction_systems[0].T.value_si == T
        assert rmg1.reaction_systems[0].P.value_si == P
        output_concentrations = {k.label: v for k, v in rmg1.reaction_systems[0].initial_concentrations.items()}
        assert pytest.approx(output_concentrations) == initialConcentrations
        for term in rmg1.reaction_systems[0].termination:
            if hasattr(term, 'time'):
                assert term.time.value_si == termination_time
        assert rmg1.solvent == solvent

        # clean up
        import os
        os.remove(liquid_output_file)

    @pytest.mark.skip(reason="Slow test that runs a full RMG job")
    def test_write_liquid_and_run(self):
        """
        Test that we can write liquid reactor input file and then run RMG without errors
        """
        import os
        import shutil

        liquid_input_file = '../../../examples/rmg/liquid_phase/input.py'
        new_run_dir = 'temp_liquid_run'
        os.makedirs(new_run_dir, exist_ok=True)
        liquid_output_file = os.path.join(new_run_dir, 'temp_liquid_input.py')

        rmg = RMG()
        inp.read_input_file(liquid_input_file, rmg)
        inp.save_input_file(liquid_output_file, rmg)

        # run RMG with the new input file
        import subprocess
        subprocess.run(['python', '../../../rmg.py', '-t', '00:00:01:30', liquid_output_file], check=True)

        # clean up
        shutil.rmtree(new_run_dir)

    @pytest.mark.skip(reason="Slow because it has to compile Julia")
    def test_write_constantVIdealGasReactor(self):
        """
        Test that we can write constant volume ideal gas reactor input file and read it back in with the same values.
        """

        rms_constant_V_input_file = '../../../examples/rmg/rms_constant_V/input.py'
        rms_constant_V_output_file = 'temp_rms_constant_V_input.py'

        rmg = RMG()
        inp.read_input_file(rms_constant_V_input_file, rmg)

        # read a bunch of values in from input file to check they are the same after writing
        T = rmg.reaction_systems[0].T.value_si
        P = rmg.reaction_systems[0].P.value_si
        tf = rmg.reaction_systems[0].tf
        init_conditions = {k: v for k, v in rmg.reaction_systems[0].initial_conditions.items()}

        termination_species = rmg.reaction_systems[0].terminations[0][0].label
        termination_conversion = rmg.reaction_systems[0].terminations[0][1]
        termination_time = rmg.reaction_systems[0].terminations[1].time

        inp.save_input_file(rms_constant_V_output_file, rmg)
        # read it back in and confirm all the values match
        rmg1 = RMG()
        inp.read_input_file(rms_constant_V_output_file, rmg1)
        assert rmg1.reaction_systems[0].T.value_si == T
        assert rmg1.reaction_systems[0].P.value_si == P
        assert rmg1.reaction_systems[0].tf == tf

        new_init_conditions = {k: v for k, v in rmg1.reaction_systems[0].initial_conditions.items()}
        assert pytest.approx(new_init_conditions) == init_conditions

        assert rmg1.reaction_systems[0].terminations[0][0].label == termination_species
        assert rmg1.reaction_systems[0].terminations[0][1] == termination_conversion
        assert rmg1.reaction_systems[0].terminations[1].time == termination_time

        # clean up
        import os
        os.remove(rms_constant_V_output_file)

    @pytest.mark.skip(reason="Slow test that runs a full RMG job")
    def test_write_constantVIdealGasReactor_and_run(self):
        """
        Test that we can write constant volume ideal gas reactor input file and then run RMG without errors
        """
        import os
        import shutil

        constant_V_input_file = '../../../examples/rmg/rms_constant_V/input.py'
        new_run_dir = 'temp_constant_V_run'
        os.makedirs(new_run_dir, exist_ok=True)
        constant_V_output_file = os.path.join(new_run_dir, 'temp_constant_V_input.py')

        rmg = RMG()
        inp.read_input_file(constant_V_input_file, rmg)
        inp.save_input_file(constant_V_output_file, rmg)

        # run RMG with the new input file
        import subprocess
        subprocess.run(['python', '../../../rmg.py', '-t', '00:00:01:30', constant_V_output_file], check=True)

        # clean up
        shutil.rmtree(new_run_dir)

    @pytest.mark.skip(reason="Slow because it has to compile Julia")
    def test_write_constantTPdealGasReactor(self):
        """
        Test that we can write constant TP ideal gas reactor input file and read it back in with the same values.
        """

        rms_constant_TP_input_file = '../../../examples/rmg/nox_transitory_edge/input.py'
        rms_constant_TP_output_file = 'temp_constant_TP_input.py'

        rmg = RMG()
        inp.read_input_file(rms_constant_TP_input_file, rmg)

        # read a bunch of values in from input file to check they are the same after writing
        T = rmg.reaction_systems[0].T.value_si
        P = rmg.reaction_systems[0].P.value_si
        tf = rmg.reaction_systems[0].tf
        init_conditions = {k: v for k, v in rmg.reaction_systems[0].initial_conditions.items()}

        termination_species = rmg.reaction_systems[0].terminations[0][0].label
        termination_conversion = rmg.reaction_systems[0].terminations[0][1]
        termination_time = rmg.reaction_systems[0].terminations[1].time

        inp.save_input_file(rms_constant_TP_output_file, rmg)
        # read it back in and confirm all the values match
        rmg1 = RMG()
        inp.read_input_file(rms_constant_TP_output_file, rmg1)
        assert rmg1.reaction_systems[0].T.value_si == T
        assert rmg1.reaction_systems[0].P.value_si == P
        assert rmg1.reaction_systems[0].tf == tf

        new_init_conditions = {k: v for k, v in rmg1.reaction_systems[0].initial_conditions.items()}
        assert pytest.approx(new_init_conditions, rel=1e-4) == init_conditions

        assert rmg1.reaction_systems[0].terminations[0][0].label == termination_species
        assert rmg1.reaction_systems[0].terminations[0][1] == termination_conversion
        assert rmg1.reaction_systems[0].terminations[1].time == termination_time

        # clean up
        import os
        os.remove(rms_constant_TP_output_file)

    @pytest.mark.skip(reason="Slow test that runs a full RMG job")
    def test_write_constantTPIdealGasReactor_and_run(self):
        """
        Test that we can write constant TP ideal gas reactor input file and then run RMG without errors
        """
        import os
        import shutil

        constant_TP_input_file = '../../../examples/rmg/nox_transitory_edge/input.py'
        new_run_dir = 'temp_constant_TP_run'
        os.makedirs(new_run_dir, exist_ok=True)
        constant_TP_output_file = os.path.join(new_run_dir, 'temp_constant_TP_input.py')

        rmg = RMG()
        inp.read_input_file(constant_TP_input_file, rmg)
        inp.save_input_file(constant_TP_output_file, rmg)

        # run RMG with the new input file
        import subprocess
        subprocess.run(['python', '../../../rmg.py', '-t', '00:00:01:30', constant_TP_output_file], check=True)

        # clean up
        shutil.rmtree(new_run_dir)
