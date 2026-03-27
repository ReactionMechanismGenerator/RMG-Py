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

"""Tests for the automated library selection module."""

import os
import unittest
from unittest.mock import MagicMock

from rmgpy import settings
from rmgpy.data.libraries import (
    AUTO,
    PAH_LIBS,
    ChemistryProfile,
    detect_chemistry,
    determine_chemistry_sets,
    determine_kinetics_families,
    expand_chemistry_sets,
    load_recommended_yml,
    merge_with_user_libraries,
    merge_with_user_reaction_libraries,
    resolve_auto_kinetics_families,
)
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import Species


def make_species(smiles, reactive=True):
    """Helper to create a Species from SMILES."""
    mol = Molecule().from_smiles(smiles)
    spec = Species(molecule=[mol])
    spec.reactive = reactive
    return spec


class MockReactor:
    """Mock reactor with a Quantity-like T attribute."""
    def __init__(self, T_value_si):
        self.T = MagicMock()
        self.T.value_si = T_value_si


class MockLiquidReactor(MockReactor):
    """Mock that also passes isinstance checks for LiquidReactor."""
    pass


class MockSurfaceReactor(MockReactor):
    """Mock that also passes isinstance checks for SurfaceReactor."""
    pass


class TestDetectChemistry(unittest.TestCase):
    """Tests for the detect_chemistry function."""

    def test_simple_hydrocarbon(self):
        """C/H species at low temperature, no solvent."""
        species = [make_species('[CH4]'), make_species('[H][H]')]
        # Methane SMILES: C, H2: [H][H]
        species = [make_species('C'), make_species('[H][H]')]
        reactor = MockReactor(500.0)
        profile = detect_chemistry(species, [reactor], solvent=None)
        self.assertIn('C', profile.elements_present)
        self.assertIn('H', profile.elements_present)
        self.assertTrue(profile.has_carbon)
        self.assertFalse(profile.has_nitrogen)
        self.assertFalse(profile.has_sulfur)
        self.assertFalse(profile.has_oxygen)
        self.assertFalse(profile.has_liquid)
        self.assertFalse(profile.has_surface)
        self.assertAlmostEqual(profile.max_temperature, 500.0)

    def test_oxygenated_fuel(self):
        """Ethanol — contains O but no explicit O2."""
        species = [make_species('CCO'), make_species('N#N')]
        reactor = MockReactor(1000.0)
        profile = detect_chemistry(species, [reactor], solvent=None)
        self.assertTrue(profile.has_oxygen)
        self.assertTrue(profile.has_carbon)
        self.assertTrue(profile.has_nitrogen)
        self.assertAlmostEqual(profile.max_temperature, 1000.0)

    def test_O2_sets_has_oxygen(self):
        """O2 species should set has_oxygen via element detection."""
        species = [make_species('C'), make_species('[O][O]')]
        reactor = MockReactor(1500.0)
        profile = detect_chemistry(species, [reactor], solvent=None)
        self.assertTrue(profile.has_oxygen)

    def test_nitrogen_species(self):
        """NH3 should trigger nitrogen detection."""
        species = [make_species('N')]
        reactor = MockReactor(500.0)
        profile = detect_chemistry(species, [reactor], solvent=None)
        self.assertTrue(profile.has_nitrogen)

    def test_sulfur_species(self):
        """H2S should trigger sulfur detection."""
        species = [make_species('S')]  # H2S
        reactor = MockReactor(500.0)
        profile = detect_chemistry(species, [reactor], solvent=None)
        self.assertTrue(profile.has_sulfur)

    def test_solvent_sets_liquid(self):
        """Specifying a solvent should set has_liquid."""
        species = [make_species('C')]
        reactor = MockReactor(300.0)
        profile = detect_chemistry(species, [reactor], solvent='water')
        self.assertTrue(profile.has_liquid)

    def test_multiple_reactors_max_temperature(self):
        """Max temperature across reactors should be used."""
        species = [make_species('C')]
        reactor1 = MockReactor(500.0)
        reactor2 = MockReactor(1200.0)
        profile = detect_chemistry(species, [reactor1, reactor2], solvent=None)
        self.assertAlmostEqual(profile.max_temperature, 1200.0)

    def test_halogen_detection(self):
        """Chloromethane should trigger halogens."""
        species = [make_species('CCl')]
        reactor = MockReactor(500.0)
        profile = detect_chemistry(species, [reactor], solvent=None)
        self.assertTrue(profile.has_halogens)
        self.assertIn('Cl', profile.elements_present)


class TestDetermineChemistrySets(unittest.TestCase):
    """Tests for the determine_chemistry_sets function."""

    def test_primary_always_included(self):
        """Primary should always be in the returned sets."""
        profile = ChemistryProfile()
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertEqual(sets, ['primary'])

    def test_nitrogen_trigger(self):
        profile = ChemistryProfile(has_nitrogen=True)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('nitrogen', sets)

    def test_sulfur_trigger(self):
        profile = ChemistryProfile(has_sulfur=True)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('sulfur', sets)

    def test_combustion_trigger_with_oxygen(self):
        """Oxygen-containing species should trigger combustion."""
        profile = ChemistryProfile(has_oxygen=True)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('combustion', sets)

    def test_ch_pyrolysis_core_high_T_carbon(self):
        """C + high T should trigger CH_pyrolysis_core."""
        profile = ChemistryProfile(has_carbon=True, max_temperature=900.0)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('CH_pyrolysis_core', sets)

    def test_ch_pyrolysis_core_not_triggered_low_T(self):
        """C + low T should NOT trigger CH_pyrolysis_core."""
        profile = ChemistryProfile(has_carbon=True, max_temperature=500.0)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertNotIn('CH_pyrolysis_core', sets)

    def test_pah_no_oxygen_auto_included(self):
        """C + high T + no O -> PAH_formation auto-included."""
        profile = ChemistryProfile(has_carbon=True, max_temperature=1000.0)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('PAH_formation', sets)
        self.assertIn('CH_pyrolysis_core', sets)

    def test_pah_with_oxygen_not_auto_included(self):
        """C + high T + O present -> PAH_formation NOT auto-included."""
        profile = ChemistryProfile(has_carbon=True, has_oxygen=True, max_temperature=1000.0)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('CH_pyrolysis_core', sets)
        self.assertNotIn('PAH_formation', sets)

    def test_pah_with_oxygen_and_keyword(self):
        """C + high T + O + <PAH_libs> keyword -> PAH_formation included."""
        profile = ChemistryProfile(has_carbon=True, has_oxygen=True, max_temperature=1000.0)
        sets = determine_chemistry_sets(profile, pah_libs_requested=True)
        self.assertIn('CH_pyrolysis_core', sets)
        self.assertIn('PAH_formation', sets)

    def test_oxygenated_fuel_pyrolysis_with_keyword(self):
        """Oxygenated fuel at high T with <PAH_libs> gets combustion + core + PAH."""
        profile = ChemistryProfile(has_oxygen=True, has_carbon=True, max_temperature=1000.0)
        sets = determine_chemistry_sets(profile, pah_libs_requested=True)
        self.assertIn('combustion', sets)
        self.assertIn('CH_pyrolysis_core', sets)
        self.assertIn('PAH_formation', sets)

    def test_liquid_oxidation_trigger(self):
        """Liquid + oxygen should trigger liquid_oxidation."""
        profile = ChemistryProfile(has_liquid=True, has_oxygen=True)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('liquid_oxidation', sets)

    def test_surface_trigger(self):
        profile = ChemistryProfile(has_surface=True)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('surface', sets)

    def test_halogens_trigger(self):
        profile = ChemistryProfile(has_halogens=True)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('halogens', sets)

    def test_lithium_trigger(self):
        profile = ChemistryProfile(has_metal=True)
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn('lithium', sets)

    def test_full_combo(self):
        """All flags on."""
        profile = ChemistryProfile(
            has_nitrogen=True,
            has_sulfur=True,
            has_oxygen=True,
            has_carbon=True,
            has_halogens=True,
            has_metal=True,
            has_surface=True,
            has_liquid=True,
            max_temperature=1200.0,
        )
        sets = determine_chemistry_sets(profile, pah_libs_requested=True)
        for expected in ['primary', 'nitrogen', 'sulfur', 'combustion',
                         'CH_pyrolysis_core', 'PAH_formation', 'liquid_oxidation',
                         'surface', 'halogens', 'lithium']:
            self.assertIn(expected, sets)


class TestDetermineKineticsFamilies(unittest.TestCase):
    """Tests for the determine_kinetics_families function."""

    def test_default_always(self):
        profile = ChemistryProfile()
        sets = determine_kinetics_families(profile)
        self.assertEqual(sets, ['default'])

    def test_liquid_peroxide(self):
        profile = ChemistryProfile(has_liquid=True, has_oxygen=True)
        sets = determine_kinetics_families(profile)
        self.assertIn('liquid_peroxide', sets)

    def test_surface_families(self):
        profile = ChemistryProfile(has_surface=True)
        sets = determine_kinetics_families(profile)
        self.assertIn('surface', sets)

    def test_halogen_families(self):
        profile = ChemistryProfile(has_halogens=True)
        sets = determine_kinetics_families(profile)
        self.assertIn('halogens', sets)

    def test_electrochem_families(self):
        profile = ChemistryProfile(has_metal=True)
        sets = determine_kinetics_families(profile)
        self.assertIn('electrochem', sets)


class TestExpandChemistrySets(unittest.TestCase):
    """Tests for the expand_chemistry_sets function."""

    def setUp(self):
        self.db_dir = settings['database.directory']
        yml_path = os.path.join(self.db_dir, 'recommended_libraries.yml')
        if not os.path.isfile(yml_path):
            self.skipTest('recommended_libraries.yml not found in RMG-database')
        self.recommended = load_recommended_yml(self.db_dir)

    def test_primary_expansion(self):
        thermo, kinetics, transport, seeds = expand_chemistry_sets(self.recommended, ['primary'])
        self.assertIn('primaryThermoLibrary', thermo)
        self.assertIn('BurkeH2O2', thermo)
        self.assertEqual(kinetics, [])  # primaryH2O2 is a seed, not a kinetics lib
        self.assertIn('primaryH2O2', seeds)
        self.assertIn('PrimaryTransportLibrary', transport)

    def test_nitrogen_expansion(self):
        thermo, kinetics, transport, seeds = expand_chemistry_sets(self.recommended, ['nitrogen'])
        self.assertIn('NH3', thermo)
        self.assertIn('primaryNitrogenLibrary', kinetics)

    def test_ch_pyrolysis_core_expansion(self):
        thermo, kinetics, transport, seeds = expand_chemistry_sets(
            self.recommended, ['CH_pyrolysis_core']
        )
        self.assertIn('CurranPentane', thermo)
        self.assertIn('C2H2_init', kinetics)

    def test_pah_formation_expansion(self):
        thermo, kinetics, transport, seeds = expand_chemistry_sets(
            self.recommended, ['PAH_formation']
        )
        self.assertIn('naphthalene_H', thermo)
        self.assertIn('Mebel_C6H5_C2H2', kinetics)
        self.assertIn('Aromatics_high_pressure/C10H9_1', kinetics)

    def test_no_duplicates_across_sets(self):
        """Expanding multiple sets should not produce duplicates."""
        thermo, kinetics, transport, seeds = expand_chemistry_sets(
            self.recommended, ['primary', 'nitrogen', 'combustion', 'CH_pyrolysis_core']
        )
        self.assertEqual(len(thermo), len(set(thermo)))
        self.assertEqual(len(kinetics), len(set(kinetics)))
        self.assertEqual(len(seeds), len(set(seeds)))


class TestMergeWithUserLibraries(unittest.TestCase):
    """Tests for the merge logic."""

    def test_auto_string_replaces(self):
        result = merge_with_user_libraries('auto', ['lib1', 'lib2'])
        self.assertEqual(result, ['lib1', 'lib2'])

    def test_no_auto_in_list(self):
        result = merge_with_user_libraries(['userA', 'userB'], ['lib1', 'lib2'])
        self.assertEqual(result, ['userA', 'userB'])

    def test_auto_token_in_list_middle(self):
        result = merge_with_user_libraries(['userA', 'auto', 'userB'], ['lib1', 'lib2'])
        self.assertEqual(result, ['userA', 'lib1', 'lib2', 'userB'])

    def test_auto_token_at_start(self):
        result = merge_with_user_libraries(['auto', 'userB'], ['lib1', 'lib2'])
        self.assertEqual(result, ['lib1', 'lib2', 'userB'])

    def test_auto_token_at_end(self):
        result = merge_with_user_libraries(['userA', 'auto'], ['lib1', 'lib2'])
        self.assertEqual(result, ['userA', 'lib1', 'lib2'])

    def test_dedup_user_and_auto(self):
        """User-specified lib should not be duplicated from auto."""
        result = merge_with_user_libraries(['userA', 'auto'], ['userA', 'lib1', 'lib2'])
        self.assertEqual(result, ['userA', 'lib1', 'lib2'])

    def test_none_passthrough(self):
        result = merge_with_user_libraries(None, ['lib1'])
        self.assertIsNone(result)

    def test_empty_list_passthrough(self):
        result = merge_with_user_libraries([], ['lib1'])
        self.assertEqual(result, [])

    def test_pah_libs_token_stripped(self):
        """<PAH_libs> should be silently removed from the list."""
        result = merge_with_user_libraries(['userA', '<PAH_libs>', 'auto'], ['lib1'])
        self.assertEqual(result, ['userA', 'lib1'])
        self.assertNotIn('<PAH_libs>', result)

    def test_pah_libs_only_no_auto(self):
        """<PAH_libs> without 'auto' — stripped, user libs returned as-is."""
        result = merge_with_user_libraries(['userA', '<PAH_libs>'], ['lib1'])
        self.assertEqual(result, ['userA'])


class TestMergeWithUserReactionLibraries(unittest.TestCase):
    """Tests for reaction library merge (tuple format)."""

    def test_auto_string(self):
        result = merge_with_user_reaction_libraries('auto', ['lib1', 'lib2'])
        self.assertEqual(result, [('lib1', False), ('lib2', False)])

    def test_auto_token_in_list(self):
        user = [('userA', True), 'auto']
        result = merge_with_user_reaction_libraries(user, ['lib1', 'lib2'])
        self.assertEqual(result, [('userA', True), ('lib1', False), ('lib2', False)])

    def test_dedup_tuple_format(self):
        user = [('lib1', True), 'auto']
        result = merge_with_user_reaction_libraries(user, ['lib1', 'lib2'])
        # lib1 already specified by user, should not duplicate
        self.assertEqual(result, [('lib1', True), ('lib2', False)])

    def test_pah_libs_stripped_from_reaction_libs(self):
        user = [('userA', True), '<PAH_libs>', 'auto']
        result = merge_with_user_reaction_libraries(user, ['lib1'])
        self.assertEqual(result, [('userA', True), ('lib1', False)])


class TestResolveAutoKineticsFamilies(unittest.TestCase):
    """Tests for the kinetics families resolution from recommended.py."""

    def setUp(self):
        self.db_dir = settings['database.directory']
        rec_path = os.path.join(self.db_dir, 'kinetics', 'families', 'recommended.py')
        if not os.path.isfile(rec_path):
            self.skipTest('recommended.py not found in RMG-database')

    def test_default_set(self):
        families = resolve_auto_kinetics_families(['default'], self.db_dir)
        self.assertIn('H_Abstraction', families)
        self.assertIn('R_Recombination', families)

    def test_default_plus_surface(self):
        families = resolve_auto_kinetics_families(['default', 'surface'], self.db_dir)
        self.assertIn('H_Abstraction', families)
        self.assertIn('Surface_Adsorption_Single', families)

    def test_no_duplicates(self):
        families = resolve_auto_kinetics_families(['default', 'surface', 'halogens'], self.db_dir)
        self.assertEqual(len(families), len(set(families)))


if __name__ == '__main__':
    unittest.main()
