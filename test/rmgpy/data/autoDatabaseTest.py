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

"""Tests for the automated library and family selection module (auto_database.py)."""

import logging
import os
import unittest

from rmgpy import settings
from rmgpy.data.auto_database import (
    AUTO,
    PAH_LIBS,
    ChemistryProfile,
    ChemistrySet,
    FamilySet,
    _get_reactor_max_temperature,
    _has_pah_libs_keyword,
    _log_lib_list,
    to_reaction_library_tuples,
    auto_select_libraries,
    detect_chemistry,
    determine_chemistry_sets,
    determine_kinetics_families,
    expand_chemistry_sets,
    load_recommended_yml,
    merge_with_user_libraries,
    resolve_auto_kinetics_families,
)
from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.rmg.main import RMG
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.surface import SurfaceReactor
from rmgpy.species import Species


def _simple_reactor(T_K: float) -> SimpleReactor:
    spc = Species().from_smiles('C')
    return SimpleReactor(T=Quantity(T_K, 'K'), P=Quantity(1e5, 'Pa'),
                         initial_mole_fractions={spc: 1.0}, n_sims=1, termination=[])


def _liquid_reactor(T_K: float) -> LiquidReactor:
    spc = Species().from_smiles('C')
    return LiquidReactor(T=Quantity(T_K, 'K'), initial_concentrations={spc: 0.1},
                         n_sims=1, termination=[])


def _surface_reactor(T_K: float) -> SurfaceReactor:
    gas_spc = Species().from_smiles('[H][H]')
    x_mol = Molecule().from_adjacency_list('1 X u0 p0 c0')
    surf_spc = Species(molecule=[x_mol])
    return SurfaceReactor(T=Quantity(T_K, 'K'), P_initial=Quantity(1e5, 'Pa'),
                          initial_gas_mole_fractions={gas_spc: 1.0},
                          initial_surface_coverages={surf_spc: 1.0},
                          surface_volume_ratio=Quantity(1e1, 'm^-1'),
                          surface_site_density=Quantity(2.72e-9, 'mol/cm^2'),
                          n_sims=1, termination=[])


def _rmg_with_auto(species_smiles, T_K, solvent=None, pah=False):
    """Build an RMG object with 'auto' library settings."""
    rmg = RMG()
    rmg.initial_species = [Species().from_smiles(s) for s in species_smiles]
    rmg.reaction_systems = [_simple_reactor(T_K)]
    rmg.solvent = solvent
    rmg.database_directory = settings['database.directory']
    rmg.thermo_libraries = AUTO
    rmg.reaction_libraries = [AUTO, PAH_LIBS] if pah else AUTO
    rmg.seed_mechanisms = AUTO
    rmg.transport_libraries = AUTO
    rmg.kinetics_families = AUTO
    rmg.reaction_libraries_output_edge = set()
    return rmg


class TestGetReactorMaxTemperature(unittest.TestCase):

    def test_simple_reactor(self):
        self.assertAlmostEqual(_get_reactor_max_temperature(_simple_reactor(1200.0)), 1200.0)

    def test_liquid_reactor(self):
        self.assertAlmostEqual(_get_reactor_max_temperature(_liquid_reactor(350.0)), 350.0)

    def test_surface_reactor(self):
        self.assertAlmostEqual(_get_reactor_max_temperature(_surface_reactor(600.0)), 600.0)

    def test_temperature_range(self):
        spc = Species().from_smiles('C')
        reactor = SimpleReactor(T=[Quantity(800, 'K'), Quantity(1500, 'K')],
                                P=Quantity(1e5, 'Pa'),
                                initial_mole_fractions={spc: 1.0}, n_sims=1, termination=[])
        self.assertAlmostEqual(_get_reactor_max_temperature(reactor), 1500.0)

    def test_no_T_attribute(self):
        self.assertIsNone(_get_reactor_max_temperature(object()))


class TestHasPahLibsKeyword(unittest.TestCase):

    def test_in_thermo(self):
        rmg = RMG()
        rmg.thermo_libraries = ['auto', PAH_LIBS]
        rmg.reaction_libraries = []
        rmg.seed_mechanisms = []
        rmg.transport_libraries = []
        self.assertTrue(_has_pah_libs_keyword(rmg))

    def test_in_reaction_libraries(self):
        rmg = RMG()
        rmg.thermo_libraries = []
        rmg.reaction_libraries = ['myLib', '<PAH_libs>']
        rmg.seed_mechanisms = []
        rmg.transport_libraries = []
        self.assertTrue(_has_pah_libs_keyword(rmg))

    def test_absent(self):
        rmg = RMG()
        rmg.thermo_libraries = ['auto']
        rmg.reaction_libraries = ['auto']
        rmg.seed_mechanisms = ['auto']
        rmg.transport_libraries = ['auto']
        self.assertFalse(_has_pah_libs_keyword(rmg))

    def test_string_field_not_iterated(self):
        rmg = RMG()
        rmg.thermo_libraries = 'auto'
        rmg.reaction_libraries = 'auto'
        rmg.seed_mechanisms = 'auto'
        rmg.transport_libraries = 'auto'
        self.assertFalse(_has_pah_libs_keyword(rmg))


class TestDetectChemistry(unittest.TestCase):

    def test_hydrocarbon(self):
        species = [Species().from_smiles('C'), Species().from_smiles('[H][H]')]
        profile = detect_chemistry(species, [_simple_reactor(500.0)], solvent=None)
        self.assertTrue(profile.has_carbon)
        self.assertFalse(profile.has_nitrogen)
        self.assertFalse(profile.has_oxygen)
        self.assertAlmostEqual(profile.max_temperature, 500.0)

    def test_oxygenated_fuel(self):
        species = [Species().from_smiles('CCO'), Species().from_smiles('N#N')]
        profile = detect_chemistry(species, [_simple_reactor(1000.0)], solvent=None)
        self.assertTrue(profile.has_oxygen)
        self.assertTrue(profile.has_carbon)
        self.assertTrue(profile.has_nitrogen)

    def test_O2(self):
        species = [Species().from_smiles('C'), Species().from_smiles('[O][O]')]
        profile = detect_chemistry(species, [_simple_reactor(1500.0)], solvent=None)
        self.assertTrue(profile.has_oxygen)

    def test_nitrogen(self):
        profile = detect_chemistry([Species().from_smiles('N')], [_simple_reactor(500.0)], solvent=None)
        self.assertTrue(profile.has_nitrogen)

    def test_sulfur(self):
        profile = detect_chemistry([Species().from_smiles('S')], [_simple_reactor(500.0)], solvent=None)
        self.assertTrue(profile.has_sulfur)

    def test_halogens(self):
        profile = detect_chemistry([Species().from_smiles('CCl')], [_simple_reactor(500.0)], solvent=None)
        self.assertTrue(profile.has_halogens)

    def test_liquid_reactor(self):
        profile = detect_chemistry([Species().from_smiles('C')], [_liquid_reactor(300.0)], solvent=None)
        self.assertTrue(profile.has_liquid)

    def test_surface_reactor(self):
        profile = detect_chemistry([Species().from_smiles('[H][H]')], [_surface_reactor(600.0)], solvent=None)
        self.assertTrue(profile.has_surface)

    def test_solvent(self):
        profile = detect_chemistry([Species().from_smiles('C')], [_simple_reactor(300.0)], solvent='water')
        self.assertTrue(profile.has_liquid)

    def test_max_T_multiple_reactors(self):
        profile = detect_chemistry([Species().from_smiles('C')],
                                   [_simple_reactor(500.0), _simple_reactor(1200.0)], solvent=None)
        self.assertAlmostEqual(profile.max_temperature, 1200.0)

    def test_surface_species(self):
        x_mol = Molecule().from_adjacency_list('1 X u0 p0 c0')
        profile = detect_chemistry([Species(molecule=[x_mol])], [_simple_reactor(500.0)], solvent=None)
        self.assertTrue(profile.has_surface)

    def test_mixed_reactors(self):
        profile = detect_chemistry([Species().from_smiles('C')],
                                   [_simple_reactor(800.0), _liquid_reactor(350.0)], solvent=None)
        self.assertTrue(profile.has_liquid)
        self.assertAlmostEqual(profile.max_temperature, 800.0)

    def test_no_reactors(self):
        profile = detect_chemistry([Species().from_smiles('CCO')], [], solvent=None)
        self.assertTrue(profile.has_oxygen)
        self.assertAlmostEqual(profile.max_temperature, 0.0)

    def test_nonreactive_species_ignored(self):
        reactive = Species().from_smiles('C')
        bath_gas = Species(reactive=False).from_smiles('N#N')
        profile = detect_chemistry([reactive, bath_gas], [_simple_reactor(500.0)], solvent=None)
        self.assertTrue(profile.has_carbon)
        self.assertFalse(profile.has_nitrogen)


class TestDetermineChemistrySets(unittest.TestCase):

    def test_primary_always(self):
        self.assertEqual(determine_chemistry_sets(ChemistryProfile(), False), [ChemistrySet.PRIMARY])

    def test_nitrogen(self):
        self.assertIn(ChemistrySet.NITROGEN, determine_chemistry_sets(ChemistryProfile(has_nitrogen=True), False))

    def test_sulfur(self):
        self.assertIn(ChemistrySet.SULFUR, determine_chemistry_sets(ChemistryProfile(has_sulfur=True), False))

    def test_oxidation(self):
        self.assertIn(ChemistrySet.OXIDATION, determine_chemistry_sets(ChemistryProfile(has_oxygen=True), False))

    def test_ch_pyrolysis_core(self):
        self.assertIn(ChemistrySet.CH_PYROLYSIS_CORE,
                      determine_chemistry_sets(ChemistryProfile(has_carbon=True, max_temperature=900.0), False))

    def test_ch_pyrolysis_core_low_T(self):
        self.assertNotIn(ChemistrySet.CH_PYROLYSIS_CORE,
                         determine_chemistry_sets(ChemistryProfile(has_carbon=True, max_temperature=500.0), False))

    def test_pah_no_oxygen(self):
        self.assertIn(ChemistrySet.PAH_FORMATION,
                      determine_chemistry_sets(ChemistryProfile(has_carbon=True, max_temperature=1000.0), False))

    def test_pah_blocked_by_oxygen(self):
        self.assertNotIn(ChemistrySet.PAH_FORMATION,
                         determine_chemistry_sets(
                             ChemistryProfile(has_carbon=True, has_oxygen=True, max_temperature=1000.0), False))

    def test_pah_with_keyword(self):
        self.assertIn(ChemistrySet.PAH_FORMATION,
                      determine_chemistry_sets(ChemistryProfile(has_carbon=True, has_oxygen=True, max_temperature=1000.0), True))

    def test_liquid_oxidation(self):
        self.assertIn(ChemistrySet.LIQUID_OXIDATION,
                      determine_chemistry_sets(ChemistryProfile(has_liquid=True, has_oxygen=True), False))

    def test_surface(self):
        sets = determine_chemistry_sets(ChemistryProfile(has_surface=True), False)
        self.assertIn(ChemistrySet.SURFACE, sets)
        self.assertNotIn(ChemistrySet.SURFACE_NITROGEN, sets)

    def test_surface_nitrogen(self):
        sets = determine_chemistry_sets(ChemistryProfile(has_surface=True, has_nitrogen=True), False)
        self.assertIn(ChemistrySet.SURFACE, sets)
        self.assertIn(ChemistrySet.SURFACE_NITROGEN, sets)

    def test_halogens(self):
        self.assertIn(ChemistrySet.HALOGENS,
                      determine_chemistry_sets(ChemistryProfile(has_halogens=True), False))

    def test_electrochem(self):
        self.assertIn(ChemistrySet.ELECTROCHEM,
                      determine_chemistry_sets(ChemistryProfile(has_electrochem=True), False))

    def test_ordering(self):
        profile = ChemistryProfile(has_nitrogen=True, has_sulfur=True, has_oxygen=True,
                                   has_carbon=True, max_temperature=1200.0)
        sets = determine_chemistry_sets(profile, False)
        self.assertEqual(sets[0], ChemistrySet.PRIMARY)
        self.assertLess(sets.index(ChemistrySet.NITROGEN), sets.index(ChemistrySet.SULFUR))

    def test_full_combo(self):
        profile = ChemistryProfile(
            has_nitrogen=True, has_sulfur=True, has_oxygen=True, has_carbon=True,
            has_halogens=True, has_electrochem=True, has_surface=True, has_liquid=True,
            max_temperature=1200.0)
        sets = determine_chemistry_sets(profile, True)
        for expected in ChemistrySet:
            self.assertIn(expected, sets)


class TestDetermineKineticsFamilies(unittest.TestCase):

    def test_default(self):
        self.assertEqual(determine_kinetics_families(ChemistryProfile()), [FamilySet.DEFAULT])

    def test_ch_pyrolysis(self):
        self.assertIn(FamilySet.CH_PYROLYSIS,
                      determine_kinetics_families(ChemistryProfile(has_carbon=True, max_temperature=1000.0)))

    def test_ch_pyrolysis_low_T(self):
        self.assertNotIn(FamilySet.CH_PYROLYSIS,
                         determine_kinetics_families(ChemistryProfile(has_carbon=True, max_temperature=500.0)))

    def test_liquid_peroxide(self):
        self.assertIn(FamilySet.LIQUID_PEROXIDE,
                      determine_kinetics_families(ChemistryProfile(has_liquid=True, has_oxygen=True)))

    def test_surface(self):
        self.assertIn(FamilySet.SURFACE,
                      determine_kinetics_families(ChemistryProfile(has_surface=True)))

    def test_halogens(self):
        self.assertIn(FamilySet.HALOGENS,
                      determine_kinetics_families(ChemistryProfile(has_halogens=True)))

    def test_electrochem(self):
        self.assertIn(FamilySet.ELECTROCHEM,
                      determine_kinetics_families(ChemistryProfile(has_electrochem=True)))


class TestLoadRecommendedYml(unittest.TestCase):

    def setUp(self):
        self.db_dir = settings['database.directory']
        if not os.path.isfile(os.path.join(self.db_dir, 'recommended_libraries.yml')):
            self.skipTest('recommended_libraries.yml not found')

    def test_loads(self):
        data = load_recommended_yml(self.db_dir)
        self.assertIsInstance(data, dict)
        self.assertIn('primary', data)

    def test_all_chemistry_sets_present(self):
        data = load_recommended_yml(self.db_dir)
        for s in ChemistrySet:
            self.assertIn(s.value, data)

    def test_required_keys(self):
        data = load_recommended_yml(self.db_dir)
        for name, entry in data.items():
            self.assertIn('thermo', entry, f'{name} missing thermo')
            self.assertIn('kinetics', entry, f'{name} missing kinetics')

    def test_invalid_path(self):
        with self.assertRaises(Exception):
            load_recommended_yml('/nonexistent')


class TestExpandChemistrySets(unittest.TestCase):

    def setUp(self):
        self.db_dir = settings['database.directory']
        if not os.path.isfile(os.path.join(self.db_dir, 'recommended_libraries.yml')):
            self.skipTest('recommended_libraries.yml not found')
        self.rec = load_recommended_yml(self.db_dir)

    def test_primary(self):
        thermo, _, transport, seeds = expand_chemistry_sets(self.rec, ['primary'])
        self.assertIn('primaryThermoLibrary', thermo)
        self.assertIn('primaryH2O2', seeds)
        self.assertIn('PrimaryTransportLibrary', transport)

    def test_nitrogen(self):
        thermo, kinetics, _, _ = expand_chemistry_sets(self.rec, ['nitrogen'])
        self.assertIn('NH3', thermo)
        self.assertIn('primaryNitrogenLibrary', kinetics)

    def test_ch_pyrolysis_core(self):
        thermo, kinetics, _, _ = expand_chemistry_sets(self.rec, ['CH_pyrolysis_core'])
        self.assertIn('Klippenstein_Glarborg2016', thermo)
        self.assertIn('C2H2_init', kinetics)

    def test_pah_formation(self):
        thermo, kinetics, _, _ = expand_chemistry_sets(self.rec, ['PAH_formation'])
        self.assertIn('naphthalene_H', thermo)
        self.assertIn('Mebel_C6H5_C2H2', kinetics)

    def test_surface_nitrogen(self):
        _, kinetics, _, _ = expand_chemistry_sets(self.rec, ['surface_nitrogen'])
        self.assertIn('Surface/Ammonia/Schneider_Pt111', kinetics)
        self.assertIn('Surface/DOC/Nitrogen', kinetics)

    def test_no_duplicates(self):
        thermo, kinetics, _, seeds = expand_chemistry_sets(
            self.rec, ['primary', 'nitrogen', 'oxidation', 'CH_pyrolysis_core'])
        self.assertEqual(len(thermo), len(set(thermo)))
        self.assertEqual(len(kinetics), len(set(kinetics)))

    def test_primary_always_first(self):
        thermo, _, _, _ = expand_chemistry_sets(self.rec, ['nitrogen', 'primary'])
        self.assertEqual(thermo[0], 'primaryThermoLibrary')

    def test_all_primary_before_other(self):
        primary_thermo = self.rec['primary']['thermo']
        thermo, _, _, _ = expand_chemistry_sets(self.rec, ['primary', 'oxidation', 'CH_pyrolysis_core'])
        last_primary = max(thermo.index(lib) for lib in primary_thermo)
        non_primary = [lib for lib in thermo if lib not in primary_thermo]
        if non_primary:
            first_other = min(thermo.index(lib) for lib in non_primary)
            self.assertLess(last_primary, first_other)

    def test_invalid_set(self):
        with self.assertRaises(Exception):
            expand_chemistry_sets(self.rec, ['nonexistent'])

    def test_empty(self):
        thermo, kinetics, transport, seeds = expand_chemistry_sets(self.rec, [])
        self.assertEqual(thermo, [])


class TestMergeWithUserLibraries(unittest.TestCase):

    def test_auto_string(self):
        self.assertEqual(merge_with_user_libraries('auto', ['a', 'b']), ['a', 'b'])

    def test_no_auto(self):
        self.assertEqual(merge_with_user_libraries(['x', 'y'], ['a']), ['x', 'y'])

    def test_auto_middle(self):
        self.assertEqual(merge_with_user_libraries(['x', 'auto', 'y'], ['a', 'b']), ['x', 'a', 'b', 'y'])

    def test_auto_start(self):
        self.assertEqual(merge_with_user_libraries(['auto', 'y'], ['a', 'b']), ['a', 'b', 'y'])

    def test_auto_end(self):
        self.assertEqual(merge_with_user_libraries(['x', 'auto'], ['a']), ['x', 'a'])

    def test_dedup(self):
        self.assertEqual(merge_with_user_libraries(['x', 'auto'], ['x', 'a']), ['x', 'a'])

    def test_none(self):
        self.assertIsNone(merge_with_user_libraries(None, ['a']))

    def test_empty(self):
        self.assertEqual(merge_with_user_libraries([], ['a']), [])

    def test_pah_stripped(self):
        self.assertEqual(merge_with_user_libraries(['x', PAH_LIBS, 'auto'], ['a']), ['x', 'a'])

    def test_pah_no_auto(self):
        self.assertEqual(merge_with_user_libraries(['x', PAH_LIBS], ['a']), ['x'])

    def test_auto_empty_libs(self):
        self.assertEqual(merge_with_user_libraries('auto', []), [])


class TestToReactionLibraryTuples(unittest.TestCase):

    def test_all_false(self):
        self.assertEqual(to_reaction_library_tuples(['a', 'b'], set()), [('a', False), ('b', False)])

    def test_some_true(self):
        self.assertEqual(to_reaction_library_tuples(['a', 'b', 'c'], {'b'}), [('a', False), ('b', True), ('c', False)])

    def test_empty(self):
        self.assertEqual(to_reaction_library_tuples([], set()), [])

    def test_all_true(self):
        self.assertEqual(to_reaction_library_tuples(['a', 'b'], {'a', 'b'}), [('a', True), ('b', True)])


class TestResolveAutoKineticsFamilies(unittest.TestCase):

    def setUp(self):
        self.db_dir = settings['database.directory']
        if not os.path.isfile(os.path.join(self.db_dir, 'kinetics', 'families', 'recommended.py')):
            self.skipTest('recommended.py not found')

    def test_default(self):
        families = resolve_auto_kinetics_families([FamilySet.DEFAULT], self.db_dir)
        self.assertIn('H_Abstraction', families)
        self.assertIn('R_Recombination', families)

    def test_default_plus_surface(self):
        families = resolve_auto_kinetics_families([FamilySet.DEFAULT, FamilySet.SURFACE], self.db_dir)
        self.assertIn('Surface_Adsorption_Single', families)

    def test_no_duplicates(self):
        families = resolve_auto_kinetics_families(
            [FamilySet.DEFAULT, FamilySet.SURFACE, FamilySet.HALOGENS], self.db_dir)
        self.assertEqual(len(families), len(set(families)))

    def test_invalid_set(self):
        with self.assertRaises(Exception):
            resolve_auto_kinetics_families(['nonexistent'], self.db_dir)

    def test_invalid_path(self):
        with self.assertRaises(Exception):
            resolve_auto_kinetics_families([FamilySet.DEFAULT], '/nonexistent')

    def test_returns_strings(self):
        for f in resolve_auto_kinetics_families([FamilySet.DEFAULT], self.db_dir):
            self.assertIsInstance(f, str)


class TestLogLibList(unittest.TestCase):

    def test_none(self):
        with self.assertLogs(level=logging.INFO) as cm:
            _log_lib_list('Test', None)
        self.assertIn('(none)', cm.output[0])

    def test_empty(self):
        with self.assertLogs(level=logging.INFO) as cm:
            _log_lib_list('Test', [])
        self.assertIn('(none)', cm.output[0])

    def test_short(self):
        with self.assertLogs(level=logging.INFO) as cm:
            _log_lib_list('Libs', ['a', 'b'])
        self.assertIn('Libs (2)', cm.output[0])

    def test_wrapping(self):
        libs = [f'very_long_library_name_{i}' for i in range(10)]
        with self.assertLogs(level=logging.INFO) as cm:
            _log_lib_list('Libs', libs, width=60)
        self.assertIn('\n', cm.output[0])


class TestAutoSelectLibraries(unittest.TestCase):

    def setUp(self):
        if not os.path.isfile(os.path.join(settings['database.directory'], 'recommended_libraries.yml')):
            self.skipTest('recommended_libraries.yml not found')

    def test_noop_without_auto(self):
        rmg = RMG()
        rmg.initial_species = [Species().from_smiles('C')]
        rmg.reaction_systems = [_simple_reactor(500.0)]
        rmg.solvent = None
        rmg.database_directory = settings['database.directory']
        rmg.thermo_libraries = ['myLib']
        rmg.reaction_libraries = ['myKinLib']
        rmg.seed_mechanisms = []
        rmg.transport_libraries = None
        rmg.kinetics_families = 'default'
        auto_select_libraries(rmg)
        self.assertEqual(rmg.thermo_libraries, ['myLib'])

    def test_methane_oxidation(self):
        rmg = _rmg_with_auto(['C', '[O][O]'], 1500.0)
        auto_select_libraries(rmg)
        self.assertIn('primaryThermoLibrary', rmg.thermo_libraries)
        self.assertIn('FFCM1(-)', rmg.thermo_libraries)
        self.assertIn('H_Abstraction', rmg.kinetics_families)

    def test_ethane_pyrolysis_includes_pah(self):
        rmg = _rmg_with_auto(['CC'], 1200.0)
        auto_select_libraries(rmg)
        self.assertIn('naphthalene_H', rmg.thermo_libraries)
        self.assertIn('Mebel_C6H5_C2H2', rmg.reaction_libraries)

    def test_ethanol_no_pah(self):
        rmg = _rmg_with_auto(['CCO'], 1000.0)
        auto_select_libraries(rmg)
        self.assertNotIn('naphthalene_H', rmg.thermo_libraries)
        self.assertIn('FFCM1(-)', rmg.thermo_libraries)

    def test_ethanol_with_pah(self):
        rmg = _rmg_with_auto(['CCO'], 1000.0, pah=True)
        auto_select_libraries(rmg)
        self.assertIn('naphthalene_H', rmg.thermo_libraries)

    def test_families_resolved(self):
        rmg = _rmg_with_auto(['C'], 500.0)
        auto_select_libraries(rmg)
        self.assertIsInstance(rmg.kinetics_families, list)
        self.assertIn('H_Abstraction', rmg.kinetics_families)

    def test_families_auto_with_exclusion(self):
        """['!H_Abstraction', 'auto'] should resolve auto families minus H_Abstraction."""
        rmg = _rmg_with_auto(['C'], 500.0)
        rmg.kinetics_families = ['!H_Abstraction', AUTO]
        auto_select_libraries(rmg)
        self.assertIsInstance(rmg.kinetics_families, list)
        self.assertNotIn('H_Abstraction', rmg.kinetics_families)
        self.assertIn('R_Recombination', rmg.kinetics_families)

    def test_families_auto_with_multiple_exclusions(self):
        rmg = _rmg_with_auto(['C'], 500.0)
        rmg.kinetics_families = ['!H_Abstraction', '!Disproportionation', AUTO]
        auto_select_libraries(rmg)
        self.assertNotIn('H_Abstraction', rmg.kinetics_families)
        self.assertNotIn('Disproportionation', rmg.kinetics_families)
        self.assertIn('R_Recombination', rmg.kinetics_families)

    def test_families_auto_with_addition(self):
        """['auto', 'MyCustomFamily'] should include auto families + the custom one."""
        rmg = _rmg_with_auto(['C'], 500.0)
        rmg.kinetics_families = [AUTO, 'MyCustomFamily']
        auto_select_libraries(rmg)
        self.assertIn('H_Abstraction', rmg.kinetics_families)
        self.assertIn('MyCustomFamily', rmg.kinetics_families)

    def test_families_auto_with_exclusion_and_addition(self):
        """['!H_Abstraction', 'auto', 'MyCustomFamily'] should work."""
        rmg = _rmg_with_auto(['C'], 500.0)
        rmg.kinetics_families = ['!H_Abstraction', AUTO, 'MyCustomFamily']
        auto_select_libraries(rmg)
        self.assertNotIn('H_Abstraction', rmg.kinetics_families)
        self.assertIn('R_Recombination', rmg.kinetics_families)
        self.assertIn('MyCustomFamily', rmg.kinetics_families)


class TestEndToEnd(unittest.TestCase):

    def test_methane_oxidation(self):
        species = [Species().from_smiles('C'), Species().from_smiles('[O][O]')]
        profile = detect_chemistry(species, [_simple_reactor(1500.0)], solvent=None)
        sets = determine_chemistry_sets(profile, False)
        self.assertIn(ChemistrySet.OXIDATION, sets)
        self.assertIn(ChemistrySet.CH_PYROLYSIS_CORE, sets)
        self.assertNotIn(ChemistrySet.PAH_FORMATION, sets)

    def test_ethane_pyrolysis(self):
        species = [Species().from_smiles('CC')]
        profile = detect_chemistry(species, [_simple_reactor(1200.0)], solvent=None)
        sets = determine_chemistry_sets(profile, False)
        self.assertIn(ChemistrySet.PAH_FORMATION, sets)
        self.assertNotIn(ChemistrySet.OXIDATION, sets)

    def test_liquid_pentane_oxidation(self):
        species = [Species().from_smiles('CCCCC'), Species().from_smiles('[O][O]')]
        profile = detect_chemistry(species, [_liquid_reactor(400.0)], solvent='pentane')
        sets = determine_chemistry_sets(profile, False)
        self.assertIn(ChemistrySet.LIQUID_OXIDATION, sets)
        self.assertIn(FamilySet.LIQUID_PEROXIDE, determine_kinetics_families(profile))

    def test_nitrogen_sulfur(self):
        species = [Species().from_smiles('N'), Species().from_smiles('S')]
        profile = detect_chemistry(species, [_simple_reactor(500.0)], solvent=None)
        sets = determine_chemistry_sets(profile, False)
        self.assertIn(ChemistrySet.NITROGEN, sets)
        self.assertIn(ChemistrySet.SULFUR, sets)
        self.assertNotIn(ChemistrySet.OXIDATION, sets)

    def test_inert_with_nitrogen(self):
        species = [Species().from_smiles('[Ar]'), Species().from_smiles('N#N')]
        profile = detect_chemistry(species, [_simple_reactor(300.0)], solvent=None)
        sets = determine_chemistry_sets(profile, False)
        self.assertEqual(sets, [ChemistrySet.PRIMARY, ChemistrySet.NITROGEN])

    def test_oxyfuel_with_inert_bath_gas(self):
        """
        Formic acid (OC=O) with non-reactive N2 bath gas at 1000 K.
        N2 is inert so nitrogen should NOT be triggered.
        Oxygen is in the fuel so oxidation + CH_pyrolysis_core should trigger.
        PAH_formation should NOT trigger (oxygen present, no <PAH_libs>).
        """
        formic_acid = Species().from_smiles('OC=O')
        n2_bath = Species(reactive=False).from_smiles('N#N')
        profile = detect_chemistry([formic_acid, n2_bath], [_simple_reactor(1000.0)], solvent=None)
        self.assertTrue(profile.has_oxygen)
        self.assertTrue(profile.has_carbon)
        self.assertFalse(profile.has_nitrogen)  # N2 is non-reactive, should be ignored
        sets = determine_chemistry_sets(profile, pah_libs_requested=False)
        self.assertIn(ChemistrySet.PRIMARY, sets)
        self.assertIn(ChemistrySet.OXIDATION, sets)
        self.assertIn(ChemistrySet.CH_PYROLYSIS_CORE, sets)
        self.assertNotIn(ChemistrySet.PAH_FORMATION, sets)
        self.assertNotIn(ChemistrySet.NITROGEN, sets)


if __name__ == '__main__':
    unittest.main()
