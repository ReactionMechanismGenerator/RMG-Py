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
Automated kinetics library, kinetics family, thermo library and transport library selection for RMG
based on the chemistry detected in the input species and reactor conditions.
"""

import logging
import os
from dataclasses import dataclass, field
from typing import Any, List, Optional, Set

import yaml

from rmgpy.exceptions import InputError
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.solver.surface import SurfaceReactor

try:
    from rmgpy.rmg.reactionmechanismsimulator_reactors import (
        ConstantTLiquidSurfaceReactor as RMSLiqSurf,
        ConstantTVLiquidReactor as RMSLiq,
    )
except ImportError:
    RMSLiqSurf = None
    RMSLiq = None

# Values used in input files to request auto-selection
AUTO = 'auto'
PAH_LIBS = '<PAH_libs>'

# Temperature threshold (in K) above which CH pyrolysis libraries are included
CH_PYROLYSIS_T_THRESHOLD = 800.0

# Elements that trigger the halogens chemistry set
HALOGEN_ELEMENTS = {'F', 'Cl', 'Br', 'I'}
# Elements that trigger the metal/electrochem chemistry set (currently only Li)
METAL_ELEMENTS = {'Li'}


@dataclass
class ChemistryProfile:
    """Detected chemistry characteristics from initial species and reactor conditions."""
    elements_present: Set[str] = field(default_factory=set)
    has_nitrogen: bool = False
    has_sulfur: bool = False
    has_oxygen: bool = False
    has_carbon: bool = False
    has_halogens: bool = False
    has_metal: bool = False
    has_surface: bool = False
    has_liquid: bool = False
    max_temperature: float = 0.0


def detect_chemistry(initial_species: list,
                     reaction_systems: list,
                     solvent: Optional[str],
                     ) -> ChemistryProfile:
    """
    Analyze initial species and reactor conditions to build a ChemistryProfile.

    Args:
        initial_species: list of Species objects from the RMG input.
        reaction_systems: list of reactor system objects.
        solvent: solvent name string, or None for gas phase.

    Returns:
        ChemistryProfile with detected characteristics.
    """
    profile = ChemistryProfile()

    for spc in initial_species:
        mol = spc.molecule[0]
        element_counts = mol.get_element_count()
        profile.elements_present.update(element_counts.keys())
        if mol.is_surface_site() or 'X' in element_counts:
            profile.has_surface = True

    profile.has_nitrogen = 'N' in profile.elements_present
    profile.has_sulfur = 'S' in profile.elements_present
    profile.has_oxygen = 'O' in profile.elements_present
    profile.has_carbon = 'C' in profile.elements_present
    profile.has_halogens = bool(profile.elements_present & HALOGEN_ELEMENTS)  # bool(set intersection)
    profile.has_metal = bool(profile.elements_present & METAL_ELEMENTS)

    for reactor in reaction_systems:
        if isinstance(reactor, LiquidReactor) or RMSLiq is not None and isinstance(reactor, (RMSLiq, RMSLiqSurf)):
            profile.has_liquid = True
        elif isinstance(reactor, SurfaceReactor) or RMSLiq is not None and isinstance(reactor, RMSLiqSurf):
            profile.has_surface = True

        T = _get_reactor_max_temperature(reactor)
        if T is not None and T > profile.max_temperature:
            profile.max_temperature = T

    if solvent is not None:
        profile.has_liquid = True

    return profile


def _get_reactor_max_temperature(reactor: Any) -> Optional[float]:
    """
    Extract the maximum temperature (in K) from a reactor system.
    Handles both single-value and range temperature specifications.

    Returns:
        float temperature in K, or None if it cannot be determined.
    """
    T = getattr(reactor, 'T', None)
    if T is None:
        # RMS reactors store T in initial_conditions
        ic = getattr(reactor, 'initial_conditions', None)
        if ic is not None:
            if isinstance(ic, dict):
                if 'T' in ic:
                    return float(ic['T'])
                # Nested dict for liquid-surface reactors
                for phase_key in ('liquid', 'gas', 'Default'):
                    if isinstance(ic.get(phase_key), dict) and 'T' in ic[phase_key]:
                        return float(ic[phase_key]['T'])
        return None

    if isinstance(T, list):
        # Temperature range: return the max
        return max(t.value_si for t in T)
    else:
        return T.value_si


def _has_pah_libs_keyword(rmg: Any) -> bool:
    """Check whether the user included the <PAH_libs> keyword in any library field."""
    for attr in ('thermo_libraries', 'reaction_libraries', 'seed_mechanisms',
                 'transport_libraries'):
        val = getattr(rmg, attr, None)
        if isinstance(val, list):
            for item in val:
                name = item[0] if isinstance(item, tuple) else item
                if name == PAH_LIBS:
                    return True
    return False


def determine_chemistry_sets(profile: ChemistryProfile,
                             pah_libs_requested: bool = False,
                             ) -> List[str]:
    """
    Determine which chemistry sets to activate based on the detected profile.

    CH pyrolysis logic:
      - CH_pyrolysis_core is always added when C present AND T >= 800 K.
      - PAH_formation is added when:
          (a) C + T >= 800 K + no O in species (pure C/H pyrolysis), OR
          (b) C + T >= 800 K + <PAH_libs> keyword requested by user.

    Args:
        profile: ChemistryProfile instance.
        pah_libs_requested: bool, True if user included <PAH_libs> keyword.

    Returns:
        List of chemistry set names in priority order.
    """
    sets = ['primary']  # always included

    if profile.has_nitrogen:
        sets.append('nitrogen')

    if profile.has_sulfur:
        sets.append('sulfur')

    if profile.has_oxygen:
        sets.append('combustion')

    high_T_carbon = profile.has_carbon and profile.max_temperature >= CH_PYROLYSIS_T_THRESHOLD

    if high_T_carbon:
        sets.append('CH_pyrolysis_core')

        if not profile.has_oxygen or pah_libs_requested:
            sets.append('PAH_formation')

    if profile.has_liquid and profile.has_oxygen:
        sets.append('liquid_oxidation')

    if profile.has_surface:
        sets.append('surface')

    if profile.has_halogens:
        sets.append('halogens')

    if profile.has_metal:
        sets.append('lithium')

    return sets


def determine_kinetics_families(profile: ChemistryProfile) -> List[str]:
    """
    Determine which kinetics family sets to activate based on the detected profile.

    These correspond to named sets in RMG-database/input/kinetics/families/recommended.py.

    Args:
        profile: ChemistryProfile instance.

    Returns:
        List of family set names to combine.
    """
    family_sets = ['default']  # always included

    if profile.has_liquid and profile.has_oxygen:
        family_sets.append('liquid_peroxide')

    if profile.has_surface:
        family_sets.append('surface')

    if profile.has_halogens:
        family_sets.append('halogens')

    if profile.has_metal:
        family_sets.append('electrochem')

    return family_sets


def load_recommended_yml(database_directory: str) -> dict:
    """
    Load the recommended_libraries.yml file from the RMG database.

    Args:
        database_directory: path to the RMG database 'input' directory.

    Returns:
        dict parsed from YAML.
    """
    yml_path = os.path.join(database_directory, 'recommended_libraries.yml')
    if not os.path.isfile(yml_path):
        raise InputError(
            f"Could not find recommended_libraries.yml at {yml_path}. "
            f"This file is required for 'auto' library selection."
        )
    with open(yml_path, 'r') as f:
        return yaml.safe_load(f)


def expand_chemistry_sets(recommended_data: dict,
                          set_names: List[str],
                          ) -> tuple:
    """
    Expand named chemistry sets into concrete library lists.

    Args:
        recommended_data: dict from recommended.yml.
        set_names: list of chemistry set names to expand.

    Returns:
        Tuple of (thermo_libraries, kinetics_libraries, transport_libraries, seed_mechanisms)
        where each is a list of library name strings.
    """
    thermo = []
    kinetics = []
    transport = []
    seeds = []

    for set_name in set_names:
        if set_name not in recommended_data:
            raise InputError(
                f"Chemistry set '{set_name}' not found in recommended.yml. "
                f"Available sets: {list(recommended_data.keys())}"
            )
        set_data = recommended_data[set_name]

        for entry in set_data.get('thermo', []):
            name = entry if isinstance(entry, str) else entry['name']
            if name not in thermo:
                thermo.append(name)

        for entry in set_data.get('kinetics', []):
            if isinstance(entry, str):
                if entry not in kinetics:
                    kinetics.append(entry)
            elif isinstance(entry, dict):
                name = entry['name']
                if entry.get('seed', False):
                    if name not in seeds:
                        seeds.append(name)
                else:
                    if name not in kinetics:
                        kinetics.append(name)

        for entry in set_data.get('transport', []):
            name = entry if isinstance(entry, str) else entry['name']
            if name not in transport:
                transport.append(name)

    return thermo, kinetics, transport, seeds


def merge_with_user_libraries(user_spec: Any, auto_libs: List[str]) -> Any:
    """
    Merge user-specified libraries with auto-selected libraries,
    respecting the position of the 'auto' token. <PAH_libs> tokens
    are silently removed (they've already been used as a signal).

    Args:
        user_spec: the user's library specification. Can be:
            - 'auto' (string): fully replace with auto_libs
            - list containing 'auto' token: replace token in-place with auto_libs
            - list without 'auto': return as-is (with <PAH_libs> stripped)
            - None or []: return as-is
        auto_libs: list of auto-selected library names.

    Returns:
        Resolved list of library names.
    """
    if user_spec == AUTO:
        return list(auto_libs)

    if not isinstance(user_spec, list):
        return user_spec

    # Collect all user-specified library names (excluding special tokens)
    user_lib_names = set()
    for item in user_spec:
        if item not in (AUTO, PAH_LIBS):
            name = item[0] if isinstance(item, tuple) else item
            user_lib_names.add(name)

    # Filter auto libs to exclude any already specified by user
    filtered_auto = [lib for lib in auto_libs if lib not in user_lib_names]

    # Replace tokens in-place
    result = []
    for item in user_spec:
        if item == AUTO:
            result.extend(filtered_auto)
        elif item == PAH_LIBS:
            continue  # silently consumed
        else:
            result.append(item)

    return result


def merge_with_user_reaction_libraries(user_spec: Any, auto_kinetics_libs: List[str]) -> Any:
    """
    Merge user-specified reaction libraries with auto-selected ones.

    Reaction libraries are stored as tuples (name, bool) in RMG.
    The user_spec is already converted to tuples by input.py.

    Args:
        user_spec: list of (name, bool) tuples, possibly containing 'auto' string.
        auto_kinetics_libs: list of library name strings from auto-selection.

    Returns:
        Resolved list of (name, bool) tuples.
    """
    if user_spec == AUTO:
        return [(name, False) for name in auto_kinetics_libs]

    if not isinstance(user_spec, list):
        return user_spec

    # Collect user-specified library names (excluding special tokens)
    user_lib_names = set()
    for item in user_spec:
        if item in (AUTO, PAH_LIBS):
            continue
        name = item[0] if isinstance(item, tuple) else item
        user_lib_names.add(name)

    filtered_auto = [(name, False) for name in auto_kinetics_libs if name not in user_lib_names]

    result = []
    for item in user_spec:
        if item == AUTO:
            result.extend(filtered_auto)
        elif item == PAH_LIBS:
            continue  # silently consumed
        else:
            result.append(item)

    return result


def resolve_auto_kinetics_families(family_set_names: List[str],
                                   database_directory: str,
                                   ) -> List[str]:
    """
    Resolve 'auto' kinetics families by combining the named sets from recommended.py.

    Reads the recommended.py file from the database and combines all requested sets
    into a single list of family names.

    Args:
        family_set_names: list of family set names (e.g., ['default', 'surface']).
        database_directory: path to the RMG database root.

    Returns:
        List of family name strings.
    """
    recommended_path = os.path.join(database_directory, 'kinetics', 'families', 'recommended.py')
    if not os.path.isfile(recommended_path):
        raise InputError(
            f"Could not find recommended.py at {recommended_path}. "
            f"This file is required for 'auto' kinetics families selection."
        )

    # Execute the recommended.py file to get the family sets
    local_context = {}
    global_context = {}
    with open(recommended_path, 'r') as f:
        exec(f.read(), global_context, local_context)

    combined = []
    for set_name in family_set_names:
        if set_name not in local_context:
            raise InputError(
                f"Kinetics family set '{set_name}' not found in recommended.py. "
                f"Available sets: {[k for k in local_context if not k.startswith('_')]}"
            )
        family_set = local_context[set_name]
        for family in family_set:
            if family not in combined:
                combined.append(family)

    return combined


def auto_select_libraries(rmg: Any) -> None:
    """
    Main entry point for auto library selection.

    Inspects the RMG object's initial_species, reaction_systems, and solvent
    to detect the chemistry, then resolves any 'auto' tokens in the library
    specifications.

    Modifies the rmg object in-place.

    Args:
        rmg: the RMG job object with populated initial_species, reaction_systems, etc.
    """
    # Check if any field uses 'auto' or '<PAH_libs>'
    has_special = False
    for attr in ('thermo_libraries', 'reaction_libraries', 'seed_mechanisms',
                 'transport_libraries', 'kinetics_families'):
        val = getattr(rmg, attr, None)
        if val in (AUTO, PAH_LIBS):
            has_special = True
            break
        if isinstance(val, list):
            for item in val:
                name = item[0] if isinstance(item, tuple) else item
                if name in (AUTO, PAH_LIBS):
                    has_special = True
                    break
            if has_special:
                break

    if not has_special:
        return

    # Check for <PAH_libs> keyword before we strip it
    pah_libs_requested = _has_pah_libs_keyword(rmg)

    # Detect chemistry
    profile = detect_chemistry(rmg.initial_species, rmg.reaction_systems, rmg.solvent)

    logging.info('')
    logging.info('=' * 80)
    logging.info('Automatic library selection')
    logging.info('=' * 80)
    logging.info(f'Elements detected: {sorted(profile.elements_present)}')
    logging.info(f'Max temperature: {profile.max_temperature:.1f} K')
    logging.info(f'Has nitrogen: {profile.has_nitrogen}')
    logging.info(f'Has sulfur: {profile.has_sulfur}')
    logging.info(f'Has oxygen: {profile.has_oxygen}')
    logging.info(f'Has carbon: {profile.has_carbon}')
    logging.info(f'Has halogens: {profile.has_halogens}')
    logging.info(f'Has lithium: {profile.has_metal}')
    logging.info(f'Has surface: {profile.has_surface}')
    logging.info(f'Has liquid: {profile.has_liquid}')
    logging.info(f'PAH_libs requested: {pah_libs_requested}')

    # Determine chemistry sets
    set_names = determine_chemistry_sets(profile, pah_libs_requested)
    logging.info(f'Triggered chemistry sets: {set_names}')

    # Load and expand recommended_libraries.yml
    recommended_data = load_recommended_yml(rmg.database_directory)
    auto_thermo, auto_kinetics, auto_transport, auto_seeds = expand_chemistry_sets(
        recommended_data, set_names
    )

    # Resolve thermo libraries
    rmg.thermo_libraries = merge_with_user_libraries(rmg.thermo_libraries, auto_thermo)
    logging.info(f'Thermo libraries: {rmg.thermo_libraries}')

    # Resolve reaction libraries
    rmg.reaction_libraries = merge_with_user_reaction_libraries(
        rmg.reaction_libraries, auto_kinetics
    )
    logging.info(f'Reaction libraries: {[name for name, _ in rmg.reaction_libraries]}')

    # Resolve seed mechanisms
    rmg.seed_mechanisms = merge_with_user_libraries(rmg.seed_mechanisms, auto_seeds)
    logging.info(f'Seed mechanisms: {rmg.seed_mechanisms}')

    # Resolve transport libraries
    rmg.transport_libraries = merge_with_user_libraries(rmg.transport_libraries, auto_transport)
    logging.info(f'Transport libraries: {rmg.transport_libraries}')

    # Resolve kinetics families
    if rmg.kinetics_families == AUTO:
        family_set_names = determine_kinetics_families(profile)
        logging.info(f'Triggered kinetics family sets: {family_set_names}')
        rmg.kinetics_families = resolve_auto_kinetics_families(
            family_set_names, rmg.database_directory
        )
        logging.info(f'Kinetics families ({len(rmg.kinetics_families)} total): {rmg.kinetics_families}')

    logging.info('='*80)
    logging.info('')
