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

import logging

from rmgpy.species import Species


# ---------------------------------------------------------------------------
# Polymer-constraint warn-once state
# ---------------------------------------------------------------------------

_unbounded_polymer_warned = False


def reset_polymer_warning():
    """Reset the once-per-run unbounded-polymer warning flag. Called from RMG.initialize."""
    global _unbounded_polymer_warned
    _unbounded_polymer_warned = False


def _warn_unbounded_polymer_once():
    global _unbounded_polymer_warned
    if not _unbounded_polymer_warned:
        logging.warning(
            "Polymer/proxy species bypass generatedSpeciesConstraints because "
            "generatePolymerConstraints is not set; polymer proxy chemistry is unbounded. "
            "Add a generatePolymerConstraints block to bound polymer-family reaction generation."
        )
        _unbounded_polymer_warned = True


def is_polymer_constraint_member(obj):
    """
    True if `obj` (a Species, Polymer, or Molecule) should be bounded by polymer
    constraints rather than the global generatedSpeciesConstraints. Uses the strongest
    signal available on the object: `is_polymer` (durable, Species/Polymer only) or
    `is_polymer_proxy` (the only signal present on a raw Molecule). Duck-typed to avoid a
    Polymer import / circular dependency.
    """
    return bool(getattr(obj, 'is_polymer', False) or getattr(obj, 'is_polymer_proxy', False))


def _normalize(species):
    """Return the Molecule to constrain: molecule[0] for a Species, else the object itself."""
    if isinstance(species, Species):
        return species.molecule[0]
    return species


def _evaluate_constraints(struct, constraints):
    """Return a reason string if `struct` violates `constraints`, else False. Pure: no get_input."""
    explicitly_allowed_molecules = constraints.get('explicitlyAllowedMolecules', [])
    for molecule in explicitly_allowed_molecules:
        if struct.is_isomorphic(molecule):
            return False

    max_carbon_atoms = constraints.get('maximumCarbonAtoms', -1)
    if max_carbon_atoms != -1:
        if struct.get_num_atoms('C') > max_carbon_atoms:
            return f"Exceeded maximumCarbonAtoms: {struct.get_num_atoms('C')} > {max_carbon_atoms}"

    max_oxygen_atoms = constraints.get('maximumOxygenAtoms', -1)
    if max_oxygen_atoms != -1:
        if struct.get_num_atoms('O') > max_oxygen_atoms:
            return f"Exceeded maximumOxygenAtoms: {struct.get_num_atoms('O')} > {max_oxygen_atoms}"

    max_nitrogen_atoms = constraints.get('maximumNitrogenAtoms', -1)
    if max_nitrogen_atoms != -1:
        if struct.get_num_atoms('N') > max_nitrogen_atoms:
            return f"Exceeded maximumNitrogenAtoms: {struct.get_num_atoms('N')} > {max_nitrogen_atoms}"

    max_silicon_atoms = constraints.get('maximumSiliconAtoms', -1)
    if max_silicon_atoms != -1:
        if struct.get_num_atoms('Si') > max_silicon_atoms:
            return f"Exceeded maximumSiliconAtoms: {struct.get_num_atoms('Si')} > {max_silicon_atoms}"

    max_sulfur_atoms = constraints.get('maximumSulfurAtoms', -1)
    if max_sulfur_atoms != -1:
        if struct.get_num_atoms('S') > max_sulfur_atoms:
            return f"Exceeded maximumSulfurAtoms: {struct.get_num_atoms('S')} > {max_sulfur_atoms}"

    max_heavy_atoms = constraints.get('maximumHeavyAtoms', -1)
    if max_heavy_atoms != -1:
        heavy_atoms = struct.get_num_atoms() - struct.get_num_atoms('H')
        if heavy_atoms > max_heavy_atoms:
            return f"Exceeded maximumHeavyAtoms: {heavy_atoms} > {max_heavy_atoms}"

    max_surface_sites = constraints.get('maximumSurfaceSites', -1)
    if max_surface_sites != -1:
        if struct.get_num_atoms('X') > max_surface_sites:
            return f"Exceeded maximumSurfaceSites: {struct.get_num_atoms('X')} > {max_surface_sites}"

    max_surface_bond_order = constraints.get('maximumSurfaceBondOrder', -1)
    if max_surface_bond_order != -1:
        for site in struct.get_surface_sites():
            if site.get_total_bond_order() > max_surface_bond_order:
                return f"Exceeded maximumSurfaceBondOrder at site: {site.get_total_bond_order()} > {max_surface_bond_order}"

    max_radicals = constraints.get('maximumRadicalElectrons', -1)
    if max_radicals != -1:
        if struct.get_radical_count() > max_radicals:
            return f"Exceeded maximumRadicalElectrons: {struct.get_radical_count()} > {max_radicals}"

    max_carbenes = constraints.get('maximumSingletCarbenes', 1)
    if max_carbenes != -1:
        if struct.get_singlet_carbene_count() > max_carbenes:
            return f"Exceeded maximumSingletCarbenes: {struct.get_singlet_carbene_count()} > {max_carbenes}"

    max_carbene_radicals = constraints.get('maximumCarbeneRadicals', 0)
    if max_carbene_radicals != -1:
        if struct.get_singlet_carbene_count() > 0 and struct.get_radical_count() > max_carbene_radicals:
            return f"Exceeded maximumCarbeneRadicals: {struct.get_radical_count()} > {max_carbene_radicals}"

    return False


def pass_cutting_threshold(species):
    """
    Pass in either a `Species` or `Molecule` object and checks whether it passes
    the speciesCuttingThreshold set by the user. The default value is 20. If yes,
    returns `True` for passing cutting threshold.
    """

    from rmgpy.rmg.input import get_input

    try:
        species_constraints = get_input('species_constraints')
    except Exception:
        logging.debug('Species constraints could not be found.')
        species_constraints = {}

    if isinstance(species, Species):
        struct = species.molecule[0]
    else:
        # expects a molecule here
        struct = species

    min_cutting_size = species_constraints.get('speciesCuttingThreshold', 20)
    if struct.get_element_count()['C'] >= min_cutting_size:
        return True

    return False

def fails_species_constraints(species):
    """
    Pass in either a `Species` or `Molecule` object and checks whether it passes
    the speciesConstraints set by the user.  If the species fails constraints, returns
    a string `reason` describing which constraint failed. If all constraints pass, returns `False`.
    """
    from rmgpy.rmg.input import get_input

    if is_polymer_constraint_member(species):
        try:
            polymer_constraints = get_input('polymer_constraints')
        except Exception:
            logging.debug('Polymer constraints could not be found.')
            polymer_constraints = None
        if polymer_constraints is None:
            _warn_unbounded_polymer_once()
            return False
        return _evaluate_constraints(_normalize(species), polymer_constraints)

    try:
        species_constraints = get_input('species_constraints')
    except Exception:
        logging.debug('Species constraints could not be found.')
        species_constraints = {}

    return _evaluate_constraints(_normalize(species), species_constraints)
