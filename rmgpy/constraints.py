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
    the speciesConstraints set by the user.  If not, returns `True` for failing speciesConstraints.
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

    explicitly_allowed_molecules = species_constraints.get('explicitlyAllowedMolecules', [])
    for molecule in explicitly_allowed_molecules:
        if struct.is_isomorphic(molecule):
            return False

    max_carbon_atoms = species_constraints.get('maximumCarbonAtoms', -1)
    if max_carbon_atoms != -1:
        if struct.get_num_atoms('C') > max_carbon_atoms:
            return True

    max_oxygen_atoms = species_constraints.get('maximumOxygenAtoms', -1)
    if max_oxygen_atoms != -1:
        if struct.get_num_atoms('O') > max_oxygen_atoms:
            return True

    max_nitrogen_atoms = species_constraints.get('maximumNitrogenAtoms', -1)
    if max_nitrogen_atoms != -1:
        if struct.get_num_atoms('N') > max_nitrogen_atoms:
            return True

    max_silicon_atoms = species_constraints.get('maximumSiliconAtoms', -1)
    if max_silicon_atoms != -1:
        if struct.get_num_atoms('Si') > max_silicon_atoms:
            return True

    max_sulfur_atoms = species_constraints.get('maximumSulfurAtoms', -1)
    if max_sulfur_atoms != -1:
        if struct.get_num_atoms('S') > max_sulfur_atoms:
            return True

    max_heavy_atoms = species_constraints.get('maximumHeavyAtoms', -1)
    if max_heavy_atoms != -1:
        if struct.get_num_atoms() - struct.get_num_atoms('H') > max_heavy_atoms:
            return True

    max_surface_sites = species_constraints.get('maximumSurfaceSites', -1)
    if max_surface_sites != -1:
        if struct.get_num_atoms('X') > max_surface_sites:
            return True

    max_surface_bond_order = species_constraints.get('maximumSurfaceBondOrder', -1)
    if max_surface_bond_order != -1:
        for site in struct.get_surface_sites():
            if site.get_total_bond_order() > max_surface_bond_order:
                return True

    max_radicals = species_constraints.get('maximumRadicalElectrons', -1)
    if max_radicals != -1:
        if struct.get_radical_count() > max_radicals:
            return True

    max_carbenes = species_constraints.get('maximumSingletCarbenes', 1)
    if max_radicals != -1:
        if struct.get_singlet_carbene_count() > max_carbenes:
            return True

    max_carbene_radicals = species_constraints.get('maximumCarbeneRadicals', 1) #changing this from 0 to 1 (PFAS chemistry involves carbene radicals)
    if max_carbene_radicals != -1:
        if struct.get_singlet_carbene_count() > 0 and struct.get_radical_count() > max_carbene_radicals:
            return True

    return False
