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
This module contains functions for writing of Cantera input files.
"""

from typing import Union, TYPE_CHECKING

import os
import shutil
import logging
import yaml

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import (
    Arrhenius, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius,
    Chebyshev, Troe, Lindemann, ThirdBody,
)
from rmgpy.reaction import Reaction
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.util import make_output_subdirectory
import rmgpy.constants as constants

if TYPE_CHECKING:
    from rmgpy.species import Species
    from rmgpy.molecule.molecule import Molecule


SYMBOL_BY_NUMBER = {0: 'e', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na',
                    12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc',
                    22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga',
                    32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb',
                    42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb',
                    52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm',
                    62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
                    72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl',
                    82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa',
                    92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md',
                    102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
                    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}
NUMBER_BY_SYMBOL = {value: key for key, value in SYMBOL_BY_NUMBER.items()}


class CanteraWriter(object):
    """
    This class listens to a RMG subject and writes a Cantera YAML file
    with the current state of the RMG model at every iteration.
    """

    def __init__(self, output_directory=''):
        self.output_directory = output_directory
        make_output_subdirectory(output_directory, 'cantera')

    def update(self, rmg):
        """
        Called whenever the RMG subject notifies listeners.
        """
        save_cantera_files(rmg)


def save_cantera_files(rmg):
    """
    Save the current reaction model to a set of Cantera YAML files.

    Creates:
      1. chem{N}.yaml (where N is num species)
      2. chem.yaml (latest copy)
    """
    # Ensure subdirectory exists
    cantera_dir = os.path.join(rmg.output_directory, 'cantera')
    if not os.path.exists(cantera_dir):
        os.mkdir(cantera_dir)
    # -------------------------------------------------------------------------
    # 1. Save Core Model
    # -------------------------------------------------------------------------
    num_species = len(rmg.reaction_model.core.species)

    # Define paths
    this_cantera_path = os.path.join(rmg.output_directory, 'cantera',
                                     'chem{0:04d}.yaml'.format(num_species))
    latest_cantera_path = os.path.join(rmg.output_directory, 'cantera', 'chem.yaml')

    logging.info(f"Saving current model core to Cantera file: {this_cantera_path}")

    # Write the YAML file
    save_cantera_model(rmg.reaction_model.core, this_cantera_path)

    # Copy to 'chem.yaml' (The latest file)
    if os.path.exists(latest_cantera_path):
        os.unlink(latest_cantera_path)
    shutil.copy2(this_cantera_path, latest_cantera_path)

    # -------------------------------------------------------------------------
    # 2. Save Edge Model (Optional, matching ChemkinWriter logic)
    # -------------------------------------------------------------------------
    if rmg.save_edge_species:
        logging.info('Saving current model core and edge to Cantera file...')

        this_edge_path = os.path.join(rmg.output_directory, 'cantera',
                                      'chem_edge{0:04d}.yaml'.format(num_species))
        latest_edge_path = os.path.join(rmg.output_directory, 'cantera', 'chem_edge.yaml')

        # Combine core and edge
        # Note: We create a temporary object or just pass list concatenations
        # Creating a simple container object to pass to save_cantera_model
        class MixedModel:
            def __init__(self, species, reactions):
                self.species = species
                self.reactions = reactions

        edge_model = MixedModel(
            rmg.reaction_model.core.species + rmg.reaction_model.edge.species,
            rmg.reaction_model.core.reactions + rmg.reaction_model.edge.reactions
        )

        save_cantera_model(edge_model, this_edge_path)

        if os.path.exists(latest_edge_path):
            os.unlink(latest_edge_path)
        shutil.copy2(this_edge_path, latest_edge_path)


def save_cantera_model(model_container, path):
    """
    Internal helper to generate the dictionary and write the YAML file.
    model_container must have .species and .reactions attributes (lists).
    """
    species_list = model_container.species
    reaction_list = model_container.reactions

    is_plasma = False
    for sp in species_list:
        if sp.is_electron():
            is_plasma = True
            break

    # Generate Data
    yaml_data = generate_cantera_data(species_list, reaction_list, is_plasma=is_plasma)

    # Write
    with open(path, 'w') as f:
        # sort_keys=False ensures 'units' comes first, then 'phases', etc.
        yaml.dump(yaml_data, f, sort_keys=False, default_flow_style=None)


def generate_cantera_data(species_list, reaction_list, is_plasma=False, search_for_additional_elements=False):
    """
    Converts RMG objects into a dictionary structure compatible with Cantera YAML.
    """
    # --- 1. Header & Units ---
    # We output everything in SI units.
    data = {
        'description': 'RMG-Py Generated Mechanism',
        'generator': 'RMG-Py CanteraWriter',
        'cantera-version': '3.1',
        'units': {
            'length': 'm',
            'time': 's',
            'quantity': 'mol',
            'activation-energy': 'J/mol'
        }
    }

    # --- 2. Phase Definition ---
    base_elements = ['H', 'C', 'O', 'N', 'Ne', 'Ar', 'He', 'Si', 'S', 'F', 'Cl', 'Br', 'I', 'E', 'Li', 'Na', 'K', 'Mg', 'Ca']
    elements_set = set(base_elements)

    if search_for_additional_elements:
        for sp in species_list:
            if sp.molecule and len(sp.molecule) > 0:
                if sp.is_electron:
                    elements_set.add('E')
                    is_plasma = True
                else:
                    for elem in sp.molecule[0].get_element_count().keys():
                        elements_set.add(elem)

    phase_def = {
        'name': 'gas',
        'thermo': 'plasma' if is_plasma else 'ideal-gas',
        'elements': sorted(list(elements_set)),
        'species': [get_label(sp, species_list) for sp in species_list],
        'kinetics': 'gas',
        'reactions': 'all'
    }

    if is_plasma:
        # Plasma specific phase settings
        phase_def['transport'] = 'ionized-gas'
        phase_def['electron-energy-distribution'] = {
            'type': 'isotropic',
            'shape-factor': 2.0,  # Maxwellian default
            'mean-electron-energy': 1.0  # Placeholder eV
        }
    else:
        phase_def['transport'] = 'mixture-averaged'

    data['phases'] = [phase_def]

    # --- 3. Species Definitions ---
    species_data = []
    for sp in species_list:
        species_data.append(species_to_dict(sp, species_list))
    data['species'] = species_data

    # --- 4. Reaction Definitions ---
    # Note: Flatten list to handle MultiKinetics (duplicates) which return lists
    reaction_data = []
    for rxn in reaction_list:
        entries = reaction_to_dict_list(rxn, species_list)  # Returns a LIST of dicts
        if entries:
            reaction_data.extend(entries)
    data['reactions'] = reaction_data

    return data


def species_to_dict(species, species_list):
    """Convert an RMG Species object to a Cantera YAML dictionary."""

    notes = list()
    try:
        notes.append(species.to_smiles())
    except:
        pass

    # Composition
    mol = species.molecule[0]
    atom_dict = dict(mol.get_element_count())

    # Calculate 'E' based on net charge: E = Z - charge
    Z_mol = sum(NUMBER_BY_SYMBOL[atom] * count for atom, count in mol.get_element_count().items())
    charge = mol.get_net_charge()
    atom_dict['E'] = Z_mol - charge

    # Sort composition by atomic number
    atom_dict = {k: atom_dict[k] for k in sorted(atom_dict.keys(), key=lambda x: NUMBER_BY_SYMBOL.get(x, 999))}

    # Thermo (NASA7)
    thermo_data = species.get_thermo_data()

    # Sort polynomials by Tmin
    sorted_polys = sorted(thermo_data.polynomials, key=lambda p: p.Tmin.value_si)

    polys = []
    for poly in sorted_polys:
        polys.append({
            'T-range': [poly.Tmin.value_si, poly.Tmax.value_si],
            'data': poly.coeffs.tolist()  # a0..a6
        })

    # Build the base dictionary
    species_entry = {
        'name': get_label(species, species_list),
        'composition': atom_dict,
        'thermo': {
            'model': 'NASA7',
            'temperature-ranges': [sorted_polys[0].Tmin.value_si, sorted_polys[0].Tmax.value_si,
                                   sorted_polys[1].Tmax.value_si],
            'data': [polys[0]['data'], polys[1]['data']]
        },
    }

    # Transport (if available)
    if species.transport_data:
        td = species.transport_data

        # Robustly handle optional parameters
        dipole = 0.0
        if td.dipoleMoment is not None:
            dipole = td.dipoleMoment.value_si * 1e21 / constants.c  # Debye

        polarizability = 0.0
        if hasattr(td, 'polarizability') and td.polarizability is not None:
            polarizability = td.polarizability.value_si * 1e30  # Angstrom^3

        rot_relax = 0.0
        if hasattr(td, 'rotrelaxcollnum') and td.rotrelaxcollnum is not None:
            rot_relax = td.rotrelaxcollnum

        species_entry['transport'] = {
            'model': 'gas',
            'geometry': 'atom' if td.shapeIndex == 0 else 'linear' if td.shapeIndex == 1 else 'nonlinear',
            'well-depth': td.epsilon.value_si / constants.R,
            'diameter': td.sigma.value_si,
            'dipole': dipole,
            'rotational-relaxation': rot_relax
        }

    if species.thermo and species.thermo.comment:
        # Clean up newlines for cleaner YAML appearance
        clean_comment = species.thermo.comment.replace('\n', '; ').strip()
        notes.append(f"Thermo Source: {clean_comment}")

    if species.transport_data and species.transport_data.comment:
        notes.append(f"Transport Source: {species.transport_data.comment.strip()}")

    if notes:
        species_entry['note'] = " | ".join(notes)

    return species_entry


def reaction_to_dict_list(reaction, species_list=None):
    """
    Convert an RMG Reaction object to a LIST of Cantera YAML dictionaries.
    Returns a list because MultiKinetics (duplicates) map to multiple YAML entries.
    """
    # Check for MultiKinetics (duplicates grouped in one RMG object)
    if isinstance(reaction.kinetics, (MultiArrhenius, MultiPDepArrhenius)):
        entries = []
        # kin.arrhenius is a list of sub-kinetics
        sub_kinetics_list = reaction.kinetics.arrhenius

        for sub_kin in sub_kinetics_list:
            # Create a temporary reaction wrapper for the sub-kinetic
            sub_rxn = Reaction(
                reactants=reaction.reactants,
                products=reaction.products,
                reversible=reaction.reversible,
                kinetics=sub_kin,
                duplicate=reaction.duplicate  # Propagate duplicate flag
            )
            # Recursively call (should return a list of 1)
            sub_result = reaction_to_dict_list(sub_rxn, species_list)
            if sub_result:
                entries.extend(sub_result)
        return entries

    # --- Single Kinetics Logic ---

    kin = reaction.kinetics

    # 1. Determine Equation String Components
    reactants_str = " + ".join([get_label(r, species_list) for r in reaction.reactants])
    products_str = " + ".join([get_label(p, species_list) for p in reaction.products])

    # Handle Third Body suffixes (Required by Cantera for these types)
    suffix = ""
    if isinstance(kin, (ThirdBody, Lindemann, Troe)):
        if hasattr(reaction, 'specific_collider') and reaction.specific_collider:
            suffix = " + " + get_label(reaction.specific_collider, species_list)
        else:
            suffix = " (+ M)"

    arrow = " <=> "

    # Assemble Equation
    equation = reactants_str + suffix + arrow + products_str + suffix

    entry = {'equation': equation}

    # Write duplicate flag if present
    if reaction.duplicate:
        entry['duplicate'] = True

    # --- Kinetics Serialization ---

    if isinstance(kin, Arrhenius):
        entry['rate-constant'] = {'A': kin.A.value_si, 'b': kin.n.value_si, 'Ea': kin.Ea.value_si}

    elif isinstance(kin, Chebyshev):
        entry['type'] = 'Chebyshev'
        entry['temperature-range'] = [kin.Tmin.value_si, kin.Tmax.value_si]
        entry['pressure-range'] = [kin.Pmin.value_si, kin.Pmax.value_si]
        entry['data'] = kin.coeffs.value_si.tolist()

    elif isinstance(kin, ThirdBody):
        entry['type'] = 'three-body'
        entry['rate-constant'] = {
            'A': kin.arrheniusLow.A.value_si,
            'b': kin.arrheniusLow.n.value_si,
            'Ea': kin.arrheniusLow.Ea.value_si
        }
        entry['efficiencies'] = {lbl: v for m, v in kin.efficiencies.items() if
                                 (lbl := get_label(m, species_list)) is not None}

    elif isinstance(kin, Troe):
        entry['type'] = 'falloff'
        entry['high-P-rate-constant'] = {
            'A': kin.arrheniusHigh.A.value_si,
            'b': kin.arrheniusHigh.n.value_si,
            'Ea': kin.arrheniusHigh.Ea.value_si
        }
        entry['low-P-rate-constant'] = {
            'A': kin.arrheniusLow.A.value_si,
            'b': kin.arrheniusLow.n.value_si,
            'Ea': kin.arrheniusLow.Ea.value_si
        }
        troe_p = {'A': kin.alpha, 'T3': kin.T3.value_si, 'T1': kin.T1.value_si}
        if kin.T2:
            troe_p['T2'] = kin.T2.value_si
        entry['Troe'] = troe_p
        entry['efficiencies'] = {lbl: v for m, v in kin.efficiencies.items() if
                                 (lbl := get_label(m, species_list)) is not None}

    elif isinstance(kin, Lindemann):
        entry['type'] = 'falloff'
        entry['high-P-rate-constant'] = {
            'A': kin.arrheniusHigh.A.value_si,
            'b': kin.arrheniusHigh.n.value_si,
            'Ea': kin.arrheniusHigh.Ea.value_si
        }
        entry['low-P-rate-constant'] = {
            'A': kin.arrheniusLow.A.value_si,
            'b': kin.arrheniusLow.n.value_si,
            'Ea': kin.arrheniusLow.Ea.value_si
        }
        entry['efficiencies'] = {lbl: v for m, v in kin.efficiencies.items() if
                                 (lbl := get_label(m, species_list)) is not None}

    elif isinstance(kin, MultiArrhenius):
        entries = []
        for sub_kin in kin.arrhenius:
            # Create a temporary wrapper reaction for the sub-kinetic
            sub_rxn = Reaction(
                reactants=reaction.reactants,
                products=reaction.products,
                reversible=reaction.reversible,
                kinetics=sub_kin,
                duplicate=True  # MultiArrhenius always implies duplicates
            )
            # Recursively handle the sub-reaction
            entries.extend(reaction_to_dict_list(sub_rxn, species_list))
        return entries

    elif isinstance(kin, PDepArrhenius):
        # Check if any pressure point uses MultiArrhenius (sum of rates)
        has_multi = any(isinstance(arr, MultiArrhenius) for arr in kin.arrhenius)

        if has_multi:
            # We must split this complex PDep into multiple "duplicate" Cantera entries.
            # 1. Determine the maximum "depth" (max number of Arrhenius terms at any pressure)
            max_terms = 0
            for arr in kin.arrhenius:
                if isinstance(arr, MultiArrhenius):
                    max_terms = max(max_terms, len(arr.arrhenius))
                else:
                    max_terms = max(max_terms, 1)

            entries = []

            # 2. Create one YAML entry per "channel" (i = 0, 1, 2...)
            for i in range(max_terms):
                sub_entry = entry.copy()
                sub_entry['type'] = 'pressure-dependent-Arrhenius'
                sub_entry['duplicate'] = True

                rates = []
                for P, arr in zip(kin.pressures.value_si, kin.arrhenius):
                    current_arr = None

                    # Logic to extract the i-th Arrhenius term at this pressure
                    if isinstance(arr, MultiArrhenius):
                        if i < len(arr.arrhenius):
                            current_arr = arr.arrhenius[i]
                    elif isinstance(arr, Arrhenius):
                        if i == 0:
                            current_arr = arr

                    if current_arr:
                        rates.append({
                            'P': P,
                            'A': current_arr.A.value_si,
                            'b': current_arr.n.value_si,
                            'Ea': current_arr.Ea.value_si
                        })
                    else:
                        # If this channel has no rate at this pressure (e.g. P1 has 2 terms, P2 has 1),
                        # Cantera requires a value for interpolation. Use a negligible rate (A=0).
                        rates.append({'P': P, 'A': 0.0, 'b': 0.0, 'Ea': 0.0})

                sub_entry['rate-constants'] = rates
                entries.append(sub_entry)

            return entries

        else:
            # Standard Case: Simple Arrhenius at every pressure
            entry['type'] = 'pressure-dependent-Arrhenius'
            rates = []
            for P, arr in zip(kin.pressures.value_si, kin.arrhenius):
                rates.append({
                    'P': P,
                    'A': arr.A.value_si,
                    'b': arr.n.value_si,
                    'Ea': arr.Ea.value_si
                })
            entry['rate-constants'] = rates

    else:
        logging.warning(f"Skipping reaction {equation}: Unknown kinetics type {type(kin)}")
        return []

    note_parts = list()
    # A. Reaction Source (Provenance)
    if isinstance(reaction, TemplateReaction):
        note_parts.append(f"Source: Template family {reaction.family}")
    elif isinstance(reaction, LibraryReaction):
        note_parts.append(f"Source: Library {reaction.library}")
    elif isinstance(reaction, PDepReaction):
        note_parts.append(f"Source: PDep Network #{reaction.network.index}")
    elif isinstance(reaction, Reaction):
        note_parts.append(f"Source: P{reaction.kinetics.comment}")

    # B. Kinetics Comments (e.g. "Matched node 1234", "Flux pairs...", etc)
    if hasattr(kin, 'comment') and kin.comment:
        # Clean up newlines to keep the YAML one-line note clean
        clean_comment = kin.comment.replace('\n', '; ').strip()
        if clean_comment:
            note_parts.append(clean_comment)

    # C. Specific Collider info (if not obvious in equation)
    if reaction.specific_collider:
        note_parts.append(f"Specific collider: {reaction.specific_collider.label}")

    if note_parts:
        entry['note'] = " | ".join(note_parts)

    return [entry]


def get_label(obj: Union['Species', 'Molecule'], species_list: list['Species']):
    if species_list:
        for sp in species_list:
            if sp.is_isomorphic(obj):
                return f'{sp.label}({sp.index})' if sp.index > 0 else sp.label
    return None
