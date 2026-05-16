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
try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper
import yaml

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import (
    Arrhenius, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius,
    Chebyshev, Troe, Lindemann, ThirdBody,
    StickingCoefficient, SurfaceArrhenius,
)
from rmgpy.reaction import Reaction
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.util import make_output_subdirectory
import rmgpy.constants as constants

from rmgpy.species import Species
if TYPE_CHECKING:
    from rmgpy.molecule.molecule import Molecule

SYMBOL_BY_NUMBER = {0: 'e', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
                    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
                    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
                    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
                    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
                    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
                    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
                    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
                    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
                    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
                    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt',
                    110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}
NUMBER_BY_SYMBOL = {value: key for key, value in SYMBOL_BY_NUMBER.items()}


class CanteraWriter2(object):
    """
    This class listens to a RMG subject and writes a Cantera YAML file
    with the current state of the RMG model at every iteration.
    """

    def __init__(self, output_directory='', config=None):
        self.output_directory = output_directory
        self.config = config
        make_output_subdirectory(output_directory, 'cantera2')

    def update(self, rmg):
        """
        Called whenever the RMG subject notifies listeners.
        """
        if self.config is not None and not self.config.should_write(
                rmg.reaction_model.iteration_num, rmg.is_final_save):
            return
        save_cantera_files(rmg, config=self.config)


def save_cantera_files(rmg, config=None):
    """
    Save the current reaction model to a set of Cantera YAML files.

    Creates:
      1. chem{N}.yaml (where N is num species)
      2. chem.yaml (latest copy)
      3. chem_annotated.yaml (if verbose_comments is True)
    """
    verbose = config.verbose_comments if (config and config.verbose_comments is not None) else rmg.verbose_comments
    save_edge = config.save_edge if (config and config.save_edge is not None) else rmg.save_edge_species

    # Ensure subdirectory exists
    cantera_dir = os.path.join(rmg.output_directory, 'cantera2')
    if not os.path.exists(cantera_dir):
        os.mkdir(cantera_dir)

    try:
        site_density = rmg.surface_site_density.value_si
    except (AttributeError, KeyError, TypeError):
        site_density = None

    # -------------------------------------------------------------------------
    # 1. Save Core Model
    # -------------------------------------------------------------------------
    num_species = len(rmg.reaction_model.core.species)

    # Define paths
    this_cantera_path = os.path.join(cantera_dir,
                                     'chem{0:04d}.yaml'.format(num_species))
    latest_cantera_path = os.path.join(cantera_dir, 'chem.yaml')

    logging.info(f"Saving current model core to Cantera file: {this_cantera_path}")

    # Write the YAML file (non-verbose)
    save_cantera_model(rmg.reaction_model.core, this_cantera_path, site_density=site_density,
                       verbose=False)

    # Copy to 'chem.yaml' (The latest file)
    if os.path.exists(latest_cantera_path):
        os.unlink(latest_cantera_path)
    shutil.copy2(this_cantera_path, latest_cantera_path)

    # Write annotated file if verbose_comments is requested
    if verbose:
        annotated_path = os.path.join(cantera_dir, 'chem_annotated.yaml')
        logging.info(f"Saving annotated Cantera file: {annotated_path}")
        save_cantera_model(rmg.reaction_model.core, annotated_path, site_density=site_density,
                           verbose=True)

    # -------------------------------------------------------------------------
    # 2. Save Edge Model (Optional, matching ChemkinWriter logic)
    # -------------------------------------------------------------------------
    if save_edge:
        from rmgpy.rmg.model import ReactionModel
        logging.info('Saving current model core and edge to Cantera file...')

        this_edge_path = os.path.join(cantera_dir,
                                      'chem_edge{0:04d}.yaml'.format(num_species))
        latest_edge_path = os.path.join(cantera_dir, 'chem_edge.yaml')

        edge_model = ReactionModel(
            species=rmg.reaction_model.core.species + rmg.reaction_model.edge.species,
            reactions=rmg.reaction_model.core.reactions + rmg.reaction_model.edge.reactions,
        )

        save_cantera_model(edge_model, this_edge_path, site_density=site_density, verbose=False)

        if os.path.exists(latest_edge_path):
            os.unlink(latest_edge_path)
        shutil.copy2(this_edge_path, latest_edge_path)

        if verbose:
            annotated_edge_path = os.path.join(cantera_dir, 'chem_edge_annotated.yaml')
            logging.info(f"Saving annotated edge Cantera file: {annotated_edge_path}")
            save_cantera_model(edge_model, annotated_edge_path, site_density=site_density,
                               verbose=True)


def save_cantera_model(model_container, path, site_density=None, verbose=False):
    """
    Internal helper to generate the dictionary and write the YAML file.
    model_container must be a :class:`rmgpy.rmg.model.ReactionModel` (or duck-typed
    equivalent with .species, .reactions, and .get_elements()).
    If verbose=True, species/reaction notes (SMILES, source, kinetics
    comments) are included in the output.
    """
    species_list = model_container.species
    reaction_list = model_container.reactions
    elements_in_use = model_container.get_elements()

    is_plasma = False
    for sp in species_list:
        if sp.is_electron():
            is_plasma = True
            break

    # Generate Data
    yaml_data = generate_cantera_data(species_list, reaction_list,
                                      elements_in_use=elements_in_use,
                                      is_plasma=is_plasma,
                                      site_density=site_density, verbose=verbose)

    # Write
    with open(path, 'w') as f:
        # sort_keys=False ensures 'units' comes first, then 'phases', etc.
        yaml.dump(yaml_data, f, Dumper=Dumper, sort_keys=False, default_flow_style=None)

def get_elements_lists(elements_in_use):
    """
    Returns custom element definitions and the elements list for phases.

    elements_in_use is a set of :class:`Element` singletons (typically from
    :meth:`rmgpy.rmg.model.ReactionModel.get_elements`). Only those elements
    are emitted; isotopes (D, T, CI, OI) and X are added only when present in
    the set. The plasma pseudo-element 'E' is added separately by the caller
    when ``is_plasma`` is true.
    """
    from rmgpy.molecule.element import H, C, O, N, Ne, Ar, He, Si, S, F, Cl, Br, I, D, T, C13, O18, X
    builtin_elements = [(H, 'H'), (C, 'C'), (O, 'O'), (N, 'N'), (Ne, 'Ne'), (Ar, 'Ar'),
                        (He, 'He'), (Si, 'Si'), (S, 'S'), (F, 'F'), (Cl, 'Cl'), (Br, 'Br'), (I, 'I')]
    elements_list = [symbol for element, symbol in builtin_elements if element in elements_in_use]
    custom_elements = []
    for isotope in (D, T, C13, O18):
        if isotope in elements_in_use:
            mass = 1000 * isotope.mass
            custom_elements.append({'symbol': isotope.chemkin_name, 'atomic-weight': mass})
            elements_list.append(isotope.chemkin_name)
    if X in elements_in_use:
        elements_list.append('X')
        custom_elements.append({'symbol': 'X', 'atomic-weight': 195.083})
    return custom_elements, elements_list

def generate_cantera_data(species_list,
                          reaction_list,
                          elements_in_use=None,
                          is_plasma=False,
                          site_density=None,
                          verbose=False,
                          ):
    """
    Converts RMG objects into a dictionary structure compatible with Cantera YAML.
    If verbose=True, species/reaction notes are included (SMILES, source, kinetics comments).

    elements_in_use is a set of :class:`Element` singletons (typically from
    :meth:`rmgpy.rmg.model.ReactionModel.get_elements`) used to size the
    'elements' block and the per-phase elements lists. Defaults to an empty
    set if None.
    """
    if elements_in_use is None:
        elements_in_use = set()
    # --- 1. Header & Units ---
    # We output everything in SI units.
    try:
        from rmgpy.rmg.main import RMG
        git_head, _ = RMG.get_git_commit(None, os.path.dirname(__file__))
        git_head = " (git commit: {0})".format(git_head[:7])
    except Exception:
        git_head = ''

    data = {
        'description': 'RMG-Py Generated Mechanism',
        'generator': f'RMG-Py CanteraWriter2 at {__file__}{git_head}',
        'cantera-version': '3.1',
        'units': {
            'length': 'm',
            'time': 's',
            'quantity': 'mol',
            'activation-energy': 'J/mol'
        }
    }

    # Sort species list by index
    sorted_species = sorted(species_list, key=lambda species: species.index)

    # --- 2. Phase Segregation (Gas vs Surface) ---
    gas_species, surface_species, gas_reactions, surface_reactions = list(), list(), list(), list()

    for spc in sorted_species:
        if spc.contains_surface_site():
            surface_species.append(spc)
        else:
            gas_species.append(spc)

    for rxn in reaction_list:
        if rxn.is_surface_reaction():
            surface_reactions.append(rxn)
        else:
            gas_reactions.append(rxn)

    # --- 3. Phase Definitions ---
    custom_elements, all_elements = get_elements_lists(elements_in_use)

    if custom_elements:
        data['elements'] = custom_elements
    elements_set = set(all_elements)
    if is_plasma:
        elements_set.add('E')

    phases = list()

    gas_phase_def = {
        'name': 'gas',
        'thermo': 'plasma' if is_plasma else 'ideal-gas',
        'elements': sorted(list(elements_set)),
        'species': [get_label(spc, species_list) for spc in gas_species],
        'kinetics': 'gas',
    }

    if surface_species:
        gas_phase_def['reactions'] = ['gas-reactions']
    if is_plasma:
        gas_phase_def['transport'] = 'ionized-gas'
        # Plasma specific defaults
        gas_phase_def['electron-energy-distribution'] = {
            'type': 'isotropic',
            'shape-factor': 2.0,
            'mean-electron-energy': 1.0
        }
    else:
        gas_phase_def['transport'] = 'mixture-averaged'

    gas_phase_def['state'] = {'T': 300.0, 'P': '1 atm'}
    phases.append(gas_phase_def)

    if surface_species:
        default_site_density = 2.5e-5  # mol/m^2

        has_coverage_dependence = any(
            hasattr(sp.thermo, 'thermo_coverage_dependence') and sp.thermo.thermo_coverage_dependence
            for sp in surface_species
        )

        surface_phase_def = {
            'name': 'surface',
            'thermo': 'coverage-dependent-surface' if has_coverage_dependence else 'ideal-surface',
            'adjacent-phases': ['gas'],
            'elements': sorted(list(elements_set)),
            'species': [get_label(sp, species_list) for sp in surface_species],
            'kinetics': 'surface',
            'reactions': ['surface-reactions'],
            'site-density': site_density or default_site_density,
            'state': {'T': 300.0, 'P': '1 atm'},
        }
        if has_coverage_dependence:
            surface_phase_def['reference-state-coverage'] = 0.11
        phases.append(surface_phase_def)

    data['phases'] = phases

    species_data = list()
    for sp in species_list:
        species_data.append(species_to_dict(sp, species_list, verbose=verbose))
    data['species'] = species_data

    # Build separate reaction lists for each phase if there are two phases
    gas_reaction_data = list()
    for rxn in gas_reactions:
        entries = reaction_to_dict_list(rxn, species_list, verbose=verbose)
        if entries:
            gas_reaction_data.extend(entries)

    if surface_species:
        data['gas-reactions'] = gas_reaction_data
    else:
        data['reactions'] = gas_reaction_data

    if surface_reactions:
        surface_reaction_data = list()
        for rxn in surface_reactions:
            entries = reaction_to_dict_list(rxn, species_list, verbose=verbose)
            if entries:
                surface_reaction_data.extend(entries)
        data['surface-reactions'] = surface_reaction_data

    return data


def species_to_dict(species, species_list, verbose=False):
    """Convert an RMG Species object to a Cantera YAML dictionary.
    If verbose=True, species notes (SMILES, thermo/transport comments) are included.
    """

    notes = list()
    if verbose:
        try:
            notes.append(species.to_smiles())
        except:
            pass

    # Composition
    mol = species.molecule[0]
    atom_dict = dict(mol.get_element_count())

    # Number of electrons 'E'
    # The special pseudo-element E is used in representing charged species, where it specifies
    # the net number of electrons compared to the number needed to form a neutral species.
    # That is, negatively charged ions will have E > 0, while positively charged ions will have E < 0.
    # https://cantera.org/3.1/userguide/creating-mechanisms.html#elemental-composition
    charge = mol.get_net_charge()
    if 'E' not in atom_dict and charge != 0:
        atom_dict['E'] = -charge

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

    num_sites = mol.number_of_surface_sites()
    if num_sites > 1:
        species_entry['sites'] = num_sites

    # Transport (if available) - Only relevant for gas phase usually
    if species.transport_data and not species.contains_surface_site():
        td = species.transport_data

        transport_dict = {
            'model': 'gas',
            'geometry': 'atom' if td.shapeIndex == 0 else 'linear' if td.shapeIndex == 1 else 'nonlinear',
            'well-depth': td.epsilon.value_si / constants.R, # Kelvin
            'diameter': td.sigma.value_si * 1e10,  # Angstroms
        }
        if td.dipoleMoment and td.dipoleMoment.value_si != 0.0:
            transport_dict['dipole'] = td.dipoleMoment.value_si * 1e21 * constants.c  # Debye
        if getattr(td, 'polarizability', None) and td.polarizability.value_si != 0.0:
            transport_dict['polarizability'] = td.polarizability.value_si * 1e30  # Angstrom^3
        if getattr(td, 'rotrelaxcollnum', None) and td.rotrelaxcollnum != 0.0:
            transport_dict['rotational-relaxation'] = td.rotrelaxcollnum
        if verbose and td.comment:
            transport_dict['note'] = td.comment.strip()
        species_entry['transport'] = transport_dict

    if verbose and species.thermo and species.thermo.comment:
        clean_comment = species.thermo.comment.replace('\n', '; ').strip()
        species_entry['thermo']['note'] = clean_comment

    if notes:
        species_entry['note'] = " | ".join(notes)

    # Add coverage-dependencies if this surface species has coverage-dependent thermo
    if (species.contains_surface_site() and
            hasattr(thermo_data, 'thermo_coverage_dependence') and
            thermo_data.thermo_coverage_dependence):
        from rmgpy.molecule.molecule import Molecule
        cov_deps = {}
        for adj_list, parameters in thermo_data.thermo_coverage_dependence.items():
            mol = Molecule().from_adjacency_list(adj_list)
            for sp in species_list:
                if sp.is_isomorphic(mol, strict=False):
                    dep_label = get_label(sp, species_list)
                    if dep_label:
                        cov_deps[dep_label] = {
                            'units': {'energy': 'J', 'quantity': 'mol'},
                            'enthalpy-coefficients': [v.value_si for v in parameters['enthalpy-coefficients']],
                            'entropy-coefficients': [v.value_si for v in parameters['entropy-coefficients']],
                        }
                    break
        if cov_deps:
            species_entry['coverage-dependencies'] = cov_deps

    return species_entry


def reaction_to_dict_list(reaction, species_list=None, verbose=False):
    """
    Convert an RMG Reaction object to a LIST of Cantera YAML dictionaries.
    If verbose=True, a 'note' field is added with source and kinetics comment.
    """
    # Check for MultiKinetics (duplicates grouped in one RMG object)
    if isinstance(reaction.kinetics, (MultiArrhenius, MultiPDepArrhenius)):
        entries = []
        sub_kinetics_list = reaction.kinetics.arrhenius

        for sub_kin in sub_kinetics_list:
            sub_rxn = Reaction(
                reactants=reaction.reactants,
                products=reaction.products,
                reversible=reaction.reversible,
                kinetics=sub_kin,
                duplicate=True
            )
            sub_result = reaction_to_dict_list(sub_rxn, species_list, verbose=verbose)
            if sub_result:
                entries.extend(sub_result)
        return entries

    kin = reaction.kinetics

    # Generate equation string
    equation = get_reaction_equation(reaction, species_list)
    entry = {'equation': equation}

    if reaction.duplicate:
        entry['duplicate'] = True

    # --- Kinetics Serialization ---

    # 1. Surface Kinetics
    if isinstance(kin, StickingCoefficient):
        entry['type'] = 'sticking-Arrhenius'
        entry['sticking-coefficient'] = {'A': kin.A.value_si, 'b': kin.n.value_si, 'Ea': kin.Ea.value_si}

    elif isinstance(kin, SurfaceArrhenius):
        entry['type'] = 'interface-Arrhenius'
        entry['rate-constant'] = {'A': kin.A.value_si, 'b': kin.n.value_si, 'Ea': kin.Ea.value_si}

    # 2. Gas Kinetics
    elif isinstance(kin, Arrhenius):
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

    elif isinstance(kin, PDepArrhenius):
        # A MultiArrhenius entry comes from a chemkin PLOG block with duplicate
        # pressures; expand it into one rate-constants entry per inner Arrhenius
        # at the shared pressure. Cantera's pressure-dependent-Arrhenius sums
        # duplicate-pressure entries at evaluation, matching chemkin semantics.
        # Splitting into separate `duplicate: true` reactions would be wrong:
        # that interpolates each component independently, then sums, which
        # gives different rates at intermediate pressures.
        entry['type'] = 'pressure-dependent-Arrhenius'
        rates = []
        for P, arr in zip(kin.pressures.value_si, kin.arrhenius):
            sub_arrhenius = arr.arrhenius if isinstance(arr, MultiArrhenius) else [arr]
            for sub in sub_arrhenius:
                rates.append({
                    'P': float(P),
                    'A': sub.A.value_si,
                    'b': sub.n.value_si,
                    'Ea': sub.Ea.value_si
                })
        entry['rate-constants'] = rates

    else:
        logging.warning(f"Skipping reaction {equation}: Unknown kinetics type {type(kin)}")
        return []

    # --- Coverage Dependencies ---
    if hasattr(kin, 'coverage_dependence') and kin.coverage_dependence:
        cov_deps = {}
        for sp, cov_params in kin.coverage_dependence.items():
            sp_label = get_label(sp, species_list)
            if sp_label:
                # Cantera YAML expects { a: ..., m: ..., E: ... }
                cov_deps[sp_label] = {
                    'a': cov_params['a'].value_si,
                    'm': cov_params['m'].value_si,
                    'E': cov_params['E'].value_si
                }
        if cov_deps:
            entry['coverage-dependencies'] = cov_deps

    # --- Metadata / Notes (only when verbose) ---
    if verbose:
        note_parts = list()
        if isinstance(reaction, TemplateReaction):
            note_parts.append(f"Source: Template family {reaction.family}")
        elif isinstance(reaction, LibraryReaction):
            note_parts.append(f"Source: Library {reaction.library}")
        elif isinstance(reaction, PDepReaction):
            note_parts.append(f"Source: PDep Network #{reaction.network.index}")

        if hasattr(kin, 'comment') and kin.comment:
            clean_comment = kin.comment.replace('\n', '; ').strip()
            if clean_comment:
                note_parts.append(clean_comment)

        if reaction.specific_collider:
            note_parts.append(f"Specific collider: {reaction.specific_collider.label}")

        if note_parts:
            entry['note'] = " | ".join(note_parts)

    return [entry]


def get_reaction_equation(reaction, species_list):
    """Helper to build reaction string"""
    reactants_str = " + ".join([get_label(r, species_list) for r in reaction.reactants])
    products_str = " + ".join([get_label(p, species_list) for p in reaction.products])

    suffix = ""
    kin = reaction.kinetics
    collider = getattr(reaction, 'specific_collider', None)
    if isinstance(kin, ThirdBody):
        # Real three-body reaction: M (or specific collider) participates as a
        # reactant on both sides without parentheses.
        m_label = get_label(collider, species_list) if collider else "M"
        suffix = " + " + m_label
    elif isinstance(kin, (Lindemann, Troe)):
        # Pressure-dependent falloff: M acts as a chaperone, written in parens.
        m_label = get_label(collider, species_list) if collider else "M"
        suffix = " (+" + m_label + ")"

    arrow = " <=> " if reaction.reversible else " => "
    return reactants_str + suffix + arrow + products_str + suffix


def get_label(obj: Union['Species', 'Molecule'], species_list: list['Species']):
    if isinstance(obj, Species):
        return f'{obj.label}({obj.index})' if obj.index > 0 else obj.label

    if species_list:
        for sp in species_list:
            if sp.is_isomorphic(obj):
                return f'{sp.label}({sp.index})' if sp.index > 0 else sp.label
    return None
