#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2024 Prof. William H. Green (whgreen@mit.edu),           #
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
This file defines functions for outputting the RMG generated mechanism to
a yaml file that can be read by Cantera
"""


import os
import shutil
import yaml
import logging

from rmgpy.species import Species
from rmgpy.kinetics.arrhenius import (
    MultiArrhenius,
    MultiPDepArrhenius,
)
from rmgpy.kinetics.falloff import ThirdBody
from rmgpy.util import make_output_subdirectory
from datetime import datetime
from rmgpy.chemkin import get_species_identifier


def _convert_anymap_to_dict(obj):
    """
    Recursively convert Cantera AnyMap objects to regular Python dicts.
    
    Cantera's input_data property returns dicts containing AnyMap objects,
    which are Cython extension types that cannot be serialized by YAML.
    This function recursively converts all AnyMaps to plain dicts.
    
    Args:
        obj: Any object (dict, list, AnyMap, or primitive type)
        
    Returns:
        The object with all AnyMaps converted to dicts
    """
    try:
        from cantera._utils import AnyMap
    except ImportError:
        # If Cantera is not available or doesn't have AnyMap, just return the object
        return obj
    
    if isinstance(obj, AnyMap):
        # Convert AnyMap to dict and recursively process values
        return {k: _convert_anymap_to_dict(v) for k, v in dict(obj).items()}
    elif isinstance(obj, dict):
        # Recursively process dict values
        return {k: _convert_anymap_to_dict(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        # Recursively process list/tuple elements
        return type(obj)(_convert_anymap_to_dict(item) for item in obj)
    else:
        # Return primitive types as-is
        return obj


def write_cantera(
    spcs,
    rxns,
    surface_site_density=None,
    solvent=None,
    solvent_data=None,
    path="chem.yml",
):
    """
    Writes yaml file depending on the type of system (gas-phase, catalysis).
    Writes beginning lines of yaml file, then uses yaml.dump(result_dict) to write species/reactions info.
    """

    # intro to file will change depending on the presence of surface species
    is_surface = False
    for spc in spcs:
        if spc.contains_surface_site():
            is_surface = True
    if is_surface:
        result_dict = get_mech_dict_surface(
            spcs, rxns, solvent=solvent, solvent_data=solvent_data
        )
        phases_block = get_phases_with_surface(
            spcs, surface_site_density
        )
    else:
        result_dict = get_mech_dict_nonsurface(
            spcs, rxns, solvent=solvent, solvent_data=solvent_data
        )
        phases_block = get_phases_gas_only(spcs)

    with open(path, "w") as f:
        # generator line
        f.write("generator: RMG\n")

        # datetime object containing current date and time
        now = datetime.now()
        dt_string = now.strftime("%a, %d %b %Y %H:%M:%S")
        f.write(f"date: {dt_string}\n")

        # units line
        f.write(
            "\nunits: {length: cm, time: s, quantity: mol, activation-energy: kcal/mol}\n\n"
        )

        f.write(phases_block)

        f.write(ELEMENTS_BLOCK)

        yaml.dump(result_dict, stream=f, sort_keys=False, default_flow_style=None, width=80)

def get_elements_block():
    """
    Returns the 'elements' section, and elements list for a phase
    """
    from rmgpy.molecule.element import get_element
    elements_list = ['H', 'C', 'O', 'N', 'Ne', 'Ar', 'He', 'Si', 'S',
                'F', 'Cl', 'Br', 'I']
    isotopes = (('H', 2), ('H', 3), ('C', 13),('O', 18))
    elements_block_list = ['', 'elements:']
    for symbol, isotope in isotopes:
        element = get_element(symbol, isotope=isotope)
        chemkin_name = element.chemkin_name
        mass = 1000 * element.mass
        elements_block_list.append(f"- symbol: {chemkin_name}\n  atomic-weight: {mass:f}")
        elements_list.append(chemkin_name)
    # Surface sites
    elements_list.append('X')
    elements_block_list.append("- symbol: X\n  atomic-weight: 195.083\n\n")
    elements_block = '\n'.join(elements_block_list)
    elements_line = f"elements: [{', '.join(elements_list)}]"
    return elements_block, elements_line
# For now this is not dynamic, and includes everything, so we just evaluate it 
# once and use it for all files.
ELEMENTS_BLOCK, ELEMENTS_LINE = get_elements_block()


def get_phases_gas_only(spcs):
    """
    Returns 'phases' sections for a file
    with only gas-phase species/reactions.
    """
    sorted_species = sorted(spcs, key=lambda spcs: spcs.index)
    species_to_write = [get_species_identifier(spec) for spec in sorted_species]
    # make sure species with "[" or "]" is in quotes
    species_to_write = [
        f"'{s}'" if "[" in s or "{" in s or "]" in s or "}" in s else s
        for s in species_to_write
    ]
    phases_block = f"""
phases:
- name: gas
  thermo: ideal-gas
  {ELEMENTS_LINE}
  species: [{', '.join(species_to_write)}]
  kinetics: gas
  transport: mixture-averaged
  state: {{T: 300.0, P: 1 atm}}
"""
    return phases_block


def get_phases_with_surface(spcs, surface_site_density):
    """
    Yaml files with surface species begin with the following blocks of text,
    which includes TWO phases instead of just one.
    Returns 'phases' sections.
    """
    surface_species = []
    gas_species = []
    for spc in spcs:

        if spc.contains_surface_site():
            surface_species.append(spc)
        else:
            gas_species.append(spc)

    sorted_surface_species = sorted(
        surface_species, key=lambda surface_species: surface_species.index
    )

    surface_species_to_write = [
        get_species_identifier(s) for s in sorted_surface_species
    ]

    # make sure species with "[" or "]" is in quotes
    surface_species_to_write = [
        f"'{s}'" if "[" in s or "{" in s or "]" in s or "}" in s else s
        for s in surface_species_to_write
    ]

    sorted_gas_species = sorted(gas_species, key=lambda gas_species: gas_species.index)
    gas_species_to_write = [get_species_identifier(s) for s in sorted_gas_species]

    # make sure species with "[" or "]" is in quotes
    gas_species_to_write = [
        f"'{s}'" if "[" in s or "{" in s or "]" in s or "}" in s else s
        for s in gas_species_to_write
    ]

    phases_block = f"""
phases:
- name: gas
  thermo: ideal-gas
  {ELEMENTS_LINE}
  species: [{', '.join(gas_species_to_write)}]
  kinetics: gas
  reactions: [gas-reactions]
  transport: mixture-averaged
  state: {{T: 300.0, P: 1 atm}}

- name: surface
  thermo: ideal-surface
  adjacent-phases: [gas]
  {ELEMENTS_LINE}
  species: [{', '.join(surface_species_to_write)}]
  kinetics: surface
  reactions: [site0-reactions]
  site-density: {surface_site_density * 1e-4 }
"""
    # surface_site_density * 1e-4 #in units of mol/cm^2

    return phases_block


def get_mech_dict_surface(spcs, rxns, solvent="solvent", solvent_data=None):
    """
    For systems with surface species/reactions.
    Adds 'species', 'gas-reactions', and 'site0-reactions' to result_dict.
    """
    gas_rxns = []
    surface_rxns = []
    for rxn in rxns:
        if rxn.is_surface_reaction():
            surface_rxns.append(rxn)
        else:
            gas_rxns.append(rxn)

    names = [x.label for x in spcs]
    for i, name in enumerate(names):  # fix duplicate names
        if names.count(name) > 1:
            names[i] += "-" + str(names.count(name))

    result_dict = dict()
    result_dict["species"] = [species_to_dict(x) for x in spcs]

    # separate gas and surface reactions

    gas_reactions = []
    for rmg_rxn in gas_rxns:
        gas_reactions.extend(reaction_to_dicts(rmg_rxn, spcs))
    result_dict["gas-reactions"] = gas_reactions

    surface_reactions = []
    for rmg_rxn in surface_rxns:
        surface_reactions.extend(reaction_to_dicts(rmg_rxn, spcs))
    result_dict["site0-reactions"] = surface_reactions

    return result_dict


def get_mech_dict_nonsurface(spcs, rxns, solvent="solvent", solvent_data=None):
    """
    For gas-phase systems.
    Adds 'species' and 'reactions' to result_dict.
    """
    names = [x.label for x in spcs]
    for i, name in enumerate(names):  # fix duplicate names
        if names.count(name) > 1:
            names[i] += "-" + str(names.count(name))

    result_dict = dict()
    result_dict["species"] = [species_to_dict(x) for x in spcs]

    reactions = []
    for rmg_rxn in rxns:
        reactions.extend(reaction_to_dicts(rmg_rxn, spcs))
    result_dict["reactions"] = reactions

    return result_dict


def _get_A_conversion_factor(n_reactants):
    """
    Get the conversion factor for the pre-exponential factor A from
    Cantera's SI default units to the declared YAML units
    (length: cm, quantity: mol).

    Cantera's input_data returns A in SI units (m, kmol, s).
    The YAML file declares units: {length: cm, quantity: mol}.

    The conversion depends on the reaction order (number of reactant
    molecules), NOT on rate_coeff_units (which is Units(0.0) for
    reactions created programmatically via to_cantera()).

    For rate constant units [length^(3*(n-1)) / quantity^(n-1) / time]:
      length: m -> cm  => multiply by (1e2)^(3*(n-1)) = 1e(6*(n-1))
      quantity: kmol -> mol => divide by (1e3)^(n-1) = 1e(3*(n-1))
      Combined: 1e(6*(n-1)) / 1e(3*(n-1)) = 1e(3*(n-1))

    Conversion factors by reaction order:
      - Unimolecular  (n=1): 1e0  = 1
      - Bimolecular   (n=2): 1e3  = 1000
      - Termolecular  (n=3): 1e6  = 1000000
    """
    order = max(n_reactants - 1, 0)
    return 10.0 ** (3 * order)


# Conversion factor for activation energy: J/kmol -> kcal/mol
_EA_CONVERSION_FACTOR = 1.0 / 4184000.0  # 4184 J/kcal * 1000 mol/kmol


def _convert_rate_constant_units(rate_dict, A_factor):
    """
    Convert a rate-constant dictionary {A, b, Ea} from Cantera SI defaults
    (m, kmol, J/kmol) to declared YAML units (cm, mol, kcal/mol).
    Modifies the dictionary in place.
    """
    if 'A' in rate_dict:
        rate_dict['A'] = rate_dict['A'] * A_factor
    if 'Ea' in rate_dict:
        rate_dict['Ea'] = rate_dict['Ea'] * _EA_CONVERSION_FACTOR


def _convert_reaction_data_units(reaction_data, n_reactants):
    """
    Convert all rate parameters in a reaction_data dict from Cantera SI
    defaults to the declared YAML units (cm, mol, kcal/mol).

    Handles simple Arrhenius (rate-constant), three-body, and
    falloff (high-P-rate-constant, low-P-rate-constant) reactions.

    n_reactants is the number of reactant molecules in the RMG reaction,
    used to determine the A conversion factor.
    """
    A_factor = _get_A_conversion_factor(n_reactants)

    if 'rate-constant' in reaction_data:
        _convert_rate_constant_units(reaction_data['rate-constant'], A_factor)
    if 'high-P-rate-constant' in reaction_data:
        _convert_rate_constant_units(reaction_data['high-P-rate-constant'], A_factor)
    if 'low-P-rate-constant' in reaction_data:
        # Low-P limit is one order higher in concentration than high-P
        low_P_A_factor = _get_A_conversion_factor(n_reactants + 1)
        _convert_rate_constant_units(reaction_data['low-P-rate-constant'], low_P_A_factor)


def reaction_to_dicts(obj, spcs):
    """
    Takes an RMG reaction object (obj), returns a list of dictionaries
    for YAML properties. For most reaction objects the list will be of
    length 1, but a MultiArrhenius or MultiPDepArrhenius will be longer.

    The returned dictionaries have rate parameters converted from Cantera's
    SI default units (m, kmol, J/kmol) to the declared YAML units
    (cm, mol, kcal/mol) so the YAML file is self-consistent.
    """

    reaction_list = []
    if isinstance(obj.kinetics, MultiArrhenius) or isinstance(
        obj.kinetics, MultiPDepArrhenius
    ):
        list_of_cantera_reactions = obj.to_cantera(use_chemkin_identifier=True)
    else:
        list_of_cantera_reactions = [obj.to_cantera(use_chemkin_identifier=True)]

    # Count reactant molecules from the RMG reaction object.
    # This is used to determine the A conversion factor since
    # rate_coeff_units is Units(0.0) for programmatically-created reactions.
    n_reactants = len(obj.reactants)

    # For three-body reactions (+ M), the third body M acts as an
    # additional reactant for unit purposes, so increment n_reactants.
    if isinstance(obj.kinetics, ThirdBody):
        n_reactants += 1

    for reaction in list_of_cantera_reactions:
        reaction_data = reaction.input_data
        _convert_reaction_data_units(reaction_data, n_reactants)
        efficiencies = getattr(obj.kinetics, "efficiencies", {})
        if efficiencies:
            reaction_data["efficiencies"] = {
                spcs[i].to_chemkin(): float(val)
                for i, val in enumerate(
                    obj.kinetics.get_effective_collider_efficiencies(spcs)
                )
                if val != 1
            }
        # Convert any AnyMap objects to regular dicts before appending
        reaction_data = _convert_anymap_to_dict(reaction_data)
        reaction_list.append(reaction_data)

    return reaction_list


def species_to_dict(species):
    """
    Takes an RMG species object, returns a list of dictionaries
    for YAML properties. Also adds in the number of surface sites
    ('sites') to dictionary.
    """
    if not isinstance(species, Species):
        raise TypeError("species object must be an RMG Species")

    cantera_species = species.to_cantera(use_chemkin_identifier=True)
    species_data = cantera_species.input_data

    try:
        transport_comment = species.transport_data.comment
        if transport_comment:
            species_data["transport"]["note"] = transport_comment
    except AttributeError:
        pass

    if "size" in species_data:
        sites = species_data["size"]
        species_data.pop("size", None)
        species_data["sites"] = sites

    # Convert any AnyMap objects to regular dicts before returning
    species_data = _convert_anymap_to_dict(species_data)

     # returns composition, name, thermo, and transport, and note
    return species_data


class CanteraWriter(object):
    """
    This class listens to a RMG subject
    and writes an YAML file with the current state of the RMG model,
    to a yaml subfolder.


    A new instance of the class can be appended to a subject as follows:

    rmg = ...
    listener = CanteraWriter(outputDirectory)
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)

    """

    def __init__(self, output_directory=""):
        super(CanteraWriter, self).__init__()
        self.output_directory = output_directory
        self.output_subdirectory = os.path.join(self.output_directory, "cantera")
        make_output_subdirectory(output_directory, "cantera")

    def update(self, rmg):

        this_output_path = os.path.join(self.output_subdirectory,
                                        f"chem{len(rmg.reaction_model.core.species):04d}.yaml")
        latest_output_path = os.path.join(self.output_subdirectory, 'chem.yaml')

        logging.info(f"Saving current model core to Cantera file: {this_output_path}")

        solvent_data = None
        if rmg.solvent:
            solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)

        surface_site_density = None
        if rmg.reaction_model.surface_site_density:
            surface_site_density = rmg.reaction_model.surface_site_density.value_si

        write_cantera(
            rmg.reaction_model.core.species,
            rmg.reaction_model.core.reactions,
            surface_site_density=surface_site_density,
            solvent=rmg.solvent,
            solvent_data=solvent_data,
            path=this_output_path
        )
        # Update the latest output path
        shutil.copy2(this_output_path, latest_output_path)
