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
try:
    from yaml import CDumper as _BaseDumper
except ImportError:
    from yaml import Dumper as _BaseDumper
import yaml
import logging


class Dumper(_BaseDumper):
    """
    YAML Dumper that emits multi-line strings using the literal block style
    ('|') instead of the default double-quoted form with embedded '\\n'.
    Used to render the multi-line 'note' field on reactions in a way that
    mirrors ck2yaml's output.
    """


def _multiline_str_representer(dumper, data):
    if "\n" in data:
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
    return dumper.represent_scalar("tag:yaml.org,2002:str", data)


Dumper.add_representer(str, _multiline_str_representer)

from rmgpy.species import Species
from rmgpy.kinetics.arrhenius import (
    MultiArrhenius,
    MultiPDepArrhenius,
)
from rmgpy.kinetics.falloff import Lindemann, ThirdBody, Troe
from rmgpy.kinetics.model import PDepKineticsModel
from rmgpy.util import make_output_subdirectory
from datetime import datetime
from rmgpy.chemkin import get_species_identifier
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.rmg.pdep import PDepReaction


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
    elements_in_use,
    surface_site_density=None,
    solvent=None,
    solvent_data=None,
    path="chem.yml",
    verbose=False,
):
    """
    Writes yaml file depending on the type of system (gas-phase, catalysis).
    Writes beginning lines of yaml file, then uses yaml.dump(result_dict) to write species/reactions info.
    If verbose=True, species and reaction notes (SMILES, source, kinetics comment) are included.

    elements_in_use is a set of :class:`Element` singletons. Only those elements
    are listed in the YAML 'elements' block and 'phases.elements' lines.
    """

    try:
        from rmgpy.rmg.main import RMG
        git_head, _ = RMG.get_git_commit(None, os.path.dirname(__file__))
        git_head = " (git commit: {0})".format(git_head[:7])
    except Exception:
        git_head = ''

    elements_block, elements_line = get_elements_block(elements_in_use)

    # intro to file will change depending on the presence of surface species
    is_surface = False
    for spc in spcs:
        if spc.contains_surface_site():
            is_surface = True
    if is_surface:
        has_coverage_dependence = any(
            hasattr(spc.thermo, 'thermo_coverage_dependence') and spc.thermo.thermo_coverage_dependence
            for spc in spcs if spc.contains_surface_site()
        )
        result_dict = get_mech_dict_surface(
            spcs, rxns, solvent=solvent, solvent_data=solvent_data, verbose=verbose
        )
        phases_block = get_phases_with_surface(
            spcs, surface_site_density, elements_line, has_coverage_dependence=has_coverage_dependence
        )
    else:
        result_dict = get_mech_dict_nonsurface(
            spcs, rxns, solvent=solvent, solvent_data=solvent_data, verbose=verbose
        )
        phases_block = get_phases_gas_only(spcs, elements_line)

    with open(path, "w") as f:
        # generator line
        generator = f"RMG-Py CanteraWriter1 at {__file__}{git_head}"
        f.write(f'generator: "{generator}"\n')

        # datetime object containing current date and time
        now = datetime.now()
        dt_string = now.strftime("%a, %d %b %Y %H:%M:%S")
        f.write(f"date: {dt_string}\n")

        # units line
        f.write(
            "\nunits: {length: m, time: s, quantity: kmol, activation-energy: J/kmol}\n\n"
        )

        f.write(phases_block)

        f.write(elements_block)

        yaml.dump(result_dict, stream=f, Dumper=Dumper, sort_keys=False, default_flow_style=None, width=80)

def get_elements_block(elements_in_use):
    """
    Returns the 'elements' section, and elements list for a phase.

    elements_in_use is a set of :class:`Element` singletons (e.g. the ones returned by
    :meth:`rmgpy.rmg.model.ReactionModel.get_elements`). Only elements present
    in the set are emitted; isotopes (D, T, CI, OI) and the surface site X are
    written to the elements block only when actually used.
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
    # Only emit the top-level 'elements:' block when there are non-builtin entries
    if custom_elements:
        elements_block = '\nelements:\n' + '\n'.join([f"- symbol: {e['symbol']}\n  atomic-weight: {e['atomic-weight']:f}" for e in custom_elements]) + '\n\n'
    else:
        elements_block = ''
    elements_line = f"elements: [{', '.join(elements_list)}]"
    return elements_block, elements_line


def get_phases_gas_only(spcs, elements_line):
    """
    Returns 'phases' sections for a file
    with only gas-phase species/reactions.

    elements_line is the pre-formatted ``elements: [...]`` string from
    :func:`get_elements_block`.
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
  {elements_line}
  species: [{', '.join(species_to_write)}]
  kinetics: gas
  transport: mixture-averaged
  state: {{T: 300.0, P: 1 atm}}
"""
    return phases_block


def get_phases_with_surface(spcs, surface_site_density, elements_line, has_coverage_dependence=False):
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

    surface_thermo = 'coverage-dependent-surface' if has_coverage_dependence else 'ideal-surface'
    reference_state_line = '\n  reference-state-coverage: 0.11' if has_coverage_dependence else ''

    phases_block = f"""
phases:
- name: gas
  thermo: ideal-gas
  {elements_line}
  species: [{', '.join(gas_species_to_write)}]
  kinetics: gas
  reactions: [gas-reactions]
  transport: mixture-averaged
  state: {{T: 300.0, P: 1 atm}}

- name: surface
  thermo: {surface_thermo}{reference_state_line}
  adjacent-phases: [gas]
  {elements_line}
  species: [{', '.join(surface_species_to_write)}]
  kinetics: surface
  reactions: [site0-reactions]
  site-density: {surface_site_density * 1e-4 }
  state: {{T: 300.0, P: 1 atm}}
"""
    # surface_site_density * 1e-4 #in units of mol/cm^2

    return phases_block


def get_mech_dict_surface(spcs, rxns, solvent="solvent", solvent_data=None, verbose=False):
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
    result_dict["species"] = [species_to_dict(x, all_species=spcs, verbose=verbose) for x in spcs]

    # separate gas and surface reactions

    gas_reactions = []
    for rmg_rxn in gas_rxns:
        gas_reactions.extend(reaction_to_dicts(rmg_rxn, spcs, verbose=verbose))
    result_dict["gas-reactions"] = gas_reactions

    surface_reactions = []
    for rmg_rxn in surface_rxns:
        surface_reactions.extend(reaction_to_dicts(rmg_rxn, spcs, verbose=verbose))
    result_dict["site0-reactions"] = surface_reactions

    return result_dict


def get_mech_dict_nonsurface(spcs, rxns, solvent="solvent", solvent_data=None, verbose=False):
    """
    For gas-phase systems.
    Adds 'species' and 'reactions' to result_dict.
    """
    names = [x.label for x in spcs]
    for i, name in enumerate(names):  # fix duplicate names
        if names.count(name) > 1:
            names[i] += "-" + str(names.count(name))

    result_dict = dict()
    result_dict["species"] = [species_to_dict(x, verbose=verbose) for x in spcs]

    reactions = []
    for rmg_rxn in rxns:
        reactions.extend(reaction_to_dicts(rmg_rxn, spcs, verbose=verbose))
    result_dict["reactions"] = reactions

    return result_dict


def _build_equation_string(obj):
    """
    Build the reaction equation string preserving the order of reactants and
    products as stored on the RMG Reaction object. Cantera's input_data sorts
    reactant/product maps internally, which loses the source ordering (e.g.
    'H + O2' instead of 'O2 + H' for HO2 formation). Match the equation
    convention used by ck2yaml/CanteraWriter2: stoichiometry coefficients are
    not collapsed, third-body M (or specific collider) is appended without
    parentheses, and falloff colliders are written as '(+M)'.
    """
    reactants = " + ".join(r.to_chemkin() for r in obj.reactants)
    products = " + ".join(p.to_chemkin() for p in obj.products)

    suffix = ""
    kin = obj.kinetics
    collider = getattr(obj, "specific_collider", None)
    if isinstance(kin, ThirdBody) and not isinstance(kin, (Lindemann, Troe)):
        m_label = collider.to_chemkin() if collider else "M"
        suffix = " + " + m_label
    elif isinstance(kin, (Lindemann, Troe)):
        m_label = collider.to_chemkin() if collider else "M"
        suffix = " (+" + m_label + ")"

    arrow = " <=> " if obj.reversible else " => "
    return reactants + suffix + arrow + products + suffix


def reaction_to_dicts(obj, spcs, verbose=False):
    """
    Takes an RMG reaction object (obj), returns a list of dictionaries
    for YAML properties. For most reaction objects the list will be of
    length 1, but a MultiArrhenius or MultiPDepArrhenius will be longer.
    If verbose=True, a 'note' field is added with source and kinetics comment.
    """

    reaction_list = []
    if isinstance(obj.kinetics, MultiArrhenius) or isinstance(
        obj.kinetics, MultiPDepArrhenius
    ):
        list_of_cantera_reactions = obj.to_cantera(use_chemkin_identifier=True)
    else:
        list_of_cantera_reactions = [obj.to_cantera(use_chemkin_identifier=True)]

    is_third_body = isinstance(obj.kinetics, PDepKineticsModel)

    rmg_equation = _build_equation_string(obj)

    for reaction in list_of_cantera_reactions:
        reaction_data = reaction.input_data
        # Cantera reorders reactant and product species (e.g. it writes
        # 'H + O2' even when the RMG reaction has them in the order O2, H),
        # and collapses repeated species into stoichiometric coefficients.
        # Overwrite with an equation built from obj.reactants/products to
        # preserve the source ordering and match ck2yaml's output style.
        reaction_data["equation"] = rmg_equation
        # Cantera's input_data omits 'type: three-body' for plain ThirdBody
        # reactions (only Lindemann/Troe falloff get a 'type' field). Add it
        # explicitly so the YAML matches what ck2yaml emits for the same input.
        if isinstance(obj.kinetics, ThirdBody) and "type" not in reaction_data:
            new_data = {"equation": reaction_data["equation"], "type": "three-body"}
            for k, v in reaction_data.items():
                if k != "equation":
                    new_data[k] = v
            reaction_data = new_data
        efficiencies = getattr(obj.kinetics, "efficiencies", {})
        if efficiencies:
            reaction_data["efficiencies"] = {
                spcs[i].to_chemkin(): float(val)
                for i, val in enumerate(
                    obj.kinetics.get_effective_collider_efficiencies(spcs)
                )
                if val != 1
            }
        elif not is_third_body:
            # Cantera's API misidentifies a species that appears on both sides
            # of a reaction (e.g. vacantX) as a third-body collider
            # when there are three or more species on one side, producing a
            # spurious 'efficiencies' entry in input_data.
            # see https://github.com/Cantera/cantera/issues/2115
            reaction_data.pop("efficiencies", None)
        # Convert any AnyMap objects to regular dicts before appending
        reaction_data = _convert_anymap_to_dict(reaction_data)

        if verbose:
            note_lines = []
            if isinstance(obj, TemplateReaction):
                note_lines.append(f"Template reaction: {obj.family}")
            elif isinstance(obj, LibraryReaction):
                note_lines.append(f"Library reaction: {obj.library}")
            elif isinstance(obj, PDepReaction):
                note_lines.append(f"PDep reaction: {obj.network}")
            if obj.specific_collider is not None:
                note_lines.append(f"Specific collider: {obj.specific_collider.label}")
            if obj.kinetics.comment:
                # Preserve the original line structure of the kinetics
                # comment (one line per source line) and right-strip each
                # line: PyYAML refuses the '|' literal block style for any
                # value whose lines have trailing whitespace.
                for line in obj.kinetics.comment.strip("\n").split("\n"):
                    note_lines.append(line.rstrip())
            if note_lines:
                # Trailing '\n' keeps PyYAML's literal block style ('|'
                # rather than '|-') so the rendered note ends with a
                # newline, matching ck2yaml's output.
                reaction_data["note"] = "\n".join(line.rstrip() for line in note_lines) + "\n"

        reaction_list.append(reaction_data)

    return reaction_list


def species_to_dict(species, all_species=None, verbose=False):
    """
    Takes an RMG species object, returns a dictionary of YAML properties.
    Also adds in the number of surface sites ('sites') to the dictionary.

    all_species: if provided, coverage-dependent thermo is resolved and
    attached to the Cantera species object before serialisation, so it
    appears in the returned dict automatically.
    If verbose=True, species SMILES and thermo/transport comments are included.
    """
    if not isinstance(species, Species):
        raise TypeError("species object must be an RMG Species")

    cantera_species = species.to_cantera(use_chemkin_identifier=True, all_species=all_species)
    species_data = cantera_species.input_data

    if verbose:
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

    if verbose:
        try:
            smiles = species.to_smiles()
            if smiles:
                species_data["note"] = smiles
        except Exception:
            pass
        if species.thermo and species.thermo.comment:
            clean_comment = species.thermo.comment.replace('\n', '; ').strip()
            if clean_comment:
                if "thermo" in species_data and isinstance(species_data["thermo"], dict):
                    species_data["thermo"]["note"] = clean_comment

    # returns composition, name, thermo, and transport, and note
    return species_data


class CanteraWriter1(object):
    """
    This class listens to a RMG subject
    and writes an YAML file with the current state of the RMG model,
    to a yaml subfolder.


    A new instance of the class can be appended to a subject as follows:

    rmg = ...
    listener = CanteraWriter1(outputDirectory)
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)

    """

    def __init__(self, output_directory="", config=None):
        super(CanteraWriter1, self).__init__()
        self.output_directory = output_directory
        self.config = config
        self.output_subdirectory = os.path.join(self.output_directory, "cantera1")
        make_output_subdirectory(output_directory, "cantera1")

    def update(self, rmg):
        if self.config is not None and not self.config.should_write(
                rmg.reaction_model.iteration_num, rmg.is_final_save):
            return
        verbose = self.config.verbose_comments if (self.config and self.config.verbose_comments is not None) else rmg.verbose_comments
        save_edge = self.config.save_edge if (self.config and self.config.save_edge is not None) else rmg.save_edge_species

        num_species = len(rmg.reaction_model.core.species)
        this_output_path = os.path.join(self.output_subdirectory,
                                        f"chem{num_species:04d}.yaml")
        latest_output_path = os.path.join(self.output_subdirectory, 'chem.yaml')

        logging.info(f"Saving current model core to Cantera file: {this_output_path}")

        solvent_data = None
        if rmg.solvent:
            solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)

        surface_site_density = None
        if rmg.reaction_model.surface_site_density:
            surface_site_density = rmg.reaction_model.surface_site_density.value_si

        core_elements = rmg.reaction_model.core.get_elements()

        write_cantera(
            rmg.reaction_model.core.species,
            rmg.reaction_model.core.reactions,
            elements_in_use=core_elements,
            surface_site_density=surface_site_density,
            solvent=rmg.solvent,
            solvent_data=solvent_data,
            path=this_output_path,
        )
        shutil.copy2(this_output_path, latest_output_path)

        if verbose:
            annotated_path = os.path.join(self.output_subdirectory, 'chem_annotated.yaml')
            logging.info(f"Saving annotated Cantera file: {annotated_path}")
            write_cantera(
                rmg.reaction_model.core.species,
                rmg.reaction_model.core.reactions,
                elements_in_use=core_elements,
                surface_site_density=surface_site_density,
                solvent=rmg.solvent,
                solvent_data=solvent_data,
                path=annotated_path,
                verbose=True,
            )

        if save_edge:
            from rmgpy.rmg.model import ReactionModel
            logging.info('Saving current model core and edge to Cantera file...')
            edge_species = rmg.reaction_model.core.species + rmg.reaction_model.edge.species
            edge_reactions = rmg.reaction_model.core.reactions + rmg.reaction_model.edge.reactions
            edge_elements = ReactionModel(species=edge_species, reactions=edge_reactions).get_elements()

            this_edge_path = os.path.join(self.output_subdirectory,
                                          f"chem_edge{num_species:04d}.yaml")
            latest_edge_path = os.path.join(self.output_subdirectory, 'chem_edge.yaml')

            write_cantera(
                edge_species,
                edge_reactions,
                elements_in_use=edge_elements,
                surface_site_density=surface_site_density,
                solvent=rmg.solvent,
                solvent_data=solvent_data,
                path=this_edge_path,
            )
            shutil.copy2(this_edge_path, latest_edge_path)

            if verbose:
                annotated_edge_path = os.path.join(self.output_subdirectory,
                                                    'chem_edge_annotated.yaml')
                logging.info(f"Saving annotated edge Cantera file: {annotated_edge_path}")
                write_cantera(
                    edge_species,
                    edge_reactions,
                    elements_in_use=edge_elements,
                    surface_site_density=surface_site_density,
                    solvent=rmg.solvent,
                    solvent_data=solvent_data,
                    path=annotated_edge_path,
                    verbose=True,
                )
