import os
import yaml

from rmgpy.species import Species
from rmgpy.kinetics.arrhenius import (
    MultiArrhenius,
    MultiPDepArrhenius,
)
from rmgpy.util import make_output_subdirectory
from datetime import datetime
from rmgpy.chemkin import get_species_identifier


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
        block1, block2, block3, block4 = write_surface_species(
            spcs, rxns, surface_site_density
        )
    else:
        # get_mech_dict writes yaml files without creating separate
        result_dict = get_mech_dict_nonsurface(
            spcs, rxns, solvent=solvent, solvent_data=solvent_data
        )
        block1, block2, block3, block4 = write_nonsurface_species(spcs)

    with open(path, "w") as f:

        # generator line
        f.write("generator: RMG\n")

        # datetime object containing current date and time
        now = datetime.now()
        dt_string = now.strftime("%a, %d %b %Y %H:%M:%S")
        f.write(f"date: {dt_string}\n")

        # units line
        f.write(
            "\nunits: {length: cm, time: s, quantity: mol, activation-energy: kcal/mol}\n"
        )
        f.write("\n")

        #'phases' line (below)

        f.write(block1)
        f.write(block2)
        f.write(block3)
        f.write(block4)
        yaml.dump(result_dict, stream=f, sort_keys=False)


def write_nonsurface_species(spcs):
    """
    Yaml files without surface species begin with the following blocks of text.
    """

    #'phases' line (below)
    block1 = """
phases:
- name: gas
  thermo: ideal-gas
  elements: [H, D, T, C, Ci, O, Oi, N, Ne, Ar, He, Si, S, F, Cl, Br, I, X]"""

    #'species' section in phases section

    sorted_species = sorted(spcs, key=lambda spcs: spcs.index)
    species_to_write = [get_species_identifier(spec) for spec in sorted_species]

    block2 = f"""
  species: [{', '.join(species_to_write)}]
  kinetics: gas"""

    block3 = """
  transport: mixture-averaged
  state: {T: 300.0, P: 1 atm}
        """

    block4 = """
elements:
- symbol: Ci
  atomic-weight: 13.003
- symbol: D
  atomic-weight: 2.014
- symbol: Oi
  atomic-weight: 17.999
- symbol: T
  atomic-weight: 3.016
- symbol: X
  atomic-weight: 195.083

"""

    return block1, block2, block3, block4


def write_surface_species(spcs, rxns, surface_site_density):
    """
    Yaml files with surface species begin with the following blocks of text, 
    which includes TWO phases instead of just one.
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
        get_species_identifier(surface_species)
        for surface_species in sorted_surface_species
    ]

    sorted_gas_species = sorted(gas_species, key=lambda gas_species: gas_species.index)
    gas_species_to_write = [
        get_species_identifier(gas_species) for gas_species in sorted_gas_species
    ]

    # gas part
    block1 = f"""
phases:
- name: gas
  thermo: ideal-gas
  elements: [H, D, T, C, Ci, O, Oi, N, Ne, Ar, He, Si, S, F, Cl, Br, I, X]
  species: [{', '.join(gas_species_to_write)}]
  kinetics: gas
  reactions: [gas_reactions]"""

    block2 = """
  transport: mixture-averaged
  state: {T: 300.0, P: 1 atm}"""

    block3 = f""" 
- name: {surface_species[0].smiles.replace("[","").replace("]","")}_surface
  thermo: ideal-surface
  adjacent-phases: [gas]
  elements: [H, D, T, C, Ci, O, Oi, N, Ne, Ar, He, Si, S, F, Cl, Br, I, X]
  species: [{', '.join(surface_species_to_write)}]
  kinetics: surface
  reactions: [surface_reactions]     
  site-density: {surface_site_density * 1e-4 }
        """

    # surface_site_density * 1e-4 #in units of mol/cm^2

    block4 = """
elements:
- symbol: Ci
  atomic-weight: 13.003
- symbol: D
  atomic-weight: 2.014
- symbol: Oi
  atomic-weight: 17.999
- symbol: T
  atomic-weight: 3.016
- symbol: X
  atomic-weight: 195.083

"""
    return block1, block2, block3, block4

def get_radicals(spc):
    if (
        spc.molecule[0].to_smiles() == "[O][O]"
    ):  # treat oxygen as stable to improve radical analysis
        return 0
    else:
        return spc.molecule[0].multiplicity - 1


def get_mech_dict_surface(spcs, rxns, solvent="solvent", solvent_data=None):
    """
    For systems with surface species/reactions. 
    Adds 'species', 'gas-reactions', and 'surface-reactions' to result_dict.
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
    result_dict["species"] = [obj_to_dict(x, spcs, names=names) for x in spcs]

    # separate gas and surface reactions

    gas_reactions = []
    for rmg_rxn in gas_rxns:
        gas_reactions.extend(reaction_to_dicts(rmg_rxn, spcs))
    result_dict["gas_reactions"] = gas_reactions

    surface_reactions = []
    for rmg_rxn in surface_rxns:
        surface_reactions.extend(reaction_to_dicts(rmg_rxn, spcs))
    result_dict["surface_reactions"] = surface_reactions

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
    result_dict["species"] = [obj_to_dict(x, spcs, names=names) for x in spcs]

    reactions = []
    for rmg_rxn in rxns:
        reactions.extend(reaction_to_dicts(rmg_rxn, spcs))
    result_dict["reactions"] = reactions

    return result_dict


def reaction_to_dicts(obj, spcs):
    """
    Takes an RMG reaction object (obj), returns a list of dictionaries
    for YAML properties. For most reaction objects the list will be of
    length 1, but a MultiArrhenius or MultiPDepArrhenius will be longer.
    """

    reaction_list = []
    if isinstance(obj.kinetics, MultiArrhenius) or isinstance(
        obj.kinetics, MultiPDepArrhenius
    ):
        list_of_cantera_reactions = obj.to_cantera(use_chemkin_identifier=True)
    else:
        list_of_cantera_reactions = [obj.to_cantera(use_chemkin_identifier=True)]

    for reaction in list_of_cantera_reactions:
        reaction_data = reaction.input_data
        efficiencies = getattr(obj.kinetics, "efficiencies", {})
        if efficiencies:
            reaction_data["efficiencies"] = {
                spcs[i].to_chemkin(): float(val)
                for i, val in enumerate(
                    obj.kinetics.get_effective_collider_efficiencies(spcs)
                )
                if val != 1
            }
        reaction_list.append(reaction_data)

    return reaction_list


def obj_to_dict(obj, spcs, names=None, label="solvent"):
    """
    Takes an RMG species object (obj), returns a list of dictionaries
    for YAML properties. Also adds in the number of surface sites 
    ('sites') to dictionary. 
    """

    result_dict = dict()

    if isinstance(obj, Species):
        s = obj.to_cantera(use_chemkin_identifier=True)
        species_data = s.input_data
        try:
            result_dict["note"] = obj.transport_data.comment
        except:
            pass
        if "size" in species_data:
            sites = species_data["size"]
            species_data.pop("size", None)
            species_data["sites"] = sites
        species_data.update(result_dict)
        return (
            species_data  # returns composition, name, thermo, and transport, and note
        )


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
        make_output_subdirectory(output_directory, "cantera")

    def update(self, rmg):
        solvent_data = None
        if rmg.solvent:
            solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
        surface_site_density = None
        NoneType = type(None)  # types.py no longer has NoneType, says slack.overview
        if not type(rmg.reaction_model.surface_site_density) == NoneType:
            surface_site_density = rmg.reaction_model.surface_site_density.value_si
        write_cantera(
            rmg.reaction_model.core.species,
            rmg.reaction_model.core.reactions,
            surface_site_density=surface_site_density,
            solvent=rmg.solvent,
            solvent_data=solvent_data,
            path=os.path.join(self.output_directory, "cantera", "chem{}.yaml").format(
                len(rmg.reaction_model.core.species)
            ),
        )
