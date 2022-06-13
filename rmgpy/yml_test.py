
import os
import yaml

from rmgpy.chemkin import load_chemkin_file
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.kinetics.arrhenius import Arrhenius, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius
from rmgpy.kinetics.falloff import Troe, ThirdBody, Lindemann
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.data.solvation import SolventData
from rmgpy.kinetics.surface import StickingCoefficient
from rmgpy.util import make_output_subdirectory     #checked




#make final chem.yml file with newest model
def convert_chemkin_to_yml(chemkin_path, dictionary_path=None, output="chem.yml"):
    if dictionary_path:
        spcs, rxns = load_chemkin_file(chemkin_path, dictionary_path=dictionary_path)
    else:
        spcs, rxns = load_chemkin_file(chemkin_path)
    write_yml(spcs, rxns, path=output)




def write_yml(spcs, rxns, solvent=None, solvent_data=None, path="chem.yml"):
    result_dict = get_mech_dict(spcs, rxns, solvent=solvent, solvent_data=solvent_data)
    with open(path, 'w') as f:
        yaml.dump(result_dict, stream=f)



def get_radicals(spc):
    if spc.molecule[0].to_smiles() == "[O][O]":  # treat oxygen as stable to improve radical analysis
        return 0
    else:
        return spc.molecule[0].multiplicity-1



###will have to edit this section!!!
def get_mech_dict(spcs, rxns, solvent='solvent', solvent_data=None):
    names = [x.label for x in spcs]
    for i,name in enumerate(names): #fix duplicate names
        if names.count(name) > 1:
            names[i] += "-"+str(names.count(name))

    is_surface = False
    for spc in spcs:
        if spc.contains_surface_site():
            is_surface = True
            break
    if not is_surface:
        result_dict = dict()
        result_dict["Units"] = dict()
        result_dict["Phases"] = [dict()]
        #result_dict["Phases"][0]["name"] = "phase"  #why is this here if it is redefined below?
        result_dict["Phases"][0]["Species"] = [obj_to_dict(x, spcs, names=names) for x in spcs]
        result_dict["Reactions"] = [obj_to_dict(x, spcs, names=names) for x in rxns]
        if solvent_data:
            result_dict["Solvents"] = [obj_to_dict(solvent_data, spcs, names=names, label=solvent)]
        return result_dict
    else:
        result_dict = dict()
        result_dict["Units"] = dict()
        result_dict["Phases"] = [dict(),dict()]
        result_dict["Phases"][0]["name"] = "gas"
        result_dict["Phases"][1]["name"] = "surface"
        result_dict["Interfaces"] = [dict()]
        result_dict["Phases"][0]["Species"] = [obj_to_dict(x,spcs,names=names) for x in spcs if not x.contains_surface_site()]
        result_dict["Phases"][1]["Species"] = [obj_to_dict(x,spcs,names=names) for x in spcs if x.contains_surface_site()]
        result_dict["Reactions"] = [obj_to_dict(x, spcs,names=names) for x in rxns]
        if solvent_data:
            result_dict["Solvents"] = [obj_to_dict(solvent_data, spcs, names=names, label=solvent)]
        return result_dict


def obj_to_dict(obj, spcs, names=None, label="solvent"):
    result_dict = dict()
    if isinstance(obj, Species):
        result_dict["name"] = names[spcs.index(obj)]
        result_dict["type"] = "Species"
        if obj.contains_surface_site():
            result_dict["adjlist"] = obj.molecule[0].to_adjacency_list()
        result_dict["smiles"] = obj.molecule[0].to_smiles()
        result_dict["thermo"] = obj_to_dict(obj.thermo, spcs)
        result_dict["radicalelectrons"] = get_radicals(obj)
        if obj.liquid_volumetric_mass_transfer_coefficient_data:
            result_dict["liquidvolumetricmasstransfercoefficient"] = dict()
            result_dict["liquidvolumetricmasstransfercoefficient"]["type"] = "TemperatureDependentLiquidVolumetricMassTransferCoefficient"
            result_dict["liquidvolumetricmasstransfercoefficient"]["Ts"] = obj.liquid_volumetric_mass_transfer_coefficient_data.Ts
            result_dict["liquidvolumetricmasstransfercoefficient"]["kLAs"] = obj.liquid_volumetric_mass_transfer_coefficient_data.kLAs
        if obj.henry_law_constant_data:
            result_dict["henrylawconstant"] = dict()
            result_dict["henrylawconstant"]["type"] = "TemperatureDependentHenryLawConstant"
            result_dict["henrylawconstant"]["Ts"] = obj.henry_law_constant_data.Ts
            result_dict["henrylawconstant"]["kHs"] = obj.henry_law_constant_data.kHs
    elif isinstance(obj, NASA):
        result_dict["polys"] = [obj_to_dict(k, spcs) for k in obj.polynomials]
        result_dict["type"] = "NASA"
    elif isinstance(obj, NASAPolynomial):
        result_dict["type"] = "NASApolynomial"
        result_dict["coefs"] = obj.coeffs.tolist()
        result_dict["Tmax"] = obj.Tmax.value_si
        result_dict["Tmin"] = obj.Tmin.value_si
    elif isinstance(obj, Reaction):
        result_dict["reactants"] = [names[spcs.index(x)] for x in obj.reactants]
        result_dict["products"] = [names[spcs.index(x)] for x in obj.products]
        result_dict["kinetics"] = obj_to_dict(obj.kinetics, spcs, names)
        result_dict["type"] = "ElementaryReaction"
        result_dict["radicalchange"] = sum([get_radicals(x) for x in obj.products]) - \
                                       sum([get_radicals(x) for x in obj.reactants])
    elif isinstance(obj, Arrhenius):
        obj.change_t0(1.0)
        result_dict["type"] = "Arrhenius"
        result_dict["A"] = obj.A.value_si
        result_dict["Ea"] = obj.Ea.value_si
        result_dict["n"] = obj.n.value_si
    elif isinstance(obj, StickingCoefficient):
        obj.change_t0(1.0)
        result_dict["type"] = "StickingCoefficient"
        result_dict["A"] = obj.A.value_si
        result_dict["Ea"] = obj.Ea.value_si
        result_dict["n"] = obj.n.value_si
    elif isinstance(obj, PDepArrhenius):
        result_dict["type"] = "PdepArrhenius"
        result_dict["Ps"] = obj.pressures.value_si.tolist()
        result_dict["arrs"] = [obj_to_dict(x, spcs) for x in obj.arrhenius]
    elif isinstance(obj, MultiArrhenius):
        result_dict["type"] = "MultiArrhenius"
        result_dict["arrs"] = [obj_to_dict(x, spcs) for x in obj.arrhenius]
    elif isinstance(obj, MultiPDepArrhenius):
        result_dict["type"] = "MultiPdepArrhenius"
        result_dict["parrs"] = [obj_to_dict(x, spcs) for x in obj.arrhenius]
    elif isinstance(obj, ThirdBody):
        result_dict["type"] = "ThirdBody"
        result_dict["arr"] = obj_to_dict(obj.arrheniusLow, spcs)
        result_dict["efficiencies"] = {spcs[i].label: float(val)
                                       for i, val in enumerate(obj.get_effective_collider_efficiencies(spcs)) if val != 1}
    elif isinstance(obj, Lindemann):
        result_dict["type"] = "Lindemann"
        result_dict["arrhigh"] = obj_to_dict(obj.arrheniusHigh, spcs)
        result_dict["arrlow"] = obj_to_dict(obj.arrheniusLow, spcs)
        result_dict["efficiencies"] = {spcs[i].label: float(val)
                                       for i, val in enumerate(obj.get_effective_collider_efficiencies(spcs)) if val != 1}
    elif isinstance(obj, Troe):
        result_dict["type"] = "Troe"
        result_dict["arrhigh"] = obj_to_dict(obj.arrheniusHigh, spcs)
        result_dict["arrlow"] = obj_to_dict(obj.arrheniusLow, spcs)
        result_dict["efficiencies"] = {spcs[i].label: float(val)
                                       for i, val in enumerate(obj.get_effective_collider_efficiencies(spcs)) if val != 1}
        result_dict["a"] = obj.alpha
        result_dict["T1"] = obj.T1.value_si
        if obj.T2:
            result_dict["T2"] = obj.T2.value_si
        else:
            result_dict["T2"] = 0.0
        result_dict["T3"] = obj.T3.value_si
    elif isinstance(obj, Chebyshev):
        result_dict["type"] = "Chebyshev"
        result_dict["coefs"] = obj.coeffs.value_si.tolist()
        result_dict["Tmin"] = obj.Tmin.value_si
        result_dict["Tmax"] = obj.Tmax.value_si
        result_dict["Pmin"] = obj.Pmin.value_si
        result_dict["Pmax"] = obj.Pmax.value_si
    elif isinstance(obj, Wilhoit):
        result_dict["type"] = "Wilhoit"
        result_dict["coefs"] = [obj.a0, obj.a1, obj.a2, obj.a3]
        result_dict["Cp0"] = obj.Cp0.value_si
        result_dict["Cpinf"] = obj.CpInf.value_si
        result_dict["H0"] = obj.H0.value_si
        result_dict["S0"] = obj.S0.value_si
        result_dict["B"] = obj.B.value_si
    elif isinstance(obj, SolventData):
        result_dict["type"] = "Solvent"
        result_dict["name"] = label
        viscosity = dict()
        viscosity["type"] = "RiedelViscosity"
        viscosity["A"] = float(obj.A)
        viscosity["B"] = float(obj.B)
        viscosity["C"] = float(obj.C)
        viscosity["D"] = float(obj.D)
        viscosity["E"] = float(obj.E)
        result_dict["mu"] = viscosity
    elif obj is None:
        return None
    else:
        raise ValueError("Object of type {} does not have a defined conversion to "
                         "ReactionMechanismSimulator format".format(type(obj)))
    return result_dict






class YAMLWriter(object):
    """
    This class listens to a RMG subject
    and writes an YAML file with the current state of the RMG model,
    to a yaml subfolder.


    A new instance of the class can be appended to a subject as follows:

    rmg = ...
    listener = YAMLWriter(outputDirectory)
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)

    """
    def __init__(self, output_directory=''):
        super(YAMLWriter, self).__init__()
        self.output_directory = output_directory
        make_output_subdirectory(output_directory, 'yaml')

    def update(self, rmg):
        solvent_data = None
        if rmg.solvent:
            solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
        write_yml(rmg.reaction_model.core.species, rmg.reaction_model.core.reactions, solvent=rmg.solvent, solvent_data=solvent_data,
                  path=os.path.join(self.output_directory, 'yaml', 'chem{}.yaml').format(len(rmg.reaction_model.core.species)))