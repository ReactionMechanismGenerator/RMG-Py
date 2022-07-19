# script to export mechanism uncertainty
# Sevy Harris
# 11/10/2021
# Steps to run this:
# 1.    Export your mechanism as a cantera yaml (needs at least Cantera 2.5)
# 2.    Specify the input chemkin file             'chem'
# 3.    Specify the input species dictionary file  'spec'
# 4.    Specify the input cantera yaml             'mechanism_yaml'
# 5.    Specify the output uncertainty yaml        'output_path'
# 6.    Adjust the database settings in uncertainty.load_database to match your RMG run
# 7.    Make sure you have the right RMG-database branch


import json
import yaml
import cantera as ct
from rmgpy.kinetics import Arrhenius, ArrheniusEP
from rmgpy.tools.uncertainty import Uncertainty
import time


def tokens_equal(RA, RB, PA, PB):
    """
    Function to check if the reactants of A match the reactants of B
    and also if the products match
    assume tokens are stripped and sorted
    """
    if len(RA) != len(RB):
        return False
    if len(PA) != len(PB):
        return False
    for i in range(0, len(RA)):
        if RA[i] != RB[i]:
            return False
    for i in range(0, len(PA)):
        if PA[i] != PB[i]:
            return False
    return True


def get_rmg_rxn_from_yaml_rxn(rmg_rxns, ct_rxn):
    """
    Function to convert from the cantera yaml equation to an RMG reaction
    rmg_rxns is a list of RMG Reactions
    ct_rxn is a dictionary, but could probably be changed to be just the reaction string
    returns a single RMG reaction from the list or None if no match is found
    """
    for rmg_rxn in rmg_rxns:
        if str(rmg_rxn) == ct_rxn['equation']:
            return rmg_rxn

    ct_str = ct_rxn['equation'].replace(' ', '')
    ct_str = ct_str.replace('(+M)', '+M')
    ct_tokens = ct_str.split('<=>')
    if len(ct_tokens) < 2:
        ct_tokens = ct_str.split('=>')
    if len(ct_tokens) < 2:
        ct_tokens = ct_str.split('<=')
    ct_reactants = ct_tokens[0].split("+")
    ct_products = ct_tokens[1].split("+")
    ct_reactants.sort()
    ct_products.sort()
    if 'M' in ct_reactants and 'M' in ct_products:
        ct_reactants.remove('M')
        ct_products.remove('M')
    for rmg_rxn in rmg_rxns:
        rmg_str = str(rmg_rxn).replace(' ', '')
        rmg_tokens = rmg_str.split('<=>')
        if len(rmg_tokens) < 2:
            rmg_tokens = rmg_str.split('=>')
        if len(rmg_tokens) < 2:
            rmg_tokens = rmg_str.split('<=')
        rmg_reactants = rmg_tokens[0].split("+")
        rmg_products = rmg_tokens[1].split("+")
        rmg_reactants.sort()
        rmg_products.sort()

        if tokens_equal(rmg_products, ct_products, rmg_reactants, ct_reactants):
            return rmg_rxn


def unpack_sensitivity(long_desc):
    start_str = 'sensitivities = '
    if start_str not in long_desc:
        return []
    start_index = long_desc.find(start_str) + len(start_str)
    sensitivities_str = long_desc[start_index:].replace("'", '"')
    sensitivities_str = sensitivities_str.replace("nan", '"-9999999"')
    sensitivities_str = sensitivities_str.replace('name', 'training_rxn_name')
    return json.loads(sensitivities_str)


def get_kinetics_dict(data):
    kinetics_dict = {}
    if isinstance(data, Arrhenius):
        ct_data = data.to_cantera_kinetics()
        kinetics_dict['A'] = ct_data.pre_exponential_factor
        kinetics_dict['Ea'] = ct_data.activation_energy
        kinetics_dict['n'] = ct_data.temperature_exponent
    elif isinstance(data, ArrheniusEP):
        kinetics_dict['A'] = data.A.value_si
        kinetics_dict['E0'] = data.E0.value_si
        kinetics_dict['n'] = data.n.value_si
        kinetics_dict['alpha'] = data.alpha.value_si
    else:
        print('Reaction type not yet implement', type(data))
    return kinetics_dict


tic = time.time()

# specify input files
chem = '/home/moon/autoscience/uncertainty/partial_nheptane/nheptane_mech_main/chem_annotated.inp'
spec = '/home/moon/autoscience/uncertainty/partial_nheptane/nheptane_mech_main/species_dictionary.txt'
mechanism_yaml = '/home/moon/autoscience/uncertainty/partial_nheptane/nheptane_mech_main/nheptane.yaml'

# specify the output path
output_path = '/home/moon/autoscience/uncertainty/partial_nheptane/nheptane_uncertainty2.yaml'

# Load the same database used to generate the model
print('Loading database...')
uncertainty = Uncertainty(output_directory='./temp/uncertainty')
uncertainty.load_model(chem, spec)
uncertainty.load_database(
    thermo_libraries=['BurkeH2O2', 'CurranPentane', 'FFCM1(-)', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'CBS_QB3_1dHR'],
    kinetics_families='default',
    reaction_libraries=['CurranPentane', 'FFCM1(-)', 'combustion_core/version5'],
    kinetics_depositories=['training'],
)

# Get the different kinetic and thermo sources
uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties()

# Load the mechanism yaml as a cantera object and as a regular yaml
gas1 = ct.Solution(mechanism_yaml)
with open(mechanism_yaml, 'r') as file:
    mech = yaml.load(file, Loader=yaml.FullLoader)

# Load the template reaction maps so we can look up training reactions
auto_gen_families = {}
for family_name in uncertainty.database.kinetics.families.keys():
    if uncertainty.database.kinetics.families[family_name].auto_generated and family_name not in auto_gen_families.keys():
        auto_gen_families[family_name] = uncertainty.database.kinetics.families[family_name].rules.get_entries()
        auto_gen_families[f'{family_name}_labels'] = [entry.label for entry in uncertainty.database.kinetics.families[family_name].rules.get_entries()]
        auto_gen_families[f'{family_name}_rxn_map'] = uncertainty.database.kinetics.families[family_name].get_reaction_matches(
            thermo_database=uncertainty.database.thermo,
            remove_degeneracy=True,
            get_reverse=True,
            exact_matches_only=False,
            fix_labels=True)


# Fill out the kinetics part of the uncertainty yaml
print('Computing kinetics uncertainty...')
for i, reaction in enumerate(mech['reactions']):
    rmg_rxn = get_rmg_rxn_from_yaml_rxn(uncertainty.reaction_list, reaction)
    rxn_uncertainty = {}
    rxn_uncertainty['value'] = float(uncertainty.kinetic_input_uncertainties[uncertainty.reaction_list.index(rmg_rxn)])
    # Note the source- assume only one source per reaction, or else this might overwrite existing reaction
    for src_key in uncertainty.reaction_sources_dict[rmg_rxn].keys():
        rxn_uncertainty['source'] = src_key
    if src_key == 'Rate Rules':
        family_name = str(rmg_rxn.family)
        rxn_uncertainty['family'] = family_name
        src = uncertainty.reaction_sources_dict[rmg_rxn]['Rate Rules'][1]
        new_trees = (src['node'] != '')
        if new_trees:
            rxn_uncertainty['rule_node'] = src['node']
            rule_index = auto_gen_families[family_name + '_labels'].index(src['node'])
            rule = auto_gen_families[family_name][rule_index]
            derivatives = unpack_sensitivity(rule.long_desc)

            # check that the sensitivities correspond to the correct training reactions
            if len(auto_gen_families[family_name + '_rxn_map'][src['node']]) != len(derivatives):
                raise ValueError(f"There is a mismatch between the number of sensitivities and the number of training reactions for {src['node']} from {family_name}")

            rxn_uncertainty['sensitivities'] = derivatives
            for j in range(0, len(derivatives)):
                t_rxn = auto_gen_families[family_name + '_rxn_map'][src['node']][j]
                ct_t_rxn = t_rxn.to_cantera()
                derivatives[j]['training_name'] = ct_t_rxn.equation
                derivatives[j]['training_rate'] = get_kinetics_dict(t_rxn.kinetics)

        # handle old trees
        else:
            uncertainties = []
            info = uncertainty.reaction_sources_dict[rmg_rxn]['Rate Rules']

            family_name = info[0]
            rxn_uncertainty['family'] = family_name
            rxn_uncertainty['exact'] = info[1]['exact']

            for rule_entry, train_entry, weight in info[1]['training']:
                rate = train_entry.data.to_cantera_kinetics()
                uncertainties.append({'training_name': train_entry.label,
                                      'training_rate': get_kinetics_dict(train_entry.data),
                                      'weight': weight,
                                      'training_desc': train_entry.long_desc})
            for rule_entry, weight in info[1]['rules']:
                uncertainties.append({'rule_name': rule_entry.label,
                                      'rule_rate': get_kinetics_dict(rule_entry.data),
                                      'weight': weight,
                                      'rule_desc': rule_entry.long_desc})
            rxn_uncertainty['train_rxn_uncertainty'] = uncertainties

    mech['reactions'][i]['uncertainty'] = rxn_uncertainty


# Fill out the thermo part of the uncertainty yaml
print('Computing thermo uncertainty...')
for i, species in enumerate(mech['species']):
    species_uncertainty = {}
    uncertainty_entry = uncertainty.species_sources_dict[uncertainty.species_list[i]]
    if 'Library' in uncertainty_entry.keys():
        species_uncertainty['source'] = 'Library'
    elif 'GAV' in uncertainty_entry.keys():
        group_entry_uncertainties = []
        species_uncertainty['source'] = 'GAV'
        for key in uncertainty_entry['GAV'].keys():
            for entry, freq in uncertainty_entry['GAV'][key]:
                group_entry_uncertainties.append({'family_name': key,
                                                  'entry_name': entry.label,
                                                  'entry_index': entry.index,
                                                  'entry_freq': freq})
        species_uncertainty['group_entries'] = group_entry_uncertainties

    else:
        print('Missing species uncertainty source')

    species_uncertainty['value'] = float(uncertainty.thermo_input_uncertainties[i])

    mech['species'][i]['uncertainty'] = species_uncertainty

# save the result
print(f'Saving uncertainty yaml...')
with open(output_path, 'w') as file:
    output = yaml.dump(mech, file)

toc = time.time()
print(f'Completed in {toc-tic} seconds')
