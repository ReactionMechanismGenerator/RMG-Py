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
This module contains functionality for parsing Arkane input files.
"""

import logging
import os.path

import numpy as np

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.rmg import get_db
from rmgpy.exceptions import InputError, DatabaseError
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.kinetics.model import PDepKineticsModel, TunnelingModel
from rmgpy.kinetics.tunneling import Wigner, Eckart
from rmgpy.molecule import Molecule
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.pdep.configuration import Configuration
from rmgpy.pdep.network import Network
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species, TransitionState
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.torsion import HinderedRotor, FreeRotor
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.thermodata import ThermoData
from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.kinetics.uncertainties import RateUncertainty
from rmgpy.transport import TransportData
from rmgpy.data.solvation import SoluteTSData
from rmgpy.util import as_list

from arkane.common import is_pdep
from arkane.encorr.ae import AEJob
from arkane.encorr.bac import BACJob
from arkane.encorr.corr import assign_frequency_scale_factor
from arkane.explorer import ExplorerJob
from arkane.kinetics import KineticsJob
from arkane.modelchem import LOT, LevelOfTheory, CompositeLevelOfTheory, model_chem_to_lot
from arkane.pdep import PressureDependenceJob
from arkane.statmech import StatMechJob
from arkane.thermo import ThermoJob

################################################################################


species_dict, transition_state_dict, reaction_dict, network_dict = dict(), dict(), dict(), dict()
job_list = list()


def database(thermoLibraries=None, transportLibraries=None, reactionLibraries=None, frequenciesLibraries=None,
             kineticsFamilies='default', kineticsDepositories='default', kineticsEstimator='rate rules'):
    """Load the RMG database"""
    thermo_libraries = as_list(thermoLibraries, default=[])
    transport_libraries = as_list(transportLibraries, default=None)
    reaction_libraries = as_list(reactionLibraries, default=[])

    database_directory = settings['database.directory']

    if kineticsDepositories == 'default':
        kinetics_depositories = ['training']
    elif kineticsDepositories == 'all':
        kinetics_depositories = None
    else:
        if not isinstance(kineticsDepositories, list):
            raise InputError(
                "kinetics_depositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe'].")
        kinetics_depositories = kineticsDepositories

    if kineticsFamilies in ('default', 'all', 'none'):
        kinetics_families = kineticsFamilies
    else:
        if not isinstance(kineticsFamilies, list):
            raise InputError(
                "kineticsFamilies should be either 'default', 'all', 'none', or a list of names eg. "
                "['H_Abstraction','R_Recombination'] or ['!Intra_Disproportionation'].")
        kinetics_families = kineticsFamilies

    rmg_database = get_db() or RMGDatabase()

    rmg_database.load(
        path=database_directory,
        thermo_libraries=thermo_libraries,
        transport_libraries=transport_libraries,
        reaction_libraries=reaction_libraries,
        seed_mechanisms=[],
        kinetics_families=kinetics_families,
        kinetics_depositories=kinetics_depositories,
        depository=False,  # Don't bother loading the depository information, as we don't use it
    )

    for family in rmg_database.kinetics.families.values():  # load training
        if not family.auto_generated:
            family.add_rules_from_training(thermo_database=rmg_database.thermo)
            family.fill_rules_by_averaging_up(verbose=True)


def species(label, *args, **kwargs):
    """Load a species from an input file"""
    global species_dict, job_list
    if label in species_dict:
        raise ValueError('Multiple occurrences of species with label {0!r}.'.format(label))
    logging.info('Loading species {0}...'.format(label))

    spec = Species(label=label)
    species_dict[label] = spec

    path = None
    if len(args) == 1:
        # The argument is a path to a conformer input file
        path = args[0]
        job = StatMechJob(species=spec, path=path)
        logging.debug('Added species {0} to a stat mech job.'.format(label))
        job_list.append(job)
    elif len(args) > 1:
        raise InputError('species {0} can only have two non-keyword argument '
                         'which should be the species label and the '
                         'path to a quantum file.'.format(spec.label))

    if len(kwargs) > 0:
        # The species parameters are given explicitly
        structure = None
        E0 = None
        modes = []
        spin_multiplicity = 0
        optical_isomers = 1
        molecular_weight = None
        collision_model = None
        energy_transfer_model = None
        thermo = None
        reactive = True
        for key, value in kwargs.items():
            if key == 'structure':
                structure = value
            elif key == 'E0':
                E0 = value
            elif key == 'modes':
                modes = value
            elif key == 'spinMultiplicity':
                spin_multiplicity = value
            elif key == 'opticalIsomers':
                optical_isomers = value
            elif key == 'molecularWeight':
                molecular_weight = value
            elif key == 'collisionModel':
                collision_model = value
            elif key == 'energyTransferModel':
                energy_transfer_model = value
            elif key == 'thermo':
                thermo = value
            elif key == 'reactive':
                reactive = value
            else:
                raise TypeError('species() got an unexpected keyword argument {0!r}.'.format(key))

        if structure:
            spec.molecule = [structure]
        spec.conformer = Conformer(E0=E0, modes=modes, spin_multiplicity=spin_multiplicity,
                                   optical_isomers=optical_isomers)
        if molecular_weight is not None:
            spec.molecular_weight = molecular_weight
        elif spec.molecular_weight is None and is_pdep(job_list):
            # If a structure was given, simply calling spec.molecular_weight will calculate the molecular weight
            # If one of the jobs is pdep and no molecular weight is given or calculated, raise an error
            raise ValueError("No molecularWeight was entered for species {0}. Since a structure wasn't given"
                             " as well, the molecularWeight, which is important for pressure dependent jobs,"
                             " cannot be reconstructed.".format(spec.label))
        spec.transport_data = collision_model
        spec.energy_transfer_model = energy_transfer_model
        spec.thermo = thermo
        spec.reactive = reactive

        if spec.reactive and path is None and spec.thermo is None and spec.conformer.E0 is None:
            if not spec.molecule:
                raise InputError('Neither thermo, E0, species file path, nor structure specified, cannot estimate'
                                 ' thermo properties of species {0}'.format(spec.label))
            try:
                db = get_db('thermo')
                if db is None:
                    raise DatabaseError('Thermo database is None.')
            except DatabaseError:
                logging.warning("The database isn't loaded, cannot estimate thermo for {0}. "
                                "If it is a bath gas, set reactive = False to avoid generating "
                                "thermo.".format(spec.label))
            else:
                logging.info('No E0 or thermo found, estimating thermo and E0 of species {0} using'
                             ' RMG-Database...'.format(spec.label))
                spec.thermo = db.get_thermo_data(spec)
                if spec.thermo.E0 is None:
                    th = spec.thermo.to_wilhoit()
                    spec.conformer.E0 = th.E0
                    spec.thermo.E0 = th.E0
                else:
                    spec.conformer.E0 = spec.thermo.E0

        if spec.reactive and spec.thermo and not spec.has_statmech() and structure is not None:
            # generate stat mech info if it wasn't provided before
            spec.generate_statmech()

        if not energy_transfer_model:
            # default to RMG's method of generating energy_transfer_model
            spec.generate_energy_transfer_model()

    return spec


def transitionState(label, *args, **kwargs):
    """Load a transition state from an input file"""
    global transition_state_dict
    if label in transition_state_dict:
        raise ValueError('Multiple occurrences of transition state with label {0!r}.'.format(label))
    logging.info('Loading transition state {0}...'.format(label))
    ts = TransitionState(label=label)
    transition_state_dict[label] = ts

    if len(args) == 1 and len(kwargs) == 0:
        # The argument is a path to a conformer input file
        path = args[0]
        job = StatMechJob(species=ts, path=path)
        job_list.append(job)

    elif len(args) == 0:
        # The species parameters are given explicitly
        E0 = None
        modes = []
        spin_multiplicity = 1
        optical_isomers = 1
        frequency = None
        for key, value in kwargs.items():
            if key == 'E0':
                E0 = value
            elif key == 'modes':
                modes = value
            elif key == 'spinMultiplicity':
                spin_multiplicity = value
            elif key == 'opticalIsomers':
                optical_isomers = value
            elif key == 'frequency':
                frequency = value
            else:
                raise TypeError('transition_state() got an unexpected keyword argument {0!r}.'.format(key))

        ts.conformer = Conformer(E0=E0, modes=modes, spin_multiplicity=spin_multiplicity,
                                 optical_isomers=optical_isomers)
        ts.frequency = frequency
    else:
        if len(args) == 0 and len(kwargs) == 0:
            raise InputError(
                'The transition_state needs to reference a quantum job file or contain kinetic information.')
        raise InputError('The transition_state can only link a quantum job or directly input information, not both.')

    return ts


def reaction(label, reactants, products, transitionState=None, kinetics=None, tunneling=''):
    """Load a reaction from an input file"""
    global reaction_dict, species_dict, transition_state_dict
    if label in reaction_dict:
        label = label + transitionState
        if label in reaction_dict:
            raise ValueError('Multiple occurrences of reaction with label {0!r}.'.format(label))
    logging.info('Loading reaction {0}...'.format(label))
    reactants = sorted([species_dict[spec] for spec in reactants])
    products = sorted([species_dict[spec] for spec in products])
    if transitionState:
        transitionState = transition_state_dict[transitionState]
    if transitionState and (tunneling == '' or tunneling is None):
        transitionState.tunneling = None
    elif tunneling.lower() == 'wigner':
        transitionState.tunneling = Wigner(frequency=None)
    elif tunneling.lower() == 'eckart':
        transitionState.tunneling = Eckart(frequency=None, E0_reac=None, E0_TS=None, E0_prod=None)

    elif transitionState and not isinstance(tunneling, TunnelingModel):
        raise ValueError('Unknown tunneling model {0!r}.'.format(tunneling))
    rxn = Reaction(label=label, reactants=reactants, products=products, transition_state=transitionState,
                   kinetics=kinetics)

    if rxn.transition_state is None and rxn.kinetics is None:
        logging.info('estimating rate of reaction {0} using RMG-database')
        if not all([m.molecule != [] for m in rxn.reactants + rxn.products]):
            raise ValueError('chemical structures of reactants and products not available for RMG estimation of '
                             'reaction {0}'.format(label))
        db = get_db('kinetics')
        rxns = db.generate_reactions_from_libraries(reactants=rxn.reactants, products=rxn.products)
        rxns = [r for r in rxns if r.elementary_high_p]

        if rxns:
            for r in rxns:
                if isinstance(rxn.kinetics, PDepKineticsModel):
                    boo = rxn.generate_high_p_limit_kinetics()
                if boo:
                    rxn = r
                    break

        if rxns == [] or not boo:
            logging.info('No library reactions tagged with elementary_high_p found for reaction {0}, generating '
                         'reactions from RMG families'.format(label))
            rxn = list(db.generate_reactions_from_families(reactants=rxn.reactants, products=rxn.products))
            model = CoreEdgeReactionModel()
            model.verbose_comments = True
            for r in rxn:
                model.apply_kinetics_to_reaction(r)

    if isinstance(rxn, Reaction):
        reaction_dict[label] = rxn
    else:
        for i in range(len(rxn)):
            reaction_dict[label + str(i)] = rxn[i]

    return rxn


def network(label, isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None):
    """Load a network from an input file"""
    global network_dict, species_dict, reaction_dict
    logging.info('Loading network {0}...'.format(label))
    isomers0 = isomers or []
    isomers = []
    for isomer in isomers0:
        if isinstance(isomer, (list, tuple)):
            raise ValueError('Only one species can be present in a unimolecular isomer.')
        isomers.append(species_dict[isomer])

    reactants0 = reactants or []
    reactants = []
    for reactant in reactants0:
        if not isinstance(reactant, (list, tuple)):
            reactant = [reactant]
        reactants.append(sorted([species_dict[spec] for spec in reactant]))

    if pathReactions is None:
        # Only add reactions that match reactants and/or isomers
        path_reactions = []
        for rxn in reaction_dict.values():
            if not rxn.is_unimolecular():
                # this reaction is not pressure dependent
                continue
            reactant_is_isomer = len(rxn.reactants) == 1 and rxn.reactants[0] in isomers
            product_is_isomer = len(rxn.products) == 1 and rxn.products[0] in isomers
            reactant_is_reactant = any(
                [frozenset(rxn.reactants) == frozenset(reactant_pair) for reactant_pair in reactants])
            product_is_reactant = any(
                [frozenset(rxn.products) == frozenset(reactant_pair) for reactant_pair in reactants])
            if reactant_is_isomer or reactant_is_reactant or product_is_isomer or product_is_reactant:
                path_reactions.append(rxn)
        logging.debug('Path reactions {} were found for network {}'.format(
                       [rxn.label for rxn in path_reactions], label))
    else:
        path_reactions_0 = pathReactions
        path_reactions = []
        for rxn in path_reactions_0:
            path_reactions.append(reaction_dict[rxn])

    if products is None:
        # Figure out which configurations are isomers, reactant channels, and product channels
        products = []
        for rxn in path_reactions:
            # Sort bimolecular configurations so that we always encounter them in the same order
            # The actual order doesn't matter, as long as it is consistent
            rxn.reactants.sort()
            rxn.products.sort()
            # All reactant configurations not already defined as reactants or 
            # isomers are assumed to be product channels
            if len(rxn.reactants) == 1 and rxn.reactants[0] not in isomers and rxn.reactants not in products:
                products.append(rxn.reactants)
            elif len(rxn.reactants) > 1 and rxn.reactants not in reactants and rxn.reactants not in products:
                products.append(rxn.reactants)
            # All product configurations not already defined as reactants or 
            # isomers are assumed to be product channels
            if len(rxn.products) == 1 and rxn.products[0] not in isomers and rxn.products not in products:
                products.append(rxn.products)
            elif len(rxn.products) > 1 and rxn.products not in reactants and rxn.products not in products:
                products.append(rxn.products)
    else:
        products0 = products or []
        products = []
        for product in products0:
            if not isinstance(product, (list, tuple)):
                product = [product]
            products.append(sorted([species_dict[spec] for spec in product]))

    isomers = [Configuration(species) for species in isomers]
    reactants = [Configuration(*species) for species in reactants]
    products = [Configuration(*species) for species in products]

    bath_gas_0 = bathGas or {}
    bath_gas = {}
    for spec, fraction in bath_gas_0.items():
        bath_gas[species_dict[spec]] = fraction

    network = Network(
        label=label,
        isomers=isomers,
        reactants=reactants,
        products=products,
        path_reactions=path_reactions,
        bath_gas=bath_gas,
    )
    network_dict[label] = network


def kinetics(label, Tmin=None, Tmax=None, Tlist=None, Tcount=0, sensitivity_conditions=None, three_params=True):
    """Generate a kinetics job"""
    global job_list, reaction_dict
    try:
        rxn = reaction_dict[label]
    except KeyError:
        raise ValueError('Unknown reaction label {0!r} for kinetics() job.'.format(label))
    job = KineticsJob(reaction=rxn, Tmin=Tmin, Tmax=Tmax, Tcount=Tcount, Tlist=Tlist,
                      sensitivity_conditions=sensitivity_conditions, three_params=three_params)
    job_list.append(job)


def statmech(label):
    """Generate a statmech job"""
    global job_list, species_dict, transition_state_dict
    if label in species_dict or label in transition_state_dict:
        for job in job_list:
            if job.species.label == label:
                break
        else:
            raise ValueError('Could not create StatMechJob for {0!r}; no path specified.'.format(label))
    else:
        raise ValueError('Unknown species or transition state label {0!r} for statmech() job.'.format(label))


def thermo(label, thermoClass):
    """Generate a thermo job"""
    global job_list, species_dict
    try:
        spec = species_dict[label]
    except KeyError:
        raise ValueError('Unknown species label {0!r} for thermo() job.'.format(label))
    job = ThermoJob(species=spec, thermo_class=thermoClass)
    job_list.append(job)


def pressureDependence(label, Tmin=None, Tmax=None, Tcount=0, Tlist=None, Pmin=None, Pmax=None, Pcount=0, Plist=None,
                       maximumGrainSize=None, minimumGrainCount=0, method=None, interpolationModel=None,
                       activeKRotor=True, activeJRotor=True, rmgmode=False, sensitivity_conditions=None, sensitivity_perturbation=(2.0,'kcal/mol')):
    """Generate a pressure dependent job"""
    global job_list, network_dict

    if isinstance(interpolationModel, str):
        interpolationModel = (interpolationModel,)

    nwk = None
    if label in list(network_dict.keys()):
        nwk = network_dict[label]

    job = PressureDependenceJob(network=nwk, Tmin=Tmin, Tmax=Tmax, Tcount=Tcount, Tlist=Tlist,
                                Pmin=Pmin, Pmax=Pmax, Pcount=Pcount, Plist=Plist,
                                maximumGrainSize=maximumGrainSize, minimumGrainCount=minimumGrainCount,
                                method=method, interpolationModel=interpolationModel,
                                activeKRotor=activeKRotor, activeJRotor=activeJRotor,
                                rmgmode=rmgmode, sensitivity_conditions=sensitivity_conditions,
                                sensitivity_perturbation=sensitivity_perturbation)
    job_list.append(job)


def explorer(source, explore_tol=0.01, energy_tol=np.inf, flux_tol=0.0, bathGas=None, maximumRadicalElectrons=np.inf):
    """Generate an explorer job"""
    global job_list, species_dict
    for job in job_list:
        if isinstance(job, PressureDependenceJob):
            pdepjob = job
            break
    else:
        raise InputError('the explorer block must occur after the pressureDependence block')

    source = [species_dict[name] for name in source]

    bath_gas = {species_dict[spec]: fraction for spec, fraction in bathGas.items()} if bathGas else None

    job = ExplorerJob(source=source, pdepjob=pdepjob, explore_tol=explore_tol,
                      energy_tol=energy_tol, flux_tol=flux_tol, bath_gas=bath_gas,
                      maximum_radical_electrons=maximumRadicalElectrons)
    job_list.append(job)


def ae(species_energies, level_of_theory=None, write_to_database=False, overwrite=False):
    """Generate an atom energy job"""
    global job_list
    job = AEJob(
        species_energies,
        level_of_theory=level_of_theory,
        write_to_database=write_to_database,
        overwrite=overwrite
    )
    job_list.append(job)


def bac(level_of_theory, bac_type='p', train_names='main', crossval_n_folds=1,
        idxs=None, exclude_idxs=None,
        exclude_elements=None, charge='all', multiplicity='all',
        weighted=False, write_to_database=False, overwrite=False,
        fit_mol_corr=True, global_opt=True, global_opt_iter=10):
    """Generate a BAC job"""
    global job_list
    job = BACJob(
        level_of_theory,
        bac_type=bac_type,
        db_names=train_names,
        crossval_n_folds=crossval_n_folds,
        idxs=idxs,
        exclude_idxs=exclude_idxs,
        exclude_elements=exclude_elements,
        charge=charge,
        multiplicity=multiplicity,
        weighted=weighted,
        write_to_database=write_to_database,
        overwrite=overwrite,
        fit_mol_corr=fit_mol_corr,
        global_opt=global_opt,
        global_opt_iter=global_opt_iter)
    job_list.append(job)


def SMILES(smiles):
    """Make a Molecule object from SMILES"""
    return Molecule().from_smiles(smiles)


def adjacencyList(adj):
    """Make a Molecule object from an adjacency list"""
    return Molecule().from_adjacency_list(adj,
                                          raise_atomtype_exception=False,
                                          raise_charge_exception=False)


def InChI(inchi):
    """Make a Molecule object from InChI"""
    return Molecule().from_inchi(inchi)


def load_necessary_databases():
    """
    loads transport and statmech databases
    """
    from rmgpy.data.statmech import StatmechDatabase
    from rmgpy.data.transport import TransportDatabase

    # only load if they are not there already.
    try:
        get_db('transport')
        get_db('statmech')
    except DatabaseError:
        logging.info("Databases not found. Making databases")
        db = RMGDatabase()
        db.statmech = StatmechDatabase()
        db.statmech.load(os.path.join(settings['database.directory'], 'statmech'))

        db.transport = TransportDatabase()
        db.transport.load(os.path.join(settings['database.directory'], 'transport'))


################################################################################


def load_input_file(path):
    """
    Load the Arkane input file located at `path` on disk, and return a list of
    the jobs defined in that file.
    """
    global species_dict, transition_state_dict, reaction_dict, network_dict, job_list

    # Clear module-level variables
    species_dict, transition_state_dict, reaction_dict, network_dict = dict(), dict(), dict(), dict()
    job_list = []

    global_context = {'__builtins__': None}
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'range': range,
        # Collision
        'TransportData': TransportData,
        'SingleExponentialDown': SingleExponentialDown,
        # Kinetics
        'Arrhenius': Arrhenius,
        'RateUncertainty': RateUncertainty,
        # Statistical mechanics
        'IdealGasTranslation': IdealGasTranslation,
        'LinearRotor': LinearRotor,
        'NonlinearRotor': NonlinearRotor,
        'KRotor': KRotor,
        'SphericalTopRotor': SphericalTopRotor,
        'HarmonicOscillator': HarmonicOscillator,
        'HinderedRotor': HinderedRotor,
        'FreeRotor': FreeRotor,
        # Thermo
        'ThermoData': ThermoData,
        'Wilhoit': Wilhoit,
        'NASA': NASA,
        'NASAPolynomial': NASAPolynomial,
        # Functions
        'reaction': reaction,
        'species': species,
        'transitionState': transitionState,
        'network': network,
        'database': database,
        # Jobs
        'kinetics': kinetics,
        'SoluteTSData': SoluteTSData,
        'statmech': statmech,
        'thermo': thermo,
        'pressureDependence': pressureDependence,
        'explorer': explorer,
        'bac': bac,
        'ae': ae,
        # Miscellaneous
        'SMILES': SMILES,
        'adjacencyList': adjacencyList,
        'InChI': InChI,
        'LevelOfTheory': LevelOfTheory,
        'CompositeLevelOfTheory': CompositeLevelOfTheory,
    }

    load_necessary_databases()

    with open(path, 'r') as f:
        content = f.read()
        try:
            exec(content, global_context, local_context)
        except (NameError, TypeError, SyntaxError) as e:
            logging.error('The input file {0!r} was invalid:'.format(path))
            line_number = e.__traceback__.tb_next.tb_lineno
            logging.error(f'Error occurred at or near line {line_number} of {path}.')
            lines = content.splitlines()
            logging.error(f'Line: {lines[line_number - 1]}')
            raise

    model_chemistry = local_context.get('modelChemistry', None)
    level_of_theory = process_model_chemistry(model_chemistry)
    if isinstance(model_chemistry, LOT):
        model_chemistry = model_chemistry.to_model_chem()

    author = local_context.get('author', '')
    if 'frequencyScaleFactor' in local_context:
        frequency_scale_factor = local_context.get('frequencyScaleFactor')
    else:
        logging.debug('Tying to assign a frequencyScaleFactor according to the '
                      'level of theory {0}'.format(level_of_theory))
        frequency_scale_factor = assign_frequency_scale_factor(level_of_theory)
    use_hindered_rotors = local_context.get('useHinderedRotors', True)
    use_atom_corrections = local_context.get('useAtomCorrections', True)
    use_bond_corrections = local_context.get('useBondCorrections', False)
    bac_type = local_context.get('bondCorrectionType', 'p')
    use_isodesmic_reactions = local_context.get('useIsodesmicReactions', False)
    reference_sets = local_context.get('referenceSets', None)
    atom_energies = local_context.get('atomEnergies', None)

    directory = os.path.dirname(path)

    for rxn in reaction_dict.values():
        rxn.elementary_high_p = True

    for job in job_list:
        if isinstance(job, StatMechJob):
            job.path = os.path.join(directory, job.path)
            job.level_of_theory = level_of_theory
            job.frequencyScaleFactor = frequency_scale_factor
            job.includeHinderedRotors = use_hindered_rotors
            job.applyAtomEnergyCorrections = use_atom_corrections
            job.applyBondEnergyCorrections = use_bond_corrections
            job.bondEnergyCorrectionType = bac_type
            job.atomEnergies = atom_energies
            job.useIsodesmicReactions = use_isodesmic_reactions
            job.referenceSets = reference_sets
        if isinstance(job, ThermoJob):
            job.arkane_species.author = author
            job.arkane_species.level_of_theory = level_of_theory
            job.arkane_species.model_chemistry = model_chemistry
            job.arkane_species.frequency_scale_factor = frequency_scale_factor
            job.arkane_species.use_hindered_rotors = use_hindered_rotors
            job.arkane_species.use_bond_corrections = use_bond_corrections
            if atom_energies is not None:
                job.arkane_species.atom_energies = atom_energies

    return job_list, reaction_dict, species_dict, transition_state_dict, network_dict, level_of_theory


def process_model_chemistry(model_chemistry):
    """Process the model chemistry string representation

    Args:
        model_chemistry (str, unicode): A representation of the model chemistry in an sp//freq format
                                        e.g., 'CCSD(T)-F12a/aug-cc-pVTZ//B3LYP/6-311++G(3df,3pd)',
                                        or a composite method, e.g. 'CBS-QB3'.

    Returns:
        LevelOfTheory, CompositeLevelOfTheory: The level of theory
    """
    if not model_chemistry:
        return model_chemistry
    if isinstance(model_chemistry, LOT):
        return model_chemistry
    if model_chemistry.count('//') > 1:
        raise InputError('The model chemistry seems wrong. It should either be a composite method (like CBS-QB3) '
                         'or of the form sp//geometry, e.g., CCSD(T)-F12a/aug-cc-pVTZ//B3LYP/6-311++G(3df,3pd), '
                         'and should not contain more than one appearance of "//".\n'
                         'Got: {0}'.format(model_chemistry))
    return model_chem_to_lot(model_chemistry)
