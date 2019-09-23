#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
from rmgpy.transport import TransportData

from arkane.common import is_pdep
from arkane.explorer import ExplorerJob
from arkane.kinetics import KineticsJob
from arkane.pdep import PressureDependenceJob
from arkane.statmech import StatMechJob, assign_frequency_scale_factor
from arkane.thermo import ThermoJob

################################################################################

speciesDict = {}
transitionStateDict = {}
reactionDict = {}
networkDict = {}
jobList = []

################################################################################


def database(thermoLibraries=None, transportLibraries=None, reactionLibraries=None, frequenciesLibraries=None,
             kineticsFamilies='default', kineticsDepositories='default', kineticsEstimator='rate rules'):
    """Load the RMG database"""
    if isinstance(thermoLibraries, str):
        thermoLibraries = [thermoLibraries]
    if isinstance(transportLibraries, str):
        transportLibraries = [transportLibraries]
    if isinstance(reactionLibraries, str):
        reactionLibraries = [reactionLibraries]
    if isinstance(frequenciesLibraries, str):
        frequenciesLibraries = [frequenciesLibraries]

    databaseDirectory = settings['database.directory']
    thermoLibraries = thermoLibraries or []
    transportLibraries = transportLibraries
    reactionLibraries = reactionLibraries or []
    kineticsEstimator = kineticsEstimator

    if kineticsDepositories == 'default':
        kineticsDepositories = ['training']
    elif kineticsDepositories == 'all':
        kineticsDepositories = None
    else:
        if not isinstance(kineticsDepositories, list):
            raise InputError(
                "kineticsDepositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe'].")
        kineticsDepositories = kineticsDepositories

    if kineticsFamilies in ('default', 'all', 'none'):
        kineticsFamilies = kineticsFamilies
    else:
        if not isinstance(kineticsFamilies, list):
            raise InputError(
                "kineticsFamilies should be either 'default', 'all', 'none', or a list of names eg. "
                "['H_Abstraction','R_Recombination'] or ['!Intra_Disproportionation'].")
        kineticsFamilies = kineticsFamilies

    rmg_database = get_db() or RMGDatabase()

    rmg_database.load(
        path=databaseDirectory,
        thermoLibraries=thermoLibraries,
        transportLibraries=transportLibraries,
        reactionLibraries=reactionLibraries,
        seedMechanisms=[],
        kineticsFamilies=kineticsFamilies,
        kineticsDepositories=kineticsDepositories,
        depository=False,  # Don't bother loading the depository information, as we don't use it
    )

    for family in rmg_database.kinetics.families.values():  # load training
        family.add_rules_from_training(thermoDatabase=rmg_database.thermo)

    for family in rmg_database.kinetics.families.values():
        family.fill_rules_by_averaging_up(verbose=True)


def species(label, *args, **kwargs):
    """Load a species from an input file"""
    global speciesDict, jobList
    if label in speciesDict:
        raise ValueError('Multiple occurrences of species with label {0!r}.'.format(label))
    logging.info('Loading species {0}...'.format(label))

    spec = Species(label=label)
    speciesDict[label] = spec

    path = None
    if len(args) == 1:
        # The argument is a path to a conformer input file
        path = args[0]
        job = StatMechJob(species=spec, path=path)
        logging.debug('Added species {0} to a stat mech job.'.format(label))
        jobList.append(job)
    elif len(args) > 1:
        raise InputError('species {0} can only have two non-keyword argument '
                         'which should be the species label and the '
                         'path to a quantum file.'.format(spec.label))

    if len(kwargs) > 0:
        # The species parameters are given explicitly
        structure = None
        E0 = None
        modes = []
        spinMultiplicity = 0
        opticalIsomers = 1
        molecularWeight = None
        collisionModel = None
        energyTransferModel = None
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
                spinMultiplicity = value
            elif key == 'opticalIsomers':
                opticalIsomers = value
            elif key == 'molecularWeight':
                molecularWeight = value
            elif key == 'collisionModel':
                collisionModel = value
            elif key == 'energyTransferModel':
                energyTransferModel = value
            elif key == 'thermo':
                thermo = value
            elif key == 'reactive':
                reactive = value
            else:
                raise TypeError('species() got an unexpected keyword argument {0!r}.'.format(key))

        if structure:
            spec.molecule = [structure]
        spec.conformer = Conformer(E0=E0, modes=modes, spin_multiplicity=spinMultiplicity, optical_isomers=opticalIsomers)
        if molecularWeight is not None:
            spec.molecularWeight = molecularWeight
        elif spec.molecularWeight is None and is_pdep(jobList):
            # If a structure was given, simply calling spec.molecularWeight will calculate the molecular weight
            # If one of the jobs is pdep and no molecular weight is given or calculated, raise an error
            raise ValueError("No molecularWeight was entered for species {0}. Since a structure wasn't given"
                             " as well, the molecularWeight, which is important for pressure dependent jobs,"
                             " cannot be reconstructed.".format(spec.label))
        spec.transportData = collisionModel
        spec.energyTransferModel = energyTransferModel
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

        if not energyTransferModel:
            # default to RMG's method of generating energyTransferModel
            spec.generate_energy_transfer_model()

    return spec


def transitionState(label, *args, **kwargs):
    """Load a transition state from an input file"""
    global transitionStateDict
    if label in transitionStateDict:
        raise ValueError('Multiple occurrences of transition state with label {0!r}.'.format(label))
    logging.info('Loading transition state {0}...'.format(label))
    ts = TransitionState(label=label)
    transitionStateDict[label] = ts

    if len(args) == 1 and len(kwargs) == 0:
        # The argument is a path to a conformer input file
        path = args[0]
        job = StatMechJob(species=ts, path=path)
        jobList.append(job)

    elif len(args) == 0:
        # The species parameters are given explicitly
        E0 = None
        modes = []
        spinMultiplicity = 1
        opticalIsomers = 1
        frequency = None
        for key, value in kwargs.items():
            if key == 'E0':
                E0 = value
            elif key == 'modes':
                modes = value
            elif key == 'spinMultiplicity':
                spinMultiplicity = value
            elif key == 'opticalIsomers':
                opticalIsomers = value
            elif key == 'frequency':
                frequency = value
            else:
                raise TypeError('transition_state() got an unexpected keyword argument {0!r}.'.format(key))

        ts.conformer = Conformer(E0=E0, modes=modes, spin_multiplicity=spinMultiplicity, optical_isomers=opticalIsomers)
        ts.frequency = frequency
    else:
        if len(args) == 0 and len(kwargs) == 0:
            raise InputError(
                'The transition_state needs to reference a quantum job file or contain kinetic information.')
        raise InputError('The transition_state can only link a quantum job or directly input information, not both.')

    return ts


def reaction(label, reactants, products, transitionState=None, kinetics=None, tunneling=''):
    """Load a reaction from an input file"""
    global reactionDict, speciesDict, transitionStateDict
    if label in reactionDict:
        label = label + transitionState
        if label in reactionDict:
            raise ValueError('Multiple occurrences of reaction with label {0!r}.'.format(label))
    logging.info('Loading reaction {0}...'.format(label))
    reactants = sorted([speciesDict[spec] for spec in reactants])
    products = sorted([speciesDict[spec] for spec in products])
    if transitionState:
        transitionState = transitionStateDict[transitionState]
    if tunneling.lower() == 'wigner':
        transitionState.tunneling = Wigner(frequency=None)
    elif tunneling.lower() == 'eckart':
        transitionState.tunneling = Eckart(frequency=None, E0_reac=None, E0_TS=None, E0_prod=None)
    elif transitionState and (tunneling == '' or tunneling is None):
        transitionState.tunneling = None
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
            model.verboseComments = True
            for r in rxn:
                model.apply_kinetics_to_reaction(r)

    if isinstance(rxn, Reaction):
        reactionDict[label] = rxn
    else:
        for i in range(len(rxn)):
            reactionDict[label + str(i)] = rxn[i]

    return rxn


def network(label, isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None):
    """Load a network from an input file"""
    global networkDict, speciesDict, reactionDict
    logging.info('Loading network {0}...'.format(label))
    isomers0 = isomers or []
    isomers = []
    for isomer in isomers0:
        if isinstance(isomer, (list, tuple)):
            raise ValueError('Only one species can be present in a unimolecular isomer.')
        isomers.append(speciesDict[isomer])

    reactants0 = reactants or []
    reactants = []
    for reactant in reactants0:
        if not isinstance(reactant, (list, tuple)):
            reactant = [reactant]
        reactants.append(sorted([speciesDict[spec] for spec in reactant]))

    if pathReactions is None:
        # Only add reactions that match reactants and/or isomers
        pathReactions = []
        for rxn in reactionDict.values():
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
                pathReactions.append(rxn)
        logging.debug('Path reactions {} were found for network {}'.format([rxn.label for rxn in pathReactions], label))
    else:
        pathReactions0 = pathReactions
        pathReactions = []
        for rxn in pathReactions0:
            pathReactions.append(reactionDict[rxn])

    if products is None:
        # Figure out which configurations are isomers, reactant channels, and product channels
        products = []
        for rxn in pathReactions:
            # Sort bimolecular configurations so that we always encounter them in the
            # same order
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
            products.append(sorted([speciesDict[spec] for spec in product]))

    isomers = [Configuration(species) for species in isomers]
    reactants = [Configuration(*species) for species in reactants]
    products = [Configuration(*species) for species in products]

    bathGas0 = bathGas or {}
    bathGas = {}
    for spec, fraction in bathGas0.items():
        bathGas[speciesDict[spec]] = fraction

    network = Network(
        label=label,
        isomers=isomers,
        reactants=reactants,
        products=products,
        path_reactions=pathReactions,
        bath_gas=bathGas,
    )
    networkDict[label] = network


def kinetics(label, Tmin=None, Tmax=None, Tlist=None, Tcount=0, sensitivity_conditions=None):
    """Generate a kinetics job"""
    global jobList, reactionDict
    try:
        rxn = reactionDict[label]
    except KeyError:
        raise ValueError('Unknown reaction label {0!r} for kinetics() job.'.format(label))
    job = KineticsJob(reaction=rxn, Tmin=Tmin, Tmax=Tmax, Tcount=Tcount, Tlist=Tlist,
                      sensitivity_conditions=sensitivity_conditions)
    jobList.append(job)


def statmech(label):
    """Generate a statmech job"""
    global jobList, speciesDict, transitionStateDict
    if label in speciesDict or label in transitionStateDict:
        for job in jobList:
            if job.species.label == label:
                break
        else:
            raise ValueError('Could not create StatMechJob for {0!r}; no path specified.'.format(label))
    else:
        raise ValueError('Unknown species or transition state label {0!r} for statmech() job.'.format(label))


def thermo(label, thermoClass):
    """Generate a thermo job"""
    global jobList, speciesDict
    try:
        spec = speciesDict[label]
    except KeyError:
        raise ValueError('Unknown species label {0!r} for thermo() job.'.format(label))
    job = ThermoJob(species=spec, thermoClass=thermoClass)
    jobList.append(job)


def pressureDependence(label, Tmin=None, Tmax=None, Tcount=0, Tlist=None, Pmin=None, Pmax=None, Pcount=0, Plist=None,
                       maximumGrainSize=None, minimumGrainCount=0, method=None, interpolationModel=None,
                       activeKRotor=True, activeJRotor=True, rmgmode=False, sensitivity_conditions=None):
    """Generate a pressure dependent job"""
    global jobList, networkDict

    if isinstance(interpolationModel, str):
        interpolationModel = (interpolationModel,)

    nwk = None
    if label in list(networkDict.keys()):
        nwk = networkDict[label]

    job = PressureDependenceJob(network=nwk, Tmin=Tmin, Tmax=Tmax, Tcount=Tcount, Tlist=Tlist,
                                Pmin=Pmin, Pmax=Pmax, Pcount=Pcount, Plist=Plist,
                                maximumGrainSize=maximumGrainSize, minimumGrainCount=minimumGrainCount,
                                method=method, interpolationModel=interpolationModel,
                                activeKRotor=activeKRotor, activeJRotor=activeJRotor,
                                rmgmode=rmgmode, sensitivity_conditions=sensitivity_conditions)
    jobList.append(job)


def explorer(source, explore_tol=0.01, energy_tol=np.inf, flux_tol=0.0, bathGas=None, maximumRadicalElectrons=np.inf):
    """Generate an explorer job"""
    global jobList, speciesDict
    for job in jobList:
        if isinstance(job, PressureDependenceJob):
            pdepjob = job
            break
    else:
        raise InputError('the explorer block must occur after the pressureDependence block')

    source = [speciesDict[name] for name in source]

    if bathGas:
        bathGas0 = bathGas or {}
        bathGas = {}
        for spec, fraction in bathGas0.items():
            bathGas[speciesDict[spec]] = fraction

    job = ExplorerJob(source=source, pdepjob=pdepjob, explore_tol=explore_tol,
                      energy_tol=energy_tol, flux_tol=flux_tol, bathGas=bathGas,
                      maximumRadicalElectrons=maximumRadicalElectrons)
    jobList.append(job)


def SMILES(smiles):
    """Make a Molecule object from SMILES"""
    return Molecule().from_smiles(smiles)


def adjacencyList(adj):
    """Make a Molecule object from an adjacency list"""
    return Molecule().from_adjacency_list(adj)


def InChI(inchi):
    """Make a Molecule object from InChI"""
    return Molecule().from_inchi(inchi)


def loadNecessaryDatabases():
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


def loadInputFile(path):
    """
    Load the Arkane input file located at `path` on disk, and return a list of
    the jobs defined in that file.
    """
    global speciesDict, transitionStateDict, reactionDict, networkDict, jobList

    # Clear module-level variables
    speciesDict = {}
    transitionStateDict = {}
    reactionDict = {}
    networkDict = {}
    jobList = []

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
        'statmech': statmech,
        'thermo': thermo,
        'pressureDependence': pressureDependence,
        'explorer': explorer,
        # Miscellaneous
        'SMILES': SMILES,
        'adjacencyList': adjacencyList,
        'InChI': InChI,
    }

    loadNecessaryDatabases()

    with open(path, 'r') as f:
        try:
            exec(f.read(), global_context, local_context)
        except (NameError, TypeError, SyntaxError):
            logging.error('The input file {0!r} was invalid:'.format(path))
            raise

    model_chemistry = local_context.get('modelChemistry', '').lower()
    sp_level, freq_level = process_model_chemistry(model_chemistry)

    author = local_context.get('author', '')
    if 'frequencyScaleFactor' in local_context:
        frequency_scale_factor = local_context.get('frequencyScaleFactor')
    else:
        logging.debug('Tying to assign a frequencyScaleFactor according to the frequency '
                      'level of theory {0}'.format(freq_level))
        frequency_scale_factor = assign_frequency_scale_factor(freq_level)
    use_hindered_rotors = local_context.get('useHinderedRotors', True)
    use_atom_corrections = local_context.get('useAtomCorrections', True)
    use_bond_corrections = local_context.get('useBondCorrections', False)
    bac_type = local_context.get('bondCorrectionType', 'p')
    atom_energies = local_context.get('atomEnergies', None)

    directory = os.path.dirname(path)

    for rxn in reactionDict.values():
        rxn.elementary_high_p = True

    for job in jobList:
        if isinstance(job, StatMechJob):
            job.path = os.path.join(directory, job.path)
            job.modelChemistry = sp_level
            job.frequencyScaleFactor = frequency_scale_factor
            job.includeHinderedRotors = use_hindered_rotors
            job.applyAtomEnergyCorrections = use_atom_corrections
            job.applyBondEnergyCorrections = use_bond_corrections
            job.bondEnergyCorrectionType = bac_type
            job.atomEnergies = atom_energies
        if isinstance(job, ThermoJob):
            job.arkane_species.author = author
            job.arkane_species.level_of_theory = model_chemistry
            job.arkane_species.frequency_scale_factor = frequency_scale_factor
            job.arkane_species.use_hindered_rotors = use_hindered_rotors
            job.arkane_species.use_bond_corrections = use_bond_corrections
            if atom_energies is not None:
                job.arkane_species.atom_energies = atom_energies

    return jobList, reactionDict, speciesDict, transitionStateDict, networkDict


def process_model_chemistry(model_chemistry):
    """Process the model chemistry string representation

    Args:
        model_chemistry (str, unicode): A representation of the model chemistry in an sp//freq format
                                        e.g., 'CCSD(T)-F12a/aug-cc-pVTZ//B3LYP/6-311++G(3df,3pd)',
                                        or a composite method, e.g. 'CBS-QB3'.

    Returns:
        str, unicode: The single point energy level of theory
        str, unicode: The frequency level of theory
    """
    if model_chemistry.count('//') > 1:
        raise InputError('The model chemistry seems wrong. It should either be a composite method (like CBS-QB3) '
                         'or of the form sp//geometry, e.g., CCSD(T)-F12a/aug-cc-pVTZ//B3LYP/6-311++G(3df,3pd), '
                         'and should not contain more than one appearance of "//".\n'
                         'Got: {0}'.format(model_chemistry))
    elif '//' in model_chemistry:
        # assume this is an sp//freq format, split
        sp_level, freq_level = model_chemistry.split('//')
    else:
        # assume the sp and freq levels are the same, assign the model chemistry to both
        # (this could also be a composite method, and we'll expect the same behavior)
        sp_level = freq_level = model_chemistry
    return sp_level, freq_level
