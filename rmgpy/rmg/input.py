#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
import os
from copy import deepcopy

import numpy as np

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.exceptions import DatabaseError, InputError
from rmgpy.molecule import Molecule
from rmgpy.molecule.group import Group
from rmgpy.quantity import Quantity, Energy, RateCoefficient, SurfaceConcentration
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.termination import TerminationTime, TerminationConversion, TerminationRateRatio
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.solver.mbSampled import MBSampledReactor
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.surface import SurfaceReactor
from rmgpy.util import as_list
from rmgpy.data.surface import MetalDatabase
from rmgpy.rmg.reactors import Reactor, ConstantVIdealGasReactor, ConstantTLiquidSurfaceReactor, ConstantTVLiquidReactor, ConstantTPIdealGasReactor
from rmgpy.data.vaporLiquidMassTransfer import liquidVolumetricMassTransferCoefficientPowerLaw
from rmgpy.molecule.fragment import Fragment

################################################################################

rmg = None
species_dict = {}


def database(
        thermoLibraries=None,
        transportLibraries=None,
        reactionLibraries=None,
        frequenciesLibraries=None,
        seedMechanisms=None,
        kineticsFamilies='default',
        kineticsDepositories='default',
        kineticsEstimator='rate rules',
):
    # This function just stores the information about the database to be loaded
    # We don't actually load the database until after we're finished reading
    # the input file
    rmg.database_directory = settings['database.directory']
    rmg.thermo_libraries = as_list(thermoLibraries, default=[])
    rmg.transport_libraries = as_list(transportLibraries, default=None)

    # Modify reaction library list such that all entries are tuples
    reaction_libraries = as_list(reactionLibraries, default=[])
    rmg.reaction_libraries = [(name, False) if not isinstance(name, tuple) else name for name in reaction_libraries]

    rmg.seed_mechanisms = as_list(seedMechanisms, default=[])
    rmg.statmech_libraries = as_list(frequenciesLibraries, default=[])
    rmg.kinetics_estimator = kineticsEstimator

    if kineticsDepositories == 'default':
        rmg.kinetics_depositories = ['training']
    elif kineticsDepositories == 'all':
        rmg.kinetics_depositories = None
    else:
        if not isinstance(kineticsDepositories, list):
            raise InputError("kinetics_depositories should be either 'default', 'all', or a list of names eg. "
                             "['training','PrIMe'].")
        rmg.kinetics_depositories = kineticsDepositories

    if kineticsFamilies in ('default', 'all', 'none'):
        rmg.kinetics_families = kineticsFamilies
    else:
        if not isinstance(kineticsFamilies, list):
            raise InputError("kineticsFamilies should be either 'default', 'all', 'none', or a list of names eg. "
                             "['H_Abstraction','R_Recombination'] or ['!Intra_Disproportionation'].")
        rmg.kinetics_families = kineticsFamilies


def catalyst_properties(bindingEnergies=None,
                        surfaceSiteDensity=None,
                        metal=None,
                        coverageDependence=False):
    """
    Specify the properties of the catalyst.
    Binding energies of C,H,O,N atoms, and the surface site density.
    Metal is the label of a metal in the surface metal library.
    Defaults to Pt(111) if not specified.
    coverageDependence defaults to False if not specified. If True then coverage dependent kinetics from libraries are used.
    """
    # Normally we wouldn't load the database until after reading the input file,
    # but we need to load the metal surfaces library to validate the input.
    metal_db = MetalDatabase()
    metal_db.load(os.path.join(settings['database.directory'], 'surface'))

    if metal and (bindingEnergies or surfaceSiteDensity):
        raise InputError("In catalyst_properties section you should only specify a 'metal' shortcut " 
                         "or the surfaceSiteDensity and bindingEnergies, but not both.")

    if metal:
        try:
            logging.info("Using catalyst surface properties from metal %r.", metal)
            rmg.binding_energies = metal_db.get_binding_energies(metal)
            rmg.surface_site_density = metal_db.get_surface_site_density(metal)
        except DatabaseError:
            logging.error('Metal %r missing from surface library. Please specify both metal and facet.', metal)
            raise
    else: # metal not specified
        if bindingEnergies is None:
            rmg.binding_energies = metal_db.get_binding_energies("Pt111")
            logging.info("Using default binding energies, Pt(111)")
        else:
            rmg.binding_energies = convert_binding_energies(bindingEnergies)
    

        if surfaceSiteDensity is None:
            rmg.surface_site_density = metal_db.get_surface_site_density("Pt111")
            logging.info("Using default surface site density, Pt(111)")
        else:
            rmg.surface_site_density = SurfaceConcentration(*surfaceSiteDensity)

    logging.info("Using binding energies:\n%r", rmg.binding_energies)
    logging.info("Using surface site density: %r", rmg.surface_site_density)

    if coverageDependence:
        logging.info("Coverage dependence is turned ON")
    else:
        logging.info("Coverage dependence is turned OFF")
    rmg.coverage_dependence = coverageDependence

def convert_binding_energies(binding_energies):
    """
    Process binding_energies dictionary from the input file
    
    It converts the values into Energy quantities, and checks that 
    all elements C,H,O, and N are present.

    :param binding_energies: a dictionary of element symbol: binding energy pairs
    :return: the processed and checked dictionary
    """
    if not isinstance(binding_energies, dict):
        raise InputError("bindingEnergies should be None (for default) or a dict.")
    new_dict = {}
    for element in 'CHON':
        try:
            new_dict[element] = Energy(binding_energies[element])
        except KeyError:
            logging.error('Element {} missing from binding energies dictionary'.format(element))
            raise
    return new_dict


def species(label, structure, reactive=True, cut=False, size_threshold=None):
    logging.debug('Found {0} species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive',
                                                         label,
                                                         structure.to_smiles()))

    if '+' in label:
        raise InputError('species {0} label cannot include a + sign'.format(label))

    if cut:
        mol_to_frag[label] = {} # key:original molecule label, value:created fragment label
        cut_frag_list = structure.cut_molecule(size_threshold=size_threshold)
        logging.info('The original molecule {0} is divided into several fragments:'.format(label))
        for initial_frag in cut_frag_list:
            frag_label = initial_frag.to_smiles()

            logging.info(frag_label)

            if frag_label not in mol_to_frag[label]:
                mol_to_frag[label][frag_label] = 1
            else:
                mol_to_frag[label][frag_label] += 1

            spec, is_new = rmg.reaction_model.make_new_species(initial_frag, label=frag_label, reactive=reactive)
            # if duplicated fragment is created, combine and add to initialMoleFractions
            if not is_new:
                continue

            rmg.initial_species.append(spec)
            species_dict[frag_label] = spec
    else:
        try:
            spec, is_new = rmg.reaction_model.make_new_species(structure, label=label, reactive=reactive)
        except:
            logging.error(f'Error when reading species "{label}" from input file.')
            raise
        if not is_new:
            raise InputError("Species {0} is a duplicate of {1}. Species in input file must be unique".format(label,
                                                                                                            spec.label))

        rmg.initial_species.append(spec)
        species_dict[label] = spec

def forbidden(label, structure):

    if '+' in label:
        raise InputError('forbidden species {0} label cannot include a + sign'.format(label))
    rmg.forbidden_structures.append(Entry(label=label,item=structure))

def smarts(string):
    return Molecule().from_smarts(string)


def smiles(string):
    return Molecule().from_smiles(string)


def inchi(string):
    return Molecule().from_inchi(string)


def adjacency_list(string):
    return Molecule().from_adjacency_list(string)

def adjacency_list_group(string):
    return Group().from_adjacency_list(string)

def react(tups):
    if not isinstance(tups, list):
        raise InputError("React takes a list of tuples of species strings.")
    for item in tups:
        if not isinstance(item, tuple):
            raise InputError("React takes a list of tuples of species strings.")
        for it in item:
            if not isinstance(it, str):
                raise InputError("React takes a list of tuples of species strings.")
    rmg.init_react_tuples = tups

def fragment_adj(string):
    return Fragment().from_adjacency_list(string)

def fragment_smiles(string):
    return Fragment().from_smiles_like_string(string)

# Reaction systems
def simple_reactor(temperature,
                   pressure,
                   initialMoleFractions,
                   nSims=6,
                   terminationConversion=None,
                   terminationTime=None,
                   terminationRateRatio=None,
                   balanceSpecies=None,
                   sensitivity=None,
                   sensitivityThreshold=1e-3,
                   sensitivityTemperature=None,
                   sensitivityPressure=None,
                   sensitivityMoleFractions=None,
                   constantSpecies=None):
    logging.debug('Found SimpleReactor reaction system')

    try:
        if mol_to_frag:
            # calculate total as denominator
            total = float(0)
            for initial_mol, value in initialMoleFractions.items():
                if initial_mol not in mol_to_frag: # there might be other species which set cut to be False but in initialMoleFractions
                    mol_to_frag[initial_mol] = {initial_mol: 1}
                total += value * sum(mol_to_frag[initial_mol].values())
            for key, frag_dict in mol_to_frag.items():
                # if not perform cutting, no need to modify initialMoleFractions
                # only 1 species in system and it does not get cut, then no need to modify initialMoleFractions
                if len(mol_to_frag.keys())==1 and len(frag_dict.keys()) == 1 and key == list(frag_dict)[0]:
                    continue
                # probably this species does not get cut but other speceis are cut,
                # so it still requires to re-calculate and normalize molar fraction
                elif len(frag_dict.keys()) == 1 and key == list(frag_dict)[0]:
                    initialMoleFractions[key] = initialMoleFractions[key] * frag_dict[key] / total
                    continue
                for frag_label, number in frag_dict.items():
                    if frag_label in initialMoleFractions:
                        initialMoleFractions[frag_label] += initialMoleFractions[key] * number / total
                    else:
                        initialMoleFractions[frag_label] = initialMoleFractions[key] * number / total
                del initialMoleFractions[key]
            logging.info('After cutting, new compositions:')
            for spec, molfrac in initialMoleFractions.items():
                logging.info('{0} = {1}'.format(spec, molfrac))
    except NameError:
        pass

    for key, value in initialMoleFractions.items():
        if not isinstance(value, list):
            initialMoleFractions[key] = float(value)
            if value < 0:
                raise InputError('Initial mole fractions cannot be negative.')
        else:
            if len(value) != 2:
                raise InputError("Initial mole fraction values must either be a number or a list with 2 entries")
            initialMoleFractions[key] = [float(value[0]), float(value[1])]
            if value[0] < 0 or value[1] < 0:
                raise InputError('Initial mole fractions cannot be negative.')
            elif value[1] < value[0]:
                raise InputError('Initial mole fraction range out of order: {0}'.format(key))

    if not isinstance(temperature, list):
        T = Quantity(temperature)
    else:
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]

    if not isinstance(pressure, list):
        P = Quantity(pressure)
    else:
        if len(pressure) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        P = [Quantity(p) for p in pressure]

    if not isinstance(temperature, list) and not isinstance(pressure, list) and all(
            [not isinstance(x, list) for x in initialMoleFractions.values()]):
        nSims = 1

    # normalize mole fractions if not using a mole fraction range
    if all([not isinstance(x, list) for x in initialMoleFractions.values()]):
        total_initial_moles = sum(initialMoleFractions.values())
        if total_initial_moles != 1:
            logging.warning('Initial mole fractions do not sum to one; normalizing.')
            logging.info('')
            logging.info('Original composition:')
            for spec, molfrac in initialMoleFractions.items():
                logging.info('{0} = {1}'.format(spec, molfrac))
            for spec in initialMoleFractions:
                initialMoleFractions[spec] /= total_initial_moles
            logging.info('')
            logging.info('Normalized mole fractions:')
            for spec, molfrac in initialMoleFractions.items():
                logging.info('{0} = {1}'.format(spec, molfrac))
            logging.info('')

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            # check whether key is in species_dict (not if molecule cut provided)
            if spec not in species_dict:
                # select first fragment as species used for terminationConversion
                repr_frag = sorted(mol_to_frag[spec].keys())[0]
                termination.append(TerminationConversion(species_dict[repr_frag], conv))
            else:
                termination.append(TerminationConversion(species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    sensitive_species = []
    if sensitivity:
        if sensitivity != 'all':
            if isinstance(sensitivity, str):
                sensitivity = [sensitivity]
            for spec in sensitivity:
                sensitive_species.append(species_dict[spec])

        else:
            sensitive_species.append('all')

    #Check the constant species exist
    if constantSpecies is not None:
        logging.debug('  Generation with constant species:')
        for const_spc in constantSpecies:
            logging.debug("  {0}".format(const_spc))
            if const_spc not in species_dict:
                raise InputError('Species {0} not found in the input file'.format(const_spc))
    
    if not isinstance(T, list):
        sensitivityTemperature = T
    if not isinstance(P, list):
        sensitivityPressure = P
    if not any([isinstance(x, list) for x in initialMoleFractions.values()]):
        sensitivityMoleFractions = deepcopy(initialMoleFractions)
    if sensitivityMoleFractions is None or sensitivityTemperature is None or sensitivityPressure is None:
        sens_conditions = None
    else:
        sens_conditions = sensitivityMoleFractions
        sens_conditions['T'] = Quantity(sensitivityTemperature).value_si
        sens_conditions['P'] = Quantity(sensitivityPressure).value_si

    system = SimpleReactor(T, P, initialMoleFractions, nSims, termination, sensitive_species, sensitivityThreshold, sens_conditions, constantSpecies)
    rmg.reaction_systems.append(system)

    assert balanceSpecies is None or isinstance(balanceSpecies, str), 'balanceSpecies should be the string corresponding to a single species'
    rmg.balance_species = balanceSpecies
    if balanceSpecies:  # check that the balanceSpecies can't be taken to zero
        total = 0.0
        for key, item in initialMoleFractions.items():
            if key == balanceSpecies:
                assert not isinstance(item, list), 'balanceSpecies must not have a defined range'
                xbspcs = item
            if isinstance(item, list):
                total += item[1] - item[0]

        if total > xbspcs:
            raise ValueError('The sum of the differences in the ranged mole fractions is greater than the mole '
                             'fraction of the balance species, this would require the balanceSpecies mole fraction to '
                             'be negative in some cases which is not allowed, either reduce the maximum mole fractions '
                             'or dont use balanceSpecies')

def constant_V_ideal_gas_reactor(temperature,
                   pressure,
                   initialMoleFractions,
                   terminationConversion=None,
                   terminationTime=None,
                   terminationRateRatio=None,
                   balanceSpecies=None):
    logging.debug('Found ConstantVIdealGasReactor reaction system')

    for key, value in initialMoleFractions.items():
        if not isinstance(value, list):
            initialMoleFractions[key] = float(value)
            if value < 0:
                raise InputError('Initial mole fractions cannot be negative.')
        else:
            if len(value) != 2:
                raise InputError("Initial mole fraction values must either be a number or a list with 2 entries")
            initialMoleFractions[key] = [float(value[0]), float(value[1])]
            if value[0] < 0 or value[1] < 0:
                raise InputError('Initial mole fractions cannot be negative.')
            elif value[1] < value[0]:
                raise InputError('Initial mole fraction range out of order: {0}'.format(key))

    if not isinstance(temperature, list):
        T = Quantity(temperature).value_si
    else:
        raise InputError("Condition ranges not supported for this reaction type")
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]

    if not isinstance(pressure, list):
        P = Quantity(pressure).value_si
    else:
        raise InputError("Condition ranges not supported for this reaction type")
        if len(pressure) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        P = [Quantity(p) for p in pressure]

    if not isinstance(temperature, list) and not isinstance(pressure, list) and all(
            [not isinstance(x, list) for x in initialMoleFractions.values()]):
        nSims = 1

    # normalize mole fractions if not using a mole fraction range
    if all([not isinstance(x, list) for x in initialMoleFractions.values()]):
        total_initial_moles = sum(initialMoleFractions.values())
        if total_initial_moles != 1:
            logging.warning('Initial mole fractions do not sum to one; normalizing.')
            logging.info('')
            logging.info('Original composition:')
            for spec, molfrac in initialMoleFractions.items():
                logging.info('{0} = {1}'.format(spec, molfrac))
            for spec in initialMoleFractions:
                initialMoleFractions[spec] /= total_initial_moles
            logging.info('')
            logging.info('Normalized mole fractions:')
            for spec, molfrac in initialMoleFractions.items():
                logging.info('{0} = {1}'.format(spec, molfrac))
            logging.info('')

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append((species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))
    
    initial_cond = initialMoleFractions
    initial_cond["T"] = T 
    initial_cond["P"] = P
    system = ConstantVIdealGasReactor(rmg.reaction_model.core.phase_system,rmg.reaction_model.edge.phase_system,initial_cond,termination)
    system.T = Quantity(T)
    system.P = Quantity(P)
    system.Trange = None
    system.Prange = None
    system.sensitive_species = []
    rmg.reaction_systems.append(system)

    assert balanceSpecies is None or isinstance(balanceSpecies, str), 'balanceSpecies should be the string corresponding to a single species'
    rmg.balance_species = balanceSpecies
    if balanceSpecies:  # check that the balanceSpecies can't be taken to zero
        total = 0.0
        for key, item in initialMoleFractions.items():
            if key == balanceSpecies:
                assert not isinstance(item, list), 'balanceSpecies must not have a defined range'
                xbspcs = item
            if isinstance(item, list):
                total += item[1] - item[0]

        if total > xbspcs:
            raise ValueError('The sum of the differences in the ranged mole fractions is greater than the mole '
                             'fraction of the balance species, this would require the balanceSpecies mole fraction to '
                             'be negative in some cases which is not allowed, either reduce the maximum mole fractions '
                             'or dont use balanceSpecies')

def constant_TP_ideal_gas_reactor(temperature,
                   pressure,
                   initialMoleFractions,
                   terminationConversion=None,
                   terminationTime=None,
                   terminationRateRatio=None,
                   balanceSpecies=None):
    logging.debug('Found ConstantTPIdealGasReactor reaction system')

    for key, value in initialMoleFractions.items():
        if not isinstance(value, list):
            initialMoleFractions[key] = float(value)
            if value < 0:
                raise InputError('Initial mole fractions cannot be negative.')
        else:
            if len(value) != 2:
                raise InputError("Initial mole fraction values must either be a number or a list with 2 entries")
            initialMoleFractions[key] = [float(value[0]), float(value[1])]
            if value[0] < 0 or value[1] < 0:
                raise InputError('Initial mole fractions cannot be negative.')
            elif value[1] < value[0]:
                raise InputError('Initial mole fraction range out of order: {0}'.format(key))

    if not isinstance(temperature, list):
        T = Quantity(temperature).value_si
    else:
        raise InputError("Condition ranges not supported for this reaction type")
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]

    if not isinstance(pressure, list):
        P = Quantity(pressure).value_si
    else:
        raise InputError("Condition ranges not supported for this reaction type")
        if len(pressure) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        P = [Quantity(p) for p in pressure]

    if not isinstance(temperature, list) and not isinstance(pressure, list) and all(
            [not isinstance(x, list) for x in initialMoleFractions.values()]):
        nSims = 1

    # normalize mole fractions if not using a mole fraction range
    if all([not isinstance(x, list) for x in initialMoleFractions.values()]):
        total_initial_moles = sum(initialMoleFractions.values())
        if total_initial_moles != 1:
            logging.warning('Initial mole fractions do not sum to one; normalizing.')
            logging.info('')
            logging.info('Original composition:')
            for spec, molfrac in initialMoleFractions.items():
                logging.info('{0} = {1}'.format(spec, molfrac))
            for spec in initialMoleFractions:
                initialMoleFractions[spec] /= total_initial_moles
            logging.info('')
            logging.info('Normalized mole fractions:')
            for spec, molfrac in initialMoleFractions.items():
                logging.info('{0} = {1}'.format(spec, molfrac))
            logging.info('')

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append((species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    initial_cond = initialMoleFractions
    initial_cond["T"] = T
    initial_cond["P"] = P
    system = ConstantTPIdealGasReactor(rmg.reaction_model.core.phase_system,rmg.reaction_model.edge.phase_system,initial_cond,termination)
    system.T = Quantity(T)
    system.P = Quantity(P)
    system.Trange = None
    system.Prange = None
    system.sensitive_species = []
    rmg.reaction_systems.append(system)

    assert balanceSpecies is None or isinstance(balanceSpecies, str), 'balanceSpecies should be the string corresponding to a single species'
    rmg.balance_species = balanceSpecies
    if balanceSpecies:  # check that the balanceSpecies can't be taken to zero
        total = 0.0
        for key, item in initialMoleFractions.items():
            if key == balanceSpecies:
                assert not isinstance(item, list), 'balanceSpecies must not have a defined range'
                xbspcs = item
            if isinstance(item, list):
                total += item[1] - item[0]

        if total > xbspcs:
            raise ValueError('The sum of the differences in the ranged mole fractions is greater than the mole '
                             'fraction of the balance species, this would require the balanceSpecies mole fraction to '
                             'be negative in some cases which is not allowed, either reduce the maximum mole fractions '
                             'or dont use balanceSpecies')

def liquid_cat_reactor(temperature,
                   initialConcentrations,
                   initialSurfaceCoverages,
                   surfaceVolumeRatio,
                   potential=None,
                   terminationConversion=None,
                   terminationTime=None,
                   terminationRateRatio=None,
                   constantSpecies=[]):
    for spec, conc in initialConcentrations.items():
        if not isinstance(conc, list):
            concentration = Quantity(conc)
            # check the dimensions are ok
            # convert to mol/m^3 (or something numerically nice? or must it be SI)
            initialConcentrations[spec] = concentration.value_si
        else:
            if len(conc) != 2:
                raise InputError("Concentration values must either be in the form of (number,units) or a list with 2 "
                                 "entries of the same format")
            initialConcentrations[spec] = [Quantity(conc[0]), Quantity(conc[1])]

    if not isinstance(temperature, list) and all([not isinstance(x, list) for x in initialConcentrations.values()]):
        nSims = 1

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append(TerminationConversion(species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    if constantSpecies is not None:
        logging.debug('  Generation with constant species:')
        for const_spc in constantSpecies:
            logging.debug("  {0}".format(const_spc))
            if const_spc not in species_dict:
                raise InputError('Species {0} not found in the input file'.format(const_spc))

    if not isinstance(temperature, list):
        T = Quantity(temperature).value_si
    else:
        raise InputError("Condition ranges not supported for this reaction type")
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]

    if not isinstance(temperature, list) and all(
            [not isinstance(x, list) for x in initialConcentrations.values()]) and all(
            [not isinstance(x, list) for x in initialSurfaceCoverages.values()]
            ):
        nSims = 1

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append((species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    initialCondLiq = dict()
    V = 1.0
    A = V*Quantity(surfaceVolumeRatio).value_si
    for key,item in initialConcentrations.items():
        initialCondLiq[key] = item*V
    initialCondLiq["T"] = T
    initialCondLiq["V"] = V
    initialCondSurf = dict()
    for key,item in initialSurfaceCoverages.items():
        initialCondSurf[key] = item*rmg.surface_site_density.value_si*A
    initialCondSurf["T"] = T
    initialCondSurf["A"] = A
    if potential:
        initialCondSurf["phi"] = Quantity(potential).value_si
    system = ConstantTLiquidSurfaceReactor(rmg.reaction_model.core.phase_system,
                                           rmg.reaction_model.edge.phase_system,
                                           {"liquid":initialCondLiq,"surface":initialCondSurf},termination,constantSpecies)
    system.T = Quantity(T)
    system.Trange = None
    system.sensitive_species = []
    rmg.reaction_systems.append(system)

def constant_T_V_liquid_reactor(temperature,
                    initialConcentrations,
                    liquidVolume=None,
                    residenceTime=None,
                    inletVolumetricFlowRate=None,
                    outletVolumetricFlowRate=None,
                    inletConcentrations=None,
                    vaporPressure=None,
                    vaporMoleFractions=None,
                    terminationConversion=None,
                    terminationTime=None,
                    terminationRateRatio=None,
                    constantSpecies=[]):


    ################################################# check input ########################################################

    if not isinstance(temperature, list):
        T = Quantity(temperature).value_si
    else:
        raise InputError("Condition ranges not supported for this reaction type")
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]

    for spec, conc in initialConcentrations.items():
        if not isinstance(conc, list):
            concentration = Quantity(conc)
            initialConcentrations[spec] = concentration.value_si
        else:
            raise InputError("Condition ranges not supported for this reaction type")
            if len(conc) != 2:
                raise InputError("Concentration values must either be in the form of (number,units) or a list with 2 "
                                 "entries of the same format")
            initialConcentrations[spec] = [Quantity(conc[0]), Quantity(conc[1])]

    V = 1.0
    if liquidVolume:
        V = Quantity(liquidVolume).value_si
    logging.debug(f'  Generation with liquid volume {V}')

    inlet_volumetric_flow_rate = None
    outlet_volumetric_flow_rate = None
    residence_time = None
    if residenceTime:
        residence_time = Quantity(residenceTime).value_si
        inlet_volumetric_flow_rate = V/residence_time
        outlet_volumetric_flow_rate = V/residence_time
        logging.debug(f'  Generation with residence time {residence_time} s. Calculated inlet/outlet volumetric flow rate is {inlet_volumetric_flow_rate} m^3/s')

    if inletVolumetricFlowRate:
        if inlet_volumetric_flow_rate:
            if inlet_volumetric_flow_rate != inletVolumetricFlowRate:
                raise InputError('  Inconsistent residence time and inlet volumetric flow rate')

        inlet_volumetric_flow_rate = Quantity(inletVolumetricFlowRate).value_si

    if outletVolumetricFlowRate:
        if outlet_volumetric_flow_rate:
            if outlet_volumetric_flow_rate != outletVolumetricFlowRate:
                raise InputError('  Inconsistent residence time and inlet volumetric flow rate')
            
        outlet_volumetric_flow_rate = Quantity(outletVolumetricFlowRate).value_si

    if inletConcentrations:
        if not inlet_volumetric_flow_rate:
            raise InputError("Inlet volumetric flow rate or residence time must be provided when inlet concentrations are specified")

        for spec, conc in inletConcentrations.items():
            if not isinstance(conc, list):
                concentration = Quantity(conc)
                inletConcentrations[spec] = concentration.value_si
            else:
                raise InputError("Condition ranges not supported for this reaction type")
                if len(conc) != 2:
                    raise InputError("Concentration values must either be in the form of (number,units) or a list with 2 "
                                    "entries of the same format")
                inletConcentrations[spec] = [Quantity(conc[0]), Quantity(conc[1])]
    else:
        if residence_time or inlet_volumetric_flow_rate:
            inletConcentrations = initialConcentrations
            logging.debug(f'  Inlet concentrations not provided. Using initial concentrations as inlet concentrations.')

    vapor_pressure = None
    if vaporPressure or vaporMoleFractions:

        if not (vaporPressure and vaporMoleFractions):
            raise ImportError("Vapor pressure, vapor mole fractions and liquid volumetric mass transfer coefficient power law"
                              "must be provided to simulate evaporation/condensation with kLA.")

        if len(vaporPressure) != 2:
            raise InputError("Vapor pressure value must be in the form of (number, units).")
        vapor_pressure = Quantity(vaporPressure).value_si
        logging.debug(f'  Generation with vapor pressure {vapor_pressure} Pa')

        for spec in vaporMoleFractions:
            vaporMoleFractions[spec] = float(vaporMoleFractions[spec])

        total_vapor_moles = sum(vaporMoleFractions.values())
        if total_vapor_moles != 1:
            logging.warning('Initial mole fractions do not sum to one; normalizing.')
            logging.info('')
            logging.info('Original composition:')
            for spec, molfrac in vaporMoleFractions.items():
                logging.info(f"{spec} = {molfrac}")
            for spec in vaporMoleFractions:
                vaporMoleFractions[spec] /= total_vapor_moles
            logging.info('')
            logging.info('Normalized mole fractions:')
            for spec, molfrac in vaporMoleFractions.items():
                logging.info(f"{spec} = {molfrac}")

    if not isinstance(temperature, list) and all([not isinstance(x, list) for x in initialConcentrations.values()]):
        nSims = 1

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append(TerminationConversion(species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    if constantSpecies is not None:
        logging.debug('  Generation with constant species:')
        for const_spc in constantSpecies:
            logging.debug("  {0}".format(const_spc))
            if const_spc not in species_dict:
                raise InputError('Species {0} not found in the input file'.format(const_spc))

    ############################################### process inputs ##############################################

    initial_conditions = dict()
    for key, item in initialConcentrations.items():
        initial_conditions[key] = item*V
    initial_conditions["T"] = T
    initial_conditions["V"] = V

    inlet_conditions = dict()
    if inletConcentrations:
        total_molar_flow_rate = 0
        for key, item in inletConcentrations.items():
            inlet_conditions[key] = item*inlet_volumetric_flow_rate
            total_molar_flow_rate += inlet_conditions[key]
        for key, item in inlet_conditions.items():
            inlet_conditions[key] = item/total_molar_flow_rate #molar fraction for each species
        inlet_conditions["T"] = initial_conditions["T"]
        inlet_conditions["P"] = 1e8
        inlet_conditions["F"] = total_molar_flow_rate

    outlet_conditions = dict()
    if outlet_volumetric_flow_rate:
        outlet_conditions["Vout"] = outlet_volumetric_flow_rate

    evap_cond_conditions = dict()
    if vaporMoleFractions:
        for key, item in vaporMoleFractions.items():
            evap_cond_conditions[key] = item
        evap_cond_conditions["P"] = vapor_pressure
        evap_cond_conditions["T"] = initial_conditions["T"]
    
    system = ConstantTVLiquidReactor(rmg.reaction_model.core.phase_system,
                                           rmg.reaction_model.edge.phase_system,
                                           initial_conditions,
                                           termination,
                                           constantSpecies,
                                           inlet_conditions,
                                           outlet_conditions,
                                           evap_cond_conditions,
                                           )
    system.T = Quantity(T)
    system.Trange = None
    system.sensitive_species = []
    rmg.reaction_systems.append(system)

# Reaction systems
def liquid_reactor(temperature,
                   initialConcentrations,
                   terminationConversion=None,
                   nSims=4,
                   terminationTime=None,
                   terminationRateRatio=None,
                   sensitivity=None,
                   sensitivityThreshold=1e-3,
                   sensitivityTemperature=None,
                   sensitivityConcentrations=None,
                   constantSpecies=None):
    logging.debug('Found LiquidReactor reaction system')

    if not isinstance(temperature, list):
        T = Quantity(temperature)
    else:
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list '
                             'with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]

    for spec, conc in initialConcentrations.items():
        if not isinstance(conc, list):
            concentration = Quantity(conc)
            # check the dimensions are ok
            # convert to mol/m^3 (or something numerically nice? or must it be SI)
            initialConcentrations[spec] = concentration.value_si
        else:
            if len(conc) != 2:
                raise InputError("Concentration values must either be in the form of (number,units) or a list with 2 "
                                 "entries of the same format")
            initialConcentrations[spec] = [Quantity(conc[0]), Quantity(conc[1])]

    if not isinstance(temperature, list) and all([not isinstance(x, list) for x in initialConcentrations.values()]):
        nSims = 1

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append(TerminationConversion(species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    sensitive_species = []
    if sensitivity:
        for spec in sensitivity:
            sensitive_species.append(species_dict[spec])

    # chatelak: check the constant species exist
    if constantSpecies is not None:
        logging.debug('  Generation with constant species:')
        for const_spc in constantSpecies:
            logging.debug("  {0}".format(const_spc))
            if const_spc not in species_dict:
                raise InputError('Species {0} not found in the input file'.format(const_spc))

    if not isinstance(T, list):
        sensitivityTemperature = T
    if not any([isinstance(x, list) for x in initialConcentrations.values()]):
        sensitivityConcentrations = initialConcentrations
    if sensitivityConcentrations is None or sensitivityTemperature is None:
        sens_conditions = None
    else:
        sens_conditions = sensitivityConcentrations
        sens_conditions['T'] = Quantity(sensitivityTemperature).value_si

    system = LiquidReactor(T, initialConcentrations, nSims, termination, sensitive_species, sensitivityThreshold,
                           sens_conditions, constantSpecies)
    rmg.reaction_systems.append(system)


# Reaction systems
def surface_reactor(temperature,
                    initialPressure,
                    initialGasMoleFractions,
                    initialSurfaceCoverages,
                    surfaceVolumeRatio,
                    nSims=4,
                    terminationConversion=None,
                    terminationTime=None,
                    terminationRateRatio=None,
                    sensitivity=None,
                    sensitivityThreshold=1e-3):
    logging.debug('Found SurfaceReactor reaction system')

    for value in list(initialGasMoleFractions.values()):
        if value < 0:
            raise InputError('Initial mole fractions cannot be negative.')
    total_initial_moles = sum(initialGasMoleFractions.values())
    if total_initial_moles != 1:
        logging.warning('Initial gas mole fractions do not sum to one; renormalizing.')
        logging.debug('')
        logging.debug('Original composition:')
        for spec, molfrac in initialGasMoleFractions.items():
            logging.debug("{0} = {1}".format(spec, molfrac))
        for spec in initialGasMoleFractions:
            initialGasMoleFractions[spec] /= total_initial_moles
        logging.info('')
        logging.debug('Normalized mole fractions:')
        for spec, molfrac in initialGasMoleFractions.items():
            logging.debug("{0} = {1}".format(spec, molfrac))

    if not isinstance(temperature, list):
        T = Quantity(temperature)
    else:
        if len(temperature) != 2:
            raise InputError('Temperature ranges can either be in the form of (number,units) or a list with 2 entries '
                             'of the same format')
        T = [Quantity(t) for t in temperature]

    if not isinstance(initialPressure, list):
        P_initial = Quantity(initialPressure)
    else:
        if len(initialPressure) != 2:
            raise InputError('Initial pressure ranges can either be in the form ''of (number,units) or a list with '
                             '2 entries of the same format')
        P_initial = [Quantity(p) for p in initialPressure]

    if not isinstance(temperature, list) and not isinstance(initialPressure, list):
        nSims = 1
    if any([isinstance(x, list) for x in initialGasMoleFractions.values()]) or \
            any([isinstance(x, list) for x in initialSurfaceCoverages.values()]):
        raise NotImplementedError("Can't do ranges on species concentrations for surface reactors yet.")

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append(TerminationConversion(species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    sensitive_species = []
    if sensitivity:
        for spec in sensitivity:
            sensitive_species.append(species_dict[spec])
    if not isinstance(T, list):
        sensitivityTemperature = T
    if not isinstance(initialPressure, list):
        sensitivityPressure = initialPressure
    sens_conditions = None
    if sensitivity:
        raise NotImplementedError("Can't currently do sensitivity with surface reactors.")
        # The problem is inside base.pyx it reads the dictionary 'sensConditions'
        # and guesses whether they're all concentrations (liquid reactor) or
        # mole fractions (simple reactor). In fact, some may be surface coverages.

    system = SurfaceReactor(T=T,
                            P_initial=P_initial,
                            initial_gas_mole_fractions=initialGasMoleFractions,
                            initial_surface_coverages=initialSurfaceCoverages,
                            surface_volume_ratio=surfaceVolumeRatio,
                            surface_site_density=rmg.surface_site_density,
                            n_sims=nSims,
                            termination=termination,
                            sensitive_species=sensitive_species,
                            sensitivity_threshold=sensitivityThreshold,
                            sens_conditions=sens_conditions,
                            coverage_dependence=rmg.coverage_dependence)
    rmg.reaction_systems.append(system)
    system.log_initial_conditions(number=len(rmg.reaction_systems))


# Reaction systems
def mb_sampled_reactor(temperature,
                       pressure,
                       initialMoleFractions,
                       mbsamplingRate,
                       terminationConversion=None,
                       terminationTime=None,
                       sensitivity=None,
                       sensitivityThreshold=1e-3,
                       constantSpecies=None,
                       ):
    logging.debug('Found MBSampledReactor reaction system')

    for value in initialMoleFractions.values():
        if value < 0:
            raise InputError('Initial mole fractions cannot be negative.')

    for spec in initialMoleFractions:
        initialMoleFractions[spec] = float(initialMoleFractions[spec])

    total_initial_moles = sum(initialMoleFractions.values())
    if total_initial_moles != 1:
        logging.warning('Initial mole fractions do not sum to one; normalizing.')
        logging.info('')
        logging.info('Original composition:')
        for spec, molfrac in initialMoleFractions.items():
            logging.info("{0} = {1}".format(spec, molfrac))
        for spec in initialMoleFractions:
            initialMoleFractions[spec] /= total_initial_moles
        logging.info('')
        logging.info('Normalized mole fractions:')
        for spec, molfrac in initialMoleFractions.items():
            logging.info("{0} = {1}".format(spec, molfrac))

    T = Quantity(temperature)
    P = Quantity(pressure)

    k_sampling = RateCoefficient(mbsamplingRate, 's^-1')

    constant_species_list = []

    for spec in constantSpecies:
        constant_species_list.append(species_dict[spec])

    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.items():
            termination.append(TerminationConversion(species_dict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if len(termination) == 0:
        raise InputError(
            'No termination conditions specified for reaction system #{0}.'.format(len(rmg.reaction_systems) + 2))

    sensitive_species = []
    if sensitivity:
        if isinstance(sensitivity, str):
            sensitivity = [sensitivity]
        for spec in sensitivity:
            sensitive_species.append(species_dict[spec])
    system = MBSampledReactor(T, P, initialMoleFractions, k_sampling, constant_species_list, termination,
                              sensitive_species, sensitivityThreshold)
    rmg.reaction_systems.append(system)


def simulator(atol, rtol, sens_atol=1e-6, sens_rtol=1e-4):
    rmg.simulator_settings_list.append(SimulatorSettings(atol, rtol, sens_atol, sens_rtol))


def solvation(solvent):
    # If solvation module in input file, set the RMG solvent variable
    if not isinstance(solvent, str):
        raise InputError("solvent should be a string like 'water'")
    rmg.solvent = solvent


def model(toleranceMoveToCore=None, toleranceMoveEdgeReactionToCore=np.inf, toleranceKeepInEdge=0.0,
          toleranceInterruptSimulation=1.0,
          toleranceMoveEdgeReactionToSurface=np.inf, toleranceMoveSurfaceSpeciesToCore=np.inf,
          toleranceMoveSurfaceReactionToCore=np.inf,
          toleranceMoveEdgeReactionToSurfaceInterrupt=None,
          toleranceMoveEdgeReactionToCoreInterrupt=None, maximumEdgeSpecies=1000000, minCoreSizeForPrune=50,
          minSpeciesExistIterationsForPrune=2, filterReactions=False, filterThreshold=1e8,
          ignoreOverallFluxCriterion=False,
          maxNumSpecies=None, maxNumObjsPerIter=1, terminateAtMaxObjects=False,
          toleranceThermoKeepSpeciesInEdge=np.inf, dynamicsTimeScale=(0.0, 'sec'),
          toleranceBranchReactionToCore=0.0, branchingIndex=0.5, branchingRatioMax=1.0,
          toleranceTransitoryDict={}, transitoryStepPeriod=20):
    """
    How to generate the model. `toleranceMoveToCore` must be specified. 
    toleranceMoveReactionToCore and toleranceReactionInterruptSimulation refers to an additional criterion for forcing an edge reaction to be included in the core
    by default this criterion is turned off
    Other parameters are optional and control the pruning.
    ignoreOverallFluxCriterion=True will cause the toleranceMoveToCore to be only applied
    to the pressure dependent network expansion and not movement of species from edge to core
    """
    if toleranceMoveToCore is None:
        raise InputError("You must provide a toleranceMoveToCore value. It should be less than or equal to "
                         "toleranceInterruptSimulation which is currently {0}".format(toleranceInterruptSimulation))
    if toleranceMoveToCore > toleranceInterruptSimulation:
        raise InputError("toleranceMoveToCore must be less than or equal to toleranceInterruptSimulation, which is "
                         "currently {0}".format(toleranceInterruptSimulation))

    rmg.model_settings_list.append(
        ModelSettings(
            tol_move_to_core=toleranceMoveToCore,
            tol_move_edge_rxn_to_core=toleranceMoveEdgeReactionToCore,
            tol_keep_in_edge=toleranceKeepInEdge,
            tol_interrupt_simulation=toleranceInterruptSimulation,
            tol_move_edge_rxn_to_surface=toleranceMoveEdgeReactionToSurface,
            tol_move_surface_spc_to_core=toleranceMoveSurfaceSpeciesToCore,
            tol_move_surface_rxn_to_core=toleranceMoveSurfaceReactionToCore,
            tol_move_edge_rxn_to_surface_interrupt=toleranceMoveEdgeReactionToSurfaceInterrupt,
            tol_move_edge_rxn_to_core_interrupt=toleranceMoveEdgeReactionToCoreInterrupt,
            maximum_edge_species=maximumEdgeSpecies,
            min_core_size_for_prune=minCoreSizeForPrune,
            min_species_exist_iterations_for_prune=minSpeciesExistIterationsForPrune,
            filter_reactions=filterReactions,
            filter_threshold=filterThreshold,
            ignore_overall_flux_criterion=ignoreOverallFluxCriterion,
            max_num_species=maxNumSpecies,
            max_num_objects_per_iter=maxNumObjsPerIter,
            terminate_at_max_objects=terminateAtMaxObjects,
            thermo_tol_keep_spc_in_edge=toleranceThermoKeepSpeciesInEdge,
            dynamics_time_scale=Quantity(dynamicsTimeScale),
            tol_branch_rxn_to_core=toleranceBranchReactionToCore,
            branching_index=branchingIndex,
            branching_ratio_max=branchingRatioMax,
            transitory_tol_dict=toleranceTransitoryDict,
            transitory_step_period=transitoryStepPeriod,
        )
    )


def quantum_mechanics(
        software,
        method,
        fileStore=None,
        scratchDirectory=None,
        onlyCyclics=False,
        maxRadicalNumber=0,
):
    from rmgpy.qm.main import QMCalculator
    rmg.quantum_mechanics = QMCalculator(
        software=software,
        method=method,
        fileStore=fileStore,
        scratchDirectory=scratchDirectory,
        onlyCyclics=onlyCyclics,
        maxRadicalNumber=maxRadicalNumber,
    )


def ml_estimator(thermo=True,
                 name='main',
                 minHeavyAtoms=1,
                 maxHeavyAtoms=None,
                 minCarbonAtoms=0,
                 maxCarbonAtoms=None,
                 minOxygenAtoms=0,
                 maxOxygenAtoms=None,
                 minNitrogenAtoms=0,
                 maxNitrogenAtoms=None,
                 onlyCyclics=False,
                 onlyHeterocyclics=False,
                 minCycleOverlap=0,
                 H298UncertaintyCutoff=(3.0, 'kcal/mol'),
                 S298UncertaintyCutoff=(2.0, 'cal/(mol*K)'),
                 CpUncertaintyCutoff=(2.0, 'cal/(mol*K)')):
    from rmgpy.ml.estimator import MLEstimator

    # Currently only support thermo
    if thermo:
        models_path = os.path.join(settings['database.directory'], 'thermo', 'ml', name)
        if not os.path.exists(models_path):
            raise InputError('Cannot find ML models folder {}'.format(models_path))
        hf298_path = os.path.join(models_path, 'hf298')
        s298_cp_path = os.path.join(models_path, 's298_cp')
        rmg.ml_estimator = MLEstimator(hf298_path, s298_cp_path)

        uncertainty_cutoffs = dict(
            H298=Quantity(*H298UncertaintyCutoff),
            S298=Quantity(*S298UncertaintyCutoff),
            Cp=Quantity(*CpUncertaintyCutoff)
        )
        rmg.ml_settings = dict(
            min_heavy_atoms=minHeavyAtoms,
            max_heavy_atoms=maxHeavyAtoms,
            min_carbon_atoms=minCarbonAtoms,
            max_carbon_atoms=maxCarbonAtoms,
            min_oxygen_atoms=minOxygenAtoms,
            max_oxygen_atoms=maxOxygenAtoms,
            min_nitrogen_atoms=minNitrogenAtoms,
            max_nitrogen_atoms=maxNitrogenAtoms,
            only_cyclics=onlyCyclics,
            only_heterocyclics=onlyHeterocyclics,
            min_cycle_overlap=minCycleOverlap,
            uncertainty_cutoffs=uncertainty_cutoffs,
        )

    # Shows warning when onlyCyclics is False and onlyHeterocyclics is True
    if minCycleOverlap > 0 and not onlyCyclics and not onlyHeterocyclics:
        logging.warning('"onlyCyclics" should be True when "minCycleOverlap" is greater than zero. '
                        'Machine learning estimator is restricted to only cyclic species thermo with the specified '
                        'minimum cycle overlap')
    elif minCycleOverlap > 0 and not onlyCyclics and onlyHeterocyclics:
        logging.warning('"onlyCyclics" should be True when "onlyHeterocyclics" is True and "minCycleOverlap" is '
                        'greater than zero. Machine learning estimator is restricted to only heterocyclic species '
                        'thermo with the specified minimum cycle overlap')
    elif onlyHeterocyclics and not onlyCyclics:
        logging.warning('"onlyCyclics" should be True when "onlyHeterocyclics" is True. '
                        'Machine learning estimator is restricted to only heterocyclic species thermo')


def pressure_dependence(
        method,
        temperatures,
        pressures,
        maximumGrainSize=0.0,
        minimumNumberOfGrains=0,
        interpolation=None,
        maximumAtoms=None,
):
    from arkane.pdep import PressureDependenceJob

    # Setting the pressureDependence attribute to non-None enables pressure dependence
    rmg.pressure_dependence = PressureDependenceJob(network=None)

    # Process method
    rmg.pressure_dependence.method = method

    # Process interpolation model
    if isinstance(interpolation, str):
        interpolation = (interpolation,)
    if interpolation[0].lower() not in ("chebyshev", "pdeparrhenius"):
        raise InputError("Interpolation model must be set to either 'Chebyshev' or 'PDepArrhenius'.")
    rmg.pressure_dependence.interpolation_model = interpolation

    # Process temperatures
    Tmin, Tmax, Tunits, Tcount = temperatures
    rmg.pressure_dependence.Tmin = Quantity(Tmin, Tunits)
    rmg.pressure_dependence.Tmax = Quantity(Tmax, Tunits)
    rmg.pressure_dependence.Tcount = Tcount
    rmg.pressure_dependence.generate_T_list()

    # Process pressures
    Pmin, Pmax, Punits, Pcount = pressures
    rmg.pressure_dependence.Pmin = Quantity(Pmin, Punits)
    rmg.pressure_dependence.Pmax = Quantity(Pmax, Punits)
    rmg.pressure_dependence.Pcount = Pcount
    rmg.pressure_dependence.generate_P_list()

    # Process grain size and count
    rmg.pressure_dependence.maximum_grain_size = Quantity(maximumGrainSize)
    rmg.pressure_dependence.minimum_grain_count = minimumNumberOfGrains

    # Process maximum atoms
    rmg.pressure_dependence.maximum_atoms = maximumAtoms

    rmg.pressure_dependence.active_j_rotor = True
    rmg.pressure_dependence.active_k_rotor = True
    rmg.pressure_dependence.rmgmode = True


def options(name='Seed', generateSeedEachIteration=True, saveSeedToDatabase=False, units='si', saveRestartPeriod=None,
            generateOutputHTML=False, generatePlots=False, saveSimulationProfiles=False, verboseComments=False,
            saveEdgeSpecies=False, keepIrreversible=False, trimolecularProductReversible=True, wallTime='00:00:00:00',
            saveSeedModulus=-1):
    if saveRestartPeriod:
        logging.warning("`saveRestartPeriod` flag was set in the input file, but this feature has been removed. Please "
                        "remove this line from the input file. This will throw an error after RMG-Py 3.1. For "
                        "restarting an RMG job see the documentation for restarting from a seed mechanism at "
                        "http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/input.html#restarting-from-a-seed-mechanism")

    rmg.name = name
    rmg.generate_seed_each_iteration = generateSeedEachIteration
    rmg.save_seed_to_database = saveSeedToDatabase
    rmg.units = units
    if generateOutputHTML:
        logging.warning('Generate Output HTML option was turned on. Note that this will slow down model generation.')
    rmg.generate_output_html = generateOutputHTML
    rmg.generate_plots = generatePlots
    rmg.save_simulation_profiles = saveSimulationProfiles
    rmg.verbose_comments = verboseComments
    if saveEdgeSpecies:
        logging.warning(
            'Edge species saving was turned on. This will slow down model generation for large simulations.')
    rmg.save_edge_species = saveEdgeSpecies
    rmg.keep_irreversible = keepIrreversible
    rmg.trimolecular_product_reversible = trimolecularProductReversible
    rmg.walltime = wallTime
    rmg.save_seed_modulus = saveSeedModulus


def generated_species_constraints(**kwargs):
    valid_constraints = [
        'allowed',
        'maximumCarbonAtoms',
        'maximumOxygenAtoms',
        'maximumNitrogenAtoms',
        'maximumSiliconAtoms',
        'maximumSulfurAtoms',
        'maximumSurfaceSites',
        'maximumHeavyAtoms',
        'maximumRadicalElectrons',
        'maximumSingletCarbenes',
        'maximumCarbeneRadicals',
        'allowSingletO2',
        'speciesCuttingThreshold',
    ]

    for key, value in kwargs.items():
        if key not in valid_constraints:
            raise InputError('Invalid generated species constraint {0!r}.'.format(key))

        rmg.species_constraints[key] = value


def thermo_central_database(host,
                            port,
                            username,
                            password,
                            application):
    from rmgpy.data.thermo import ThermoCentralDatabaseInterface
    rmg.thermo_central_database = ThermoCentralDatabaseInterface(host,
                                                                 port,
                                                                 username,
                                                                 password,
                                                                 application)


def uncertainty(localAnalysis=False, globalAnalysis=False, uncorrelated=True, correlated=True,
                localNumber=10, globalNumber=5, terminationTime=None,
                pceRunTime=1800, pceErrorTol=None, pceMaxEvals=None, logx=True):
    if not localAnalysis and globalAnalysis:
        logging.info('Enabling local uncertainty analysis as prerequisite for running global uncertainty analysis.')

    rmg.uncertainty = {
        'local': localAnalysis if not globalAnalysis else True,  # Must run local before global
        'global': globalAnalysis,
        'uncorrelated': uncorrelated,
        'correlated': correlated,
        'localnum': localNumber,
        'globalnum': globalNumber,
        'time': Quantity(terminationTime) if terminationTime else terminationTime,
        'pcetime': pceRunTime,
        'pcetol': pceErrorTol,
        'pceevals': pceMaxEvals,
        'logx': logx,
    }


def restart_from_seed(path=None, coreSeed=None, edgeSeed=None, filters=None, speciesMap=None):
    parent_dir = os.path.dirname(rmg.input_file)
    rmg.restart = True
    doc_link = 'http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/input.html#restarting-from-a-seed-mechanism.'

    if path:
        if any((coreSeed, edgeSeed, filters, speciesMap)):
            raise InputError('For restarting an RMG job from a seed mechanism, either the path to the RMG generated '
                             'seed mechanism should be given as `path`, or the path for each of the required files '
                             'should be explicitly given, but not both. Please take one approach or the other. For '
                             'further information see the RMG documentation on restarting from a seed mechanism at '
                             '{0}.'.format(doc_link))

        if not os.path.isabs(path):
            path = os.path.join(parent_dir, path)

        if not os.path.exists(path):
            raise ValueError('Unable to find the path to the restart seed folder. {0} does not exist'.format(path))

        # Try to find the paths for all of the required modules
        rmg.core_seed_path = os.path.join(path, 'seed')
        rmg.edge_seed_path = os.path.join(path, 'seed_edge')
        rmg.filters_path = os.path.join(path, 'filters', 'filters.h5')
        rmg.species_map_path = os.path.join(path, 'filters', 'species_map.yml')

    else:  # The user has specified each of the paths individually
        rmg.core_seed_path = coreSeed
        rmg.edge_seed_path = edgeSeed
        rmg.filters_path = filters
        rmg.species_map_path = speciesMap

    rmg_paths = [rmg.core_seed_path, rmg.edge_seed_path, rmg.filters_path, rmg.species_map_path]
    path_errors = [filePath for filePath in rmg_paths if not os.path.exists(filePath)]

    if path_errors:
        if path:
            raise InputError('Could not find one or more of the required files/directories for restarting from a seed '
                             'mechanism: {0}. Try specifying the file paths individually. See the RMG documentation '
                             'at {1} for more information'.format(path_errors, doc_link))
        else:
            raise InputError('Could not find one or more of the required files/directories for restarting from a seed '
                             'mechanism: {0}. See the RMG documentation at {1} for more information'.format(path_errors,
                                                                                                            doc_link))

def liquid_volumetric_mass_transfer_coefficient_power_law(prefactor=(0,"1/s"), diffusionCoefficientPower=0, solventViscosityPower=0, solventDensityPower=0):

    rmg.liquid_volumetric_mass_transfer_coefficient_power_law = liquidVolumetricMassTransferCoefficientPowerLaw(prefactor=Quantity(prefactor).value_si,
                                                                                                                diffusion_coefficient_power=diffusionCoefficientPower,
                                                                                                                solvent_viscosity_power=solventViscosityPower,
                                                                                                                solvent_density_power=solventDensityPower)

################################################################################

def set_global_rmg(rmg0):
    """
    sets the global variable rmg to rmg0. This is used to allow for unittesting
    of above methods
    """
    global rmg
    rmg = rmg0


def read_input_file(path, rmg0):
    """
    Read an RMG input file at `path` on disk into the :class:`RMG` object 
    `rmg`.
    """
    global rmg, species_dict, mol_to_frag

    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise

    logging.info('Reading input file "{0}"...'.format(full_path))
    logging.info(f.read())
    f.seek(0)  # return to beginning of file

    set_global_rmg(rmg0)
    rmg.reaction_model = CoreEdgeReactionModel()
    rmg.initial_species = []
    rmg.reaction_systems = []
    species_dict = {}
    mol_to_frag = {}

    global_context = {'__builtins__': None}
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'database': database,
        'catalystProperties': catalyst_properties,
        'species': species,
        'forbidden': forbidden,
        'SMARTS': smarts,
        'fragment_adj': fragment_adj,
        'fragment_SMILES': fragment_smiles,
        'SMILES': smiles,
        'InChI': inchi,
        'adjacencyList': adjacency_list,
        'adjacencyListGroup': adjacency_list_group,
        'react': react,
        'simpleReactor': simple_reactor,
        'constantVIdealGasReactor' : constant_V_ideal_gas_reactor,
        'constantTPIdealGasReactor' : constant_TP_ideal_gas_reactor,
        'liquidSurfaceReactor' : liquid_cat_reactor,
        'constantTVLiquidReactor': constant_T_V_liquid_reactor,
        'liquidReactor': liquid_reactor,
        'surfaceReactor': surface_reactor,
        'mbsampledReactor': mb_sampled_reactor,
        'simulator': simulator,
        'solvation': solvation,
        'liquidVolumetricMassTransferCoefficientPowerLaw': liquid_volumetric_mass_transfer_coefficient_power_law,
        'model': model,
        'quantumMechanics': quantum_mechanics,
        'mlEstimator': ml_estimator,
        'pressureDependence': pressure_dependence,
        'options': options,
        'generatedSpeciesConstraints': generated_species_constraints,
        'thermoCentralDatabase': thermo_central_database,
        'uncertainty': uncertainty,
        'restartFromSeed': restart_from_seed,
    }

    thermo_libraries = rmg0.thermo_libraries if isinstance(rmg0.thermo_libraries, list) else None
    reaction_libraries = rmg0.reaction_libraries if isinstance(rmg0.reaction_libraries, list) else None

    try:
        exec(f.read(), global_context, local_context)
    except (NameError, TypeError, SyntaxError) as e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    if thermo_libraries is not None:
        rmg0.thermo_libraries.extend(thermo_libraries)
    if reaction_libraries is not None:
        rmg0.reaction_libraries.extend(reaction_libraries)

    rmg.species_constraints['explicitlyAllowedMolecules'] = []

    # convert keys from species names into species objects.
    for reaction_system in rmg.reaction_systems:
        if not isinstance(reaction_system, Reactor):
            reaction_system.convert_initial_keys_to_species_objects(species_dict)

    if rmg.quantum_mechanics:
        rmg.quantum_mechanics.set_default_output_directory(rmg.output_directory)
        rmg.quantum_mechanics.initialize()

    logging.info('')


################################################################################

def read_thermo_input_file(path, rmg0):
    """
    Read an thermo estimation input file at `path` on disk into the :class:`RMG` object 
    `rmg`.
    """

    global rmg, species_dict

    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise

    logging.info('Reading input file "{0}"...'.format(full_path))

    rmg = rmg0
    rmg.reaction_model = CoreEdgeReactionModel()
    rmg.initial_species = []
    rmg.reaction_systems = []
    species_dict = {}

    global_context = {'__builtins__': None}
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'database': database,
        'catalystProperties': catalyst_properties,
        'species': species,
        'SMARTS': smarts,
        'SMILES': smiles,
        'InChI': inchi,
        'solvation': solvation,
        'adjacencyList': adjacency_list,
        'quantumMechanics': quantum_mechanics,
        'mlEstimator': ml_estimator,
    }

    try:
        exec(f.read(), global_context, local_context)
    except (NameError, TypeError, SyntaxError) as e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    if rmg.quantum_mechanics:
        rmg.quantum_mechanics.set_default_output_directory(rmg.output_directory)
        rmg.quantum_mechanics.initialize()

    logging.info('')


################################################################################

def save_input_file(path, rmg):
    """
    Save an RMG input file at `path` on disk from the :class:`RMG` object 
    `rmg`.
    """

    f = open(path, 'w')

    # Databases
    f.write('database(\n')
    # f.write('    "{0}",\n'.format(rmg.database_directory))
    f.write('    thermoLibraries = {0!r},\n'.format(rmg.thermo_libraries))
    f.write('    reactionLibraries = {0!r},\n'.format(rmg.reaction_libraries))
    f.write('    seedMechanisms = {0!r},\n'.format(rmg.seed_mechanisms))
    f.write('    kinetics_depositories = {0!r},\n'.format(rmg.kinetics_depositories))
    f.write('    kineticsFamilies = {0!r},\n'.format(rmg.kinetics_families))
    f.write('    kineticsEstimator = {0!r},\n'.format(rmg.kinetics_estimator))
    f.write(')\n\n')

    if rmg.surfaceSiteDenisty or rmg.binding_energies:
        f.write('catalystProperties(\n')
        if rmg.surfaceSiteDenisty:
            f.write('    surface_site_density = {0!r},'.format(rmg.surface_site_density))
        if rmg.binding_energies:
            f.write('    binding_energies = {0!r},'.format(rmg.binding_energies))
        f.write(')\n\n')

    # Species
    for spcs in rmg.initial_species:
        f.write('species(\n')
        f.write('    label = "{0}",\n'.format(spcs.label))
        f.write('    reactive = {0},\n'.format(spcs.reactive))
        f.write('    structure = adjacencyList(\n')
        f.write('"""\n')
        f.write(spcs.molecule[0].to_adjacency_list())
        f.write('"""),\n')
        f.write(')\n\n')

    # Reaction systems
    for system in rmg.reaction_systems:
        if rmg.solvent:
            f.write('liquidReactor(\n')
            f.write('    temperature = ({0:g},"{1!s}"),\n'.format(system.T.value, system.T.units))
            f.write('    initialConcentrations={\n')
            for spcs, conc in system.initial_concentrations.items():
                f.write('        "{0!s}": ({1:g},"{2!s}"),\n'.format(spcs.label, conc.value, conc.units))
        else:
            f.write('simpleReactor(\n')
            f.write('    temperature = ({0:g},"{1!s}"),\n'.format(system.T.value, system.T.units))
            # Convert the pressure from SI pascal units to bar here
            # Do something more fancy later for converting to user's desired units for both T and P..
            f.write('    pressure = ({0:g},"{1!s}"),\n'.format(system.P.value, system.P.units))
            f.write('    initialMoleFractions={\n')
            for spcs, molfrac in system.initial_mole_fractions.items():
                f.write('        "{0!s}": {1:g},\n'.format(spcs.label, molfrac))
        f.write('    },\n')

        # Termination criteria
        conversions = ''
        for term in system.termination:
            if isinstance(term, TerminationTime):
                f.write('    terminationTime = ({0:g},"{1!s}"),\n'.format(term.time.value, term.time.units))

            else:
                conversions += '        "{0:s}": {1:g},\n'.format(term.species.label, term.conversion)
        if conversions:
            f.write('    terminationConversion = {\n')
            f.write(conversions)
            f.write('    },\n')

        # Sensitivity analysis
        if system.sensitive_species:
            sensitivity = []
            for item in system.sensitive_species:
                sensitivity.append(item.label)
            f.write('    sensitivity = {0},\n'.format(sensitivity))
            f.write('    sensitivityThreshold = {0},\n'.format(system.sensitivity_threshold))

        f.write(')\n\n')

    if rmg.solvent:
        f.write("solvation(\n    solvent = '{0!s}'\n)\n\n".format(rmg.solvent))

    # Simulator tolerances
    f.write('simulator(\n')
    f.write('    atol = {0:g},\n'.format(rmg.simulator_settings_list[0].atol))
    f.write('    rtol = {0:g},\n'.format(rmg.simulator_settings_list[0].rtol))
    f.write('    sens_atol = {0:g},\n'.format(rmg.simulator_settings_list[0].sens_atol))
    f.write('    sens_rtol = {0:g},\n'.format(rmg.simulator_settings_list[0].sens_rtol))
    f.write(')\n\n')

    # Model
    f.write('model(\n')
    f.write('    toleranceMoveToCore = {0:g},\n'.format(rmg.model_settings_list[0].tol_move_to_core))
    f.write('    toleranceKeepInEdge = {0:g},\n'.format(rmg.model_settings_list[0].tol_keep_in_edge))
    f.write('    toleranceInterruptSimulation = {0:g},\n'.format(rmg.model_settings_list[0].tol_interrupt_simulation))
    f.write('    maximumEdgeSpecies = {0:d},\n'.format(rmg.model_settings_list[0].maximum_edge_species))
    f.write('    minCoreSizeForPrune = {0:d},\n'.format(rmg.model_settings_list[0].min_core_size_for_prune))
    f.write('    minSpeciesExistIterationsForPrune = {0:d},\n'.format(rmg.model_settings_list[0].min_species_exist_iterations_for_prune))
    f.write('    filterReactions = {0:d},\n'.format(rmg.model_settings_list[0].filter_reactions))
    f.write('    filterThreshold = {0:g},\n'.format(rmg.model_settings_list[0].filter_threshold))
    f.write(')\n\n')

    # Pressure Dependence
    if rmg.pressure_dependence:
        f.write('pressureDependence(\n')
        f.write('    method = {0!r},\n'.format(rmg.pressure_dependence.method))
        f.write('    maximumGrainSize = ({0:g},"{1!s}"),\n'.format(rmg.pressure_dependence.grain_size.value,
                                                                   rmg.pressure_dependence.grain_size.units))
        f.write('    minimumNumberOfGrains = {0},\n'.format(rmg.pressure_dependence.grain_count))
        f.write('    temperatures = ({0:g},{1:g},"{2!s}",{3:d}),\n'.format(
            rmg.pressure_dependence.Tmin.value,
            rmg.pressure_dependence.Tmax.value,
            rmg.pressure_dependence.Tmax.units,
            rmg.pressure_dependence.Tcount,
        ))
        f.write('    pressures = ({0:g},{1:g},"{2!s}",{3:d}),\n'.format(
            rmg.pressure_dependence.Pmin.value,
            rmg.pressure_dependence.Pmax.value,
            rmg.pressure_dependence.Pmax.units,
            rmg.pressure_dependence.Pcount,
        ))
        f.write('    interpolation = {0},\n'.format(rmg.pressure_dependence.interpolation_model))
        f.write('    maximumAtoms = {0}, \n'.format(rmg.pressure_dependence.maximum_atoms))
        f.write(')\n\n')

    # Quantum Mechanics
    if rmg.quantum_mechanics:
        f.write('quantumMechanics(\n')
        f.write('    software = {0!r},\n'.format(rmg.quantum_mechanics.settings.software))
        f.write('    method = {0!r},\n'.format(rmg.quantum_mechanics.settings.method))
        # Split paths created by QMSettings
        if rmg.quantum_mechanics.settings.fileStore:
            f.write('    fileStore = {0!r},\n'.format(os.path.split(rmg.quantum_mechanics.settings.fileStore)[0]))
        else:
            f.write('    fileStore = None,\n')
        if rmg.quantum_mechanics.settings.scratchDirectory:
            f.write('    scratchDirectory = {0!r},\n'.format(
                os.path.split(rmg.quantum_mechanics.settings.scratchDirectory)[0]))
        else:
            f.write('    scratchDirectory = None,\n')
        f.write('    onlyCyclics = {0},\n'.format(rmg.quantum_mechanics.settings.onlyCyclics))
        f.write('    maxRadicalNumber = {0},\n'.format(rmg.quantum_mechanics.settings.maxRadicalNumber))
        f.write(')\n\n')

    # Species Constraints
    if rmg.species_constraints:
        f.write('generatedSpeciesConstraints(\n')
        for constraint, value in sorted(list(rmg.species_constraints.items()), key=lambda constraint: constraint[0]):
            if value is not None:
                f.write('    {0} = {1},\n'.format(constraint, value))
        f.write(')\n\n')

    # Options
    f.write('options(\n')
    f.write('    units = "{0}",\n'.format(rmg.units))
    f.write('    generateOutputHTML = {0},\n'.format(rmg.generate_output_html))
    f.write('    generatePlots = {0},\n'.format(rmg.generate_plots))
    f.write('    saveSimulationProfiles = {0},\n'.format(rmg.save_simulation_profiles))
    f.write('    saveEdgeSpecies = {0},\n'.format(rmg.save_edge_species))
    f.write('    keepIrreversible = {0},\n'.format(rmg.keep_irreversible))
    f.write('    trimolecularProductReversible = {0},\n'.format(rmg.trimolecular_product_reversible))
    f.write('    verboseComments = {0},\n'.format(rmg.verbose_comments))
    f.write('    wallTime = {0},\n'.format(rmg.walltime))
    f.write(')\n\n')

    f.close()


def get_input(name):
    """
    Returns the RMG input object that corresponds
    to the parameter name.

    First, the module level is queried. If this variable
    is empty, the broadcasted variables are queried.
    """
    global rmg

    if rmg:
        if name == 'species_constraints':
            return rmg.species_constraints
        elif name == 'quantum_mechanics':
            return rmg.quantum_mechanics
        elif name == 'ml_estimator':
            return rmg.ml_estimator, rmg.ml_settings
        elif name == 'thermo_central_database':
            return rmg.thermo_central_database
        else:
            raise Exception('Unrecognized keyword: {}'.format(name))

    raise Exception('Could not get variable with name: {}'.format(name))
