#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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
import warnings
import quantities
import os
import numpy
from copy import deepcopy

from rmgpy import settings

from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.solver.base import TerminationTime, TerminationConversion, TerminationRateRatio
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from model import CoreEdgeReactionModel

from rmgpy.scoop_framework.util import broadcast, get
from rmgpy.exceptions import InputError
################################################################################

rmg = None
speciesDict = {}

def database(
             thermoLibraries = None,
             transportLibraries = None,
             reactionLibraries = None,
             frequenciesLibraries = None,
             seedMechanisms = None,
             kineticsFamilies = 'default',
             kineticsDepositories = 'default',
             kineticsEstimator = 'rate rules',
             ):
    # This function just stores the information about the database to be loaded
    # We don't actually load the database until after we're finished reading
    # the input file
    if isinstance(thermoLibraries, str): thermoLibraries = [thermoLibraries]
    if isinstance(transportLibraries, str): transportLibraries = [transportLibraries]
    if isinstance(reactionLibraries, str): reactionLibraries = [reactionLibraries]
    if isinstance(seedMechanisms, str): seedMechanisms = [seedMechanisms]
    if isinstance(frequenciesLibraries, str): frequenciesLibraries = [frequenciesLibraries]
    rmg.databaseDirectory = settings['database.directory']
    rmg.thermoLibraries = thermoLibraries or []
    rmg.transportLibraries = transportLibraries
    # Modify reactionLibraries if the user didn't specify tuple input
    if reactionLibraries:
        index = 0
        while index < len(reactionLibraries):
            if isinstance(reactionLibraries[index],tuple):
                pass
            elif isinstance(reactionLibraries[index],str):
                reactionLibraries[index] = (reactionLibraries[index], False)
            else:
                raise TypeError('reaction libraries must be input as tuples or strings')
            index += 1
    rmg.reactionLibraries = reactionLibraries or []
    rmg.seedMechanisms = seedMechanisms or []
    rmg.statmechLibraries = frequenciesLibraries or []
    rmg.kineticsEstimator = kineticsEstimator
    if kineticsDepositories == 'default':
        rmg.kineticsDepositories = ['training']
    elif kineticsDepositories == 'all':
        rmg.kineticsDepositories = None
    else:
        if not isinstance(kineticsDepositories,list):
            raise InputError("kineticsDepositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe'].")
        rmg.kineticsDepositories = kineticsDepositories

    if kineticsFamilies in ('default', 'all', 'none'):
        rmg.kineticsFamilies = kineticsFamilies
    else:
        if not isinstance(kineticsFamilies,list):
            raise InputError("kineticsFamilies should be either 'default', 'all', 'none', or a list of names eg. ['H_Abstraction','R_Recombination'] or ['!Intra_Disproportionation'].")
        rmg.kineticsFamilies = kineticsFamilies

def species(label, structure, reactive=True):
    logging.debug('Found {0} species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
    
    if '+' in label:
        raise InputError('species {0} label cannot include a + sign'.format(label))
    
    spec, isNew = rmg.reactionModel.makeNewSpecies(structure, label=label, reactive=reactive)
    if not isNew:
        raise InputError("Species {0} is a duplicate of {1}. Species in input file must be unique".format(label,spec.label))
    # Force RMG to add the species to edge first, prior to where it is added to the core, in case it is found in 
    # any reaction libraries along the way
    rmg.reactionModel.addSpeciesToEdge(spec)
    rmg.initialSpecies.append(spec)
    speciesDict[label] = spec
    
def SMARTS(string):
    return Molecule().fromSMARTS(string)

def SMILES(string):
    return Molecule().fromSMILES(string)

def InChI(string):
    return Molecule().fromInChI(string)

def adjacencyList(string):
    return Molecule().fromAdjacencyList(string)

# Reaction systems
def simpleReactor(temperature,
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
                  ):
    logging.debug('Found SimpleReactor reaction system')
    
    for key,value in initialMoleFractions.iteritems():
        if not isinstance(value,list):
            initialMoleFractions[key] = float(value)
            if value < 0:
                raise InputError('Initial mole fractions cannot be negative.')
        else:
            if len(value) != 2:
                raise InputError("Initial mole fraction values must either be a number or a list with 2 entries")
            initialMoleFractions[key] = [float(value[0]),float(value[1])]
            if value[0] < 0 or value[1] < 0:
                raise InputError('Initial mole fractions cannot be negative.')
            elif value[1] < value[0]:
                raise InputError('Initial mole fraction range out of order: {0}'.format(key))
    
    if not isinstance(temperature,list):
        T = Quantity(temperature)
    else:
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]
        
    if not isinstance(pressure,list):
        P = Quantity(pressure)
    else:
        if len(pressure) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list with 2 entries of the same format')
        P = [Quantity(p) for p in pressure]
        
        
    if not isinstance(temperature,list) and not isinstance(pressure,list) and all([not isinstance(x,list) for x in initialMoleFractions.values()]):
        nSims=1
    
    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.iteritems():
            termination.append(TerminationConversion(speciesDict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reactionSystems)+2))
    
    sensitiveSpecies = []
    if sensitivity:
        if sensitivity != 'all':
            if isinstance(sensitivity, str): sensitivity = [sensitivity]
            for spec in sensitivity:
                sensitiveSpecies.append(speciesDict[spec])

        else:
            sensitiveSpecies.append('all')
    
    if not isinstance(T,list):
        sensitivityTemperature = T
    if not isinstance(P,list):
        sensitivityPressure = P
    if not any([isinstance(x,list) for x in initialMoleFractions.itervalues()]):
        sensitivityMoleFractions = deepcopy(initialMoleFractions)
    if sensitivityMoleFractions is None or sensitivityTemperature is None or sensitivityPressure is None:
        sensConditions = None
    else:
        sensConditions = sensitivityMoleFractions
        sensConditions['T'] = Quantity(sensitivityTemperature).value_si
        sensConditions['P'] = Quantity(sensitivityPressure).value_si

    
    system = SimpleReactor(T, P, initialMoleFractions, nSims, termination, sensitiveSpecies, sensitivityThreshold,sensConditions)
    rmg.reactionSystems.append(system)
    
    assert balanceSpecies is None or isinstance(balanceSpecies,str), 'balanceSpecies should be the string corresponding to a single species'
    rmg.balanceSpecies = balanceSpecies
    if balanceSpecies: #check that the balanceSpecies can't be taken to zero
        total = 0.0
        for key,item in initialMoleFractions.iteritems():
            if key == balanceSpecies:
                assert not isinstance(item,list), 'balanceSpecies must not have a defined range'
                xbspcs = item
            if isinstance(item,list):
                total += item[1]-item[0]

        if total > xbspcs:
            raise ValueError('The sum of the differences in the ranged mole fractions is greater than the mole fraction of the balance species, this would require the balanceSpecies mole fraction to be negative in some cases which is not allowed, either reduce the maximum mole fractions or dont use balanceSpecies')

# Reaction systems
def liquidReactor(temperature,
                  initialConcentrations,
                  terminationConversion=None,
                  nSims = 4,
                  terminationTime=None,
                  terminationRateRatio=None,
                  sensitivity=None,
                  sensitivityThreshold=1e-3,
                  sensitivityTemperature=None,
                  sensitivityConcentrations=None,
                  constantSpecies=None):
    
    logging.debug('Found LiquidReactor reaction system')
    
    if not isinstance(temperature,list):
        T = Quantity(temperature)
    else:
        if len(temperature) != 2:
            raise InputError('Temperature and pressure ranges can either be in the form of (number,units) or a list with 2 entries of the same format')
        T = [Quantity(t) for t in temperature]
    
    for spec,conc in initialConcentrations.iteritems():
        if not isinstance(conc,list):
            concentration = Quantity(conc)
            # check the dimensions are ok
            # convert to mol/m^3 (or something numerically nice? or must it be SI)
            initialConcentrations[spec] = concentration.value_si
        else:
            if len(conc) != 2:
                raise InputError("Concentration values must either be in the form of (number,units) or a list with 2 entries of the same format")
            initialConcentrations[spec] = [Quantity(conc[0]),Quantity(conc[1])]
            
    if not isinstance(temperature,list) and all([not isinstance(x,list) for x in initialConcentrations.itervalues()]):
        nSims=1
        
    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.iteritems():
            termination.append(TerminationConversion(speciesDict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if terminationRateRatio is not None:
        termination.append(TerminationRateRatio(terminationRateRatio))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reactionSystems)+2))
    
    sensitiveSpecies = []
    if sensitivity:
        for spec in sensitivity:
            sensitiveSpecies.append(speciesDict[spec])
    
    ##chatelak: check the constant species exist
    if constantSpecies is not None:
        logging.debug('  Generation with constant species:')
        for constantSpecie in constantSpecies:
            logging.debug("  {0}".format(constantSpecie))
            if not speciesDict.has_key(constantSpecie):
                raise InputError('Species {0} not found in the input file'.format(constantSpecie))
    
    if not isinstance(T,list):
        sensitivityTemperature = T
    if not any([isinstance(x,list) for x in initialConcentrations.itervalues()]):
        sensitivityConcentrations = initialConcentrations
    if sensitivityConcentrations is None or sensitivityTemperature is None:
        sensConditions = None
    else:
        sensConditions = sensitivityConcentrations
        sensConditions['T'] = Quantity(sensitivityTemperature).value_si
        
    system = LiquidReactor(T, initialConcentrations, nSims, termination, sensitiveSpecies, sensitivityThreshold, sensConditions, constantSpecies)
    rmg.reactionSystems.append(system)
    
def simulator(atol, rtol, sens_atol=1e-6, sens_rtol=1e-4):
    rmg.simulatorSettingsList.append(SimulatorSettings(atol, rtol, sens_atol, sens_rtol))
    
def solvation(solvent):
    # If solvation module in input file, set the RMG solvent variable
    if not isinstance(solvent,str):
        raise InputError("solvent should be a string like 'water'")
    rmg.solvent = solvent

def model(toleranceMoveToCore=None, toleranceMoveEdgeReactionToCore=numpy.inf,toleranceKeepInEdge=0.0, toleranceInterruptSimulation=1.0, 
          toleranceMoveEdgeReactionToSurface=numpy.inf, toleranceMoveSurfaceSpeciesToCore=numpy.inf, toleranceMoveSurfaceReactionToCore=numpy.inf,
          toleranceMoveEdgeReactionToSurfaceInterrupt=None,
          toleranceMoveEdgeReactionToCoreInterrupt=None, maximumEdgeSpecies=1000000, minCoreSizeForPrune=50, 
          minSpeciesExistIterationsForPrune=2, filterReactions=False, filterThreshold=1e8, ignoreOverallFluxCriterion=False,
          maxNumSpecies=None,maxNumObjsPerIter=1,terminateAtMaxObjects=False,toleranceThermoKeepSpeciesInEdge=numpy.inf,dynamicsTimeScale=(0.0,'sec'),
          toleranceBranchReactionToCore=0.0, branchingIndex=0.5, branchingRatioMax=1.0):
    """
    How to generate the model. `toleranceMoveToCore` must be specified. 
    toleranceMoveReactionToCore and toleranceReactionInterruptSimulation refers to an additional criterion for forcing an edge reaction to be included in the core
    by default this criterion is turned off
    Other parameters are optional and control the pruning.
    ignoreOverallFluxCriterion=True will cause the toleranceMoveToCore to be only applied
    to the pressure dependent network expansion and not movement of species from edge to core
    """
    if toleranceMoveToCore is None:
        raise InputError("You must provide a toleranceMoveToCore value. It should be less than or equal to toleranceInterruptSimulation which is currently {0}".format(toleranceInterruptSimulation))
    if toleranceMoveToCore > toleranceInterruptSimulation:
        raise InputError("toleranceMoveToCore must be less than or equal to toleranceInterruptSimulation, which is currently {0}".format(toleranceInterruptSimulation))
    
    rmg.modelSettingsList.append(
        ModelSettings(
            toleranceMoveToCore=toleranceMoveToCore,
            toleranceMoveEdgeReactionToCore=toleranceMoveEdgeReactionToCore,
            toleranceKeepInEdge=toleranceKeepInEdge,
            toleranceInterruptSimulation=toleranceInterruptSimulation,
            toleranceMoveEdgeReactionToSurface=toleranceMoveEdgeReactionToSurface,
            toleranceMoveSurfaceSpeciesToCore=toleranceMoveSurfaceSpeciesToCore,
            toleranceMoveSurfaceReactionToCore=toleranceMoveSurfaceReactionToCore,
            toleranceMoveEdgeReactionToSurfaceInterrupt=toleranceMoveEdgeReactionToSurfaceInterrupt,
            toleranceMoveEdgeReactionToCoreInterrupt=toleranceMoveEdgeReactionToCoreInterrupt,
            maximumEdgeSpecies=maximumEdgeSpecies,
            minCoreSizeForPrune=minCoreSizeForPrune,
            minSpeciesExistIterationsForPrune=minSpeciesExistIterationsForPrune,
            filterReactions=filterReactions,
            filterThreshold=filterThreshold,
            ignoreOverallFluxCriterion=ignoreOverallFluxCriterion,
            maxNumSpecies=maxNumSpecies,
            maxNumObjsPerIter=maxNumObjsPerIter,
            terminateAtMaxObjects=terminateAtMaxObjects,
            toleranceThermoKeepSpeciesInEdge=toleranceThermoKeepSpeciesInEdge,
            dynamicsTimeScale=Quantity(dynamicsTimeScale),
            toleranceBranchReactionToCore=toleranceBranchReactionToCore,
            branchingIndex=branchingIndex,
            branchingRatioMax=branchingRatioMax,
        )
    )
    
def quantumMechanics(
                    software,
                    method,
                    fileStore = None,
                    scratchDirectory = None,
                    onlyCyclics = False,
                    maxRadicalNumber = 0,
                    ):
    from rmgpy.qm.main import QMCalculator
    rmg.quantumMechanics = QMCalculator(software = software,
                                        method = method,
                                        fileStore = fileStore,
                                        scratchDirectory = scratchDirectory,
                                        onlyCyclics = onlyCyclics,
                                        maxRadicalNumber = maxRadicalNumber,
                                        )

def mlEstimator(thermo=True,
                name='main',
                minHeavyAtoms=1,
                maxHeavyAtoms=None,
                H298UncertaintyCutoff=(3.0, 'kcal/mol'),
                S298UncertaintyCutoff=(2.0, 'cal/(mol*K)'),
                CpUncertaintyCutoff=(2.0, 'cal/(mol*K)')):
    from rmgpy.ml.estimator import MLEstimator

    # Currently only support thermo
    if thermo:
        models_path = os.path.join(settings['database.directory'], 'thermo', 'ml', name)
        if not os.path.exists(models_path):
            raise InputError('Cannot find ML models folder {}'.format(models_path))
        H298_path = os.path.join(models_path, 'H298')
        S298_path = os.path.join(models_path, 'S298')
        Cp_path = os.path.join(models_path, 'Cp')
        rmg.ml_estimator = MLEstimator(H298_path, S298_path, Cp_path)

        uncertainty_cutoffs = dict(
            H298=Quantity(*H298UncertaintyCutoff),
            S298=Quantity(*S298UncertaintyCutoff),
            Cp=Quantity(*CpUncertaintyCutoff)
        )
        rmg.ml_settings = dict(
            min_heavy_atoms=minHeavyAtoms,
            max_heavy_atoms=maxHeavyAtoms,
            uncertainty_cutoffs=uncertainty_cutoffs,
        )
                    

def pressureDependence(
                       method,
                       temperatures,
                       pressures,
                       maximumGrainSize = 0.0,
                       minimumNumberOfGrains = 0,
                       interpolation = None,
                       maximumAtoms=None,
                       ):

    from rmgpy.cantherm.pdep import PressureDependenceJob
    
    # Setting the pressureDependence attribute to non-None enables pressure dependence
    rmg.pressureDependence = PressureDependenceJob(network=None)
    
    # Process method
    rmg.pressureDependence.method = method
    
    # Process interpolation model
    if isinstance(interpolation, str):
        interpolation = (interpolation,)
    if interpolation[0].lower() not in ("chebyshev","pdeparrhenius"):
        raise InputError("Interpolation model must be set to either 'Chebyshev' or 'PDepArrhenius'.")
    rmg.pressureDependence.interpolationModel = interpolation

    # Process temperatures
    Tmin, Tmax, Tunits, Tcount = temperatures
    rmg.pressureDependence.Tmin = Quantity(Tmin, Tunits)
    rmg.pressureDependence.Tmax = Quantity(Tmax, Tunits)
    rmg.pressureDependence.Tcount = Tcount
    rmg.pressureDependence.generateTemperatureList()
    
    # Process pressures
    Pmin, Pmax, Punits, Pcount = pressures
    rmg.pressureDependence.Pmin = Quantity(Pmin, Punits)
    rmg.pressureDependence.Pmax = Quantity(Pmax, Punits)
    rmg.pressureDependence.Pcount = Pcount
    rmg.pressureDependence.generatePressureList()
    
    # Process grain size and count
    rmg.pressureDependence.maximumGrainSize = Quantity(maximumGrainSize)
    rmg.pressureDependence.minimumGrainCount = minimumNumberOfGrains
    
    # Process maximum atoms
    rmg.pressureDependence.maximumAtoms = maximumAtoms
    
    rmg.pressureDependence.activeJRotor = True
    rmg.pressureDependence.activeKRotor = True
    rmg.pressureDependence.rmgmode = True

def options(name='Seed', generateSeedEachIteration=False, saveSeedToDatabase=False, units='si', saveRestartPeriod=None, 
            generateOutputHTML=False, generatePlots=False, saveSimulationProfiles=False, verboseComments=False, 
            saveEdgeSpecies=False, keepIrreversible=False, trimolecularProductReversible=True, wallTime='00:00:00:00'):
    rmg.name = name
    rmg.generateSeedEachIteration=generateSeedEachIteration
    rmg.saveSeedToDatabase=saveSeedToDatabase
    rmg.units = units
    rmg.saveRestartPeriod = Quantity(saveRestartPeriod) if saveRestartPeriod else None
    if generateOutputHTML:
        logging.warning('Generate Output HTML option was turned on. Note that this will slow down model generation.')
    rmg.generateOutputHTML = generateOutputHTML 
    rmg.generatePlots = generatePlots
    rmg.saveSimulationProfiles = saveSimulationProfiles
    rmg.verboseComments = verboseComments
    if saveEdgeSpecies:
        logging.warning('Edge species saving was turned on. This will slow down model generation for large simulations.')
    rmg.saveEdgeSpecies = saveEdgeSpecies
    rmg.keepIrreversible = keepIrreversible
    rmg.trimolecularProductReversible = trimolecularProductReversible
    rmg.wallTime = wallTime

def generatedSpeciesConstraints(**kwargs):

    validConstraints = [
        'allowed',
        'maximumCarbonAtoms',
        'maximumOxygenAtoms',
        'maximumNitrogenAtoms',
        'maximumSiliconAtoms',
        'maximumSulfurAtoms',
        'maximumHeavyAtoms',
        'maximumRadicalElectrons',
        'maximumSingletCarbenes',
        'maximumCarbeneRadicals',
        'allowSingletO2',
    ]

    for key, value in kwargs.items():
        if key not in validConstraints:
            raise InputError('Invalid generated species constraint {0!r}.'.format(key))
        
        rmg.speciesConstraints[key] = value

def thermoCentralDatabase(host,
                        port,
                        username,
                        password,
                        application):
    
    from rmgpy.data.thermo import ThermoCentralDatabaseInterface
    rmg.thermoCentralDatabase = ThermoCentralDatabaseInterface(host,
                                                            port,
                                                            username,
                                                            password,
                                                            application)
                    

################################################################################

def setGlobalRMG(rmg0):
    """
    sets the global variable rmg to rmg0. This is used to allow for unittesting
    of above methods
    """
    global rmg
    rmg = rmg0

def readInputFile(path, rmg0):
    """
    Read an RMG input file at `path` on disk into the :class:`RMG` object 
    `rmg`.
    """
    global rmg, speciesDict
    
    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise e

    logging.info('Reading input file "{0}"...'.format(full_path))
    logging.info(f.read())
    f.seek(0)# return to beginning of file

    setGlobalRMG(rmg0)
    rmg.reactionModel = CoreEdgeReactionModel()
    rmg.initialSpecies = []
    rmg.reactionSystems = []
    speciesDict = {}
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'database': database,
        'species': species,
        'SMARTS': SMARTS,
        'SMILES': SMILES,
        'InChI': InChI,
        'adjacencyList': adjacencyList,
        'simpleReactor': simpleReactor,
        'liquidReactor': liquidReactor,
        'simulator': simulator,
        'solvation': solvation,
        'model': model,
        'quantumMechanics': quantumMechanics,
        'mlEstimator': mlEstimator,
        'pressureDependence': pressureDependence,
        'options': options,
        'generatedSpeciesConstraints': generatedSpeciesConstraints,
        'thermoCentralDatabase': thermoCentralDatabase
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()
    
    rmg.speciesConstraints['explicitlyAllowedMolecules'] = []         
    broadcast(rmg.speciesConstraints, 'speciesConstraints')

    # convert keys from species names into species objects.
    for reactionSystem in rmg.reactionSystems:
        reactionSystem.convertInitialKeysToSpeciesObjects(speciesDict)

    if rmg.quantumMechanics:
        rmg.quantumMechanics.setDefaultOutputDirectory(rmg.outputDirectory)
        rmg.quantumMechanics.initialize()
    broadcast(rmg.quantumMechanics, 'quantumMechanics')

    logging.info('')
    
################################################################################

def readThermoInputFile(path, rmg0):
    """
    Read an thermo estimation input file at `path` on disk into the :class:`RMG` object 
    `rmg`.
    """

    global rmg, speciesDict
    
    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise e

    logging.info('Reading input file "{0}"...'.format(full_path))

    rmg = rmg0
    rmg.reactionModel = CoreEdgeReactionModel()
    rmg.initialSpecies = []
    rmg.reactionSystems = []
    speciesDict = {}
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'database': database,
        'species': species,
        'SMARTS': SMARTS,
        'SMILES': SMILES,
        'InChI': InChI,
        'solvation': solvation,
        'adjacencyList': adjacencyList,
        'quantumMechanics': quantumMechanics,
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    if rmg.quantumMechanics:
        rmg.quantumMechanics.setDefaultOutputDirectory(rmg.outputDirectory)
        rmg.quantumMechanics.initialize()
    broadcast(rmg.quantumMechanics, 'quantumMechanics')
    
    logging.info('')    

################################################################################

def saveInputFile(path, rmg):
    """
    Save an RMG input file at `path` on disk from the :class:`RMG` object 
    `rmg`.
    """

    f = open(path, 'w')

    # Databases
    f.write('database(\n')
    #f.write('    "{0}",\n'.format(rmg.databaseDirectory))
    f.write('    thermoLibraries = {0!r},\n'.format(rmg.thermoLibraries))
    f.write('    reactionLibraries = {0!r},\n'.format(rmg.reactionLibraries))
    f.write('    seedMechanisms = {0!r},\n'.format(rmg.seedMechanisms))
    f.write('    kineticsDepositories = {0!r},\n'.format(rmg.kineticsDepositories))
    f.write('    kineticsFamilies = {0!r},\n'.format(rmg.kineticsFamilies))
    f.write('    kineticsEstimator = {0!r},\n'.format(rmg.kineticsEstimator))
    f.write(')\n\n')

    # Species
    for species in rmg.initialSpecies:
        f.write('species(\n')
        f.write('    label = "{0}",\n'.format(species.label))
        f.write('    reactive = {0},\n'.format(species.reactive))
        f.write('    structure = adjacencyList(\n')
        f.write('"""\n')
        f.write(species.molecule[0].toAdjacencyList())
        f.write('"""),\n')
        f.write(')\n\n')
    
    # Reaction systems
    for system in rmg.reactionSystems:
        if rmg.solvent:
            f.write('liquidReactor(\n')
            f.write('    temperature = ({0:g},"{1!s}"),\n'.format(system.T.getValue(),system.T.units))
            f.write('    initialConcentrations={\n')
            for species, conc in system.initialConcentrations.iteritems():
                f.write('        "{0!s}": ({1:g},"{2!s}"),\n'.format(species.label,conc.getValue(),conc.units))
        else:
            f.write('simpleReactor(\n')
            f.write('    temperature = ({0:g},"{1!s}"),\n'.format(system.T.getValue(),system.T.units))
            # Convert the pressure from SI pascal units to bar here
            # Do something more fancy later for converting to user's desired units for both T and P..
            f.write('    pressure = ({0:g},"{1!s}"),\n'.format(system.P.getValue(),system.P.units))
            f.write('    initialMoleFractions={\n')
            for species, molfrac in system.initialMoleFractions.iteritems():
                f.write('        "{0!s}": {1:g},\n'.format(species.label, molfrac))
        f.write('    },\n')               
        
        # Termination criteria
        conversions = ''
        for term in system.termination:
            if isinstance(term, TerminationTime):
                f.write('    terminationTime = ({0:g},"{1!s}"),\n'.format(term.time.getValue(),term.time.units))
                
            else:
                conversions += '        "{0:s}": {1:g},\n'.format(term.species.label, term.conversion)
        if conversions:        
            f.write('    terminationConversion = {\n')
            f.write(conversions)
            f.write('    },\n')
        
        # Sensitivity analysis
        if system.sensitiveSpecies:
            sensitivity = []
            for item in system.sensitiveSpecies:
                sensitivity.append(item.label)
            f.write('    sensitivity = {0},\n'.format(sensitivity))
            f.write('    sensitivityThreshold = {0},\n'.format(system.sensitivityThreshold))
        
        f.write(')\n\n')
    
    if rmg.solvent:
        f.write("solvation(\n    solvent = '{0!s}'\n)\n\n".format(rmg.solvent))
        
    # Simulator tolerances
    f.write('simulator(\n')
    f.write('    atol = {0:g},\n'.format(rmg.absoluteTolerance))
    f.write('    rtol = {0:g},\n'.format(rmg.relativeTolerance))
    f.write('    sens_atol = {0:g},\n'.format(rmg.sensitivityAbsoluteTolerance))
    f.write('    sens_rtol = {0:g},\n'.format(rmg.sensitivityRelativeTolerance))
    f.write(')\n\n')

    # Model
    f.write('model(\n')
    f.write('    toleranceMoveToCore = {0:g},\n'.format(rmg.fluxToleranceMoveToCore))
    f.write('    toleranceKeepInEdge = {0:g},\n'.format(rmg.fluxToleranceKeepInEdge))
    f.write('    toleranceInterruptSimulation = {0:g},\n'.format(rmg.fluxToleranceInterrupt))
    f.write('    maximumEdgeSpecies = {0:d},\n'.format(rmg.maximumEdgeSpecies))
    f.write('    minCoreSizeForPrune = {0:d},\n'.format(rmg.minCoreSizeForPrune))
    f.write('    minSpeciesExistIterationsForPrune = {0:d},\n'.format(rmg.minSpeciesExistIterationsForPrune))
    f.write('    filterReactions = {0:d},\n'.format(rmg.filterReactions))
    f.write('    filterThreshold = {0:g},\n'.format(rmg.filterThreshold))
    f.write(')\n\n')

    # Pressure Dependence
    if rmg.pressureDependence:
        f.write('pressureDependence(\n')
        f.write('    method = {0!r},\n'.format(rmg.pressureDependence.method))
        f.write('    maximumGrainSize = ({0:g},"{1!s}"),\n'.format(rmg.pressureDependence.grainSize.getValue(),rmg.pressureDependence.grainSize.units))
        f.write('    minimumNumberOfGrains = {0},\n'.format(rmg.pressureDependence.grainCount))
        f.write('    temperatures = ({0:g},{1:g},"{2!s}",{3:d}),\n'.format(
            rmg.pressureDependence.Tmin.getValue(),
            rmg.pressureDependence.Tmax.getValue(),
            rmg.pressureDependence.Tmax.units,
            rmg.pressureDependence.Tcount,
        ))
        f.write('    pressures = ({0:g},{1:g},"{2!s}",{3:d}),\n'.format(
            rmg.pressureDependence.Pmin.getValue(),
            rmg.pressureDependence.Pmax.getValue(),
            rmg.pressureDependence.Pmax.units,
            rmg.pressureDependence.Pcount,
        ))
        f.write('    interpolation = {0},\n'.format(rmg.pressureDependence.interpolationModel))     
        f.write('    maximumAtoms = {0}, \n'.format(rmg.pressureDependence.maximumAtoms))
        f.write(')\n\n')
    
    # Quantum Mechanics
    if rmg.quantumMechanics:
        f.write('quantumMechanics(\n')
        f.write('    software = {0!r},\n'.format(rmg.quantumMechanics.settings.software))
        f.write('    method = {0!r},\n'.format(rmg.quantumMechanics.settings.method))
        # Split paths created by QMSettings
        if rmg.quantumMechanics.settings.fileStore:
            f.write('    fileStore = {0!r},\n'.format(os.path.split(rmg.quantumMechanics.settings.fileStore)[0]))
        else:
            f.write('    fileStore = None,\n')
        if rmg.quantumMechanics.settings.scratchDirectory:
            f.write('    scratchDirectory = {0!r},\n'.format(os.path.split(rmg.quantumMechanics.settings.scratchDirectory)[0]))
        else:
            f.write('    scratchDirectory = None,\n')
        f.write('    onlyCyclics = {0},\n'.format(rmg.quantumMechanics.settings.onlyCyclics))
        f.write('    maxRadicalNumber = {0},\n'.format(rmg.quantumMechanics.settings.maxRadicalNumber))
        f.write(')\n\n')
    
    # Species Constraints
    if rmg.speciesConstraints:
        f.write('generatedSpeciesConstraints(\n')
        for constraint, value in sorted(rmg.speciesConstraints.items(), key=lambda constraint: constraint[0]):
            if value is not None: f.write('    {0} = {1},\n'.format(constraint,value))
        f.write(')\n\n')
    
    # Options
    f.write('options(\n')
    f.write('    units = "{0}",\n'.format(rmg.units))
    if rmg.saveRestartPeriod:
        warnings.warn("The option saveRestartPeriod is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        f.write('    saveRestartPeriod = ({0},"{1}"),\n'.format(rmg.saveRestartPeriod.getValue(), rmg.saveRestartPeriod.units))
    else:
        f.write('    saveRestartPeriod = None,\n')
    f.write('    generateOutputHTML = {0},\n'.format(rmg.generateOutputHTML))
    f.write('    generatePlots = {0},\n'.format(rmg.generatePlots))
    f.write('    saveSimulationProfiles = {0},\n'.format(rmg.saveSimulationProfiles))
    f.write('    saveEdgeSpecies = {0},\n'.format(rmg.saveEdgeSpecies))
    f.write('    keepIrreversible = {0},\n'.format(rmg.keepIrreversible))
    f.write('    trimolecularProductReversible = {0},\n'.format(rmg.trimolecularProductReversible))
    f.write('    verboseComments = {0},\n'.format(rmg.verboseComments))
    f.write('    wallTime = {0},\n'.format(rmg.wallTime))
    f.write(')\n\n')
    
    f.close()

def getInput(name):
    """
    Returns the RMG input object that corresponds
    to the parameter name.

    First, the module level is queried. If this variable
    is empty, the broadcasted variables are queried.
    """
    global rmg

    if rmg:
        if name == 'speciesConstraints':
            return rmg.speciesConstraints
        elif name == 'quantumMechanics':
            return rmg.quantumMechanics
        elif name == 'MLEstimator':
            return rmg.ml_estimator, rmg.ml_settings
        elif name == 'thermoCentralDatabase':
            return rmg.thermoCentralDatabase
        else:
            raise Exception('Unrecognized keyword: {}'.format(name))
    else:
        try:
            obj = get(name)
            if obj:
                return obj
            else:
                raise Exception
        except Exception, e:
            logging.debug("Did not find a way to obtain the variable for {}.".format(name))
            raise e

    raise Exception('Could not get variable with name: {}'.format(name))
