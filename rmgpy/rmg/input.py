#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import logging
import quantities
import os

from rmgpy import settings

from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.data.rmg import RMGDatabase
from rmgpy.quantity import Quantity
from rmgpy.solver.base import TerminationTime, TerminationConversion
from rmgpy.solver.simple import SimpleReactor

from model import CoreEdgeReactionModel

################################################################################

class InputError(Exception): pass

################################################################################

rmg = None
speciesDict = {}

def database(
             thermoLibraries = None,
             reactionLibraries = None,
             frequenciesLibraries = None,
             seedMechanisms = None,
             kineticsFamilies = 'default',
             kineticsDepositories = 'default',
             kineticsEstimator = 'group additivity',
             ):
    # This function just stores the information about the database to be loaded
    # We don't actually load the database until after we're finished reading
    # the input file
    if isinstance(thermoLibraries, str): thermoLibraries = [thermoLibraries]
    if isinstance(reactionLibraries, str): reactionLibraries = [reactionLibraries]
    if isinstance(seedMechanisms, str): seedMechanisms = [seedMechanisms]
    if isinstance(frequenciesLibraries, str): frequenciesLibraries = [frequenciesLibraries]
    rmg.databaseDirectory = settings['database.directory']
    rmg.thermoLibraries = thermoLibraries or []
    rmg.reactionLibraries = reactionLibraries or []
    rmg.seedMechanisms = seedMechanisms or []
    rmg.statmechLibraries = frequenciesLibraries or []
    rmg.kineticsEstimator = kineticsEstimator
    if kineticsDepositories == 'default':
        rmg.kineticsDepositories = ['training']
    elif kineticsDepositories == 'all':
        rmg.kineticsDepositories = None
    else:
        assert isinstance(kineticsDepositories,list), "kineticsDepositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe']."
        rmg.kineticsDepositories = kineticsDepositories
    if kineticsFamilies == 'default':
        pass
    elif kineticsFamilies == 'all':
        pass
    else:
        assert isinstance(kineticsFamilies,list), "kineticsFamilies should be either 'default', 'all', or a list of names eg. ['H_Abstraction','!Intra_Disproportionation']."
        rmg.kineticsFamilies = kineticsFamilies

def species(label, structure, reactive=True):
    logging.debug('Found {0} species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
    spec, isNew = rmg.reactionModel.makeNewSpecies(structure, label=label, reactive=reactive)
    assert isNew, "Species {0} is a duplicate of {1}. Species in input file must be unique".format(label,spec.label)
    rmg.initialSpecies.append(spec)
    speciesDict[label] = spec
    
def CML(string):
    return Molecule().fromCML(string)

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
                  terminationConversion=None,
                  terminationTime=None,
                  sensitivity=None,
                  sensitivityThreshold=1e-3
                  ):
    logging.debug('Found SimpleReactor reaction system')

    if sum(initialMoleFractions.values()) != 1:
        logging.warning('Initial mole fractions do not sum to one; renormalizing.')
        for spec in initialMoleFractions:
            initialMoleFractions[spec] /= sum(initialMoleFractions.values())

    T = Quantity(temperature)
    P = Quantity(pressure)
    
    termination = []
    if terminationConversion is not None:
        for spec, conv in terminationConversion.iteritems():
            termination.append(TerminationConversion(speciesDict[spec], conv))
    if terminationTime is not None:
        termination.append(TerminationTime(Quantity(terminationTime)))
    if len(termination) == 0:
        raise InputError('No termination conditions specified for reaction system #{0}.'.format(len(rmg.reactionSystems)+2))
    
    sensitivitySpecies = []
    if sensitivity:
        for spec in sensitivity:
            sensitivitySpecies.append(speciesDict[spec])
    system = SimpleReactor(T, P, initialMoleFractions, termination, sensitivitySpecies, sensitivityThreshold)
    rmg.reactionSystems.append(system)

def simulator(atol, rtol):
    rmg.absoluteTolerance = atol
    rmg.relativeTolerance = rtol

def model(toleranceMoveToCore, toleranceKeepInEdge=0.0, toleranceInterruptSimulation=1.0, maximumEdgeSpecies=None):
    rmg.fluxToleranceKeepInEdge = toleranceKeepInEdge
    rmg.fluxToleranceMoveToCore = toleranceMoveToCore
    rmg.fluxToleranceInterrupt = toleranceInterruptSimulation
    rmg.maximumEdgeSpecies = maximumEdgeSpecies

def quantumMechanics(
                    software,
                    fileStore = None,
                    scratchDirectory = None,
                    onlyCyclics = False,
                    maxRadicalNumber = 0,
                    ):
    from rmgpy.qm.main import QMCalculator
    rmg.quantumMechanics = QMCalculator()
    rmg.quantumMechanics.settings.software = software
    rmg.quantumMechanics.settings.fileStore = fileStore
    rmg.quantumMechanics.settings.scratchDirectory = scratchDirectory
    rmg.quantumMechanics.settings.onlyCyclics = onlyCyclics
    rmg.quantumMechanics.settings.maxRadicalNumber = maxRadicalNumber
                    

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

def options(units='si', saveRestartPeriod=None, drawMolecules=False, generatePlots=False, saveConcentrationProfiles=False, verboseComments=False):
    rmg.units = units
    rmg.saveRestartPeriod = Quantity(saveRestartPeriod) if saveRestartPeriod else None
    rmg.drawMolecules = drawMolecules
    rmg.generatePlots = generatePlots
    rmg.saveConcentrationProfiles = saveConcentrationProfiles
    rmg.verboseComments = verboseComments

def generatedSpeciesConstraints(**kwargs):
    validConstraints = [
        'maximumCarbonAtoms',
        'maximumHydrogenAtoms',
        'maximumOxygenAtoms',
        'maximumNitrogenAtoms',
        'maximumSiliconAtoms',
        'maximumSulfurAtoms',
        'maximumHeavyAtoms',
        'maximumRadicalElectrons',
    ]
    constraints = {}
    for key, value in kwargs.items():
        if key not in validConstraints:
            raise InputError('Invalid generated species constraint {0!r}.'.format(key))
        rmg.reactionGenerationOptions[key] = value

################################################################################

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
        'CML': CML,
        'SMILES': SMILES,
        'InChI': InChI,
        'adjacencyList': adjacencyList,
        'simpleReactor': simpleReactor,
        'simulator': simulator,
        'model': model,
        'quantumMechanics': quantumMechanics,
        'pressureDependence': pressureDependence,
        'options': options,
        'generatedSpeciesConstraints': generatedSpeciesConstraints,
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    for reactionSystem in rmg.reactionSystems:
        initialMoleFractions = {}
        for label, moleFrac in reactionSystem.initialMoleFractions.iteritems():
            initialMoleFractions[speciesDict[label]] = moleFrac
        reactionSystem.initialMoleFractions = initialMoleFractions

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
        f.write('simpleReactor(\n')
        f.write('    temperature = ({0:g},"{1!s}"),\n'.format(system.T.getValueInGivenUnits(),system.T.units))
        # Convert the pressure from SI pascal units to bar here
        # Do something more fancy later for converting to user's desired units for both T and P..
        f.write('    pressure = ({0:g},"{1!s}"),\n'.format(system.P.getValueInGivenUnits(),system.P.units))
        f.write('    initialMoleFractions={\n')
        for species, molfrac in system.initialMoleFractions.iteritems():
            f.write('        "{0!s}": {1:g},\n'.format(species.label, molfrac))
        f.write('    },\n')               
        
        # Termination criteria
        conversions = ''
        for term in system.termination:
            if isinstance(term, TerminationTime):
                f.write('    terminationTime = ({0:g},"{1!s}"),\n'.format(term.time.getValueInGivenUnits(),term.time.units))
                
            else:
                conversions += '        "{0:s}": {1:g},\n'.format(term.species.label, term.conversion)
        if conversions:        
            f.write('    terminationConversion = {\n')
            f.write(conversions)
            f.write('    },\n')
        
        # Sensitivity analysis
        if system.sensitivity:
            f.write('    sensitivity = {0},\n'.format(system.sensitivity))       
        if system.sensitivityThreshold:
            f.write('    sensitivityThreshold = {0},\n'.format(system.sensitivity))      
        
        f.write(')\n\n')
        
        
    # Simulator tolerances
    f.write('simulator(\n')
    f.write('    atol = {0:g},\n'.format(rmg.absoluteTolerance))
    f.write('    rtol = {0:g},\n'.format(rmg.relativeTolerance))
    f.write(')\n\n')

    # Model
    f.write('model(\n')
    f.write('    toleranceKeepInEdge = {0:g},\n'.format(rmg.fluxToleranceKeepInEdge))
    f.write('    toleranceMoveToCore = {0:g},\n'.format(rmg.fluxToleranceMoveToCore))
    f.write('    toleranceInterruptSimulation = {0:g},\n'.format(rmg.fluxToleranceInterrupt))
    f.write('    maximumEdgeSpecies = {0:d},\n'.format(rmg.maximumEdgeSpecies))
    f.write(')\n\n')

    # Pressure Dependence
    if rmg.pressureDependence:
        f.write('pressureDependence(\n')
        f.write('    method = "{0!s}",\n'.format(rmg.pressureDependence.method))
        f.write('    maximumGrainSize = ({0:g},"{1!s}"),\n'.format(rmg.pressureDependence.grainSize.getValueInGivenUnits(),rmg.pressureDependence.grainSize.units))
        f.write('    minimumNumberOfGrains = {0},\n'.format(rmg.pressureDependence.grainCount))
        f.write('    temperatures = ({0:g},{1:g},"{2!s}",{3:d}),\n'.format(
            rmg.pressureDependence.Tmin.getValueInGivenUnits(),
            rmg.pressureDependence.Tmax.getValueInGivenUnits(),
            rmg.pressureDependence.Tmax.units,
            rmg.pressureDependence.Tcount,
        ))
        f.write('    pressures = ({0:g},{1:g},"{2!s}",{3:d}),\n'.format(
            rmg.pressureDependence.Pmin.getValueInGivenUnits(),
            rmg.pressureDependence.Pmax.getValueInGivenUnits(),
            rmg.pressureDependence.Pmax.units,
            rmg.pressureDependence.Pcount,
        ))
        f.write('    interpolation = {0},\n'.format(rmg.pressureDependence.model))
        f.write(')\n\n')
        
    # Options
    f.write('options(\n')
    f.write('    units = "{0}",\n'.format(rmg.units))
    if rmg.saveRestartPeriod:
        f.write('    saveRestartPeriod = ({0},"{1}"),\n'.format(rmg.saveRestartPeriod.getValueInGivenUnits(), rmg.saveRestartPeriod.units))
    else:
        f.write('    saveRestartPeriod = None,\n')
    f.write('    drawMolecules = {0},\n'.format(rmg.drawMolecules))
    f.write('    generatePlots = {0},\n'.format(rmg.generatePlots))
    f.write('    saveConcentrationProfiles = {0},\n'.format(rmg.saveConcentrationProfiles))
    f.write('    verboseComments = {0},\n'.format(rmg.verboseComments))
    f.write(')\n\n')
        
    f.close()
