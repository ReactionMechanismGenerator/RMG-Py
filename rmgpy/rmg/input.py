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

from rmgpy.molecule import Molecule

from rmgpy.data.rmg import RMGDatabase

from rmgpy.solver.base import TerminationTime, TerminationConversion
from rmgpy.solver.simple import SimpleReactor

from model import *

################################################################################

class InputError(Exception): pass

################################################################################

rmg = None
speciesDict = {}

def database(path, thermoLibraries=None, reactionLibraries=None, frequenciesLibraries=None, seedMechanisms=None, kineticsDepositories='default'):
    # This function just stores the information about the database to be loaded
    # We don't actually load the database until after we're finished reading
    # the input file
    if isinstance(thermoLibraries, str): thermoLibraries = [thermoLibraries]
    if isinstance(reactionLibraries, str): reactionLibraries = [reactionLibraries]
    if isinstance(seedMechanisms, str): seedMechanisms = [seedMechanisms]
    if isinstance(frequenciesLibraries, str): frequenciesLibraries = [frequenciesLibraries]
    rmg.databaseDirectory = os.path.abspath(os.path.expandvars(path))
    rmg.thermoLibraries = thermoLibraries or []
    rmg.reactionLibraries = reactionLibraries or []
    rmg.seedMechanisms = seedMechanisms or []
    rmg.statmechLibraries = frequenciesLibraries or []
    if kineticsDepositories == 'default':
        rmg.kineticsDepositories = ['training']
    elif kineticsDepositories == 'all':
        rmg.kineticsDepositories = None
    else:
        assert isinstance(kineticsDepositories,list), "kineticsDepositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe']."
        rmg.kineticsDepositories = kineticsDepositories

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
def simpleReactor(temperature, pressure, initialMoleFractions):
    logging.debug('Found SimpleReactor reaction system')

    if sum(initialMoleFractions.values()) != 1:
        logging.warning('Initial mole fractions do not sum to one; renormalizing.')
        for spec in initialMoleFractions:
            initialMoleFractions[spec] /= sum(initialMoleFractions.values())

    T = Quantity(temperature)
    P = Quantity(pressure)
    
    system = SimpleReactor(T, P, initialMoleFractions)
    rmg.reactionSystems.append(system)

def termination(conversion=None, time=None):
    rmg.termination = []
    if conversion is not None:
        for spec, conv in conversion.iteritems():
            rmg.termination.append(TerminationConversion(speciesDict[spec], conv))
    if time is not None:
        rmg.termination.append(TerminationTime(Quantity(time)))

def simulator(atol, rtol):
    rmg.absoluteTolerance = atol
    rmg.relativeTolerance = rtol

def model(toleranceMoveToCore, toleranceKeepInEdge=0.0, toleranceInterruptSimulation=1.0, maximumEdgeSpecies=None):
    rmg.fluxToleranceKeepInEdge = toleranceKeepInEdge
    rmg.fluxToleranceMoveToCore = toleranceMoveToCore
    rmg.fluxToleranceInterrupt = toleranceInterruptSimulation
    rmg.maximumEdgeSpecies = maximumEdgeSpecies

def pressureDependence(method, temperatures, pressures, maximumGrainSize=0.0, minimumNumberOfGrains=0, interpolation=None):

    from rmgpy.measure.input import getTemperaturesForModel, getPressuresForModel
    from rmgpy.measure.main import MEASURE
    
    # Setting the pressureDependence attribute to non-None enables pressure dependence
    rmg.pressureDependence = MEASURE()
    
    # Process method
    rmg.pressureDependence.method = method
    
    # Process temperatures
    Tmin, Tmax, T_units, Tcount = temperatures
    rmg.pressureDependence.Tmin = Quantity(Tmin, T_units)
    rmg.pressureDependence.Tmax = Quantity(Tmax, T_units)
    rmg.pressureDependence.Tcount = Tcount
    Tlist = getTemperaturesForModel(interpolation, rmg.pressureDependence.Tmin.value, rmg.pressureDependence.Tmax.value, rmg.pressureDependence.Tcount)
    rmg.pressureDependence.Tlist = Quantity(Tlist,"K")
    
    # Process pressures
    Pmin, Pmax, P_units, Pcount = pressures
    rmg.pressureDependence.Pmin = Quantity(Pmin, P_units)
    rmg.pressureDependence.Pmax = Quantity(Pmax, P_units)
    rmg.pressureDependence.Pcount = Pcount
    Plist = getPressuresForModel(interpolation, rmg.pressureDependence.Pmin.value, rmg.pressureDependence.Pmax.value, rmg.pressureDependence.Pcount)
    rmg.pressureDependence.Plist = Quantity(Plist,"Pa")
    
    # Process grain size and count
    rmg.pressureDependence.grainSize = Quantity(maximumGrainSize)
    rmg.pressureDependence.grainCount = minimumNumberOfGrains

    # Process interpolation model
    rmg.pressureDependence.model = interpolation

def options(units='si', saveRestartPeriod=None, drawMolecules=False, generatePlots=False, saveConcentrationProfiles=True):
    rmg.units = units
    rmg.saveRestartPeriod = Quantity(saveRestartPeriod) if saveRestartPeriod else None
    rmg.drawMolecules = drawMolecules
    rmg.generatePlots = generatePlots
    rmg.saveConcentrationProfiles = saveConcentrationProfiles

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
        'termination': termination,
        'simulator': simulator,
        'model': model,
        'pressureDependence': pressureDependence,
        'options': options,
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
    f.write('    "{0}",\n'.format(rmg.databaseDirectory))
    f.write('    thermoLibraries = {0!r},\n'.format(rmg.thermoLibraries))
    f.write('    reactionLibraries = {0!r},\n'.format(rmg.reactionLibraries))
    f.write('    seedMechanisms = {0!r},\n'.format(rmg.seedMechanisms))
    f.write('    kineticsDepositories = {0!r},\n'.format(rmg.kineticsDepositories))
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
        f.write(')\n\n')

    # Termination
    f.write('termination(\n')
    for term in rmg.termination:
        if isinstance(term,TerminationTime):
            f.write('    time = ({0:g},"{1!s}"),\n'.format(term.time.getValueInGivenUnits(),term.time.units))
        if isinstance(term,TerminationConversion):
            f.write('    conversion = {\n')
            f.write('        "{0:s}": {1:g},\n'.format(term.species.label, term.conversion))
            f.write('    }\n')
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
    f.write(')\n\n')
        
    f.close()
