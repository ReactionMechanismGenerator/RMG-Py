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

speciesDict = {}
databases = {}
reactionSystems = []
reactionModel = None

def database(path, thermoLibraries=None, reactionLibraries=None, frequenciesLibraries=None, seedMechanisms=None):
    global databases
    # This function just stores the information about the database to be loaded
    # We don't actually load the database until after we're finished reading
    # the input file
    if isinstance(thermoLibraries, str): thermoLibraries = [thermoLibraries]
    if isinstance(reactionLibraries, str): reactionLibraries = [reactionLibraries]
    if isinstance(seedMechanisms, str): seedMechanisms = [seedMechanisms]
    if isinstance(frequenciesLibraries, str): frequenciesLibraries = [frequenciesLibraries]
    databases['path'] = os.path.abspath(path)
    databases['thermoLibraries'] = thermoLibraries or []
    databases['reactionLibraries'] = reactionLibraries or []
    databases['seedMechanisms'] = seedMechanisms or []
    databases['frequenciesLibraries'] = frequenciesLibraries or []

def species(label, structure, reactive=True):
    global speciesDict, reactionModel
    logging.debug('Found {0} species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
    speciesDict[label], isNew = reactionModel.makeNewSpecies(label=label, molecule=structure, reactive=reactive)

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
    global reactionSystems
    logging.debug('Found SimpleReactor reaction system')

    if sum(initialMoleFractions.values()) != 1:
        logging.warning('Initial mole fractions do not sum to one; renormalizing.')
        for spec in initialMoleFractions:
            initialMoleFractions[spec] /= sum(initialMoleFractions.values())

    T = processQuantity(temperature)[0]
    P = processQuantity(pressure)[0]
    
    system = SimpleReactor(T, P, initialMoleFractions)
    reactionSystems.append(system)

def termination(conversion=None, time=None):
    global reactionModel, speciesDict
    reactionModel.termination = []
    if conversion is not None:
        for spec, conv in conversion.iteritems():
            reactionModel.termination.append(TerminationConversion(speciesDict[spec], conv))
    if time is not None:
        reactionModel.termination.append(TerminationTime(processQuantity(time)[0]))

def simulator(atol, rtol):
    global reactionModel
    reactionModel.absoluteTolerance = atol
    reactionModel.relativeTolerance = rtol

def model(toleranceMoveToCore, toleranceKeepInEdge=0.0, toleranceInterruptSimulation=1.0, maximumEdgeSpecies=None):
    global reactionModel
    reactionModel.fluxToleranceKeepInEdge = toleranceKeepInEdge
    reactionModel.fluxToleranceMoveToCore = toleranceMoveToCore
    reactionModel.fluxToleranceInterrupt = toleranceInterruptSimulation
    reactionModel.maximumEdgeSpecies = maximumEdgeSpecies

def pressureDependence(method, temperatures, pressures, minimumGrainSize=0.0, minimumNumberOfGrains=0, interpolation=None):

    from rmgpy.measure.input import getTemperaturesForModel, getPressuresForModel
    # Process temperatures
    Tmin, Tmin_units, Tmax, Tmax_units, Tcount = temperatures
    Tmin = processQuantity((Tmin, Tmin_units))[0]
    Tmax = processQuantity((Tmax, Tmax_units))[0]
    Tlist = getTemperaturesForModel(interpolation, Tmin, Tmax, Tcount)

    # Process pressures
    Pmin, Pmin_units, Pmax, Pmax_units, Pcount = pressures
    Pmin = processQuantity((Pmin, Pmin_units))[0]
    Pmax = processQuantity((Pmax, Pmax_units))[0]
    Plist = getPressuresForModel(interpolation, Pmin, Pmax, Pcount)
    
    # Process grain size
    minimumGrainSize = processQuantity(minimumGrainSize)[0]

    # Save settings (setting this to non-None enables pressure dependence)
    settings.pressureDependence = (method, Tmin, Tmax, Tlist, Pmin, Pmax, Plist, minimumGrainSize, minimumNumberOfGrains, interpolation)
    
def options(units='si', saveRestart=False, drawMolecules=False, generatePlots=False):
    settings.saveRestart = saveRestart

################################################################################

def processQuantity(quantity):
    if isinstance(quantity, tuple) or isinstance(quantity, list):
        value, units = quantity
    else:
        value = quantity; units = ''
    newUnits = str(quantities.Quantity(1.0, units).simplified.units).split()[1]
    if isinstance(value, tuple) or isinstance(value, list):
        return [float(quantities.Quantity(v, units).simplified) for v in value], newUnits
    else:
        return float(quantities.Quantity(value, units).simplified), newUnits

################################################################################

def readInputFile(path):

    global speciesDict, reactionModel, databases, reactionSystems

    try:
        f = open(path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(path))
        logging.info('Check that the file exists and that you have read access.')
        return

    logging.info('Reading input file "{0}"...'.format(path))

    reactionModel = CoreEdgeReactionModel()

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
        logging.error('The input file "{0}" was invalid:'.format(path))
        logging.exception(e)
        raise
    finally:
        f.close()

    logging.info('')

    # Load databases
    rmgDatabase = RMGDatabase()
    rmgDatabase.load(
        path = databases['path'],
        thermoLibraries = databases['thermoLibraries'],
        reactionLibraries = databases['reactionLibraries'],
        seedMechanisms = databases['seedMechanisms'],
        #frequenciesLibraries = databases['frequenciesLibraries'],
        depository = False, # Don't bother loading the depository information, as we don't use it
    )
    seedMechanisms = databases['seedMechanisms']
    logging.info("")
    
    speciesList = speciesDict.values()
    speciesList.sort(cmp=lambda x, y: x.index - y.index)

    for reactionSystem in reactionSystems:
        initialMoleFractions = {}
        for label, moleFrac in reactionSystem.initialMoleFractions.iteritems():
            initialMoleFractions[speciesDict[label]] = moleFrac
    reactionSystem.initialMoleFractions = initialMoleFractions

    return reactionModel, speciesList, reactionSystems, rmgDatabase, seedMechanisms

