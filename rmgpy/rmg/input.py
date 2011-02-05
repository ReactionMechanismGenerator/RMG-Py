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

from rmgpy.chem.molecule import Molecule

from rmgpy.data import getDatabaseDirectory
from rmgpy.data.thermo import loadThermoDatabase
from rmgpy.data.kinetics import loadKineticsDatabase
from rmgpy.data.states import loadFrequencyDatabase

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

def database(thermo_groups, kinetics_groups, thermo_libraries=None, 
  kinetics_libraries=None, reaction_libraries=None, seed_mechanisms=None,
  frequencies_groups=None, frequencies_libraries=None):
    global databases
    if isinstance(thermo_groups, str): thermo_groups = [thermo_groups]
    if isinstance(kinetics_groups, str): kinetics_groups = [kinetics_groups]
    if isinstance(thermo_libraries, str): thermo_libraries = [thermo_libraries]
    if isinstance(kinetics_libraries, str): kinetics_libraries = [kinetics_libraries]
    if isinstance(reaction_libraries, str): reaction_libraries = [reaction_libraries]
    if isinstance(seed_mechanisms, str): seed_mechanisms = [seed_mechanisms]
    if isinstance(frequencies_groups, str): frequencies_groups = [frequencies_groups]
    if isinstance(frequencies_libraries, str): frequencies_libraries = [frequencies_libraries]
    databases['thermo_groups'] = thermo_groups or []
    databases['kinetics_groups'] = kinetics_groups or []
    databases['thermo_libraries'] = thermo_libraries or []
    databases['kinetics_libraries'] = kinetics_libraries or []
    databases['reaction_libraries'] = reaction_libraries or []
    databases['seed_mechanisms'] = seed_mechanisms or []
    databases['frequencies_groups'] = frequencies_groups or []
    databases['frequencies_libraries'] = frequencies_libraries or []

def species(label, structure, reactive=True):
    global speciesDict, reactionModel
    logging.debug('Found %s species "%s" (%s)' % ('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
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
    pass

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
        logging.error('The input file "%s" could not be opened.' % path)
        logging.info('Check that the file exists and that you have read access.')
        return

    logging.info('Reading input file "%s"...' % path)

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
        logging.error('The input file "%s" was invalid:' % path)
        logging.exception(e)
        raise
    finally:
        f.close()

    logging.info('')

    # Load databases
    for d in databases['thermo_libraries']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadThermoDatabase(path, group=False, old=True)
    for d in databases['thermo_groups']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadThermoDatabase(path, group=True, old=True)
    for d in databases['seed_mechanisms']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadKineticsDatabase(path, group=False, old=True, seedMechanism=True)
    for d in databases['reaction_libraries']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadKineticsDatabase(path, group=False, old=True, reactionLibrary=True)
    for d in databases['kinetics_libraries']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadKineticsDatabase(path, group=False, old=True)
    for d in databases['kinetics_groups']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadKineticsDatabase(path, group=True, old=True)
    for d in databases['frequencies_libraries']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadFrequencyDatabase(path, group=False, old=True)
    for d in databases['frequencies_groups']:
        path = os.path.join(getDatabaseDirectory(), d)
        loadFrequencyDatabase(path, group=True, old=True)

    speciesList = speciesDict.values()
    speciesList.sort(cmp=lambda x, y: x.index - y.index)

    for reactionSystem in reactionSystems:
        initialMoleFractions = {}
        for label, moleFrac in reactionSystem.initialMoleFractions.iteritems():
            initialMoleFractions[speciesDict[label]] = moleFrac
    reactionSystem.initialMoleFractions = initialMoleFractions

    return reactionModel, speciesList, reactionSystems

