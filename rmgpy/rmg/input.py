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

def database(path, thermoLibraries=None, reactionLibraries=None, frequenciesLibraries=None, seedMechanisms=None):
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

def species(label, structure, reactive=True):
    logging.debug('Found {0} species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
    spec, isNew = rmg.reactionModel.makeNewSpecies(structure, label=label, reactive=reactive)
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

    T = Quantity(temperature).value
    P = Quantity(pressure).value
    
    system = SimpleReactor(T, P, initialMoleFractions)
    rmg.reactionSystems.append(system)

def termination(conversion=None, time=None):
    rmg.termination = []
    if conversion is not None:
        for spec, conv in conversion.iteritems():
            rmg.termination.append(TerminationConversion(speciesDict[spec], conv))
    if time is not None:
        rmg.termination.append(TerminationTime(Quantity(time).value))

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
    Tmin, Tmin_units, Tmax, Tmax_units, Tcount = temperatures
    rmg.pressureDependence.Tmin = Quantity(Tmin, Tmin_units)
    rmg.pressureDependence.Tmax = Quantity(Tmax, Tmax_units)
    rmg.pressureDependence.Tcount = Tcount
    Tlist = getTemperaturesForModel(interpolation, rmg.pressureDependence.Tmin.value, rmg.pressureDependence.Tmax.value, rmg.pressureDependence.Tcount)
    rmg.pressureDependence.Tlist = Quantity(Tlist,"K")
    
    # Process pressures
    Pmin, Pmin_units, Pmax, Pmax_units, Pcount = pressures
    rmg.pressureDependence.Pmin = Quantity(Pmin, Pmin_units)
    rmg.pressureDependence.Pmax = Quantity(Pmax, Pmax_units)
    rmg.pressureDependence.Pcount = Pcount
    Plist = getPressuresForModel(interpolation, rmg.pressureDependence.Pmin.value, rmg.pressureDependence.Pmax.value, rmg.pressureDependence.Pcount)
    rmg.pressureDependence.Plist = Quantity(Plist,"Pa")
    
    # Process grain size and count
    rmg.pressureDependence.grainSize = Quantity(maximumGrainSize)
    rmg.pressureDependence.grainCount = minimumNumberOfGrains
   
    # Process interpolation model
    rmg.pressureDependence.model = interpolation

def options(units='si', saveRestart=False, drawMolecules=False, generatePlots=False):
    rmg.units = units
    rmg.saveRestart = saveRestart
    rmg.drawMolecules = drawMolecules
    rmg.generatePlots = generatePlots

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

class InputFile():
    """
    A class for storing the information in a RMG-Py input file as well as loading
    and saving it.

    ======================= =========== ========================================
    Attribute               Type        Description
    ======================= =========== ========================================
    `reactionModel`		        rmgpy.rmg.model.CoreEdgeReactionModel instance
    `reactionSystems`			rmgpy.solver.simple.SimpleReactor object
    `speciesList `          ``list`     list of species objects
    `databases`             ``dict``    dictionary of thermo, reaction, seedmech, etc. libraries
    `pdepSettings`          ``dict``    dictionary of pressure dependence settings
    `runSettings`              ``dict``    run settings
    ======================= =========== ========================================

   """

    def __init__(self,reactionModel=None, reactionSystems = [], speciesList= '', databases={}, pdepSettings={},runSettings={}):
        self.reactionModel = reactionModel
        self.reactionSystems = reactionSystems
        self.speciesList = speciesList
        self.databases = databases
        self.pdepSettings = pdepSettings
        self.runSettings = runSettings

    def load(self, path):
        """
        Load input.py file as InputFile object.
        """
        global reactionModel, reactionSystems, speciesDict, databases, pdepSettings, runSettings
        speciesDict = {}
        databases = {}
        pdepSettings = {}
        runSettings = {}
        reactionSystems = []
        reactionModel = None
        
        full_path = os.path.abspath(os.path.expandvars(path))
        try:
            f = open(full_path)
        except IOError, e:
            logging.error('The input file "{0}" could not be opened.'.format(full_path))
            logging.info('Check that the file exists and that you have read access.')
            raise e

        logging.info('Reading input file "{0}"...'.format(full_path))

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
            logging.error('The input file "{0}" was invalid:'.format(full_path))
            logging.exception(e)
            raise
        finally:
            f.close()


        speciesList = speciesDict.values()
        speciesList.sort(cmp=lambda x, y: x.index - y.index)

        for reactionSystem in reactionSystems:
            initialMoleFractions = {}
            for label, moleFrac in reactionSystem.initialMoleFractions.iteritems():
                initialMoleFractions[speciesDict[label]] = moleFrac
            reactionSystem.initialMoleFractions = initialMoleFractions

        self.reactionModel = reactionModel
        self.reactionSystems = reactionSystems
        self.speciesList = speciesList
        self.databases = databases
        self.pdepSettings = pdepSettings
        self.runSettings = runSettings
        return self

    def save(self, path):
        """
        Saves InputFile object as a input.py file.
        """

        f = open(path, 'w')

        # Databases
        f.write("database(\n")
	f.write("\t'{}',\n".format(self.databases['path']))
	f.write("\tthermoLibraries = {},\n".format(self.databases['thermoLibraries']))
	f.write("\treactionLibraries = {},\n".format(self.databases['reactionLibraries']))
	f.write("\tseedMechanisms = {},\n".format(self.databases['seedMechanisms']))
	f.write(")\n\n")

	# Species
	for species in self.speciesList:
		f.write("species(\n")
		f.write("\tlabel = '{}',\n".format(species.label))
		f.write("\treactive = {},\n".format(species.reactive))
                adjlist = species.molecule[0].toAdjacencyList()
		f.write("\tstructure = adjacencyList(\n")
                f.write('\t\t"""\n')
                for line in adjlist.splitlines():
                    f.write("\t\t{}\n".format(line))
                f.write('\t\t"""),\n')

		f.write(")\n\n")

	# Reaction systems
	for system in self.reactionSystems:
		f.write("simpleReactor(\n")
		f.write("\ttemperature = ({}, 'K'),\n".format(system.T))
		# Convert the pressure from SI pascal units to bar here
		# Do something more fancy later for converting to user's desired units for both T and P..
		f.write("\tpressure = ({}, 'bar'),\n".format(system.P/1e5))

		f.write("\tinitialMoleFractions={\n")
		for species, molfrac in system.initialMoleFractions.iteritems():
			f.write('\t\t"{}": {},\n'.format(species.label, molfrac))
		f.write("\t},\n")

		f.write(")\n\n")

	# Termination
	f.write("termination(\n")
	for term in self.reactionModel.termination:
		if isinstance(term,TerminationTime):
			f.write("\ttime = ({}, 's'),\n".format(term.time))
		if isinstance(term,TerminationConversion):
			f.write("\tconversion = {\n")
			f.write("'{}': {},\n".format(term.species.label, term.conversion))
			f.write("\t}\n")
	f.write(")\n\n")

	# Simulator tolerances
	f.write("simulator(\n")
	f.write("\tatol = {},\n".format(self.reactionModel.absoluteTolerance))
	f.write("\trtol = {},\n".format(self.reactionModel.relativeTolerance))
	f.write(")\n\n")

	# Model
	f.write("model(\n")
	f.write("\ttoleranceKeepInEdge = {},\n".format(self.reactionModel.fluxToleranceKeepInEdge))
    	f.write("\ttoleranceMoveToCore = {},\n".format(self.reactionModel.fluxToleranceMoveToCore))
    	f.write("\ttoleranceInterruptSimulation = {},\n".format(self.reactionModel.fluxToleranceInterrupt))
    	f.write("\tmaximumEdgeSpecies = {},\n".format(self.reactionModel.maximumEdgeSpecies))
	f.write(")\n\n")

	# Pressure Dependence
        if pdepSettings:
            f.write("pressureDependence(\n")
            f.write("\tmethod = '{}',\n".format(self.pdepSettings['method']))
            # change from J/mol to kJ/mol
            f.write("\tminimumGrainSize = ({},'kJ/mol'),\n".format(self.pdepSettings['minimumGrainSize']/1e3))
            f.write("\tminimumNumberOfGrains = {},\n".format(self.pdepSettings['minimumNumberOfGrains']))
            f.write("\ttemperatures = ({},'K',{},'K',{}),\n".format(self.pdepSettings['Tmin'],self.pdepSettings['Tmax'],self.pdepSettings['Tcount']))
            # convert pressure to bar here
            f.write("\tpressures = ({},'bar',{},'bar',{}),\n".format(self.pdepSettings['Pmin']/1e5,self.pdepSettings['Pmax']/1e5,self.pdepSettings['Pcount']))
            f.write("\tinterpolation = {},\n".format(self.pdepSettings['interpolation']))
            f.write(")\n\n")
        
	# Options
        f.write("options(\n")
        for property, value in self.runSettings.iteritems():
            f.write("\t{} = '{}',\n".format(property,value))
        f.write(")\n")
        
        f.close()
