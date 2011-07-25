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


def database(path, thermoLibraries=None, reactionLibraries=None, frequenciesLibraries=None, seedMechanisms=None):
    global databases
    # This function just stores the information about the database to be loaded
    # We don't actually load the database until after we're finished reading
    # the input file
    if isinstance(thermoLibraries, str): thermoLibraries = [thermoLibraries]
    if isinstance(reactionLibraries, str): reactionLibraries = [reactionLibraries]
    if isinstance(seedMechanisms, str): seedMechanisms = [seedMechanisms]
    if isinstance(frequenciesLibraries, str): frequenciesLibraries = [frequenciesLibraries]
    databases['path'] = os.path.abspath(os.path.expandvars(path))
    databases['thermoLibraries'] = thermoLibraries or []
    databases['reactionLibraries'] = reactionLibraries or []
    databases['seedMechanisms'] = seedMechanisms or []
    databases['frequenciesLibraries'] = frequenciesLibraries or []

def species(label, structure, reactive=True):
    global speciesDict, reactionModel
    logging.debug('Found {0} species "{1}" ({2})'.format('reactive' if reactive else 'nonreactive', label, structure.toSMILES()))
    speciesDict[label], isNew = reactionModel.makeNewSpecies(structure, label=label, reactive=reactive)

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

#    from rmgpy.measure.input import getTemperaturesForModel, getPressuresForModel
#    # Process temperatures
#    Tmin, Tmin_units, Tmax, Tmax_units, Tcount = temperatures
#    Tmin = processQuantity((Tmin, Tmin_units))[0]
#    Tmax = processQuantity((Tmax, Tmax_units))[0]
#    Tlist = getTemperaturesForModel(interpolation, Tmin, Tmax, Tcount)
#
#    # Process pressures
#    Pmin, Pmin_units, Pmax, Pmax_units, Pcount = pressures
#    Pmin = processQuantity((Pmin, Pmin_units))[0]
#    Pmax = processQuantity((Pmax, Pmax_units))[0]
#    Plist = getPressuresForModel(interpolation, Pmin, Pmax, Pcount)
#
#    # Process grain size
#    minimumGrainSize = processQuantity(minimumGrainSize)[0]
#
#    # Save settings (setting this to non-None enables pressure dependence)
#    settings.pressureDependence = (method, Tmin, Tmax, Tlist, Pmin, Pmax, Plist, minimumGrainSize, minimumNumberOfGrains, interpolation)

    global pdepSettings
    from rmgpy.measure.input import getTemperaturesForModel, getPressuresForModel

    pdepSettings['method']= method
    # Process temperatures
    Tmin, Tmin_units, Tmax, Tmax_units, Tcount = temperatures
    pdepSettings['Tmin'] = processQuantity((Tmin, Tmin_units))[0]
    pdepSettings['Tmax'] = processQuantity((Tmax, Tmax_units))[0]
    pdepSettings['Tlist'] = getTemperaturesForModel(interpolation, Tmin, Tmax, Tcount)
    pdepSettings['Tcount'] = Tcount

    # Process pressures
    Pmin, Pmin_units, Pmax, Pmax_units, Pcount = pressures
    pdepSettings['Pmin'] = processQuantity((Pmin, Pmin_units))[0]
    pdepSettings['Pmax'] = processQuantity((Pmax, Pmax_units))[0]
    pdepSettings['Plist'] = getPressuresForModel(interpolation, Pmin, Pmax, Pcount)
    pdepSettings['Pcount'] = Pcount
    # Process grain size
    pdepSettings['minimumGrainSize'] = processQuantity(minimumGrainSize)[0]

    # Save settings (setting this to non-None enables pressure dependence)
    pdepSettings['minimumNumberOfGrains'] = minimumNumberOfGrains
    pdepSettings['interpolation'] = interpolation
    
    settings.pressureDependence = True

def options(units='si', saveRestart=False, drawMolecules=False, generatePlots=False):
    # settings.saveRestart = saveRestart
    
    # currently drawMolecules, generatePlots don't actually work yet...
    global runSettings
    runSettings['units']=units
    runSettings['saveRestart']=saveRestart
    runSettings['drawMolecules']=drawMolecules
    runSettings['generatePlots']=generatePlots

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
