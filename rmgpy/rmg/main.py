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

"""
This module contains the main execution functionality for Reaction Mechanism
Generator (RMG).
"""

import sys
import warnings
import time
import logging
import os
import shutil

import numpy as np
import gc
import copy
from copy import deepcopy
from scipy.optimize import brute
from cantera import ck2cti

from rmgpy.rmg.settings import ModelSettings
from rmgpy.constraints import failsSpeciesConstraints
from rmgpy.molecule import Molecule
from rmgpy.solver.base import TerminationTime, TerminationConversion
from rmgpy.solver.simple import SimpleReactor
from rmgpy.data.rmg import RMGDatabase
from rmgpy.exceptions import ForbiddenStructureException, DatabaseError, CoreError
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.kinetics.family import KineticsFamily, TemplateReaction
from rmgpy.rmg.pdep import PDepReaction

from rmgpy.data.thermo import ThermoLibrary
from rmgpy.data.base import Entry
from rmgpy import settings

from rmgpy.kinetics.diffusionLimited import diffusionLimiter

from model import Species, CoreEdgeReactionModel
from rmgpy.reaction import Reaction
from pdep import PDepNetwork
import rmgpy.util as util

from rmgpy.chemkin import ChemkinWriter
from rmgpy.rmg.output import OutputHTMLWriter
from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
from rmgpy.restart import RestartWriter
from rmgpy.qm.main import QMDatabaseWriter
from rmgpy.stats import ExecutionStatsWriter
from rmgpy.thermo.thermoengine import submit
from rmgpy.tools.simulate import plot_sensitivity
################################################################################

solvent = None

class RMG(util.Subject):
    """
    A representation of a Reaction Mechanism Generator (RMG) job. The 
    attributes are:
    
    =================================== ================================================
    Attribute                           Description
    =================================== ================================================
    `inputFile`                         The path to the input file
    ----------------------------------- ------------------------------------------------
    `databaseDirectory`                 The directory containing the RMG database
    `thermoLibraries`                   The thermodynamics libraries to load
    `reactionLibraries`                 The kinetics libraries to load
    `statmechLibraries`                 The statistical mechanics libraries to load
    `seedMechanisms`                    The seed mechanisms included in the model
    `kineticsFamilies`                  The kinetics families to use for reaction generation
    `kineticsDepositories`              The kinetics depositories to use for looking up kinetics in each family
    `kineticsEstimator`                 The method to use to estimate kinetics: 'group additivity' or 'rate rules'
    `solvent`                           If solvation estimates are required, the name of the solvent.
    ----------------------------------- ------------------------------------------------
    `reactionModel`                     The core-edge reaction model generated by this job
    `reactionSystems`                   A list of the reaction systems used in this job
    `database`                          The RMG database used in this job
    ----------------------------------- ------------------------------------------------
    `modelSettingsList`                 List of ModelSettings objects containing information related to how to manage species/reaction movement
    `simulatorSettingsList`             List of SimulatorSettings objects containing information on how to run simulations
    `trimolecular`                      ``True`` to consider reactions between three species (i.e., if trimolecular reaction families are present)
    `unimolecularThreshold`             Array of flags indicating whether a species is above the unimolecular reaction threshold
    `bimolecularThreshold`              Array of flags indicating whether two species are above the bimolecular reaction threshold
    `trimolecularThreshold`             Array of flags indicating whether three species are above the trimolecular reaction threshold
    `unimolecularReact`                 Array of flags indicating whether a species should react unimolecularly in the enlarge step
    `bimolecularReact`                  Array of flags indicating whether two species should react in the enlarge step
    `trimolecularReact`                 Array of flags indicating whether three species should react in the enlarge step
    `termination`                       A list of termination targets (i.e :class:`TerminationTime` and :class:`TerminationConversion` objects)
    `speciesConstraints`                Dictates the maximum number of atoms, carbons, electrons, etc. generated by RMG
    ----------------------------------- ------------------------------------------------
    `outputDirectory`                   The directory used to save output files
    `verbosity`                         The level of logging verbosity for console output
    `loadRestart`                       ``True`` if restarting a previous job, ``False`` otherwise
    `saveRestartPeriod`                 The time period to periodically save a restart file (:class:`Quantity`), or ``None`` for never.
    `units`                             The unit system to use to save output files (currently must be 'si')
    `generateOutputHTML`                ``True`` to draw pictures of the species and reactions, saving a visualized model in an output HTML file.  ``False`` otherwise
    `generatePlots`                     ``True`` to generate plots of the job execution statistics after each iteration, ``False`` otherwise
    `verboseComments`                   ``True`` to keep the verbose comments for database estimates, ``False`` otherwise
    `saveEdgeSpecies`                   ``True`` to save chemkin and HTML files of the edge species, ``False`` otherwise
    `keepIrreversible`                  ``True`` to keep ireversibility of library reactions as is ('<=>' or '=>'). ``False`` (default) to force all library reactions to be reversible ('<=>')
    `trimolecularProductReversible`     ``True`` (default) to allow families with trimolecular products to react in the reverse direction, ``False`` otherwise
    `pressureDependence`                Whether to process unimolecular (pressure-dependent) reaction networks
    `quantumMechanics`                  Whether to apply quantum mechanical calculations instead of group additivity to certain molecular types.
    `ml_estimator`                      To use thermo estimation with machine learning
    `ml_settings`                       Settings for ML estimation
    `wallTime`                          The maximum amount of CPU time in the form DD:HH:MM:SS to expend on this job; used to stop gracefully so we can still get profiling information
    `kineticsdatastore`                 ``True`` if storing details of each kinetic database entry in text file, ``False`` otherwise
    ----------------------------------- ------------------------------------------------
    `initializationTime`                The time at which the job was initiated, in seconds since the epoch (i.e. from time.time())
    `done`                              Whether the job has completed (there is nothing new to add)
    =================================== ================================================
    
    """
    
    def __init__(self, inputFile=None, outputDirectory=None):
        super(RMG, self).__init__()
        self.inputFile = inputFile
        self.outputDirectory = outputDirectory
        self.clear()
        self.modelSettingsList = []
        self.simulatorSettingsList = []
        self.Tmin = 0.0
        self.Tmax = 0.0
        self.Pmin = 0.0
        self.Pmax = 0.0
    
    def clear(self):
        """
        Clear all loaded information about the job (except the file paths).
        """
        self.databaseDirectory = None
        self.thermoLibraries = None
        self.transportLibraries = None
        self.reactionLibraries = None
        self.statmechLibraries = None
        self.seedMechanisms = None
        self.kineticsFamilies = None
        self.kineticsDepositories = None
        self.kineticsEstimator = 'group additivity'
        self.solvent = None
        self.diffusionLimiter = None
        
        self.reactionModel = None
        self.reactionSystems = None
        self.database = None
        self.reactionSystem = None
        
        self.modelSettingsList = []
        self.simulatorSettingsList = []
        self.balanceSpecies = None
        
        self.filterReactions=False
        self.trimolecular = False
        self.unimolecularReact = None
        self.bimolecularReact = None
        self.trimolecularReact = None
        self.unimolecularThreshold = None
        self.bimolecularThreshold = None
        self.trimolecularThreshold = None
        self.termination = []
        
        self.done = False
        self.verbosity = logging.INFO
        self.loadRestart = None
        self.saveRestartPeriod = None
        self.units = 'si'
        self.generateOutputHTML = None
        self.generatePlots = None
        self.saveSimulationProfiles = None
        self.verboseComments = None
        self.saveEdgeSpecies = None
        self.keepIrreversible = None
        self.trimolecularProductReversible = None
        self.pressureDependence = None
        self.quantumMechanics = None
        self.ml_estimator = None
        self.ml_settings = None
        self.speciesConstraints = {}
        self.wallTime = '00:00:00:00'
        self.initializationTime = 0
        self.kineticsdatastore = None
        
        self.name = 'Seed'
        self.generateSeedEachIteration = True
        self.saveSeedToDatabase = False

        self.thermoCentralDatabase = None

        self.execTime = []
    
    def loadInput(self, path=None):
        """
        Load an RMG job from the input file located at `inputFile`, or
        from the `inputFile` attribute if not given as a parameter.
        """
        from input import readInputFile
        if path is None: path = self.inputFile
        readInputFile(path, self)
        self.reactionModel.kineticsEstimator = self.kineticsEstimator
        # If the output directory is not yet set, then set it to the same
        # directory as the input file by default
        if not self.outputDirectory:
            self.outputDirectory = os.path.dirname(path)
        if self.pressureDependence:
            self.pressureDependence.outputFile = self.outputDirectory
            self.reactionModel.pressureDependence = self.pressureDependence
        if self.solvent:
            self.reactionModel.solventName = self.solvent

        self.reactionModel.verboseComments = self.verboseComments
        self.reactionModel.saveEdgeSpecies = self.saveEdgeSpecies
        
        if self.quantumMechanics:
            self.reactionModel.quantumMechanics = self.quantumMechanics
            
    def loadThermoInput(self, path=None):
        """
        Load an Thermo Estimation job from a thermo input file located at `inputFile`, or
        from the `inputFile` attribute if not given as a parameter.
        """
        from input import readThermoInputFile
        if path is None: path = self.inputFile
        if not self.outputDirectory:
            self.outputDirectory = os.path.dirname(path)
        readThermoInputFile(path, self)
        
        if self.quantumMechanics:
            self.reactionModel.quantumMechanics = self.quantumMechanics
        
    def checkInput(self):
        """
        Check for a few common mistakes in the input file.
        """
        if self.pressureDependence:
            for index, reactionSystem in enumerate(self.reactionSystems):
                if reactionSystem.T:
                    logging.info(reactionSystem.T)
                    assert (reactionSystem.T.value_si < self.pressureDependence.Tmax.value_si), "Reaction system T is above pressureDependence range."
                    assert (reactionSystem.T.value_si > self.pressureDependence.Tmin.value_si), "Reaction system T is below pressureDependence range."
                else:
                    assert (reactionSystem.Trange[1].value_si < self.pressureDependence.Tmax.value_si), "Reaction system T is above pressureDependence range."
                    assert (reactionSystem.Trange[0].value_si > self.pressureDependence.Tmin.value_si), "Reaction system T is below pressureDependence range."
                if reactionSystem.P:
                    assert (reactionSystem.P.value_si < self.pressureDependence.Pmax.value_si), "Reaction system P is above pressureDependence range."
                    assert (reactionSystem.P.value_si > self.pressureDependence.Pmin.value_si), "Reaction system P is below pressureDependence range."
                else:
                    assert (reactionSystem.Prange[1].value_si < self.pressureDependence.Pmax.value_si), "Reaction system P is above pressureDependence range."
                    assert (reactionSystem.Prange[0].value_si > self.pressureDependence.Pmin.value_si), "Reaction system P is below pressureDependence range."
                

            assert any([not s.reactive for s in reactionSystem.initialMoleFractions.keys()]), \
                "Pressure Dependence calculations require at least one inert (nonreacting) species for the bath gas."

    def checkLibraries(self):
        """
        Check unwanted use of libraries:
        Liquid phase libraries in Gas phase simulation.
        Loading a Liquid phase library obtained in another solvent than the one defined in the input file.
        Other checks can be added here.
        """
        #Liquid phase simulation checks
        if self.solvent:
            #check thermo librairies
            for libIter in self.database.thermo.libraries.iterkeys():
                if self.database.thermo.libraries[libIter].solvent:
                    if not self.solvent ==  self.database.thermo.libraries[libIter].solvent:
                        raise DatabaseError('''Thermo library "{2}" was obtained in "{1}" and cannot be used with this liquid phase simulation in "{0}"
                        '''.format(self.solvent, self.database.thermo.libraries[libIter].solvent, self.database.thermo.libraries[libIter].name))   
            #Check kinetic librairies
            for libIter in self.database.kinetics.libraries.iterkeys():
                if self.database.kinetics.libraries[libIter].solvent:
                    if not self.solvent ==  self.database.kinetics.libraries[libIter].solvent:
                        raise DatabaseError('''Kinetics library "{2}" was obtained in "{1}" and cannot be used with this liquid phase simulation in "{0}"
                        '''.format(self.solvent, self.database.kinetics.libraries[libIter].solvent, self.database.kinetics.libraries[libIter].name))
        #Gas phase simulation checks
        else:
            #check thermo librairies
            for libIter in self.database.thermo.libraries.iterkeys():
                if self.database.thermo.libraries[libIter].solvent:
                    raise DatabaseError('''Thermo library "{1}" was obtained in "{0}" solvent and cannot be used in gas phase simulation
                    '''.format(self.database.thermo.libraries[libIter].solvent, self.database.thermo.libraries[libIter].name))   
            #Check kinetic librairies
            for libIter in self.database.kinetics.libraries.iterkeys():
                if self.database.kinetics.libraries[libIter].solvent:
                    raise DatabaseError('''Kinetics library "{1}" was obtained in "{0}" solvent and cannot be used in gas phase simulation
                    '''.format(self.database.kinetics.libraries[libIter].solvent, self.database.kinetics.libraries[libIter].name))
    
    def saveInput(self, path=None):
        """
        Save an RMG job to the input file located at `path`, or
        from the `outputFile` attribute if not given as a parameter.
        """
        from input import saveInputFile
        if path is None: path = self.outputFile
        saveInputFile(path, self)
        
    def loadDatabase(self):
        
        self.database = RMGDatabase()
        self.database.load(
            path = self.databaseDirectory,
            thermoLibraries = self.thermoLibraries,
            transportLibraries = self.transportLibraries,
            reactionLibraries = [library for library, option in self.reactionLibraries],
            seedMechanisms = self.seedMechanisms,
            kineticsFamilies = self.kineticsFamilies,
            kineticsDepositories = self.kineticsDepositories,
            #frequenciesLibraries = self.statmechLibraries,
            depository = False, # Don't bother loading the depository information, as we don't use it
        )

        # Turn off reversibility for families with three products if desired
        if not self.trimolecularProductReversible:
            for family in self.database.kinetics.families.itervalues():
                if len(family.forwardTemplate.products) > 2:
                    family.reversible = False
                    family.reverseTemplate = None
                    family.reverseRecipe = None
                    family.reverse = None

        # Determine if trimolecular families are present
        for family in self.database.kinetics.families.itervalues():
            if len(family.forwardTemplate.reactants) > 2:
                logging.info('Trimolecular reactions are turned on')
                self.trimolecular = True
                break
        # Only check products if we want to react them
        if not self.trimolecular and self.trimolecularProductReversible:
            for family in self.database.kinetics.families.itervalues():
                if len(family.forwardTemplate.products) > 2:
                    logging.info('Trimolecular reactions are turned on')
                    self.trimolecular = True
                    break

        #check libraries
        self.checkLibraries()
        
        #set global variable solvent
        if self.solvent:
            global solvent
            solvent=self.solvent
        
        if self.kineticsEstimator == 'rate rules':
            if '!training' not in self.kineticsDepositories:
                logging.info('Adding rate rules from training set in kinetics families...')
                # Temporarily remove species constraints for the training reactions
                copySpeciesConstraints=copy.copy(self.speciesConstraints)
                self.speciesConstraints={}
                for family in self.database.kinetics.families.values():
                    family.addKineticsRulesFromTrainingSet(thermoDatabase=self.database.thermo)

                    #If requested by the user, write a text file for each kinetics family detailing the source of each entry
                    if self.kineticsdatastore:
                        logging.info('Writing sources of kinetic entries in family {0} to text file'.format(family.label))
                        path = os.path.join(self.outputDirectory, 'kinetics_database', family.label + '.txt')
                        with open(path, 'w') as f:
                            for template_label, entries in family.rules.entries.iteritems():
                                f.write("Template [{0}] uses the {1} following source(s):\n".format(template_label,str(len(entries))))
                                for entry_index, entry in enumerate(entries):
                                    f.write(str(entry_index+1) + ". " + entry.shortDesc + "\n" + entry.longDesc + "\n")
                                f.write('\n')
                            f.write('\n')

                self.speciesConstraints=copySpeciesConstraints
            else:
                logging.info('Training set explicitly not added to rate rules in kinetics families...')
            logging.info('Filling in rate rules in kinetics families by averaging...')
            for family in self.database.kinetics.families.values():
                family.fillKineticsRulesByAveragingUp(verbose=self.verboseComments)
    
    def initialize(self, **kwargs):
        """
        Initialize an RMG job using the command-line arguments `args` as returned
        by the :mod:`argparse` package.
        """
    
        # Save initialization time
        self.initializationTime = time.time()
    
        # Log start timestamp
        logging.info('RMG execution initiated at ' + time.asctime() + '\n')
    
        # Print out RMG header
        self.logHeader()
        
        try:
            restart = kwargs['restart']
        except KeyError:
            restart = False

        if restart:
            if not os.path.exists(os.path.join(self.outputDirectory,'restart.pkl')):
                logging.error("Could not find restart file (restart.pkl). Please run without --restart option.")
                raise Exception("No restart file")
            
        # Read input file
        self.loadInput(self.inputFile)

        # Check input file 
        self.checkInput()
        
        #Properly set filterReactions to initialize flags properly
        if len(self.modelSettingsList) > 0:
            self.filterReactions = self.modelSettingsList[0].filterReactions
        
        # See if memory profiling package is available
        try:
            import psutil
        except ImportError:
            logging.info('Optional package dependency "psutil" not found; memory profiling information will not be saved.')
    
        
        # Make output subdirectories
        util.makeOutputSubdirectory(self.outputDirectory, 'pdep')
        util.makeOutputSubdirectory(self.outputDirectory, 'solver')
        util.makeOutputSubdirectory(self.outputDirectory, 'kinetics_database')

        # Specifies if details of kinetic database entries should be stored according to user
        try:
            self.kineticsdatastore = kwargs['kineticsdatastore']
        except KeyError:
            self.kineticsdatastore = False

        # Load databases
        self.loadDatabase()

        # Set trimolecular reactant flags of reaction systems
        if self.trimolecular:
            for reactionSystem in self.reactionSystems:
                reactionSystem.trimolecular = True
        
        # Do all liquid-phase startup things:
        if self.solvent:
            solventData = self.database.solvation.getSolventData(self.solvent)
            diffusionLimiter.enable(solventData, self.database.solvation)
            logging.info("Setting solvent data for {0}".format(self.solvent))

            # Set solvent viscosity for reaction filtering
            for reactionSystem in self.reactionSystems:
                reactionSystem.viscosity = solventData.getSolventViscosity(reactionSystem.T.value_si)

        try:
            self.wallTime = kwargs['walltime']
        except KeyError:
            pass

        data = self.wallTime.split(':')
        if not len(data) == 4:
            raise ValueError('Invalid format for wall time {0}; should be DD:HH:MM:SS.'.format(self.wallTime))
        self.wallTime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3]) + 86400 * int(data[-4])

        # Initialize reaction model
        if restart:
            self.initializeRestartRun(os.path.join(self.outputDirectory,'restart.pkl'))
        else:
    
            # Seed mechanisms: add species and reactions from seed mechanism
            # DON'T generate any more reactions for the seed species at this time
            for seedMechanism in self.seedMechanisms:
                self.reactionModel.addSeedMechanismToCore(seedMechanism, react=False)

            # Reaction libraries: add species and reactions from reaction library to the edge so
            # that RMG can find them if their rates are large enough
            for library, option in self.reactionLibraries:
                self.reactionModel.addReactionLibraryToEdge(library)
                
            # Also always add in a few bath gases (since RMG-Java does)
            for label, smiles in [('Ar','[Ar]'), ('He','[He]'), ('Ne','[Ne]'), ('N2','N#N')]:
                molecule = Molecule().fromSMILES(smiles)
                spec, isNew = self.reactionModel.makeNewSpecies(molecule, label=label, reactive=False)
                if isNew:
                    self.initialSpecies.append(spec)
            
            # Perform species constraints and forbidden species checks on input species
            for spec in self.initialSpecies:
                if self.database.forbiddenStructures.isMoleculeForbidden(spec.molecule[0]):
                    if 'allowed' in self.speciesConstraints and 'input species' in self.speciesConstraints['allowed']:
                        logging.warning('Input species {0} is globally forbidden.  It will behave as an inert unless found in a seed mechanism or reaction library.'.format(spec.label))
                    else:
                        raise ForbiddenStructureException("Input species {0} is globally forbidden. You may explicitly allow it, but it will remain inert unless found in a seed mechanism or reaction library.".format(spec.label))
                if failsSpeciesConstraints(spec):
                    if 'allowed' in self.speciesConstraints and 'input species' in self.speciesConstraints['allowed']:
                        self.speciesConstraints['explicitlyAllowedMolecules'].append(spec.molecule[0])
                        pass
                    else:
                        raise ForbiddenStructureException("Species constraints forbids input species {0}. Please reformulate constraints, remove the species, or explicitly allow it.".format(spec.label))

            # For liquidReactor, checks whether the solvent is listed as one of the initial species.
            if self.solvent:
                solventStructure = self.database.solvation.getSolventStructure(self.solvent)
                self.database.solvation.checkSolventinInitialSpecies(self,solventStructure)

            #Check to see if user has input Singlet O2 into their input file or libraries
            #This constraint is special in that we only want to check it once in the input instead of every time a species is made
            if 'allowSingletO2' in self.speciesConstraints and self.speciesConstraints['allowSingletO2']:
                pass
            else:
                #Here we get a list of all species that from the user input
                allInputtedSpecies=[spec for spec in self.initialSpecies]
                #Because no iterations have taken place, the only things in the core are from seed mechanisms
                allInputtedSpecies.extend(self.reactionModel.core.species)
                #Because no iterations have taken place, the only things in the edge are from reaction libraries
                allInputtedSpecies.extend(self.reactionModel.edge.species)
                
                O2Singlet=Molecule().fromSMILES('O=O')
                for spec in allInputtedSpecies:
                    if spec.isIsomorphic(O2Singlet):
                        raise ForbiddenStructureException("""Species constraints forbids input species {0}
                        RMG expects the triplet form of oxygen for correct usage in reaction families. Please change your input to SMILES='[O][O]'
                        If you actually want to use the singlet state, set the allowSingletO2=True inside of the Species Constraints block in your input file.
                        """.format(spec.label))

            for spec in self.initialSpecies:
                submit(spec,self.solvent)
                
            # Add nonreactive species (e.g. bath gases) to core first
            # This is necessary so that the PDep algorithm can identify the bath gas            
            for spec in self.initialSpecies:
                if not spec.reactive:
                    self.reactionModel.enlarge(spec)
            for spec in self.initialSpecies:
                if spec.reactive:
                    self.reactionModel.enlarge(spec)
            
            #chatelak: store constant SPC indices in the reactor attributes if any constant SPC provided in the input file
            #advantages to write it here: this is run only once (as species indexes does not change over the generation)
            if self.solvent is not None:
                for index, reactionSystem in enumerate(self.reactionSystems):
                    if reactionSystem.constSPCNames is not None: #if no constant species provided do nothing
                        reactionSystem.get_constSPCIndices(self.reactionModel.core.species)  ##call the function to identify indices in the solver         
                                  
            self.initializeReactionThresholdAndReactFlags()


        self.reactionModel.initializeIndexSpeciesDict()
            

    def register_listeners(self):
        """
        Attaches listener classes depending on the options 
        found in the RMG input file.
        """

        self.attach(ChemkinWriter(self.outputDirectory))

        if self.generateOutputHTML:
            self.attach(OutputHTMLWriter(self.outputDirectory))

        if self.saveRestartPeriod:
            warnings.warn("The option saveRestartPeriod is no longer supported and may be"
                          " removed in version 2.3.", DeprecationWarning)
            self.attach(RestartWriter()) 

        if self.quantumMechanics:
            self.attach(QMDatabaseWriter()) 

        self.attach(ExecutionStatsWriter(self.outputDirectory))

        if self.saveSimulationProfiles:

            for index, reactionSystem in enumerate(self.reactionSystems):
                reactionSystem.attach(SimulationProfileWriter(
                    self.outputDirectory, index, self.reactionModel.core.species))   
                reactionSystem.attach(SimulationProfilePlotter(
                    self.outputDirectory, index, self.reactionModel.core.species))  
        

    def execute(self, **kwargs):
        """
        Execute an RMG job using the command-line arguments `args` as returned
        by the :mod:`argparse` package.
        """
    
        self.initialize(**kwargs)

        # register listeners
        self.register_listeners()

        self.done = False

        # determine min and max values for T and P (don't determine P values for liquid reactors)
        self.Tmin = min([r_sys.Trange[0] if r_sys.Trange else r_sys.T for r_sys in self.reactionSystems]).value_si
        self.Tmax = max([r_sys.Trange[1] if r_sys.Trange else r_sys.T for r_sys in self.reactionSystems]).value_si
        try:
            self.Pmin = min([x.Prange[0] if x.Prange else x.P for x in self.reactionSystems]).value_si
            self.Pmax = max([x.Prange[1] if x.Prange else x.P for x in self.reactionSystems]).value_si
        except AttributeError:
            # For LiquidReactor, Pmin and Pmax remain with the default value of `None`
            pass
        
        self.rmg_memories = []

        logging.info('Initialization complete. Starting model generation.\n')

        # Initiate first reaction discovery step after adding all core species
        for index, reactionSystem in enumerate(self.reactionSystems):
            # Initialize memory object to track conditions for ranged reactors
            self.rmg_memories.append(RMG_Memory(reactionSystem, self.balanceSpecies))
            self.rmg_memories[index].generate_cond()
            log_conditions(self.rmg_memories, index)

            # Update react flags
            if self.filterReactions:
                # Run the reaction system to update threshold and react flags
                reactionSystem.initializeModel(
                    coreSpecies=self.reactionModel.core.species,
                    coreReactions=self.reactionModel.core.reactions,
                    edgeSpecies=[],
                    edgeReactions=[],
                    pdepNetworks=self.reactionModel.networkList,
                    atol=self.simulatorSettingsList[0].atol,
                    rtol=self.simulatorSettingsList[0].rtol,
                    filterReactions=True,
                    conditions = self.rmg_memories[index].get_cond(),
                )

                self.updateReactionThresholdAndReactFlags(
                    rxnSysUnimolecularThreshold=reactionSystem.unimolecularThreshold,
                    rxnSysBimolecularThreshold=reactionSystem.bimolecularThreshold,
                    rxnSysTrimolecularThreshold=reactionSystem.trimolecularThreshold,
                )

                logging.info('Generating initial reactions for reaction system {0}...'.format(index + 1))
            else:
                # If we're not filtering reactions, then we only need to react
                # the first reaction system since they share the same core
                if index > 0:
                    continue
                logging.info('Generating initial reactions...')

            # React core species to enlarge edge
            self.reactionModel.enlarge(reactEdge=True,
                unimolecularReact=self.unimolecularReact,
                bimolecularReact=self.bimolecularReact,
                trimolecularReact=self.trimolecularReact)

        if not np.isinf(self.modelSettingsList[0].toleranceThermoKeepSpeciesInEdge):
            self.reactionModel.setThermodynamicFilteringParameters(
                self.Tmax,
                toleranceThermoKeepSpeciesInEdge=self.modelSettingsList[0].toleranceThermoKeepSpeciesInEdge,
                minCoreSizeForPrune=self.modelSettingsList[0].minCoreSizeForPrune,
                maximumEdgeSpecies=self.modelSettingsList[0].maximumEdgeSpecies,
                reactionSystems=self.reactionSystems
            )

        if not np.isinf(self.modelSettingsList[0].toleranceThermoKeepSpeciesInEdge):
            self.reactionModel.thermoFilterDown(maximumEdgeSpecies=self.modelSettingsList[0].maximumEdgeSpecies)
        
        logging.info('Completed initial enlarge edge step.\n')
        
        self.saveEverything()
        
        if self.generateSeedEachIteration:
            self.makeSeedMech(firstTime=True)

        maxNumSpcsHit = False #default
        
        for q,modelSettings in enumerate(self.modelSettingsList):
            if len(self.simulatorSettingsList) > 1: 
                simulatorSettings = self.simulatorSettingsList[q]
            else: #if they only provide one input for simulator use that everytime
                simulatorSettings = self.simulatorSettingsList[0]

            self.filterReactions = modelSettings.filterReactions

            logging.info('Beginning model generation stage {0}...\n'.format(q+1))
            
            self.done = False

            # Main RMG loop
            while not self.done:
                
                self.reactionModel.iterationNum += 1
                self.done = True
                
                allTerminated = True
                numCoreSpecies = len(self.reactionModel.core.species)
                
                prunableSpecies = self.reactionModel.edge.species[:]
                prunableNetworks = self.reactionModel.networkList[:]
                
                for index, reactionSystem in enumerate(self.reactionSystems):
                    
                    reactionSystem.prunableSpecies = prunableSpecies   #these lines reset pruning for a new cycle
                    reactionSystem.prunableNetworks = prunableNetworks
                    reactionSystem.reset_max_edge_species_rate_ratios() 
                    
                    for p in xrange(reactionSystem.nSims):
                        reactorDone = True
                        objectsToEnlarge = []
                        self.reactionSystem = reactionSystem
                        # Conduct simulation
                        logging.info('Conducting simulation of reaction system %s...' % (index+1))
                        prune = True
                        
                        self.reactionModel.adjustSurface()
                        
                        if numCoreSpecies < modelSettings.minCoreSizeForPrune:
                            # Turn pruning off if we haven't reached minimum core size.
                            prune = False
                            
                        try: terminated,resurrected,obj,newSurfaceSpecies,newSurfaceReactions,t,x = reactionSystem.simulate(
                            coreSpecies = self.reactionModel.core.species,
                            coreReactions = self.reactionModel.core.reactions,
                            edgeSpecies = self.reactionModel.edge.species,
                            edgeReactions = self.reactionModel.edge.reactions,
                            surfaceSpecies = self.reactionModel.surface.species,
                            surfaceReactions = self.reactionModel.surface.reactions,
                            pdepNetworks = self.reactionModel.networkList,
                            prune = prune,
                            modelSettings=modelSettings,
                            simulatorSettings = simulatorSettings,
                            conditions = self.rmg_memories[index].get_cond()
                        )
                        except:
                            logging.error("Model core reactions:")
                            if len(self.reactionModel.core.reactions) > 5:
                                logging.error("Too many to print in detail")
                            else:
                                from rmgpy.cantherm.output import prettify
                                logging.error(prettify(repr(self.reactionModel.core.reactions)))
                            if self.generateSeedEachIteration:
                                self.makeSeedMech()
                            else:
                                self.makeSeedMech(firstTime=True)
                            raise
                        
                        self.rmg_memories[index].add_t_conv_N(t,x,len(obj))
                        self.rmg_memories[index].generate_cond()
                        log_conditions(self.rmg_memories,index)
                        
                        if self.generateSeedEachIteration:
                            self.makeSeedMech()
                            
                        reactorDone = self.reactionModel.addNewSurfaceObjects(obj,newSurfaceSpecies,newSurfaceReactions,reactionSystem)
                        
                        allTerminated = allTerminated and terminated
                        logging.info('')
                            
                        # If simulation is invalid, note which species should be added to
                        # the core
                        if obj != [] and not (obj is None):
                            objectsToEnlarge = self.processToSpeciesNetworks(obj)
    
                            reactorDone = False
                        # Enlarge objects identified by the simulation for enlarging
                        # These should be Species or Network objects
                        logging.info('')
    
                        objectsToEnlarge = list(set(objectsToEnlarge))
    
                        # Add objects to enlarge to the core first
                        for objectToEnlarge in objectsToEnlarge:
                            self.reactionModel.enlarge(objectToEnlarge)
                            
                        if modelSettings.filterReactions:
                            # Run a raw simulation to get updated reaction system threshold values
                            # Run with the same conditions as with pruning off
                            tempModelSettings = deepcopy(modelSettings)
                            tempModelSettings.fluxToleranceKeepInEdge = 0
                            if not resurrected:
                                try:
                                    reactionSystem.simulate(
                                        coreSpecies = self.reactionModel.core.species,
                                        coreReactions = self.reactionModel.core.reactions,
                                        edgeSpecies = [],
                                        edgeReactions = [],
                                        surfaceSpecies = self.reactionModel.surface.species,
                                        surfaceReactions = self.reactionModel.surface.reactions,
                                        pdepNetworks = self.reactionModel.networkList,
                                        modelSettings = tempModelSettings,
                                        simulatorSettings = simulatorSettings,
                                        conditions = self.rmg_memories[index].get_cond()
                                    )
                                except:
                                    self.updateReactionThresholdAndReactFlags(
                                        rxnSysUnimolecularThreshold = reactionSystem.unimolecularThreshold,
                                        rxnSysBimolecularThreshold = reactionSystem.bimolecularThreshold,
                                        rxnSysTrimolecularThreshold = reactionSystem.trimolecularThreshold,
                                        skipUpdate=True)
                                    logging.warn('Reaction thresholds/flags for Reaction System {0} was not updated due to simulation failure'.format(index+1))
                                else:
                                    self.updateReactionThresholdAndReactFlags(
                                        rxnSysUnimolecularThreshold = reactionSystem.unimolecularThreshold,
                                        rxnSysBimolecularThreshold = reactionSystem.bimolecularThreshold,
                                        rxnSysTrimolecularThreshold = reactionSystem.trimolecularThreshold
                                    )
                            else:
                                self.updateReactionThresholdAndReactFlags(
                                    rxnSysUnimolecularThreshold = reactionSystem.unimolecularThreshold,
                                    rxnSysBimolecularThreshold = reactionSystem.bimolecularThreshold,
                                    rxnSysTrimolecularThreshold = reactionSystem.trimolecularThreshold,
                                    skipUpdate = True
                                )
                                logging.warn('Reaction thresholds/flags for Reaction System {0} was not updated due to resurrection'.format(index+1))

                            logging.info('')
                        else:
                            self.updateReactionThresholdAndReactFlags()

                        if not np.isinf(modelSettings.toleranceThermoKeepSpeciesInEdge):
                            self.reactionModel.setThermodynamicFilteringParameters(self.Tmax, toleranceThermoKeepSpeciesInEdge=modelSettings.toleranceThermoKeepSpeciesInEdge,
                                                              minCoreSizeForPrune=modelSettings.minCoreSizeForPrune, 
                                                              maximumEdgeSpecies=modelSettings.maximumEdgeSpecies,
                                                              reactionSystems=self.reactionSystems)
        
                        oldEdgeSize = len(self.reactionModel.edge.reactions)
                        oldCoreSize = len(self.reactionModel.core.reactions)
                        self.reactionModel.enlarge(reactEdge=True, 
                                unimolecularReact=self.unimolecularReact, 
                                bimolecularReact=self.bimolecularReact,
                                trimolecularReact=self.trimolecularReact)
                            
                        if oldEdgeSize != len(self.reactionModel.edge.reactions) or oldCoreSize != len(self.reactionModel.core.reactions):
                            reactorDone = False
                            
                        if not np.isinf(self.modelSettingsList[0].toleranceThermoKeepSpeciesInEdge):
                            self.reactionModel.thermoFilterDown(maximumEdgeSpecies=modelSettings.maximumEdgeSpecies)
                        
                        maxNumSpcsHit = len(self.reactionModel.core.species) >= modelSettings.maxNumSpecies

                        self.saveEverything()

                        if maxNumSpcsHit:  # breaks the nSims loop
                            # self.done is still True, which will break the while loop
                            break

                        if not reactorDone:
                            self.done = False
                        
                    if maxNumSpcsHit:  # breaks the reactionSystems loop
                        break

                if not self.done: # There is something that needs exploring/enlarging

                    # If we reached our termination conditions, then try to prune
                    # species from the edge
                    if allTerminated and modelSettings.fluxToleranceKeepInEdge>0.0:
                        logging.info('Attempting to prune...')
                        self.reactionModel.prune(self.reactionSystems, modelSettings.fluxToleranceKeepInEdge, modelSettings.fluxToleranceMoveToCore, modelSettings.maximumEdgeSpecies, modelSettings.minSpeciesExistIterationsForPrune)
                        # Perform garbage collection after pruning
                        collected = gc.collect()
                        logging.info('Garbage collector: collected %d objects.' % (collected))
    
                # Consider stopping gracefully if the next iteration might take us
                # past the wall time
                if self.wallTime > 0 and len(self.execTime) > 1:
                    t = self.execTime[-1]
                    dt = self.execTime[-1] - self.execTime[-2]
                    if t + 3 * dt > self.wallTime:
                        logging.info('MODEL GENERATION TERMINATED')
                        logging.info('')
                        logging.info('There is not enough time to complete the next iteration before the wall time is reached.')
                        logging.info('The output model may be incomplete.')
                        logging.info('')
                        coreSpec, coreReac, edgeSpec, edgeReac = self.reactionModel.getModelSize()
                        logging.info('The current model core has %s species and %s reactions' % (coreSpec, coreReac))
                        logging.info('The current model edge has %s species and %s reactions' % (edgeSpec, edgeReac))
                        return
                    
            if maxNumSpcsHit: #resets maxNumSpcsHit and continues the settings for loop
                logging.info('The maximum number of species ({0}) has been hit, Exiting stage {1} ...'.format(modelSettings.maxNumSpecies,q+1))
                maxNumSpcsHit = False
                continue
        
        if not self.generateSeedEachIteration:
            self.makeSeedMech(firstTime=True)
            
        # Run sensitivity analysis post-model generation if sensitivity analysis is on
        for index, reactionSystem in enumerate(self.reactionSystems):
            
            if reactionSystem.sensitiveSpecies and reactionSystem.sensConditions:
                logging.info('Conducting sensitivity analysis of reaction system %s...' % (index+1))

                if reactionSystem.sensitiveSpecies == ['all']:
                    reactionSystem.sensitiveSpecies = self.reactionModel.core.species
                    
                sensWorksheet = []
                for spec in reactionSystem.sensitiveSpecies:
                    csvfilePath = os.path.join(self.outputDirectory, 'solver', 'sensitivity_{0}_SPC_{1}.csv'.format(index+1, spec.index))
                    sensWorksheet.append(csvfilePath)
                
                terminated, resurrected,obj, surfaceSpecies, surfaceReactions,t,x = reactionSystem.simulate(
                    coreSpecies = self.reactionModel.core.species,
                    coreReactions = self.reactionModel.core.reactions,
                    edgeSpecies = self.reactionModel.edge.species,
                    edgeReactions = self.reactionModel.edge.reactions,
                    surfaceSpecies = [],
                    surfaceReactions = [],
                    pdepNetworks = self.reactionModel.networkList,
                    sensitivity = True,
                    sensWorksheet = sensWorksheet,
                    modelSettings = ModelSettings(toleranceMoveToCore=1e8,toleranceInterruptSimulation=1e8),
                    simulatorSettings = self.simulatorSettingsList[-1],
                    conditions = reactionSystem.sensConditions,
                )
                
                plot_sensitivity(self.outputDirectory, index, reactionSystem.sensitiveSpecies)

        # generate Cantera files chem.cti & chem_annotated.cti in a designated `cantera` output folder
        try:
            self.generateCanteraFiles(os.path.join(self.outputDirectory, 'chemkin', 'chem.inp'))
            self.generateCanteraFiles(os.path.join(self.outputDirectory, 'chemkin', 'chem_annotated.inp'))
        except EnvironmentError:
            logging.error('Could not generate Cantera files due to EnvironmentError. Check read\write privileges in output directory.')
        
        self.check_model()
        
        # Write output file
        logging.info('')
        logging.info('MODEL GENERATION COMPLETED')
        logging.info('')
        coreSpec, coreReac, edgeSpec, edgeReac = self.reactionModel.getModelSize()
        logging.info('The final model core has %s species and %s reactions' % (coreSpec, coreReac))
        logging.info('The final model edge has %s species and %s reactions' % (edgeSpec, edgeReac))
        
        self.finish()
    
    def check_model(self):
        """
        Run checks on the RMG model
        """
        logging.info('Performing final model checks...')

        # Check that no two species in core or edge are isomorphic
        for i, spc in enumerate(self.reactionModel.core.species):
            for j in xrange(i):
                spc2 = self.reactionModel.core.species[j]
                if spc.isIsomorphic(spc2):
                    raise CoreError(
                        'Although the model has completed, species {0} is isomorphic to species {1} in the core. '
                        'Please open an issue on GitHub with the following output:'
                        '\n{2}\n{3}'.format(spc.label, spc2.label, spc.toAdjacencyList(), spc2.toAdjacencyList())
                    )

        for i, spc in enumerate(self.reactionModel.edge.species):
            for j in xrange(i):
                spc2 = self.reactionModel.edge.species[j]
                if spc.isIsomorphic(spc2):
                    logging.warning(
                        'Species {0} is isomorphic to species {1} in the edge. This does not affect '
                        'the generated model. If you would like to report this to help make RMG better '
                        'please open a GitHub issue with the following output:'
                        '\n{2}\n{3}'.format(spc.label, spc2.label, spc.toAdjacencyList(), spc2.toAdjacencyList())
                    )

        # Check all core reactions (in both directions) for collision limit violation
        violators = []
        num_rxn_violators = 0
        for rxn in self.reactionModel.core.reactions:
            violator_list = rxn.check_collision_limit_violation(t_min=self.Tmin, t_max=self.Tmax,
                                                            p_min=self.Pmin, p_max=self.Pmax)
            if violator_list:
                violators.extend(violator_list)
                num_rxn_violators += 1
        # Whether or not violators were found, rename 'collision_rate_violators.log' if it exists
        new_file = os.path.join(self.outputDirectory, 'collision_rate_violators.log')
        old_file = os.path.join(self.outputDirectory, 'collision_rate_violators_OLD.log')
        if os.path.isfile(new_file):
            # If there are no violators, yet the violators log exists (probably from a previous run
            # in the same folder), rename it.
            if os.path.isfile(old_file):
                os.remove(old_file)
            os.rename(new_file, old_file)
        if violators:
            logging.info("\n")
            logging.warning("{0} CORE reactions violate the collision rate limit!"
                            "\nSee the 'collision_rate_violators.log' for details.\n\n".format(num_rxn_violators))
            with open(new_file, 'w') as violators_f:
                violators_f.write('*** Collision rate limit violators report ***\n'
                                  '"Violation factor" is the ratio of the rate coefficient to the collision limit'
                                  ' rate at the relevant conditions\n\n')
                for violator in violators:
                    rxn_string = str(violator[0])
                    kinetics = violator[0].kinetics
                    comment=''
                    if isinstance(violator[0], TemplateReaction):
                        comment = violator[0].kinetics.comment
                        violator[0].kinetics.comment = ''  # the comment is printed better when outside of the object
                    if isinstance(violator[0], LibraryReaction):
                        comment = 'Kinetic library: {0}'.format(violator[0].library)
                    if isinstance(violator[0], PDepReaction):
                        comment = 'Network #{0}'.format(violator[0].network)
                    direction = violator[1]
                    ratio = violator[2]
                    condition = violator[3]
                    violators_f.write('{0}\n{1}\n{2}\nDirection: {3}\nViolation factor: {4:.2f}\n'
                                      'Violation condition: {5}\n\n'.format(
                                        rxn_string, kinetics, comment, direction, ratio, condition))
                    if isinstance(violator[0], TemplateReaction):
                        # although this is the end of the run, restore the original comment
                        violator[0].kinetics.comment = comment
        else:
            logging.info("No collision rate violators found.")

    def makeSeedMech(self,firstTime=False):
        """
        causes RMG to make a seed mechanism out of the current chem_annotated.inp and species_dictionary.txt
        this seed mechanism is outputted in a seed folder within the run directory and automatically
        added to as the (or replaces the current) 'Seed' thermo and kinetics libraries in database
        
        if run with firstTime=True it will change self.name to be unique within the thermo/kinetics libraries
        by adding integers to the end of the name to prevent overwritting
        """
        
        logging.info('Making seed mechanism...')
        
        name = self.name
        
        if self.saveSeedToDatabase and firstTime: #make sure don't overwrite current libraries
            thermoNames = self.database.thermo.libraries.keys()
            kineticsNames = self.database.kinetics.libraries.keys()
                
            if name in thermoNames or name in kineticsNames: 
                q = 1
                while name+str(q) in thermoNames or name+str(q) in kineticsNames:
                    q += 1
                self.name = name + str(q)
        
        seedDir = os.path.join(self.outputDirectory,'seed')
        
        if firstTime and not os.path.exists(seedDir): #if seed directory does not exist make it
            os.mkdir(seedDir)
        else:
            shutil.rmtree(seedDir) #otherwise delete the old seed and make a new directory
            os.mkdir(seedDir)
            
        speciesList = self.reactionModel.core.species
        reactionList = self.reactionModel.core.reactions
        edgeSpeciesList = self.reactionModel.edge.species
        edgeReactionList = self.reactionModel.edge.reactions
        
        # Make species labels independent
        oldLabels = self.makeSpeciesLabelsIndependent(speciesList)
        edgeOldLabels = self.makeSpeciesLabelsIndependent(edgeSpeciesList)
        
    
        # load kinetics library entries                    
        kineticsLibrary = KineticsLibrary(name=name,autoGenerated=True)
        kineticsLibrary.entries = {}
        for i in range(len(reactionList)):
            reaction = reactionList[i]        
            entry = Entry(
                    index = i+1,
                    label = reaction.toLabeledStr(),
                    item = reaction,
                    data = reaction.kinetics,
                )
            
            if 'rate rule' in reaction.kinetics.comment:
                entry.longDesc = reaction.kinetics.comment
            elif hasattr(reaction,'library') and reaction.library:
                entry.longDesc = 'Originally from reaction library: ' + reaction.library + "\n" + reaction.kinetics.comment
            else:
                entry.longDesc = reaction.kinetics.comment
            
            kineticsLibrary.entries[i+1] = entry
        
        # load kinetics library entries                    
        edgeKineticsLibrary = KineticsLibrary(name=name+'_edge',autoGenerated=True)
        edgeKineticsLibrary.entries = {}
        for i,reaction in enumerate(edgeReactionList):       
            entry = Entry(
                    index = i+1,
                    label = reaction.toLabeledStr(),
                    item = reaction,
                    data = reaction.kinetics,
                )
            try:
                entry.longDesc = 'Originally from reaction library: ' + reaction.library + "\n" + reaction.kinetics.comment
            except AttributeError:
                entry.longDesc = reaction.kinetics.comment
            edgeKineticsLibrary.entries[i+1] = entry
        
        #save in database
        if self.saveSeedToDatabase:
            databaseDirectory = settings['database.directory']
            try:
                os.makedirs(os.path.join(databaseDirectory, 'kinetics', 'libraries',name))
            except:
                pass
            kineticsLibrary.save(os.path.join(databaseDirectory, 'kinetics', 'libraries', name, 'reactions.py'))
            kineticsLibrary.saveDictionary(os.path.join(databaseDirectory, 'kinetics', 'libraries', name, 'dictionary.txt'))
            
            try:
                os.makedirs(os.path.join(databaseDirectory, 'kinetics', 'libraries',name+'_edge'))
            except:
                pass
            edgeKineticsLibrary.save(os.path.join(databaseDirectory, 'kinetics', 'libraries', name+'_edge', 'reactions.py'))
            edgeKineticsLibrary.saveDictionary(os.path.join(databaseDirectory, 'kinetics', 'libraries', name+'_edge', 'dictionary.txt'))

        #save in output directory
        kineticsLibrary.save(os.path.join(seedDir, name, 'reactions.py'))
        kineticsLibrary.saveDictionary(os.path.join(seedDir, name, 'dictionary.txt'))
        
        edgeKineticsLibrary.save(os.path.join(seedDir, name+'_edge', 'reactions.py'))
        edgeKineticsLibrary.saveDictionary(os.path.join(seedDir, name+'_edge', 'dictionary.txt'))
        
        #change labels back so species aren't renamed
        for i,label in enumerate(oldLabels):
            speciesList[i].label = label
        
        for i,label in enumerate(edgeOldLabels):
            edgeSpeciesList[i].label = label
            
    def makeSpeciesLabelsIndependent(self, species):
        """
        This method looks at the core species labels and makes sure none of them conflict
        If a conflict occurs, the second occurance will have '-2' added
        returns a list of the old labels
        """
        oldLabels = []
        labels = set()
        for spec in species:
            oldLabels.append(spec.label)
            duplicate_index = 1
            if '+' in spec.label:
                L = spec.molecule[0].getFormula()
            else:
                L = spec.label
            potential_label = L
            while potential_label in labels:
                duplicate_index += 1
                potential_label = L + '-{}'.format(duplicate_index)

            spec.label = potential_label
            labels.add(potential_label)
            
            
        return oldLabels
    
    ################################################################################
    def processToSpeciesNetworks(self,obj):
        """
        breaks down the objects returned by simulate into Species and PDepNetwork
        components
        """
        
        if isinstance(obj, PDepNetwork):
            out = [self.processPdepNetworks(obj)]
            return out
        elif isinstance(obj, Species):
            return [obj]
        elif isinstance(obj,Reaction):
            return list(self.processReactionsToSpecies(obj))
        elif isinstance(obj,list): #list of species
            rspcs = self.processReactionsToSpecies([k for k in obj if isinstance(k,Reaction)])
            spcs = {k for k in obj if isinstance(k,Species)} | rspcs
            nworks,pspcs = self.processPdepNetworks([k for k in obj if isinstance(k,PDepNetwork)])
            spcs = list(spcs-pspcs) #avoid duplicate species
            return spcs+nworks
        else:
            raise TypeError("improper call, obj input was incorrect")

    def processPdepNetworks(self,obj):
        """
        properly processes PDepNetwork objects and lists of PDepNetwork objects returned from simulate
        """
        reactionSystem = self.reactionSystem
        if isinstance(obj, PDepNetwork):
            # Determine which species in that network has the highest leak rate
            # We do this here because we need a temperature and pressure
            # Store the maximum leak species along with the associated network
            ob = (obj, obj.getMaximumLeakSpecies(reactionSystem.T.value_si, reactionSystem.P.value_si))
            return ob
        elif isinstance(obj,list):
            spcs = [ob.getMaximumLeakSpecies(reactionSystem.T.value_si, reactionSystem.P.value_si) for ob in obj]
            nworks = [(obj[i],spcs[i]) for i in xrange(len(obj))]
            return nworks,set(spcs)
        else:
            raise TypeError("improper call, obj input was incorrect")
            
    def processReactionsToSpecies(self,obj):
        """
        properly processes Reaction objects and lists of Reaction objects returned from simulate
        """
        coreSpecies = self.reactionModel.core.species
        filterFcn = lambda x: not ((x in coreSpecies)) #remove species already in core
        if isinstance(obj,Reaction):
            potentialSpcs = obj.reactants+obj.products
            potentialSpcs = filter(filterFcn,potentialSpcs)
        elif isinstance(obj,list) or isinstance(obj,set):
            potentialSpcs = set()
            for ob in obj:
                potentialSpcs = potentialSpcs | set(ob.reactants+ob.products)
            potentialSpcs = {sp for sp in potentialSpcs if filterFcn(sp)}
        else:
            raise TypeError("improper call, obj input was incorrect")
        return potentialSpcs

    def generateCanteraFiles(self, chemkinFile, **kwargs):
        """
        Convert a chemkin mechanism chem.inp file to a cantera mechanism file chem.cti
        and save it in the cantera directory
        """
        transportFile = os.path.join(os.path.dirname(chemkinFile), 'tran.dat')
        fileName = os.path.splitext(os.path.basename(chemkinFile))[0] + '.cti'
        outName = os.path.join(self.outputDirectory, 'cantera', fileName)
        canteraDir = os.path.dirname(outName)
        try:
            os.makedirs(canteraDir)
        except OSError:
            if not os.path.isdir(canteraDir):
                raise
        if os.path.exists(outName):
            os.remove(outName)
        parser = ck2cti.Parser()
        parser.convertMech(chemkinFile, transportFile=transportFile, outName=outName, quiet=True, permissive=True, **kwargs)

    def initializeReactionThresholdAndReactFlags(self):
        numCoreSpecies = len(self.reactionModel.core.species)
        if self.filterReactions:
            self.unimolecularReact = np.zeros((numCoreSpecies),bool)
            self.bimolecularReact = np.zeros((numCoreSpecies, numCoreSpecies),bool)
            self.unimolecularThreshold = np.zeros((numCoreSpecies),bool)
            self.bimolecularThreshold = np.zeros((numCoreSpecies, numCoreSpecies),bool)
            if self.trimolecular:
                self.trimolecularReact = np.zeros((numCoreSpecies, numCoreSpecies, numCoreSpecies), bool)
                self.trimolecularThreshold = np.zeros((numCoreSpecies, numCoreSpecies, numCoreSpecies), bool)
        else:
            # By default, react everything
            self.unimolecularReact = np.ones((numCoreSpecies),bool)
            self.bimolecularReact = np.ones((numCoreSpecies, numCoreSpecies),bool)
            if self.trimolecular:
                self.trimolecularReact = np.ones((numCoreSpecies, numCoreSpecies, numCoreSpecies),bool)
            # No need to initialize reaction threshold arrays in this case
    
    def updateReactionThresholdAndReactFlags(self,
                                             rxnSysUnimolecularThreshold=None,
                                             rxnSysBimolecularThreshold=None,
                                             rxnSysTrimolecularThreshold=None,
                                             skipUpdate=False):
        """
        updates the length and boolean value of the unimolecular and bimolecular react and threshold flags
        """
        numCoreSpecies = len(self.reactionModel.core.species)
        prevNumCoreSpecies = len(self.unimolecularReact)
        new_core_species = numCoreSpecies > prevNumCoreSpecies

        # Always reset the react arrays from prior iterations
        self.unimolecularReact = np.zeros((numCoreSpecies), bool)
        self.bimolecularReact = np.zeros((numCoreSpecies, numCoreSpecies), bool)
        if self.trimolecular:
            self.trimolecularReact = np.zeros((numCoreSpecies, numCoreSpecies, numCoreSpecies), bool)

        if self.filterReactions:
            if new_core_species:
                # Expand the threshold arrays if there were new core species added
                unimolecularThreshold = np.zeros((numCoreSpecies), bool)
                bimolecularThreshold = np.zeros((numCoreSpecies, numCoreSpecies), bool)

                # Broadcast original thresholds
                unimolecularThreshold[:prevNumCoreSpecies] = self.unimolecularThreshold
                bimolecularThreshold[:prevNumCoreSpecies,:prevNumCoreSpecies] = self.bimolecularThreshold
                self.unimolecularThreshold = unimolecularThreshold
                self.bimolecularThreshold = bimolecularThreshold

                if self.trimolecular:
                    trimolecularThreshold = np.zeros((numCoreSpecies, numCoreSpecies, numCoreSpecies), bool)
                    trimolecularThreshold[:prevNumCoreSpecies,
                                          :prevNumCoreSpecies,
                                          :prevNumCoreSpecies] = self.trimolecularThreshold
                    self.trimolecularThreshold = trimolecularThreshold
                
            if skipUpdate:
                return
            
            # Always update the react and threshold arrays
            for i in xrange(numCoreSpecies):
                if not self.unimolecularThreshold[i] and rxnSysUnimolecularThreshold[i]:
                    # We've shifted from not reacting to reacting
                    self.unimolecularReact[i] = True
                    self.unimolecularThreshold[i] = True

            for i in xrange(numCoreSpecies):
                for j in xrange(i, numCoreSpecies):
                    if not self.bimolecularThreshold[i,j] and rxnSysBimolecularThreshold[i,j]:
                        # We've shifted from not reacting to reacting
                        self.bimolecularReact[i,j] = True
                        self.bimolecularThreshold[i,j] = True

            if self.trimolecular:
                for i in xrange(numCoreSpecies):
                    for j in xrange(i, numCoreSpecies):
                        for k in xrange(j, numCoreSpecies):
                            if not self.trimolecularThreshold[i,j,k] and rxnSysTrimolecularThreshold[i,j,k]:
                                # We've shifted from not reacting to reacting
                                self.trimolecularReact[i,j,k] = True
                                self.trimolecularThreshold[i,j,k] = True
        else:
            # We are not filtering reactions
            if new_core_species:
                # React all the new core species unimolecularly
                for i in xrange(prevNumCoreSpecies, numCoreSpecies):
                    self.unimolecularReact[i] = True
                
                # React all the new core species with all the core species bimolecularly
                for i in xrange(numCoreSpecies):
                    for j in xrange(prevNumCoreSpecies,numCoreSpecies):
                        self.bimolecularReact[i,j] = True

                # React all the new core species with all bimolecular combinations trimolecularly
                if self.trimolecular:
                    for i in xrange(numCoreSpecies):
                        for j in xrange(numCoreSpecies):
                            for k in xrange(prevNumCoreSpecies, numCoreSpecies):
                                self.trimolecularReact[i,j,k] = True

        
    def saveEverything(self):
        """
        Saves the output HTML, the Chemkin file, and the Restart file (if appropriate).
        
        The restart file is only saved if self.saveRestartPeriod or self.done.
        """
        # If the user specifies it, add unused reaction library reactions to
        # an additional output species and reaction list which is written to the ouput HTML
        # file as well as the chemkin file
        
        if self.reactionLibraries:
            # First initialize the outputReactionList and outputSpeciesList to empty
            self.reactionModel.outputSpeciesList = []
            self.reactionModel.outputReactionList = []
            for library, option in self.reactionLibraries:
                if option:
                    self.reactionModel.addReactionLibraryToOutput(library)
        
        self.execTime.append(time.time() - self.initializationTime)

        # Notify registered listeners:
        self.notify()
            
    def finish(self):
        """
        Complete the model generation.
        """
        # Log end timestamp
        logging.info('')
        logging.info('RMG execution terminated at ' + time.asctime())
    
    def getGitCommit(self, modulePath):
        import subprocess
        if os.path.exists(os.path.join(modulePath,'..','.git')):
            try:
                return subprocess.check_output(['git', 'log',
                                                '--format=%H%n%cd', '-1'],
                                                cwd=modulePath).splitlines()
            except:
                return '', ''
        else:
            return '', ''
    
    def logHeader(self, level=logging.INFO):
        """
        Output a header containing identifying information about RMG to the log.
        """
        from rmgpy import __version__, getPath
        logging.log(level, '#########################################################')
        logging.log(level, '# RMG-Py - Reaction Mechanism Generator in Python       #')
        logging.log(level, '# Version: {0:44s} #'.format(__version__))
        logging.log(level, '# Authors: RMG Developers (rmg_dev@mit.edu)             #')
        logging.log(level, '# P.I.s:   William H. Green (whgreen@mit.edu)           #')
        logging.log(level, '#          Richard H. West (r.west@neu.edu)             #')
        logging.log(level, '# Website: http://reactionmechanismgenerator.github.io/ #')
        logging.log(level, '#########################################################\n')
    
        # Extract git commit from RMG-Py
        head, date = self.getGitCommit(getPath())
        if head != '' and date != '':
            logging.log(level, 'The current git HEAD for RMG-Py is:')
            logging.log(level, '\t%s' % head)
            logging.log(level, '\t%s' % date)
            logging.log(level, '')
        else:
            # If we cannot get git info, try checking if it is a conda package instead:
            condaPackage = get_condaPackage('rmg')
            if condaPackage != '':
                logging.log(level, 'The current anaconda package for RMG-Py is:')
                logging.log(level, condaPackage)
                logging.log(level,'')
                
        databaseHead, databaseDate = self.getGitCommit(settings['database.directory'])
        if databaseHead !='' and databaseDate !='':
            logging.log(level, 'The current git HEAD for RMG-database is:')
            logging.log(level, '\t%s' % databaseHead)
            logging.log(level, '\t%s' % databaseDate)
            logging.log(level, '')
        else:
            databaseCondaPackage=get_condaPackage('rmgdatabase')
            if databaseCondaPackage != '':
                logging.log(level, 'The current anaconda package for RMG-database is:')
                logging.log(level, databaseCondaPackage)
                logging.log(level,'')

    def initializeRestartRun(self, path):

        from rmgpy.rmg.model import getFamilyLibraryObject

        # read restart file
        self.loadRestartFile(path)

        # A few things still point to the species in the input file, so update
        # those to point to the equivalent species loaded from the restart file
    
        # The termination conversions still point to the old species
        from rmgpy.solver.base import TerminationConversion
        for reactionSystem in self.reactionSystems:
            for term in reactionSystem.termination:
                if isinstance(term, TerminationConversion):
                    term.species, isNew = self.reactionModel.makeNewSpecies(term.species.molecule[0], term.species.label, term.species.reactive)
    
        # The initial mole fractions in the reaction systems still point to the old species
        for reactionSystem in self.reactionSystems:
            initialMoleFractions = {}
            for spec0, moleFrac in reactionSystem.initialMoleFractions.iteritems():
                spec, isNew = self.reactionModel.makeNewSpecies(spec0.molecule[0], spec0.label, spec0.reactive)
                initialMoleFractions[spec] = moleFrac
            reactionSystem.initialMoleFractions = initialMoleFractions
    
        # The reactions and reactionDict still point to the old reaction families
        reactionDict = {}
        for family0_label in self.reactionModel.reactionDict:
    
            # Find the equivalent library or family in the newly-loaded kinetics database
            family_label = None
            family0_obj = getFamilyLibraryObject(family0_label)
            if isinstance(family0_obj, KineticsLibrary):
                for label, database in self.database.kinetics.libraries.iteritems():
                    if database.label == family0_label:
                        family_label = database.label
                        break
            elif isinstance(family0_obj, KineticsFamily):
                for label, database in self.database.kinetics.families.iteritems():
                    if database.label == family0_label:
                        family_label = database.label
                        break    
            else:
                import pdb; pdb.set_trace()
            if family_label is None:
                raise Exception("Unable to find matching reaction family for %s" % family0_label)
    
            # Update each affected reaction to point to that new family
            # Also use that new family in a duplicate reactionDict
            reactionDict[family_label] = {}
            for reactant1 in self.reactionModel.reactionDict[family0_label]:
                reactionDict[family_label][reactant1] = {}
                for reactant2 in self.reactionModel.reactionDict[family0_label][reactant1]:
                    reactionDict[family_label][reactant1][reactant2] = []
                    if isinstance(family0_obj, KineticsLibrary):
                        for rxn in self.reactionModel.reactionDict[family0_label][reactant1][reactant2]:
                            assert isinstance(rxn, LibraryReaction)
                            rxn.library = family_label
                            reactionDict[family_label][reactant1][reactant2].append(rxn)
                    elif isinstance(family0_obj, KineticsFamily):
                        for rxn in self.reactionModel.reactionDict[family0_label][reactant1][reactant2]:
                            assert isinstance(rxn, TemplateReaction)
                            rxn.family_label = family_label
                            reactionDict[family_label][reactant1][reactant2].append(rxn)
        
        self.reactionModel.reactionDict = reactionDict
    
    def loadRestartFile(self, path):
        """
        Load a restart file at `path` on disk.
        """
        import cPickle
    
        # Unpickle the reaction model from the specified restart file
        logging.info('Loading previous restart file...')
        f = open(path, 'rb')
        rmg_restart = cPickle.load(f)
        f.close()

        self.reactionModel = rmg_restart.reactionModel
        self.unimolecularReact = rmg_restart.unimolecularReact
        self.bimolecularReact = rmg_restart.bimolecularReact
        self.trimolecularReact = rmg_restart.trimolecularReact
        if self.filterReactions:
            self.unimolecularThreshold = rmg_restart.unimolecularThreshold
            self.bimolecularThreshold = rmg_restart.bimolecularThreshold
            self.trimolecularThreshold = rmg_restart.trimolecularThreshold
        
    def loadRMGJavaInput(self, path):
        """
        Load an RMG-Java job from the input file located at `inputFile`, or
        from the `inputFile` attribute if not given as a parameter.
        """
        warnings.warn("The RMG-Java input is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        # NOTE: This function is currently incomplete!
        # It only loads a subset of the available information.
    
        self.reactionModel = CoreEdgeReactionModel()
        self.initialSpecies = []
        self.reactionSystems = []
    
        Tlist = []; Plist = []; concentrationList = []; speciesDict = {}
        termination = []
        
        with open(path, 'r') as f:
            line = self.readMeaningfulLineJava(f)
            while line != '':
                
                
                if line.startswith('TemperatureModel:'):
                    tokens = line.split()
                    units = tokens[2][1:-1]
                    assert units in ['C', 'F', 'K']
                    if units == 'C':
                        Tlist = [float(T)+273.15 for T in tokens[3:]]
                    elif units == 'F':
                        Tlist = [(float(T)+459.67)*5./9. for T in tokens[3:]]
                    else:
                        Tlist = [float(T) for T in tokens[3:]]
                
                elif line.startswith('PressureModel:'):
                    tokens = line.split()
                    units = tokens[2][1:-1]
                    assert units in ['atm', 'bar', 'Pa', 'torr']
                    if units == 'atm':
                        Plist = [float(P)*101325. for P in tokens[3:]]
                    elif units == 'bar':
                        Plist = [float(P)*100000. for P in tokens[3:]]
                    elif units == 'torr':
                        Plist = [float(P)/760.*101325. for P in tokens[3:]]
                    else:
                        Plist = [float(P) for P in tokens[3:]]
                        
                elif line.startswith('InitialStatus:'):
                    label = ''; concentrations = []; adjlist = ''
                    
                    line = self.readMeaningfulLineJava(f)
                    while line != 'END':
                        
                        if line == '' and label != '':
                            species = Species(label=label, molecule=[Molecule().fromAdjacencyList(adjlist)])
                            self.initialSpecies.append(species)
                            speciesDict[label] = species
                            concentrationList.append(concentrations)
                            label = ''; concentrations = []; adjlist = ''
                        
                        elif line != '' and label == '':
                            tokens = line.split()
                            label = tokens[0]
                            units = tokens[1][1:-1]
                            if tokens[-1] in ['Unreactive', 'ConstantConcentration']:
                                tokens.pop(-1)
                            assert units in ['mol/cm3', 'mol/m3', 'mol/l']
                            if units == 'mol/cm3':
                                concentrations = [float(C)*1.0e6 for C in tokens[2:]]
                            elif units == 'mol/l':
                                concentrations = [float(C)*1.0e3 for C in tokens[2:]]
                            else:
                                concentrations = [float(C) for C in tokens[2:]]
                        
                        elif line != '':
                            adjlist += line + '\n'
                        
                        line = f.readline().strip()
                        if '//' in line: line = line[0:line.index('//')]
                        
                elif line.startswith('InertGas:'):
                    
                    line = self.readMeaningfulLineJava(f)
                    while line != 'END':
                        
                        tokens = line.split()
                        label = tokens[0]
                        assert label in ['N2', 'Ar', 'He', 'Ne']
                        if label == 'Ne':
                            smiles = '[Ne]'
                        elif label == 'Ar':
                            smiles = '[Ar]'
                        elif label == 'He':
                            smiles = '[He]'
                        else:
                            smiles = 'N#N'
                        units = tokens[1][1:-1]
                        assert units in ['mol/cm3', 'mol/m3', 'mol/l']
                        if units == 'mol/cm3':
                            concentrations = [float(C)*1.0e6 for C in tokens[2:]]
                        elif units == 'mol/l':
                            concentrations = [float(C)*1.0e3 for C in tokens[2:]]
                        else:
                            concentrations = [float(C) for C in tokens[2:]]
                        
                        species = Species(label=label, reactive=False, molecule=[Molecule().fromSMILES(smiles)])
                        self.initialSpecies.append(species)
                        speciesDict[label] = species
                        concentrationList.append(concentrations)
                            
                        line = self.readMeaningfulLineJava(f)
                
                elif line.startswith('FinishController:'):
                    
                    # First meaningful line is a termination time or conversion
                    line = self.readMeaningfulLineJava(f)
                    tokens = line.split()
                    if tokens[2].lower() == 'conversion:':
                        label = tokens[3]
                        conversion = float(tokens[4])
                        termination.append(TerminationConversion(spec=speciesDict[label], conv=conversion))
                    elif tokens[2].lower() == 'reactiontime:':
                        time = float(tokens[3])
                        units = tokens[4][1:-1]
                        assert units in ['sec', 'min', 'hr', 'day']
                        if units == 'min':
                            time *= 60.
                        elif units == 'hr':
                            time *= 60. * 60.
                        elif units == 'day':
                            time *= 60. * 60. * 24.
                        termination.append(TerminationTime(time=time))
                            
                    # Second meaningful line is the error tolerance
                    # We're not doing anything with this information yet!
                    line = self.readMeaningfulLineJava(f)
                
                line = self.readMeaningfulLineJava(f)
        
        assert len(Tlist) > 0
        assert len(Plist) > 0
        concentrationList = np.array(concentrationList)
        assert concentrationList.shape[1] > 0  # An arbitrary number of concentrations is acceptable, and should be run for each reactor system 
        
        # Make a reaction system for each (T,P) combination
        for T in Tlist:
            for P in Plist:
                for i in range(concentrationList.shape[1]):
                    concentrations = concentrationList[:,i]
                    totalConc = np.sum(concentrations)
                    initialMoleFractions = dict([(self.initialSpecies[i], concentrations[i] / totalConc) for i in range(len(self.initialSpecies))])
                    reactionSystem = SimpleReactor(T, P, initialMoleFractions=initialMoleFractions, termination=termination)
                    self.reactionSystems.append(reactionSystem)
    
    def readMeaningfulLineJava(self, f):
        """
        Read a meaningful line from an RMG-Java condition file object `f`,
        returning the line with any comments removed.
        """
        warnings.warn("The RMG-Java input is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        line = f.readline()
        if line != '':
            line = line.strip()
            if '//' in line: line = line[0:line.index('//')]
            while line == '':
                line = f.readline()
                if line == '': break
                line = line.strip()
                if '//' in line: line = line[0:line.index('//')]
        return line
    
################################################################################

def initializeLog(verbose, log_file_name):
    """
    Set up a logger for RMG to use to print output to stdout. The
    `verbose` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(verbose)

    # Create console handler and set level to debug; send everything to stdout
    # rather than stderr
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(verbose)

    logging.addLevelName(logging.CRITICAL, 'Critical: ')
    logging.addLevelName(logging.ERROR, 'Error: ')
    logging.addLevelName(logging.WARNING, 'Warning: ')
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, '')
    logging.addLevelName(1, '')

    # Create formatter and add to console handler
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    #formatter = Formatter('%(message)s', '%Y-%m-%d %H:%M:%S')
    formatter = logging.Formatter('%(levelname)s%(message)s')
    ch.setFormatter(formatter)

    # create file handler
    if os.path.exists(log_file_name):
        backup = os.path.join(log_file_name[:-7], 'RMG_backup.log')
        if os.path.exists(backup):
            logging.info("Removing old "+backup)
            os.remove(backup)
        logging.info('Moving {0} to {1}\n'.format(log_file_name, backup))
        shutil.move(log_file_name, backup)
    fh = logging.FileHandler(filename=log_file_name) #, backupCount=3)
    fh.setLevel(min(logging.DEBUG,verbose)) # always at least VERBOSE in the file
    fh.setFormatter(formatter)
    # notice that STDERR does not get saved to the log file
    # so errors from underlying libraries (eg. openbabel) etc. that report
    # on stderr will not be logged to disk.

    # remove old handlers!
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Add console and file handlers to logger
    logger.addHandler(ch)
    logger.addHandler(fh)

################################################################################
class RMG_Memory:
    """
    class for remembering RMG simulations
    and determining what simulation to run next
    """
    def __init__(self,reactionSystem,bspc):
        self.Ranges = dict()
        
        if hasattr(reactionSystem,'Trange') and isinstance(reactionSystem.Trange, list):
            Trange = reactionSystem.Trange
            self.Ranges['T'] = [T.value_si for T in Trange]
        if hasattr(reactionSystem,'Prange') and isinstance(reactionSystem.Prange, list):
            Prange = reactionSystem.Prange
            self.Ranges['P'] = [np.log(P.value_si) for P in Prange]
        if hasattr(reactionSystem,'initialMoleFractions'):
            if bspc:
                self.initialMoleFractions = deepcopy(reactionSystem.initialMoleFractions)
                self.balanceSpecies = [x for x in self.initialMoleFractions.keys() if x.label == bspc][0]  #find the balance species
            for key,value in reactionSystem.initialMoleFractions.iteritems():
                assert key != 'T' and key != 'P', 'naming a species T or P is forbidden'
                if isinstance(value, list):
                    self.Ranges[key] = value
        if hasattr(reactionSystem,'initialConcentrations'):
            for key,value in reactionSystem.initialConcentrations.iteritems():
                assert key != 'T' and key != 'P', 'naming a species T or P is forbidden'
                if isinstance(value, list):
                    self.Ranges[key] = [v.value_si for v in value]
                    
        for term in reactionSystem.termination:
            if isinstance(term, TerminationTime):
                self.tmax = term.time.value_si
        
        self.reactionSystem = reactionSystem
        self.conditionList = []
        self.scaledConditionList = []
        self.ts = []
        self.convs = []
        self.Ns = []
        self.randState = np.random.RandomState(1)
        
    def add_t_conv_N(self,t,conv,N):
        """
        adds the completion time and conversion and the number of objects added 
        from a given run to the memory
        """
        if hasattr(self,'tmax'):
            self.ts.append(t/self.tmax)
        self.convs.append(conv)
        self.Ns.append(N)
        
    def get_cond(self):
        """
        Returns the condition being run
        """
        if self.Ranges == dict():
            return None
        else:
            return self.conditionList[-1]
    
    def calculate_cond(self,obj,Ndims,Ns=20):
        """
        Weighted Stochastic Grid Sampling algorithm 
        obj is evaluated at a grid of points and the evaluations are normalized 
        and then sampled randomly based on their normalized value
        then a random step of length 1/(2*Ns) is taken from that point to give a final condition point
        if this process were to impact runtime under some conditions you could decrease the value of Ns to speed it up
        """
        bounds = tuple((0.0,1.0) for k in xrange(Ndims))
        x0,fval,grid,Jout = brute(obj,bounds,Ns=Ns,full_output=True,finish=None) #run brute just to easily get the evaluations at each grid point (we don't care about the optimal value)
        Jout += abs(Jout.min(tuple(xrange(Ndims)))) #shifts Jout positive so tot is positive
        tot = np.sum(Jout,axis=tuple(xrange(len(Jout.shape))))
        Jout /= tot #normalize Jout
        n = self.randState.uniform(0,1,1)[0] #draw a random number between 0 and 1
        s = 0.0
        for indexes in np.ndenumerate(Jout): #choose a coordinate such that grid[indexes] is choosen with probability Jout[indexes]
            s += Jout[indexes[0]]
            if s > n:
                break
        if len(bounds) != 1:
            yf = np.array([grid[i][indexes[0]] for i in xrange(len(grid))])
        else:
            yf = np.array([grid[indexes[0]] for i in xrange(len(grid))])
        
        step = self.randState.uniform(0,1,len(Jout.shape)) #take a step in a random direction in a length between 0 and 1/(2*Ns)
        step /= step.sum()
        mag = self.randState.uniform(0,1,1)[0]

        yf += step*mag*np.sqrt(2)/(2.0*Ns)
        
        return yf
    
    def generate_cond(self):
        """
        find the next condition to run at by solving an optimization problem
        this optimization problem maximizes distance from prior conditions weighted more if they are more recent
        and maximizes number of objects added
        the resulting condition is added to the end of conditionList
        """
        if self.conditionList == []:
            self.conditionList.append({key:value[0] for key,value in self.Ranges.iteritems()})
            self.scaledConditionList.append({key:0.0 for key,value in self.Ranges.iteritems()})
        elif len(self.conditionList[0]) == 0:
            pass
        else:
            ykey = self.conditionList[0].keys()
            Ns = self.Ns
            def obj(y):
                boo = y.shape == tuple()
                vec = []
                N = len(self.conditionList)
                for i,cond in enumerate(self.scaledConditionList):
                    for j,key in enumerate(ykey):
                        if not boo:
                            vec.append(10.0*N/((N-i)*(Ns[i]+1))*abs(y[j]-cond[key])**0.3) 
                        else:
                            vec.append(10.0*N/((N-i)*(Ns[i]+1))*abs(y-cond[key])**0.3) 
                return -np.array(vec).sum()

            yf = self.calculate_cond(obj,len(ykey))
            
            scaledNewCond = {ykey[i]:yf[i] for i in xrange(len(ykey))}
            newCond = {yk:yf[i]*(self.Ranges[yk][1]-self.Ranges[yk][0])+self.Ranges[yk][0] for i,yk in enumerate(ykey) }
            if 'P' in newCond.keys():
                newCond['P'] = np.exp(newCond['P'])
            
            if hasattr(self,'initialMoleFractions'):
                for key in self.initialMoleFractions.keys():
                    if not isinstance(self.initialMoleFractions[key],list):
                        newCond[key] = self.initialMoleFractions[key]
                total = sum([val for key,val in newCond.iteritems() if key != 'T' and key != 'P'])
                if self.balanceSpecies is None:
                    for key,val in newCond.iteritems():
                        if key != 'T' and key != 'P':
                            newCond[key] = val/total
                else:
                    newCond[self.balanceSpecies] = self.initialMoleFractions[self.balanceSpecies] + 1.0 - total

            self.conditionList.append(newCond)
            self.scaledConditionList.append(scaledNewCond)
        return 

def log_conditions(RMG_Memories,index):
    """
    log newly generated reactor conditions
    """
    if RMG_Memories[index].get_cond() is not None:
        s = 'conditions choosen for reactor {0} were: '.format(index)
        for key,item in RMG_Memories[index].get_cond().iteritems():
            if key == 'T':
                s += 'T = {0} K, '.format(item)
            elif key == 'P':
                s += 'P = {0} bar, '.format(item/1.0e5)
            else:
                s += key.label + ' = {0}, '.format(item)
        
        logging.info(s)
    
class Tee:
    """A simple tee to create a stream which prints to many streams.
    
    This is used to report the profiling statistics to both the log file
    and the standard output.
    """
    def __init__(self, *fileobjects):
        self.fileobjects=fileobjects
    def write(self, string):
        for fileobject in self.fileobjects:
            fileobject.write(string)
            
def get_condaPackage(module):
    """
    Check the version of any conda package
    """
    import subprocess
    try:
        lines = subprocess.check_output(['conda', 'list', '-f', module]).splitlines()
        
        packages=[]
        # Strip comments
        for line in lines:
            if line[:1]=='#':
                pass
            else:
                packages.append(line)
                
        return '\n'.join(packages)
    except:
        return ''

def processProfileStats(stats_file, log_file):
    import pstats
    out_stream = Tee(sys.stdout,open(log_file,'a')) # print to screen AND append to RMG.log
    print >>out_stream, "="*80
    print >>out_stream, "Profiling Data".center(80)
    print >>out_stream, "="*80
    stats = pstats.Stats(stats_file,stream=out_stream)
    stats.strip_dirs()
    print >>out_stream, "Sorted by internal time"
    stats.sort_stats('time')
    stats.print_stats(25)
    stats.print_callers(25)
    print >>out_stream, "Sorted by cumulative time"
    stats.sort_stats('cumulative')
    stats.print_stats(25)
    stats.print_callers(25)
    stats.print_callees(25)

def makeProfileGraph(stats_file):
    """
    Uses gprof2dot to create a graphviz dot file of the profiling information.
    
    This requires the gprof2dot package available via `pip install gprof2dot`.
    Render the result using the program 'dot' via a command like
    `dot -Tps2 input.dot -o output.ps2`.
    
    Rendering the ps2 file to pdf requires an external pdf converter
    `ps2pdf output.ps2` which produces a `output.ps2.pdf` file.
    """
    try:
        from gprof2dot import PstatsParser, DotWriter, SAMPLES, themes
    except ImportError:
        logging.warning('Trouble importing from package gprof2dot. Unable to create a graph of the profile statistics.')
        logging.warning('Try getting the latest version with something like `pip install --upgrade gprof2dot`.')
        return
    import subprocess
    
    #create an Options class to mimic optparser output as much as possible:
    class Options:
        pass
    
    options = Options()
    options.node_thres = 0.8
    options.edge_thres = 0.1
    options.strip = False
    options.show_samples = False
    options.root = ""
    options.leaf = ""
    options.wrap = True
    
    theme = themes['color'] # bw color gray pink
    theme.fontname = "ArialMT" # default "Arial" leads to PostScript warnings in dot (on Mac OS)
    parser = PstatsParser(stats_file)
    profile = parser.parse()
    
    dot_file = stats_file + '.dot'
    output = open(dot_file,'wt')
    dot = DotWriter(output)
    dot.strip = options.strip
    dot.wrap = options.wrap
    
    if options.show_samples:
        dot.show_function_events.append(SAMPLES)
    
    profile = profile
    profile.prune(options.node_thres/100.0, options.edge_thres/100.0)

    if options.root:
        rootId = profile.getFunctionId(options.root)
        if not rootId:
            sys.stderr.write('root node ' + options.root + ' not found (might already be pruned : try -e0 -n0 flags)\n')
            sys.exit(1)
        profile.prune_root(rootId)
    if options.leaf:
        leafId = profile.getFunctionId(options.leaf)
        if not leafId:
            sys.stderr.write('leaf node ' + options.leaf + ' not found (maybe already pruned : try -e0 -n0 flags)\n')
            sys.exit(1)
        profile.prune_leaf(leafId)

    dot.graph(profile, theme)

    output.close()
    
    try:
        subprocess.check_call(['dot', '-Tps2', dot_file, '-o', '{0}.ps2'.format(dot_file)])
    except subprocess.CalledProcessError:
        logging.error("Error returned by 'dot' when generating graph of the profile statistics.")
        logging.info("To try it yourself:\n     dot -Tps2 {0} -o {0}.ps2".format(dot_file))
    except OSError:
        logging.error("Couldn't run 'dot' to create graph of profile statistics. Check graphviz is installed properly and on your path.")
        logging.info("Once you've got it, try:\n     dot -Tps2 {0} -o {0}.ps2".format(dot_file))
    
    try:
        subprocess.check_call(['ps2pdf', '{0}.ps2'.format(dot_file), '{0}.pdf'.format(dot_file)])
    except OSError:
        logging.error("Couldn't run 'ps2pdf' to create pdf graph of profile statistics. Check that ps2pdf converter is installed.")
        logging.info("Once you've got it, try:\n     pd2pdf {0}.ps2 {0}.pdf".format(dot_file))    
    else:
        logging.info("Graph of profile statistics saved to: \n {0}.pdf".format(dot_file))

