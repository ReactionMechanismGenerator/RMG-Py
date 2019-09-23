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
This module contains the main execution functionality for Reaction Mechanism
Generator (RMG).
"""
from __future__ import print_function

import copy
import gc
import logging
import os
import resource
import shutil
import sys
import time
import warnings
from copy import deepcopy

import h5py
import numpy as np
import psutil
import yaml
from cantera import ck2cti
from scipy.optimize import brute

import rmgpy.util as util
from rmgpy.rmg.model import Species, CoreEdgeReactionModel
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy import settings
from rmgpy.chemkin import ChemkinWriter
from rmgpy.constraints import failsSpeciesConstraints
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.exceptions import ForbiddenStructureException, DatabaseError, CoreError
from rmgpy.kinetics.diffusionLimited import diffusionLimiter
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMDatabaseWriter
from rmgpy.reaction import Reaction
from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
from rmgpy.rmg.output import OutputHTMLWriter
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.rmg.settings import ModelSettings
from rmgpy.solver.base import TerminationTime, TerminationConversion
from rmgpy.solver.simple import SimpleReactor
from rmgpy.stats import ExecutionStatsWriter
from rmgpy.thermo.thermoengine import submit
from rmgpy.tools.plot import plot_sensitivity
from rmgpy.tools.uncertainty import Uncertainty, process_local_results
from rmgpy.yml import RMSWriter

################################################################################

# This module uses the HDF5 data format, which can cause problems on files systems that use NFS (common for network
# mounted file systems. The following sets an environment variable that prevents file locking that would otherwise
# cause a problem for NFS.
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

solvent = None

# Maximum number of user defined processors
maxproc = 1


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
        self.surfaceSiteDensity = None
        self.bindingEnergies = None

        self.reactionModel = None
        self.reactionSystems = None
        self.database = None
        self.reactionSystem = None

        self.modelSettingsList = []
        self.simulatorSettingsList = []
        self.balanceSpecies = None

        self.filterReactions = False
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
        self.restart = False
        self.coreSeedPath = None
        self.edgeSeedPath = None
        self.filtersPath = None
        self.speciesMapPath = None

        self.name = 'Seed'
        self.generateSeedEachIteration = True
        self.saveSeedToDatabase = False

        self.thermoCentralDatabase = None
        self.uncertainty = None

        self.execTime = []

    def loadInput(self, path=None):
        """
        Load an RMG job from the input file located at `inputFile`, or
        from the `inputFile` attribute if not given as a parameter.
        """
        from rmgpy.rmg.input import readInputFile
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

        if self.surfaceSiteDensity:
            self.reactionModel.surfaceSiteDensity = self.surfaceSiteDensity

        self.reactionModel.verboseComments = self.verboseComments
        self.reactionModel.saveEdgeSpecies = self.saveEdgeSpecies

        if self.quantumMechanics:
            self.reactionModel.quantumMechanics = self.quantumMechanics

    def loadThermoInput(self, path=None):
        """
        Load an Thermo Estimation job from a thermo input file located at `inputFile`, or
        from the `inputFile` attribute if not given as a parameter.
        """
        from rmgpy.rmg.input import readThermoInputFile
        if path is None:
            path = self.inputFile
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
        # Liquid phase simulation checks
        if self.solvent:
            # check thermo librairies
            for libIter in self.database.thermo.libraries.keys():
                if self.database.thermo.libraries[libIter].solvent:
                    if not self.solvent == self.database.thermo.libraries[libIter].solvent:
                        raise DatabaseError("Thermo library '{2}' was obtained in '{1}' and cannot be used with this "
                                            "liquid phase simulation in '{0}' "
                                            .format(self.solvent,
                                                    self.database.thermo.libraries[libIter].solvent,
                                                    self.database.thermo.libraries[libIter].name))
            # Check kinetic librairies
            for libIter in self.database.kinetics.libraries.keys():
                if self.database.kinetics.libraries[libIter].solvent:
                    if not self.solvent == self.database.kinetics.libraries[libIter].solvent:
                        raise DatabaseError("Kinetics library '{2}' was obtained in '{1}' and cannot be used with this "
                                            "liquid phase simulation in '{0}'"
                                            .format(self.solvent,
                                                    self.database.kinetics.libraries[libIter].solvent,
                                                    self.database.kinetics.libraries[libIter].name))
        # Gas phase simulation checks
        else:
            # check thermo librairies
            for libIter in self.database.thermo.libraries.keys():
                if self.database.thermo.libraries[libIter].solvent:
                    raise DatabaseError("Thermo library '{1}' was obtained in '{0}' solvent and cannot be used in gas "
                                        "phase simulation"
                                        .format(self.database.thermo.libraries[libIter].solvent,
                                                self.database.thermo.libraries[libIter].name))
                    # Check kinetic librairies
            for libIter in self.database.kinetics.libraries.keys():
                if self.database.kinetics.libraries[libIter].solvent:
                    raise DatabaseError("Kinetics library '{1}' was obtained in '{0}' solvent and cannot be used in "
                                        "gas phase simulation"
                                        .format(self.database.kinetics.libraries[libIter].solvent,
                                                self.database.kinetics.libraries[libIter].name))

    def saveInput(self, path=None):
        """
        Save an RMG job to the input file located at `path`, or
        from the `outputFile` attribute if not given as a parameter.
        """
        from rmgpy.rmg.input import saveInputFile
        if path is None:
            path = self.outputFile
        saveInputFile(path, self)

    def loadDatabase(self):

        self.database = RMGDatabase()
        self.database.load(
            path=self.databaseDirectory,
            thermoLibraries=self.thermoLibraries,
            transportLibraries=self.transportLibraries,
            reactionLibraries=[library for library, option in self.reactionLibraries],
            seedMechanisms=self.seedMechanisms,
            kineticsFamilies=self.kineticsFamilies,
            kineticsDepositories=self.kineticsDepositories,
            # frequenciesLibraries = self.statmechLibraries,
            depository=False,  # Don't bother loading the depository information, as we don't use it
        )

        # Turn off reversibility for families with three products if desired
        if not self.trimolecularProductReversible:
            for family in self.database.kinetics.families.values():
                if len(family.forwardTemplate.products) > 2:
                    family.reversible = False
                    family.reverseTemplate = None
                    family.reverseRecipe = None
                    family.reverse = None

        # Determine if trimolecular families are present
        for family in self.database.kinetics.families.values():
            if len(family.forwardTemplate.reactants) > 2:
                logging.info('Trimolecular reactions are turned on')
                self.trimolecular = True
                break
        # Only check products if we want to react them
        if not self.trimolecular and self.trimolecularProductReversible:
            for family in self.database.kinetics.families.values():
                if len(family.forwardTemplate.products) > 2:
                    logging.info('Trimolecular reactions are turned on')
                    self.trimolecular = True
                    break

        # check libraries
        self.checkLibraries()

        if self.bindingEnergies:
            self.database.thermo.set_delta_atomic_adsorption_energies(self.bindingEnergies)

        # set global variable solvent
        if self.solvent:
            global solvent
            solvent = self.solvent

        if self.kineticsEstimator == 'rate rules':
            if '!training' not in self.kineticsDepositories:
                logging.info('Adding rate rules from training set in kinetics families...')
                # Temporarily remove species constraints for the training reactions
                copy_species_constraints = copy.copy(self.speciesConstraints)
                self.speciesConstraints = {}
                for family in self.database.kinetics.families.values():
                    if not family.autoGenerated:
                        family.add_rules_from_training(thermoDatabase=self.database.thermo)

                    # If requested by the user, write a text file for each kinetics family detailing the source of each entry
                    if self.kineticsdatastore:
                        logging.info(
                            'Writing sources of kinetic entries in family {0} to text file'.format(family.label))
                        path = os.path.join(self.outputDirectory, 'kinetics_database', family.label + '.txt')
                        with open(path, 'w') as f:
                            for template_label, entries in family.rules.entries.items():
                                f.write("Template [{0}] uses the {1} following source(s):\n".format(template_label,
                                                                                                    str(len(entries))))
                                for entry_index, entry in enumerate(entries):
                                    f.write(str(entry_index+1) + ". " + entry.shortDesc + "\n" + entry.longDesc + "\n")
                                f.write('\n')
                            f.write('\n')

                self.speciesConstraints = copy_species_constraints
            else:
                logging.info('Training set explicitly not added to rate rules in kinetics families...')
            logging.info('Filling in rate rules in kinetics families by averaging...')
            for family in self.database.kinetics.families.values():
                if not family.autoGenerated:
                    family.fill_rules_by_averaging_up(verbose=self.verboseComments)

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

        # Read input file
        self.loadInput(self.inputFile)

        if kwargs.get('restart', ''):
            import rmgpy.rmg.input
            rmgpy.rmg.input.restartFromSeed(path=kwargs['restart'])

        # Check input file 
        self.checkInput()

        # Properly set filterReactions to initialize flags properly
        if len(self.modelSettingsList) > 0:
            self.filterReactions = self.modelSettingsList[0].filterReactions

        # Make output subdirectories
        util.makeOutputSubdirectory(self.outputDirectory, 'pdep')
        util.makeOutputSubdirectory(self.outputDirectory, 'solver')
        util.makeOutputSubdirectory(self.outputDirectory, 'kinetics_database')

        # Specifies if details of kinetic database entries should be stored according to user
        try:
            self.kineticsdatastore = kwargs['kineticsdatastore']
        except KeyError:
            self.kineticsdatastore = False

        global maxproc
        try:
            maxproc = kwargs['maxproc']
        except KeyError:
            pass

        if maxproc > psutil.cpu_count():
            raise ValueError("""Invalid input for user defined maximum number of processes {0}; 
            should be an integer and smaller or equal to your available number of 
            processors {1}""".format(maxproc, psutil.cpu_count()))

        # Load databases
        self.loadDatabase()

        # Load restart seed mechanism (if specified)
        if self.restart:
            # Copy the restart files to a separate folder so that the job does not overwrite it
            restart_dir = os.path.join(self.outputDirectory, 'previous_restart')
            core_restart = os.path.join(restart_dir, 'restart')
            edge_restart = os.path.join(restart_dir, 'restart_edge')
            filters_restart = os.path.join(restart_dir, 'filters')
            util.makeOutputSubdirectory(self.outputDirectory, 'previous_restart')
            shutil.copytree(self.coreSeedPath, core_restart)
            shutil.copytree(self.edgeSeedPath, edge_restart)
            os.mkdir(filters_restart)
            shutil.copyfile(self.filtersPath, os.path.join(filters_restart, 'filters.h5'))
            shutil.copyfile(self.speciesMapPath, os.path.join(filters_restart, 'species_map.yml'))

            # Load the seed mechanism to get the core and edge species
            self.database.kinetics.load_libraries(restart_dir, libraries=['restart', 'restart_edge'])
            self.seedMechanisms.append('restart')
            self.reactionLibraries.append(('restart_edge', False))

        # Set trimolecular reactant flags of reaction systems
        if self.trimolecular:
            for reactionSystem in self.reactionSystems:
                reactionSystem.trimolecular = True

        # Do all liquid-phase startup things:
        if self.solvent:
            solvent_data = self.database.solvation.get_solvent_data(self.solvent)
            diffusionLimiter.enable(solvent_data, self.database.solvation)
            logging.info("Setting solvent data for {0}".format(self.solvent))

            # Set solvent viscosity for reaction filtering
            for reactionSystem in self.reactionSystems:
                reactionSystem.viscosity = solvent_data.get_solvent_viscosity(reactionSystem.T.value_si)

        try:
            self.wallTime = kwargs['walltime']
        except KeyError:
            pass

        data = self.wallTime.split(':')
        if not len(data) == 4:
            raise ValueError('Invalid format for wall time {0}; should be DD:HH:MM:SS.'.format(self.wallTime))
        self.wallTime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3]) + 86400 * int(data[-4])

        # Initialize reaction model

        # Seed mechanisms: add species and reactions from seed mechanism
        # DON'T generate any more reactions for the seed species at this time
        for seedMechanism in self.seedMechanisms:
            self.reactionModel.addSeedMechanismToCore(seedMechanism, react=False)

        # Reaction libraries: add species and reactions from reaction library to the edge so
        # that RMG can find them if their rates are large enough
        for library, option in self.reactionLibraries:
            self.reactionModel.addReactionLibraryToEdge(library)

        # Also always add in a few bath gases (since RMG-Java does)
        for label, smiles in [('Ar', '[Ar]'), ('He', '[He]'), ('Ne', '[Ne]'), ('N2', 'N#N')]:
            molecule = Molecule().from_smiles(smiles)
            spec, is_new = self.reactionModel.makeNewSpecies(molecule, label=label, reactive=False)
            if is_new:
                self.initialSpecies.append(spec)

        # Perform species constraints and forbidden species checks on input species
        for spec in self.initialSpecies:
            if self.database.forbiddenStructures.is_molecule_forbidden(spec.molecule[0]):
                if 'allowed' in self.speciesConstraints and 'input species' in self.speciesConstraints['allowed']:
                    logging.warning('Input species {0} is globally forbidden.  It will behave as an inert unless found '
                                    'in a seed mechanism or reaction library.'.format(spec.label))
                else:
                    raise ForbiddenStructureException("Input species {0} is globally forbidden. You may explicitly "
                                                      "allow it, but it will remain inert unless found in a seed "
                                                      "mechanism or reaction library.".format(spec.label))
            if failsSpeciesConstraints(spec):
                if 'allowed' in self.speciesConstraints and 'input species' in self.speciesConstraints['allowed']:
                    self.speciesConstraints['explicitlyAllowedMolecules'].append(spec.molecule[0])
                else:
                    raise ForbiddenStructureException("Species constraints forbids input species {0}. Please "
                                                      "reformulate constraints, remove the species, or explicitly "
                                                      "allow it.".format(spec.label))

        # For liquidReactor, checks whether the solvent is listed as one of the initial species.
        if self.solvent:
            solvent_structure_list = self.database.solvation.get_solvent_structure(self.solvent)
            for spc in solvent_structure_list:
                self.database.solvation.check_solvent_in_initial_species(self, spc)

        # Check to see if user has input Singlet O2 into their input file or libraries
        # This constraint is special in that we only want to check it once in the input instead of every time a species is made
        if 'allowSingletO2' in self.speciesConstraints and self.speciesConstraints['allowSingletO2']:
            pass
        else:
            # Here we get a list of all species that from the user input
            all_inputted_species = [spec for spec in self.initialSpecies]
            # Because no iterations have taken place, the only things in the core are from seed mechanisms
            all_inputted_species.extend(self.reactionModel.core.species)
            # Because no iterations have taken place, the only things in the edge are from reaction libraries
            all_inputted_species.extend(self.reactionModel.edge.species)

            O2Singlet = Molecule().from_smiles('O=O')
            for spec in all_inputted_species:
                if spec.is_isomorphic(O2Singlet):
                    raise ForbiddenStructureException("Species constraints forbids input species {0} RMG expects the "
                                                      "triplet form of oxygen for correct usage in reaction families. "
                                                      "Please change your input to SMILES='[O][O]' If you actually "
                                                      "want to use the singlet state, set the allowSingletO2=True "
                                                      "inside of the Species Constraints block in your input file."
                                                      .format(spec.label))

        for spec in self.initialSpecies:
            submit(spec, self.solvent)

        # Add nonreactive species (e.g. bath gases) to core first
        # This is necessary so that the PDep algorithm can identify the bath gas
        for spec in self.initialSpecies:
            if not spec.reactive:
                self.reactionModel.enlarge(spec)
        for spec in self.initialSpecies:
            if spec.reactive:
                self.reactionModel.enlarge(spec)

        # chatelak: store constant SPC indices in the reactor attributes if any constant SPC provided in the input file
        # advantages to write it here: this is run only once (as species indexes does not change over the generation)
        if self.solvent is not None:
            for index, reactionSystem in enumerate(self.reactionSystems):
                if reactionSystem.constSPCNames is not None:  # if no constant species provided do nothing
                    reactionSystem.get_constSPCIndices(
                        self.reactionModel.core.species)  # call the function to identify indices in the solver

        self.initializeReactionThresholdAndReactFlags()
        self.reactionModel.initializeIndexSpeciesDict()

    def register_listeners(self):
        """
        Attaches listener classes depending on the options 
        found in the RMG input file.
        """

        self.attach(ChemkinWriter(self.outputDirectory))
        self.attach(RMSWriter(self.outputDirectory))

        if self.generateOutputHTML:
            self.attach(OutputHTMLWriter(self.outputDirectory))

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
        self.Tmin = min([x.Trange[0].value_si if x.Trange else x.T.value_si for x in self.reactionSystems])
        self.Tmax = max([x.Trange[1].value_si if x.Trange else x.T.value_si for x in self.reactionSystems])
        try:
            self.Pmin = min([x.Prange[0].value_si if x.Prange else x.P.value_si for x in self.reactionSystems])
            self.Pmax = max([x.Prange[1].value_si if x.Prange else x.P.value_si for x in self.reactionSystems])
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
                    conditions=self.rmg_memories[index].get_cond(),
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

        max_num_spcs_hit = False  # default

        for q, modelSettings in enumerate(self.modelSettingsList):
            if len(self.simulatorSettingsList) > 1:
                simulator_settings = self.simulatorSettingsList[q]
            else:  # if they only provide one input for simulator use that everytime
                simulator_settings = self.simulatorSettingsList[0]

            self.filterReactions = modelSettings.filterReactions

            logging.info('Beginning model generation stage {0}...\n'.format(q + 1))

            self.done = False

            # Main RMG loop
            while not self.done:
                if self.generateSeedEachIteration:
                    self.makeSeedMech()

                self.reactionModel.iterationNum += 1
                self.done = True

                all_terminated = True
                num_core_species = len(self.reactionModel.core.species)

                prunable_species = self.reactionModel.edge.species[:]
                prunable_networks = self.reactionModel.networkList[:]

                for index, reactionSystem in enumerate(self.reactionSystems):

                    reactionSystem.prunableSpecies = prunable_species  # these lines reset pruning for a new cycle
                    reactionSystem.prunableNetworks = prunable_networks
                    reactionSystem.reset_max_edge_species_rate_ratios()

                    for p in range(reactionSystem.nSims):
                        reactor_done = True
                        objects_to_enlarge = []
                        self.reactionSystem = reactionSystem
                        # Conduct simulation
                        logging.info('Conducting simulation of reaction system %s...' % (index + 1))
                        prune = True

                        self.reactionModel.adjustSurface()

                        if num_core_species < modelSettings.minCoreSizeForPrune:
                            # Turn pruning off if we haven't reached minimum core size.
                            prune = False

                        try:
                            terminated, resurrected, obj, new_surface_species, new_surface_reactions, t, x = reactionSystem.simulate(
                                coreSpecies=self.reactionModel.core.species,
                                coreReactions=self.reactionModel.core.reactions,
                                edgeSpecies=self.reactionModel.edge.species,
                                edgeReactions=self.reactionModel.edge.reactions,
                                surfaceSpecies=self.reactionModel.surface.species,
                                surfaceReactions=self.reactionModel.surface.reactions,
                                pdepNetworks=self.reactionModel.networkList,
                                prune=prune,
                                modelSettings=modelSettings,
                                simulatorSettings=simulator_settings,
                                conditions=self.rmg_memories[index].get_cond()
                            )
                        except:
                            logging.error("Model core reactions:")
                            if len(self.reactionModel.core.reactions) > 5:
                                logging.error("Too many to print in detail")
                            else:
                                from arkane.output import prettify
                                logging.error(prettify(repr(self.reactionModel.core.reactions)))
                            if not self.generateSeedEachIteration:  # Then we haven't saved the seed mechanism yet
                                self.makeSeedMech(firstTime=True)  # Just in case the user wants to restart from this
                            raise

                        self.rmg_memories[index].add_t_conv_N(t, x, len(obj))
                        self.rmg_memories[index].generate_cond()
                        log_conditions(self.rmg_memories, index)

                        reactor_done = self.reactionModel.addNewSurfaceObjects(obj, new_surface_species,
                                                                              new_surface_reactions, reactionSystem)

                        all_terminated = all_terminated and terminated
                        logging.info('')

                        # If simulation is invalid, note which species should be added to
                        # the core
                        if obj != [] and not (obj is None):
                            objects_to_enlarge = self.processToSpeciesNetworks(obj)

                            reactor_done = False
                        # Enlarge objects identified by the simulation for enlarging
                        # These should be Species or Network objects
                        logging.info('')

                        objects_to_enlarge = list(set(objects_to_enlarge))

                        # Add objects to enlarge to the core first
                        for objectToEnlarge in objects_to_enlarge:
                            self.reactionModel.enlarge(objectToEnlarge)

                        if modelSettings.filterReactions:
                            # Run a raw simulation to get updated reaction system threshold values
                            # Run with the same conditions as with pruning off
                            temp_model_settings = deepcopy(modelSettings)
                            temp_model_settings.fluxToleranceKeepInEdge = 0
                            if not resurrected:
                                try:
                                    reactionSystem.simulate(
                                        coreSpecies=self.reactionModel.core.species,
                                        coreReactions=self.reactionModel.core.reactions,
                                        edgeSpecies=[],
                                        edgeReactions=[],
                                        surfaceSpecies=self.reactionModel.surface.species,
                                        surfaceReactions=self.reactionModel.surface.reactions,
                                        pdepNetworks=self.reactionModel.networkList,
                                        modelSettings=temp_model_settings,
                                        simulatorSettings=simulator_settings,
                                        conditions=self.rmg_memories[index].get_cond()
                                    )
                                except:
                                    self.updateReactionThresholdAndReactFlags(
                                        rxnSysUnimolecularThreshold=reactionSystem.unimolecularThreshold,
                                        rxnSysBimolecularThreshold=reactionSystem.bimolecularThreshold,
                                        rxnSysTrimolecularThreshold=reactionSystem.trimolecularThreshold,
                                        skipUpdate=True)
                                    logging.warning('Reaction thresholds/flags for Reaction System {0} was not updated '
                                                    'due to simulation failure'.format(index + 1))
                                else:
                                    self.updateReactionThresholdAndReactFlags(
                                        rxnSysUnimolecularThreshold=reactionSystem.unimolecularThreshold,
                                        rxnSysBimolecularThreshold=reactionSystem.bimolecularThreshold,
                                        rxnSysTrimolecularThreshold=reactionSystem.trimolecularThreshold
                                    )
                            else:
                                self.updateReactionThresholdAndReactFlags(
                                    rxnSysUnimolecularThreshold=reactionSystem.unimolecularThreshold,
                                    rxnSysBimolecularThreshold=reactionSystem.bimolecularThreshold,
                                    rxnSysTrimolecularThreshold=reactionSystem.trimolecularThreshold,
                                    skipUpdate=True
                                )
                                logging.warning('Reaction thresholds/flags for Reaction System {0} was not updated due '
                                                'to resurrection'.format(index + 1))

                            logging.info('')
                        else:
                            self.updateReactionThresholdAndReactFlags()

                        if not np.isinf(modelSettings.toleranceThermoKeepSpeciesInEdge):
                            self.reactionModel.setThermodynamicFilteringParameters(self.Tmax,
                                                                                   toleranceThermoKeepSpeciesInEdge=modelSettings.toleranceThermoKeepSpeciesInEdge,
                                                                                   minCoreSizeForPrune=modelSettings.minCoreSizeForPrune,
                                                                                   maximumEdgeSpecies=modelSettings.maximumEdgeSpecies,
                                                                                   reactionSystems=self.reactionSystems
                            )

                        old_edge_size = len(self.reactionModel.edge.reactions)
                        old_core_size = len(self.reactionModel.core.reactions)
                        self.reactionModel.enlarge(reactEdge=True,
                                                   unimolecularReact=self.unimolecularReact,
                                                   bimolecularReact=self.bimolecularReact,
                                                   trimolecularReact=self.trimolecularReact)

                        if old_edge_size != len(self.reactionModel.edge.reactions) or old_core_size != len(
                                self.reactionModel.core.reactions):
                            reactor_done = False

                        if not np.isinf(self.modelSettingsList[0].toleranceThermoKeepSpeciesInEdge):
                            self.reactionModel.thermoFilterDown(maximumEdgeSpecies=modelSettings.maximumEdgeSpecies)

                        max_num_spcs_hit = len(self.reactionModel.core.species) >= modelSettings.maxNumSpecies

                        self.saveEverything()

                        if max_num_spcs_hit:  # breaks the nSims loop
                            # self.done is still True, which will break the while loop
                            break

                        if not reactor_done:
                            self.done = False

                    if max_num_spcs_hit:  # breaks the reactionSystems loop
                        break

                if not self.done:  # There is something that needs exploring/enlarging

                    # If we reached our termination conditions, then try to prune
                    # species from the edge
                    if all_terminated and modelSettings.fluxToleranceKeepInEdge > 0.0:
                        logging.info('Attempting to prune...')
                        self.reactionModel.prune(self.reactionSystems, modelSettings.fluxToleranceKeepInEdge,
                                                 modelSettings.fluxToleranceMoveToCore,
                                                 modelSettings.maximumEdgeSpecies,
                                                 modelSettings.minSpeciesExistIterationsForPrune)
                        # Perform garbage collection after pruning
                        collected = gc.collect()
                        logging.info('Garbage collector: collected %d objects.' % collected)

                # Consider stopping gracefully if the next iteration might take us
                # past the wall time
                if self.wallTime > 0 and len(self.execTime) > 1:
                    t = self.execTime[-1]
                    dt = self.execTime[-1] - self.execTime[-2]
                    if t + 3 * dt > self.wallTime:
                        logging.info('MODEL GENERATION TERMINATED')
                        logging.info('')
                        logging.info(
                            'There is not enough time to complete the next iteration before the wall time is reached.')
                        logging.info('The output model may be incomplete.')
                        logging.info('')
                        core_spec, core_reac, edge_spec, edge_reac = self.reactionModel.getModelSize()
                        logging.info('The current model core has %s species and %s reactions' % (core_spec, core_reac))
                        logging.info('The current model edge has %s species and %s reactions' % (edge_spec, edge_reac))
                        return

            if max_num_spcs_hit:  # resets maxNumSpcsHit and continues the settings for loop
                logging.info('The maximum number of species ({0}) has been hit, Exiting stage {1} ...'.format(
                    modelSettings.maxNumSpecies, q + 1))
                max_num_spcs_hit = False
                continue

        # Save the final seed mechanism
        if self.generateSeedEachIteration:
            self.makeSeedMech()
        else:
            self.makeSeedMech(firstTime=True)

        self.run_model_analysis()

        # generate Cantera files chem.cti & chem_annotated.cti in a designated `cantera` output folder
        try:
            if any([s.contains_surface_site() for s in self.reactionModel.core.species]):
                self.generateCanteraFiles(os.path.join(self.outputDirectory, 'chemkin', 'chem-gas.inp'),
                                          surfaceFile=(
                                              os.path.join(self.outputDirectory, 'chemkin', 'chem-surface.inp')))
                self.generateCanteraFiles(os.path.join(self.outputDirectory, 'chemkin', 'chem_annotated-gas.inp'),
                                          surfaceFile=(os.path.join(self.outputDirectory, 'chemkin',
                                                                    'chem_annotated-surface.inp')))
            else:  # gas phase only
                self.generateCanteraFiles(os.path.join(self.outputDirectory, 'chemkin', 'chem.inp'))
                self.generateCanteraFiles(os.path.join(self.outputDirectory, 'chemkin', 'chem_annotated.inp'))
        except EnvironmentError:
            logging.exception('Could not generate Cantera files due to EnvironmentError. Check read\write privileges '
                              'in output directory.')
        except Exception:
            logging.exception('Could not generate Cantera files for some reason.')

        self.check_model()
        # Write output file
        logging.info('')
        logging.info('MODEL GENERATION COMPLETED')
        logging.info('')
        core_spec, core_reac, edge_spec, edge_reac = self.reactionModel.getModelSize()
        logging.info('The final model core has %s species and %s reactions' % (core_spec, core_reac))
        logging.info('The final model edge has %s species and %s reactions' % (edge_spec, edge_reac))

        self.finish()

    def run_model_analysis(self, number=10):
        """
        Run sensitivity and uncertainty analysis if requested.
        """
        # Run sensitivity analysis post-model generation if sensitivity analysis is on
        for index, reactionSystem in enumerate(self.reactionSystems):

            if reactionSystem.sensitiveSpecies and reactionSystem.sensConditions:
                logging.info('Conducting sensitivity analysis of reaction system {0}...'.format(index + 1))

                if reactionSystem.sensitiveSpecies == ['all']:
                    reactionSystem.sensitiveSpecies = self.reactionModel.core.species

                sens_worksheet = []
                for spec in reactionSystem.sensitiveSpecies:
                    csvfile_path = os.path.join(self.outputDirectory, 'solver',
                                               'sensitivity_{0}_SPC_{1}.csv'.format(index + 1, spec.index))
                    sens_worksheet.append(csvfile_path)

                terminated, resurrected, obj, surface_species, surface_reactions, t, x = reactionSystem.simulate(
                    coreSpecies=self.reactionModel.core.species,
                    coreReactions=self.reactionModel.core.reactions,
                    edgeSpecies=self.reactionModel.edge.species,
                    edgeReactions=self.reactionModel.edge.reactions,
                    surfaceSpecies=[],
                    surfaceReactions=[],
                    pdepNetworks=self.reactionModel.networkList,
                    sensitivity=True,
                    sensWorksheet=sens_worksheet,
                    modelSettings=ModelSettings(toleranceMoveToCore=1e8, toleranceInterruptSimulation=1e8),
                    simulatorSettings=self.simulatorSettingsList[-1],
                    conditions=reactionSystem.sensConditions,
                )

                plot_sensitivity(self.outputDirectory, index, reactionSystem.sensitiveSpecies, number=number)

        self.run_uncertainty_analysis()

    def run_uncertainty_analysis(self):
        """
        Run uncertainty analysis if proper settings are available.
        """
        if self.uncertainty is not None and self.uncertainty['global']:
            try:
                import libmuqModelling
            except ImportError:
                logging.error('Unable to import MUQ. Skipping global uncertainty analysis.')
                self.uncertainty['global'] = False
            else:
                import re
                import random
                from rmgpy.tools.canteraModel import Cantera
                from rmgpy.tools.muq import ReactorPCEFactory

        if self.uncertainty is not None and self.uncertainty['local']:
            correlation = []
            if self.uncertainty['uncorrelated']: correlation.append(False)
            if self.uncertainty['correlated']: correlation.append(True)

            # Set up Uncertainty object
            uncertainty = Uncertainty(outputDirectory=self.outputDirectory)
            uncertainty.database = self.database
            uncertainty.speciesList, uncertainty.reactionList = self.reactionModel.core.species, self.reactionModel.core.reactions
            uncertainty.extractSourcesFromModel()

            # Reload reaction families with verbose comments if necessary
            if not self.verboseComments:
                logging.info('Reloading kinetics families with verbose comments for uncertainty analysis...')
                self.database.kinetics.load_families(os.path.join(self.databaseDirectory, 'kinetics', 'families'),
                                                     self.kineticsFamilies, self.kineticsDepositories)
                # Temporarily remove species constraints for the training reactions
                self.speciesConstraints, speciesConstraintsCopy = {}, self.speciesConstraints
                for family in self.database.kinetics.families.values():
                    family.add_rules_from_training(thermoDatabase=self.database.thermo)
                    family.fill_rules_by_averaging_up(verbose=True)
                self.speciesConstraints = speciesConstraintsCopy

            for correlated in correlation:
                uncertainty.assignParameterUncertainties(correlated=correlated)

                for index, reactionSystem in enumerate(self.reactionSystems):
                    if reactionSystem.sensitiveSpecies and reactionSystem.sensConditions:
                        logging.info('Conducting {0}correlated local uncertainty analysis for '
                                     'reaction system {1}...\n'.format('un' if not correlated else '', index + 1))
                        results = uncertainty.localAnalysis(reactionSystem.sensitiveSpecies,
                                                            reactionSystemIndex=index,
                                                            correlated=correlated,
                                                            number=self.uncertainty['localnum'])
                        logging.info('Local uncertainty analysis results for reaction system {0}:\n'.format(index + 1))
                        local_result, local_result_str = process_local_results(results, reactionSystem.sensitiveSpecies,
                                                                               number=self.uncertainty['localnum'])
                        logging.info(local_result_str)

                        if self.uncertainty['global']:
                            logging.info('Conducting {0}correlated global uncertainty analysis for '
                                         'reaction system {1}...'.format('un' if not correlated else '', index + 1))
                            # Get simulation conditions
                            for criteria in reactionSystem.termination:
                                if isinstance(criteria, TerminationTime):
                                    time_criteria = ([criteria.time.value], criteria.time.units)
                                    break
                            else:
                                time_criteria = self.uncertainty['time']
                            Tlist = ([reactionSystem.sensConditions['T']], 'K')
                            Plist = ([reactionSystem.sensConditions['P']], 'Pa')
                            molFracList = [reactionSystem.sensConditions.copy()]
                            del molFracList[0]['T']
                            del molFracList[0]['P']

                            # Set up Cantera reactor
                            job = Cantera(speciesList=uncertainty.speciesList, reactionList=uncertainty.reactionList,
                                          outputDirectory=os.path.join(self.outputDirectory, 'global_uncertainty'))
                            job.loadModel()
                            job.generateConditions(
                                reactorTypeList=['IdealGasConstPressureTemperatureReactor'],
                                reactionTimeList=time_criteria,
                                molFracList=molFracList,
                                Tlist=Tlist,
                                Plist=Plist,
                            )

                            # Extract uncertain parameters from local analysis
                            k_params = []
                            g_params = []
                            for spc in reactionSystem.sensitiveSpecies:
                                _, reaction_c, thermo_c = local_result[spc]
                                for label, _, _ in reaction_c[:self.uncertainty['globalnum']]:
                                    if correlated:
                                        k_param = label
                                    else:
                                        # For uncorrelated, we need the reaction index
                                        k_index = label.split(':')[0]  # Looks like 'k1234: A+B=C+D'
                                        k_param = int(k_index[1:])
                                    if k_param not in k_params:
                                        k_params.append(k_param)
                                for label, _, _ in thermo_c[:self.uncertainty['globalnum']]:
                                    if correlated:
                                        g_param = label
                                    else:
                                        # For uncorrelated, we need the species index
                                        match = re.search(r'dG\[\S+\((\d+)\)\]', label)
                                        g_param = int(match.group(1))
                                    if g_param not in g_params:
                                        g_params.append(g_param)

                            reactor_pce_factory = ReactorPCEFactory(
                                cantera=job,
                                outputSpeciesList=reactionSystem.sensitiveSpecies,
                                kParams=k_params,
                                kUncertainty=uncertainty.kineticInputUncertainties,
                                gParams=g_params,
                                gUncertainty=uncertainty.thermoInputUncertainties,
                                correlated=correlated,
                                logx=self.uncertainty['logx'],
                            )

                            logging.info('Generating PCEs...')
                            reactor_pce_factory.generatePCE(runTime=self.uncertainty['pcetime'])

                            # Try a test point to see how well the PCE performs
                            reactor_pce_factory.compareOutput(
                                [random.uniform(-1.0, 1.0) for i in range(len(k_params) + len(g_params))])

                            # Analyze results and save statistics
                            reactor_pce_factory.analyzeResults()
                    else:
                        logging.info('Unable to run uncertainty analysis. Must specify sensitivity analysis options in '
                                     'reactor options.')

    def check_model(self):
        """
        Run checks on the RMG model
        """
        logging.info('Performing final model checks...')

        # Check that no two species in core or edge are isomorphic
        for i, spc in enumerate(self.reactionModel.core.species):
            for j in range(i):
                spc2 = self.reactionModel.core.species[j]
                if spc.is_isomorphic(spc2):
                    raise CoreError(
                        'Although the model has completed, species {0} is isomorphic to species {1} in the core. '
                        'Please open an issue on GitHub with the following output:'
                        '\n{2}\n{3}'.format(spc.label, spc2.label, spc.to_adjacency_list(), spc2.to_adjacency_list())
                    )

        for i, spc in enumerate(self.reactionModel.edge.species):
            for j in range(i):
                spc2 = self.reactionModel.edge.species[j]
                if spc.is_isomorphic(spc2):
                    logging.warning(
                        'Species {0} is isomorphic to species {1} in the edge. This does not affect '
                        'the generated model. If you would like to report this to help make RMG better '
                        'please open a GitHub issue with the following output:'
                        '\n{2}\n{3}'.format(spc.label, spc2.label, spc.to_adjacency_list(), spc2.to_adjacency_list())
                    )

        # Check all core reactions (in both directions) for collision limit violation
        violators = []
        num_rxn_violators = 0
        for rxn in self.reactionModel.core.reactions:
            if rxn.is_surface_reaction():
                # Don't check collision limits for surface reactions.
                continue
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
                    comment = ''
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
                                      'Violation condition: {5}\n\n'.format(rxn_string, kinetics, comment, direction,
                                                                            ratio, condition))
                    if isinstance(violator[0], TemplateReaction):
                        # although this is the end of the run, restore the original comment
                        violator[0].kinetics.comment = comment
        else:
            logging.info("No collision rate violators found.")

    def makeSeedMech(self, firstTime=False):
        """
        causes RMG to make a seed mechanism out of the current chem_annotated.inp and species_dictionary.txt
        this seed mechanism is outputted in a seed folder within the run directory and automatically
        added to as the (or replaces the current) 'Seed' thermo and kinetics libraries in database
        
        if run with firstTime=True it will change self.name to be unique within the thermo/kinetics libraries
        by adding integers to the end of the name to prevent overwritting

        This also writes the filter tensors to the `filters` sub-folder for restarting an RMG job from a seed mechanism
        """

        logging.info('Making seed mechanism...')

        name = self.name

        if self.saveSeedToDatabase and firstTime:  # make sure don't overwrite current libraries
            thermo_names = list(self.database.thermo.libraries.keys())
            kinetics_names = list(self.database.kinetics.libraries.keys())

            if name in thermo_names or name in kinetics_names:
                q = 1
                while name + str(q) in thermo_names or name + str(q) in kinetics_names:
                    q += 1
                self.name = name + str(q)

        seed_dir = os.path.join(self.outputDirectory, 'seed')
        filter_dir = os.path.join(seed_dir, 'filters')
        temp_seed_dir = os.path.join(self.outputDirectory, 'seed_tmp')

        if firstTime:
            if os.path.exists(seed_dir):  # This is a seed from a previous RMG run. Delete it
                shutil.rmtree(seed_dir)
        else:  # This is a seed from the previous iteration. Move it to a temporary directory in case we run into errors
            os.rename(seed_dir, os.path.join(temp_seed_dir))

        # Now that we have either deleted or moved the seed mechanism folder, create a new one
        os.mkdir(seed_dir)

        try:
            species_list = self.reactionModel.core.species
            reaction_list = self.reactionModel.core.reactions
            edge_species_list = self.reactionModel.edge.species
            edge_reaction_list = self.reactionModel.edge.reactions

            # Make species labels independent
            old_labels = self.makeSpeciesLabelsIndependent(species_list)
            edge_old_labels = self.makeSpeciesLabelsIndependent(edge_species_list)

            # load kinetics library entries
            kinetics_library = KineticsLibrary(name=name, autoGenerated=True)
            kinetics_library.entries = {}
            for i in range(len(reaction_list)):
                reaction = reaction_list[i]
                entry = Entry(
                    index=i + 1,
                    label=reaction.toLabeledStr(),
                    item=reaction,
                    data=reaction.kinetics,
                )

                if 'rate rule' in reaction.kinetics.comment:
                    entry.longDesc = reaction.kinetics.comment
                elif hasattr(reaction, 'library') and reaction.library:
                    entry.longDesc = 'Originally from reaction library: ' + reaction.library + "\n" + reaction.kinetics.comment
                else:
                    entry.longDesc = reaction.kinetics.comment

                kinetics_library.entries[i + 1] = entry

            # load kinetics library entries
            edge_kinetics_library = KineticsLibrary(name=name + '_edge', autoGenerated=True)
            edge_kinetics_library.entries = {}
            for i, reaction in enumerate(edge_reaction_list):
                entry = Entry(
                    index=i + 1,
                    label=reaction.toLabeledStr(),
                    item=reaction,
                    data=reaction.kinetics,
                )
                try:
                    entry.longDesc = 'Originally from reaction library: ' + reaction.library + "\n" + reaction.kinetics.comment
                except AttributeError:
                    entry.longDesc = reaction.kinetics.comment
                edge_kinetics_library.entries[i + 1] = entry

            # save in database
            if self.saveSeedToDatabase:
                database_directory = settings['database.directory']
                try:
                    os.makedirs(os.path.join(database_directory, 'kinetics', 'libraries', name))
                except:
                    pass
                kinetics_library.save(os.path.join(database_directory, 'kinetics', 'libraries', name, 'reactions.py'))
                kinetics_library.save_dictionary(
                    os.path.join(database_directory, 'kinetics', 'libraries', name, 'dictionary.txt'))

                try:
                    os.makedirs(os.path.join(database_directory, 'kinetics', 'libraries', name + '_edge'))
                except:
                    pass
                edge_kinetics_library.save(
                    os.path.join(database_directory, 'kinetics', 'libraries', name + '_edge', 'reactions.py'))
                edge_kinetics_library.save_dictionary(
                    os.path.join(database_directory, 'kinetics', 'libraries', name + '_edge', 'dictionary.txt'))

            # save in output directory
            # Rename for the output directory, as these names should not be dynamic
            kinetics_library.name = 'seed'
            kinetics_library.save(os.path.join(seed_dir, 'seed', 'reactions.py'))
            kinetics_library.save_dictionary(os.path.join(seed_dir, 'seed', 'dictionary.txt'))

            edge_kinetics_library.name = 'seed_edge'
            edge_kinetics_library.save(os.path.join(seed_dir, 'seed_edge', 'reactions.py'))
            edge_kinetics_library.save_dictionary(os.path.join(seed_dir, 'seed_edge', 'dictionary.txt'))

            # Save the filter tensors
            if not os.path.exists(filter_dir):
                os.mkdir(filter_dir)
            with h5py.File(os.path.join(filter_dir, 'filters.h5'), 'w') as f:
                if self.unimolecularThreshold is not None:
                    f.create_dataset('unimolecularThreshold', data=self.unimolecularThreshold)
                if self.bimolecularThreshold is not None:
                    f.create_dataset('bimolecularThreshold', data=self.bimolecularThreshold)
                if self.trimolecularThreshold is not None:
                    f.create_dataset('trimolecularThreshold', data=self.trimolecularThreshold)

            # Save a map of species indices
            spcs_map = [spc.molecule[0].to_adjacency_list() for spc in self.reactionModel.core.species]

            with open(os.path.join(filter_dir, 'species_map.yml'), 'w') as f:
                yaml.dump(data=spcs_map, stream=f)

            # Generate a file for restarting from a seed mechanism if this is not a restart job
            if firstTime and (not self.restart):
                with open(os.path.join(self.outputDirectory, 'restart_from_seed.py'), 'w') as f:
                    f.write('restartFromSeed(path=\'seed\')\n\n')
                    with open(self.inputFile, 'r') as inputFile:
                        f.write(''.join(inputFile.readlines()))

            # Finally, delete the seed mechanism from the previous iteration (if it exists)
            if os.path.exists(temp_seed_dir):
                shutil.rmtree(temp_seed_dir)

        except Exception as e:
            # Move the seed mechanism from the previous iteration (if it exists) back
            if os.path.exists(temp_seed_dir):
                shutil.rmtree(seed_dir)  # Delete the bad save of the current seed mechanism
                os.rename(temp_seed_dir, seed_dir)
                logging.error('Error in writing the seed mechanism for the current iteration. The seed mechanism from '
                              'the previous iteration has been restored')
            raise e

        # change labels back so species aren't renamed
        for i, label in enumerate(old_labels):
            species_list[i].label = label

        for i, label in enumerate(edge_old_labels):
            edge_species_list[i].label = label

    def makeSpeciesLabelsIndependent(self, species):
        """
        This method looks at the core species labels and makes sure none of them conflict
        If a conflict occurs, the second occurance will have '-2' added
        returns a list of the old labels
        """
        old_labels = []
        labels = set()
        for spec in species:
            old_labels.append(spec.label)
            duplicate_index = 1
            if '+' in spec.label:
                L = spec.molecule[0].get_formula()
            else:
                L = spec.label
            potential_label = L
            while potential_label in labels:
                duplicate_index += 1
                potential_label = L + '-{}'.format(duplicate_index)

            spec.label = potential_label
            labels.add(potential_label)

        return old_labels

    ################################################################################
    def processToSpeciesNetworks(self, obj):
        """
        breaks down the objects returned by simulate into Species and PDepNetwork
        components
        """

        if isinstance(obj, PDepNetwork):
            out = [self.processPdepNetworks(obj)]
            return out
        elif isinstance(obj, Species):
            return [obj]
        elif isinstance(obj, Reaction):
            return list(self.processReactionsToSpecies(obj))
        elif isinstance(obj, list):  # list of species
            rspcs = self.processReactionsToSpecies([k for k in obj if isinstance(k, Reaction)])
            spcs = {k for k in obj if isinstance(k, Species)} | rspcs
            nworks, pspcs = self.processPdepNetworks([k for k in obj if isinstance(k, PDepNetwork)])
            spcs = list(spcs - pspcs)  # avoid duplicate species
            return spcs + nworks
        else:
            raise TypeError("improper call, obj input was incorrect")

    def processPdepNetworks(self, obj):
        """
        properly processes PDepNetwork objects and lists of PDepNetwork objects returned from simulate
        """
        reaction_system = self.reactionSystem
        if isinstance(obj, PDepNetwork):
            # Determine which species in that network has the highest leak rate
            # We do this here because we need a temperature and pressure
            # Store the maximum leak species along with the associated network
            ob = (obj, obj.getMaximumLeakSpecies(reaction_system.T.value_si, reaction_system.P.value_si))
            return ob
        elif isinstance(obj, list):
            spcs = [ob.getMaximumLeakSpecies(reaction_system.T.value_si, reaction_system.P.value_si) for ob in obj]
            nworks = [(obj[i], spcs[i]) for i in range(len(obj))]
            return nworks, set(spcs)
        else:
            raise TypeError("improper call, obj input was incorrect")

    def processReactionsToSpecies(self, obj):
        """
        properly processes Reaction objects and lists of Reaction objects returned from simulate
        """
        core_species = self.reactionModel.core.species
        filter_fcn = lambda x: not ((x in core_species))  # remove species already in core
        if isinstance(obj, Reaction):
            potential_spcs = obj.reactants + obj.products
            potential_spcs = list(filter(filter_fcn, potential_spcs))
        elif isinstance(obj, list) or isinstance(obj, set):
            potential_spcs = set()
            for ob in obj:
                potential_spcs = potential_spcs | set(ob.reactants + ob.products)
            potential_spcs = {sp for sp in potential_spcs if filter_fcn(sp)}
        else:
            raise TypeError("improper call, obj input was incorrect")
        return potential_spcs

    def generateCanteraFiles(self, chemkinFile, **kwargs):
        """
        Convert a chemkin mechanism chem.inp file to a cantera mechanism file chem.cti
        and save it in the cantera directory
        """
        transport_file = os.path.join(os.path.dirname(chemkinFile), 'tran.dat')
        file_name = os.path.splitext(os.path.basename(chemkinFile))[0] + '.cti'
        out_name = os.path.join(self.outputDirectory, 'cantera', file_name)
        if 'surfaceFile' in kwargs:
            out_name = out_name.replace('-gas.', '.')
        cantera_dir = os.path.dirname(out_name)
        try:
            os.makedirs(cantera_dir)
        except OSError:
            if not os.path.isdir(cantera_dir):
                raise
        if os.path.exists(out_name):
            os.remove(out_name)
        parser = ck2cti.Parser()
        try:
            parser.convertMech(chemkinFile, transportFile=transport_file, outName=out_name, quiet=True, permissive=True,
                               **kwargs)
        except ck2cti.InputParseError:
            logging.exception("Error converting to Cantera format.")
            logging.info("Trying again without transport data file.")
            parser.convertMech(chemkinFile, outName=out_name, quiet=True, permissive=True, **kwargs)

    def initializeReactionThresholdAndReactFlags(self):
        num_core_species = len(self.reactionModel.core.species)

        # Initialize everything to react by default, but we will handle the restart and filtering case immediately after
        self.unimolecularReact = np.ones(num_core_species, bool)
        self.bimolecularReact = np.ones((num_core_species, num_core_species), bool)
        if self.trimolecular:
            self.trimolecularReact = np.ones((num_core_species, num_core_species, num_core_species), bool)

        if self.filterReactions or self.restart:  # Otherwise no need to initialize thresholds or fix react flags
            self.unimolecularThreshold = np.zeros(num_core_species, bool)
            self.bimolecularThreshold = np.zeros((num_core_species, num_core_species), bool)
            if self.trimolecular:
                self.trimolecularThreshold = np.zeros((num_core_species, num_core_species, num_core_species), bool)

            if self.restart:
                # Load in the restart mapping
                with open(os.path.join(self.speciesMapPath), 'r') as f:
                    restart_species_list = yaml.safe_load(stream=f)

                num_restart_spcs = len(restart_species_list)
                restart_species_list = [Species().from_adjacency_list(adj_list) for adj_list in restart_species_list]

                # Load in the restart filter tensors
                with h5py.File(self.filtersPath, 'r') as f:
                    try:
                        unimolecular_threshold_restart = f.get('unimolecularThreshold').value
                        bimolecular_threshold_restart = f.get('bimolecularThreshold').value
                        if self.trimolecular:
                            trimolecular_threshold_restart = f.get('trimolecularThreshold').value

                        # Expand Thresholds to match number of species in the current model.
                        # Note that we are about to reorder the core species to match the order in the restart seed
                        # mechanism, so we only need to broadcast to the indices up to numRestartSpcs. Any indices after
                        # this are additional species that should have `False` for their threshold
                        unimolecular_threshold = np.zeros(num_core_species, bool)
                        unimolecular_threshold[:num_restart_spcs] = unimolecular_threshold_restart

                        bimolecular_threshold = np.zeros((num_core_species, num_core_species), bool)
                        bimolecular_threshold[:num_restart_spcs, :num_restart_spcs] = bimolecular_threshold_restart

                        if self.trimolecular:
                            trimolecular_threshold = np.zeros((num_core_species, num_core_species, num_core_species), bool)
                            trimolecular_threshold[:num_restart_spcs, :num_restart_spcs, :num_restart_spcs] = \
                                trimolecular_threshold_restart

                        filters_found = True

                    except KeyError:  # If we can't find the filters then this is because filtering was not used
                        logging.warning('No filters were found in file {0}. This is to be expected if the restart '
                                        'files specified are from an RMG job without reaction filtering. Therefore, '
                                        'RMG will assume that all of the core species from the restart core seed have '
                                        'already been reacted and will not react them again. Additional species added '
                                        'to the input file but not in the restart core seed WILL be reacted with the '
                                        'rest of the core.'.format(self.filtersPath))

                        filters_found = False

                # Reorder the core species to match the indices of the restart filter tensors
                reordered_core_species = []
                for spc in restart_species_list:
                    for j, oldCoreSpc in enumerate(self.reactionModel.core.species):
                        if oldCoreSpc.is_isomorphic(spc, strict=False):
                            reordered_core_species.append(self.reactionModel.core.species.pop(j))
                            break
                    else:
                        raise RuntimeError('Species {0} was defined in the restart file, but was not included in the'
                                           'core.'.format(spc))

                # Append the remaining species left in the core to the very end
                self.reactionModel.core.species = reordered_core_species + self.reactionModel.core.species

                # If we are restarting we must be able to handle all four possible combinations of whether or not
                # filtering was used in the restart job or if filtering is used in the current job. To summarize:
                #
                # If filtering is being used in this job, we only need to set the threshold flags properly because the
                # react flags are updated based on this before the first enlarge. Otherwise, we need to explicitly set
                # the react flags properly because the thresholds will not be used.
                #
                # If filtering was used in the job we are restarting from, we can load these values in for all of the
                # species present in that job to help us set the current threshold or react flags.
                # Otherwise, if filtering was not used in the previous job then we can assume that all species
                # present in the job we are restarting from have already reacted, and thus either set their
                # filters to be True or their react flags to be false.

                # Fill in the filter tensors.
                if filters_found:  # Add in the filter data where we have it.
                    # Note that additional species that have been defined in the input file but were not
                    # present in the restart core seed will be reacted.

                    if self.filterReactions:  # Filling in the filter thresholds will suffice
                        # Fill in the newly initialized filter tensors
                        self.unimolecularThreshold = unimolecular_threshold
                        self.bimolecularThreshold = bimolecular_threshold
                        if self.trimolecular:
                            self.trimolecularThreshold = trimolecular_threshold

                    else:  # We must set the react flags instead. If it was `True` in the threshold, it should not react
                        self.unimolecularReact = np.logical_not(unimolecular_threshold)
                        self.bimolecularReact = np.logical_not(bimolecular_threshold)
                        if self.trimolecular:
                            self.trimolecularReact = np.logical_not(trimolecular_threshold)

                else:  # Assume that all species found in the restart core seed have already been reacted
                    if self.filterReactions:  # Filling in the filter thresholds will suffice
                        self.unimolecularThreshold[:num_restart_spcs] = True
                        self.bimolecularThreshold[:num_restart_spcs, :num_restart_spcs] = True
                        if self.trimolecular:
                            self.trimolecularThreshold[:num_restart_spcs, :num_restart_spcs, :num_restart_spcs] = True

                    else:  # We must set the react flags instead.
                        # Don't react any species that were present in the restart core seed
                        self.unimolecularReact[:num_restart_spcs] = False
                        self.bimolecularReact[:num_restart_spcs, :num_restart_spcs] = False
                        if self.trimolecular:
                            self.trimolecularReact[:num_restart_spcs, :num_restart_spcs, :num_restart_spcs] = False

    def updateReactionThresholdAndReactFlags(self,
                                             rxnSysUnimolecularThreshold=None,
                                             rxnSysBimolecularThreshold=None,
                                             rxnSysTrimolecularThreshold=None,
                                             skipUpdate=False):
        """
        updates the length and boolean value of the unimolecular and bimolecular react and threshold flags
        """
        num_core_species = len(self.reactionModel.core.species)
        prev_num_core_species = len(self.unimolecularReact)
        new_core_species = num_core_species > prev_num_core_species

        # Always reset the react arrays from prior iterations
        self.unimolecularReact = np.zeros((num_core_species), bool)
        self.bimolecularReact = np.zeros((num_core_species, num_core_species), bool)
        if self.trimolecular:
            self.trimolecularReact = np.zeros((num_core_species, num_core_species, num_core_species), bool)

        if self.filterReactions:
            if new_core_species:
                # Expand the threshold arrays if there were new core species added
                unimolecular_threshold = np.zeros((num_core_species), bool)
                bimolecular_threshold = np.zeros((num_core_species, num_core_species), bool)

                # Broadcast original thresholds
                unimolecular_threshold[:prev_num_core_species] = self.unimolecularThreshold
                bimolecular_threshold[:prev_num_core_species, :prev_num_core_species] = self.bimolecularThreshold
                self.unimolecularThreshold = unimolecular_threshold
                self.bimolecularThreshold = bimolecular_threshold

                if self.trimolecular:
                    trimolecular_threshold = np.zeros((num_core_species, num_core_species, num_core_species), bool)
                    trimolecular_threshold[:prev_num_core_species, :prev_num_core_species, :prev_num_core_species] = self.trimolecularThreshold
                    self.trimolecularThreshold = trimolecular_threshold

            if skipUpdate:
                return

            # Always update the react and threshold arrays
            for i in range(num_core_species):
                if not self.unimolecularThreshold[i] and rxnSysUnimolecularThreshold[i]:
                    # We've shifted from not reacting to reacting
                    self.unimolecularReact[i] = True
                    self.unimolecularThreshold[i] = True

            for i in range(num_core_species):
                for j in range(i, num_core_species):
                    if not self.bimolecularThreshold[i, j] and rxnSysBimolecularThreshold[i, j]:
                        # We've shifted from not reacting to reacting
                        self.bimolecularReact[i, j] = True
                        self.bimolecularThreshold[i, j] = True

            if self.trimolecular:
                for i in range(num_core_species):
                    for j in range(i, num_core_species):
                        for k in range(j, num_core_species):
                            if not self.trimolecularThreshold[i, j, k] and rxnSysTrimolecularThreshold[i, j, k]:
                                # We've shifted from not reacting to reacting
                                self.trimolecularReact[i, j, k] = True
                                self.trimolecularThreshold[i, j, k] = True
        else:
            # We are not filtering reactions
            if new_core_species:
                # React all the new core species unimolecularly
                for i in range(prev_num_core_species, num_core_species):
                    self.unimolecularReact[i] = True

                # React all the new core species with all the core species bimolecularly
                for i in range(num_core_species):
                    for j in range(prev_num_core_species, num_core_species):
                        self.bimolecularReact[i, j] = True

                # React all the new core species with all bimolecular combinations trimolecularly
                if self.trimolecular:
                    for i in range(num_core_species):
                        for j in range(num_core_species):
                            for k in range(prev_num_core_species, num_core_species):
                                self.trimolecularReact[i, j, k] = True

    def saveEverything(self):
        """
        Saves the output HTML and the Chemkin file
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
        # Print neural network-generated quote
        import datetime
        import textwrap
        try:
            from textgenrnn.quotes import get_quote
        except ImportError:
            pass
        else:
            quote = '"' + get_quote() + '"'
            logging.info('')
            logging.info(textwrap.fill(quote, subsequent_indent=' '))
            logging.info('             ---Quote-generating neural network, {}'.format(
                datetime.datetime.now().strftime("%B %Y")
            ))

        # Log end timestamp
        logging.info('')
        logging.info('RMG execution terminated at ' + time.asctime())

    def getGitCommit(self, modulePath):
        import subprocess
        if os.path.exists(os.path.join(modulePath, '..', '.git')):
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
            conda_package = get_condaPackage('rmg')
            if conda_package != '':
                logging.log(level, 'The current anaconda package for RMG-Py is:')
                logging.log(level, conda_package)
                logging.log(level, '')

        database_head, database_date = self.getGitCommit(settings['database.directory'])
        if database_head != '' and database_date != '':
            logging.log(level, 'The current git HEAD for RMG-database is:')
            logging.log(level, '\t%s' % database_head)
            logging.log(level, '\t%s' % database_date)
            logging.log(level, '')
        else:
            database_conda_package = get_condaPackage('rmgdatabase')
            if database_conda_package != '':
                logging.log(level, 'The current anaconda package for RMG-database is:')
                logging.log(level, database_conda_package)
                logging.log(level, '')

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

        T_list = []
        P_list = []
        concentration_list = []
        species_dict = {}
        termination = []

        with open(path, 'r') as f:
            line = self.readMeaningfulLineJava(f)
            while line != '':

                if line.startswith('TemperatureModel:'):
                    tokens = line.split()
                    units = tokens[2][1:-1]
                    assert units in ['C', 'F', 'K']
                    if units == 'C':
                        T_list = [float(T) + 273.15 for T in tokens[3:]]
                    elif units == 'F':
                        T_list = [(float(T) + 459.67) * 5. / 9. for T in tokens[3:]]
                    else:
                        T_list = [float(T) for T in tokens[3:]]

                elif line.startswith('PressureModel:'):
                    tokens = line.split()
                    units = tokens[2][1:-1]
                    assert units in ['atm', 'bar', 'Pa', 'torr']
                    if units == 'atm':
                        P_list = [float(P) * 101325. for P in tokens[3:]]
                    elif units == 'bar':
                        P_list = [float(P) * 100000. for P in tokens[3:]]
                    elif units == 'torr':
                        P_list = [float(P) / 760. * 101325. for P in tokens[3:]]
                    else:
                        P_list = [float(P) for P in tokens[3:]]

                elif line.startswith('InitialStatus:'):
                    label = ''
                    concentrations = []
                    adjlist = ''

                    line = self.readMeaningfulLineJava(f)
                    while line != 'END':

                        if line == '' and label != '':
                            species = Species(label=label, molecule=[Molecule().from_adjacency_list(adjlist)])
                            self.initialSpecies.append(species)
                            species_dict[label] = species
                            concentration_list.append(concentrations)
                            label = ''
                            concentrations = []
                            adjlist = ''

                        elif line != '' and label == '':
                            tokens = line.split()
                            label = tokens[0]
                            units = tokens[1][1:-1]
                            if tokens[-1] in ['Unreactive', 'ConstantConcentration']:
                                tokens.pop(-1)
                            assert units in ['mol/cm3', 'mol/m3', 'mol/l']
                            if units == 'mol/cm3':
                                concentrations = [float(C) * 1.0e6 for C in tokens[2:]]
                            elif units == 'mol/l':
                                concentrations = [float(C) * 1.0e3 for C in tokens[2:]]
                            else:
                                concentrations = [float(C) for C in tokens[2:]]

                        elif line != '':
                            adjlist += line + '\n'

                        line = f.readline().strip()
                        if '//' in line:
                            line = line[0:line.index('//')]

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
                            concentrations = [float(C) * 1.0e6 for C in tokens[2:]]
                        elif units == 'mol/l':
                            concentrations = [float(C) * 1.0e3 for C in tokens[2:]]
                        else:
                            concentrations = [float(C) for C in tokens[2:]]

                        species = Species(label=label, reactive=False, molecule=[Molecule().from_smiles(smiles)])
                        self.initialSpecies.append(species)
                        species_dict[label] = species
                        concentration_list.append(concentrations)

                        line = self.readMeaningfulLineJava(f)

                elif line.startswith('FinishController:'):

                    # First meaningful line is a termination time or conversion
                    line = self.readMeaningfulLineJava(f)
                    tokens = line.split()
                    if tokens[2].lower() == 'conversion:':
                        label = tokens[3]
                        conversion = float(tokens[4])
                        termination.append(TerminationConversion(spec=species_dict[label], conv=conversion))
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

        assert len(T_list) > 0
        assert len(P_list) > 0
        concentration_list = np.array(concentration_list)
        # An arbitrary number of concentrations is acceptable, and should be run for each reactor system
        if not concentration_list.shape[1] > 0:
            raise AssertionError()

        # Make a reaction system for each (T,P) combination
        for T in T_list:
            for P in P_list:
                for i in range(concentration_list.shape[1]):
                    concentrations = concentration_list[:, i]
                    total_conc = np.sum(concentrations)
                    initial_mole_fractions = dict([(self.initialSpecies[i], concentrations[i] / total_conc) for i in
                                                 range(len(self.initialSpecies))])
                    reaction_system = SimpleReactor(T, P, initialMoleFractions=initial_mole_fractions,
                                                   termination=termination)
                    self.reactionSystems.append(reaction_system)

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
            if '//' in line:
                line = line[0:line.index('//')]
            while line == '':
                line = f.readline()
                if line == '':
                    break
                line = line.strip()
                if '//' in line:
                    line = line[0:line.index('//')]
        return line


################################################################################

def determine_procnum_from_RAM():
    """
    Get available RAM (GB)and procnum dependent on OS.
    """
    if sys.platform.startswith('linux'):
        # linux
        memory_available = psutil.virtual_memory().free / (1000.0 ** 3)
        memory_use = psutil.Process(os.getpid()).memory_info()[0] / (1000.0 ** 3)
        tmp = divmod(memory_available, memory_use)
        tmp2 = min(maxproc, tmp[0])
        procnum = max(1, int(tmp2))
    elif sys.platform == "darwin":
        # OS X
        memory_available = psutil.virtual_memory().available / (1000.0 ** 3)
        memory_use = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1000.0 ** 3)
        tmp = divmod(memory_available, memory_use)
        tmp2 = min(maxproc, tmp[0])
        procnum = max(1, int(tmp2))
    else:
        # Everything else
        procnum = 1

    # Return the maximal number of processes for multiprocessing
    return procnum


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
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    # formatter = Formatter('%(message)s', '%Y-%m-%d %H:%M:%S')
    formatter = logging.Formatter('%(levelname)s%(message)s')
    ch.setFormatter(formatter)

    # create file handler
    if os.path.exists(log_file_name):
        name, ext = os.path.splitext(log_file_name)
        backup = name + '_backup' + ext
        if os.path.exists(backup):
            logging.info("Removing old " + backup)
            os.remove(backup)
        logging.info('Moving {0} to {1}\n'.format(log_file_name, backup))
        shutil.move(log_file_name, backup)
    fh = logging.FileHandler(filename=log_file_name)  # , backupCount=3)
    fh.setLevel(min(logging.DEBUG, verbose))  # always at least VERBOSE in the file
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
class RMG_Memory(object):
    """
    class for remembering RMG simulations
    and determining what simulation to run next
    """

    def __init__(self, reactionSystem, bspc):
        self.Ranges = dict()

        if hasattr(reactionSystem, 'Trange') and isinstance(reactionSystem.Trange, list):
            T_range = reactionSystem.Trange
            self.Ranges['T'] = [T.value_si for T in T_range]
        if hasattr(reactionSystem, 'Prange') and isinstance(reactionSystem.Prange, list):
            P_range = reactionSystem.Prange
            self.Ranges['P'] = [np.log(P.value_si) for P in P_range]
        if hasattr(reactionSystem, 'initialMoleFractions'):
            if bspc:
                self.initialMoleFractions = deepcopy(reactionSystem.initialMoleFractions)
                self.balanceSpecies = [x for x in self.initialMoleFractions.keys() if x.label == bspc][
                    0]  # find the balance species
            for key, value in reactionSystem.initialMoleFractions.items():
                assert key != 'T' and key != 'P', 'naming a species T or P is forbidden'
                if isinstance(value, list):
                    self.Ranges[key] = value
        if hasattr(reactionSystem, 'initialConcentrations'):
            for key, value in reactionSystem.initialConcentrations.items():
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

    def add_t_conv_N(self, t, conv, N):
        """
        adds the completion time and conversion and the number of objects added 
        from a given run to the memory
        """
        if hasattr(self, 'tmax'):
            self.ts.append(t / self.tmax)
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

    def calculate_cond(self, obj, Ndims, Ns=20):
        """
        Weighted Stochastic Grid Sampling algorithm 
        obj is evaluated at a grid of points and the evaluations are normalized 
        and then sampled randomly based on their normalized value
        then a random step of length 1/(2*Ns) is taken from that point to give a final condition point
        if this process were to impact runtime under some conditions you could decrease the value of Ns to speed it up
        """
        bounds = tuple((0.0, 1.0) for k in range(Ndims))
        x0, fval, grid, Jout = brute(obj, bounds, Ns=Ns, full_output=True,
                                     finish=None)  # run brute just to easily get the evaluations at each grid point (we don't care about the optimal value)
        Jout += abs(Jout.min(tuple(range(Ndims))))  # shifts Jout positive so tot is positive
        tot = np.sum(Jout, axis=tuple(range(len(Jout.shape))))
        Jout /= tot  # normalize Jout
        n = self.randState.uniform(0, 1, 1)[0]  # draw a random number between 0 and 1
        s = 0.0
        for indexes in np.ndenumerate(Jout):  # choose a coordinate such that grid[indexes] is choosen with probability Jout[indexes]
            s += Jout[indexes[0]]
            if s > n:
                break
        if len(bounds) != 1:
            yf = np.array([grid[i][indexes[0]] for i in range(len(grid))])
        else:
            yf = np.array([grid[indexes[0]] for i in range(len(grid))])

        step = self.randState.uniform(0, 1, len(Jout.shape))  # take a step in a random direction in a length between 0 and 1/(2*Ns)
        step /= step.sum()
        mag = self.randState.uniform(0, 1, 1)[0]

        yf += step * mag * np.sqrt(2) / (2.0 * Ns)

        return yf

    def generate_cond(self):
        """
        find the next condition to run at by solving an optimization problem
        this optimization problem maximizes distance from prior conditions weighted more if they are more recent
        and maximizes number of objects added
        the resulting condition is added to the end of conditionList
        """
        if self.conditionList == []:
            self.conditionList.append({key: value[0] for key, value in self.Ranges.items()})
            self.scaledConditionList.append({key: 0.0 for key, value in self.Ranges.items()})
        elif len(self.conditionList[0]) == 0:
            pass
        else:
            ykey = list(self.conditionList[0].keys())
            Ns = self.Ns

            def obj(y):
                boo = y.shape == tuple()
                vec = []
                N = len(self.conditionList)
                for i, cond in enumerate(self.scaledConditionList):
                    for j, key in enumerate(ykey):
                        if not boo:
                            vec.append(10.0 * N / ((N - i) * (Ns[i] + 1)) * abs(y[j] - cond[key]) ** 0.3)
                        else:
                            vec.append(10.0 * N / ((N - i) * (Ns[i] + 1)) * abs(y - cond[key]) ** 0.3)
                return -np.array(vec).sum()

            yf = self.calculate_cond(obj, len(ykey))

            scaled_new_cond = {ykey[i]: yf[i] for i in range(len(ykey))}
            new_cond = {yk: yf[i] * (self.Ranges[yk][1] - self.Ranges[yk][0]) + self.Ranges[yk][0] for i, yk in
                       enumerate(ykey)}
            if 'P' in list(new_cond.keys()):
                new_cond['P'] = np.exp(new_cond['P'])

            if hasattr(self, 'initialMoleFractions'):
                for key in self.initialMoleFractions.keys():
                    if not isinstance(self.initialMoleFractions[key], list):
                        new_cond[key] = self.initialMoleFractions[key]
                total = sum([val for key, val in new_cond.items() if key != 'T' and key != 'P'])
                if self.balanceSpecies is None:
                    for key, val in new_cond.items():
                        if key != 'T' and key != 'P':
                            new_cond[key] = val / total
                else:
                    new_cond[self.balanceSpecies] = self.initialMoleFractions[self.balanceSpecies] + 1.0 - total

            self.conditionList.append(new_cond)
            self.scaledConditionList.append(scaled_new_cond)
        return


def log_conditions(RMG_Memories, index):
    """
    log newly generated reactor conditions
    """
    if RMG_Memories[index].get_cond() is not None:
        s = 'conditions choosen for reactor {0} were: '.format(index)
        for key, item in RMG_Memories[index].get_cond().items():
            if key == 'T':
                s += 'T = {0} K, '.format(item)
            elif key == 'P':
                s += 'P = {0} bar, '.format(item / 1.0e5)
            else:
                s += key.label + ' = {0}, '.format(item)

        logging.info(s)


class Tee(object):
    """A simple tee to create a stream which prints to many streams.
    
    This is used to report the profiling statistics to both the log file
    and the standard output.
    """

    def __init__(self, *fileobjects):
        self.fileobjects = fileobjects

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

        packages = []
        # Strip comments
        for line in lines:
            if line[:1] == '#':
                pass
            else:
                packages.append(line)

        return '\n'.join(packages)
    except:
        return ''


def processProfileStats(stats_file, log_file):
    import pstats
    out_stream = Tee(sys.stdout, open(log_file, 'a'))  # print to screen AND append to RMG.log
    print("=" * 80, file=out_stream)
    print("Profiling Data".center(80), file=out_stream)
    print("=" * 80, file=out_stream)
    stats = pstats.Stats(stats_file, stream=out_stream)
    stats.strip_dirs()
    print("Sorted by internal time", file=out_stream)
    stats.sort_stats('time')
    stats.print_stats(25)
    stats.print_callers(25)
    print("Sorted by cumulative time", file=out_stream)
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

    # create an Options class to mimic optparser output as much as possible:
    class Options(object):
        pass

    options = Options()
    options.node_thres = 0.8
    options.edge_thres = 0.1
    options.strip = False
    options.show_samples = False
    options.root = ""
    options.leaf = ""
    options.wrap = True

    theme = themes['color']  # bw color gray pink
    theme.fontname = "ArialMT"  # default "Arial" leads to PostScript warnings in dot (on Mac OS)
    parser = PstatsParser(stats_file)
    profile = parser.parse()

    dot_file = stats_file + '.dot'
    output = open(dot_file, 'wt')
    dot = DotWriter(output)
    dot.strip = options.strip
    dot.wrap = options.wrap

    if options.show_samples:
        dot.show_function_events.append(SAMPLES)

    profile = profile
    profile.prune(options.node_thres / 100.0, options.edge_thres / 100.0)

    if options.root:
        root_id = profile.getFunctionId(options.root)
        if not root_id:
            sys.stderr.write('root node ' + options.root + ' not found (might already be pruned : try -E0 -n0 flags)\n')
            sys.exit(1)
        profile.prune_root(root_id)
    if options.leaf:
        leaf_id = profile.getFunctionId(options.leaf)
        if not leaf_id:
            sys.stderr.write('leaf node ' + options.leaf + ' not found (maybe already pruned : try -E0 -n0 flags)\n')
            sys.exit(1)
        profile.prune_leaf(leaf_id)

    dot.graph(profile, theme)

    output.close()

    try:
        subprocess.check_call(['dot', '-Tps2', dot_file, '-o', '{0}.ps2'.format(dot_file)])
    except subprocess.CalledProcessError:
        logging.error("Error returned by 'dot' when generating graph of the profile statistics.")
        logging.info("To try it yourself:\n     dot -Tps2 {0} -o {0}.ps2".format(dot_file))
    except OSError:
        logging.error("Couldn't run 'dot' to create graph of profile statistics. Check graphviz is installed properly "
                      "and on your path.")
        logging.info("Once you've got it, try:\n     dot -Tps2 {0} -o {0}.ps2".format(dot_file))

    try:
        subprocess.check_call(['ps2pdf', '{0}.ps2'.format(dot_file), '{0}.pdf'.format(dot_file)])
    except OSError:
        logging.error("Couldn't run 'ps2pdf' to create pdf graph of profile statistics. Check that ps2pdf converter "
                      "is installed.")
        logging.info("Once you've got it, try:\n     pd2pdf {0}.ps2 {0}.pdf".format(dot_file))
    else:
        logging.info("Graph of profile statistics saved to: \n {0}.pdf".format(dot_file))
