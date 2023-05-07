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

"""
This module contains the main execution functionality for Reaction Mechanism
Generator (RMG).
"""

import copy
import gc
import logging
import marshal
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
from cantera import ck2yaml
from scipy.optimize import brute

import rmgpy.util as util
from rmgpy.rmg.model import Species, CoreEdgeReactionModel
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy import settings
from rmgpy.chemkin import ChemkinWriter
from rmgpy.constraints import fails_species_constraints
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.exceptions import ForbiddenStructureException, DatabaseError, CoreError, InputError
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.data.vaporLiquidMassTransfer import vapor_liquid_mass_transfer
from rmgpy.kinetics import ThirdBody
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
from rmgpy.rmg.reactors import Reactor

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

    ========================================================== ================================================
    Attribute                                                  Description
    ========================================================== ================================================
    `input_file`                                               The path to the input file
    `profiler`                                                 A cProfile.Profile object for time profiling RMG
    ---------------------------------------------------------- ------------------------------------------------
    `database_directory`                                       The directory containing the RMG database
    `thermo_libraries`                                         The thermodynamics libraries to load
    `reaction_libraries`                                       The kinetics libraries to load
    `statmech_libraries`                                       The statistical mechanics libraries to load
    `seed_mechanisms`                                          The seed mechanisms included in the model
    `kinetics_families`                                        The kinetics families to use for reaction generation
    `kinetics_depositories`                                    The kinetics depositories to use for looking up kinetics in each family
    `kinetics_estimator`                                       The method to use to estimate kinetics: 'group additivity' or 'rate rules'
    `solvent`                                                  If solvation estimates are required, the name of the solvent.
    `liquid_volumetric_mass_transfer_coefficient_power_law`    If kLA estimates are required, the coefficients for kLA power law
    ---------------------------------------------------------- ------------------------------------------------
    `reaction_model`                                           The core-edge reaction model generated by this job
    `reaction_systems`                                         A list of the reaction systems used in this job
    `database`                                                 The RMG database used in this job
    ---------------------------------------------------------- ------------------------------------------------
    `model_settings_list`                                      List of ModelSettings objects containing information related to how to manage species/reaction movement
    `simulator_settings_list`                                  List of SimulatorSettings objects containing information on how to run simulations
    `init_react_tuples`                                        List of name tuples of species to react at beginning of run
    `trimolecular`                                             ``True`` to consider reactions between three species (i.e., if trimolecular reaction families are present)
    `unimolecular_threshold`                                   Array of flags indicating whether a species is above the unimolecular reaction threshold
    `bimolecular_threshold`                                    Array of flags indicating whether two species are above the bimolecular reaction threshold
    `trimolecular_threshold`                                   Array of flags indicating whether three species are above the trimolecular reaction threshold
    `unimolecular_react`                                       Array of flags indicating whether a species should react unimolecularly in the enlarge step
    `bimolecular_react`                                        Array of flags indicating whether two species should react in the enlarge step
    `trimolecular_react`                                       Array of flags indicating whether three species should react in the enlarge step
    `termination`                                              A list of termination targets (i.e :class:`TerminationTime` and :class:`TerminationConversion` objects)
    `species_constraints`                                      Dictates the maximum number of atoms, carbons, electrons, etc. generated by RMG
    ---------------------------------------------------------- ------------------------------------------------
    `output_directory`                                         The directory used to save output files
    `verbosity`                                                The level of logging verbosity for console output
    `units`                                                    The unit system to use to save output files (currently must be 'si')
    `generate_output_html`                                     ``True`` to draw pictures of the species and reactions, saving a visualized model in an output HTML file.  ``False`` otherwise
    `generate_plots`                                           ``True`` to generate plots of the job execution statistics after each iteration, ``False`` otherwise
    `verbose_comments`                                         ``True`` to keep the verbose comments for database estimates, ``False`` otherwise
    `save_edge_species`                                        ``True`` to save chemkin and HTML files of the edge species, ``False`` otherwise
    `keep_irreversible`                                        ``True`` to keep ireversibility of library reactions as is ('<=>' or '=>'). ``False`` (default) to force all library reactions to be reversible ('<=>')
    `trimolecular_product_reversible`                          ``True`` (default) to allow families with trimolecular products to react in the reverse direction, ``False`` otherwise
    `pressure_dependence`                                      Whether to process unimolecular (pressure-dependent) reaction networks
    `quantum_mechanics`                                        Whether to apply quantum mechanical calculations instead of group additivity to certain molecular types.
    `ml_estimator`                                             To use thermo estimation with machine learning
    `ml_settings`                                              Settings for ML estimation
    `walltime`                                                 The maximum amount of CPU time in the form DD:HH:MM:SS to expend on this job; used to stop gracefully so we can still get profiling information
    `max_iterations`                                           The maximum number of RMG iterations allowed, after which the job will terminate
    `kinetics_datastore`                                       ``True`` if storing details of each kinetic database entry in text file, ``False`` otherwise
    ---------------------------------------------------------- ------------------------------------------------
    `initialization_time`                                      The time at which the job was initiated, in seconds since the epoch (i.e. from time.time())
    `done`                                                     Whether the job has completed (there is nothing new to add)
    ========================================================== ================================================

    """

    def __init__(self, input_file=None, output_directory=None, profiler=None, stats_file=None):
        super(RMG, self).__init__()
        self.input_file = input_file
        self.output_directory = output_directory
        self.profiler = profiler
        self.clear()
        self.model_settings_list = []
        self.simulator_settings_list = []
        self.max_iterations = None
        self.Tmin = 0.0
        self.Tmax = 0.0
        self.Pmin = 0.0
        self.Pmax = 0.0
        self.database = None

    def clear(self):
        """
        Clear all loaded information about the job (except the file paths).
        """
        self.database_directory = None
        self.thermo_libraries = None
        self.transport_libraries = None
        self.reaction_libraries = None
        self.statmech_libraries = None
        self.seed_mechanisms = None
        self.kinetics_families = None
        self.kinetics_depositories = None
        self.kinetics_estimator = 'group additivity'
        self.solvent = None
        self.diffusion_limiter = None
        self.surface_site_density = None
        self.binding_energies = None
        self.coverage_dependence = False
        self.forbidden_structures = []

        self.reaction_model = None
        self.reaction_systems = None
        self.database = None
        self.reaction_system = None

        self.model_settings_list = []
        self.simulator_settings_list = []
        self.balance_species = None

        self.filter_reactions = False
        self.init_react_tuples = []
        self.trimolecular = False
        self.unimolecular_react = None
        self.bimolecular_react = None
        self.trimolecular_react = None
        self.unimolecular_threshold = None
        self.bimolecular_threshold = None
        self.trimolecular_threshold = None
        self.termination = []

        self.done = False
        self.verbosity = logging.INFO
        self.units = 'si'
        self.generate_output_html = None
        self.generate_plots = None
        self.save_simulation_profiles = None
        self.verbose_comments = None
        self.save_edge_species = None
        self.keep_irreversible = None
        self.trimolecular_product_reversible = None
        self.pressure_dependence = None
        self.quantum_mechanics = None
        self.ml_estimator = None
        self.ml_settings = None
        self.species_constraints = {}
        self.walltime = '00:00:00:00'
        self.save_seed_modulus = -1
        self.max_iterations = None
        self.initialization_time = 0
        self.kinetics_datastore = None
        self.restart = False
        self.core_seed_path = None
        self.edge_seed_path = None
        self.filters_path = None
        self.species_map_path = None

        self.name = 'Seed'
        self.generate_seed_each_iteration = True
        self.save_seed_to_database = False

        self.thermo_central_database = None
        self.uncertainty = None

        self.exec_time = []
        self.liquid_volumetric_mass_transfer_coefficient_power_law = None

    def load_input(self, path=None):
        """
        Load an RMG job from the input file located at `input_file`, or
        from the `input_file` attribute if not given as a parameter.
        """
        from rmgpy.rmg.input import read_input_file
        if path is None:
            path = self.input_file
        read_input_file(path, self)
        self.reaction_model.kinetics_estimator = self.kinetics_estimator
        # If the output directory is not yet set, then set it to the same
        # directory as the input file by default
        if not self.output_directory:
            self.output_directory = os.path.dirname(path)
        if self.pressure_dependence:
            self.pressure_dependence.output_file = self.output_directory
            self.reaction_model.pressure_dependence = self.pressure_dependence
        if self.solvent:
            self.reaction_model.solvent_name = self.solvent

        if self.surface_site_density:
            self.reaction_model.surface_site_density = self.surface_site_density
            self.reaction_model.core.phase_system.phases["Surface"].site_density = self.surface_site_density.value_si
            self.reaction_model.edge.phase_system.phases["Surface"].site_density = self.surface_site_density.value_si
        self.reaction_model.coverage_dependence = self.coverage_dependence
            
        self.reaction_model.verbose_comments = self.verbose_comments
        self.reaction_model.save_edge_species = self.save_edge_species

        if self.quantum_mechanics:
            self.reaction_model.quantum_mechanics = self.quantum_mechanics

        for reaction_system in self.reaction_systems:
            self.reaction_model.reaction_systems.append(reaction_system)

    def load_thermo_input(self, path=None):
        """
        Load an Thermo Estimation job from a thermo input file located at `input_file`, or
        from the `input_file` attribute if not given as a parameter.
        """
        from rmgpy.rmg.input import read_thermo_input_file
        if path is None:
            path = self.input_file
        if not self.output_directory:
            self.output_directory = os.path.dirname(path)
        read_thermo_input_file(path, self)

        if self.quantum_mechanics:
            self.reaction_model.quantum_mechanics = self.quantum_mechanics

    def check_input(self):
        """
        Check for a few common mistakes in the input file.
        """
        if self.pressure_dependence:
            for index, reaction_system in enumerate(self.reaction_systems):
                if reaction_system.T:
                    logging.info(reaction_system.T)
                    assert (reaction_system.T.value_si < self.pressure_dependence.Tmax.value_si), "Reaction system T is above pressure_dependence range."
                    assert (reaction_system.T.value_si > self.pressure_dependence.Tmin.value_si), "Reaction system T is below pressure_dependence range."
                else:
                    assert (reaction_system.Trange[1].value_si < self.pressure_dependence.Tmax.value_si), "Reaction system T is above pressure_dependence range."
                    assert (reaction_system.Trange[0].value_si > self.pressure_dependence.Tmin.value_si), "Reaction system T is below pressure_dependence range."
                if reaction_system.P:
                    assert (reaction_system.P.value_si < self.pressure_dependence.Pmax.value_si), "Reaction system P is above pressure_dependence range."
                    assert (reaction_system.P.value_si > self.pressure_dependence.Pmin.value_si), "Reaction system P is below pressure_dependence range."
                else:
                    assert (reaction_system.Prange[1].value_si < self.pressure_dependence.Pmax.value_si), "Reaction system P is above pressure_dependence range."
                    assert (reaction_system.Prange[0].value_si > self.pressure_dependence.Pmin.value_si), "Reaction system P is below pressure_dependence range."

            assert any([not s.reactive for s in reaction_system.initial_mole_fractions.keys()]), \
                "Pressure Dependence calculations require at least one inert (nonreacting) species for the bath gas."

    def check_libraries(self):
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

    def save_input(self, path=None):
        """
        Save an RMG job to the input file located at `path`.
        """
        from rmgpy.rmg.input import save_input_file
        save_input_file(path, self)

    def load_database(self):

        self.database = RMGDatabase()
        self.database.load(
            path=self.database_directory,
            thermo_libraries=self.thermo_libraries,
            transport_libraries=self.transport_libraries,
            reaction_libraries=[library for library, option in self.reaction_libraries],
            seed_mechanisms=self.seed_mechanisms,
            kinetics_families=self.kinetics_families,
            kinetics_depositories=self.kinetics_depositories,
            statmech_libraries = self.statmech_libraries,
            depository=False,  # Don't bother loading the depository information, as we don't use it
        )

        # Turn off reversibility for families with three products if desired
        if not self.trimolecular_product_reversible:
            for family in self.database.kinetics.families.values():
                if len(family.forward_template.products) > 2:
                    family.reversible = False
                    family.reverse_template = None
                    family.reverse_recipe = None
                    family.reverse = None

        # Determine if trimolecular families are present
        for family in self.database.kinetics.families.values():
            if len(family.forward_template.reactants) > 2:
                logging.info('Trimolecular reactions are turned on')
                self.trimolecular = True
                break
        # Only check products if we want to react them
        if not self.trimolecular and self.trimolecular_product_reversible:
            for family in self.database.kinetics.families.values():
                if len(family.forward_template.products) > 2:
                    logging.info('Trimolecular reactions are turned on')
                    self.trimolecular = True
                    break

        # check libraries
        self.check_libraries()

        # set global binding energies variable
        if self.binding_energies:
            self.database.thermo.set_binding_energies(self.binding_energies)

        # set global variable solvent
        if self.solvent:
            global solvent
            solvent = self.solvent

        # add any forbidden structures in the input file to the forbidden structures database
        for forbidden_structure_entry in self.forbidden_structures:
            label = forbidden_structure_entry.label
            if label in self.database.forbidden_structures.entries:
                raise InputError("""
        Forbidden structure {0} label is already in the forbidden structure database.
        Please choose a different label for this structure that is not already in
        {1}""".format(label,os.path.join(self.database_directory,'forbiddenStructures.py')))
            logging.info('Adding {0} to the forbidden structures database...'.format(label))
            self.database.forbidden_structures.entries[label] = forbidden_structure_entry

        if self.kinetics_estimator == 'rate rules':
            if '!training' not in self.kinetics_depositories:
                logging.info('Adding rate rules from training set in kinetics families...')
                # Temporarily remove species constraints for the training reactions
                copy_species_constraints = copy.copy(self.species_constraints)
                self.species_constraints = {}
                for family in self.database.kinetics.families.values():
                    if not family.auto_generated:
                        family.add_rules_from_training(thermo_database=self.database.thermo)

                    # If requested by the user, write a text file for each kinetics family detailing the source of each entry
                    if self.kinetics_datastore:
                        logging.info(
                            'Writing sources of kinetic entries in family {0} to text file'.format(family.label))
                        path = os.path.join(self.output_directory, 'kinetics_database', family.label + '.txt')
                        with open(path, 'w') as f:
                            for template_label, entries in family.rules.entries.items():
                                f.write("Template [{0}] uses the {1} following source(s):\n".format(template_label,
                                                                                                    str(len(entries))))
                                for entry_index, entry in enumerate(entries):
                                    f.write(str(entry_index+1) + ". " + entry.short_desc + "\n" + entry.long_desc + "\n")
                                f.write('\n')
                            f.write('\n')

                self.species_constraints = copy_species_constraints
            else:
                logging.info('Training set explicitly not added to rate rules in kinetics families...')
            logging.info('Filling in rate rules in kinetics families by averaging...')
            for family in self.database.kinetics.families.values():
                if not family.auto_generated:
                    family.fill_rules_by_averaging_up(verbose=self.verbose_comments)

    def initialize(self, **kwargs):
        """
        Initialize an RMG job using the command-line arguments `args` as returned
        by the :mod:`argparse` package.
        """

        # Save initialization time
        self.initialization_time = time.time()

        # Log start timestamp
        logging.info('RMG execution initiated at ' + time.asctime() + '\n')

        # Print out RMG header
        self.log_header()

        # Read input file
        self.load_input(self.input_file)

        if kwargs.get('restart', ''):
            import rmgpy.rmg.input
            rmgpy.rmg.input.restart_from_seed(path=kwargs['restart'])

        # Check input file
        self.check_input()

        # Properly set filter_reactions to initialize flags properly
        if len(self.model_settings_list) > 0:
            self.filter_reactions = self.model_settings_list[0].filter_reactions

        # Make output subdirectories
        util.make_output_subdirectory(self.output_directory, 'pdep')
        util.make_output_subdirectory(self.output_directory, 'solver')
        util.make_output_subdirectory(self.output_directory, 'kinetics_database')

        # Specifies if details of kinetic database entries should be stored according to user
        try:
            self.kinetics_datastore = kwargs['kinetics_datastore']
        except KeyError:
            self.kinetics_datastore = False

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
        self.load_database()

        for spec in self.initial_species:
            self.reaction_model.add_species_to_edge(spec)

        for reaction_system in self.reaction_systems:
            if isinstance(reaction_system, Reactor):
                reaction_system.finish_termination_criteria()

        # Load restart seed mechanism (if specified)
        if self.restart:
            # Copy the restart files to a separate folder so that the job does not overwrite it
            restart_dir = os.path.join(self.output_directory, 'previous_restart')
            core_restart = os.path.join(restart_dir, 'restart')
            edge_restart = os.path.join(restart_dir, 'restart_edge')
            filters_restart = os.path.join(restart_dir, 'filters')
            util.make_output_subdirectory(self.output_directory, 'previous_restart')
            shutil.copytree(self.core_seed_path, core_restart)
            shutil.copytree(self.edge_seed_path, edge_restart)
            os.mkdir(filters_restart)
            shutil.copyfile(self.filters_path, os.path.join(filters_restart, 'filters.h5'))
            shutil.copyfile(self.species_map_path, os.path.join(filters_restart, 'species_map.yml'))

            # Load the seed mechanism to get the core and edge species
            self.database.kinetics.load_libraries(restart_dir, libraries=['restart', 'restart_edge'])
            self.seed_mechanisms.append('restart')
            self.reaction_libraries.append(('restart_edge', False))

        # Set trimolecular reactant flags of reaction systems
        if self.trimolecular:
            for reaction_system in self.reaction_systems:
                reaction_system.trimolecular = True

        # Do all liquid-phase startup things:
        if self.solvent:
            solvent_data = self.database.solvation.get_solvent_data(self.solvent)
            if not self.reaction_model.core.phase_system.in_nose:
                self.reaction_model.core.phase_system.phases["Default"].set_solvent(solvent_data)
                self.reaction_model.edge.phase_system.phases["Default"].set_solvent(solvent_data)

            diffusion_limiter.enable(solvent_data, self.database.solvation)
            logging.info("Setting solvent data for {0}".format(self.solvent))

            if self.liquid_volumetric_mass_transfer_coefficient_power_law:
                vapor_liquid_mass_transfer.enable(solvent_data, self.database.solvation, self.liquid_volumetric_mass_transfer_coefficient_power_law)
                logging.info("Setting vapor liquid mass transfer with {0} as solvent".format(self.solvent))

            # Set solvent viscosity for reaction filtering
            for reaction_system in self.reaction_systems:
                if reaction_system.T:
                    reaction_system.viscosity = solvent_data.get_solvent_viscosity(reaction_system.T.value_si)

        try:
            self.walltime = kwargs['walltime']
        except KeyError:
            pass

        try:
            self.max_iterations = kwargs['max_iterations']
        except KeyError:
            pass

        data = self.walltime.split(':')
        if not len(data) == 4:
            raise ValueError('Invalid format for wall time {0}; should be DD:HH:MM:SS.'.format(self.walltime))
        self.walltime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3]) + 86400 * int(data[-4])

        # Initialize reaction model

        # Seed mechanisms: add species and reactions from seed mechanism
        # DON'T generate any more reactions for the seed species at this time
        for seed_mechanism in self.seed_mechanisms:
            self.reaction_model.add_seed_mechanism_to_core(seed_mechanism, react=False)

        # Reaction libraries: add species and reactions from reaction library to the edge so
        # that RMG can find them if their rates are large enough
        for library, option in self.reaction_libraries:
            self.reaction_model.add_reaction_library_to_edge(library)

        # Also always add in a few bath gases (since RMG-Java does)
        for label, smiles in [('Ar', '[Ar]'), ('He', '[He]'), ('Ne', '[Ne]'), ('N2', 'N#N')]:
            molecule = Molecule().from_smiles(smiles)
            spec, is_new = self.reaction_model.make_new_species(molecule, label=label, reactive=False)
            if is_new:
                self.initial_species.append(spec)

        # Perform species constraints and forbidden species checks on input species
        for spec in self.initial_species:
            if self.database.forbidden_structures.is_molecule_forbidden(spec.molecule[0]):
                if 'allowed' in self.species_constraints and 'input species' in self.species_constraints['allowed']:
                    spec.explicitly_allowed = True
                    logging.warning("Input species {0} is globally forbidden but will be explicitly allowed "
                            " since input species are permitted by the user's species constraints".format(spec.label))
                else:
                    raise ForbiddenStructureException("Input species {0} is globally forbidden. You may explicitly "
                                                      "allow it by adding 'input species' to the `generatedSpeciesConstraints` `allowed` list.".format(spec.label))
            if fails_species_constraints(spec):
                if 'allowed' in self.species_constraints and 'input species' in self.species_constraints['allowed']:
                    self.species_constraints['explicitlyAllowedMolecules'].append(spec.molecule[0])
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
        if 'allowSingletO2' in self.species_constraints and self.species_constraints['allowSingletO2']:
            pass
        else:
            # Here we get a list of all species that from the user input
            all_inputted_species = [spec for spec in self.initial_species]
            # Because no iterations have taken place, the only things in the core are from seed mechanisms
            all_inputted_species.extend(self.reaction_model.core.species)
            # Because no iterations have taken place, the only things in the edge are from reaction libraries
            all_inputted_species.extend(self.reaction_model.edge.species)

            O2Singlet = Molecule().from_smiles('O=O')
            for spec in all_inputted_species:
                if spec.is_isomorphic(O2Singlet):
                    raise ForbiddenStructureException("Species constraints forbids input species {0} RMG expects the "
                                                      "triplet form of oxygen for correct usage in reaction families. "
                                                      "Please change your input to SMILES='[O][O]' If you actually "
                                                      "want to use the singlet state, set the allowSingletO2=True "
                                                      "inside of the Species Constraints block in your input file."
                                                      .format(spec.label))

        for spec in self.initial_species:
            submit(spec, self.solvent)
            if vapor_liquid_mass_transfer.enabled:
                spec.get_liquid_volumetric_mass_transfer_coefficient_data()
                spec.get_henry_law_constant_data()

        # Add nonreactive species (e.g. bath gases) to core first
        # This is necessary so that the PDep algorithm can identify the bath gas
        for spec in self.initial_species:
            if not spec.reactive:
                self.reaction_model.enlarge(spec)
        for spec in self.initial_species:
            if spec.reactive:
                self.reaction_model.enlarge(spec)

        # chatelak: store constant SPC indices in the reactor attributes if any constant SPC provided in the input file
        # advantages to write it here: this is run only once (as species indexes does not change over the generation)
        if self.solvent is not None:
            for index, reaction_system in enumerate(self.reaction_systems):
                if not isinstance(reaction_system, Reactor) and reaction_system.const_spc_names is not None:  # if no constant species provided do nothing
                    reaction_system.get_const_spc_indices(
                        self.reaction_model.core.species)  # call the function to identify indices in the solver

        self.initialize_reaction_threshold_and_react_flags()
        if self.filter_reactions and self.init_react_tuples:
            self.react_init_tuples()
        self.reaction_model.initialize_index_species_dict()

        self.initialize_seed_mech()

    def register_listeners(self):
        """
        Attaches listener classes depending on the options
        found in the RMG input file.
        """

        self.attach(ChemkinWriter(self.output_directory))
        self.attach(RMSWriter(self.output_directory))

        if self.generate_output_html:
            self.attach(OutputHTMLWriter(self.output_directory))

        if self.quantum_mechanics:
            self.attach(QMDatabaseWriter())

        self.attach(ExecutionStatsWriter(self.output_directory))

        if self.save_simulation_profiles:

            for index, reaction_system in enumerate(self.reaction_systems):
                if isinstance(reaction_system, Reactor):
                    typ = type(reaction_system)
                    raise InputError(f"save_simulation_profiles=True not compatible with reactor of type {typ}")
                reaction_system.attach(SimulationProfileWriter(
                    self.output_directory, index, self.reaction_model.core.species))
                reaction_system.attach(SimulationProfilePlotter(
                    self.output_directory, index, self.reaction_model.core.species))

    def execute(self, initialize=True, **kwargs):
        """
        Execute an RMG job using the command-line arguments `args` as returned
        by the :mod:`argparse` package.
        ``initialize`` is a ``bool`` type flag used to determine whether to call self.initialize()
        """

        if initialize:
            self.initialize(**kwargs)

        # register listeners
        self.register_listeners()

        self.done = False

        # determine min and max values for T and P (don't determine P values for liquid reactors)
        self.Tmin = min([x.Trange[0].value_si if x.Trange else x.T.value_si for x in self.reaction_systems])
        self.Tmax = max([x.Trange[1].value_si if x.Trange else x.T.value_si for x in self.reaction_systems])
        try:
            self.Pmin = min([x.Prange[0].value_si if hasattr(x, 'Prange') and x.Prange else x.P.value_si for x in self.reaction_systems])
            self.Pmax = max([x.Prange[1].value_si if hasattr(x, 'Prange') and x.Prange else x.P.value_si for x in self.reaction_systems])
        except AttributeError:
            pass

        self.rmg_memories = []

        logging.info('Initialization complete. Starting model generation.\n')

        # Initiate first reaction discovery step after adding all core species
        for index, reaction_system in enumerate(self.reaction_systems):
            # Initialize memory object to track conditions for ranged reactors
            self.rmg_memories.append(RMG_Memory(reaction_system, self.balance_species))
            self.rmg_memories[index].generate_cond()
            log_conditions(self.rmg_memories, index)

            # Update react flags
            if self.filter_reactions:
                # Run the reaction system to update threshold and react flags
                if isinstance(reaction_system, Reactor):
                    self.update_reaction_threshold_and_react_flags(
                        rxn_sys_unimol_threshold=np.zeros((len(self.reaction_model.core.species),),bool),
                        rxn_sys_bimol_threshold=np.zeros((len(self.reaction_model.core.species),len(self.reaction_model.core.species)),bool),
                        rxn_sys_trimol_threshold=np.zeros((len(self.reaction_model.core.species),len(self.reaction_model.core.species),len(self.reaction_model.core.species)),bool),
                    )

                else:
                    reaction_system.initialize_model(
                        core_species=self.reaction_model.core.species,
                        core_reactions=self.reaction_model.core.reactions,
                        edge_species=[],
                        edge_reactions=[],
                        pdep_networks=self.reaction_model.network_list,
                        atol=self.simulator_settings_list[0].atol,
                        rtol=self.simulator_settings_list[0].rtol,
                        filter_reactions=True,
                        conditions=self.rmg_memories[index].get_cond(),
                    )

                    self.update_reaction_threshold_and_react_flags(
                        rxn_sys_unimol_threshold=reaction_system.unimolecular_threshold,
                        rxn_sys_bimol_threshold=reaction_system.bimolecular_threshold,
                        rxn_sys_trimol_threshold=reaction_system.trimolecular_threshold,
                    )

                logging.info('Generating initial reactions for reaction system {0}...'.format(index + 1))
            else:
                # If we're not filtering reactions, then we only need to react
                # the first reaction system since they share the same core
                if index > 0:
                    continue
                logging.info('Generating initial reactions...')

            # React core species to enlarge edge
            self.reaction_model.enlarge(react_edge=True,
                                        unimolecular_react=self.unimolecular_react,
                                        bimolecular_react=self.bimolecular_react,
                                        trimolecular_react=self.trimolecular_react)

        if not np.isinf(self.model_settings_list[0].thermo_tol_keep_spc_in_edge):
            self.reaction_model.set_thermodynamic_filtering_parameters(
                self.Tmax,
                thermo_tol_keep_spc_in_edge=self.model_settings_list[0].thermo_tol_keep_spc_in_edge,
                min_core_size_for_prune=self.model_settings_list[0].min_core_size_for_prune,
                maximum_edge_species=self.model_settings_list[0].maximum_edge_species,
                reaction_systems=self.reaction_systems
            )

        if not np.isinf(self.model_settings_list[0].thermo_tol_keep_spc_in_edge):
            self.reaction_model.thermo_filter_down(maximum_edge_species=self.model_settings_list[0].maximum_edge_species)

        logging.info('Completed initial enlarge edge step.\n')

        self.save_everything()

        if self.generate_seed_each_iteration:
            self.make_seed_mech()

        max_num_spcs_hit = False  # default

        for q, model_settings in enumerate(self.model_settings_list):
            if len(self.simulator_settings_list) > 1:
                simulator_settings = self.simulator_settings_list[q]
            else:  # if they only provide one input for simulator use that everytime
                simulator_settings = self.simulator_settings_list[0]

            self.filter_reactions = model_settings.filter_reactions

            logging.info('Beginning model generation stage {0}...\n'.format(q + 1))

            self.done = False

            # Main RMG loop
            while not self.done:
                # iteration number starts at 0. Increment it before entering make_seed_mech
                self.reaction_model.iteration_num += 1
                self.done = True

                if self.generate_seed_each_iteration:
                    self.make_seed_mech()

                all_terminated = True
                num_core_species = len(self.reaction_model.core.species)

                prunable_species = self.reaction_model.edge.species[:]
                prunable_networks = self.reaction_model.network_list[:]

                for index, reaction_system in enumerate(self.reaction_systems):

                    reaction_system.prunable_species = prunable_species  # these lines reset pruning for a new cycle
                    reaction_system.prunable_networks = prunable_networks
                    reaction_system.reset_max_edge_species_rate_ratios()

                    for p in range(reaction_system.n_sims):
                        reactor_done = True
                        objects_to_enlarge = []

                        conditions = self.rmg_memories[index].get_cond()
                        if conditions and self.solvent:
                            T = conditions['T']
                            # Set solvent viscosity
                            solvent_data = self.database.solvation.get_solvent_data(self.solvent)
                            reaction_system.viscosity = solvent_data.get_solvent_viscosity(T)

                        self.reaction_system = reaction_system
                        # Conduct simulation
                        logging.info('Conducting simulation of reaction system %s...' % (index + 1))
                        prune = True

                        self.reaction_model.adjust_surface()

                        if num_core_species < model_settings.min_core_size_for_prune:
                            # Turn pruning off if we haven't reached minimum core size.
                            prune = False

                        try:
                            if isinstance(reaction_system, Reactor):
                                terminated,resurrected,obj,unimolecular_threshold,bimolecular_threshold,trimolecular_threshold,max_edge_species_rate_ratios,t,x = reaction_system.simulate(model_settings=model_settings,
                                    simulator_settings=simulator_settings,
                                    conditions=self.rmg_memories[index].get_cond()
                                )
                                reaction_system.unimolecular_threshold = unimolecular_threshold
                                reaction_system.bimolecular_threshold = bimolecular_threshold
                                reaction_system.trimolecular_threshold = trimolecular_threshold
                                if hasattr(reaction_system,"max_edge_species_rate_ratios"):
                                    max_edge_species_rate_ratios_temp = np.zeros(len(max_edge_species_rate_ratios))
                                    for i in range(len(max_edge_species_rate_ratios)):
                                        if i < len(reaction_system.max_edge_species_rate_ratios):
                                            max_edge_species_rate_ratios_temp[i] = max(reaction_system.max_edge_species_rate_ratios[i],max_edge_species_rate_ratios[i])
                                        else:
                                            max_edge_species_rate_ratios_temp[i] = max_edge_species_rate_ratios[i]
                                    reaction_system.max_edge_species_rate_ratios = max_edge_species_rate_ratios_temp
                                else:
                                    reaction_system.max_edge_species_rate_ratios = max_edge_species_rate_ratios
                                new_surface_species = []
                                new_surface_reactions =  []
                                obj_temp = []
                                for item in obj:
                                    if hasattr(item,"name"):
                                        obj_temp.append(self.reaction_model.edge.phase_system.species_dict[item.name])
                                    else: #Reaction
                                        for val in item.reactants+item.products:
                                            spc = self.reaction_model.edge.phase_system.species_dict[val.name]
                                            if spc not in self.reaction_model.core.species:
                                                obj_temp.append(spc)
                                        assert len(obj_temp) > 0
                                obj = obj_temp
                            else:
                                terminated, resurrected, obj, new_surface_species, new_surface_reactions, t, x = reaction_system.simulate(
                                    core_species=self.reaction_model.core.species,
                                    core_reactions=self.reaction_model.core.reactions,
                                    edge_species=self.reaction_model.edge.species,
                                    edge_reactions=self.reaction_model.edge.reactions,
                                    surface_species=self.reaction_model.surface.species,
                                    surface_reactions=self.reaction_model.surface.reactions,
                                    pdep_networks=self.reaction_model.network_list,
                                    prune=prune,
                                    model_settings=model_settings,
                                    simulator_settings=simulator_settings,
                                    conditions=self.rmg_memories[index].get_cond()
                                )
                        except:
                            logging.error("Model core reactions:")
                            if len(self.reaction_model.core.reactions) > 5:
                                logging.error("Too many to print in detail")
                            else:
                                from arkane.output import prettify
                                logging.error(prettify(repr(self.reaction_model.core.reactions)))
                            if not self.generate_seed_each_iteration:  # Then we haven't saved the seed mechanism yet
                                self.make_seed_mech()  # Just in case the user wants to restart from this
                            raise

                        self.rmg_memories[index].add_t_conv_N(t, x, len(obj))
                        self.rmg_memories[index].generate_cond()
                        log_conditions(self.rmg_memories, index)

                        reactor_done = self.reaction_model.add_new_surface_objects(obj, new_surface_species,
                                                                                   new_surface_reactions, reaction_system)

                        all_terminated = all_terminated and terminated
                        logging.info('')

                        # If simulation is invalid, note which species should be added to
                        # the core
                        if obj != [] and not (obj is None):
                            objects_to_enlarge = self.process_to_species_networks(obj)

                            reactor_done = False
                        # Enlarge objects identified by the simulation for enlarging
                        # These should be Species or Network objects
                        logging.info('')

                        objects_to_enlarge = list(set(objects_to_enlarge))

                        # Add objects to enlarge to the core first
                        for objectToEnlarge in objects_to_enlarge:
                            self.reaction_model.enlarge(objectToEnlarge)

                        if model_settings.filter_reactions:
                            # Run a raw simulation to get updated reaction system threshold values
                            # Run with the same conditions as with pruning off
                            temp_model_settings = deepcopy(model_settings)
                            temp_model_settings.tol_keep_in_edge = 0
                            if not resurrected:
                                try:
                                    if isinstance(reaction_system, Reactor):
                                        terminated,resurrected,obj,unimolecular_threshold,bimolecular_threshold,trimolecular_threshold,max_edge_species_rate_ratios,t,x = reaction_system.simulate(model_settings=model_settings,
                                            simulator_settings=simulator_settings,
                                            conditions=self.rmg_memories[index].get_cond()
                                        )
                                        reaction_system.unimolecular_threshold = unimolecular_threshold
                                        reaction_system.bimolecular_threshold = bimolecular_threshold
                                        reaction_system.trimolecular_threshold = trimolecular_threshold
                                        if hasattr(reaction_system,"max_edge_species_rate_ratios"):
                                            max_edge_species_rate_ratios_temp = np.zeros(len(max_edge_species_rate_ratios))
                                            for i in range(len(max_edge_species_rate_ratios)):
                                                if i < len(reaction_system.max_edge_species_rate_ratios):
                                                    max_edge_species_rate_ratios_temp[i] = max(reaction_system.max_edge_species_rate_ratios[i],max_edge_species_rate_ratios[i])
                                                else:
                                                    max_edge_species_rate_ratios_temp[i] = max_edge_species_rate_ratios[i]
                                            reaction_system.max_edge_species_rate_ratios = max_edge_species_rate_ratios_temp
                                        else:
                                            reaction_system.max_edge_species_rate_ratios = max_edge_species_rate_ratios
                                    else:
                                        reaction_system.simulate(
                                            core_species=self.reaction_model.core.species,
                                            core_reactions=self.reaction_model.core.reactions,
                                            edge_species=[],
                                            edge_reactions=[],
                                            surface_species=self.reaction_model.surface.species,
                                            surface_reactions=self.reaction_model.surface.reactions,
                                            pdep_networks=self.reaction_model.network_list,
                                            model_settings=temp_model_settings,
                                            simulator_settings=simulator_settings,
                                            conditions=self.rmg_memories[index].get_cond()
                                        )
                                except:
                                    self.update_reaction_threshold_and_react_flags(
                                        rxn_sys_unimol_threshold=reaction_system.unimolecular_threshold,
                                        rxn_sys_bimol_threshold=reaction_system.bimolecular_threshold,
                                        rxn_sys_trimol_threshold=reaction_system.trimolecular_threshold,
                                        skip_update=True)
                                    logging.warning('Reaction thresholds/flags for Reaction System {0} was not updated '
                                                    'due to simulation failure'.format(index + 1))
                                else:
                                    self.update_reaction_threshold_and_react_flags(
                                        rxn_sys_unimol_threshold=reaction_system.unimolecular_threshold,
                                        rxn_sys_bimol_threshold=reaction_system.bimolecular_threshold,
                                        rxn_sys_trimol_threshold=reaction_system.trimolecular_threshold
                                    )
                            else:
                                self.update_reaction_threshold_and_react_flags(
                                    rxn_sys_unimol_threshold=reaction_system.unimolecular_threshold,
                                    rxn_sys_bimol_threshold=reaction_system.bimolecular_threshold,
                                    rxn_sys_trimol_threshold=reaction_system.trimolecular_threshold,
                                    skip_update=True
                                )
                                logging.warning('Reaction thresholds/flags for Reaction System {0} was not updated due '
                                                'to resurrection'.format(index + 1))

                            logging.info('')
                        else:
                            self.update_reaction_threshold_and_react_flags()

                        if not np.isinf(model_settings.thermo_tol_keep_spc_in_edge):
                            self.reaction_model.set_thermodynamic_filtering_parameters(
                                self.Tmax,
                                thermo_tol_keep_spc_in_edge=model_settings.thermo_tol_keep_spc_in_edge,
                                min_core_size_for_prune=model_settings.min_core_size_for_prune,
                                maximum_edge_species=model_settings.maximum_edge_species,
                                reaction_systems=self.reaction_systems
                            )

                        old_edge_size = len(self.reaction_model.edge.reactions)
                        old_core_size = len(self.reaction_model.core.reactions)
                        self.reaction_model.enlarge(react_edge=True,
                                                    unimolecular_react=self.unimolecular_react,
                                                    bimolecular_react=self.bimolecular_react,
                                                    trimolecular_react=self.trimolecular_react)

                        if old_edge_size != len(self.reaction_model.edge.reactions) or old_core_size != len(
                                self.reaction_model.core.reactions):
                            reactor_done = False

                        if not np.isinf(self.model_settings_list[0].thermo_tol_keep_spc_in_edge):
                            self.reaction_model.thermo_filter_down(maximum_edge_species=model_settings.maximum_edge_species)

                        max_num_spcs_hit = len(self.reaction_model.core.species) >= model_settings.max_num_species

                        self.save_everything()

                        if max_num_spcs_hit:  # breaks the n_sims loop
                            # self.done is still True, which will break the while loop
                            break

                        if not reactor_done:
                            self.done = False

                    if max_num_spcs_hit:  # breaks the reaction_systems loop
                        break

                if not self.done:  # There is something that needs exploring/enlarging

                    # If we reached our termination conditions, then try to prune
                    # species from the edge
                    if all_terminated and model_settings.tol_keep_in_edge > 0.0:
                        logging.info('Attempting to prune...')
                        self.reaction_model.prune(self.reaction_systems, model_settings.tol_keep_in_edge,
                                                  model_settings.tol_move_to_core,
                                                  model_settings.maximum_edge_species,
                                                  model_settings.min_species_exist_iterations_for_prune)
                        # Perform garbage collection after pruning
                        collected = gc.collect()
                        logging.info('Garbage collector: collected %d objects.' % collected)

                # Consider stopping gracefully if the next iteration might take us
                # past the wall time
                if self.walltime > 0 and len(self.exec_time) > 1:
                    t = self.exec_time[-1]
                    dt = self.exec_time[-1] - self.exec_time[-2]
                    if t + 3 * dt > self.walltime:
                        logging.info('MODEL GENERATION TERMINATED')
                        logging.info('')
                        logging.info(
                            'There is not enough time to complete the next iteration before the wall time is reached.')
                        logging.info('The output model may be incomplete.')
                        logging.info('')
                        core_spec, core_reac, edge_spec, edge_reac = self.reaction_model.get_model_size()
                        logging.info('The current model core has %s species and %s reactions' % (core_spec, core_reac))
                        logging.info('The current model edge has %s species and %s reactions' % (edge_spec, edge_reac))
                        return

                if self.max_iterations and (self.reaction_model.iteration_num >= self.max_iterations):
                    logging.info('MODEL GENERATION TERMINATED')
                    logging.info('')
                    logging.info('The maximum number of iterations of {0} has been reached'.format(self.max_iterations))
                    logging.info('The output model may be incomplete.')
                    logging.info('')
                    core_spec, core_reac, edge_spec, edge_reac = self.reaction_model.get_model_size()
                    logging.info('The current model core has %s species and %s reactions' % (core_spec, core_reac))
                    logging.info('The current model edge has %s species and %s reactions' % (edge_spec, edge_reac))
                    return

            if max_num_spcs_hit:  # resets maxNumSpcsHit and continues the settings for loop
                logging.info('The maximum number of species ({0}) has been hit, Exiting stage {1} ...'.format(
                    model_settings.max_num_species, q + 1))
                max_num_spcs_hit = False

        # Save the final seed mechanism
        self.make_seed_mech()

        self.run_model_analysis()

        # generate Cantera files chem.yaml & chem_annotated.yaml in a designated `cantera` output folder
        try:
            if any([s.contains_surface_site() for s in self.reaction_model.core.species]):
                self.generate_cantera_files(os.path.join(self.output_directory, 'chemkin', 'chem-gas.inp'),
                                            surface_file=(
                                              os.path.join(self.output_directory, 'chemkin', 'chem-surface.inp')))
                self.generate_cantera_files(os.path.join(self.output_directory, 'chemkin', 'chem_annotated-gas.inp'),
                                            surface_file=(os.path.join(self.output_directory, 'chemkin',
                                                                    'chem_annotated-surface.inp')))
            else:  # gas phase only
                self.generate_cantera_files(os.path.join(self.output_directory, 'chemkin', 'chem.inp'))
                self.generate_cantera_files(os.path.join(self.output_directory, 'chemkin', 'chem_annotated.inp'))
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
        core_spec, core_reac, edge_spec, edge_reac = self.reaction_model.get_model_size()
        logging.info('The final model core has %s species and %s reactions' % (core_spec, core_reac))
        logging.info('The final model edge has %s species and %s reactions' % (edge_spec, edge_reac))

        self.finish()

    def run_model_analysis(self, number=10):
        """
        Run sensitivity and uncertainty analysis if requested.
        """
        # Run sensitivity analysis post-model generation if sensitivity analysis is on
        for index, reaction_system in enumerate(self.reaction_systems):

            if reaction_system.sensitive_species and reaction_system.sens_conditions:
                logging.info('Conducting sensitivity analysis of reaction system {0}...'.format(index + 1))

                if reaction_system.sensitive_species == ['all']:
                    reaction_system.sensitive_species = self.reaction_model.core.species

                sens_worksheet = []
                for spec in reaction_system.sensitive_species:
                    csvfile_path = os.path.join(self.output_directory, 'solver',
                                                'sensitivity_{0}_SPC_{1}.csv'.format(index + 1, spec.index))
                    sens_worksheet.append(csvfile_path)

                terminated, resurrected, obj, surface_species, surface_reactions, t, x = reaction_system.simulate(
                    core_species=self.reaction_model.core.species,
                    core_reactions=self.reaction_model.core.reactions,
                    edge_species=self.reaction_model.edge.species,
                    edge_reactions=self.reaction_model.edge.reactions,
                    surface_species=[],
                    surface_reactions=[],
                    pdep_networks=self.reaction_model.network_list,
                    sensitivity=True,
                    sens_worksheet=sens_worksheet,
                    model_settings=ModelSettings(tol_move_to_core=1e8, tol_interrupt_simulation=1e8),
                    simulator_settings=self.simulator_settings_list[-1],
                    conditions=reaction_system.sens_conditions,
                )

                plot_sensitivity(self.output_directory, index, reaction_system.sensitive_species, number=number)

        self.run_uncertainty_analysis()

    def run_uncertainty_analysis(self):
        """
        Run uncertainty analysis if proper settings are available.
        """
        if self.uncertainty is not None and self.uncertainty['global']:
            try:
                import muq
            except ImportError:
                logging.error('Unable to import MUQ. Skipping global uncertainty analysis.')
                self.uncertainty['global'] = False
            else:
                import re
                import random
                from rmgpy.tools.canteramodel import Cantera
                from rmgpy.tools.globaluncertainty import ReactorPCEFactory

        if self.uncertainty is not None and self.uncertainty['local']:
            correlation = []
            if self.uncertainty['uncorrelated']: correlation.append(False)
            if self.uncertainty['correlated']: correlation.append(True)

            # Set up Uncertainty object
            uncertainty = Uncertainty(output_directory=self.output_directory)
            uncertainty.database = self.database
            uncertainty.species_list, uncertainty.reaction_list = self.reaction_model.core.species, self.reaction_model.core.reactions
            uncertainty.extract_sources_from_model()

            # Reload reaction families with verbose comments if necessary
            if not self.verbose_comments:
                logging.info('Reloading kinetics families with verbose comments for uncertainty analysis...')
                self.database.kinetics.load_families(os.path.join(self.database_directory, 'kinetics', 'families'),
                                                     self.kinetics_families, self.kinetics_depositories)
                # Temporarily remove species constraints for the training reactions
                self.species_constraints, speciesConstraintsCopy = {}, self.species_constraints
                for family in self.database.kinetics.families.values():
                    if not family.auto_generated:
                        family.add_rules_from_training(thermo_database=self.database.thermo)
                        family.fill_rules_by_averaging_up(verbose=True)
                self.species_constraints = speciesConstraintsCopy

            for correlated in correlation:
                uncertainty.assign_parameter_uncertainties(correlated=correlated)

                for index, reaction_system in enumerate(self.reaction_systems):
                    if reaction_system.sensitive_species and reaction_system.sens_conditions:
                        logging.info('Conducting {0}correlated local uncertainty analysis for '
                                     'reaction system {1}...\n'.format('un' if not correlated else '', index + 1))
                        results = uncertainty.local_analysis(reaction_system.sensitive_species,
                                                            reaction_system_index=index,
                                                            correlated=correlated,
                                                            number=self.uncertainty['localnum'])
                        logging.info('Local uncertainty analysis results for reaction system {0}:\n'.format(index + 1))
                        local_result, local_result_str = process_local_results(results, reaction_system.sensitive_species,
                                                                               number=self.uncertainty['localnum'])
                        logging.info(local_result_str)

                        if self.uncertainty['global']:
                            logging.info('Conducting {0}correlated global uncertainty analysis for '
                                         'reaction system {1}...'.format('un' if not correlated else '', index + 1))
                            # Get simulation conditions
                            for criteria in reaction_system.termination:
                                if isinstance(criteria, TerminationTime):
                                    time_criteria = ([criteria.time.value], criteria.time.units)
                                    break
                            else:
                                if self.uncertainty['time']:
                                    time_criteria = ([self.uncertainty['time'].value], self.uncertainty['time'].units)
                                else:
                                    raise InputError('If terminationTime was not specified in the reactor options block, it'
                                                    'must be specified in the uncertainty options block for global uncertainty'
                                                    'analysis.')


                            Tlist = ([reaction_system.sens_conditions['T']], 'K')
                            try:
                                Plist = ([reaction_system.sens_conditions['P']], 'Pa')
                            except KeyError:
                                #LiquidReactor
                                Plist = ([1e8], 'Pa')
                                logging.warning('Using 1e8 Pa as the reaction system pressure to approximate liquid phase density.')

                            mol_frac_list = [reaction_system.sens_conditions.copy()]
                            del mol_frac_list[0]['T']
                            if 'P' in mol_frac_list[0]:
                                del mol_frac_list[0]['P']

                            # Set up Cantera reactor
                            job = Cantera(species_list=uncertainty.species_list, reaction_list=uncertainty.reaction_list,
                                          output_directory=os.path.join(self.output_directory, 'global_uncertainty'))
                            job.load_model()
                            job.generate_conditions(
                                reactor_type_list=['IdealGasConstPressureTemperatureReactor'],
                                reaction_time_list=time_criteria,
                                mol_frac_list=mol_frac_list,
                                Tlist=Tlist,
                                Plist=Plist,
                            )

                            # Extract uncertain parameters from local analysis
                            k_params = []
                            g_params = []
                            for spc in reaction_system.sensitive_species:
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
                                output_species_list=reaction_system.sensitive_species,
                                k_params=k_params,
                                k_uncertainty=uncertainty.kinetic_input_uncertainties,
                                g_params=g_params,
                                g_uncertainty=uncertainty.thermo_input_uncertainties,
                                correlated=correlated,
                                logx=self.uncertainty['logx'],
                            )

                            logging.info('Generating PCEs...')
                            reactor_pce_factory.generate_pce(
                                max_adapt_time=self.uncertainty['pcetime'],
                                error_tol=self.uncertainty['pcetol'],
                                max_evals=self.uncertainty['pceevals'],
                            )

                            # Try a test point to see how well the PCE performs
                            reactor_pce_factory.compare_output(
                                [random.uniform(-1.0, 1.0) for i in range(len(k_params) + len(g_params))])

                            # Analyze results and save statistics
                            reactor_pce_factory.analyze_results()
                    else:
                        logging.info('Unable to run uncertainty analysis. Must specify sensitivity analysis options in '
                                     'reactor options.')

    def check_model(self):
        """
        Run checks on the RMG model
        """
        logging.info('Performing final model checks...')

        # Check that no two species in core or edge are isomorphic
        for i, spc in enumerate(self.reaction_model.core.species):
            for j in range(i):
                spc2 = self.reaction_model.core.species[j]
                if spc.is_isomorphic(spc2):
                    raise CoreError(
                        'Although the model has completed, species {0} is isomorphic to species {1} in the core. '
                        'Please open an issue on GitHub with the following output:'
                        '\n{2}\n{3}'.format(spc.label, spc2.label, spc.to_adjacency_list(), spc2.to_adjacency_list())
                    )

        for i, spc in enumerate(self.reaction_model.edge.species):
            for j in range(i):
                spc2 = self.reaction_model.edge.species[j]
                if spc.is_isomorphic(spc2):
                    logging.warning(
                        'Species {0} is isomorphic to species {1} in the edge. This does not affect '
                        'the generated model. If you would like to report this to help make RMG better '
                        'please open a GitHub issue with the following output:'
                        '\n{2}\n{3}'.format(spc.label, spc2.label, spc.to_adjacency_list(), spc2.to_adjacency_list())
                    )

        # Check all core reactions (in both directions) for collision limit violation
        violators = []
        for rxn in self.reaction_model.core.reactions:
            if rxn.is_surface_reaction():
                # Don't check collision limits for surface reactions.
                continue
            violator_list = rxn.check_collision_limit_violation(t_min=self.Tmin, t_max=self.Tmax,
                                                                p_min=self.Pmin, p_max=self.Pmax)
            if violator_list:
                violators.extend(violator_list)
        # Whether or not violators were found, rename 'collision_rate_violators.log' if it exists
        new_file = os.path.join(self.output_directory, 'collision_rate_violators.log')
        old_file = os.path.join(self.output_directory, 'collision_rate_violators_OLD.log')
        if os.path.isfile(new_file):
            # If there are no violators, yet the violators log exists (probably from a previous run
            # in the same folder), rename it.
            if os.path.isfile(old_file):
                os.remove(old_file)
            os.rename(new_file, old_file)
        if violators:
            logging.info("\n")
            logging.warning("{0} CORE reactions violate the collision rate limit!"
                            "\nSee the 'collision_rate_violators.log' for details.\n\n".format(len(violators)))
            with open(new_file, 'w') as violators_f:
                violators_f.write('*** Collision rate limit violators report ***\n'
                                  '"Violation factor" is the ratio of the rate coefficient to the collision limit'
                                  ' rate at the relevant conditions\n\n')
                for violator in violators:
                    if isinstance(violator[0].kinetics, ThirdBody):
                        rxn_string = violator[0].to_chemkin(self.reaction_model.core.species)
                    else:
                        rxn_string = violator[0].to_chemkin()
                    direction = violator[1]
                    ratio = violator[2]
                    condition = violator[3]
                    violators_f.write(f'{rxn_string}\n'
                                      f'Direction: {direction}\n'
                                      f'Violation factor: {ratio:.2f}\n'
                                      f'Violation condition: {condition}\n\n\n')
        else:
            logging.info("No collision rate violators found in the model's core.")

    def initialize_seed_mech(self):
        """
        Initialize the process of saving the seed mechanism by performing the following:

        1. Create the initial seed mechanism folder (the seed from a previous iterations will be deleted)
        2. Save the restart-from-seed file (unless the current job is itself a restart job)
        3. Ensure that we don't overwrite existing libraries in the database that have the same name as this job
        4. Create the previous_seeds directory to save intermediate seeds if the user gives a value for saveSeedModulus
        """
        # Make the initial seed mechanism folder
        seed_dir = os.path.join(self.output_directory, 'seed')
        if os.path.exists(seed_dir):  # This is a seed from a previous RMG run. Delete it
            shutil.rmtree(seed_dir)
        os.mkdir(seed_dir)

        # Generate a file for restarting from a seed mechanism if this is not a restart job
        if not self.restart:
            with open(os.path.join(self.output_directory, 'restart_from_seed.py'), 'w') as f:
                f.write('restartFromSeed(path=\'seed\')\n\n')
                with open(self.input_file, 'r') as input_file:
                    f.write(input_file.read())

        # Change self.name to be unique within the thermo/kinetics libraries by adding integers to the end of the
        # name to prevent overwriting
        name = self.name
        if self.save_seed_to_database:  # make sure we don't overwrite current libraries
            thermo_names = list(self.database.thermo.libraries.keys())
            kinetics_names = list(self.database.kinetics.libraries.keys())

            if name in thermo_names or name in kinetics_names:
                q = 1
                while name + str(q) in thermo_names or name + str(q) in kinetics_names:
                    q += 1
                self.name = name + str(q)

        previous_seeds_dir = os.path.join(self.output_directory, 'previous_seeds')
        if os.path.exists(previous_seeds_dir):  # These are seeds from a previous RMG run. Delete them
            shutil.rmtree(previous_seeds_dir)
        if self.save_seed_modulus != -1:
            os.makedirs(previous_seeds_dir, exist_ok=True)

    def make_seed_mech(self):
        """
        Save a seed mechanism (both core and edge) in the 'seed' sub-folder of the output directory. Additionally, save
        the filter tensors to the 'seed/filters' sub-folder so that the RMG job can be restarted from a seed mechanism.
        If `self.save_seed_to_database` is True then the seed mechanism is also saved as libraries (one each for the
        core and edge) in the RMG-database.

        Notes:
            `initialize_seed_mech` should be called one time before this function is ever called.

        """
        logging.info('Making seed mechanism...')

        name = self.name

        seed_dir = os.path.join(self.output_directory, 'seed')
        filter_dir = os.path.join(seed_dir, 'filters')
        previous_seeds_dir = os.path.join(self.output_directory, 'previous_seeds')
        temp_seed_dir = os.path.join(self.output_directory, 'seed_tmp')

        # Move the seed from the previous iteration to a temporary directory in case we run into errors

        # First, remove the temporary seed mechanism if it exists from a previous run
        if os.path.exists(temp_seed_dir):
            shutil.rmtree(temp_seed_dir)
        try:
            os.rename(seed_dir, temp_seed_dir)
        except PermissionError:  # The Windows Subsystem for Linux (WSL) can have problems with renaming
            # Try copying over the files instead. Unfortunately, this takes more time
            if os.path.exists(temp_seed_dir):  # First, delete the contents of the old folder if it exists
                shutil.rmtree(temp_seed_dir)
            shutil.copytree(seed_dir, temp_seed_dir)

            # Now remove the contents of the seed directory
            shutil.rmtree(seed_dir)

        # Now that we have moved the seed mechanism folder, create a new one
        os.mkdir(seed_dir)

        try:
            species_list = self.reaction_model.core.species
            reaction_list = self.reaction_model.core.reactions
            edge_species_list = self.reaction_model.edge.species
            edge_reaction_list = self.reaction_model.edge.reactions

            # Make species labels independent
            old_labels = self.make_species_labels_independent(species_list)
            edge_old_labels = self.make_species_labels_independent(edge_species_list)

            # load kinetics library entries
            kinetics_library = KineticsLibrary(name=name, auto_generated=True)
            kinetics_library.entries = {}
            for i in range(len(reaction_list)):
                reaction = reaction_list[i]
                entry = Entry(
                    index=i + 1,
                    label=reaction.to_labeled_str(),
                    item=reaction,
                    data=reaction.kinetics,
                )

                if 'rate rule' in reaction.kinetics.comment:
                    entry.long_desc = reaction.kinetics.comment
                elif hasattr(reaction, 'library') and reaction.library:
                    entry.long_desc = 'Originally from reaction library: ' + reaction.library + "\n" + reaction.kinetics.comment
                else:
                    entry.long_desc = reaction.kinetics.comment

                kinetics_library.entries[i + 1] = entry

            # load kinetics library entries
            edge_kinetics_library = KineticsLibrary(name=name + '_edge', auto_generated=True)
            edge_kinetics_library.entries = {}
            for i, reaction in enumerate(edge_reaction_list):
                entry = Entry(
                    index=i + 1,
                    label=reaction.to_labeled_str(),
                    item=reaction,
                    data=reaction.kinetics,
                )
                try:
                    entry.long_desc = 'Originally from reaction library: ' + reaction.library + "\n" + reaction.kinetics.comment
                except AttributeError:
                    entry.long_desc = reaction.kinetics.comment
                edge_kinetics_library.entries[i + 1] = entry

            # save in database
            if self.save_seed_to_database:
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
                if self.unimolecular_threshold is not None:
                    f.create_dataset('unimolecular_threshold', data=self.unimolecular_threshold)
                if self.bimolecular_threshold is not None:
                    f.create_dataset('bimolecular_threshold', data=self.bimolecular_threshold)
                if self.trimolecular_threshold is not None:
                    f.create_dataset('trimolecular_threshold', data=self.trimolecular_threshold)

            # Save a map of species indices
            spcs_map = [spc.molecule[0].to_adjacency_list() for spc in self.reaction_model.core.species]

            with open(os.path.join(filter_dir, 'species_map.yml'), 'w') as f:
                yaml.dump(data=spcs_map, stream=f)

            # Also, save the seed to the previous_seeds directory on specified iterations
            if self.save_seed_modulus != -1 and self.reaction_model.iteration_num % self.save_seed_modulus == 0:
                dst = os.path.join(previous_seeds_dir, 'iteration_number_{0}'.format(self.reaction_model.iteration_num))
                logging.info('Copying seed from seed directory to {0}'.format(dst))
                shutil.copytree(seed_dir, dst)

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

    def make_species_labels_independent(self, species):
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
    def process_to_species_networks(self, obj):
        """
        breaks down the objects returned by simulate into Species and PDepNetwork
        components
        """

        if isinstance(obj, PDepNetwork):
            out = [self.process_pdep_networks(obj)]
            return out
        elif isinstance(obj, Species):
            return [obj]
        elif isinstance(obj, Reaction):
            return list(self.process_reactions_to_species(obj))
        elif isinstance(obj, list):  # list of species
            rspcs = self.process_reactions_to_species([k for k in obj if isinstance(k, Reaction)])
            spcs = {k for k in obj if isinstance(k, Species)} | rspcs
            nworks, pspcs = self.process_pdep_networks([k for k in obj if isinstance(k, PDepNetwork)])
            spcs = list(spcs - pspcs)  # avoid duplicate species
            return spcs + nworks
        else:
            raise TypeError("improper call, obj input was incorrect")

    def process_pdep_networks(self, obj):
        """
        properly processes PDepNetwork objects and lists of PDepNetwork objects returned from simulate
        """
        reaction_system = self.reaction_system
        if isinstance(obj, PDepNetwork):
            # Determine which species in that network has the highest leak rate
            # We do this here because we need a temperature and pressure
            # Store the maximum leak species along with the associated network
            ob = (obj, obj.get_maximum_leak_species(reaction_system.T.value_si, reaction_system.P.value_si))
            return ob
        elif isinstance(obj, list):
            spcs = [ob.get_maximum_leak_species(reaction_system.T.value_si, reaction_system.P.value_si) for ob in obj]
            nworks = [(obj[i], spcs[i]) for i in range(len(obj))]
            return nworks, set(spcs)
        else:
            raise TypeError("improper call, obj input was incorrect")

    def process_reactions_to_species(self, obj):
        """
        properly processes Reaction objects and lists of Reaction objects returned from simulate
        """
        core_species = self.reaction_model.core.species
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

    def generate_cantera_files(self, chemkin_file, **kwargs):
        """
        Convert a chemkin mechanism chem.inp file to a cantera mechanism file chem.yaml
        and save it in the cantera directory
        """
        transport_file = os.path.join(os.path.dirname(chemkin_file), 'tran.dat')
        file_name = os.path.splitext(os.path.basename(chemkin_file))[0] + '.yaml'
        out_name = os.path.join(self.output_directory, 'cantera', file_name)
        if 'surface_file' in kwargs:
            out_name = out_name.replace('-gas.', '.')
        cantera_dir = os.path.dirname(out_name)
        try:
            os.makedirs(cantera_dir)
        except OSError:
            if not os.path.isdir(cantera_dir):
                raise
        if os.path.exists(out_name):
            os.remove(out_name)
        parser = ck2yaml.Parser()
        try:
            parser.convert_mech(chemkin_file, transport_file=transport_file, out_name=out_name, quiet=True, permissive=True,
                               **kwargs)
        except ck2yaml.InputError:
            logging.exception("Error converting to Cantera format.")
            logging.info("Trying again without transport data file.")
            parser.convert_mech(chemkin_file, out_name=out_name, quiet=True, permissive=True, **kwargs)

    def initialize_reaction_threshold_and_react_flags(self):
        num_core_species = len(self.reaction_model.core.species)

        # Initialize everything to react by default, but we will handle the restart and filtering case immediately after
        self.unimolecular_react = np.ones(num_core_species, bool)
        self.bimolecular_react = np.ones((num_core_species, num_core_species), bool)
        if self.trimolecular:
            self.trimolecular_react = np.ones((num_core_species, num_core_species, num_core_species), bool)

        if self.filter_reactions or self.restart:  # Otherwise no need to initialize thresholds or fix react flags
            self.unimolecular_threshold = np.zeros(num_core_species, bool)
            self.bimolecular_threshold = np.zeros((num_core_species, num_core_species), bool)
            if self.trimolecular:
                self.trimolecular_threshold = np.zeros((num_core_species, num_core_species, num_core_species), bool)

            if self.restart:
                # Load in the restart mapping
                with open(os.path.join(self.species_map_path), 'r') as f:
                    restart_species_list = yaml.safe_load(stream=f)

                num_restart_spcs = len(restart_species_list)
                restart_species_list = [Species().from_adjacency_list(adj_list) for adj_list in restart_species_list]

                # Load in the restart filter tensors
                with h5py.File(self.filters_path, 'r') as f:
                    if 'unimolecular_threshold' in f.keys():
                        unimolecular_threshold_restart = f.get('unimolecular_threshold')[()]
                        bimolecular_threshold_restart = f.get('bimolecular_threshold')[()]
                        if self.trimolecular:
                            trimolecular_threshold_restart = f.get('trimolecular_threshold')[()]

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

                    else:  # If we can't find the filters then this is because filtering was not used
                        logging.warning('No filters were found in file {0}. This is to be expected if the restart '
                                        'files specified are from an RMG job without reaction filtering. Therefore, '
                                        'RMG will assume that all of the core species from the restart core seed have '
                                        'already been reacted and will not react them again. Additional species added '
                                        'to the input file but not in the restart core seed WILL be reacted with the '
                                        'rest of the core.'.format(self.filters_path))

                        filters_found = False

                # Reorder the core species to match the indices of the restart filter tensors
                reordered_core_species = []
                for spc in restart_species_list:
                    for j, oldCoreSpc in enumerate(self.reaction_model.core.species):
                        if oldCoreSpc.is_isomorphic(spc, strict=False):
                            reordered_core_species.append(self.reaction_model.core.species.pop(j))
                            break
                    else:
                        raise RuntimeError('Species {0} was defined in the restart file, but was not included in the'
                                           'core.'.format(spc))

                # Append the remaining species left in the core to the very end
                self.reaction_model.core.species = reordered_core_species + self.reaction_model.core.species

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

                    if self.filter_reactions:  # Filling in the filter thresholds will suffice
                        # Fill in the newly initialized filter tensors
                        self.unimolecular_threshold = unimolecular_threshold
                        self.bimolecular_threshold = bimolecular_threshold
                        if self.trimolecular:
                            self.trimolecular_threshold = trimolecular_threshold

                    else:  # We must set the react flags instead. If it was `True` in the threshold, it should not react
                        self.unimolecular_react = np.logical_not(unimolecular_threshold)
                        self.bimolecular_react = np.logical_not(bimolecular_threshold)
                        if self.trimolecular:
                            self.trimolecular_react = np.logical_not(trimolecular_threshold)

                else:  # Assume that all species found in the restart core seed have already been reacted
                    if self.filter_reactions:  # Filling in the filter thresholds will suffice
                        self.unimolecular_threshold[:num_restart_spcs] = True
                        self.bimolecular_threshold[:num_restart_spcs, :num_restart_spcs] = True
                        if self.trimolecular:
                            self.trimolecular_threshold[:num_restart_spcs, :num_restart_spcs, :num_restart_spcs] = True

                    else:  # We must set the react flags instead.
                        # Don't react any species that were present in the restart core seed
                        self.unimolecular_react[:num_restart_spcs] = False
                        self.bimolecular_react[:num_restart_spcs, :num_restart_spcs] = False
                        if self.trimolecular:
                            self.trimolecular_react[:num_restart_spcs, :num_restart_spcs, :num_restart_spcs] = False

    def react_init_tuples(self):
        """
        Reacts tuples given in the react block
        """
        logging.info("Reacting Given Initial Tuples...")
        num_core_species = len(self.reaction_model.core.species)
        self.unimolecular_react = np.zeros((num_core_species), bool)
        self.bimolecular_react = np.zeros((num_core_species, num_core_species), bool)
        if self.trimolecular:
            self.trimolecular_react = np.zeros((num_core_species, num_core_species, num_core_species), bool)

        sts = [spc.label for spc in self.reaction_model.core.species]
        for tup in self.init_react_tuples:
            if len(tup) == 1:
                ind = sts.index(tup[0])
                if not self.unimolecular_threshold[ind]:
                    self.unimolecular_react[ind] = True
                    self.unimolecular_threshold[ind] = True
            elif len(tup) == 2:
                inds = sorted([sts.index(it) for it in tup])
                if not self.bimolecular_threshold[inds[0], inds[1]]:
                    self.bimolecular_react[inds[0], inds[1]] = True
                    self.bimolecular_threshold[inds[0], inds[1]] = True
            elif self.trimolecular and len(tup) == 3:
                inds = sorted([sts.index(it) for it in tup])
                if not self.trimolecular_threshold[inds[0], inds[1], inds[2]]:
                    self.trimolecular_react[inds[0], inds[1], inds[2]] = True
                    self.trimolecular_threshold[inds[0], inds[1], inds[2]] = True
        self.reaction_model.enlarge(react_edge=True,
                                    unimolecular_react=self.unimolecular_react,
                                    bimolecular_react=self.bimolecular_react,
                                    trimolecular_react=self.trimolecular_react)

    def update_reaction_threshold_and_react_flags(self,
                                                  rxn_sys_unimol_threshold=None,
                                                  rxn_sys_bimol_threshold=None,
                                                  rxn_sys_trimol_threshold=None,
                                                  skip_update=False):
        """
        updates the length and boolean value of the unimolecular and bimolecular react and threshold flags
        """
        num_core_species = len(self.reaction_model.core.species)
        prev_num_core_species = len(self.unimolecular_react)
        new_core_species = num_core_species > prev_num_core_species

        # Always reset the react arrays from prior iterations
        self.unimolecular_react = np.zeros((num_core_species), bool)
        self.bimolecular_react = np.zeros((num_core_species, num_core_species), bool)
        if self.trimolecular:
            self.trimolecular_react = np.zeros((num_core_species, num_core_species, num_core_species), bool)

        if self.filter_reactions:
            if new_core_species:
                # Expand the threshold arrays if there were new core species added
                unimolecular_threshold = np.zeros((num_core_species), bool)
                bimolecular_threshold = np.zeros((num_core_species, num_core_species), bool)

                # Broadcast original thresholds
                unimolecular_threshold[:prev_num_core_species] = self.unimolecular_threshold
                bimolecular_threshold[:prev_num_core_species, :prev_num_core_species] = self.bimolecular_threshold
                self.unimolecular_threshold = unimolecular_threshold
                self.bimolecular_threshold = bimolecular_threshold

                if self.trimolecular:
                    trimolecular_threshold = np.zeros((num_core_species, num_core_species, num_core_species), bool)
                    trimolecular_threshold[:prev_num_core_species, :prev_num_core_species, :prev_num_core_species] = self.trimolecular_threshold
                    self.trimolecular_threshold = trimolecular_threshold

            if skip_update:
                return

            # Always update the react and threshold arrays
            for i in range(num_core_species):
                if not self.unimolecular_threshold[i] and rxn_sys_unimol_threshold[i]:
                    # We've shifted from not reacting to reacting
                    self.unimolecular_react[i] = True
                    self.unimolecular_threshold[i] = True

            for i in range(num_core_species):
                for j in range(i, num_core_species):
                    if not self.bimolecular_threshold[i, j] and rxn_sys_bimol_threshold[i, j]:
                        # We've shifted from not reacting to reacting
                        self.bimolecular_react[i, j] = True
                        self.bimolecular_threshold[i, j] = True

            if self.trimolecular:
                for i in range(num_core_species):
                    for j in range(i, num_core_species):
                        for k in range(j, num_core_species):
                            if not self.trimolecular_threshold[i, j, k] and rxn_sys_trimol_threshold[i, j, k]:
                                # We've shifted from not reacting to reacting
                                self.trimolecular_react[i, j, k] = True
                                self.trimolecular_threshold[i, j, k] = True
        else:
            # We are not filtering reactions
            if new_core_species:
                # React all the new core species unimolecularly
                for i in range(prev_num_core_species, num_core_species):
                    self.unimolecular_react[i] = True

                # React all the new core species with all the core species bimolecularly
                for i in range(num_core_species):
                    for j in range(prev_num_core_species, num_core_species):
                        self.bimolecular_react[i, j] = True

                # React all the new core species with all bimolecular combinations trimolecularly
                if self.trimolecular:
                    for i in range(num_core_species):
                        for j in range(num_core_species):
                            for k in range(prev_num_core_species, num_core_species):
                                self.trimolecular_react[i, j, k] = True

    def save_profiler_info(self):
        if self.profiler:  # Save the profile information in case the job crashes
            with open(os.path.join(self.output_directory, 'RMG.profile'), 'wb') as f:
                self.profiler.snapshot_stats()
                marshal.dump(self.profiler.stats, f)

    def save_everything(self):
        """
        Saves the output HTML and the Chemkin file. If the job is being profiled this is saved as well.
        """
        # If the user specifies it, add unused reaction library reactions to
        # an additional output species and reaction list which is written to the output HTML
        # file as well as the chemkin file

        if self.reaction_libraries:
            # First initialize the output_reaction_list and output_species_list to empty
            self.reaction_model.output_species_list = []
            self.reaction_model.output_reaction_list = []
            for library, option in self.reaction_libraries:
                if option:
                    self.reaction_model.add_reaction_library_to_output(library)

        self.exec_time.append(time.time() - self.initialization_time)

        # Notify registered listeners:
        self.notify()

        self.save_profiler_info()

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

    def get_git_commit(self, module_path):
        import subprocess
        if os.path.exists(os.path.join(module_path, '..', '.git')):
            try:
                return subprocess.check_output(['git', 'log',
                                                '--format=%H%n%cd', '-1'],
                                               cwd=module_path).splitlines()
            except:
                return '', ''
        else:
            return '', ''

    def log_header(self, level=logging.INFO):
        """
        Output a header containing identifying information about RMG to the log.
        """
        from rmgpy import __version__, get_path
        logging.log(level, '#########################################################')
        logging.log(level, '# RMG-Py - Reaction Mechanism Generator in Python       #')
        logging.log(level, '# Version: {0:44s} #'.format(__version__))
        logging.log(level, '# Authors: RMG Developers (rmg_dev@mit.edu)             #')
        logging.log(level, '# P.I.s:   William H. Green (whgreen@mit.edu)           #')
        logging.log(level, '#          Richard H. West (r.west@neu.edu)             #')
        logging.log(level, '# Website: http://reactionmechanismgenerator.github.io/ #')
        logging.log(level, '#########################################################\n')

        # Extract git commit from RMG-Py
        head, date = self.get_git_commit(get_path())
        if head != '' and date != '':
            logging.log(level, 'The current git HEAD for RMG-Py is:')
            logging.log(level, '\t%s' % head)
            logging.log(level, '\t%s' % date)
            logging.log(level, '')
        else:
            # If we cannot get git info, try checking if it is a conda package instead:
            conda_package = get_conda_package('rmg')
            if conda_package != '':
                logging.log(level, 'The current anaconda package for RMG-Py is:')
                logging.log(level, conda_package)
                logging.log(level, '')

        database_head, database_date = self.get_git_commit(settings['database.directory'])
        if database_head != '' and database_date != '':
            logging.log(level, 'The current git HEAD for RMG-database is:')
            logging.log(level, '\t%s' % database_head)
            logging.log(level, '\t%s' % database_date)
            logging.log(level, '')
        else:
            database_conda_package = get_conda_package('rmgdatabase')
            if database_conda_package != '':
                logging.log(level, 'The current anaconda package for RMG-database is:')
                logging.log(level, database_conda_package)
                logging.log(level, '')

    def load_rmg_java_input(self, path):
        """
        Load an RMG-Java job from the input file located at `input_file`, or
        from the `input_file` attribute if not given as a parameter.
        """
        warnings.warn("The RMG-Java input is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        # NOTE: This function is currently incomplete!
        # It only loads a subset of the available information.

        self.reaction_model = CoreEdgeReactionModel()
        self.initial_species = []
        self.reaction_systems = []

        T_list = []
        P_list = []
        concentration_list = []
        species_dict = {}
        termination = []

        with open(path, 'r') as f:
            line = self.read_meaningful_line_java(f)
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

                    line = self.read_meaningful_line_java(f)
                    while line != 'END':

                        if line == '' and label != '':
                            species = Species(label=label, molecule=[Molecule().from_adjacency_list(adjlist)])
                            self.initial_species.append(species)
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

                    line = self.read_meaningful_line_java(f)
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
                        self.initial_species.append(species)
                        species_dict[label] = species
                        concentration_list.append(concentrations)

                        line = self.read_meaningful_line_java(f)

                elif line.startswith('FinishController:'):

                    # First meaningful line is a termination time or conversion
                    line = self.read_meaningful_line_java(f)
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
                    line = self.read_meaningful_line_java(f)

                line = self.read_meaningful_line_java(f)

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
                    initial_mole_fractions = dict([(self.initial_species[i], concentrations[i] / total_conc) for i in
                                                   range(len(self.initial_species))])
                    reaction_system = SimpleReactor(T, P, initial_mole_fractions=initial_mole_fractions,
                                                   termination=termination)
                    self.reaction_systems.append(reaction_system)

    def read_meaningful_line_java(self, f):
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

def determine_procnum_from_ram():
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


def initialize_log(verbose, log_file_name):
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

    def __init__(self, reaction_system, bspc):
        self.Ranges = dict()

        if hasattr(reaction_system, 'Trange') and isinstance(reaction_system.Trange, list):
            Trange = reaction_system.Trange
            self.Ranges['T'] = [T.value_si for T in Trange]
        if hasattr(reaction_system, 'Prange') and isinstance(reaction_system.Prange, list):
            Prange = reaction_system.Prange
            self.Ranges['P'] = [np.log(P.value_si) for P in Prange]
        if hasattr(reaction_system, 'initial_mole_fractions'):
            if bspc:
                self.initial_mole_fractions = deepcopy(reaction_system.initial_mole_fractions)
                self.balance_species = [x for x in self.initial_mole_fractions.keys() if x.label == bspc][
                    0]  # find the balance species
            for key, value in reaction_system.initial_mole_fractions.items():
                assert key != 'T' and key != 'P', 'naming a species T or P is forbidden'
                if isinstance(value, list):
                    self.Ranges[key] = value
        if hasattr(reaction_system, 'initial_concentrations'):
            for key, value in reaction_system.initial_concentrations.items():
                assert key != 'T' and key != 'P', 'naming a species T or P is forbidden'
                if isinstance(value, list):
                    self.Ranges[key] = [v.value_si for v in value]
        
        if isinstance(reaction_system, Reactor):
            self.tmax = reaction_system.tf
        else:
            for term in reaction_system.termination:
                if isinstance(term, TerminationTime):
                    self.tmax = term.time.value_si

        self.reaction_system = reaction_system
        self.condition_list = []
        self.scaled_condition_list = []
        self.ts = []
        self.convs = []
        self.Ns = []
        self.rand_state = np.random.RandomState(1)

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
            return self.condition_list[-1]

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
        n = self.rand_state.uniform(0, 1, 1)[0]  # draw a random number between 0 and 1
        s = 0.0
        for indexes in np.ndenumerate(Jout):  # choose a coordinate such that grid[indexes] is choosen with probability Jout[indexes]
            s += Jout[indexes[0]]
            if s > n:
                break
        if len(bounds) != 1:
            yf = np.array([grid[i][indexes[0]] for i in range(len(grid))])
        else:
            yf = np.array([grid[indexes[0]] for i in range(len(grid))])

        step = self.rand_state.uniform(0, 1, len(Jout.shape))  # take a step in a random direction in a length between 0 and 1/(2*Ns)
        step /= step.sum()
        mag = self.rand_state.uniform(0, 1, 1)[0]

        yf += step * mag * np.sqrt(2) / (2.0 * Ns)

        return yf

    def generate_cond(self):
        """
        find the next condition to run at by solving an optimization problem
        this optimization problem maximizes distance from prior conditions weighted more if they are more recent
        and maximizes number of objects added
        the resulting condition is added to the end of condition_list
        """
        if self.condition_list == []:
            self.condition_list.append({key: value[0] for key, value in self.Ranges.items()})
            self.scaled_condition_list.append({key: 0.0 for key, value in self.Ranges.items()})
        elif len(self.condition_list[0]) == 0:
            pass
        else:
            ykey = list(self.condition_list[0].keys())
            Ns = self.Ns

            def obj(y):
                boo = y.shape == tuple()
                vec = []
                N = len(self.condition_list)
                for i, cond in enumerate(self.scaled_condition_list):
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

            if hasattr(self, 'initial_mole_fractions'):
                for key in self.initial_mole_fractions.keys():
                    if not isinstance(self.initial_mole_fractions[key], list):
                        new_cond[key] = self.initial_mole_fractions[key]
                total = sum([val for key, val in new_cond.items() if key != 'T' and key != 'P'])
                if self.balance_species is None:
                    for key, val in new_cond.items():
                        if key != 'T' and key != 'P':
                            new_cond[key] = val / total
                else:
                    new_cond[self.balance_species] = self.initial_mole_fractions[self.balance_species] + 1.0 - total

            self.condition_list.append(new_cond)
            self.scaled_condition_list.append(scaled_new_cond)
        return


def log_conditions(rmg_memories, index):
    """
    log newly generated reactor conditions
    """
    if rmg_memories[index].get_cond() is not None:
        s = 'conditions choosen for reactor {0} were: '.format(index)
        for key, item in rmg_memories[index].get_cond().items():
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


def get_conda_package(module):
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


def process_profile_stats(stats_file, log_file):
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


def make_profile_graph(stats_file, force_graph_generation=False):
    """
    Uses gprof2dot to create a graphviz dot file of the profiling information.

    This requires the gprof2dot package available via `pip install gprof2dot`.
    Render the result using the program 'dot' via a command like
    `dot -Tps2 input.dot -o output.ps2`.

    Rendering the ps2 file to pdf requires an external pdf converter
    `ps2pdf output.ps2` which produces a `output.ps2.pdf` file.

    Will only generate a graph if a display is present as errors can occur otherwise. If `force_graph_generation` is
    True then the graph generation will be attempted either way
    """
    # Making the profile graph requires a display. See if one is available first
    display_found = False

    try:
        display_found = bool(os.environ['DISPLAY'])
    except KeyError:  # This means that no display was found
        pass

    if display_found or force_graph_generation:
        try:
            from gprof2dot import PstatsParser, DotWriter, SAMPLES, themes, TIME, TIME_RATIO, TOTAL_TIME, TOTAL_TIME_RATIO
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

        # Add both total time and self time in seconds to the graph output
        dot.show_function_events = [TOTAL_TIME, TOTAL_TIME_RATIO, TIME, TIME_RATIO]

        if options.show_samples:
            dot.show_function_events.append(SAMPLES)

        profile = profile
        profile.prune(options.node_thres / 100.0, options.edge_thres / 100.0, [], False)

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

    else:
        logging.warning('Could not find a display, which is required in order to generate the profile graph. This '
                        'is likely due to this job being run on a remote server without performing X11 forwarding '
                        'or running the job through a job manager like SLURM.\n\n The graph can be generated later '
                        'by running with the postprocessing flag `rmg.py -P input.py` from any directory/computer '
                        'where both the input file and RMG.profile file are located and a display is available.\n\n'
                        'Note that if the postprocessing flag is specified, this will force the graph generation '
                        'regardless of if a display was found, which could cause this program to crash or freeze.')
