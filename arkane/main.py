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
This module contains the :class:`Arkane` class, the main class used to run Arkane.
"""

import argparse
import csv
import logging
import os
import os.path
import subprocess
import sys
import time

import numpy as np

try:
    import matplotlib
    matplotlib.rc('mathtext', default='regular')
except ImportError:
    pass

from rmgpy import __version__, get_path, settings
from rmgpy.chemkin import write_elements_section
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import InputError

from arkane.common import is_pdep
from arkane.encorr.ae import AEJob
from arkane.encorr.bac import BACJob
from arkane.explorer import ExplorerJob
from arkane.input import load_input_file
from arkane.kinetics import KineticsJob
from arkane.output import save_thermo_lib, save_kinetics_lib
from arkane.pdep import PressureDependenceJob
from arkane.statmech import StatMechJob
from arkane.thermo import ThermoJob

################################################################################


class Arkane(object):
    """
    The :class:`Arkane` class represents an instance of Arkane, a tool for
    computing properties of chemical species and reactions. The attributes are:
    
    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `job_list`          A list of the jobs to execute
    `input_file`        The path of the input file defining the jobs to execute
    `output_directory`  The directory in which to write the output files
    `verbose`           The level of detail in the generated logging messages
    =================== ========================================================
    
    The output directory defaults to the same directory as the input file if
    not explicitly specified.
    
    To use this class programmatically, create an instance and set its
    attributes using either the :meth:`__init__()` method or by directly
    accessing the attributes, and then invoke the :meth:`execute()` method.
    You can also populate the attributes from the command line using the
    :meth:`parse_command_line_arguments()` method before running :meth:`execute()`.
    """

    def __init__(self, input_file=None, output_directory=None, verbose=logging.INFO, save_rmg_libraries=True):
        self.job_list = []
        self.input_file = input_file
        self.output_directory = output_directory
        self.verbose = verbose
        self.save_rmg_libraries = save_rmg_libraries

    def parse_command_line_arguments(self):
        """
        Parse the command-line arguments being passed to Arkane. This uses the
        :mod:`argparse` module, which ensures that the command-line arguments are
        sensible, parses them, and returns them.
        """

        parser = argparse.ArgumentParser(description="""
        Arkane is a Python toolkit for computing chemical reaction rates
        and other properties used in detailed kinetics models
        using various methodologies and theories.
        """)
        parser.add_argument('file', metavar='FILE', type=str, nargs=1, help='a file describing the job to execute')

        # Options for controlling the amount of information printed to the console
        # By default a moderate level of information is printed; you can either
        # ask for less (quiet), more (verbose), or much more (debug)
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-q', '--quiet', action='store_const', const=logging.WARNING, default=logging.INFO,
                           dest='verbose', help='only print warnings and errors')
        group.add_argument('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO,
                           dest='verbose', help='print more verbose output')
        group.add_argument('-d', '--debug', action='store_const', const=0, default=logging.INFO, dest='verbose',
                           help='print debug information')

        # Add options for controlling what directories files are written to
        parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
                            metavar='DIR', help='use DIR as output directory')

        # Add options for controlling generation of plots
        parser.add_argument('-p', '--no-plot', action='store_false', default=True,
                            help='prevent generating plots', dest='plot')

        args = parser.parse_args()

        # Extract the input file
        self.input_file = args.file[0]

        # Extract the log verbosity
        self.verbose = args.verbose

        # Extract the plot settings
        self.plot = args.plot

        # Determine the output directory
        # By default the directory containing the input file is used, unless an
        # alternate directory is specified using the -o flag
        if args.output_directory and os.path.isdir(args.output_directory[0]):
            self.output_directory = os.path.abspath(args.output_directory[0])
        else:
            self.output_directory = os.path.dirname(os.path.abspath(args.file[0]))

    def load_input_file(self, input_file):
        """
        Load a set of jobs from the given `input_file` on disk. Returns the
        loaded set of jobs as a list.
        """
        self.input_file = input_file
        self.job_list, self.reaction_dict, self.species_dict, self.transition_state_dict, self.network_dict, \
            self.level_of_theory = load_input_file(self.input_file)
        logging.info('')
        return self.job_list

    def execute(self):
        """
        Execute, in order, the jobs found in input file specified by the
        `input_file` attribute.
        """

        # Initialize the logging system (both to the console and to a file in the
        # output directory)
        initialize_log(self.verbose, os.path.join(self.output_directory, 'arkane.log'))

        # Print some information to the beginning of the log
        log_header()

        # Load the input file for the job
        self.job_list = self.load_input_file(self.input_file)
        logging.info('')

        # Initialize (and clear!) the output files for the job
        if self.output_directory is None:
            self.output_directory = os.path.dirname(os.path.abspath(self.input_file))
        output_file = os.path.join(self.output_directory, 'output.py')
        with open(output_file, 'w'):
            pass
        chemkin_file = os.path.join(self.output_directory, 'chem.inp')

        # write the chemkin files and run the thermo and then kinetics jobs
        with open(chemkin_file, 'w') as f:
            write_elements_section(f)

            f.write('SPECIES\n\n')

            # write each species in species block
            for job in self.job_list:
                if isinstance(job, ThermoJob):
                    f.write(job.species.to_chemkin())
                    f.write('\n')

            f.write('\nEND\n\n\n\n')
            f.write('THERM ALL\n')
            f.write('    300.000  1000.000  5000.000\n\n')

        # run thermo and statmech jobs (also writes thermo blocks to Chemkin file)
        supporting_info = []
        hindered_rotor_info = []
        bacjob_num = 1
        for job in self.job_list:
            if isinstance(job, ThermoJob):
                job.execute(output_directory=self.output_directory, plot=self.plot)
            if isinstance(job, StatMechJob):
                job.execute(output_directory=self.output_directory, plot=self.plot, pdep=is_pdep(self.job_list))
                if hasattr(job, 'supporting_info'):
                    supporting_info.append(job.supporting_info)
                if hasattr(job, 'raw_hindered_rotor_data'):
                    for hr_info in job.raw_hindered_rotor_data:
                        hindered_rotor_info.append(hr_info)
            if isinstance(job, BACJob):
                job.execute(output_directory=self.output_directory, plot=self.plot, jobnum=bacjob_num)
                bacjob_num += 1
            if isinstance(job, AEJob):
                job.execute(output_file=output_file)

        with open(chemkin_file, 'a') as f:
            f.write('\n')
            f.write('END\n\n\n\n')
            f.write('REACTIONS    KCAL/MOLE   MOLES\n\n')

        if supporting_info:
            # write supporting_info.csv for statmech jobs
            supporting_info_file = os.path.join(self.output_directory, 'supporting_information.csv')
            with open(supporting_info_file, 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(['Label', 'Symmetry Number', 'Number of optical isomers', 'Symmetry Group',
                                 'Rotational constant (cm-1)',
                                 'Calculated Frequencies (unscaled and prior to projection, cm^-1)',
                                 'Electronic energy (J/mol)', 'E0 (electronic energy + ZPE, J/mol)',
                                 'E0 with atom and bond corrections (J/mol)', 'Atom XYZ coordinates (angstrom)',
                                 'T1 diagnostic', 'D1 diagnostic'])
                for row in supporting_info:
                    label = row[0]
                    rot = '-'
                    freq = '-'
                    if row[4] is not None and isinstance(row[4].rotationalConstant.value, float):
                        # diatomic species have a single rotational constant
                        rot = '{0:.2f}'.format(row[4].rotationalConstant.value)
                    elif row[4] is not None:
                        rot = ', '.join(['{0:.2f}'.format(s) for s in row[4].rotationalConstant.value])
                    if row[5] is not None:
                        freq = ''
                        if row[6] is not None:  # there is a negative frequency
                            freq = '{0:.1f}'.format(abs(row[6])) + 'i, '
                        freq += ', '.join(['{0:.1f}'.format(s) for s in row[5]])
                    atoms = ', '.join(["{0}    {1}".format(atom, "    ".join([str(c) for c in coords]))
                                       for atom, coords in zip(row[10], row[11])])
                    writer.writerow([label, row[1], row[2], row[3], rot, freq, row[7], row[8], row[9], atoms, row[12],
                                     row[13]])
        if hindered_rotor_info:
            hr_file = os.path.join(self.output_directory, 'hindered_rotor_scan_data.csv')
            # find longest length to set column number for energies
            max_energy_length = max([len(hr[4]) for hr in hindered_rotor_info])
            with open(hr_file, 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(['species', 'rotor_number', 'symmetry', 'resolution (degrees)',
                                 'pivot_atoms', 'frozen_atoms'] +
                                ['energy (J/mol) {}'.format(i) for i in range(max_energy_length)])
                for row in hindered_rotor_info:
                    writer.writerow([row[0], row[1], row[2], row[3][1] * 180 / np.pi,
                                     row[5], row[6]] + [a for a in row[4]])
        # run kinetics and pdep jobs (also writes reaction blocks to Chemkin file)
        for job in self.job_list:
            if isinstance(job, KineticsJob):
                job.execute(output_directory=self.output_directory, plot=self.plot)
            elif isinstance(job, PressureDependenceJob) and not any([isinstance(job, ExplorerJob) for job in
                                                                     self.job_list]):
                # if there is an explorer job the pdep job will be run in the explorer job
                if job.network is None:
                    raise InputError(
                        'No network matched the label of the pressureDependence block and there is no explorer block '
                        'to generate a network')
                job.execute(output_file=output_file, plot=self.plot)
            elif isinstance(job, ExplorerJob):
                thermo_library, kinetics_library, species_list = self.get_libraries()
                job.execute(output_file=output_file, plot=self.plot, species_list=species_list,
                            thermo_library=thermo_library, kinetics_library=kinetics_library)

        with open(chemkin_file, 'a') as f:
            f.write('END\n\n')

        # Print some information to the end of the log
        log_footer()

        if self.save_rmg_libraries:
            # save RMG thermo and kinetics libraries
            species, reactions = list(), list()
            for job in self.job_list:
                if isinstance(job, ThermoJob) and len(job.species.molecule):
                    species.append(job.species)
                elif isinstance(job, KineticsJob) \
                        and all([len(species.molecule) for species in job.reaction.reactants + job.reaction.products]):
                    reactions.append(job.reaction)
                elif isinstance(job, PressureDependenceJob):
                    for reaction in job.network.path_reactions:
                        if all([len(species.molecule) for species in reaction.reactants + reaction.products]):
                            reactions.append(reaction)
            lib_path = os.path.join(self.output_directory, 'RMG_libraries')
            level_of_theory = f' using {self.level_of_theory}' if self.level_of_theory is not None else ''
            lib_long_desc = f'Calculated using Arkane v{__version__}{level_of_theory}.'
            save_thermo_lib(species_list=species, path=lib_path, name='thermo', lib_long_desc=lib_long_desc)
            save_kinetics_lib(rxn_list=reactions, path=lib_path, name='kinetics', lib_long_desc=lib_long_desc)

    def get_libraries(self):
        """Get RMG kinetics and thermo libraries"""
        name = 'kineticsjobs'

        species_list = list(self.species_dict.values())
        reaction_list = list(self.reaction_dict.values())

        # remove duplicate species
        for rxn in reaction_list:
            for i, rspc in enumerate(rxn.reactants):
                for spc in species_list:
                    if spc.is_isomorphic(rspc):
                        rxn.reactants[i] = spc
                        break
            for i, rspc in enumerate(rxn.products):
                for spc in species_list:
                    if spc.is_isomorphic(rspc):
                        rxn.products[i] = spc
                        break
        del_inds = []
        for i, spc1 in enumerate(species_list):
            for j, spc2 in enumerate(species_list):
                if j > i and spc1.is_isomorphic(spc2):
                    del_inds.append(j)

        for j in sorted(del_inds)[::-1]:
            del species_list[j]

        thermo_library = ThermoLibrary(name=name)
        for i, species in enumerate(species_list):
            if species.thermo:
                thermo_library.load_entry(index=i + 1,
                                          label=species.label,
                                          molecule=species.molecule[0].to_adjacency_list(),
                                          thermo=species.thermo,
                                          shortDesc=species.thermo.comment)
            else:
                logging.warning(
                    'Species {0} did not contain any thermo data and was omitted from the thermo library.'.format(
                        str(species)))

        # load kinetics library entries                    
        kinetics_library = KineticsLibrary(name=name, auto_generated=True)
        kinetics_library.entries = {}
        for i, reaction in enumerate(reaction_list):
            entry = Entry(
                index=i + 1,
                label=reaction.to_labeled_str(),
                item=reaction,
                data=reaction.kinetics)

            if reaction.kinetics is not None:
                if hasattr(reaction, 'library') and reaction.library:
                    entry.long_desc = 'Originally from reaction library: ' + \
                                      reaction.library + "\n" + reaction.kinetics.comment
                else:
                    entry.long_desc = reaction.kinetics.comment

            kinetics_library.entries[i + 1] = entry

        kinetics_library.label = name

        return thermo_library, kinetics_library, species_list


def initialize_log(verbose=logging.INFO, log_file=None):
    """
    Set up a logger for Arkane to use to print output to stdout. The
    `verbose` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(verbose)

    # Use custom level names for cleaner log output
    logging.addLevelName(logging.CRITICAL, 'Critical: ')
    logging.addLevelName(logging.ERROR, 'Error: ')
    logging.addLevelName(logging.WARNING, 'Warning: ')
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, '')
    logging.addLevelName(0, '')

    # Create formatter and add to handlers
    formatter = logging.Formatter('%(levelname)s%(message)s')

    # Remove old handlers before adding ours
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Create console handler; send everything to stdout rather than stderr
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(verbose)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Create file handler; always be at least verbose in the file
    if log_file:
        fh = logging.FileHandler(filename=log_file)
        fh.setLevel(min(logging.DEBUG, verbose))
        fh.setFormatter(formatter)
        logger.addHandler(fh)


def get_git_commit(path):
    """
    Get the recent git commit to be logged.
    """
    head, date = '', ''
    if os.path.exists(os.path.join(path, '..', '.git')):
        try:
            head, date = subprocess.check_output(['git', 'log', '--format=%H%n%cd', '-1'], cwd=path).splitlines()
            head, date = head.decode(), date.decode()
        except (subprocess.CalledProcessError, OSError):
            return head, date
    return head, date


def get_conda_package(module):
    """
    Check the version of any conda package.
    """
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


def log_header(level=logging.INFO):
    """
    Output a header containing identifying information about Arkane to the log.
    """
    logging.log(level, 'Arkane execution initiated at {0}'.format(time.asctime()))
    logging.log(level, '')

    logging.log(level, '################################################################')
    logging.log(level, '#                                                              #')
    logging.log(level, '# Automated Reaction Kinetics and Network Exploration (Arkane) #')
    logging.log(level, '#                                                              #')
    logging.log(level, '#   Version: {0:49s} #'.format(__version__))
    logging.log(level, '#   Authors: RMG Developers (rmg_dev@mit.edu)                  #')
    logging.log(level, '#   P.I.s:   William H. Green (whgreen@mit.edu)                #')
    logging.log(level, '#            Richard H. West (r.west@neu.edu)                  #')
    logging.log(level, '#   Website: http://reactionmechanismgenerator.github.io/      #')
    logging.log(level, '#                                                              #')
    logging.log(level, '################################################################')
    logging.log(level, '')

    # Extract git commit from RMG-Py
    head, date = get_git_commit(get_path())
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

    # Extract git commit from RMG-database
    database_head, database_date = get_git_commit(settings['database.directory'])
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


def log_footer(level=logging.INFO):
    """
    Output a footer to the log.
    """
    logging.log(level, '')
    logging.log(level, 'Arkane execution terminated at {0}'.format(time.asctime()))
