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
This module contains the :class:`Arkane` class, the main class used to run Arkane.
"""

import os
import os.path
import sys
import logging
import argparse
import time
import csv
try:
    import matplotlib
    matplotlib.rc('mathtext', default='regular')
except ImportError:
    pass

from rmgpy.chemkin import writeElementsSection
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import InputError

from arkane.input import loadInputFile
from arkane.kinetics import KineticsJob
from arkane.statmech import StatMechJob
from arkane.thermo import ThermoJob
from arkane.pdep import PressureDependenceJob
from arkane.explorer import ExplorerJob

################################################################################


class Arkane:
    """
    The :class:`Arkane` class represents an instance of Arkane, a tool for
    computing properties of chemical species and reactions. The attributes are:
    
    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `jobList`           A list of the jobs to execute
    `inputFile`         The path of the input file defining the jobs to execute
    `outputDirectory`   The directory in which to write the output files
    `verbose`           The level of detail in the generated logging messages
    =================== ========================================================
    
    The output directory defaults to the same directory as the input file if
    not explicitly specified.
    
    To use this class programmatically, create an instance and set its
    attributes using either the :meth:`__init__()` method or by directly
    accessing the attributes, and then invoke the :meth:`execute()` method.
    You can also populate the attributes from the command line using the
    :meth:`parseCommandLineArguments()` method before running :meth:`execute()`.
    """
    
    def __init__(self, inputFile=None, outputDirectory=None, verbose=logging.INFO):
        self.jobList = []
        self.inputFile = inputFile
        self.outputDirectory = outputDirectory
        self.verbose = verbose
    
    def parseCommandLineArguments(self):
        """
        Parse the command-line arguments being passed to Arkane. This uses the
        :mod:`argparse` module, which ensures that the command-line arguments are
        sensible, parses them, and returns them.
        """
    
        parser = argparse.ArgumentParser(description=
        """
        Arkane is a Python toolkit for computing chemical reaction rates and other
        properties used in detailed kinetics models using various methodologies
        and theories.
        """)
        parser.add_argument('file', metavar='FILE', type=str, nargs=1,
            help='a file describing the job to execute')
    
        # Options for controlling the amount of information printed to the console
        # By default a moderate level of information is printed; you can either
        # ask for less (quiet), more (verbose), or much more (debug)
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-q', '--quiet', action='store_const', const=logging.WARNING, default=logging.INFO, dest='verbose', help='only print warnings and errors')
        group.add_argument('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO, dest='verbose', help='print more verbose output')
        group.add_argument('-d', '--debug', action='store_const', const=0, default=logging.INFO, dest='verbose', help='print debug information')
    
        # Add options for controlling what directories files are written to
        parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
            metavar='DIR', help='use DIR as output directory')

        # Add options for controlling generation of plots
        parser.add_argument('-p', '--plot', action='store_true', default=True, help='generate plots of results')

        args = parser.parse_args()
        
        # Extract the input file
        self.inputFile = args.file[0] 
        
        # Extract the log verbosity
        self.verbose = args.verbose
        
        # Extract the plot settings
        self.plot = args.plot
        
        # Determine the output directory
        # By default the directory containing the input file is used, unless an
        # alternate directory is specified using the -o flag
        if args.output_directory and os.path.isdir(args.output_directory[0]):
            self.outputDirectory = os.path.abspath(args.output_directory[0])
        else:
            self.outputDirectory = os.path.dirname(os.path.abspath(args.file[0]))
    
    def initializeLog(self, verbose=logging.INFO, logFile=None):
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
        if logFile:
            fh = logging.FileHandler(filename=logFile)
            fh.setLevel(min(logging.DEBUG,verbose))
            fh.setFormatter(formatter)
            logger.addHandler(fh)
    
    def logHeader(self, level=logging.INFO):
        """
        Output a header containing identifying information about Arkane to the log.
        """
        from rmgpy import __version__
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
    
    def logFooter(self, level=logging.INFO):
        """
        Output a footer to the log.
        """
        logging.log(level, '')
        logging.log(level, 'Arkane execution terminated at {0}'.format(time.asctime()))

    def loadInputFile(self, inputFile):
        """
        Load a set of jobs from the given `inputFile` on disk. Returns the
        loaded set of jobs as a list.
        """
        self.inputFile = inputFile
        self.jobList, self.reactionDict, self.speciesDict, self.transitionStateDict, self.networkDict = loadInputFile(self.inputFile)
        logging.info('')
        return self.jobList
        
    def execute(self):
        """
        Execute, in order, the jobs found in input file specified by the
        `inputFile` attribute.
        """
        
        # Initialize the logging system (both to the console and to a file in the
        # output directory)
        self.initializeLog(self.verbose, os.path.join(self.outputDirectory, 'arkane.log'))
        
        # Print some information to the beginning of the log
        self.logHeader()
        
        # Load the input file for the job
        self.jobList = self.loadInputFile(self.inputFile)
        logging.info('')
        
        # Initialize (and clear!) the output files for the job
        if self.outputDirectory is None:
            self.outputDirectory = os.path.dirname(os.path.abspath(self.inputFile))
        outputFile = os.path.join(self.outputDirectory, 'output.py')
        with open(outputFile, 'w') as f:
            pass
        chemkinFile = os.path.join(self.outputDirectory, 'chem.inp')

        # write the chemkin files and run the thermo and then kinetics jobs
        with open(chemkinFile, 'w') as f:
            writeElementsSection(f)
            
            f.write('SPECIES\n\n')

            # write each species in species block
            for job in self.jobList:
                if isinstance(job,ThermoJob):
                    f.write(job.species.toChemkin())
                    f.write('\n')

            f.write('\nEND\n\n\n\n')
            f.write('THERM ALL\n')
            f.write('    300.000  1000.000  5000.000\n\n')

        # run thermo and statmech jobs (also writes thermo blocks to Chemkin file)
        supporting_info = []
        for job in self.jobList:
            if isinstance(job, ThermoJob):
                job.execute(outputFile=outputFile, plot=self.plot)
            if isinstance(job, StatMechJob):
                job.execute(outputFile=outputFile, plot=self.plot)
                supporting_info.append(job.supporting_info)

        with open(chemkinFile, 'a') as f:
            f.write('\n')
            f.write('END\n\n\n\n')
            f.write('REACTIONS    KCAL/MOLE   MOLES\n\n')

        supporting_info_file = os.path.join(self.outputDirectory, 'supporting_information.csv')
        with open(supporting_info_file, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['Label','Rotational constant (cm-1)','Unscaled frequencies (cm-1)'])
            for row in supporting_info:
                label = row[0]
                rot = '-'
                freq = '-'
                if len(row) > 1:  # monoatomic species have no frequencies nor rotational constants
                    if isinstance(row[1].rotationalConstant.value, float):
                        # diatomic species have a single rotational constant
                        rot = '{0:.2f}'.format(row[1].rotationalConstant.value)
                    else:
                        rot = ', '.join(['{0:.2f}'.format(s) for s in row[1].rotationalConstant.value])
                    freq = ''
                    if len(row) == 4:
                        freq = '{0:.1f}'.format(abs(row[3])) + 'i, '
                    freq += ', '.join(['{0:.1f}'.format(s) for s in row[2]])
                writer.writerow([label, rot, freq])

        # run kinetics and pdep jobs (also writes reaction blocks to Chemkin file)
        for job in self.jobList:
            if isinstance(job,KineticsJob):
                job.execute(outputFile=outputFile, plot=self.plot)
            elif isinstance(job, PressureDependenceJob) and not any([isinstance(job,ExplorerJob) for job in self.jobList]): #if there is an explorer job the pdep job will be run in the explorer job
                if job.network is None:
                    raise InputError('No network matched the label of the pressureDependence block and there is no explorer block to generate a network')
                job.execute(outputFile=outputFile, plot=self.plot)
            elif isinstance(job, ExplorerJob):
                thermoLibrary,kineticsLibrary,speciesList = self.getLibraries()
                job.execute(outputFile=outputFile, plot=self.plot, speciesList=speciesList, thermoLibrary=thermoLibrary, kineticsLibrary=kineticsLibrary)

        with open(chemkinFile, 'a') as f:
            f.write('END\n\n')

        # Print some information to the end of the log
        self.logFooter()
    
    def getLibraries(self):

        name = 'kineticsjobs'
                
        speciesList = self.speciesDict.values()
        reactionList = self.reactionDict.values()

        # remove duplicate species
        for rxn in reactionList:
            for i,rspc in enumerate(rxn.reactants):
                for spc in speciesList:
                    if spc.isIsomorphic(rspc):
                        rxn.reactants[i] = spc
                        break
            for i,rspc in enumerate(rxn.products):
                for spc in speciesList:
                    if spc.isIsomorphic(rspc):
                        rxn.products[i] = spc
                        break
        del_inds = []
        for i,spc1 in enumerate(speciesList):
            for j,spc2 in enumerate(speciesList):
                if j>i and spc1.isIsomorphic(spc2):
                    del_inds.append(j)
        
        for j in sorted(del_inds)[::-1]:
            del speciesList[j]
            
        thermoLibrary = ThermoLibrary(name=name)
        for i,species in enumerate(speciesList): 
            if species.thermo:
                thermoLibrary.loadEntry(index = i + 1,
                                        label = species.label,
                                        molecule = species.molecule[0].toAdjacencyList(),
                                        thermo = species.thermo,
                                        shortDesc = species.thermo.comment
               )                
            else:
                logging.warning('Species {0} did not contain any thermo data and was omitted from the thermo library.'.format(str(species)))

        # load kinetics library entries                    
        kineticsLibrary = KineticsLibrary(name=name,autoGenerated=True)
        kineticsLibrary.entries = {}
        for i,reaction in enumerate(reactionList):      
            entry = Entry(
                    index = i+1,
                    label = reaction.toLabeledStr(),
                    item = reaction,
                    data = reaction.kinetics,
                )

            if reaction.kinetics is not None:
                if hasattr(reaction,'library') and reaction.library:
                    entry.longDesc = 'Originally from reaction library: ' +\
                                     reaction.library + "\n" + reaction.kinetics.comment
                else:
                    entry.longDesc = reaction.kinetics.comment
            
            kineticsLibrary.entries[i+1] = entry
        
        kineticsLibrary.label = name
        
        return thermoLibrary,kineticsLibrary,speciesList
