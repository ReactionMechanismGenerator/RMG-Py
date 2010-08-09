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

"""
RMG is an automatic chemical mechanism generator. It is awesomely awesome.
"""

import os.path
import sys
import argparse
import logging
import time

import rmg.settings as settings
from rmg.input import readInputFile

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to RMG Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser(description=
    """
    Reaction Mechanism Generator (RMG) is an automatic chemical reaction
    mechanism generator that constructs kinetic models composed of
    elementary chemical reaction steps using a general understanding of
    how molecules react.
    """)
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
        help='a file describing the job to execute')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')

    # Add options for controlling what directories files are written to
    parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as output directory')
    parser.add_argument('-s', '--scratch-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as scratch directory')
    parser.add_argument('-l', '--library-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as library directory')

    # Add restart option
    parser.add_argument('-r', '--restart', action='store_true', help='restart an incomplete job')

    parser.add_argument('-p', '--profile', action='store_true', help='run under cProfile to gather profiling statistics, and postprocess them if job completes')
    parser.add_argument('-P', '--postprocess', action='store_true', help='postprocess profiling statistics from previous [failed] run; does not run the simulation')

    parser.add_argument('-t', '--walltime', type=str, nargs=1, default='0',
        metavar='HH:MM:SS', help='set the maximum execution time')

    return parser.parse_args()

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
    logging.addLevelName(0, '')

    # Create formatter and add to console handler
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    #formatter = Formatter('%(message)s', '%Y-%m-%d %H:%M:%S')
    formatter = logging.Formatter('%(levelname)s%(message)s')
    ch.setFormatter(formatter)

    # create file handler
    if os.path.exists(log_file_name):
        backup_name = '_backup'.join(os.path.splitext(log_file_name))
        if os.path.exists(backup_name):
            print "Removing old file %s" % backup_name
            os.remove(backup_name)
        print "Renaming %s to %s"%(log_file_name, backup_name)
        print
        os.rename(log_file_name, backup_name)
    fh = logging.FileHandler(filename=log_file_name) #, backupCount=3)
    fh.setLevel(min(logging.DEBUG,verbose)) # always at least VERBOSE in the file
    fh.setFormatter(formatter)
    # notice that STDERR does not get saved to the log file
    # so errors from underlying libraries (eg. openbabel) etc. that report
    # on stderr will not be logged to disk.

    # remove old handlers!
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Add ch to logger
    logger.addHandler(ch)
    logger.addHandler(fh)

def logHeader(level=logging.INFO):
    """
    Output a header containing identifying information about RMG to the log.
    """

    logging.log(level, '#################################################')
    logging.log(level, '# RMG - Reaction Mechanism Generator            #')
    logging.log(level, '# Version: 0.1.0 (14 May 2009)                  #')
    logging.log(level, '# Authors: RMG Developers (rmg_dev@mit.edu)     #')
    logging.log(level, '# P.I.:    William H. Green (whgreen@mit.edu)   #')
    logging.log(level, '# Website: http://rmg.sourceforge.net/          #')
    logging.log(level, '#################################################\n')

    import os

    def getGitCommit():
        try:
            f = os.popen('git log --format="%H %n %cd" -1')
            lines = []
            for line in f: lines.append(line)
            f.close()
            head = lines[0].strip()
            date = lines[1].strip()
            return head, date
        except IndexError:
            return '', ''

    head, date = getGitCommit()
    if head != '' and date != '':
        logging.log(level, 'The current git HEAD is:')
        logging.log(level, '\t%s' % head)
        logging.log(level, '\t%s' % date)

    logging.log(level, '')

################################################################################

def execute(args):
    """
    Generate a reaction model for the set of reaction systems specified in an
    input file at location `inputFile`. Output and temporary files will be
    placed in the directories `outputDir` and `scratchDir`, respectively; if
    either of these are empty strings, the corresponding files will be placed
    in the same directory as the input file. The `libraryDir` parameter can be
    used to specify additional directories to search for RMG databases. The
    `verbose` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """

    inputFile = args.file[0]

    # Set directories
    settings.outputDirectory = args.output_directory
    settings.scratchDirectory = args.scratch_directory
    settings.libraryDirectory = args.library_directory

    # Set wall time
    if args.walltime == '0': settings.wallTime = 0
    else:
        data = args.walltime.split(':')
        if len(data) == 1:
            settings.walltime = int(data[-1])
        elif len(data) == 2:
            settings.walltime = int(data[-1]) + 60 * int(data[-2])
        elif len(data) == 3:
            settings.walltime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3])
        elif len(data) == 4:
            settings.walltime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3]) + 86400 * int(data[-4])
        else:
            raise ValueError('Invalid format for wall time; should be HH:MM:SS.')

    # Save initialization time
    settings.initializationTime = time.time()

    # Log start timestamp
    logging.info('RMG execution initiated at ' + time.asctime() + '\n')

    # Print out RMG header
    logHeader()

    # Make output subdirectories
    plotDir = os.path.join(settings.outputDirectory, 'plot')
    if os.path.exists(plotDir):
        for f in os.listdir(plotDir):
            os.remove(plotDir + '/' + f)
        os.rmdir(plotDir)
    os.mkdir(plotDir)

    specDir = os.path.join(settings.outputDirectory, 'species')
    if os.path.exists(specDir):
        for f in os.listdir(specDir):
            os.remove(specDir + '/' + f)
        os.rmdir(specDir)
    os.mkdir(specDir)

    # Read input file
    reactionModel, coreSpecies, reactionSystems = readInputFile(inputFile)

    # Initialize reaction model
    if args.restart:
        import gzip
        import cPickle
        import ctml_writer
        logging.info('Loading previous restart file...')
        f = gzip.GzipFile(os.path.join(settings.outputDirectory,'restart.pkl'), 'rb')
        species.speciesList, species.speciesCounter, reaction.reactionDict, \
            reactionModel, reactionSystems = cPickle.load(f)
        f.close()
        # Cantera stuff
        reload(ctml_writer) # ensure new empty ctml_writer._species and ._reactions lists
        for reactor in reactionSystems:
            # initialise the ctml_writer thing
            reactor.initializeCantera()
        for spec in reactionModel.core.species:
            # add species to ctml_writer._species list
            spec.toCantera()
        for rxn in reactionModel.core.reactions:
            # add reaction to ctml_writer._reactions list
            rxn.toCantera()
        #print "enter 'c' to continue"; import pdb; pdb.set_trace()
        options.restart = False # have already restarted
    else:
        for spec in coreSpecies:
            reactionModel.enlarge(spec)
#
#    # RMG execution statistics
#    coreSpeciesCount = []
#    coreReactionCount = []
#    edgeSpeciesCount = []
#    edgeReactionCount = []
#    execTime = []
#    restartSize = []
#    memoryUse = []
#
#    # Handle unimolecular (pressure dependent) reaction networks
#    if settings.unimolecularReactionNetworks:
#        reactionModel.updateUnimolecularReactionNetworks()
#        logging.info('')
#
#    # Main RMG loop
#    done = False
#    while not done:
#
#        done = True
#        objectsToEnlarge = []
#        for index, reactionSystem in enumerate(reactionSystems):
#
#            # Conduct simulation
#            logging.info('Conducting simulation of reaction system %s...' % (index+1))
#            t, y, dydt, valid, obj = reactionSystem.simulate(reactionModel)
#
#            # Postprocess results
#            logging.info('')
#            logging.info('Saving simulation results for reaction system %s...' % (index+1))
#            reactionSystem.postprocess(reactionModel, t, y, dydt, str(index+1))
#
#            # If simulation is invalid, note which species should be added to
#            # the core
#            if not valid:
#                objectsToEnlarge.append(obj)
#                done = False
#
#        if not done:
#            # Enlarge objects identified by the simulation for enlarging
#            # These should be Species or Network objects
#            logging.info('')
#            objectsToEnlarge = list(set(objectsToEnlarge))
#            for object in objectsToEnlarge:
#                reactionModel.enlarge(object)
#
#            # Handle unimolecular (pressure dependent) reaction networks
#            if settings.unimolecularReactionNetworks:
#                reactionModel.updateUnimolecularReactionNetworks()
#                logging.info('')
#
#            # Save the restart file
#            # In order to get all the references preserved, you must pickle all of
#            # the objects in one concerted dump; this also has the added benefits
#            # of using less space and running faster
#            # We also compress the restart file to save space (and lower the
#            # disk read/write time)
#            if settings.saveRestart:
#                frequency, iterations, lastSaveTime, lastSaveIteration = settings.saveRestart
#                settings.saveRestart[-2], settings.saveRestart[-1] = \
#                    saveRestartFile(frequency, iterations, lastSaveTime, lastSaveIteration,
#                    reactionModel, reactionSystems)
#
#            # Update RMG execution statistics
#            logging.info('Updating RMG execution statistics...')
#            coreSpeciesCount.append(len(reactionModel.core.species))
#            coreReactionCount.append(len(reactionModel.core.reactions))
#            edgeSpeciesCount.append(len(reactionModel.edge.species))
#            edgeReactionCount.append(len(reactionModel.edge.reactions))
#            execTime.append(time.time() - settings.initializationTime)
#            from guppy import hpy
#            hp = hpy()
#            memoryUse.append(hp.heap().size / 1.0e6)
#            logging.debug('Execution time: %s s' % (execTime[-1]))
#            logging.debug('Memory used: %s MB' % (memoryUse[-1]))
#            restartSize.append(os.path.getsize(os.path.join(settings.outputDirectory,'restart.pkl')) / 1.0e6)
#            saveExecutionStatistics(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)
#            generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)
#
#        logging.info('')
#
#        # Consider stopping gracefully if the next iteration might take us
#        # past the wall time
#        if settings.wallTime > 0 and len(execTime) > 1:
#            t = execTime[-1]
#            dt = execTime[-1] - execTime[-2]
#            if t + 2 * dt > settings.wallTime:
#                logging.info('MODEL GENERATION TERMINATED')
#                logging.info('')
#                logging.info('There is not enough time to complete the next iteration before the wall time is reached.')
#                logging.info('The output model may be incomplete.')
#                logging.info('')
#                logging.info('The current model core has %s species and %s reactions' % (len(reactionModel.core.species), len(reactionModel.core.reactions)))
#                logging.info('The current model edge has %s species and %s reactions' % (len(reactionModel.edge.species), len(reactionModel.edge.reactions)))
#                io.writeOutputFile(os.path.join(settings.outputDirectory,'output.xml'), reactionModel, reactionSystems)
#                return
#
#    # Write output file
#    logging.info('MODEL GENERATION COMPLETED')
#    logging.info('')
#    logging.info('The final model core has %s species and %s reactions' % (len(reactionModel.core.species), len(reactionModel.core.reactions)))
#    logging.info('The final model edge has %s species and %s reactions' % (len(reactionModel.edge.species), len(reactionModel.edge.reactions)))
#    io.writeOutputFile(os.path.join(settings.outputDirectory,'output.xml'), reactionModel, reactionSystems)

    # Log end timestamp
    logging.info('')
    logging.info('RMG execution terminated at ' + time.asctime())

################################################################################

def processStats(stats_file, log_file):
        out_stream = tee(sys.stdout,open(log_file,'a')) # print to screen AND append to RMG.log
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

################################################################################

if __name__ == '__main__':

    # Initialize the memory profiler
    # It works best if we do this as the very first thing
    from guppy import hpy
    hp = hpy()
    hp.heap()

    # Parse the command-line arguments (requires the argparse module)
    args = parseCommandLineArguments()

    # For output and scratch directories, if they are empty strings, set them
    # to match the input file location
    import os.path
    inputDirectory = os.path.abspath(os.path.dirname(args.file[0]))
    if args.output_directory == '':
        args.output_directory = inputDirectory
    if args.scratch_directory == '':
        args.scratch_directory = inputDirectory

    # Initialize the logging system
    level = logging.INFO
    if args.debug: level = 0
    elif args.verbose: level = logging.DEBUG
    elif args.quiet: level = logging.WARNING
    initializeLog(level, os.path.join(args.output_directory,'RMG.log'))

    if args.postprocess:
        print "Postprocessing the profiler statistics (will be appended to RMG.log)"
        args.profile = True

    if args.profile:
        import cProfile, sys, pstats, os
        global_vars = {}
        local_vars = {'args': args, 'options':options, 'rmg':rmg}
        command = """execute(args)"""
        stats_file = os.path.join(options.outputDirectory,'RMG.profile')
        print("Running under cProfile")
        if not args.postprocess:
            # actually run the program!
            cProfile.runctx(command, global_vars, local_vars, stats_file)
        # postprocess the stats
        log_file = os.path.join(options.outputDirectory,'RMG.log')
        processStats(stats_file, log_file)

    else:
        execute(args)

