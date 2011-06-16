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
import numpy
import shutil

import rmgpy.rmg.settings as settings
from rmgpy.rmg.input import readInputFile
from rmgpy.rmg.model import Species, PDepNetwork

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

def makeOutputSubdirectory(folder):
    """
    Create a subdirectory `folder` in the output directory. If the folder
    already exists (e.g. from a previous job) its contents are deleted.
    """
    dir = os.path.join(settings.outputDirectory, folder)
    if os.path.exists(dir):
        # The directory already exists, so delete it (and all its content!)
        shutil.rmtree(dir)
    os.mkdir(dir)

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
        data = args.walltime[0].split(':')
        if len(data) == 1:
            settings.wallTime = int(data[-1])
        elif len(data) == 2:
            settings.wallTime = int(data[-1]) + 60 * int(data[-2])
        elif len(data) == 3:
            settings.wallTime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3])
        elif len(data) == 4:
            settings.wallTime = int(data[-1]) + 60 * int(data[-2]) + 3600 * int(data[-3]) + 86400 * int(data[-4])
        else:
            raise ValueError('Invalid format for wall time; should be HH:MM:SS.')

    # Save initialization time
    settings.initializationTime = time.time()

    # Log start timestamp
    logging.info('RMG execution initiated at ' + time.asctime() + '\n')

    # Print out RMG header
    logHeader()

    # See if memory profiling package is available
    try:
        import os
        import psutil
    except ImportError:
        logging.info('Optional package dependency "psutil" not found; memory profiling information will not be saved.')

    # Make output subdirectories
    makeOutputSubdirectory('plot')
    makeOutputSubdirectory('species')
    makeOutputSubdirectory('pdep')
    makeOutputSubdirectory('chemkin')
    
    # Read input file
    reactionModel, coreSpecies, reactionSystems, database, seedMechanisms = readInputFile(inputFile)
    
    # Delete previous HTML file
    from rmgpy.rmg.output import saveOutputHTML
    saveOutputHTML(os.path.join(settings.outputDirectory, 'output.html'), reactionModel)
    
    # Initialize reaction model
    if args.restart:
        reactionModel = loadRestartFile(os.path.join(settings.outputDirectory,'restart.pkl.gz'))
    else:

        # Seed mechanisms: add species and reactions from seed mechanism
        # DON'T generate any more reactions for the seed species at this time
        for seedMechanism in seedMechanisms:
            reactionModel.addSeedMechanismToCore(database, seedMechanism, react=False)
        
        # Add nonreactive species (e.g. bath gases) to core first
        # This is necessary so that the PDep algorithm can identify the bath gas
        for spec in coreSpecies:
            if not spec.reactive:
                reactionModel.enlarge(spec, database)
        # Then add remaining reactive species
        for spec in coreSpecies:
            if spec.reactive:
                spec.generateThermoData(database)
        for spec in coreSpecies:
            if spec.reactive:
                reactionModel.enlarge(spec, database)

    # RMG execution statistics
    coreSpeciesCount = []
    coreReactionCount = []
    edgeSpeciesCount = []
    edgeReactionCount = []
    execTime = []
    restartSize = []
    memoryUse = []

    # Main RMG loop
    done = False
    while not done:

        done = True
        objectsToEnlarge = []
        allTerminated = True
        for index, reactionSystem in enumerate(reactionSystems):

            # Conduct simulation
            logging.info('Conducting simulation of reaction system %s...' % (index+1))
            terminated, obj = reactionSystem.simulate(
                coreSpecies = reactionModel.core.species,
                coreReactions = reactionModel.core.reactions,
                edgeSpecies = reactionModel.edge.species,
                edgeReactions = reactionModel.edge.reactions,
                toleranceKeepInEdge = reactionModel.fluxToleranceKeepInEdge,
                toleranceMoveToCore = reactionModel.fluxToleranceMoveToCore,
                toleranceInterruptSimulation = reactionModel.fluxToleranceInterrupt,
                termination = reactionModel.termination,
                pdepNetworks = reactionModel.unirxnNetworks,
            )
            allTerminated = allTerminated and terminated
            logging.info('')
            
            # If simulation is invalid, note which species should be added to
            # the core
            if obj:
                if isinstance(obj, PDepNetwork):
                    # Determine which species in that network has the highest leak rate
                    # We do this here because we need a temperature and pressure
                    # Store the maximum leak species along with the associated network
                    obj = (obj, obj.getMaximumLeakSpecies(reactionSystem.T, reactionSystem.P))
                objectsToEnlarge.append(obj)
                done = False

        if not done:

            # If we reached our termination conditions, then try to prune
            # species from the edge
            if allTerminated:
                reactionModel.prune(reactionSystems)

            # Enlarge objects identified by the simulation for enlarging
            # These should be Species or Network objects
            logging.info('')
            objectsToEnlarge = list(set(objectsToEnlarge))
            for object in objectsToEnlarge:
                reactionModel.enlarge(object, database)

        # Save the current state of the model core to a pretty HTML file
        logging.info('Saving latest model core to HTML file...')
        saveOutputHTML(os.path.join(settings.outputDirectory, 'output.html'), reactionModel)

        # Save a Chemkin file containing the current model core
        logging.info('Saving latest model core to Chemkin file...')
        this_chemkin_path = os.path.join(settings.outputDirectory, 'chemkin', 'chem%04i.inp' % len(reactionModel.core.species))
        latest_chemkin_path = os.path.join(settings.outputDirectory, 'chemkin','chem.inp')
        reactionModel.saveChemkinFile(this_chemkin_path)
        if os.path.exists(latest_chemkin_path):
            os.unlink(latest_chemkin_path)
        os.link(this_chemkin_path,latest_chemkin_path)

        # Save the restart file if desired
        if settings.saveRestart:
            saveRestartFile(os.path.join(settings.outputDirectory,'restart.pkl.gz'), reactionModel)

        # Update RMG execution statistics
        logging.info('Updating RMG execution statistics...')
        coreSpec, coreReac, edgeSpec, edgeReac = reactionModel.getModelSize()
        coreSpeciesCount.append(coreSpec)
        coreReactionCount.append(coreReac)
        edgeSpeciesCount.append(edgeSpec)
        edgeReactionCount.append(edgeReac)
        execTime.append(time.time() - settings.initializationTime)
        logging.info('    Execution time (HH:MM:SS): %s' % (time.strftime("%H:%M:%S", time.gmtime(execTime[-1]))))
        try:
            import psutil
            process = psutil.Process(os.getpid())
            rss, vms = process.get_memory_info()
            memoryUse.append(rss / 1.0e6)
            logging.info('    Memory used: %.2f MB' % (memoryUse[-1]))
        except ImportError:
            memoryUse.append(0.0)
        if os.path.exists(os.path.join(settings.outputDirectory,'restart.pkl.gz')):
            restartSize.append(os.path.getsize(os.path.join(settings.outputDirectory,'restart.pkl.gz')) / 1.0e6)
            logging.info('    Restart file size: %.2f MB' % (restartSize[-1]))
        else:
            restartSize.append(0.0)
        saveExecutionStatistics(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)
        generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount, edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize)

        logging.info('')

        # Consider stopping gracefully if the next iteration might take us
        # past the wall time
        if settings.wallTime > 0 and len(execTime) > 1:
            t = execTime[-1]
            dt = execTime[-1] - execTime[-2]
            if t + 3 * dt > settings.wallTime:
                logging.info('MODEL GENERATION TERMINATED')
                logging.info('')
                logging.info('There is not enough time to complete the next iteration before the wall time is reached.')
                logging.info('The output model may be incomplete.')
                logging.info('')
                logging.info('The current model core has %s species and %s reactions' % (len(reactionModel.core.species), len(reactionModel.core.reactions)))
                logging.info('The current model edge has %s species and %s reactions' % (len(reactionModel.edge.species), len(reactionModel.edge.reactions)))
                return

    # Write output file
    logging.info('')
    logging.info('MODEL GENERATION COMPLETED')
    logging.info('')
    coreSpec, coreReac, edgeSpec, edgeReac = reactionModel.getModelSize()
    logging.info('The final model core has %s species and %s reactions' % (coreSpec, coreReac))
    logging.info('The final model edge has %s species and %s reactions' % (edgeSpec, edgeReac))
    
    # Log end timestamp
    logging.info('')
    logging.info('RMG execution terminated at ' + time.asctime())

################################################################################

def loadRestartFile(path):
    """
    Load a restart file at `path` on disk.
    """

    import gzip
    import cPickle

    # Unpickle the reaction model from the specified restart file
    logging.info('Loading previous restart file...')
    f = gzip.GzipFile(path, 'rb')
    reactionModel = cPickle.load(f)
    f.close()

    # A few things still point to the species in the input file, so update
    # those to point to the equivalent species loaded from the restart file

    # The termination conversions still point to the old species
    from rmgpy.solver.base import TerminationConversion
    for term in reactionModel.termination:
        if isinstance(term, TerminationConversion):
            term.species, isNew = reactionModel.makeNewSpecies(term.species.molecule[0], term.species.label, term.species.reactive)

    # The initial mole fractions in the reaction systems still point to the old species
    for reactionSystem in reactionSystems:
        initialMoleFractions = {}
        for spec0, moleFrac in reactionSystem.initialMoleFractions.iteritems():
            spec, isNew = reactionModel.makeNewSpecies(spec0.molecule[0], spec0.label, spec0.reactive)
            initialMoleFractions[spec] = moleFrac
        reactionSystem.initialMoleFractions = initialMoleFractions

    # The reactions and reactionDict still point to the old reaction families
    reactionDict = {}
    for family0 in reactionModel.reactionDict:

        # Find the equivalent family in the newly-loaded database
        import rmgpy.data.kinetics
        family = None
        for kineticsDatabase in rmgpy.data.kinetics.kineticsDatabases:
            if isinstance(kineticsDatabase, rmgpy.data.kinetics.KineticsPrimaryDatabase):
                if kineticsDatabase.label == family0.label:
                    family = kineticsDatabase
                    break
            elif isinstance(kineticsDatabase, rmgpy.data.kinetics.KineticsGroupDatabase):
                for label, fam in kineticsDatabase.families.iteritems():
                    if fam.label == family0.label:
                        family = fam
                        break
        if family is None:
            raise Exception("Unable to find matching reaction family for %s" % family0.label)

        # Update each affected reaction to point to that new family
        # Also use that new family in a duplicate reactionDict
        reactionDict[family] = {}
        for reactant1 in reactionModel.reactionDict[family0]:
            reactionDict[family][reactant1] = {}
            for reactant2 in reactionModel.reactionDict[family0][reactant1]:
                reactionDict[family][reactant1][reactant2] = []
                for rxn in reactionModel.reactionDict[family0][reactant1][reactant2]:
                    rxn.family = family
                    rxn.reverse.family = family
                    reactionDict[family][reactant1][reactant2].append(rxn)

    # Return the unpickled reaction model
    return reactionModel

def saveRestartFile(path, reactionModel):
    """
    Save a restart file to `path` on disk containing the contents of the
    provided `reactionModel`.
    """
    import gzip
    import cPickle
    
    # Pickle the reaction model to the specified file
    # We also compress the restart file to save space (and lower the disk read/write time)
    logging.info('Saving restart file...')
    f = gzip.GzipFile(path, 'wb')
    cPickle.dump(reactionModel, f)
    f.close()

################################################################################

def saveExecutionStatistics(execTime, coreSpeciesCount, coreReactionCount, \
    edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize):
    """
    Save the statistics of the RMG job to an Excel spreadsheet for easy viewing
    after the run is complete. The statistics are saved to the file
    `statistics.xls` in the output directory. The ``xlwt`` package is used to
    create the spreadsheet file; if this package is not installed, no file is
    saved.
    """

    # Attempt to import the xlwt package; return if not installed
    try:
        import xlwt
    except ImportError:
        logging.warning('Package xlwt not found. Unable to save execution statistics.')
        return

    # Create workbook and sheet for statistics to be places
    workbook = xlwt.Workbook()
    sheet = workbook.add_sheet('Statistics')

    # First column is execution time
    sheet.write(0,0,'Execution time (s)')
    for i, etime in enumerate(execTime):
        sheet.write(i+1,0,etime)

    # Second column is number of core species
    sheet.write(0,1,'Core species')
    for i, count in enumerate(coreSpeciesCount):
        sheet.write(i+1,1,count)

    # Third column is number of core reactions
    sheet.write(0,2,'Core reactions')
    for i, count in enumerate(coreReactionCount):
        sheet.write(i+1,2,count)

    # Fourth column is number of edge species
    sheet.write(0,3,'Edge species')
    for i, count in enumerate(edgeSpeciesCount):
        sheet.write(i+1,3,count)

    # Fifth column is number of edge reactions
    sheet.write(0,4,'Edge reactions')
    for i, count in enumerate(edgeReactionCount):
        sheet.write(i+1,4,count)

    # Sixth column is memory used
    sheet.write(0,5,'Memory used (MB)')
    for i, memory in enumerate(memoryUse):
        sheet.write(i+1,5,memory)

    # Seventh column is restart file size
    sheet.write(0,6,'Restart file size (MB)')
    for i, memory in enumerate(restartSize):
        sheet.write(i+1,6,memory)

    # Save workbook to file
    fstr = os.path.join(settings.outputDirectory, 'statistics.xls')
    workbook.save(fstr)

################################################################################

def generateExecutionPlots(execTime, coreSpeciesCount, coreReactionCount,
    edgeSpeciesCount, edgeReactionCount, memoryUse, restartSize):
    """
    Generate a number of plots describing the statistics of the RMG job,
    including the reaction model core and edge size and memory use versus
    execution time. These will be placed in the output directory in the plot/
    folder.
    """

    # Only generate plots if that flag is turned on (in input file)
    if not settings.generatePlots:
        return

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.semilogx(execTime, coreSpeciesCount, 'o-b')
    ax1.set_xlabel('Execution time (s)')
    ax1.set_ylabel('Number of core species')
    ax2 = ax1.twinx()
    ax2.semilogx(execTime, coreReactionCount, 'o-r')
    ax2.set_ylabel('Number of core reactions')
    plt.savefig(os.path.join(settings.outputDirectory, 'plot/coreSize.svg'))
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.loglog(execTime, edgeSpeciesCount, 'o-b')
    ax1.set_xlabel('Execution time (s)')
    ax1.set_ylabel('Number of edge species')
    ax2 = ax1.twinx()
    ax2.loglog(execTime, edgeReactionCount, 'o-r')
    ax2.set_ylabel('Number of edge reactions')
    plt.savefig(os.path.join(settings.outputDirectory, 'plot/edgeSize.svg'))
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.semilogx(execTime, memoryUse, 'o-k')
    ax1.semilogx(execTime, restartSize, 'o-g')
    ax1.set_xlabel('Execution time (s)')
    ax1.set_ylabel('Memory (MB)')
    ax1.legend(['RAM', 'Restart file'], loc=2)
    plt.savefig(os.path.join(settings.outputDirectory, 'plot/memoryUse.svg'))
    plt.clf()

################################################################################

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

def processProfileStats(stats_file, log_file):
    import pstats
    out_stream = Tee(sys.stdout,open(log_file,'a')) # print to screen AND append to RMG.log
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
    `dot -Tpdf input.dot -o output.pdf`.
    """
    try:
        from gprof2dot import gprof2dot
    except ImportError:
        logging.warning('Package gprof2dot not found. Unable to create a graph of the profile statistics.')
        # `pip install gprof2dot` if you don't have it.
        return
    m = gprof2dot.Main()
    class Options:
        pass
    m.options = Options()
    m.options.node_thres = 0.5
    m.options.edge_thres = 0.1
    m.options.strip = False
    m.options.wrap = False
    m.theme = m.themes['color']
    parser = gprof2dot.PstatsParser(stats_file)
    m.profile = parser.parse()
    dot_file = stats_file + '.dot'
    m.output = open(dot_file,'wt')
    m.write_graph()
    logging.info("Now try:\n     dot -Tsvg %s -o %s.svg"%(dot_file,dot_file))
    # we could actually try this here using subprocess.Popen() or something
    # wrapped in a try: block.
    
################################################################################

if __name__ == '__main__':

    # Initialize the memory profiler
    # It works best if we do this as the very first thing
    # If the memory profiler package is not installed then carry on
    try:
        from guppy import hpy
        hp = hpy()
        hp.heap()
    except ImportError:
        pass

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
        local_vars = {'args': args, 'execute': execute}
        command = """execute(args)"""
        stats_file = os.path.join(args.output_directory,'RMG.profile')
        print("Running under cProfile")
        if not args.postprocess:
            # actually run the program!
            cProfile.runctx(command, global_vars, local_vars, stats_file)
        # postprocess the stats
        log_file = os.path.join(args.output_directory,'RMG.log')
        processProfileStats(stats_file, log_file)
        makeProfileGraph(stats_file)
        

    else:
        execute(args)

