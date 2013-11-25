#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script runs stand-alone thermo estimation using RMG for a list of species in a
thermo input file.  It generates an output.txt file containing the chemkin format
thermochemistry as well as a ThermoLibrary file containing the enthalpy, entropy, and
heat capacity data in RMG-database format.
"""

import os.path
import logging
from rmgpy.rmg.main import RMG, initializeLog, processProfileStats, makeProfileGraph
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.chemkin import writeThermoEntry
from rmgpy.rmg.model import makeThermoForSpecies
from scoop import futures,shared
import resource # to see memory usage
################################################################################
def chunks(l, n):
    """ 
    Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]
        
def runThermoEstimator(inputFile,chunkSize):
    """
    Estimate thermo for a list of species using RMG and the settings chosen inside a thermo input file.
    """
    logging.debug("Maximum memory usage:{0} MBs.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))
    rmg = RMG()

    logging.debug("RMG object created...")
    logging.debug("Maximum memory usage:{0} MBs.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))
    rmg.loadThermoInput(inputFile)
    logging.debug("Input file loaded...")
    logging.debug("Maximum memory usage:{0} MBs.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))
    
    # initialize and load the database as well as any QM settings
    rmg.loadThermoDatabase()
    logging.debug("Thermo database loaded...")
    logging.debug("Maximum memory usage:{0} MBs.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))
    if rmg.quantumMechanics:
        rmg.quantumMechanics.initialize()
        logging.debug("QM module initialized...")
        logging.debug("Maximum memory usage:{0} MBs.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))
    
    # Generate the thermo for all the species and write them to chemkin format as well as
    # ThermoLibrary format with values for H, S, and Cp's.
    output = open(os.path.join(rmg.outputDirectory, 'output.txt'),'wb')
    listOfSpecies=rmg.initialSpecies
    logging.debug("Initial species loaded...")
    logging.debug("Maximum memory usage:{0} MBs.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))

    chunkIndex=0
    shared.setConst(qmValue=rmg.reactionModel.quantumMechanics)
    for chunk in list(chunks(listOfSpecies,chunkSize)):
        # There will be no stdout from workers except the main one.
        outputList = futures.map(makeThermoForSpecies, chunk)
        if chunkIndex == 0: libraryName = 'ThermoLibrary'
        else: libraryName = 'ThermoLibrary'+ str(chunkIndex)
        library = ThermoLibrary(name=libraryName)
        for species, thermo in zip(chunk, outputList):
            logging.debug("Species {0}".format(species.label))
            species.thermo = thermo   
            library.loadEntry(
                index = len(library.entries) + 1,
                label = species.label,
                molecule = species.molecule[0].toAdjacencyList(),
                thermo = species.thermo.toThermoData(),
                shortDesc = species.thermo.comment,
            )
            output.write(writeThermoEntry(species))
            output.write('\n')
        logging.debug("Thermo library created...")
        logging.debug("Maximum memory usage:{0} MBs.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))
        library.save(os.path.join(rmg.outputDirectory, libraryName + '.py'))
        del library
        chunkIndex += 1
    output.close()
    logging.debug("runThermoEstimator is done.")
    


################################################################################

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description=
     """
     thermoEstimator.py generates thermochemical parameters based on Benson group additivity
     or quantum mechanical calculations. \n
     Generates three output files.
     RMG.log: Contains information about the process.
     output.txt: Contains string representations of the NASA model for each species, readable by Chemkin.
     ThermoLibrary.py: Thermo library that can be used in RMG simulations. Can be uploaded to RMG-database.
     """)
    parser.add_argument('input', metavar='FILE', type=str, nargs=1,
        help='Thermo input file')
    parser.add_argument('CHUNKSIZE', type=int, default=10000,nargs='?', help='''chunk size that determines number of species passed to 
         workers at once, should be larger than the number of processors. (default value is 10000)''')
    parser.add_argument('-p', '--profile', action='store_true', help='run under cProfile to gather profiling statistics, and postprocess them if job completes')
    parser.add_argument('-P', '--postprocess', action='store_true', help='postprocess profiling statistics from previous [failed] run; does not run the simulation')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    
    
    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    inputDirectory = os.path.abspath(os.path.dirname(args.input[0]))
    chunkSize = args.CHUNKSIZE
    if args.postprocess:
        print "Postprocessing the profiler statistics (will be appended to thermo.log)"
        print  "Use `dot -Tpdf RMG.profile.dot -o RMG.profile.pdf`"
        args.profile = True
    
    if args.profile:
        import cProfile, sys, pstats, os
        global_vars = {}
        local_vars = {'inputFile': inputFile,'runThermoEstimator':runThermoEstimator}
        command = """runThermoEstimator(inputFile)"""
        stats_file = 'RMG.profile'
        print("Running under cProfile")
        if not args.postprocess:
        # actually run the program!
            cProfile.runctx(command, global_vars, local_vars, stats_file)
        # postprocess the stats
        log_file = os.path.join(inputDirectory,'RMG.log')
        processProfileStats(stats_file, log_file)
        makeProfileGraph(stats_file)
        
    else:
        
        if args.debug: level = logging.DEBUG
        elif args.quiet: level = logging.WARNING
        else: level = logging.INFO
        initializeLog(level, 'RMG.log')
        logging.debug("runThermoEstimator starts...")
        runThermoEstimator(inputFile,chunkSize)
