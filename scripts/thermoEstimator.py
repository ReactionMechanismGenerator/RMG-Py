#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script runs stand-alone thermo estimation using RMG for a list of species in a
thermo input file.  It generates an output.txt file containing the chemkin format
thermochemistry as well as a ThermoLibrary file containing the enthalpy, entropy, and
heat capacity data in RMG-database format.
"""

import os.path
from rmgpy.rmg.main import RMG
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.chemkin import writeThermoEntry
from rmgpy.rmg.model import Species

################################################################################

def runThermoEstimator(inputFile):
    """
    Estimate thermo for a list of species using RMG and the settings chosen inside a thermo input file.
    """
    
    rmg = RMG()
    rmg.loadThermoInput(inputFile)
    
    # initialize and load the database as well as any QM settings
    rmg.loadDatabase()
    if rmg.quantumMechanics:
        rmg.quantumMechanics.initialize()
   
    if rmg.solvent:
        Species.solventData = rmg.database.solvation.getSolventData(rmg.solvent)
        Species.solventName = rmg.solvent
        
    # Generate the thermo for all the species and write them to chemkin format as well as
    # ThermoLibrary format with values for H, S, and Cp's.
    output = open(os.path.join(rmg.outputDirectory, 'output.txt'),'wb')
    library = ThermoLibrary(name='Thermo Estimation Library')
    for species in rmg.initialSpecies:
        species.generateThermoData(rmg.database, quantumMechanics=rmg.reactionModel.quantumMechanics)

        library.loadEntry(
            index = len(library.entries) + 1,
            label = species.label,
            molecule = species.molecule[0].toAdjacencyList(),
            thermo = species.thermo.toThermoData(),
            shortDesc = species.thermo.comment,
        )
        output.write(writeThermoEntry(species))
        output.write('\n')
    
    output.close()
    library.save(os.path.join(rmg.outputDirectory,'ThermoLibrary.py'))


################################################################################

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='Thermo input file')
    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    
    runThermoEstimator(inputFile)