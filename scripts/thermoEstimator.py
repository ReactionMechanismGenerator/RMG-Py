#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script runs stand-alone thermo estimation using RMG for a list of species in a
thermo input file.  It generates an output.txt file containing the chemkin format
thermochemistry as well as a ThermoLibrary file containing the enthalpy, entropy, and
heat capacity data in RMG-database format.
"""

import os.path
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import RMG
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.chemkin import saveChemkinFile, saveSpeciesDictionary
from rmgpy.species import Species
from rmgpy.thermo.thermoengine import submit
                     
################################################################################

def runThermoEstimator(inputFile):
    """
    Estimate thermo for a list of species using RMG and the settings chosen inside a thermo input file.
    """
    
    rmg = RMG()
    rmg.loadThermoInput(inputFile)
    
    rmg.database = RMGDatabase()
    path = os.path.join(settings['database.directory'])

    # forbidden structure loading
    rmg.database.loadThermo(os.path.join(path, 'thermo'), rmg.thermoLibraries, depository=False)
   
    if rmg.solvent:
        rmg.database.loadSolvation(os.path.join(path, 'solvation'))
        rmg.solvent.solventData = rmg.database.solvation.getSolventData(rmg.solvent.solventName)

    for species in rmg.initialSpecies:
        # Because thermoEstimator does not used rmgpy.rmg.main where solvent is set as global variable,
        # the solvent data must be passed on directly to thermoEngine. 
        submit(species, rmg.solvent)

    # library = ThermoLibrary(name='Thermo Estimation Library')
    # for spc in rmg.initialSpecies:
    #     library.loadEntry(
    #         index = len(library.entries) + 1,
    #         label = species.label,
    #         molecule = species.molecule[0].toAdjacencyList(),
    #         thermo = species.getThermoData().toThermoData(),
    #         shortDesc = species.getThermoData().comment,
    #     )
    # library.save(os.path.join(rmg.outputDirectory,'ThermoLibrary.py'))
    

    # Save the thermo data to chemkin format output files and dictionary, with no reactions    
    saveChemkinFile(os.path.join(rmg.outputDirectory, 'chem_annotated.inp'), species=rmg.initialSpecies, reactions=[])
    saveSpeciesDictionary(os.path.join(rmg.outputDirectory, 'species_dictionary.txt'), species=rmg.initialSpecies)

################################################################################

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='Thermo input file')
    args = parser.parse_args()
    
    inputFile = os.path.abspath(args.input[0])
    
    runThermoEstimator(inputFile)