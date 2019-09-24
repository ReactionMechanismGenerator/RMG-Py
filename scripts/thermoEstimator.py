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
This script runs stand-alone thermo estimation using RMG for a list of species in a
thermo input file.  It generates an output.txt file containing the chemkin format
thermochemistry as well as a ThermoLibrary file containing the enthalpy, entropy, and
heat capacity data in RMG-database format.
"""

import os.path

from rmgpy import settings
from rmgpy.chemkin import save_chemkin_file, save_species_dictionary
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import Species
from rmgpy.thermo.thermoengine import submit


################################################################################

def run_thermo_estimator(input_file, library_flag):
    """
    Estimate thermo for a list of species using RMG and the settings chosen inside a thermo input file.
    """

    rmg = RMG()
    rmg.load_thermo_input(input_file)

    rmg.database = RMGDatabase()
    path = os.path.join(settings['database.directory'])

    # forbidden structure loading
    rmg.database.load_thermo(os.path.join(path, 'thermo'), rmg.thermo_libraries, depository=False)

    if rmg.solvent:
        rmg.database.load_solvation(os.path.join(path, 'solvation'))
        Species.solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
        Species.solvent_name = rmg.solvent

    for species in rmg.initial_species:
        submit(species)

    if library_flag:
        library = ThermoLibrary(name='Thermo Estimation Library')
        for species in rmg.initial_species:
            library.load_entry(
                index=len(library.entries) + 1,
                label=species.label,
                molecule=species.molecule[0].to_adjacency_list(),
                thermo=species.get_thermo_data().to_thermo_data(),
                shortDesc=species.get_thermo_data().comment,
            )
        library.save(os.path.join(rmg.output_directory, 'ThermoLibrary.py'))

    # Save the thermo data to chemkin format output files and dictionary, with no reactions    
    save_chemkin_file(os.path.join(rmg.output_directory, 'chem_annotated.inp'), species=rmg.initial_species, reactions=[])
    save_species_dictionary(os.path.join(rmg.output_directory, 'species_dictionary.txt'), species=rmg.initial_species)


################################################################################

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
                        help='Thermo input file')
    parser.add_argument('-l', '--library', action='store_true', help='generate RMG thermo library')

    args = parser.parse_args()

    input_file = os.path.abspath(args.input[0])

    run_thermo_estimator(input_file, args.library)
