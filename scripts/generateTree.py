#!/usr/bin/env python3

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
This script cleans and generates new trees for the specified family running in parallel on the specified number of processors
`python generateTree.py familyName nprocs`
ex:  `python generateTree.py intra_H_migration 6`
Note that 6 is the maximum number of processors used currently by this script
This script will work only for the certain reaction families that have been fixed to run this script.
Running this script on other reaction families will throw error.
The documentation for this script will be ready in the future.
"""

import argparse
import logging
import os
import os.path

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import initialize_log


################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('name', metavar='NAME', type=str, nargs=1,
                        help='Family Name')

    parser.add_argument('nprocs', metavar='NPROCS', type=int, nargs=1,
                        help='Number of Processors for Parallelization')
    args = parser.parse_args()

    name = args.name[0]
    nprocs = args.nprocs[0]

    return name, nprocs


def main():
    initialize_log(logging.INFO, 'treegen.log')
    dbdir = settings['database.directory']
    family_name, nprocs = parse_arguments()
    database = RMGDatabase()
    database.load(
        path=dbdir,
        thermo_libraries=['Klippenstein_Glarborg2016', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo',
                         'primaryThermoLibrary', 'primaryNS', 'NitrogenCurran', 'NOx2018', 'FFCM1(-)',
                         'SulfurLibrary', 'SulfurGlarborgH2S'],
        transport_libraries=[],
        reaction_libraries=[],
        seed_mechanisms=[],
        kinetics_families=[family_name],
        kinetics_depositories=['training'],
        # frequenciesLibraries = self.statmechLibraries,
        depository=False,  # Don't bother loading the depository information, as we don't use it
    )
    family = database.kinetics.families[family_name]
    family.clean_tree()
    family.generate_tree(thermo_database=database.thermo, nprocs=min(4, nprocs))
    family.check_tree()
    family.regularize()
    template_rxn_map = family.get_reaction_matches(thermo_database=database.thermo, remove_degeneracy=True,
                                                   get_reverse=True, fix_labels=True)
    family.make_bm_rules_from_template_rxn_map(template_rxn_map, nprocs=min(6, nprocs))
    family.check_tree()
    family.save(os.path.join(dbdir, 'kinetics', 'families', family_name))


################################################################################

if __name__ == '__main__':
    main()
