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
This script cleans and generates new trees for the specified family running in parallel on the specified number of processors
`python generateTree.py familyName nprocs`
ex:  `python generateTree.py intra_H_migration 6`
Note that 6 is the maximum number of processors used currently by this script
"""

import os
import os.path
import argparse
import logging
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import initializeLog

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
    initializeLog(logging.INFO,'treegen.log')
    dbdir = settings['database.directory']
    familyName, nprocs = parse_arguments()
    database = RMGDatabase()
    database.load(
        path=dbdir,
        thermoLibraries=['Klippenstein_Glarborg2016', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo',
                           'primaryThermoLibrary', 'primaryNS', 'NitrogenCurran', 'NOx2018', 'FFCM1(-)',
                           'SulfurLibrary', 'SulfurGlarborgH2S'],
        transportLibraries=[],
        reactionLibraries=[],
        seedMechanisms=[],
        kineticsFamilies=[familyName],
        kineticsDepositories=['training'],
        # frequenciesLibraries = self.statmechLibraries,
        depository=False,  # Don't bother loading the depository information, as we don't use it
    )
    family = database.kinetics.families[familyName]
    family.cleanTree(database.thermo)
    family.generateTree(thermoDatabase=database.thermo, nprocs=min(4, nprocs))
    family.checkTree()
    family.regularize()
    templateRxnMap = family.getReactionMatches(thermoDatabase=database.thermo, removeDegeneracy=True, getReverse=True, fixLabels=True)
    family.makeBMRulesFromTemplateRxnMap(templateRxnMap, nprocs=min(6, nprocs))
    family.checkTree()
    family.save(os.path.join(dbdir,'kinetics','families',familyName))
################################################################################

if __name__ == '__main__':
    main()
