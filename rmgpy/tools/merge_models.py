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
This module contains a function for merging two or more Chemkin files (and
associated species dictionaries) into a single unified Chemkin file.

Use the RMG-Py/scripts/mergeModels.py script for command line use.
"""

from __future__ import print_function

import os

from rmgpy.chemkin import loadChemkinFile, saveChemkinFile, saveSpeciesDictionary, saveTransportFile
from rmgpy.rmg.model import ReactionModel

################################################################################


def execute(input_files, output_directory=None, rmg=None, add_all=False):

    if output_directory is None:
        output_directory = os.getcwd()

    chem_output = os.path.join(output_directory, 'chem_merged.inp')
    dict_output = os.path.join(output_directory, 'species_dictionary_merged.txt')
    tran_output = os.path.join(output_directory, 'tran_merged.dat')
    
    # Load the models to merge
    transport = False
    models = []
    for chem_path, dict_path, tran_path in input_files:
        if tran_path is not None:
            transport = True
        print('Loading model #{0:d}...'.format(len(models) + 1))
        model = ReactionModel()
        model.species, model.reactions = loadChemkinFile(chem_path, dict_path, transportPath=tran_path)
        models.append(model)

    final_model = ReactionModel()
    for i, model in enumerate(models):        
        print('Ignoring common species and reactions from model #{0:d}...'.format(i + 1))
        nspec0 = len(final_model.species)
        nrxn0 = len(final_model.reactions)
        final_model = final_model.merge(model)
        nspec = len(final_model.species)
        nrxn = len(final_model.reactions)

        if len(model.species) == 0:
            print('Added no species from model #{0:d} because it does not contain any species.'.format(i + 1))
        else:
            print('Added {1:d} out of {2:d} ({3:.1f}%) unique species from model '
                  '#{0:d}.'.format(i + 1, nspec - nspec0, len(model.species),
                                   (nspec - nspec0) * 100. / len(model.species)))
        if len(model.reactions) == 0:
            print('Added no reactions from model #{0:d} because it does not contain any reactions.'.format(i + 1))
        else:
            print('Added {1:d} out of {2:d} ({3:.1f}%) unique reactions from model '
                  '#{0:d}.'.format(i + 1, nrxn - nrxn0, len(model.reactions),
                                   (nrxn - nrxn0) * 100. / len(model.reactions)))

    print('The merged model has {0:d} species and {1:d} reactions'.format(len(final_model.species),
                                                                          len(final_model.reactions)))

    # Save the merged model to disk
    saveChemkinFile(chem_output, final_model.species, final_model.reactions)
    saveSpeciesDictionary(dict_output, final_model.species)
    if transport:
        saveTransportFile(tran_output, final_model.species)

    print('Merged Chemkin file saved to {0}'.format(chem_output))
    print('Merged species dictionary saved to {0}'.format(dict_output))
    if transport:
        print('Merged transport file saved to {0}'.format(tran_output))
