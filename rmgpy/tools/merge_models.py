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
This script enables the automatic merging of two or more Chemkin files (and
associated species dictionaries) into a single unified Chemkin file. Simply
pass the paths of the Chemkin files and species dictionaries on the 
command-line, e.g.

    $ python mergeModels.py --model1 /path/to/chem1.inp /path/to/species_dictionary1.txt --model2 /path/to/chem2.inp /path/to/species_dictionary2.txt

The resulting merged files are placed in ``chem.inp`` and
``species_dictionary.txt`` in the execution directory.
"""

from __future__ import division, print_function

import argparse
import os
import os.path

from rmgpy.chemkin import load_chemkin_file, save_chemkin_file, save_species_dictionary, save_transport_file
from rmgpy.rmg.model import ReactionModel


################################################################################

def parse_command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model1', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the first model to merge')
    parser.add_argument('--model2', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the second model to merge')
    parser.add_argument('--model3', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the third model to merge')
    parser.add_argument('--model4', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the fourth model to merge')
    parser.add_argument('--model5', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the fifth model to merge')

    args = parser.parse_args()
    return args


def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.
    """
    # Parse the command-line arguments (requires the argparse module)
    args = parse_command_line_arguments()

    transport = False
    input_model_files = []
    for model in [args.model1, args.model2, args.model3, args.model4, args.model5]:
        if model is None:
            continue
        if len(model) == 2:
            input_model_files.append((model[0], model[1], None))
        elif len(model) == 3:
            transport = True
            input_model_files.append((model[0], model[1], model[2]))
        else:
            raise Exception

    kwargs = {
        'wd': os.getcwd(),
        'transport': transport,
    }

    execute(input_model_files, **kwargs)


def execute(input_model_files, **kwargs):
    try:
        wd = kwargs['wd']
    except KeyError:
        wd = os.getcwd()

    transport = kwargs['transport']

    output_chemkin_file = os.path.join(wd, 'chem.inp')
    output_species_dictionary = os.path.join(wd, 'species_dictionary.txt')
    output_transport_file = os.path.join(wd, 'tran.dat') if transport else None

    models = get_models_to_merge(input_model_files)

    final_model = combine_models(models)

    # Save the merged model to disk
    save_chemkin_file(output_chemkin_file, final_model.species, final_model.reactions)
    save_species_dictionary(output_species_dictionary, final_model.species)
    if transport:
        save_transport_file(output_transport_file, final_model.species)

    print('Merged Chemkin file saved to {0}'.format(output_chemkin_file))
    print('Merged species dictionary saved to {0}'.format(output_species_dictionary))
    if transport:
        print('Merged transport file saved to {0}'.format(output_transport_file))


def get_models_to_merge(input_model_files):
    """
    Reads input file paths and creates a list of ReactionModel
    """
    models = []
    for chemkin, species_path, transport_path in input_model_files:
        print('Loading model #{0:d}...'.format(len(models) + 1))
        model = ReactionModel()
        model.species, model.reactions = load_chemkin_file(chemkin, species_path, transport_path=transport_path)
        models.append(model)
    return models


def combine_models(models):
    """
    Takes in a list of ReactionModels and and merges them into a single ReactionModel
    Reindexes species with the same label and index
    """
    final_model = ReactionModel()
    for i, model in enumerate(models):
        print('Ignoring common species and reactions from model #{0:d}...'.format(i + 1))
        nspec0 = len(final_model.species)
        nrxn0 = len(final_model.reactions)
        final_model = final_model.merge(model)
        nspec = len(final_model.species)
        nrxn = len(final_model.reactions)
        if len(model.species) > 0:
            print('Added {1:d} out of {2:d} ({3:.1f}%) unique species from model '
                  '#{0:d}.'.format(i + 1, nspec - nspec0, len(model.species),
                                   (nspec - nspec0) * 100. / len(model.species)))
        else:
            print('Added {1:d} out of {2:d} unique species from model '
                  '#{0:d}.'.format(i + 1, nspec - nspec0, len(model.species)))

        if len(model.reactions) > 0:
            print('Added {1:d} out of {2:d} ({3:.1f}%) unique reactions from model '
                  '#{0:d}.'.format(i + 1, nrxn - nrxn0, len(model.reactions),
                                   (nrxn - nrxn0) * 100. / len(model.reactions)))
        else:
            print('Added {1:d} out of {2:d} unique reactions from model '
                  '#{0:d}.'.format(i + 1, nrxn - nrxn0, len(model.reactions)))
    print('The merged model has {0:d} species and {1:d} reactions'
          ''.format(len(final_model.species), len(final_model.reactions)))

    # ensure no species with same name and index
    label_index_dict = {}
    for s in final_model.species:
        if s.label not in label_index_dict:
            label_index_dict[s.label] = [s.index]
        else:
            if s.index in label_index_dict[s.label]:
                # obtained a duplicate
                s.index = max(label_index_dict[s.label]) + 1
                print("Reindexed {0} due to dublicate labels and index".format(s.label))
            label_index_dict[s.label].append(s.index)

    return final_model
