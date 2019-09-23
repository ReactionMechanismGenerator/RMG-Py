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
This script enables the automatic renaming of species names of of two or more Chemkin files (and
associated species dictionaries) so that they use consistent, matching names.  Simply
pass the paths of the Chemkin files and species dictionaries on the 
command-line, e.g.

    $ python standardizeModelSpeciesNames.py --model1 /path/to/chem1.inp /path/to/species_dictionary1.txt --model2 /path/to/chem2.inp /path/to/species_dictionary2.txt

The resulting files are saved as ``chem1.inp`` and
``species_dictionary1.txt``, ``chem2.inp``, ``species_dictionary2.txt`` and so forth in the execution directory.
"""
from __future__ import print_function

import argparse

from rmgpy.chemkin import load_chemkin_file, save_chemkin_file, save_species_dictionary
from rmgpy.rmg.model import ReactionModel

################################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--model1', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the first model')
    parser.add_argument('--model2', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the second model')
    parser.add_argument('--model3', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the third model')
    parser.add_argument('--model4', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the fourth model')
    parser.add_argument('--model5', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the fifth model')

    args = parser.parse_args()

    transport = False
    inputModelFiles = []
    for model in [args.model1, args.model2, args.model3, args.model4, args.model5]:
        if model is None:
            continue
        if len(model) == 2:
            inputModelFiles.append((model[0], model[1], None))
        elif len(model) == 3:
            transport = True
            inputModelFiles.append((model[0], model[1], model[2]))
        else:
            raise Exception

    outputChemkinFile = 'chem.inp'
    outputSpeciesDictionary = 'species_dictionary.txt'
    outputTransportFile = 'tran.dat' if transport else None

    # Load the models to merge
    models = []
    for chemkin, speciesPath, transportPath in inputModelFiles:
        print('Loading model #{0:d}...'.format(len(models) + 1))
        model = ReactionModel()
        model.species, model.reactions = load_chemkin_file(chemkin, speciesPath, transportPath=transportPath)
        models.append(model)

    allSpecies = []
    speciesIndices = [[] for i in range(len(models))]
    for i, model in enumerate(models):
        speciesIndices[i] = []
        for j, species in enumerate(model.species):
            for index, species0 in enumerate(allSpecies):
                if species0.is_isomorphic(species):
                    speciesIndices[i].append(index)
                    break
            else:
                allSpecies.append(species)
                speciesIndices[i].append(allSpecies.index(species))
    # Reassign species names and labels according to the list of all species in all models
    # We must retain the original thermochemistry
    for i, model in enumerate(models):
        for j, species in enumerate(model.species):
            index = speciesIndices[i][j]
            species.label = allSpecies[index].label
            species.index = allSpecies[index].index

        # Resave the models    
        save_chemkin_file('chem{0}.inp'.format(i + 1), model.species, model.reactions)
        save_species_dictionary('species_dictionary{0}.txt'.format(i + 1), model.species)

    print('Saving of new models with consistent names is complete!')
