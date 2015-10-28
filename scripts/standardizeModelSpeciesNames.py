#!/usr/bin/env python
# encoding: utf-8

"""
This script enables the automatic renaming of species names of of two or more Chemkin files (and
associated species dictionaries) so that they use consistent, matching names.  Simply
pass the paths of the Chemkin files and species dictionaries on the 
command-line, e.g.

    $ python standardizeModelSpeciesNames.py --model1 /path/to/chem1.inp /path/to/species_dictionary1.txt --model2 /path/to/chem2.inp /path/to/species_dictionary2.txt

The resulting files are saved as ``chem1.inp`` and
``species_dictionary1.txt``, ``chem2.inp``, ``species_dictionary2.txt`` and so forth in the execution directory.
"""

import os.path
import argparse

from rmgpy.chemkin import loadChemkinFile, saveChemkinFile, saveSpeciesDictionary, saveTransportFile
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
        if model is None: continue
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
        print 'Loading model #{0:d}...'.format(len(models)+1)
        model = ReactionModel()
        model.species, model.reactions = loadChemkinFile(chemkin, speciesPath, transportPath=transportPath)
        models.append(model)

    allSpecies = []
    speciesIndices = [[] for i in range(len(models))]
    for i, model in enumerate(models):       
        speciesIndices[i] = []
        for j, species in enumerate(model.species):
            for index, species0 in enumerate(allSpecies):
                if species0.isIsomorphic(species):
                    speciesIndices[i].append(index)
                    break; 
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
        saveChemkinFile('chem{0}.inp'.format(i+1), model.species, model.reactions)
        saveSpeciesDictionary('species_dictionary{0}.txt'.format(i+1), model.species)
        
    print 'Saving of new models with consistent names is complete!'