#!/usr/bin/env python
# encoding: utf-8

"""
This script enables the automatic merging of two or more Chemkin files (and
associated species dictionaries) into a single unified Chemkin file. Simply
pass the paths of the Chemkin files and species dictionaries on the 
command-line, e.g.

    $ python mergeModels.py /path/to/chem1.inp /path/to/species_dictionary1.txt /path/to/chem2.inp /path/to/species_dictionary2.txt

The resulting merged files are placed in ``chem.inp`` and
``species_dictionary.txt`` in the execution directory.
"""

import os.path
import argparse

from rmgpy.chemkin import loadChemkinFile, saveChemkinFile, saveSpeciesDictionary
from rmgpy.reaction import ReactionModel

################################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('files', metavar='FILE', type=str, nargs='+',
        help='the Chemkin files and species dictionaries of each model to merge')
    
    args = parser.parse_args()
    
    outputChemkinFile = 'chem.inp'
    outputSpeciesDictionary = 'species_dictionary.txt'
    
    assert len(args.files) % 2 == 0
    
    # Load the models to merge
    models = []
    for chemkin, speciesDict in zip(args.files[0::2], args.files[1::2]):
        print 'Loading model #{0:d}...'.format(len(models)+1)
        model = ReactionModel()
        model.species, model.reactions = loadChemkinFile(chemkin, speciesDict)
        models.append(model)
    


    finalModel = ReactionModel()
    for i, model in enumerate(models):        
        print 'Ignoring common species and reactions from model #{0:d}...'.format(i+1)
        Nspec0 = len(finalModel.species)
        Nrxn0 = len(finalModel.reactions)
        finalModel = finalModel.merge(model)
        Nspec = len(finalModel.species)
        Nrxn = len(finalModel.reactions)
        print 'Added {1:d} out of {2:d} ({3:.1f}%) unique species from model #{0:d}.'.format(i+1, Nspec - Nspec0, len(model.species), (Nspec - Nspec0) * 100. / len(model.species))
        print 'Added {1:d} out of {2:d} ({3:.1f}%) unique reactions from model #{0:d}.'.format(i+1, Nrxn - Nrxn0, len(model.reactions), (Nrxn - Nrxn0) * 100. / len(model.reactions))
    
    print 'The merged model has {0:d} species and {1:d} reactions'.format(len(finalModel.species), len(finalModel.reactions))
        
    # Save the merged model to disk
    saveChemkinFile(outputChemkinFile, finalModel.species, finalModel.reactions)
    saveSpeciesDictionary(outputSpeciesDictionary, finalModel.species)

    print 'Merged Chemkin file saved to {0}'.format(outputChemkinFile)
    print 'Merged species dictionary saved to {0}'.format(outputSpeciesDictionary)
