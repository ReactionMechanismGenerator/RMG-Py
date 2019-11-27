#!/bin/bash

# merging two models
python ../../../scripts/mergeModels.py --model1 Models/minimal_chem.inp Models/minimal_species_dictionary.txt --model2 Models/superminimal_chem.inp Models/superminimal_species_dictionary.txt

# merging three models
#python ../../../scripts/mergeModels.py --model1 Models/minimal_chem.inp Models/minimal_species_dictionary.txt --model2 Models/superminimal_chem.inp Models/superminimal_species_dictionary.txt --model3 Models/minimal_ml_chem.inp Models/minimal_ml_species_dictionary.txt 
