#!/bin/bash

# standardize species names for the two models
python ../../../scripts/standardizeModelSpeciesNames.py --model1 Models/minimal_chem.inp Models/minimal_species_dictionary.txt --model2 Models/superminimal_chem.inp Models/superminimal_species_dictionary.txt
