#!/bin/bash

# Normal Output
python ../../../scripts/diffModels.py Models/DFT_chem.inp Models/DFT_species_dictionary.txt Models/Primary_chem.inp Models/Primary_species_dictionary.txt

# Do not show identical species or reactions
#python ../../../scripts/diffModels.py Models/DFT_chem.inp Models/DFT_species_dictionary.txt Models/Primary_chem.inp Models/Primary_species_dictionary.txt --diffOnly

# Only show species or reactions that are present in both models but have different values
#python ../../../scripts/diffModels.py Models/DFT_chem.inp Models/DFT_species_dictionary.txt Models/Primary_chem.inp Models/Primary_species_dictionary.txt --commonDiffOnly
