.. _standardizeModelSpeciesNames:

*******************************
Standardize Model Species Names
*******************************

This script enables the automatic renaming of species names of of two or more Chemkin files (and
associated species dictionaries) so that they use consistent, matching names.  Simply
pass the paths of the Chemkin files and species dictionaries on the 
command-line, e.g.::

    python standardizeModelSpeciesNames.py --model1 /path/to/chem1.inp /path/to/species_dictionary1.txt --model2 /path/to/chem2.inp /path/to/species_dictionary2.txt

The resulting files are saved as ``chem1.inp`` and
``species_dictionary1.txt``, ``chem2.inp``, ``species_dictionary2.txt`` and so forth (depending on
how many models you want to standardize) and will be saved in the execution directory.
