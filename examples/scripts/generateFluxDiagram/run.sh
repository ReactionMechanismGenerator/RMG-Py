#!/bin/bash

# Run generateFluxDiagram.py on a RMG input file, a Chemkin input file, and a RMG species dictionary file
python ../../../scripts/generateFluxDiagram.py input.py chem_annotated.inp species_dictionary.txt
