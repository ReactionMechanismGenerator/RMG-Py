#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script can be used to compare two RMG-generated kinetics models. To use,
pass the chem.inp and species_dictionary.txt files to the script. The syntax
is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 CHEMKIN2 SPECIESDICT2

Optionally, you may use the --thermo1 and/or --thermo2 flags to add separate
thermo chemkin files.

The optional --web flag is used for running this script through the RMG-website

With all options the syntax is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 --thermo1 THERMO1 CHEMKIN2 SPECIESDICT2 --thermo2 THERMO2 --web 
"""
import rmgpy.tools.diff_models as diff_models

################################################################################

def main():
    diff_models.main()

if __name__ == '__main__':
    main()
    
    
