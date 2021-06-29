#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This script can be used to compare two RMG-generated kinetics models. To use,
pass the chem.inp and species_dictionary.txt files to the script. The syntax
is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 CHEMKIN2 SPECIESDICT2

Optionally, you may use the --thermo1 and/or --thermo2 flags to add separate
thermo chemkin files.

The optional --web flag is used for running this script through the RMG-website

With all the above options the syntax is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 --thermo1 THERMO1 CHEMKIN2 SPECIESDICT2 --thermo2 THERMO2 --web

Further option flags:
======================= ====================================================================================
Flag                    Description
======================= ====================================================================================
--diffOnly              Only show species and reactions which are unique or have different values
--commonDiffOnly        Only show species and reactions present in BOTH models which have different values
======================= ==================================================================================== 
"""
import rmgpy.tools.diffmodels as diff_models


################################################################################

def main():
    diff_models.main()


if __name__ == '__main__':
    main()
