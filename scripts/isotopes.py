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
This script accepts one input file (e.g. input.py) with the RMG-Py model to generate,
optional parameters `--original [folder of original rmg model] ` can allow using
a starting RMG model. A special path can be added  with the argument `--output` for the
path to output the final files.
"""

import argparse
import logging
import os
import os.path

from rmgpy.exceptions import InputError
from rmgpy.rmg.main import initialize_log
from rmgpy.tools.isotopes import run


################################################################################


def parse_command_line_arguments():
    """
    Parse the command-line arguments being passed to RMG-Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='RMG input file')
    parser.add_argument('--output', type=str, nargs=1, default='', help='Output folder')
    parser.add_argument('--original', type=str, nargs=1, default='',
                        help='Location of the isotopeless mechanism')
    parser.add_argument('--maximumIsotopicAtoms', type=int, nargs=1, default=[1000000],
                        help='The maxuminum number of isotopes you allow in a specific molecule')
    parser.add_argument('--useOriginalReactions', action='store_true', default=False,
                        help='Use reactions from the original rmgpy generated chem_annotated.inp file')
    parser.add_argument('--kineticIsotopeEffect', type=str, nargs=1, default='',
                        help='Type of kinetic isotope effects to use, currently only "simple" supported.')
    args = parser.parse_args()

    return args


def main():
    args = parse_command_line_arguments()
    if args.useOriginalReactions and not args.original:
        raise InputError('Cannot use original reactions without a previously run RMG job')
    maximum_isotopic_atoms = args.maximumIsotopicAtoms[0]
    use_original_reactions = args.useOriginalReactions
    input_file = args.input
    outputdir = os.path.abspath(args.output[0]) if args.output else os.path.abspath('.')
    original = os.path.abspath(args.original[0]) if args.original else None
    kie = args.kineticIsotopeEffect[0] if args.kineticIsotopeEffect else None
    supported_kie_methods = ['simple']
    if kie not in supported_kie_methods and kie is not None:
        raise InputError('The kie input, {0}, is not one of the currently supported methods, {1}'.format(kie, supported_kie_methods))
    initialize_log(logging.INFO, os.path.join(os.getcwd(), 'RMG.log'))
    run(input_file, outputdir, original=original,
        maximum_isotopic_atoms=maximum_isotopic_atoms,
        use_original_reactions=use_original_reactions,
        kinetic_isotope_effect=kie)


if __name__ == '__main__':
    main()
