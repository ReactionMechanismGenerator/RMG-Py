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
This script enables the automatic merging of two or more Chemkin files (and
associated species dictionaries) into a single unified Chemkin file. Simply
pass the paths of the Chemkin files and species dictionaries on the 
command-line, e.g.

    $ python mergeModels.py --model1 /path/to/chem1.inp /path/to/species_dictionary1.txt --model2 /path/to/chem2.inp /path/to/species_dictionary2.txt

The resulting merged files are placed in ``chem.inp`` and
``species_dictionary.txt`` in the execution directory by default.
The output directory can be changed by passing a new path via the
``--output-directory`` argument.
"""

import argparse
import logging
import os

import rmgpy.tools.merge_models as merge_models
from rmgpy.rmg.main import RMG, initializeLog

################################################################################


def parse_command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m1', '--model1', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the first model to merge')
    parser.add_argument('-m2', '--model2', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the second model to merge')
    parser.add_argument('-m3', '--model3', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the third model to merge')
    parser.add_argument('-m4', '--model4', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the fourth model to merge')
    parser.add_argument('-m5', '--model5', metavar='FILE', type=str, nargs='+',
                        help='the Chemkin files and species dictionaries of the fifth model to merge')
    parser.add_argument('-o', '--output-directory', metavar='DIR', nargs=1,
                        help='output directory for final merged mechanism')
    parser.add_argument('-r', '--react', metavar='INPUT', nargs=1,
                        help='provide RMG input file to generate bimolecular cross reactions between models')
    parser.add_argument('-a', '--add-all', action='store_true',
                        help='add all new species generated from cross reactions')

    # RMG options
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')

    parser.add_argument('-n', '--maxproc', type=int, nargs=1, default=1,
                        help='max number of processes used during reaction generation')
    parser.add_argument('-t', '--walltime', type=str, nargs=1, default='00:00:00:00',
                        metavar='DD:HH:MM:SS', help='set the maximum execution time')

    return parser.parse_args()


def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.
    """
    args = parse_command_line_arguments()

    input_files = []
    for i, model in enumerate([args.model1, args.model2, args.model3, args.model4, args.model5]):
        if model is None:
            continue
        elif len(model) == 2:
            input_files.append((model[0], model[1], None))
        elif len(model) == 3:
            input_files.append((model[0], model[1], model[2]))
        else:
            raise ValueError('Unexpected number of arguments for model {0}. Arguments should include '
                             'Chemkin file, species dictionary, and an optional transport file.'.format(i + 1))

    if args.output_directory is not None:
        output_directory = args.output_directory[0]
    else:
        output_directory = os.getcwd()

    level = logging.INFO
    if args.debug:
        level = 0
    elif args.verbose:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING

    initializeLog(level, os.path.join(output_directory, 'merge.log'))

    if args.react is not None:
        input_file = args.react[0]

        rmg = RMG(inputFile=input_file, outputDirectory=output_directory)

        if args.walltime != '00:00:00:00':
            args.walltime = args.walltime[0]

        if args.maxproc != 1:
            args.maxproc = args.maxproc[0]

        kwargs = {
            'walltime': args.walltime,
            'maxproc': args.maxproc,
        }

        rmg.initialize(**kwargs)
    else:
        rmg = None

    merge_models.execute(input_files, output_directory=output_directory, rmg=rmg, add_all=args.add_all)


if __name__ == '__main__':
    main()
