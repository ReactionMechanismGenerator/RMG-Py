#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
This script is used to generate all the possible reactions involving a given
set of reactants, including pressure-dependent effects if desired. This is
effectively the first step in the RMG rate-based mechanism generation algorithm.

The input file is a subset of that used with regular RMG jobs. 
"""

import logging
import os.path

from rmgpy.chemkin import ChemkinWriter
from rmgpy.rmg.main import initialize_log, RMG
from rmgpy.rmg.output import OutputHTMLWriter
from rmgpy.util import parse_command_line_arguments


def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.
    """
    # Parse the command-line arguments (requires the argparse module)
    args = parse_command_line_arguments()

    # Initialize the logging system (resets the RMG.log file)
    level = logging.INFO
    if args.debug:
        level = 0
    elif args.verbose:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING

    kwargs = {
        'restart': args.restart,
        'walltime': args.walltime,
        'log': level,
        'maxproc': args.maxproc,
        'kineticsdatastore': args.kineticsdatastore
    }

    initialize_log(level, os.path.join(args.output_directory, 'RMG.log'))

    rmg = RMG(input_file=args.file, output_directory=args.output_directory)

    # Add output listeners:
    rmg.attach(ChemkinWriter(args.output_directory))
    rmg.attach(OutputHTMLWriter(args.output_directory))

    execute(rmg, **kwargs)


def execute(rmg, **kwargs):
    """

    Generates all the possible reactions involving a given
    set of reactants, including pressure-dependent effects if desired. This is
    effectively the first step in the RMG rate-based mechanism generation algorithm.

    The input file is a subset of that used with regular RMG jobs. 

    Returns an RMG object.
    """
    rmg.initialize(**kwargs)

    rmg.reaction_model.enlarge(react_edge=True,
                               unimolecular_react=rmg.unimolecular_react,
                               bimolecular_react=rmg.bimolecular_react,
                               trimolecular_react=rmg.trimolecular_react)
    # Show all core and edge species and reactions in the output
    rmg.reaction_model.output_species_list.extend(rmg.reaction_model.edge.species)
    rmg.reaction_model.output_reaction_list.extend(rmg.reaction_model.edge.reactions)

    rmg.save_everything()

    rmg.finish()

    return rmg
