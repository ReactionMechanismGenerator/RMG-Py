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
RMG is an automatic chemical mechanism generator. It is awesomely awesome.
"""

import os.path
import logging
import rmgpy

from rmgpy.rmg.main import RMG, initializeLog, processProfileStats, makeProfileGraph
from rmgpy.util import parse_command_line_arguments

################################################################################


def main():
    # Parse the command-line arguments (requires the argparse module)
    args = parse_command_line_arguments()

    if args.postprocess:
        logging.info("Postprocessing the profiler statistics (will be appended to RMG.log)")
    else:
        # Initialize the logging system (resets the RMG.log file)
        level = logging.INFO
        if args.debug:
            level = 0
        elif args.verbose:
            level = logging.DEBUG
        elif args.quiet:
            level = logging.WARNING
        initializeLog(level, os.path.join(args.output_directory, 'RMG.log'))

    logging.info(rmgpy.settings.report())

    kwargs = {
        'restart': args.restart,
        'walltime': args.walltime,
        'maxproc': args.maxproc,
        'kineticsdatastore': args.kineticsdatastore
    }

    if args.profile:
        import cProfile
        global_vars = {}
        local_vars = {
            'inputFile': args.file,
            'output_dir': args.output_directory,
            'kwargs': kwargs,
            'RMG': RMG
        }

        command = """rmg = RMG(inputFile=inputFile, outputDirectory=output_dir); rmg.execute(**kwargs)"""

        stats_file = os.path.join(args.output_directory, 'RMG.profile')
        print("Running under cProfile")
        if not args.postprocess:
            # actually run the program!
            cProfile.runctx(command, global_vars, local_vars, stats_file)
        # postprocess the stats
        log_file = os.path.join(args.output_directory, 'RMG.log')
        processProfileStats(stats_file, log_file)
        makeProfileGraph(stats_file)

    else:

        rmg = RMG(inputFile=args.file, outputDirectory=args.output_directory)
        rmg.execute(**kwargs)


################################################################################

if __name__ == '__main__':
    main()
