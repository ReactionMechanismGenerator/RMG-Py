#!/usr/bin/env python

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
RMG is an automatic chemical mechanism generator. It is awesomely awesome.
"""

import os.path
import logging

# Before importing any RMG modules, check Python version
try:
    import utilities
except ImportError:  # This is likely the binary install version, which does not include utilities
    # If this is in fact the binary version, conda has ensured that the python version is correct.
    # It is okay to skip this test here
    utilities = None
else:
    utilities.check_python()

import rmgpy
from rmgpy.rmg.main import RMG, initialize_log, process_profile_stats, make_profile_graph
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
        initialize_log(level, os.path.join(args.output_directory, 'RMG.log'))

    logging.info(rmgpy.settings.report())

    kwargs = {
        'restart': args.restart,
        'walltime': args.walltime,
        'maxproc': args.maxproc,
        'kineticsdatastore': args.kineticsdatastore,
        'max_iterations': args.maxiter,
    }

    if args.profile:
        import cProfile
        rmg_profiler = cProfile.Profile()
        stats_file = os.path.join(args.output_directory, 'RMG.profile')

        global_vars = {}
        local_vars = {
            'inputFile': args.file,
            'output_dir': args.output_directory,
            'kwargs': kwargs,
            'RMG': RMG,
            'rmg_profiler': rmg_profiler,
        }

        command = """rmg = RMG(input_file=inputFile, output_directory=output_dir, profiler=rmg_profiler); rmg.execute(**kwargs)"""

        print("Running under cProfile")
        if not args.postprocess:
            # actually run the program!
            rmg_profiler.runctx(command, global_vars, local_vars)
            rmg_profiler.dump_stats(stats_file)
        # postprocess the stats
        log_file = os.path.join(args.output_directory, 'RMG.log')
        process_profile_stats(stats_file, log_file)

        # Make the profile graph. Force graph generation regardless of display status if args.postprocess
        make_profile_graph(stats_file, force_graph_generation=args.postprocess)
        
    else:

        rmg = RMG(input_file=args.file, output_directory=args.output_directory)
        rmg.execute(**kwargs)


################################################################################

if __name__ == '__main__':
    main()
