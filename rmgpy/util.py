#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

import argparse
import logging
import os.path
import shutil
import time
from functools import wraps


class Subject(object):
    """Subject in Observer Pattern"""

    def __init__(self):
        self._observers = []

    """
    Call this method when your (self-implemented)
    observer class should start listening to the Subject
    class.

    e.g.:

    listener = YourOwnListener()
    subject.attach(listener)
    """

    def attach(self, observer):
        if not observer in self._observers:
            self._observers.append(observer)

    """
    Call this method when your (self-implemented)
    observer class should stop listening to the Subject
    class.

    e.g.:
    listener = YourOwnListener()
    subject.attach(listener)

    ...<do some work>...

    subject.detach(listener)

    """

    def detach(self, observer):
        try:
            self._observers.remove(observer)
        except ValueError:
            pass

    """
    Call this method in classes that implement
    Subject, when the data that your interested in,
    is available.

    e.g.:
    class YourClass(Subject):
        ...
        def simulate(...)
            <stuff is being done>

            self.notify()

            <continue doing other stuff>            


    Make sure that your listener class implements the update(subject)
    method!

    e.g.:

    class YourOwnListener(object):
        def __init__(self):
            self.data = []

        def update(self, subject):
            self.data.append(subject.data)

    """

    def notify(self, modifier=None):
        for observer in self._observers:
            if modifier != observer:
                observer.update(self)


def make_output_subdirectory(output_directory, folder):
    """
    Create a subdirectory `folder` in the output directory. If the folder
    already exists (e.g. from a previous job) its contents are deleted.
    """
    dirname = os.path.join(output_directory, folder)
    if os.path.exists(dirname):
        # The directory already exists, so delete it (and all its content!)
        shutil.rmtree(dirname)
    os.mkdir(dirname)


def timefn(fn):
    @wraps(fn)
    def measure_time(*args, **kwargs):
        t1 = time.time()
        result = fn(*args, **kwargs)
        t2 = time.time()
        logging.info("@timefn: {} took {:.2f} seconds".format(fn.__name__, t2 - t1))
        return result

    return measure_time


def parse_command_line_arguments(command_line_args=None):
    """
    Parse the command-line arguments being passed to RMG Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser(description='Reaction Mechanism Generator (RMG) is an automatic chemical reaction '
                                                 'mechanism generator that constructs kinetic models composed of '
                                                 'elementary chemical reaction steps using a general understanding of '
                                                 'how molecules react.')

    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
                        help='a file describing the job to execute')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')

    # Add options for controlling what directories files are written to
    parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
                        metavar='DIR', help='use DIR as output directory')

    # Add restart option
    parser.add_argument('-r', '--restart', type=str, nargs=1, metavar='path/to/seed/', help='restart RMG from a seed',
                        default='')

    parser.add_argument('-p', '--profile', action='store_true',
                        help='run under cProfile to gather profiling statistics, and postprocess them if job completes')
    parser.add_argument('-P', '--postprocess', action='store_true',
                        help='postprocess profiling statistics from previous [failed] run; does not run the simulation')

    parser.add_argument('-t', '--walltime', type=str, nargs=1, default='00:00:00:00',
                        metavar='DD:HH:MM:SS', help='set the maximum execution time')

    parser.add_argument('-i', '--maxiter', type=int, nargs=1, default=None,
                        help='set the maximum number of RMG iterations')

    # Add option to select max number of processes for reaction generation
    parser.add_argument('-n', '--maxproc', type=int, nargs=1, default=1,
                        help='max number of processes used during reaction generation')

    # Add option to output a folder that stores the details of each kinetic database entry source
    parser.add_argument('-k', '--kineticsdatastore', action='store_true',
                        help='output a folder, kinetics_database, that contains a .txt file for each reaction family '
                             'listing the source(s) for each entry')

    args = parser.parse_args(command_line_args)

    # Process args to set correct default values and format

    # For output and scratch directories, if they are empty strings, set them
    # to match the input file location
    args.file = args.file[0]

    # If walltime was specified, retrieve this string from the element 1 list
    if args.walltime != '00:00:00:00':
        args.walltime = args.walltime[0]

    if args.restart:
        args.restart = args.restart[0]

    if args.maxiter:
        args.maxiter = args.maxiter[0]

    if args.maxproc != 1:
        args.maxproc = args.maxproc[0]

    # Set directories
    input_directory = os.path.abspath(os.path.dirname(args.file))

    if args.output_directory == '':
        args.output_directory = input_directory
    # If output directory was specified, retrieve this string from the element 1 list
    else:
        args.output_directory = args.output_directory[0]

    if args.postprocess:
        args.profile = True

    return args


def as_list(item, default=None):
    """
    Wrap the given item in a list if it is not None and not already a list.

    Args:
        item: the item to be put in a list
        default (optional): a default value to return if the item is None
    """
    if isinstance(item, list):
        return item
    elif item is None:
        return default
    else:
        return [item]
