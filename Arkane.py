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
This is the main executable script for Arkane, a tool for computing chemical
reaction rates and other properties used in detailed kinetics models using
various methodologies and theories. To run Arkane, use the command ::

    $ python Arkane.py FILE

where ``FILE`` is the path to an Arkane input file describing the job to
execute. Arkane will run the specified job, writing the output to
``output.py`` and a log to both the console and to ``Arkane.log``, with both
files appearing in the same directory as the input file. Some additional
command-line arguments are available; run the command ::

    $ python Arkane.py -h

for more information.
"""

import os
import logging

from arkane.main import Arkane

arkane = Arkane()

# Parse and validate the command-line arguments
arkane.parse_command_line_arguments()

# Execute the job
arkane.execute()

try:
    import psutil

    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    logging.info('Memory used: %.2f MB' % (memory_info.rss / 1024.0 / 1024.0))
except ImportError:
    logging.info('Optional package dependency "psutil" not found; memory profiling information will not be saved.')
