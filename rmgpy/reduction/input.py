#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os.path

from rmgpy.tools.loader import loadRMGPyJob
from rmgpy.scoop_framework.util import logger as logging

def loadReductionInput(reductionFile):
    """
    Load an reduction job from the input file located at `reductionFile`
    """

    targets = None
    tolerance = -1

    full_path = os.path.abspath(os.path.expandvars(reductionFile))
    try:
        f = open(full_path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise e

    logging.info('Reading input file "{0}"...'.format(full_path))
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'targets': targets,
        'tolerance': tolerance
    }

    try:
        exec f in global_context, local_context

        targets = local_context['targets']
        tolerance = local_context['tolerance']

    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    assert targets is not None
    assert tolerance != -1

    return targets, tolerance

def load(rmgInputFile, reductionFile, chemkinFile, speciesDict):
    
    rmg = loadRMGPyJob(rmgInputFile, chemkinFile, speciesDict, generateImages=False, useChemkinNames=True)
    targets, tolerance = loadReductionInput(reductionFile)

    return rmg, targets, tolerance
