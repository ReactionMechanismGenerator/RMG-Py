#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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

import sys
import os.path

from input import load
from output import write_model
from optimization import optimize
from reduction import compute_conversion, initialize

from rmgpy.scoop_framework.util import logger as logging

def main():
    
    inputFile, reductionFile, chemkinFile, spc_dict = sys.argv[-4:]

    for f in [inputFile, reductionFile, chemkinFile, spc_dict]:
        assert os.path.isfile(f), 'Could not find {}'.format(f)

    inputDirectory = os.path.abspath(os.path.dirname(inputFile))
    output_directory = inputDirectory

    rmg, target_label, error = load(inputFile, reductionFile, chemkinFile, spc_dict)
    logging.info('Allowed error in target conversion: {0:.0f}%'.format(error * 100))

    reactionModel = rmg.reactionModel
    initialize(rmg.outputDirectory, reactionModel.core.reactions)

    atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
    index = 0
    reactionSystem = rmg.reactionSystems[index]    
    
    #compute original target conversion
    Xorig = compute_conversion(target_label, reactionModel, reactionSystem, index,\
     rmg.absoluteTolerance, rmg.relativeTolerance)

    logger.info('Observables of original model:')
    for target, observable in zip(targets, observables):
        logger.info('{}: {:.2f}%'.format(target, observable * 100))

    # optimize reduction tolerance
    tol, important_reactions = optimize(target_label, reactionModel, rmg, index, error, Xorig)
    logging.info('Optimized tolerance: {:.0E}'.format(tol))

    # plug the important reactions into the RMG object and write:
    rmg.reactionModel.core.reactions = important_reactions
    write_model(rmg)


if __name__ == '__main__':
    main()