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

import numpy as np

from reduction import reduceModel
from output import writeModel
from rmgpy.scoop_framework.util import logger as logging

def optimize(target_label, reactionModel, rmg, reactionSystemIndex, error, orig_observable):
    """
    The optimization algorithm that searches for the most reduced model that satisfies the
    applied constraints.

    The error introduced by the reduced model for a response variable
    of a target is used as the objective function.

    The optimization algorithm increments the trial tolerance from a very low value
    until the introduced error is greater than the user-provided threshold.
    """


    low = -30 
    high = 0

    """
    Tolerance to decide whether a reaction is unimportant for the formation/destruction of a species

    Tolerance is a floating point value between 0 and 1.

    A high tolerance means that many reactions will be deemed unimportant, and the reduced model will be drastically
    smaller.

    A low tolerance means that few reactions will be deemed unimportant, and the reduced model will only differ from the full
    model by a few reactions.
    """

    tol, importantReactions = \
     bisect(low, high, error, target_label, reactionModel, rmg, reactionSystemIndex, orig_observable)

    return tol, importantReactions

def computeDeviation(original, reduced, targets):
    """
    Computes the relative deviation between the observables of the
    original and reduced model.

    Assumes the observables are numpy arrays.
    """
    devs = np.abs((reduced - original) / original)

    logging.info('Deviations: '.format())
    for dev, target in zip(devs, targets):
        logging.info('Deviation for {}: {:.2f}%'.format(target, dev * 100))

    return devs

def isInvalid(devs, error):
    """
    Check if the reduced observables differ from the original
    observables more than the parameter error threshold.
    """
    invalid = np.any(devs > error)
    return invalid

def bisect(low, high, error, targets, reactionModel, rmg, reactionSystemIndex, orig_observable):
    """
    Bisect method in log space.

    Interrupt iterations when two consecutive, successful iterations differ less than a
    threshold value.
    """

    THRESHOLD = 0.05

    importantReactions = None
    final_devs = None
    old_trial = low
    while True:
        midpoint = (low + high) / 2.0
        reduced_observable, newImportantReactions = evaluate(midpoint, targets, reactionModel, rmg, reactionSystemIndex)
        
        devs = computeDeviation(orig_observable, reduced_observable, targets)

        if isInvalid(devs, error):
            high = midpoint
        else:
            if len(newImportantReactions) == 0:
                logging.error('Model reduction resulted in a model with 0 reactions.')
                logging.error('Perhaps change reactor conditions to allow for more adequate reduction. Exiting...')
                break
            low = midpoint
            importantReactions = newImportantReactions
            final_devs = devs
            writeModel(rmg, chemkin_name='chem_reduced_{}.inp'.format(len(importantReactions)))
            
        if np.abs((midpoint - old_trial) / old_trial) < THRESHOLD:
            break

        old_trial = low

    if not importantReactions:
        logging.error("Could not find a good guess...")
        importantReactions = []

    logging.info('Final deviations: '.format())
    for dev, target in zip(final_devs, targets):
        logging.info('Final deviation for {}: {:.2f}%'.format(target, dev * 100))


    return low, importantReactions

def evaluate(guess, targets, reactionModel, rmg, reactionSystemIndex):
    """
    
    """
    logging.info('Trial tolerance: {:.2E}'.format(10**guess))

    observable, newImportantReactions = reduceModel(10**guess, targets, reactionModel, rmg, reactionSystemIndex)

    return observable, newImportantReactions
