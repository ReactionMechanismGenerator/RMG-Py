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

"""
This module contains functionality for the parallel execution of RMG-Py.
"""

import sys
import traceback
import warnings
from functools import wraps

try:
    from scoop import futures
    from scoop import shared
    from scoop import logger as logging

except ImportError:
    import logging as logging
    logging.debug("Could not properly import SCOOP.")

def warnScoopStartedProperly(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        
        futures_not_loaded = 'scoop.futures' not in sys.modules

        try:
            controller_not_started = not (
                sys.modules['scoop.futures'].__dict__.get("_controller", None)
            )
        except KeyError, e:
            warnings.warn(
                    "SCOOP was not started properly.\n"
                    "Be sure to start your program with the "
                    "'-m scoop' parameter. You can find "
                    "further information in the "
                    "documentation.\n",
                    RuntimeWarning
                )
            return

        if futures_not_loaded or controller_not_started:
            warnings.warn(
                    "SCOOP was not started properly.\n"
                    "Be sure to start your program with the "
                    "'-m scoop' parameter. You can find "
                    "further information in the "
                    "documentation.\n",
                    RuntimeWarning
                )

            return

        return func(*args, **kwargs)
    return wrapper

class WorkerWrapper(object):
    """
    This class can be used to expose the exception trace of a worker
    that was running on a remote worker.

    Use it as follows:

    Wrap the function that will be running on the remote worker with the current class:

    futures.map(WorkerWrapper(f), mapped_data, ...) or 
    futures.submit(WorkerWrapper(f), mapped_data, ...)

    """
    __name__ = 'WorkerWrapper'

    def __init__(self, myfn):
        self.myfn = myfn

    def __call__(self, *args, **kwargs):
        try:
            return self.myfn(*args, **kwargs)
        except:
            type, value, tb = sys.exc_info()
            lines = traceback.format_exception(type, value, tb)
            print ''.join(lines)
            raise

@warnScoopStartedProperly
def broadcast(obj, key):
    """
    Broadcasts the object across the workers using the key parameter as the key.
    """      
    
    kwargs = {key : obj}
    try:
        if shared.getConst(key):
            logging.debug('An object with the key {} was already broadcasted.'.format(key))
        else:
            shared.setConst(**kwargs)
    except NameError, e:
        """
        Name error will be caught when the SCOOP library is not imported properly.
        """
        logging.debug('SCOOP not loaded. Not broadcasting the object {}'.format(obj))

@warnScoopStartedProperly
def get(key):    
    """
    Searches for the shared variable to retrieve identified by the 
    parameter key.
    """

    try:
        data = shared.getConst(key)
        return data
    except KeyError, e:
        logging.error('An object with the key {} could not be found.'.format(key))
        raise e
    except NameError:
        """
        Name error will be caught when the SCOOP library is not imported properly.
        """
        logging.debug('SCOOP not loaded. Not retrieving the shared object with key {}'.format(key))