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

"""
This module contains functionality for the parallel execution of RMG-Py.
"""

import sys
import traceback
import warnings
from functools import wraps

logger = None

try:
    from scoop import futures
    from scoop.futures import map, submit
    from scoop import shared
    from scoop import logger as scooplogger
    logger = scooplogger
    # logger.setLevel(20)#10 : debug, 20: info
except ImportError:
    import logging as logging
    logger = logging.getLogger()
    logger.debug("Could not properly import SCOOP.")

def warnScoopStartedProperly(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        
        futures_not_loaded = 'scoop.futures' not in sys.modules

        warnings.simplefilter('ignore', RuntimeWarning)

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
            logger.debug('An object with the key {} was already broadcasted.'.format(key))
        else:
            shared.setConst(**kwargs)
    except NameError, e:
        """
        Name error will be caught when the SCOOP library is not imported properly.
        """
        logger.debug('SCOOP not loaded. Not broadcasting the object {}'.format(obj))

@warnScoopStartedProperly
def get(key):    
    """
    Searches for the shared variable to retrieve identified by the 
    parameter key.
    """

    try:
        data = shared.getConst(key, timeout=1e-9)
        return data
    except NameError:
        """
        Name error will be caught when the SCOOP library is not imported properly.
        """
        logger.debug('SCOOP not loaded. Not retrieving the shared object with key {}'.format(key))

def map_(*args, **kwargs):
    return map(WorkerWrapper(args[0]), *args[1:], **kwargs)

def submit_(func, *args, **kwargs):
    """
    Task submission of a function.

    returns the return value of the called function, or
    when SCOOP is loaded, the future object.
    """
    try:
        task = submit(WorkerWrapper(func), *args, **kwargs)#returns immediately
        return task
    except Exception, e:
        """
        Name error will be caught when the SCOOP library is not imported properly.
        """
        logger.debug('SCOOP not loaded. Submitting serial mode.')
        return func(*args, **kwargs)
