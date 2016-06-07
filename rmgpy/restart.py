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


import os.path
import logging
import cPickle
import time

def save(rmg):
    # Save the restart file if desired
    if rmg.saveRestartPeriod or rmg.done:
        saveRestartFile( os.path.join(rmg.outputDirectory,'restart.pkl'),
                              rmg.reactionModel,
                              delay=0 if rmg.done else rmg.saveRestartPeriod.value_si
                            )

        
        
def saveRestartFile(path, reactionModel, delay=0):
    """
    Save a restart file to `path` on disk containing the contents of the
    provided `reactionModel`. The `delay` parameter is a time in seconds; if
    the restart file is not at least that old, the save is aborted. (Use the
    default value of 0 to force the restart file to be saved.)
    """
    
    # Saving of a restart file is very slow (likely due to all the Quantity objects)
    # Therefore, to save it less frequently, don't bother if the restart file is less than an hour old
    if os.path.exists(path) and time.time() - os.path.getmtime(path) < delay:
        logging.info('Not saving restart file in this iteration.')
        return
    
    # Pickle the reaction model to the specified file
    # We also compress the restart file to save space (and lower the disk read/write time)
    logging.info('Saving restart file...')
    f = open(path, 'wb')
    cPickle.dump(reactionModel, f, cPickle.HIGHEST_PROTOCOL)
    f.close()

class RestartWriter(object):
    """
    This class listens to a RMG subject
    and writes a restart file with the current state of the RMG model.


    A new instance of the class can be appended to a subject as follows:
    
    rmg = ...
    listener = RestartWriter()
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)
    
    """
    def __init__(self):
        super(RestartWriter, self).__init__()
    
    def update(self, rmg):
      	save(rmg)

        
    