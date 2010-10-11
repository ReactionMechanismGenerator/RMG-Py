#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   CanTherm
#    
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

import numpy
import logging

from chempy.kinetics import ArrheniusModel

################################################################################

def generateKineticsModel(reaction, tunneling='', plot=False):
    
    logging.info('Calculating rate coefficient for %s...' % (reaction))
    
    Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
    klist = reaction.calculateTSTRateCoefficients(Tlist, tunneling)
    arrhenius = ArrheniusModel().fitToData(Tlist, klist)
    klist2 = arrhenius.getRateCoefficients(Tlist)
    
    logging.debug('    %s' % (arrhenius))
    
    if plot:
        logging.info('Plotting thermo model for %s...' % (reaction))
        import pylab
        pylab.semilogy(1000.0 / Tlist, klist , 'ok')
        pylab.semilogy(1000.0 / Tlist, klist2, '-k')
        pylab.xlabel('1000 / Temperature (1000/K)')
        pylab.ylabel('Rate coefficient (SI units)')
        pylab.show()
    
    return arrhenius
    