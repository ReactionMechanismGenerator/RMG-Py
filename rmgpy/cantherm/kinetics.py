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

from rmgpy.kinetics import Arrhenius

################################################################################

def generateKineticsModel(reaction, tunneling='', plot=False):
    
    logging.info('Calculating rate coefficient for {0}...'.format(reaction))
    
    if len(reaction.reactants) == 1:
        kunits = 's^-1'
    elif len(reaction.reactants) == 2:
        kunits = 'm^3/(mol*s)'
    elif len(reaction.reactants) == 3:
        kunits = 'm^6/(mol^2*s)'
    else:
        kunits = ''
    
    Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
    klist = reaction.calculateTSTRateCoefficients(Tlist, tunneling)
    arrhenius = Arrhenius().fitToData(Tlist, klist, kunits)
    klist2 = arrhenius.getRateCoefficients(Tlist)
    
    reaction.kinetics = arrhenius
    
    if plot:
        logging.info('Plotting kinetics model for {0}...'.format(reaction))
        import pylab
        pylab.semilogy(1000.0 / Tlist, klist  * reaction.degeneracy, 'ok')
        pylab.semilogy(1000.0 / Tlist, klist2 * reaction.degeneracy, '-k')
        pylab.xlabel('1000 / Temperature (1000/K)')
        pylab.ylabel('Rate coefficient (SI units)')
        pylab.show()
    
################################################################################

def saveKinetics(reaction, tunneling, label, path):
    """
    Append the kinetic model generated for `reaction` with associated
    string `label` to the file located at `path` on disk.
    """
    
    f = open(path, 'a')

    Aunits = reaction.kinetics.A.units
    
    f.write('#   ======= =========== =========== =========== ===============\n')
    f.write('#   Temp.   k (TST)     Tunneling   k (TST+T)   Units\n')
    f.write('#   ======= =========== =========== =========== ===============\n')
    for T in [300,400,500,600,800,1000,1500,2000]:  
        k0 = reaction.calculateTSTRateCoefficient(T, tunneling='')
        if tunneling.lower() == 'wigner':
            kappa = reaction.calculateWignerTunnelingCorrection(T)
        elif tunneling.lower() == 'eckart':
            kappa = reaction.calculateEckartTunnelingCorrection(T)
        else:
            kappa = 1.0
        k = k0 * kappa
        f.write('#    {0:4g} K {1:11.3e} {2:11.3f} {3:11.3e} {4}\n'.format(T, k0, kappa, k, Aunits))
    f.write('#   ======= =========== =========== =========== ===============\n')
    f.write('kinetics(\n')
    f.write('    label = "{0}",\n'.format(label))
    
    # Reaction path degeneracy is INCLUDED in the kinetics itself!
    if isinstance(reaction.kinetics, Arrhenius):
        f.write('    kinetics = {0!r},\n'.format(reaction.kinetics))
    f.write('    Tmin = (300.0,"K"),\n')
    f.write('    Tmax = (2000.0,"K"),\n')
    f.write('    degeneracy = {0:d},\n'.format(reaction.degeneracy))
    f.write('    short_comment = "",\n')
    f.write('    long_comment = \n')
    f.write('"""\n')
    f.write('\n')
    f.write('""",\n')
    f.write(')\n')
    f.write('\n')
    
    f.close()
