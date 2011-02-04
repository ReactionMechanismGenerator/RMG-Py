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

from rmgpy.chem.kinetics import ArrheniusModel

################################################################################

def generateKineticsModel(reaction, tunneling='', plot=False):
    
    logging.info('Calculating rate coefficient for %s...' % (reaction))
    
    Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
    klist = reaction.calculateTSTRateCoefficients(Tlist, tunneling)
    arrhenius = ArrheniusModel().fitToData(Tlist, klist)
    klist2 = arrhenius.getRateCoefficients(Tlist)
    
    reaction.kinetics = arrhenius
    
    if plot:
        logging.info('Plotting thermo model for %s...' % (reaction))
        import pylab
        pylab.semilogy(1000.0 / Tlist, klist  * reaction.degeneracy, 'ok')
        pylab.semilogy(1000.0 / Tlist, klist2 * reaction.degeneracy, '-k')
        pylab.xlabel('1000 / Temperature (1000/K)')
        pylab.ylabel('Rate coefficient (SI units)')
        pylab.show()
    
################################################################################

def saveKinetics(reaction, label, path):
    """
    Append the kinetic model generated for `reaction` with associated
    string `label` to the file located at `path` on disk.
    """
    
    f = open(path, 'a')

    Nreac = len(reaction.reactants)
    if Nreac == 1:   Aunits = 's^-1'
    elif Nreac == 2: Aunits = 'cm^3/(mol*s)'
    elif Nreac == 3: Aunits = 'cm^6/(mol^2*s)'
    else:            Aunits = 'cm^%i/(mol^%i*s)' % (3*(Nreac-1), Nreac-1)

    f.write('# Kinetics for %s:\n' % (label)) # Note: This includes reaction-path degeneracy!
    f.write('#   k(300 K)  = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(300, 1e5) * 1.0e6**(Nreac - 1), Aunits))
    f.write('#   k(400 K)  = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(400, 1e5) * 1.0e6**(Nreac - 1), Aunits))
    f.write('#   k(500 K)  = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(500, 1e5) * 1.0e6**(Nreac - 1), Aunits))
    f.write('#   k(600 K)  = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(600, 1e5) * 1.0e6*(Nreac - 1), Aunits))
    f.write('#   k(800 K)  = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(800, 1e5) * 1.0e6**(Nreac - 1), Aunits))
    f.write('#   k(1000 K) = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(1000, 1e5) * 1.0e6**(Nreac - 1), Aunits))
    f.write('#   k(1500 K) = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(1500, 1e5) * 1.0e6**(Nreac - 1), Aunits))
    f.write('#   k(2000 K) = %12.3e %s\n' % (reaction.degeneracy * reaction.kinetics.getRateCoefficient(2000, 1e5) * 1.0e6**(Nreac - 1), Aunits))

    Nreac = len(reaction.reactants)
    if Nreac == 1:   Aunits = 's^-1'
    elif Nreac == 2: Aunits = 'm^3/(mol*s)'
    elif Nreac == 3: Aunits = 'm^6/(mol^2*s)'
    else:            Aunits = 'm^%i/(mol^%i*s)' % (3*(Nreac-1), Nreac-1)

    f.write('kinetics(\n')
    f.write('    label = "%s",\n' % label)
    
    # Reaction path degeneracy is NOT included in the model itself; you must
    # add it yourself later
    if isinstance(reaction.kinetics, ArrheniusModel):
        f.write('    kinetics = ArrheniusModel(\n')
        f.write('        A = (%g, "%s"),\n' % (reaction.kinetics.A, Aunits))
        f.write('        n = %g,\n' % (reaction.kinetics.n))
        f.write('        Ea = (%g, "kcal/mol"),\n' % (reaction.kinetics.Ea / 4184))
        f.write('        T0 = (%g, "K"),\n' % (reaction.kinetics.T0))
        f.write('    ),\n')
    
    f.write('    Tmin = 300.0,\n')
    f.write('    Tmax = 2000.0,\n')
    f.write('    degeneracy = %i,\n' % int(reaction.degeneracy))
    f.write('    short_comment = "",\n')
    f.write('    long_comment = \n')
    f.write('"""\n')
    f.write('\n')
    f.write('""",\n')
    f.write(')\n')
    f.write('\n')
    
    f.close()
