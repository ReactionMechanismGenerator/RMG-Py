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

import math
import numpy.linalg
import logging

import chempy.constants as constants
from chempy.thermo import ThermoGAModel, WilhoitModel
from chempy.ext.thermo_converter import convertWilhoitToNASA

################################################################################

def generateThermoModel(species, model, plot=False):
    """
    Generate a thermodynamics model for a given `species`. The type of `model`
    generated can be any of ``'Group additivity'``, ``'Wilhoit'``, or 
    ``'NASA'``. Set `plot` to ``True`` to see a plot of the generated
    thermodynamic values.
    """

    if model.lower() not in ['group additivity', 'wilhoit', 'nasa']:
        raise Exception('Unknown thermodynamic model "%s".' % model)

    logging.info('Generating %s thermo model for %s...' % (model, species))
    linear = species.states.modes[1].linear
    Nfreq = len(species.states.modes[2].frequencies)
    Nrotors = len(species.states.modes[3:])

    H298 = species.states.getEnthalpy(298.15) + species.E0
    S298 = species.states.getEntropy(298.15)
    
    if model.lower() == 'group additivity':
        Tdata = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)
        Cpdata = species.states.getHeatCapacities(Tdata)
        species.thermo = ThermoGAModel(Tdata=Tdata, Cpdata=Cpdata, H298=H298, S298=S298)
    else:
        Tlist = numpy.arange(10.0, 3001.0, 10.0, numpy.float64)
        Cplist = species.states.getHeatCapacities(Tlist)
        wilhoit = WilhoitModel()
        wilhoit.fitToData(Tlist, Cplist, linear, Nfreq, Nrotors, H298, S298, B0=500.0)
        if model.lower() == 'nasa':
            species.thermo = convertWilhoitToNASA(wilhoit, Tmin=10.0, Tmax=3000.0, Tint=500.0, fixedTint=False, weighting=True, continuity=3)
        else:
            species.thermo = wilhoit
            
    # Plots to compare with the states model predictions
    if plot:
        print 'Plotting thermo model for %s...' % (species)

        import pylab
        Tlist = numpy.arange(10.0, 2501.0, 10.0)
        Cplist = species.states.getHeatCapacities(Tlist)
        Cplist1 = species.thermo.getHeatCapacities(Tlist)
        Slist = species.states.getEntropies(Tlist)
        Slist1 = species.thermo.getEntropies(Tlist)
        Hlist = species.states.getEnthalpies(Tlist) + species.E0
        Hlist1 = species.thermo.getEnthalpies(Tlist)
        Glist = Hlist - Tlist * Slist
        Glist1 = species.thermo.getFreeEnergies(Tlist)

        fig = pylab.figure(figsize=(10,8))

        pylab.subplot(2,2,1)
        pylab.plot(Tlist, Cplist / 4.184, '-r', Tlist, Cplist1 / 4.184, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Heat capacity (cal/mol*K)')
        pylab.legend(['states', 'thermo'], loc=4)

        pylab.subplot(2,2,2)
        pylab.plot(Tlist, Slist / 4.184, '-r', Tlist, Slist1 / 4.184, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Entropy (cal/mol*K)')

        pylab.subplot(2,2,3)
        pylab.plot(Tlist, Hlist / 4184.0, '-r', Tlist, Hlist1 / 4184.0, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Enthalpy (kcal/mol)')

        pylab.subplot(2,2,4)
        pylab.plot(Tlist, Glist / 4184.0, '-r', Tlist, Glist1 / 4184.0, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Gibbs free energy (kcal/mol)')

        fig.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.95, wspace=0.35, hspace=0.20)
        pylab.show()

################################################################################

def saveThermoData(fstr, species, label):

    from datetime import date
    f = open(fstr, 'w')
    f.write('Thermo(\n')
    f.write('    species = \n')
    f.write('"""\n')
    f.write('%s\n' % label)
    f.write('\n')
    f.write('"""\n')
    f.write('    thermo = %s,\n' % repr(species.thermo))
    f.write('    Tmin = 0.0,\n')
    f.write('    Tmax = 9999.9,\n')
    f.write('    short_comment = "",\n')
    f.write('    long_comment = \n')
    f.write('"""\n')
    f.write('\n')
    f.write('"""\n')
    f.write('\n')
    f.write('"""\n')
    f.write('    history = ("%s","jwallen","Added to database.)",\n' % date.today().strftime('%Y-%m-%d'))
    f.write(')\n')
    f.write('\n')
    f.write('States(\n')
    f.write('    species = \n')
    f.write('"""\n')
    f.write('%s\n' % label)
    f.write('\n')
    f.write('"""\n')
    f.write('    modes = [\n')
    for mode in species.states.modes:
        f.write('        %s,\n' % repr(mode))
    f.write('    ],\n')
    f.write('    short_comment = "",\n')
    f.write('    long_comment = \n')
    f.write('"""\n')
    f.write('\n')
    f.write('"""\n')
    f.write('    history = ("%s","jwallen","Added to database.)",\n' % date.today().strftime('%Y-%m-%d'))
    f.write(')')
    f.close()
