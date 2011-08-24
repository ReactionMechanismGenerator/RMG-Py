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

from rmgpy.quantity import Quantity
from rmgpy.thermo import ThermoData, Wilhoit, MultiNASA, convertThermoModel
from rmgpy.statmech import RigidRotor

################################################################################

def generateThermoModel(species, model, plot=False):
    """
    Generate a thermodynamics model for a given `species`. The type of `model`
    generated can be any of ``'Group additivity'``, ``'Wilhoit'``, or 
    ``'NASA'``. Set `plot` to ``True`` to see a plot of the generated
    thermodynamic values.
    """

    if model.lower() not in ['group additivity', 'wilhoit', 'nasa']:
        raise Exception('Unknown thermodynamic model "{0}".'.format(model))

    logging.info('Generating {0} thermo model for {1}...'.format(model, species))
    if not any([isinstance(mode, RigidRotor) for mode in species.states.modes]):
        linear = False
        Nfreq = 0
        Nrotors = 0
    else:
        linear = species.states.modes[1].linear
        Nfreq = len(species.states.modes[2].frequencies.values)
        Nrotors = len(species.states.modes[3:])

    H298 = species.states.getEnthalpy(298.15) + species.E0.value
    S298 = species.states.getEntropy(298.15)
    
    if model.lower() == 'group additivity':
        Tdata = numpy.arange(300.0, 2001.0, 100.0, numpy.float64)
        Cpdata = species.states.getHeatCapacities(Tdata)
        species.thermo = ThermoData(Tdata=(Tdata,"K"), Cpdata=(Cpdata,"J/(mol*K)"), H298=(H298/1000.,"kJ/mol"), S298=(S298,"J/(mol*K)"))
    else:
        Tlist = numpy.arange(10.0, 3001.0, 10.0, numpy.float64)
        Cplist = species.states.getHeatCapacities(Tlist)
        wilhoit = Wilhoit()
        if Nfreq == 0 and Nrotors == 0:
            wilhoit.cp0 = Quantity(Cplist[0],"J/(mol*K)")
            wilhoit.cpInf = Quantity(Cplist[0],"J/(mol*K)")
            wilhoit.B = Quantity(500,"K")
            wilhoit.H0 = Quantity(0,"J/mol")
            wilhoit.S0 = Quantity(0,"J/(mol*K)")
            wilhoit.H0 = Quantity(species.states.getEnthalpy(298.15) - wilhoit.getEnthalpy(298.15) + species.E0.value,"J/mol")
            wilhoit.S0 = Quantity(species.states.getEntropy(298.15) - wilhoit.getEntropy(298.15),"J/(mol*K)")
        else:
            wilhoit.fitToData(Tlist, Cplist, linear, Nfreq, Nrotors, H298, S298, B0=500.0)
        if model.lower() == 'nasa':
            species.thermo = convertThermoModel(wilhoit, MultiNASA, Tmin=10.0, Tmax=3000.0, Tint=500.0)
        else:
            species.thermo = wilhoit
            
    # Plots to compare with the states model predictions
    if plot:
        print 'Plotting thermo model for {0}...'.format(species)

        import pylab
        Tlist = numpy.arange(10.0, 2501.0, 10.0)
        Cplist = species.states.getHeatCapacities(Tlist)
        Cplist1 = species.thermo.getHeatCapacities(Tlist)
        Slist = species.states.getEntropies(Tlist)
        Slist1 = species.thermo.getEntropies(Tlist)
        Hlist = species.states.getEnthalpies(Tlist) + species.E0.value
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

def saveThermo(species, label, path):
    """
    Append the thermodynamic model generated for `species` with associated
    string `label` to the file located at `path` on disk.
    """
    
    f = open(path, 'a')

    f.write('# Thermodynamics for %s:\n' % (label))
    f.write('#   H(298 K)   = {0:9.3f} kcal/mol\n'.format(species.thermo.getEnthalpy(298) / 4184.))
    f.write('#   S(298 K)   = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getEntropy(298) / 4.184))
    f.write('#   Cp(300 K)  = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getHeatCapacity(300) / 4.184))
    f.write('#   Cp(400 K)  = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getHeatCapacity(400) / 4.184))
    f.write('#   Cp(500 K)  = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getHeatCapacity(500) / 4.184))
    f.write('#   Cp(600 K)  = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getHeatCapacity(600) / 4.184))
    f.write('#   Cp(800 K)  = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getHeatCapacity(800) / 4.184))
    f.write('#   Cp(1000 K) = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getHeatCapacity(1000) / 4.184))
    f.write('#   Cp(1500 K) = {0:9.3f} cal/(mol*K)\n'.format(species.thermo.getHeatCapacity(1500) / 4.184))

    f.write('thermo(\n')
    f.write('    label = "%s",\n' % label)
    
    if isinstance(species.thermo, ThermoData):
        f.write('    thermo = ThermoData(\n')
        f.write('        Tdata = {0!r},\n'.format(species.thermo.Tdata))
        f.write('        Cpdata = {0!r},\n'.format(species.thermo.Cpdata))
        f.write('        H298 = {0!r},\n'.format(species.thermo.H298))
        f.write('        S298 = {0!r},\n'.format(species.thermo.S298))
        f.write('    ),\n')
    elif isinstance(species.thermo, Wilhoit):
        f.write('    thermo = Wilhoit(\n')
        f.write('        cp0 = {0!r},\n'.format(species.thermo.cp0))
        f.write('        cpInf = {0!r},\n'.format(species.thermo.cpInf))
        f.write('        B = {0!r},\n'.format(species.thermo.B))
        f.write('        a0 = {0:g},\n'.format(species.thermo.a0))
        f.write('        a1 = {0:g},\n'.format(species.thermo.a1))
        f.write('        a2 = {0:g},\n'.format(species.thermo.a2))
        f.write('        a3 = {0:g},\n'.format(species.thermo.a3))
        f.write('        H0 = {0!r},\n'.format(species.thermo.H0))
        f.write('        S0 = {0!r},\n'.format(species.thermo.S0))
        f.write('    ),\n')
    elif isinstance(species.thermo, MultiNASA):
        f.write('    thermo = MultiNASA(polynomials = [\n')
        for poly in species.thermo.polynomials:
            f.write('        {0!r},\n'.format(poly))
        f.write('    ]),\n')
       
    f.write('    Tmin = (0.0,"K"),\n')
    f.write('    Tmax = (3000.0,"K"),\n')
    f.write('    short_comment = "",\n')
    f.write('    long_comment = \n')
    f.write('"""\n')
    f.write('\n')
    f.write('""",\n')
    f.write(')\n')
    f.write('\n')
    
    f.close()
